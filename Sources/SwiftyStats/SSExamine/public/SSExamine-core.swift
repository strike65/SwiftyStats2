//
//  Created by VT on 28.10.25.
//  Copyright © 2025 Volker Thieme. All rights reserved.
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.
//

/// Core SSExamine container implementation with locking and metadata management.

import SwiftyStatsPrelude
import ZIPFoundation
import SwiftyBoost

/// A generic, univariate sample container with frequency tracking and sequence preservation.
///
/// SSExamine collects elements, tracks absolute and relative frequencies, and can reconstruct
/// the original insertion order (the “raw” sequence). It supports categorical and numeric data
/// and carries a level-of-measurement hint to guide downstream statistical operations.
///
/// Type parameters:
/// - SSElement: The element type for observations. Must be `Hashable` and `Comparable`
///   to support frequency maps and sorted views, and `Codable`/`Sendable` for persistence/concurrency.
/// - FP: The floating-point type used for relative frequencies and numeric parameters (e.g., alpha).
///   Must conform to `RealLike` (see SwiftyBoostPrelude) and be `Codable`/`Sendable`.
///
/// Design highlights:
/// - Reference semantics (class): mutations affect all references.
/// - Preserves insertion order via an internal index map, enabling `.original` reconstruction.
/// - Frequency maps are updated incrementally on append; derived caches can be invalidated with `needsUpdateFrequencies`.
/// - Codable: encodes metadata and the raw sequence for lossless round-tripping.
/// - Concurrency aware: all public entry points synchronise via an internal recursive lock so readers/writers can run in parallel tasks.
///
/// Thread-safety:
/// - An `NSRecursiveLock` guards state. Helpers snapshot under lock to avoid leaking locked references, and the type is marked `@unchecked Sendable` to advertise safe sharing despite inheriting from `NSObject`.
///
/// Equality semantics:
/// - Two instances are equal if their raw sequences are identical (order and values).
///
/// Copy semantics:
/// - `copy()` produces a deep copy of the sequence and metadata into a new instance.
///
/// Performance notes:
/// - `append(_:)`: amortized O(1).
/// - `itemsAsArray(sorted:)`: O(n log n) for sorted orders; O(n) for `.original` and `.none`. Results are produced from a temporary snapshot to sidestep locking exposure.
/// - `remove(_:allOccurences:)`: O(n) due to sequence rebuild; invoked under the lock for atomicity.
///
/// Examples:
/// ```swift
/// // Numeric example
/// let ex = SSExamine<Int, Double>()
/// ex.append(1)
/// ex.append(repeating: 2, item: 3) // adds 3 twice
/// let raw = ex.itemsAsArray(sorted: .original) // [1, 3, 3]
/// let freq3 = ex.frequency(3)             // 2
/// let rFreq3 = ex.relativeFrequency(3)    // 2/3
///
/// // Categorical (characters from a string)
/// let cs = CharacterSet.letters
/// let chars = try SSExamine<String, Double>(using: "abac", levelOfMeasurement: .nominal, name: "letters", characterSet: cs)
/// let sortedAsc = chars.itemsAsArray(sorted: .ascending) // ["a", "a", "b", "c"]
/// ```
///
public class SSExamine<SSElement, FP>: NSObject, NSCopying, Codable
where SSElement: Hashable & Comparable & Codable & Sendable, FP: RealLike {

    
    /// Convenience alias for the element type required by `ExamineContainer`.
    // typealias ExamineItem = SSElement
    
    /// Convenience alias for the floating-point frequency type required by `ExamineContainer`.
    // typealias Frequency = FP
    
    /// A recursive lock guarding all internal mutable state.
    internal let lock = NSRecursiveLock()
    /// Optional external object to associate with this instance (for UI or metadata).
    internal var rootObjectStorage: Any?
    /// Arbitrary identifier.
    internal var tagStorage: String?
    /// Human-readable name.
    internal var nameStorage: String?
    /// Level of measurement associated with the data (nominal/ordinal/interval/ratio).
    internal var levelOfMeasurementStorage: LevelOfMeasurement = .nominal
    /// A descriptive string useful for debugging or UI.
    internal var descriptionStringStorage: String = "SSExamine"
    /// Significance level for downstream statistical procedures.
    internal var alphaStorage: FP = FP(0.05)
    /// Tracks whether the dataset has been modified since the last stable state.
    internal var hasChangesStorage: Bool = false
    /// Heuristic flag indicating whether the elements are numeric.
    internal var isNumericStorage: Bool = true
    /// Cached total (Σ value * frequency) for numeric datasets.
    internal var totalStorage: FP? = nil
    /// Bookkeeping: last removal event (value and count).
    internal var lastItemRemoved: (item: SSElement, count: Int)? = nil
    /// Bookkeeping: last addition event (value and count).
    internal var lastItemAdded: (item: SSElement, count: Int)? = nil
    /// For each distinct element, stores the 1-based insertion indices at which it appeared.
    ///
    /// Used to reconstruct `.original` order efficiently.
    internal var sequence: Dictionary<SSElement, Array<Int>> = [:]
    /// Absolute frequencies per distinct element.
    internal var items: Dictionary<SSElement, Int> = [:]
    /// Cumulative frequency cache per distinct element (recomputed on demand).
    internal var cumulativeFrequenciesStorage: Dictionary<SSElement, Int> = [:]
    /// Cumulative relative frequency cache per distinct element (recomputed on demand).
    internal var cumulativeRelativeFrequenciesStorage: Dictionary<SSElement, FP> = [:]
    /// Indicates whether derived frequency caches require recomputation.
    internal var needsUpdateFrequencies: Bool = true
    /// Cache for ascending order of distinct items (not actively maintained in this file).
    internal var itemsAscending: Array<SSElement> = []
    /// Total number of observations (including duplicates).
    internal var count: Int = 0
    
    /// Optional root object to associate external metadata with this examine.
    public var rootObject: Any? {
        get { withLock { rootObjectStorage } }
        set { withLock { rootObjectStorage = newValue } }
    }
    /// Arbitrary tag or identifier. Defaults to a UUID when name is not provided.
    public var tag: String? {
        get { withLock { tagStorage } }
        set { withLock { tagStorage = newValue } }
    }
    /// Human-readable name for the examine instance. Defaults to `tag` if not explicitly set.
    public var name: String? {
        get { withLock { nameStorage } }
        set { withLock { nameStorage = newValue } }
    }
    /// The statistical level of measurement for the contained data.
    ///
    /// This does not constrain the API but can be used by consumers to select appropriate analyses.
    /// Defaults to `.nominal` for empty instances and when initialising from `String`.
    public var levelOfMeasurement: LevelOfMeasurement {
        get { withLock { levelOfMeasurementStorage } }
        set { withLock { levelOfMeasurementStorage = newValue } }
    }
    /// A descriptive string for debugging or UI display.
    public var descriptionString: String {
        get { withLock { descriptionStringStorage } }
        set { withLock { descriptionStringStorage = newValue } }
    }
    /// Significance level used by downstream statistical procedures (if any).
    ///
    /// Defaults to 0.05. Changing this value does not affect stored data, only subsequent analyses.
    public var alpha: FP {
        get { withLock { alphaStorage } }
        set { withLock { alphaStorage = newValue } }
    }
    /// True if there are no observations.
    public var isEmpty: Bool { withLock { count == 0 } }
    /// Flag indicating the data has been mutated since last known stable state.
    public var hasChanges: Bool {
        get { withLock { hasChangesStorage } }
        set { withLock { hasChangesStorage = newValue } }
    }
    /// Heuristic marker indicating whether the elements are numeric.
    ///
    /// This is set based on the first ingested element (or source kind) and is a best-effort hint.
    public var isNumeric: Bool {
        get { withLock { isNumericStorage } }
        set { withLock { isNumericStorage = newValue } }
    }
    
    /// A snapshot of cumulative absolute frequencies for all distinct items.
    ///
    /// The dictionary is recomputed lazily when invalidated by mutations.
    /// Access is thread-safe and will refresh internal caches as needed.
    /// Cumulative absolute frequencies F\_cum(x) over the observed support.
    ///
    /// Values are recomputed lazily under lock to satisfy `F_cum(x_j) = Σ_{i≤j} f(x_i)`
    /// in ascending key order, mirroring the empirical CDF numerator.
    public var cumulativeFrequencies: Dictionary<SSElement, Int> {
        get {
            withLock {
                let _ = updateCumulativeFrequenciesLocked()
                return cumulativeFrequenciesStorage
            }
        }
    }

    /// A snapshot of cumulative relative frequencies for all distinct items.
    ///
    /// The dictionary is recomputed lazily when invalidated by mutations.
    /// Access is thread-safe and will refresh internal caches as needed.
    public var cumulativeRelativeFrequencies: Dictionary<SSElement, FP>{
        get {
            withLock {
                let _ = updateCumulativeFrequenciesLocked()
                return cumulativeRelativeFrequenciesStorage
            }
        }
    }

   
    /// Creates an empty examine with default metadata.
    ///
    /// The instance starts with:
    /// - `alpha = 0.05`
    /// - `levelOfMeasurement = .nominal`
    /// - `isNumeric = true`
    /// - Empty sequence and frequency maps
    public override init() {
        super.init( )
        initialize()
    }
    
    /// Resets all internal state to the default empty configuration.
    ///
    /// This clears all observations, frequency maps, and caches. It preserves the default `alpha = 0.05`.
    fileprivate func initialize() {
        withLock {
            initializeLocked()
        }
    }
    
    /// Creates an examine by ingesting a source object (a `String` or an array of `SSElement`).
    ///
    /// - Parameters:
    ///   - object: Either a `String` (only if `SSElement == String`) or `[SSElement]`. For `String`,
    ///             characters are filtered by `characterSet` and appended as `String` elements.
    ///   - levelOfMeasurement: The measurement scale to associate with this data.
    ///   - name: Optional human-readable name. If `nil`, a UUID-based `tag` and `name` are created.
    ///   - characterSet: Optional filter when ingesting from `String`. Only unicode scalars contained in this set are appended.
    /// - Throws: `SSError(.missingData, ...)` if `object` is not a supported type for the current `SSElement`.
    /// - Note: When ingesting from `String`, `levelOfMeasurement` is set to `.nominal` regardless of the provided value.
    public init(using object: Any, levelOfMeasurement: LevelOfMeasurement, name: String?, characterSet: CharacterSet?) throws {
        let isStringElement = (SSElement.self is String.Type)
        guard (object is String && isStringElement) || (object is Array<SSElement>) else {
            let message = "Cannot initialize SSExamine with object of type: \(type(of: object))"
            SSLog.developerError(message)
            throw SSError(type: .missingData, file: #fileID, line: #line, function: #function, reason: message)
        }
        super.init()
        self.levelOfMeasurement = levelOfMeasurement
        self.createName(name: name)
        if let s = object as? String {
            self.levelOfMeasurement = .nominal
            _ = self.initialize(string: s, characterSet: characterSet)
        }
        else if let array = object as? Array<SSElement> {
            _ = self.initialize(using: array)
        }
    }

    /// Initialises the instance from an array of elements, updating cumulative caches.
    ///
    /// - Parameter a: Source elements.
    /// - Parameter name: Optional name of the instance
    /// - Parameter characterSet: Optional `CharacterSet` to use (only meaningful for strings).
    /// - Returns: `true` if at least one element was appended; otherwise `false`.
    /// - Thread safety: Acquires the lock.
    public init(usingArray a: Array<SSElement>, name: String?, characterSet: CharacterSet?) {
        super.init()
        createName(name: name)
        _ = self.initialize(using: a)
    }
    
    /// Creates a deep copy of the receiver, including metadata and the raw sequence.
    ///
    /// - Returns: A new `SSExamine` instance with the same elements (in the same order) and metadata.
    /// - Note: The copy preserves `hasChanges` from the source instance.
    public func copy(with zone: NSZone? = nil) -> Any {
        let snapshot = withLock {
            (
                isEmpty: count == 0,
                tag: tagStorage,
                description: descriptionStringStorage,
                name: nameStorage,
                alpha: alphaStorage,
                level: levelOfMeasurementStorage,
                hasChanges: hasChangesStorage,
                rawItems: itemsAsArrayLocked(sorted: .original)
            )
        }
        let result = SSExamine<SSElement, FP>()
        if !snapshot.isEmpty {
            result.tag = snapshot.tag
            result.descriptionString = snapshot.description
            result.name = snapshot.name
            result.alpha = snapshot.alpha
            result.levelOfMeasurement = snapshot.level
            // Bulk load to preserve order without marking as changed
            _ = result.initialize(using: snapshot.rawItems)
            // Preserve original change state (or set to false if you prefer clean copies)
            result.hasChanges = snapshot.hasChanges
        }
        return result
    }

    /// Compares two examine instances by their raw sequences.
    ///
    /// - Parameter object: Another object to compare against.
    /// - Returns: `true` if `object` is an `SSExamine<SSElement, FP>` and both instances have identical raw sequences,
    ///   or both are empty; otherwise `false`.
    public override func isEqual(_ object: Any?) -> Bool {
        guard let other = object as? SSExamine<SSElement, FP> else {
            return false
        }
        return self.itemsAsArray(sorted: .original) == other.itemsAsArray(sorted: .original)
    }
    
    /// Coding keys used for `Codable` conformance.
    internal enum codingKeys: String, CodingKey {
        case tag = "TAG"
        case name
        case descriptionString
        case alpha
        case levelOfMeasurement
        case isNumeric
        case data
    }
    
    /// Encodes metadata and the raw data sequence.
    ///
    /// The encoded payload contains:
    /// - `name`, `tag`, `descriptionString`, `alpha`, `levelOfMeasurement`, `isNumeric`
    /// - `data`: the full raw sequence `[SSElement]` (omitted if empty)
    public func encode(to encoder: any Encoder) throws {
        let snapshot = withLock {
            (
                name: nameStorage,
                tag: tagStorage,
                description: descriptionStringStorage,
                alpha: alphaStorage,
                level: levelOfMeasurementStorage,
                isNumeric: isNumericStorage,
                raw: itemsAsArrayLocked(sorted: .original)
            )
        }
        var container = encoder.container(keyedBy: codingKeys.self)
        try container.encodeIfPresent(snapshot.name, forKey: .name)
        try container.encodeIfPresent(snapshot.tag, forKey: .tag)
        try container.encodeIfPresent(snapshot.description, forKey: .descriptionString)
        try container.encode(snapshot.alpha, forKey: .alpha)
        try container.encode(snapshot.level, forKey: .levelOfMeasurement)
        try container.encode(snapshot.isNumeric, forKey: .isNumeric)
        if !snapshot.raw.isEmpty {
            try container.encode(snapshot.raw, forKey: .data)
        }
    }
    
    /// Decodes an instance from an encoded payload.
    ///
    /// Decoding strategy:
    /// - Reads metadata fields and applies them directly to storage.
    /// - If a `data` array is present, the instance is populated using the raw sequence via `initialize(using:)`.
    ///
    /// Defaults if fields are absent:
    /// - `descriptionString`: "SSExamine"
    /// - `alpha`: `.zero` (but will be reset to `0.05` if `initialize(using:)` is called due to its internal reset)
    /// - `levelOfMeasurement`: `.interval`
    /// - `isNumeric`: `true`
    ///
    /// Important:
    /// - When `data` is present, `initialize(using:)` eventually calls `initializeLocked()`, which resets `alpha` to `0.05`
    ///   and `levelOfMeasurement` to `.nominal`. This matches the test expectations and ensures a canonical post-load state.
    required public init(from decoder: Decoder) throws {
        super.init()
        let container = try decoder.container(keyedBy: codingKeys.self)
        let decodedTag = try container.decodeIfPresent(String.self, forKey: .tag)
        let decodedName = try container.decodeIfPresent(String.self, forKey: .name)
        let decodedDescription = try container.decodeIfPresent(String.self, forKey: .descriptionString) ?? "SSExamine"
        let decodedAlpha = try container.decodeIfPresent(FP.self, forKey: .alpha) ?? FP.zero
        let decodedLevel = try container.decodeIfPresent(LevelOfMeasurement.self, forKey: .levelOfMeasurement) ?? .interval
        let decodedIsNumeric = try container.decodeIfPresent(Bool.self, forKey: .isNumeric) ?? true
        let decodedData = try container.decodeIfPresent(Array<SSElement>.self, forKey: .data)
        withLock {
            tagStorage = decodedTag
            nameStorage = decodedName
            descriptionStringStorage = decodedDescription
            alphaStorage = decodedAlpha
            levelOfMeasurementStorage = decodedLevel
            isNumericStorage = decodedIsNumeric
        }
        if let data = decodedData {
            _ = initialize(using: data)
        }
    }
    
    /// Hash of the raw sequence.
    ///
    /// - Returns 0 for empty datasets.
    /// - For non-empty datasets, hashes the raw array of items to produce a stable value across equal instances.
    public  override var hash: Int {
        guard !self.isEmpty else { return 0 }
        return withLock {
            let a = itemsAsArray(sorted: .original)
            var hasher = Hasher()
            hasher.combine(a)
            return hasher.finalize()
        }
    }

}

/// Concurrency marker indicating instances can be shared across threads/tasks.
///
/// Internal synchronization is provided via a recursive lock.
/// This is marked `@unchecked` because `NSObject` is not `Sendable`,
/// but the class enforces thread-safety for its mutable state.

extension SSExamine: @unchecked Sendable {}


/// Container-style APIs for `SSExamine` that treat the sample as a counting measure.
///
/// Frequencies are derived as `f(x)` and `f(x)/n` under a lock, preserving probability
/// mass even when concurrent writers append or remove data. These utilities keep the
/// multiset view (counts) and the raw insertion order in sync so statistical estimators
/// can treat the same dataset either as a bag of observations or as an ordered sequence.
extension SSExamine: ExamineContainer {
    /// The total number of observations (including duplicates).
    public var sampleSize: Int {
        return withLock { count }
    }
    
    /// The number of distinct elements present.
    public var length: Int {
        return withLock { items.count }
    }
    
    /// Indicates whether a specific element is present at least once.
    /// - Parameter item: The element to check.
    /// - Returns: `true` if present, otherwise `false`.
    public func contains(_ item: SSElement) -> Bool {
        return withLock { items[item] != nil }
    }
    
    /// Returns the absolute frequency of an element.
    /// - Parameter item: The element whose frequency is requested.
    /// - Returns: The count of occurrences, or 0 if absent.
    public func frequency(_ item: SSElement) -> Int {
        return withLock { items[item] ?? 0 }
    }
    
    /// Returns the relative frequency of an element.
    /// - Parameter item: The element whose relative frequency is requested.
    /// - Returns: The proportion of this element in the sample, or 0 if the sample is empty or the element is absent.
    public func relativeFrequency(_ item: SSElement) -> FP {
        return withLock {
            guard let elementCount = items[item], count > 0 else { return 0 }
            return FP(elementCount) / FP(count)
        }
    }
    
    /// Appends a single observation to the sample.
    ///
    /// Updates the absolute frequency and records the insertion index to enable `.original` reconstruction.
    /// - Parameter item: The element to append.
    public func append(_ item: SSElement) {
        withLock {
            appendLocked(item)
            let _ = updateCumulativeFrequenciesLocked()
        }
    }
    
    /// Appends the same observation multiple times.
    /// - Parameters:
    ///   - count: The number of times to append. Must be > 0.
    ///   - item: The element to append.
    /// - Returns: `true` if any items were appended; otherwise `false`.
    public func append(repeating count: Int, item: SSElement) -> Bool {
        guard count > 0 else { return false }
        return withLock {
            for _ in 0..<count {
                appendLocked(item)
            }
            needsUpdateFrequencies = true
            return updateCumulativeFrequenciesLocked()
        }
    }
    
    /// Appends the contents of an array in order.
    ///
    /// - Heuristically sets `isNumeric` based on the first element if the instance is currently empty.
    /// - Parameter array: The elements to append.
    /// - Returns: `true` if any items were appended; otherwise `false`.
    public func append(contentOf array: Array<SSElement>) -> Bool {
        guard let firstElement = array.first else { return false }
        return withLock {
            if count == 0 {
                isNumericStorage = Helpers.isNumeric(firstElement)
            }
            for element in array {
                appendLocked(element)
            }
            needsUpdateFrequencies = true
            return updateCumulativeFrequenciesLocked()
        }
    }
    
    /// Appends characters of a string as `String` elements, optionally filtered by a `CharacterSet`.
    ///
    /// - Important: Requires `SSElement == String`.
    /// - Parameters:
    ///   - text: The source text.
    ///   - characterSet: Optional filter. Only scalars included in the set are appended.
    /// - Throws: `SSError(.missingData, ...)` if `SSElement` is not `String`.
    /// - Returns: `true` if any items were appended; otherwise `false`.
    public func append(text: String, characterSet: CharacterSet?) throws -> Bool {
        guard (SSElement.self is String.Type) else {
            let message = "append(text:characterSet:) requires SSElement == String"
            SSLog.statisticsError(message)
            throw SSError(type: .missingData, file: #fileID, line: #line, function: #function, reason: message)
        }
        guard !text.isEmpty else { return false }
        let additions: [SSElement]
        if let cs = characterSet {
            additions = text.unicodeScalars.compactMap { scalar in
                cs.contains(scalar) ? (String(scalar) as? SSElement) : nil
            }
        } else {
            additions = text.map { String($0) as! SSElement }
        }
        guard !additions.isEmpty else { return false }
        return withLock {
            for element in additions {
                appendLocked(element)
            }
            hasChangesStorage = true
            needsUpdateFrequencies = true
            return updateCumulativeFrequenciesLocked()
        }
    }
    
    /// Removes occurrences of an element from the sample.
    ///
    /// - Parameters:
    ///   - item: The element to remove.
    ///   - allOccurences: If `true`, removes all occurrences; if `false`, removes only the first occurrence in raw order.
    /// - Note: This operation rebuilds internal state from the filtered raw array and is O(n).
    public func remove(_ item: SSElement, allOccurences: Bool) {
        let current = itemsAsArray(sorted: .original)
        let filtered: [SSElement]
        if allOccurences {
            filtered = current.filter { $0 != item }
        } else {
            var removed = false
            filtered = current.filter { elem in
                if !removed && elem == item {
                    removed = true
                    return false
                }
                return true
            }
        }
        withLock {
            _ = initializeLocked(using: filtered)
            hasChangesStorage = true
            needsUpdateFrequencies = true
            let _ = updateCumulativeFrequenciesLocked()
        }
    }
    
    /// Removes all observations and resets caches and flags.
    public func removeAll() {
        withLock {
            initializeLocked()
            hasChangesStorage = true
            let _ = updateCumulativeFrequenciesLocked()
        }
    }

}

/// Convenience renderings and projections of the sample's items.
///
/// These methods expose the empirical distribution `F_n` in different coordinate
/// systems: raw insertion order, order statistics, or frequency-ranked categories.
/// All accessors in this extension are thread-safe and acquire the instance's lock
/// internally via `withLock` to provide a consistent snapshot of the underlying data.
/// Unless otherwise noted, operations that materialize arrays are O(n) with respect
/// to the number of observations, plus any cost for sorting if a sorted order is requested.
extension SSExamine {
    /// Renders the items as a delimited string or a multi-line list.
    ///
    /// - Parameters:
    ///   - del: The delimiter to use when `asRow` is `true`. If `nil`, an empty string is used.
    ///   - asRow: If `true`, returns a single-line row with elements joined by `del`;
    ///            if `false`, returns a multi-line string with one element per line,
    ///            optionally prefixed by `name` (if present).
    ///   - encl: If non-`nil` and non-empty, each element string is enclosed by this string.
    ///           Occurrences of `encl` inside elements are escaped by doubling the enclosure string.
    /// - Returns: A formatted string, or `nil` if there are no items to render.
    ///
    /// Example:
    /// ```swift
    /// let ex = SSExamine<String, Double>()
    /// _ = try ex.append(text: "ab", characterSet: .letters)
    /// let csv = ex.itemsAsString(delimiter: ",", asRow: true, encloseElementsBy: "\"")
    /// // csv == "a,b"
    /// ```
    public func itemsAsString(delimiter del: String?, asRow: Bool = true, encloseElementsBy encl: String? = nil) -> String? {
        let e = self.itemsAsArray(sorted: .original)
        let parts: [String] = e.map { item in
            let itemString = "\(item)"
            if let e = encl, !e.isEmpty {
                let escaped = itemString.replacingOccurrences(of: e, with: e + e)
                return e + escaped + e
            }
            else {
                return itemString
            }
        }
        var res: String
        if asRow {
            res = parts.joined(separator: del ?? "")
        }
        else {
            var lines: [String] = []
            if let n = self.name {
                lines.append(n)
            }
            lines.append(contentsOf: parts)
            res = lines.joined(separator: "\n")
        }
        if res.isEmpty {
            SSLog.statisticsError("itemsAsString produced an empty rendering for dataset '\(self.name ?? "unnamed")'")
            return nil
        }
        return res
    }
    
    /// Materialises the sample into an array according to the requested order.
    ///
    /// - Parameter sorted: The desired sort order.
    ///   - `.ascending` / `.descending`: Returns the multiset sorted by element order.
    ///   - `.none`: Returns the current multiset as-is without ordering guarantees (hash-map iteration order).
    ///   - `.original`: Reconstructs the original insertion order using internal sequence indices.
    /// - Returns: An array of all observations in the requested order.
    /// - Complexity: O(n log n) for `.ascending`/`.descending`; O(n) for `.original` and `.none`.
    public func itemsAsArray(sorted: SortOrderDataArray) -> Array<SSElement> {
        return withLock {
            itemsAsArrayLocked(sorted: sorted)
        }
    }
    
    /// Returns the set of unique items present in the sample in the requested order.
    ///
    /// - Parameter sorted: The desired order for unique items (e.g., by natural element order,
    ///   by frequency, or by insertion/original order depending on `SortOrderUniqeItems`).
    /// - Returns: An array containing each distinct `SSElement` exactly once, ordered per `sorted`.
    /// - Complexity: Typically O(n) to collect unique elements plus any sorting cost dictated by `sorted`.
    public func uniqueItems(sorted: SortOrderUniqeItems) -> Array<SSElement> {
        guard !isEmpty else { return [] }
        return withLock {
            uniqueItemsLocked(sorted: sorted)
        }
    }
    
    /// Builds a frequency table of the sample's items.
    ///
    /// - Parameter sorted: The desired order for the resulting frequency rows (e.g., by key,
    ///   by count ascending/descending, or by insertion/original order depending on `SortOrderFrequencyTable`).
    /// - Returns: An array of frequency table rows, each describing an item and its frequency-related metrics.
    /// - Complexity: O(n) to accumulate counts plus any sorting cost dictated by `sorted`.
    public func frequencyTable(sorted: SortOrderFrequencyTable) -> Array<FrequencyTableItem<SSElement, FP>> {
        return withLock {
            frequencyTableLocked(sorted: sorted)
        }
    }
    
    /// Builds a cumulative frequency table of the sample's items.
    ///
    /// - Parameter format: The desired cumulative table format, which controls the aggregation domain
    ///   (e.g., cumulative by sorted key, cumulative probability, running totals, etc.).
    /// - Returns: An array of cumulative frequency entries in the order implied by `format`.
    /// - Complexity: O(n) to accumulate plus any sorting cost implied by `format`.
    public func cumulativeFrequencyTable(format: CumulativeFrequencyTableFormat) -> Array<CumulativeFrequencyTableItem<SSElement, FP>> {
        return withLock {
            cumulativeFrequencyTableLocked(format: format)
        }
    }
    
    /// Attempts to project all items to a numeric array.
    ///
    /// If `SSElement` is numeric or can be mapped to `FP` (the floating-point type used by this instance),
    /// returns the corresponding numeric values. If any item cannot be represented as `FP`, this may return `nil`
    /// (depending on the implementation of `itemsAsNumericArrayLocked()`).
    ///
    /// - Returns: An array of numeric values corresponding to the items, or `nil` if the projection is unavailable.
    /// - Note: The returned values reflect a snapshot taken under lock at the time of the call.
    public var itemsAsNumericArray: Array<FP>? {
        return withLock {
            return itemsAsNumericArrayLocked()
        }
    }
    
    /// Accesses the item at the given index in original insertion order.
    ///
    /// - Parameter index: The zero-based index in the original order of insertion.
    /// - Returns: The element at `index` if it exists; otherwise `nil`.
    /// - Note: This subscript materializes the items in original order and then indexes into that array.
    ///         It validates bounds and emptiness before attempting access.
    subscript (index: Int) -> SSElement? {
        guard isValidIndex(index: index), !self.isEmpty else {
            SSLog.statisticsError("Index \(index) is invalid or dataset is empty for SSExamine subscript access")
            return nil
        }
        let a = self.itemsAsArray(sorted: .original)
        return a[index]
    }
    
    /// Evaluates the empirical cumulative distribution function (ECDF) at the given item.
    ///
    /// The ECDF at `x` is defined as the proportion of observations less than or equal to `x`
    /// (depending on the exact definition used internally). The result is expressed in the floating-point
    /// type `FP` associated with this instance.
    ///
    /// - Parameter item: The query value.
    /// - Returns: The ECDF value in `[0, 1]` if it can be computed; otherwise `nil` (e.g., for an empty sample).
    public func empiricalCDF(item: SSElement) -> FP? {
        return withLock { () -> FP? in
            return self.empiricalCDFLocked(item)
        }
    }
}

extension SSExamine {
    /// Convenience goodness-of-fit check for normality using the one-sample Kolmogorov-Smirnov test.
    ///
    /// - Returns: `true` when the bootstrap p-value is at least `alpha` (fail to reject normality),
    ///   `false` when normality is rejected, or `nil` when the dataset is empty or non-numeric.
    /// - Discussion: Uses the container's `alpha` and the `.normal` target from ``Inferential/HypothesisTesting/GoodnessOfFit/KSTest``.
    public var isGaussian: Bool? {
        // Grab numeric snapshot under lock without forcing the closure to return Bool.
        let data: [FP]? = withLock { self.itemsAsNumericArrayLocked() }
        guard let values = data, !values.isEmpty else {
            SSLog.statisticsError("Kolmogorov-Smirnov test requires numeric data")
            return nil
        }
        // Run KS test against Normal (N(0,1) or estimated μ, σ depending on your target choice).
        let result: GoodnessOfFit<FP>.KSTestResult = GoodnessOfFit<FP>.ksTestOneSample(data: values, testTarget: .normal)
        // Decide Gaussian-ness by comparing p-value to alpha (use container’s alpha).
        // Convention: return true if we fail to reject normality (p >= alpha).
        let alphaValue = FP(self.alpha)
        return result.pBootstrap >= alphaValue
    }
    

    /// Runs a one-sample Kolmogorov-Smirnov test against a specified target distribution.
    ///
    /// - Parameter target: The distribution to test against (normal, logistic, studentT, or custom CDF via the target enum).
    /// - Returns: The ``Inferential/HypothesisTesting/GoodnessOfFit/KSTestResult`` when the dataset is non-empty and numeric, or `nil` otherwise.
    /// - Note: Uses the bootstrap-based KS implementation and does not mutate the underlying container.
    public func testForDistribution(target: GOFTestTarget) -> GoodnessOfFit<FP>.KSTestResult? {
        let data: [FP]? = withLock { self.itemsAsNumericArrayLocked() }
        guard let values = data, !values.isEmpty else {
            SSLog.statisticsError("Kolmogorov-Smirnov test requires numeric data")
            return nil
        }
        let result = GoodnessOfFit<FP>.ksTestOneSample(data: values, testTarget: target)
        return result
    }
    
}


/// Lock-aware helpers shared across SSExamine location, dispersion, and entropy routines.
///
/// These utilities manage reference counting of cached aggregates (totals, sums of powers)
/// under a recursive lock so moment-based estimators respect algebraic identities such as
/// `Σ x = n·E[X]` even in concurrent contexts. Their design keeps numerical stability in
/// focus by preferring compensated accumulation and logarithmic summations when the data
/// spans many magnitudes.
extension SSExamine {
    /// Locked implementation of `inverseTotal`.
    ///
    /// See `inverseTotal` for behavior and edge cases.
    ///
    /// - Returns: `nil` for non-numeric datasets; `0` for empty numeric datasets; `NaN` if undefined.
    /// - Thread safety: The caller must hold the lock.
    @inline(__always)
    internal func inverseTotalLocked() -> FP? {
        guard self.isNumeric else {
            SSLog.statisticsError("inverseTotalLocked requires numeric data")
            return nil
        }
        guard !self.isEmpty else {
            // Important branch: empty dataset -> 0 for numeric, NaN for non-numeric (defensive).
            if self.isNumeric {
                return FP.zero
            }
            else {
                SSLog.statisticsError("inverseTotalLocked received an empty non-numeric dataset")
                return nil
            }
        }
        do {
            // Build terms first, fail fast on NaN, then sum by increasing |term|
            var terms = [FP]()
            terms.reserveCapacity(items.count)
            var fpValue = FP.zero
            var fpFreq = FP.zero
            for (item, freq) in self.items {
                fpValue = try RealConverter.to(item, locale: nil)
                fpFreq = FP(freq)
                // Important: division can yield ±∞ (when value == 0); allowed and handled by `sum`.
                let term = fpFreq / fpValue
                if term.isNaN {
                    // Important branch: NaN encountered (e.g., 0/0), short-circuit to NaN.
                    return FP.nan
                }
                terms.append(term)
            }
            return Helpers.sum(&terms)
        }
        catch _ {
            // Important branch: conversion failed (non-numeric); propagate as nil to signal not applicable.
            SSLog.statisticsError("inverseTotalLocked failed to convert one or more values to a numeric representation")
            return nil
        }
    }

    /// Locked implementation of `powerSum(power:)`.
    ///
    /// See `powerSum(power:)` for behavior and edge cases.
    ///
    /// Complexity:
    /// - O(k log k) for k distinct items (due to sorting inside `sum`).
    /// - Conversion of items to `BinaryFloatingPoint` is O(k).
    ///
    /// - Thread safety: The caller must hold the lock.
    @inline(__always)
    internal func powerSumLocked(power: FP) -> FP? {
        guard self.isNumeric else {
            SSLog.statisticsError("powerSumLocked requires numeric data")
            return nil
        }
        guard !self.isEmpty else {
            // Important branch: empty dataset -> 0 for numeric, NaN for non-numeric (defensive).
            if self.isNumeric {
                return FP.zero
            }
            else {
                SSLog.statisticsError("powerSumLocked received an empty non-numeric dataset")
                return nil
            }
        }
        // Fast paths
        if power == 1 {
            // Important branch: reuse cached/locked total.
            return self.total
        }
        else if power == FP.zero {
            // Important branch: Σ x^0 * f = Σ 1 * f = count (sample size).
            return FP(self.count)
        }
        do {
            // Build terms first, fail fast on NaN, then sum by increasing |term|
            var terms = [FP]()
            terms.reserveCapacity(items.count)
            var fpValue = FP.zero
            var fpFreq = FP.zero
            let p = FP(power)
            var term: FP = 0
            for (item, freq) in self.items {
                fpValue = try RealConverter.to(item, locale: nil)
                fpFreq = FP(freq)
                if p == 2 {
                    // Faster and more precise than pow(fpValue, 2)
                    term = fpValue * fpValue * fpFreq
                }
                else {
                    // Important branch: pow can yield NaN for negative base with non-integer exponent.
                    term = FP.pow(fpValue, p) * fpFreq
                }
                if term.isNaN {
                    // Important branch: if any term is NaN, the entire sum is NaN.
                    return FP.nan
                }
                terms.append(term)
            }
            return Helpers.sum(&terms)
        }
        catch _ {
            // Important branch: conversion failed; in this path we signal NaN to distinguish from non-numeric dataset.
            SSLog.statisticsError("powerSumLocked failed to convert one or more values to a numeric representation")
            return FP.nan
        }
    }
    
    /// Updates `totalStorage` with a numerically stable sum of Σ value * frequency.
    ///
    /// - For empty numeric datasets, stores `0`.
    /// - For empty non-numeric datasets, stores `NaN`.
    ///
    /// Implementation detail:
    /// - Builds terms eagerly to fail fast on `NaN`, then sums via `sum(_:)`
    ///   (sorted by magnitude with compensation).
    ///
    /// Thread-safety:
    /// - Must be called under `withLock`.
    @inline(__always)
    internal func updateTotalLocked() {
        guard self.isNumeric else {
            // Important branch: non-numeric dataset -> clear cache (nil).
            self.totalStorage = nil
            return
        }
        guard !self.isEmpty else {
            // Important branch: empty numeric dataset -> total is 0.
            if self.isNumeric {
                self.totalStorage = FP.zero
                return
            }
            else {
                // Defensive: non-numeric empty -> nil.
                self.totalStorage = nil
                return
            }
        }
        do {
            // Build terms first, fail fast on NaN, then sum by increasing |term|
            var terms = [FP]()
            terms.reserveCapacity(items.count)
            var fpValue = FP.zero
            var fpFreq = FP.zero
            for (item, freq) in self.items {
                fpValue = try RealConverter.to(item, locale: nil)
                fpFreq = FP(freq)
                let term = fpValue * fpFreq
                if term.isNaN {
                    // Important branch: NaN term (e.g., undefined arithmetic) invalidates the total.
                    self.totalStorage = FP.nan
                    return
                }
                terms.append(term)
            }
            self.totalStorage = Helpers.sum(&terms)
        }
        catch _ {
            // Important branch: conversion failed; store NaN to indicate invalid numeric state.
            self.totalStorage = FP.nan
        }
    }
    
    /// Locked implementation for `logSum`.
    ///
    /// - Returns: 0 for empty numeric datasets, -∞ if any value is zero, `nil` if any value < 0 or NaN.
    /// - Uses stable summation on logarithms of values.
    /// - Thread safety: The caller must hold the lock.
    @inline(__always)
    internal func logSumLocked() -> FP? {
        guard self.isNumeric else {
            SSLog.statisticsError("logSumLocked requires numeric data")
            return nil
        }
        if count == 0 { return 0 } // Important branch: empty dataset -> neutral element for addition is 0.
        guard let values: [FP] = self.itemsAsNumericArrayLocked() else {
            SSLog.statisticsError("logSumLocked could not project items to numeric values")
            return nil
        }
        // Build logs directly to avoid an intermediate array.
        var logs = [FP]()
        logs.reserveCapacity(values.count)
        for fpValue in values {
            if fpValue.isNaN || fpValue < 0 {
                // Important branch: log undefined for negative or NaN -> signal not applicable.
                SSLog.statisticsError("logSumLocked encountered a negative or NaN value")
                return nil
            }
            if fpValue.isZero {
                // Important branch: log(0) = -∞ -> return immediately.
                return -FP.infinity
            }
            logs.append(FP.log(fpValue))
        }
        return Helpers.sum(&logs)
    }
    
    /// Sum of logs of magnitudes: Σ log(|x_i|).
    ///
    /// - Returns nil if any value is NaN.
    /// - Returns -∞ if any value is exactly zero.
    ///
    /// Notes:
    /// - This is a helper for logarithmic product computation where sign is handled separately.
    /// - Thread safety: Pure function over the provided values.
    @inline(__always)
    internal func logSumOfMagnitudes(_ values: [FP]) -> FP? {
        if values.isEmpty { return 0 } // Important branch: empty set -> neutral element for addition is 0.
        var logs = [FP]()
        logs.reserveCapacity(values.count)
        for v in values {
            if v.isNaN {
                SSLog.statisticsError("logSumOfMagnitudes encountered NaN input")
                return nil
            }
            let a = v.magnitude
            if a.isZero { return -FP.infinity } // Important branch: log(0) = -inf.
            // log of magnitude; allow +∞ if a == +∞
            let lv = FP.log(a)
            if lv.isNaN {
                SSLog.statisticsError("logSumOfMagnitudes produced NaN while computing log(\(a))")
                return nil
            }
            logs.append(lv)
        }
        return Helpers.sum(&logs)
    }
}

extension SSExamine {
    /// Executes the given closure while holding the instance's recursive lock.
    ///
    /// Use this helper to serialize access to mutable internal state.
    ///
    /// - Parameter body: A closure to execute while the lock is held.
    /// - Returns: The value returned by `body`.
    /// - Throws: Rethrows any error thrown by `body`.
    /// - Thread safety: Acquires and releases the instance lock.
    @inline(__always)
    internal func withLock<T>(_ body: () throws -> T) rethrows -> T {
        lock.lock()
        defer { lock.unlock() }
        return try body()
    }

    /// Convenience helper to log a statistics error message and return `nil`.
    ///
    /// - Parameter message: Human-readable reason why the operation could not produce a value.
    /// - Returns: Always `nil`, allowing call-sites to stay concise.
    @inline(__always)
    internal func logNil<T>(_ message: String) -> T? {
        SSLog.statisticsError(message)
        return nil
    }
}
