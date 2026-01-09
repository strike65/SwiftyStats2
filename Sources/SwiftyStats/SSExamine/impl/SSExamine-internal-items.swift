//
//  Created by VT on 22.11.25.
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

/// Internal item management and sequence tracking for SSExamine containers.

import SwiftyStatsPrelude
import SwiftyBoost

// MARK: Items

/// Internal, lock-required helpers for mutating and materializing sample data.
///
/// These routines keep the empirical distribution `F_n` coherent by updating both
/// the counting measure (frequency dictionary) and the raw insertion order indices
/// in the same critical section. All methods in this extension that end with
/// `Locked` must be called only while the instance's lock is held. Public-facing
/// wrappers typically acquire the lock via `withLock { ... }` and then delegate
/// to these helpers.
extension SSExamine {
    /// Appends a single observation to the dataset (locked variant).
    ///
    /// - Important: The caller must already hold the lock.
    /// - Effects:
    ///   - Increments `count`.
    ///   - Updates absolute frequencies in `items`.
    ///   - Appends the new 1-based index to `sequence[item]`.
    ///   - Marks cumulative frequencies as needing update.
    ///   - Sets `hasChangesStorage = true`.
    /// - Parameter item: The observation to append.
    @inline(__always)
    internal func appendLocked(_ item: SSElement) {
        count += 1
        items[item, default: 0] += 1
        sequence[item, default: []].append(count)
        needsUpdateFrequencies = true
        hasChangesStorage = true
    }
    
    /// Resets all internal storage to a clean initial state (locked variant).
    ///
    /// - Important: The caller must already hold the lock.
    /// - Effects:
    ///   - Clears `sequence`, `items`, and cumulative frequency caches.
    ///   - Resets `needsUpdateFrequencies`, `itemsAscending`, and `count`.
    ///   - Restores default metadata:
    ///     - `descriptionStringStorage = "SSExamine instance - standard"`
    ///     - `hasChangesStorage = false`
    ///     - `isNumericStorage = true`
    ///     - `alphaStorage = 0.05`
    ///     - `levelOfMeasurementStorage = .nominal`
    @inline(__always)
    internal func initializeLocked() {
        sequence.removeAll()
        items.removeAll()
        cumulativeFrequenciesStorage.removeAll()
        cumulativeRelativeFrequenciesStorage.removeAll()
        needsUpdateFrequencies = true
        itemsAscending.removeAll()
        count = 0
        descriptionStringStorage = "SSExamine instance - standard"
        hasChangesStorage = false
        isNumericStorage = true
        alphaStorage = FP(0.05)
        levelOfMeasurementStorage = .nominal
    }
    
    /// Bulk-initializes the instance from an array of elements (locked variant).
    ///
    /// - Important: The caller must already hold the lock.
    /// - Effects:
    ///   - Resets the instance via `initializeLocked()`.
    ///   - Heuristically sets `isNumericStorage` based on the first element.
    ///   - Appends all elements via `appendLocked(_:)`.
    ///   - Marks cumulative frequencies as needing update.
    ///   - Does not mark the instance as changed (bulk load).
    /// - Parameter array: Source elements. No-op if empty.
    /// - Returns: `true` if at least one element was appended; otherwise `false`.
    /// - Complexity: O(n) to append n elements.
    @inline(__always)
    @discardableResult
    internal func initializeLocked(using array: Array<SSElement>) -> Bool {
        initializeLocked()
        guard !array.isEmpty else {
            return false
        }
        isNumericStorage = Helpers.isNumeric(array.first as Any)
        for element in array {
            appendLocked(element)
        }
        needsUpdateFrequencies = true
        // Bulk initialization should not mark the instance as changed
        hasChangesStorage = false
        return true
    }
    
    /// Bulk-initializes the instance from a `String`, optionally filtering by `CharacterSet` (locked variant).
    ///
    /// - Important:
    ///   - Requires `SSElement == String`.
    ///   - The caller must already hold the lock.
    ///   - Each character (or scalar, when `characterSet` is provided) is converted to `String`
    ///     and appended as a separate observation.
    /// - Parameters:
    ///   - string: Source text. No-op if empty.
    ///   - characterSet: If provided, only scalars contained in the set are appended.
    /// - Returns: `true` if at least one element was appended; otherwise `false`.
    /// - Complexity: O(n) over number of characters or filtered scalars.
    @inline(__always)
    @discardableResult
    internal func initializeLocked(string: String, characterSet: CharacterSet?) -> Bool {
        initializeLocked()
        guard !string.isEmpty else { return false }
        guard (SSElement.self is String.Type) else { return false }
        isNumericStorage = false
        if let cs = characterSet {
            for scalar in string.unicodeScalars where cs.contains(scalar) {
                appendLocked(String(scalar) as! SSElement)
            }
        } else {
            for character in string {
                appendLocked(String(character) as! SSElement)
            }
        }
        needsUpdateFrequencies = true
        // Bulk initialization should not mark the instance as changed
        hasChangesStorage = false
        return true
    }
    
    /// Initialises the instance from an array of elements, updating cumulative caches.
    ///
    /// - Parameter array: Source elements.
    /// - Returns: `true` if at least one element was appended; otherwise `false`.
    /// - Thread safety: Acquires the lock.
    internal func initialize(using array: Array<SSElement>) -> Bool {
        return withLock {
            initializeLocked(using: array) && updateCumulativeFrequenciesLocked()
        }
    }

    /// Initialises the instance from a `String` of characters, optionally filtered by a `CharacterSet`.
    ///
    /// - Important: This requires `SSElement == String`. Each character (or scalar) is converted to `String`
    ///   and appended as a separate observation.
    /// - Parameters:
    ///   - string: The source text. No operation is performed if empty.
    ///   - characterSet: Optional filter. If provided, only scalars contained in the set are appended.
    /// - Returns: `true` if at least one element was appended; otherwise `false`.
    /// - Thread safety: Acquires the lock.
    internal func initialize(string: String, characterSet: CharacterSet?) -> Bool {
        return withLock {
            initializeLocked(string: string, characterSet: characterSet) && updateCumulativeFrequenciesLocked()
        }
    }

    /// Creates or assigns a name and tag.
    ///
    /// - If `name` is provided, it is set and `tag` is left unchanged.
    /// - If `name` is `nil`, a new UUID is generated and assigned to both `tag` and `name`.
    /// - Parameter name: Optional human-readable name to assign.
    /// - Thread safety: Acquires the lock.
    internal func createName(name: String?) {
        withLock {
            if let provided = name {
                nameStorage = provided
            } else {
                let identifier = UUID().uuidString
                tagStorage = identifier
                nameStorage = identifier
            }
        }
    }
    
    internal func checkNumericLocked() {
        guard self.count > 0 else {
            self.isNumericStorage = false
            return
        }
        if let a = self.itemsAsNumericArrayLocked() {
            self.isNumericStorage = a.allSatisfy({ Helpers.isNumeric($0) })
        }
        else {
            self.isNumericStorage = false
        }
    }
    
    /// Returns all observations as an array, with optional sorting (locked variant).
    ///
    /// - Important: The caller must already hold the lock.
    /// - Parameter sorted: Sorting mode.
    /// - Returns:
    ///   - `.ascending` / `.descending`: Sorted by natural order of `SSElement`.
    ///   - `.none`: Unordered expansion based on `items` counts.
    ///   - `.original`: Original insertion order reconstructed from `sequence`.
    /// - Complexity: O(n log n) when sorting; O(n) otherwise.
    @inline(__always)
    internal func itemsAsArrayLocked(sorted: SortOrderDataArray) -> Array<SSElement> {
        guard count > 0 else { return [] }
        var tmp: Array<SSElement> = []
        tmp.reserveCapacity(count)
        for (item, freq) in items {
            for _ in 0..<freq {
                tmp.append(item)
            }
        }
        switch sorted {
        case .ascending:
            return tmp.sorted(by: { $0 < $1 })
        case .descending:
            return tmp.sorted(by: { $0 > $1 })
        case .none:
            return tmp
        case .original:
            var raw = Array(repeating: tmp.first!, count: count)
            for (item, seq) in sequence {
                for index in seq {
                    raw[index - 1] = item
                }
            }
            return raw
        }
    }
    
    /// Returns the unique set of observed items, with optional sorting (locked variant).
    ///
    /// - Important: The caller must already hold the lock.
    /// - Parameter sorted: Sorting mode for unique items.
    /// - Returns: The array of unique `SSElement` values.
    /// - Complexity: O(k log k) for k unique items when sorting; O(k) otherwise.
    @inline(__always)
    internal func uniqueItemsLocked(sorted: SortOrderUniqeItems) -> Array<SSElement> {
        guard !self.isEmpty else { return [] }
        var tmp: Array<SSElement> = []
        tmp.reserveCapacity(count)
        for (item, _) in items {
            tmp.append(item)
        }
        switch sorted {
        case .ascending:
            return tmp.sorted(by: { $0 < $1 } )
        case .descending:
            return tmp.sorted(by: { $0 > $1 } )
        case .none:
            return tmp
        }
    }
    
    /// Checks whether a zero-based index is valid for the raw sequence.
    ///
    /// - Parameter index: The index to validate.
    /// - Returns: `true` if `0 <= index < count`; otherwise `false`.
    @inline(__always)
    internal func isValidIndex(index: Int) -> Bool {
        return index.liesInRightOpenRange(from: 0, to: count)
    }
    
    /// Attempts to project the observations into numeric values of type `FP` (locked variant).
    ///
    /// - Returns: An array of `FP` values if all elements can be converted; otherwise `nil`.
    /// - Throws: This method internally catches any conversion errors and returns `nil` on failure.
    /// - Discussion:
    ///   - Materializes the items with `.none` ordering and maps each element via `RealConverter.to(_:locale:)`.
    ///   - The conversion semantics depend on `RealConverter` and the concrete `SSElement` type.
    @inline(__always)
    internal func itemsAsNumericArrayLocked() -> Array<FP>? {
        do {
            let array = self.itemsAsArray(sorted: .none)
            return try array.map { try RealConverter.to($0, locale: nil) as FP }
        }
        catch _ {
            return nil
        }
    }


}
