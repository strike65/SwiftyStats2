
//
//  Created by VT on 27.10.25.
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

/// Protocols describing the SSExamine API surface for containers and frequency views.

import SwiftyStatsPrelude

/// A generic container protocol for examining categorical or discrete data,
/// providing frequencies and relative frequencies over a multiset of items.
///
/// Conceptually, conformers expose a counting measure μ on a finite support S,
/// where μ(x) is the absolute frequency and μ(x)/n is the empirical probability
/// mass at x. Conforming types manage a collection of `ExamineItem` values where
/// items may appear multiple times (i.e. a bag/multiset). The protocol exposes
/// utilities to query counts, relative frequencies, and to mutate the contents.
///
/// - Generics:
///   - ExamineItem: The item type stored by the container. Must be `Hashable` for
///     counting and lookups, `Comparable` for ordered operations, and `Sendable`
///     for safe use with concurrency.
///   - FloatingPointValue: The floating-point type used for relative frequency
///     calculations (e.g. `Double`, `Float`, etc.). Must conform to `Real`,
///     `BinaryFloatingPoint`, and `Sendable`.
internal protocol ExamineContainer<ExamineItem, FloatingPointValue> where ExamineItem: Hashable & Comparable & Sendable, FloatingPointValue: RealLike {
    /// The item type stored in this container.
    associatedtype ExamineItem
    /// The floating-point type used when reporting relative frequencies.
    associatedtype FloatingPointValue

    /// The total number of observations currently stored, including duplicates.
    ///
    /// For example, if the underlying multiset contains ["A", "A", "B"],
    /// `sampleSize` is 3.
    var sampleSize: Int { get }

    /// The number of distinct items currently represented.
    ///
    /// Continuing the example ["A", "A", "B"], `length` is 2.
    var length: Int { get }

    /// A Boolean value indicating whether the container has no observations.
    var isEmpty: Bool { get }

    /// Returns whether at least one occurrence of the given item exists.
    ///
    /// - Parameter item: The item to look up.
    /// - Returns: `true` if the item occurs at least once; otherwise, `false`.
    mutating func contains(_ item: ExamineItem) -> Bool

    /// Returns the number of occurrences of the given item.
    ///
    /// - Parameter item: The item to count.
    /// - Returns: The frequency (non-negative integer) of `item`.
    mutating func frequency(_ item: ExamineItem) -> Int

    /// Returns the relative frequency (proportion) of the given item.
    ///
    /// The value is typically computed as `frequency(item) / FloatingPointValue(sampleSize)`.
    /// If `sampleSize` is zero, conformers should define a sensible behavior
    /// (e.g. return `0`).
    ///
    /// - Parameter item: The item to compute a relative frequency for.
    /// - Returns: The proportion of the item in the container.
    mutating func relativeFrequency(_ item: ExamineItem) -> FloatingPointValue

    /// Appends a single item occurrence to the container.
    ///
    /// - Parameter item: The item to append.
    mutating func append(_ item: ExamineItem)

    /// Appends the given item repeated a number of times.
    ///
    /// - Parameters:
    ///   - count: The number of times to append `item`. Implementations may treat
    ///     non-positive counts as a no-op and return `false`.
    ///   - item: The item to append.
    /// - Returns: `true` if at least one occurrence was appended; otherwise, `false`.
    mutating func append(repeating count: Int, item: ExamineItem) -> Bool

    /// Appends the contents of an array of items.
    ///
    /// - Parameter array: The items to append.
    /// - Returns: `true` if at least one item was appended; otherwise, `false`.
    mutating func append(contentOf array: Array<ExamineItem>) -> Bool

    /// Appends items parsed from a textual source.
    ///
    /// Conformers may interpret `text` and `characterSet` to derive `ExamineItem`s.
    /// For example, a character-based container may extract characters matching
    /// `characterSet`, or a token-based container may split on non-matching characters.
    ///
    /// - Parameters:
    ///   - text: The source text to parse.
    ///   - characterSet: An optional character set that influences parsing or filtering.
    /// - Returns: `true` if at least one item was appended; otherwise, `false`.
    /// - Throws: An error if parsing fails or if items cannot be constructed from `text`.
    mutating func append(text: String, characterSet: CharacterSet?) throws -> Bool

    /// Removes occurrences of a given item.
    ///
    /// - Parameters:
    ///   - item: The item to remove.
    ///   - allOccurences: If `true`, remove all occurrences of `item`; if `false`,
    ///     remove only a single occurrence (if present).
    mutating func remove(_ item: ExamineItem, allOccurences: Bool)

    /// Removes all items and resets the container to the empty state.
    mutating func removeAll()
}

/// A protocol describing a tabular container of named examine-series (columns),
/// similar to a lightweight data frame.
///
/// Conforming types manage a collection of homogeneous series (each of type `Examine`)
/// that share a common sample size (number of rows/observations).
internal protocol SSDataFrameContainer {
    /// The element/series type stored in each column.
    associatedtype Examine

    /// The number of columns currently stored.
    var nCols: Int { get }

    /// The number of rows/observations shared across columns.
    ///
    /// Conformers may require that all appended `Examine` instances agree on
    /// their sample size; this value reflects that common length.
    var totalSampleSize: Int { get  }

    /// A Boolean value indicating whether no columns are stored.
    var isEmpty: Bool { get }

    /// Appends a new column (series) to the container.
    ///
    /// - Parameters:
    ///   - examine: The series to append.
    ///   - name: An optional column name. Conformers may enforce uniqueness.
    /// - Throws: An error if the column cannot be appended (e.g. name collision
    ///   or incompatible sample size).
    mutating func append(_ examine: Examine, name: String?) throws

    /// Removes and returns the column associated with a name.
    ///
    /// - Parameter name: The name of the column to remove.
    /// - Returns: The removed series if it existed; otherwise, `nil`.
    mutating func remove(name: String!) -> Examine?

    /// Removes all columns and resets the container to the empty state.
    mutating func removeAll()
}
