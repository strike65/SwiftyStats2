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

/// Internal frequency and empirical distribution helpers for SSExamine.

import SwiftyStatsPrelude
import SwiftyBoost

// MARK: Frequencies

extension SSExamine {
    /// Updates cumulative frequency caches in a thread-safe manner.
    ///
    /// - Returns: `true` if caches are up-to-date or were successfully updated.
    /// - Thread safety: Acquires the lock.
    internal func updateCumulativeFrequencies() -> Bool {
        withLock {
            return updateCumulativeFrequenciesLocked()
        }
    }
    
    /// Rebuilds cumulative frequency and relative frequency caches if needed (locked variant).
    ///
    /// - Important:
    ///   - The caller must already hold the lock.
    ///   - A no-op if `needsUpdateFrequencies == false`.
    /// - Returns: `true` if caches are up-to-date or were successfully updated.
    /// - Complexity: O(k log k) for k unique elements (due to sorting).
    @inline(__always)
    internal func updateCumulativeFrequenciesLocked() -> Bool {
        guard needsUpdateFrequencies else { return true }
        cumulativeFrequenciesStorage.removeAll(keepingCapacity: true)
        cumulativeRelativeFrequenciesStorage.removeAll(keepingCapacity: true)
        guard count > 0 else {
            needsUpdateFrequencies = false
            return true
        }
        // Work on unique keys in ascending order while holding the lock.
        let keys = uniqueItemsLocked(sorted: .ascending)
        var runningFrequency: Int = 0
        var runningRelative: FP = .zero
        for key in keys {
            let f = items[key] ?? 0
            runningFrequency += f
            let rf: FP = FP(f) / FP(count)
            runningRelative += rf
            cumulativeFrequenciesStorage[key] = runningFrequency
            cumulativeRelativeFrequenciesStorage[key] = runningRelative
        }
        needsUpdateFrequencies = false
        return true
    }
    
    /// Builds a frequency table (absolute and relative) for all unique items (locked variant).
    ///
    /// - Important: The caller must already hold the lock.
    /// - Parameter sorted: Sorting mode for the resulting table.
    /// - Returns: An array of `FrequencyTableItem` entries.
    /// - Complexity: O(k) to build; O(k log k) when sorting.
    @inline(__always)
    internal func frequencyTableLocked(sorted: SortOrderFrequencyTable) -> Array<FrequencyTableItem<SSElement, FP>> {
        var result = Array<FrequencyTableItem<SSElement, FP>>()
        var tableItem: FrequencyTableItem<SSElement, FP>
        guard !self.isEmpty else { return result }
        let n = FP(self.sampleSize)
        for (item, freq) in self.items {
            let f = FP(freq)
            tableItem = FrequencyTableItem(item: item, relativeFrequency: f / n, absoluteFrequency: freq)
            result.append(tableItem)
        }
        switch sorted {
        case .none:
            return result
        case .valueAscending:
            return result.sorted(by: { $0.item < $1.item })
        case .valueDescending:
            return result.sorted(by: { $0.item > $1.item })
        case .frequencyAscending:
            return result.sorted(by: { $0.frequency < $1.frequency })
        case .frequencyDescending:
            return result.sorted(by: { $0.frequency > $1.frequency })
        }
    }
    
    /// Builds a cumulative frequency table according to the specified format (locked variant).
    ///
    /// - Important: The caller must already hold the lock.
    /// - Parameter format: Defines whether each observation or each unique item is represented.
    /// - Returns: An array of `CumulativeFrequencyTableItem` entries.
    /// - Complexity: O(n) for `eachItem` (expanded), O(k) for `eachUniqueItem`.
    ///
    /// Notes:
    /// - The cumulative counts are computed in ascending order of item values.
    /// - `.eachItem` expands entries proportionally to their absolute frequencies;
    ///   `.eachUniqueItem` lists each distinct item once.
    @inline(__always)
    internal func cumulativeFrequencyTableLocked(format: CumulativeFrequencyTableFormat) -> Array<CumulativeFrequencyTableItem<SSElement, FP>> {
        var tableItem: CumulativeFrequencyTableItem<SSElement, FP>
        var runningRelativeFrequency: FP = 0
        var runningAbsoluteFrequency: Int = 0
        var result = Array<CumulativeFrequencyTableItem<SSElement, FP>>()
        guard !self.isEmpty else { return result }
        let frequencyTable = frequencyTableLocked(sorted: .valueAscending)
        switch format {
        case .eachItem:
            for fItem: FrequencyTableItem<SSElement, FP> in frequencyTable {
                runningAbsoluteFrequency += fItem.frequency
                runningRelativeFrequency += fItem.relativeFrequency
                // Expand each absolute count into individual cumulative entries.
                for _ in (runningAbsoluteFrequency - fItem.frequency)..<(runningAbsoluteFrequency - 1) {
                    tableItem = CumulativeFrequencyTableItem<SSElement, FP>(item: fItem.item, cumulativeCount: runningAbsoluteFrequency, cumulativeRelativeCount: runningRelativeFrequency)
                    result.append(tableItem)
                }
            }
        case .eachUniqueItem:
            for fItem:FrequencyTableItem<SSElement, FP> in frequencyTable where fItem.frequency > 1 {
                runningAbsoluteFrequency += fItem.frequency
                runningRelativeFrequency += fItem.relativeFrequency
                tableItem = CumulativeFrequencyTableItem<SSElement, FP>(item: fItem.item, cumulativeCount: runningAbsoluteFrequency, cumulativeRelativeCount: runningRelativeFrequency)
                    result.append(tableItem)
            }
        }
        return result
    }
    

    
    /// Evaluates the empirical cumulative distribution function (ECDF) at a given item (locked variant).
    ///
    /// - Important: The caller must already hold the lock.
    /// - Parameter item: The query value for which to evaluate the ECDF.
    /// - Returns: The ECDF value as `FP` in `[0, 1]` when it can be determined; `nil` if the sample is empty.
    ///
    /// Discussion:
    /// - If `item` is strictly smaller than the current minimum, returns `0`.
    /// - If `item` is strictly greater than the current maximum, returns `1`.
    /// - Otherwise, attempts to return the stored cumulative relative frequency for `item`;
    ///   if not present (i.e., `item` not observed), returns the cumulative relative frequency
    ///   of the greatest observed value strictly less than `item`.
    @inline(__always)
    internal func empiricalCDFLocked(_ item: SSElement) -> FP? {
        guard !self.isEmpty else { return nil }
        _ = updateCumulativeFrequenciesLocked()
        guard let minItem = self.smallest, let maxItem = self.largest else { return nil }

        if item < minItem {
            return FP.zero
        }
        if item >= maxItem {
            return FP.one
        }

        if let exact = cumulativeRelativeFrequencies[item] {
            return exact
        }

        var result: FP = 0
        let orderedKeys = uniqueItemsLocked(sorted: .ascending)
        for key in orderedKeys where key < item {
            if let cumulative = cumulativeRelativeFrequencies[key] {
                result = cumulative
            }
        }
        return result
    }
}
