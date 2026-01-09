//
//  Created by VT on 04.12.25.
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

/// Ranking utilities and helpers for ordered statistical operations.
///
/// - Includes: stable grouped sorting via ``GroupedData``, average ranks with
///   tie-aware handling, per-group rank aggregates, and per-block tie metadata
///   (`TieBlock`) capturing start index, length, value, and exact `t^3 - t`
///   correction.

import SwiftyStatsPrelude

/// A helper container for paired values and group labels, with a stable sort that
/// maintains the association between values and their groups.
///
/// Use this to ensure your values and groups are sorted consistently before
/// computing ranks with ``Rank``.
///
/// - Generics:
///   - Value: The comparable value type to sort.
///   - Group: The associated group identifier type to be permuted in lockstep.
///
/// Example:
/// ```swift
/// let gd = GroupedData(data: [3, 1, 2, 2], groups: ["A", "B", "A", "B"])
/// let (sortedData, sortedGroups) = gd.sorted()
/// // sortedData  == [1, 2, 2, 3]
/// // sortedGroups == ["B", "A", "B", "A"]
/// ```
public struct GroupedData<Value: Comparable, Group> {
    /// The values to be sorted.
    public var data: [Value]
    /// The group labels associated with each value.
    public var groups: [Group]
    
    /// Creates a grouped dataset.
    /// - Parameters:
    ///   - data: Values to be processed.
    ///   - groups: Group labels aligned with `data`.
    public init(data: [Value], groups: [Group]) {
        self.data = data
        self.groups = groups
    }
    
    struct Item<V: Comparable, G> {
        let value: V
        let group: G
        let index: Int
    }
    
    /// Returns a stable, nondecreasing sort of `data` along with the
    /// correspondingly permuted `groups`.
    ///
    /// If two values are equal, the element with the smaller original index
    /// comes first, preserving stability.
    ///
    /// - Returns: A tuple of `(sortedData, sortedGroups)`.
    public func sorted() -> (sortedData: [Value], sortedGroups: [Group]) {
        precondition(self.data.count == self.groups.count)
        let n = self.data.count
        if n <= 1 {
            return (self.data, self.groups)
        }
        var items: [Item<Value, Group>] = []
        items.reserveCapacity(n)
        for i in 0..<n {
            items.append(Item(value: self.data[i], group: self.groups[i], index: i))
        }
        items.sort { lhs, rhs in
            if lhs.value == rhs.value {
                return lhs.index < rhs.index
            }
            else {
                return lhs.value < rhs.value
            }
        }
        var sortedData: [Value] = []
        var sortedGroups: [Group] = []
        sortedData.reserveCapacity(n)
        sortedGroups.reserveCapacity(n)
        for item in items {
            sortedData.append(item.value)
            sortedGroups.append(item.group)
        }
        return (sortedData, sortedGroups)
    }
}


    /// A generic rank computation utility for already-sorted data, optionally grouped,
    /// with support for tie-aware average ranks, per-group summaries, and detailed
    /// per-block tie corrections.
///
/// `Rank` computes:
/// - Per-observation average ranks (with ties receiving the mean of their rank positions).
/// - Tie-correction terms (t^3 − t) for each tie block, useful in nonparametric tests
///   such as Wilcoxon–Mann–Whitney or Kruskal–Wallis.
/// - Optional group aggregations: sum of ranks, mean ranks, and sample sizes per group.
///
/// - Important: `sortedData` must be sorted in nondecreasing order. If you have unsorted data,
///   sort it first (e.g., with ``GroupedData`` if you also have groups).
///
/// - Generics:
///   - Value: Comparable, Hashable value type that is already sorted nondecreasingly.
///   - Group: A hashable group identifier (e.g., `Int`, `String`).
///   - FP: A floating-point–like type used for ranks and tie corrections (see `RealLike`).
///
/// Example:
/// ```swift
/// // Example data with ties and two groups
/// let rawData = [7.0, 3.0, 5.0, 5.0, 10.0, 3.0]
/// let rawGroups = ["A", "B", "A", "B", "A", "B"]
///
/// // Sort data and groups together
/// let gd = GroupedData(data: rawData, groups: rawGroups)
/// let (sortedData, sortedGroups) = gd.sorted()
///
/// // Compute ranks (using Double as FP here; Double conforms to RealLike in this project)
/// let rankInfo = Rank<Double, String, Double>(sortedData: sortedData, sortedGroups: sortedGroups)
///
/// // Per-observation ranks (aligned with sortedData)
/// let ranks = rankInfo.ranks
///
/// // Group summaries
/// let levels = rankInfo.groupLevels
/// let sums = rankInfo.sumOfRanksPerGroup
/// let means = rankInfo.meanRanksPerGroup
/// let sizes = rankInfo.sampleSizesPerGroup
///
/// // Tie information
/// let tieBlocks = rankInfo.numberOfTies
/// let tieCorrections = rankInfo.tieCorrections
/// let totalCorrection = rankInfo.totalTieCorrection
/// ```
///
/// Ungrouped usage:
/// ```swift
/// let data = [1, 2, 2, 4, 9]  // already sorted
/// let r = Rank<Int, Int, Double>(sortedData: data) // Group is unused; FP is Double
/// // r.ranks == [1, 2.5, 2.5, 4, 5]
/// // r.sumOfRanksPerGroup == [15.0]  // single overall group
/// // r.sampleSizesPerGroup == [5]
/// ```
public struct Rank<Value: Comparable & Hashable & Sendable & Codable,
                   Group: Hashable & Sendable & Codable,
                   FP: RealLike> {
    
    /// Average ranks for each observation, aligned with the input order of `sortedData`.
    public let ranks: [FP]
    
    /// Group labels for each observation, aligned with `sortedData` if provided, otherwise `nil`.
    public let groups: [Group]?
    
    /// Unique group levels in order of their first appearance in the input.
    ///
    /// This defines the indexing for ``sumOfRanksPerGroup``, ``meanRanksPerGroup``,
    /// and ``sampleSizesPerGroup``.
    public let groupLevels: [Group]
    
    /// Sum of ranks per group, aligned with ``groupLevels``.
    public let sumOfRanksPerGroup: [FP]
    
    /// Mean (average) rank per group, aligned with ``groupLevels``.
    public let meanRanksPerGroup: [FP]
    
    /// Sample size (count of observations) per group, aligned with ``groupLevels``.
    public let sampleSizesPerGroup: [Int]

    /// Tie-correction terms for each tie block, defined as t^3 − t where t is the tie frequency.
    public let tieCorrections: [FP]
    
    /// Detailed metadata for each tie block.
    public struct TieBlock: Sendable, Codable {
        /// Start index in `sortedData` where the tie block begins.
        public let startIndex: Int
        /// Number of tied observations in this block.
        public let length: Int
        /// Exact tie-correction term `t^3 − t` for this block.
        public let correction: FP
        /// The tied value shared by this block.
        public let value: Value
        
        /// Creates a tie-block descriptor for a contiguous block of equal values.
        ///
        /// - Parameters:
        ///   - startIndex: Start index in `sortedData` where the tie block begins.
        ///   - length: Number of tied observations in this block.
        ///   - correction: Tie correction term `t^3 - t` for `t = length`.
        ///   - value: The tied value shared by this block.
        public init(startIndex: Int, length: Int, correction: FP, value: Value) {
            self.startIndex = startIndex
            self.length = length
            self.correction = correction
            self.value = value
        }
    }
    
    /// Per-block tie information (one entry per tie block in `sortedData`).
    public let tieBlocks: [TieBlock]
    
    /// Number of tie blocks discovered in `sortedData`.
    public let numberOfTies: Int
    
    /// Sum of all tie-correction terms, useful in variance adjustments for rank-based tests.
    public var totalTieCorrection: FP {
        tieCorrections.reduce(FP.zero, +)
    }
    
    /// Creates rank information for already sorted data, optionally grouped.
    ///
    /// - Important: `sortedData` must be sorted in nondecreasing order. If you also have
    ///   group labels, they must be passed as `sortedGroups` in the same order.
    ///
    /// - Parameters:
    ///   - sortedData: Data sorted in nondecreasing order.
    ///   - sortedGroups: Optional group labels, same length and ordering as `sortedData`.
    public init(sortedData: [Value], sortedGroups: [Group]? = nil) {
        precondition(sortedGroups == nil || sortedGroups!.count == sortedData.count,
                     "Data and groups must have the same length.")
        
        // Compute ranks and tie corrections on the sorted data.
        let (computedRanks, tieBlocks) = Rank.computeRanks(forSorted: sortedData)
        self.ranks = computedRanks
        self.tieBlocks = tieBlocks
        self.tieCorrections = tieBlocks.map { $0.correction }
        self.numberOfTies = tieBlocks.count
        
        // Handle group aggregation.
        if let groups = sortedGroups {
            self.groups = groups
            
            // Build group index map: group value -> 0-based index in groupLevels.
            var groupIndex: [Group: Int] = [:]
            var levels: [Group] = []
            levels.reserveCapacity(groups.count)
            
            for g in groups {
                if groupIndex[g] == nil {
                    let idx = levels.count
                    levels.append(g)
                    groupIndex[g] = idx
                }
            }
            
            self.groupLevels = levels
            
            let gCount = levels.count
            var sumRanks = Array(repeating: FP.zero, count: gCount)
            var sampleSizes = Array(repeating: 0, count: gCount)
            
            // Single pass: accumulate rank sums and sample sizes per group.
            for (k, g) in groups.enumerated() {
                guard let idx = groupIndex[g] else {
                    // This should never happen if the map was built correctly.
                    continue
                }
                sumRanks[idx] += computedRanks[k]
                sampleSizes[idx] += 1
            }
            
            self.sumOfRanksPerGroup = sumRanks
            self.sampleSizesPerGroup = sampleSizes
            
            var meanRanks = Array(repeating: FP.zero, count: gCount)
            for i in 0..<gCount {
                if sampleSizes[i] > 0 {
                    meanRanks[i] = sumRanks[i] / FP(sampleSizes[i])
                } else {
                    meanRanks[i] = FP.zero
                }
            }
            self.meanRanksPerGroup = meanRanks
            
        } else {
            // Ungrouped case: treat all observations as one group.
            self.groups = nil
            self.groupLevels = []
            
            let totalSum = computedRanks.reduce(FP.zero, +)
            let n = computedRanks.count
            
            self.sumOfRanksPerGroup = [totalSum]
            self.sampleSizesPerGroup = [n]
            self.meanRanksPerGroup = n > 0 ? [totalSum / FP(n)] : [FP.zero]
        }
    }
    
    /// Returns the rank sum for group `g`.
    ///
    /// - Parameter g: A value from ``groupLevels``.
    /// - Returns: The rank sum for `g`, or `nan` when `g` is not present.
    public func sumOfRanks(g: Group) -> FP {
        if !self.groupLevels.isEmpty {
            if let idx = self.groupLevels.firstIndex(of: g) {
                return self.sumOfRanksPerGroup[idx]
            }
            else {
                return .nan
            }
        }
        else {
            return self.sumOfRanksPerGroup.first ?? .nan
        }
    }

    /// Returns the mean rank for group `g`.
    ///
    /// - Parameter g: A value from ``groupLevels``.
    /// - Returns: The mean rank for `g`, or `nan` when `g` is not present.
    public func meanRank(g: Group) -> FP {
        if !self.groupLevels.isEmpty {
            if let idx = self.groupLevels.firstIndex(of: g) {
                return self.meanRanksPerGroup[idx]
            }
            else {
                return .nan
            }
        }
        else {
            return self.meanRanksPerGroup.first ?? .nan
        }
    }
    
    /// Returns the sample size (count) for group `g`.
    ///
    /// - Parameter g: A value from ``groupLevels``.
    /// - Returns: The sample size for `g` as a floating-point value, or `nan` when `g` is not present.
    public func sampleSize(g: Group) -> FP {
        if !self.groupLevels.isEmpty {
            if let idx = self.groupLevels.firstIndex(of: g) {
                return FP(self.sampleSizesPerGroup[idx])
            }
            else {
                return .nan
            }
        }
        else {
            let v = self.sampleSizesPerGroup.first ?? 0
            return FP(v)
        }
    }


    
    
    // MARK: - Internal rank computation
    
    /// Computes average ranks for already sorted data and tie-block metadata.
    ///
    /// - Important: `data` must be sorted in nondecreasing order.
    ///
    /// - Parameter data: The sorted values to rank.
    /// - Returns: A tuple containing:
    ///   - `ranks`: average ranks per observation (same length as `data`).
    ///   - `tieBlocks`: metadata for each tie block including exact corrections.
    private static func computeRanks(forSorted data: [Value])
        -> (ranks: [FP], tieBlocks: [TieBlock])
    {
        let n = data.count
        if n == 0 {
            return ([], [])
        }
        
        var ranks: [FP] = []
        ranks.reserveCapacity(n)
        
        var tieBlocks: [TieBlock] = []
        tieBlocks.reserveCapacity(8)
        
        var i = 0
        while i < n {
            let start = i
            let value = data[start]
            
            // Advance j while equal values (tie block).
            var j = start + 1
            while j < n && data[j] == value {
                j += 1
            }
            let freq = j - start
            
            if freq == 1 {
                // Single observation: rank = position + 1.
                let rank = FP(start + 1)
                ranks.append(rank)
            } else {
                // Tie block of length freq.
                let correctionInt = freq * freq * freq - freq // t^3 - t in integer space
                let correction = FP(exactly: correctionInt) ?? FP(correctionInt)
                tieBlocks.append(TieBlock(startIndex: start,
                                          length: freq,
                                          correction: correction,
                                          value: value))
                
                // Average rank for positions (start+1)...(start+freq)
                // mean = start + (freq + 1) / 2
                let numerator = FP(2 * start + freq + 1)   // 2*start + freq + 1
                let avgRank = numerator / FP(2)
                
                for _ in 0..<freq {
                    ranks.append(avgRank)
                }
            }
            
            i = j
        }
        
        return (ranks, tieBlocks)
    }
}
