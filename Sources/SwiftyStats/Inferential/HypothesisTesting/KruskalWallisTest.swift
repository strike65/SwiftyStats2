//
//  Created by VT on 13.12.25.
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
import SwiftyStatsPrelude

extension Inferential.HypothesisTesting.NonParametricTests {
    
    /// Performs the Kruskal-Wallis one-way analysis of variance by ranks.
    ///
    /// This non-parametric test compares `k >= 2` independent groups using the pooled ranks.
    /// The (uncorrected) test statistic is
    /// `H = (12 / (N * (N + 1))) * sum_i (R_i^2 / n_i) - 3 * (N + 1)`,
    /// where `R_i` is the rank sum of group `i` and `N` is the total sample size.
    ///
    /// When ties are present, a tie-corrected statistic is reported using
    /// `Hc = H / (1 - sum(t^3 - t) / (N^3 - N))`, where `t` denotes a tie block size.
    ///
    /// - Parameters:
    ///   - data: One array per group. Each group must contain at least two observations.
    ///   - alpha: Significance level in `(0, 1)` used to compute the chi-square critical value.
    /// - Returns: A ``KruskalWallisHTestResult`` or `nil` when the input is degenerate (e.g. fewer than two valid groups).
    /// - Throws: ``SSError`` if the chi-square distribution cannot be constructed or evaluated.
    public static func kruskalWallisTest(data: [[T]], alpha: T) throws -> KruskalWallisHTestResult<T>? {
        guard data.count > 1 else {
            return nil
        }
        var groups: [Int] = []
        var a1: [T] = []
        var k: Int = 1
        var N: T = 0
        for a in data {
            guard a.count > 1 else {
                return nil
            }
            groups.append(contentsOf: Array<Int>.init(repeating: k, count: a.count))
            k += 1
            N += T(a.count)
            a1.append(contentsOf: a)
        }
        let gd = GroupedData<T, Int>(data: a1, groups: groups)
        let (sortedData, sortedGroups) = gd.sorted()
        let rankInfo = Rank<T, Int, T>(sortedData: sortedData, sortedGroups: sortedGroups)
        var sumOfRanks: T = 0
        var t1: T, t2: T, t3: T
        sumOfRanks = rankInfo.sumOfRanksPerGroup.reduce(.zero, +)
        t1 = N + .one
        t2 = N * t1
        t3 = t2 / .two
        if sumOfRanks != t3 {
            return nil
        }
        var rankSum: [T] = []
        for i in 1...rankInfo.groupLevels.count {
            rankSum.append(T.pow(rankInfo.sumOfRanks(g: i), 2) / rankInfo.sampleSize(g: i))
        }
        let sum: T = Helpers.sum(&rankSum)
        t1 = N * (N + .one)
        t2 = T(3) * (N + .one)
        t3 = T(12) / t1
        let H: T = t3 * sum - t2
        let df: T = T(rankInfo.groupLevels.count) - .one
        var ts: T = .zero
        ts = rankInfo.tieCorrections.reduce(.zero) { (acc, x) in
            return acc + x
        }
        ts = .one - (ts / T.pow(N, 3.0) - N)
        let Hc: T = ts != .zero ? H / ts : .nan
        let dist: Distribution.ChiSquared<T> = try .init(degreesOfFreedom: df)
        let p: T = .one - (try dist.cdf(H))
        let cv: T = try dist.quantile(.one - alpha)
        return KruskalWallisHTestResult<T>.init(H_value: H,
                                                H_value_corrected: Hc,
                                                pValue: p,
                                                nGroups: rankInfo.groupLevels.count,
                                                df: Int(df),
                                                nObservations: Int(N),
                                                meanRanks: rankInfo.meanRanksPerGroup,
                                                sumRanks: rankInfo.sumOfRanksPerGroup,
                                                cv: cv,
                                                nTies: rankInfo.numberOfTies,
                                                alpha: alpha)
    }
    
    /// The results of the H test
    public struct KruskalWallisHTestResult<FP: RealLike>: CustomStringConvertible, Codable {
        /// Chi
        public var H_value: FP?
        /// Chi square corrected for ties
        public var H_value_corrected: FP?
        /// one sided p-value
        public var pValue: FP?
        /// number of Groups
        public var nGroups: Int?
        /// Degrees of Freedom
        public var df: Int?
        /// number of observations
        public var nObservations: Int?
        /// array of mean ranks per group
        public var meanRanks: Array<FP>?
        /// array of rank sums per group
        public var sumRanks: Array<FP>?
        /// critical value at alpha
        public var cv: FP?
        /// Number of ties
        public var nTies: Int?
        /// alpha
        public var alpha: FP?
        /// Returns a description
        public var description: String {
            get {
                var descr = String()
                if let chi = self.H_value, let chic = self.H_value_corrected, let p = self.pValue, let ng = self.nGroups, let sdf = self.df, let no = self.nObservations, let mr = self.meanRanks, let sr = self.sumRanks, let scv = self.cv, let a = self.alpha, let nt = self.nTies {
                    descr.append("KRUSKAL-WALLIS H TEST\n")
                    descr.append("***************** ***\n")
                    descr.append("p-value: \(niceNumber(p))\n")
                    descr.append("critical value: \(niceNumber(scv))\n")
                    descr.append("Chi square: \(niceNumber(chi))\n")
                    descr.append("Chi square corrected for ties: \(niceNumber(chic))\n")
                    descr.append("count of groups: \(ng)\n")
                    descr.append("degrees of freedom: \(sdf)\n")
                    descr.append("count of observations: \(no)\n")
                    descr.append("count of ties: \(nt)\n")
                    descr.append("mean of ranks: [\(mr.map { niceNumber($0) }.joined(separator:", " ))]\n")
                    descr.append("sum of ranks: [\(sr.map { niceNumber($0) }.joined(separator:", "))]\n")
                    descr.append("alpha: \(niceNumber(a))\n")
                }
                return descr
            }
        }
    }
}
