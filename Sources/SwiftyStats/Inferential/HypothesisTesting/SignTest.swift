//
//  Created by VT on 06.12.25.
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
    /// Paired sign test comparing two matched samples.
    ///
    /// Computes the counts of positive/negative paired differences, discards ties,
    /// and returns:
    /// - two-sided exact p-value via the Binomial(n, 0.5) tail for non-tied pairs,
    /// - two-sided normal approximation with continuity correction,
    /// - counts of positive, negative, and tied pairs.
    ///
    /// - Parameters:
    ///   - set1: First paired sample container.
    ///   - set2: Second paired sample container (must match `set1` size and order).
    /// - Returns: ``SignTestRestult`` with exact and approximate p-values when inputs are valid; otherwise `nil` if sizes mismatch or data are non-numeric.
    public static func signTest<S: Comparable & Hashable & Codable & Sendable>(set1: SSExamine<S, T>, set2: SSExamine<S, T>) throws -> SignTestRestult<T>? {
        guard set1.sampleSize >= 2, set2.sampleSize >= 2, set1.sampleSize == set2.sampleSize, set1.isNumeric, set2.isNumeric else { return nil}
        let s1 = set1.itemsAsArray(sorted: .original)
        let s2 = set2.itemsAsArray(sorted: .original)
        var np: Int = 0
        var nn: Int = 0
        var nties: Int = 0
        var i: Int = 0
        while i < s1.count {
            if s2[i] > s1[i] {
                np += 1
            } else if s2[i] < s1[i] {
                nn += 1
            } else {
                nties += 1
            }
            i += 1
        }
        let effectiveSample = np + nn
        guard effectiveSample > 0 else {
            return SignTestRestult(pValueExact: T.one, pValueApprox: nil, nPosDiff: np, nNegDiff: nn, nTies: nties, total: effectiveSample, zStat: nil)
        }
        guard let pExactTwoSided = try Parametric<T>.binomialTest(successes: np, numberOfTrials: effectiveSample, p0: .half, alternative: .twoSided) else {
            return nil
        }
        let nnpnp: T = T(effectiveSample)
        let diff = T(np) - (.half * nnpnp)
        let continuity: T = diff >= .zero ? .half : -.half
        let std = T.sqrt(nnpnp * T(0.25))
        let z: T = (diff - continuity) / std
        let dist: Distribution.Normal<T> = try .init()
        let absZ = z < .zero ? -z : z
        let pAsympOneSided: T = T.one - (try dist.cdf(absZ))
        let pAsymp: T = min(T.one, T.two * pAsympOneSided)
        return SignTestRestult(pValueExact: pExactTwoSided, pValueApprox: pAsymp, nPosDiff: np, nNegDiff: nn, nTies: nties, total: effectiveSample, zStat: z)
    }
    
    /// Convenience overload accepting raw arrays for the paired sign test.
    ///
    /// - Parameters:
    ///   - set1: First paired sample.
    ///   - set2: Second paired sample.
    /// - Returns: ``SignTestRestult`` when the paired sign test can be evaluated; otherwise `nil` if sizes mismatch or data are non-numeric.
    public static func signTest<S: Comparable & Hashable & Codable & Sendable>(set1: [S], set2: [S]) throws -> SignTestRestult<T>? {
        let ex1: SSExamine<S, T> = .init(usingArray: set1, name: "set1", characterSet: nil)
        let ex2: SSExamine<S, T> = .init(usingArray: set2, name: "set2", characterSet: nil)
        return try signTest(set1: ex1, set2: ex2)
    }
    
    /// One-sample sign test against a median of zero.
    ///
    /// Interprets each observation's sign (positive, negative, zero) relative to `0`,
    /// discards ties, and returns the exact binomial two-sided p-value along with the
    /// continuity-corrected normal approximation.
    ///
    /// - Parameter data: Sample container to test; values must be numeric.
    /// - Returns: ``SignTestRestult`` with counts and p-values; `nil` if the sample is too small or non-numeric.
    public static func signTest<S: Comparable & Hashable & Codable & Sendable>(data: SSExamine<S, T>) throws -> SignTestRestult<T>? {
        guard data.sampleSize >= 2, data.isNumeric else { return nil}
        guard let s1 = data.itemsAsNumericArray else { return nil }
        var np: Int = 0
        var nn: Int = 0
        var nties: Int = 0
        var i: Int = 0
        while i < s1.count {
            if s1[i] > .zero {
                np += 1
            } else if s1[i] < .zero {
                nn += 1
            } else {
                nties += 1
            }
            i += 1
        }
        let effectiveSample = np + nn
        guard effectiveSample > 0 else {
            return SignTestRestult(pValueExact: T.one, pValueApprox: nil, nPosDiff: np, nNegDiff: nn, nTies: nties, total: effectiveSample, zStat: nil)
        }
        guard let pExactTwoSided = try Parametric<T>.binomialTest(successes: np, numberOfTrials: effectiveSample, p0: .half, alternative: .twoSided) else {
            return nil
        }
        let nnpnp: T = T(effectiveSample)
        let diff = T(np) - (.half * nnpnp)
        let continuity: T = diff >= .zero ? .half : -.half
        let std = T.sqrt(nnpnp * T.quarter)
        let z: T = (diff - continuity) / std
        let dist: Distribution.Normal<T> = try .init()
        let absZ = z < .zero ? -z : z
        let pAsympOneSided: T = T.one - (try dist.cdf(absZ))
        let pAsymp: T = min(T.one, T.two * pAsympOneSided)
        return SignTestRestult(pValueExact: pExactTwoSided, pValueApprox: pAsymp, nPosDiff: np, nNegDiff: nn, nTies: nties, total: effectiveSample, zStat: z)
    }
    
    /// Convenience overload for the one-sample sign test from a raw array.
    ///
    /// - Parameter data: Sample to test against a median of zero.
    /// - Returns: ``SignTestRestult`` mirroring the SSExamine-based overload; `nil` when the test cannot be evaluated.
    public static func signTest<S: Comparable & Hashable & Codable & Sendable>(data: Array<S>) throws -> SignTestRestult<T>? {
        let ex: SSExamine<S, T> = .init(usingArray: data, name: nil, characterSet: nil)
        return try signTest(data: ex)
    }
    
    /// Sign test results
    public struct SignTestRestult<FP: RealLike>: CustomStringConvertible, Codable {
        /// Exact two-sided p-value from the Binomial(n, 0.5) tail (ties removed).
        public var pValueExact: FP?
        /// Two-sided normal approximation with continuity correction.
        public var pValueApprox: FP?
        /// Count of positive paired differences (or positive signs vs 0).
        public var nPosDiff: Int?
        /// Count of negative paired differences (or negative signs vs 0).
        public var nNegDiff: Int?
        /// Count of ties (differences equal to 0).
        public var nTies: Int?
        /// Effective sample size used by the test (n = nPosDiff + nNegDiff).
        public var total: Int?
        /// z statistic used for the normal approximation (if available).
        public var zStat: FP?
        /// Human-readable multiline summary of results.
        public var description: String {
            get {
                var descr = String()
                if let pe = self.pValueExact, let pa = self.pValueApprox, let np = self.nPosDiff, let nn = self.nNegDiff, let nt = self.nTies, let n = self.total, let z = self.zStat {
                    descr.append("SIGN TEST\n")
                    descr.append("*********\n")
                    descr.append("p-value exact: \(niceNumber(pe))\n")
                    descr.append("p-value asymp: \(niceNumber(pa))\n")
                    descr.append("z: \(niceNumber(z))\n")
                    descr.append("count of \"+\": \(np)\n")
                    descr.append("count of \"-\": \(nn)\n")
                    descr.append("count of ties: \(nt)\n")
                    descr.append("n: \(n)\n")
                    descr.append("H1: true median is equal to 0\n")
                }
                return descr
            }
        }
    }
}
