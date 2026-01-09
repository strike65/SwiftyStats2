//
//  Created by VT on 05.12.25.
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
    
   
    /// Mann–Whitney U test (Wilcoxon rank-sum) for two independent samples.
    ///
    /// - Parameters:
    ///   - group1: First sample.
    ///   - group2: Second sample.
    ///   - continuityCorrection: Apply a ±0.5 continuity correction to the normal
    ///     approximation (exact p-values are unaffected).
    /// - Returns: Test result (U, W, z, exact/asymptotic p-values, effect size) or
    ///   `nil` when either sample is empty.
    public static func mannWhitneyUTest<S: Codable & Sendable & Comparable & Hashable>(group1: SSExamine<S, T>, group2: SSExamine<S, T>, continuityCorrection: Bool = false) throws -> MannWhitneyUTestResult<T>? {
        guard group1.sampleSize >= 1, group2.sampleSize >= 1 else {
            return nil
        }
        let groups: [String] = Array<String>.init(repeating: "group 1", count: group1.sampleSize) + Array<String>.init(repeating: "group 2", count: group2.sampleSize)
        let data = group1.itemsAsArray(sorted: .original) + group2.itemsAsArray(sorted: .original)
        let gd = GroupedData<S, String>(data: data, groups: groups)
        let (sortedData, sortedGroups) = gd.sorted()
        let rankInfo = Rank<S, String, T>(sortedData: sortedData, sortedGroups: sortedGroups)
        let sumRanksGroup1: T = rankInfo.sumOfRanks(g: "group 1")
        let sumRanksGroup2: T = rankInfo.sumOfRanks(g: "group 2")
        let n1: T = T(rankInfo.sampleSize(g: "group 1"))
        let n2: T = T(rankInfo.sampleSize(g: "group 2"))
        let mn: T = n1 * n2
        let S: T = n1 + n2
        let U1: T = mn + n1 * (n1 + .one) / .two - sumRanksGroup1
        let U2: T = mn + n2 * (n2 + .one) / .two - sumRanksGroup2
        if (U1 + U2) != mn {
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Wrong")
        }
        let dist: Distribution.Normal<T> = try .init(mean: .zero, sd: .one)
        let U: T = U1 > mn / .two ? mn - U1 : U1
        let meanU: T = mn / .two
        let tieAdjustment: T = rankInfo.totalTieCorrection
        let variance:T = (mn / T(12)) * ((S + T.one) - tieAdjustment / (S * (S - T.one)))
        let sd = T.sqrt(T.maximum(0, variance))
        let diff = U - meanU
        let continuityShift: T = continuityCorrection ? (diff >= .zero ? -(.half) : (.half)) : .zero
        let z = variance.isZero ? .zero : (diff + continuityShift) / sd
        let lowerTail = try dist.cdf(z)
        let pValueAsymptotic = variance.isZero ? (U == meanU ? .one : .zero ) : T.minimum(.one - lowerTail, lowerTail)
        var pValueExact: T = .nan
        let term3 = mn + T.minimum(n1, n2)
        if rankInfo.numberOfTies == 0 && mn < 400 && term3 <= 220 {
            pValueExact = try MMannWhithneyDistribution.pValue(U: Int(U), n1: Int(n1), n2: Int(n2))
        }
        let W = rankInfo.sumOfRanks(g: "group 1")
        let res: MannWhitneyUTestResult<T> = MannWhitneyUTestResult<T>(sumRanks1: rankInfo.sumOfRanks(g: "group 1"),
                                                                       sumRanks2: rankInfo.sumOfRanks(g: "group 2"),
                                                                       meanRank1: rankInfo.meanRank(g: "group 1"),
                                                                       meanRank2: rankInfo.meanRank(g: "group 2"),
                                                                       sumTiedRanks: rankInfo.totalTieCorrection,
                                                                       nTies: rankInfo.numberOfTies,
                                                                       U: U,
                                                                       W: W,
                                                                       zStat: z.isNaN ? .zero : z,
                                                                       pValueAsymTwoSided: T.minimum(.one, pValueAsymptotic * .two),
                                                                       pValueExactTwoSided: pValueExact.isNaN ? pValueExact : T.minimum(.one, pValueExact * .two),
                                                                       pValueAsym: pValueAsymptotic,
                                                                       pValueExact: pValueExact,
                                                                       effectSize: (!z.isNaN && !variance.isZero) ? z / T.sqrt(S) : .zero,
                                                                       continuityCorrected: continuityCorrection)
        return res
    }
    
    /// Convenience overload accepting raw arrays.
    public static func mannWhitneyUTest<S: Codable & Sendable & Comparable & Hashable>(group1: [S], group2: [S], continuityCorrection: Bool = false) throws -> MannWhitneyUTestResult<T>? {
        let g1: SSExamine<S, T> = .init(usingArray: group1, name: "group 1", characterSet: nil)
        let g2: SSExamine<S, T> = .init(usingArray: group2, name: "group 2", characterSet: nil)
        return try mannWhitneyUTest(group1: g1, group2: g2, continuityCorrection: continuityCorrection)
    }
    
    /// Sample size planning for Mann–Whitney U given target power.
    ///
    /// - Parameters:
    ///   - effectSize: Probability that a random draw from group1 exceeds group2 (Pr[X > Y]).
    ///   - alpha: Significance level.
    ///   - alternative: Test alternative (`.twoSided`/`.less`/`.greater`).
    ///   - power: Desired power (1 - β).
    /// - Returns: Tuple of (power, n1, n2) using equal allocation (c = 0.5), or `nil` when inputs are zeroed.
    public static func nMannWhitney(effectSize: T, alpha: T, alternative: Alternative = .twoSided, power: T) throws -> (power: T, n1: Int, n2: Int)? {
        return try nPowerU(effectSize: effectSize, alpha: alpha, alternative: alternative, power: power)
    }

    /// Power calculation for Mann–Whitney U given total sample size.
    ///
    /// - Parameters:
    ///   - effectSize: Probability that a random draw from group1 exceeds group2 (Pr[X > Y]).
    ///   - alpha: Significance level.
    ///   - alternative: Test alternative (`.twoSided`/`.less`/`.greater`).
    ///   - N: Total sample size (split 50/50).
    /// - Returns: Tuple of (power, n1, n2) using equal allocation, or `nil` when inputs are zeroed.
    public static func powerMannWhitney(effectSize: T, alpha: T, alternative: Alternative = .twoSided, N: Int) throws -> (power: T, n1: Int, n2: Int)? {
        return try nPowerU(effectSize: effectSize, alpha: alpha, alternative: alternative, N: N)
    }

    
    fileprivate static func nPowerU(effectSize: T, alpha: T = .half, alternative: Alternative = .twoSided, power: T = 0, N: Int = 0) throws -> (power: T, n1: Int, n2: Int)? {
        guard N != 0 || power != 0 else {
            return nil
        }
        var sigLevel: T = .zero
        if alternative != .twoSided {
            sigLevel = .two * alpha
        }
        else {
            sigLevel = alpha
        }
        let c: T = .half
        let nd: Distribution.Normal<T> = try .init(mean: .zero, sd: .one)
        var zBeta: T
        var zAlpha2: T
        var num: T, den: T
        var NN: T = .zero
        var n1: Int = 0, n2: Int = 0
        var PP: T = .zero
        if N == 0 {
            zBeta = try nd.quantile(.one - power)
            zAlpha2 = try nd.quantile(sigLevel / 2)
            num = T.pow(zBeta + zAlpha2, 2)
            den = T(12) * c * (.one - c) * T.pow(effectSize - .half, 2)
            NN = num / den
            n1 = Int((NN * c).rounded(.up))
            n2 = Int((NN * (1 - c)).rounded(.up))
            PP = power
        }
        else if power == 0 {
            let sigma0 = T.sqrt(.one / (T(12) * c * (.one - c)))
            zAlpha2 = try nd.quantile(sigLevel / 2)
            let delta: T = (effectSize - 0.5).magnitude
            let z: T = (T.sqrt(T(N)) * delta + zAlpha2 * sigma0) / sigma0
            PP = try nd.cdf(z)
            n1 = Int( (T(N) * c).rounded(.up) )
            n2 = N - n1
        }
        return (PP, n1, n2)
    }
    
    internal struct MMannWhithneyDistribution {
        private static func pmf(n1: Int, n2: Int) throws -> [T] {
            precondition(n1 > 0 && n2 > 0)
            let m = min(n1, n2)
            let n = max(n1, n2)
            let maxU = m * n
            
            // dynamic programming on counts N(i,j,u)
            // We store only two "layers" in j: prev (for j-1) and curr (for j)
            var prev = Array(repeating: Array<T>(repeating: .zero, count: maxU + 1), count: m + 1)
            for i in 0...m {
                prev[i][0] = .one
            }
            
            if n == 0 {
                return prev[m]
            }
            
            // Recurrence over j = 1,...,n
            // N(i,j,u) = N(i, j-1, u) + N(i-1, j, u - j)
            for j in 1...n {
                var curr = Array(repeating: Array<T>(repeating: .zero, count: maxU + 1), count: m + 1)
                curr[0][0] = .one
                for i in 1...m {
                    let maxUij = i * j
                    for u in 0...maxUij {
                        var value = prev[i][u]
                        if u >= j {
                            value += curr[i - 1][u - j]
                        }
                        curr[i][u] = value
                    }
                }
                prev = curr
            }
            let totalComb = T(try SpecialFunctions.binomial_coeff(UInt32(n + m), UInt32(m)))
            var pmf: [T] = prev[m]
            let invTotal = T.one / totalComb
            for u in 0...maxU {
                pmf[u] *= invTotal
            }
            return pmf
        }
        
        private static func cdf(n1: Int, n2: Int) throws -> [T] {
            let p = try pmf(n1: n1, n2: n2)
            var cdf = p
            var running: T = .zero
            for u in 0..<p.count {
                running += p[u]
                cdf[u] = running
            }
            if let last = cdf.last, (last - .one).magnitude > T.ulpOfOne {
                cdf[cdf.count - 1] = .one
            }
            return cdf
        }
        
        /// Exact tail probability for a Mann-Whitney `U` statistic.
        ///
        /// The distribution is computed exactly via dynamic programming over the number of
        /// interleavings of ranks (equivalently: the Wilcoxon rank-sum distribution).
        ///
        /// - Parameters:
        ///   - U: Observed `U` statistic, constrained to `0 ... n1 * n2`.
        ///   - n1: Sample size of the first group.
        ///   - n2: Sample size of the second group.
        ///   - alternative: Tail definition to apply (`less`, `greater`, or `twoSided`).
        /// - Returns: The requested exact p-value under the null hypothesis.
        /// - Throws: ``SSError`` if binomial coefficients cannot be computed.
        public static func pValue(U: Int, n1: Int, n2: Int, alternative: Alternative = .twoSided) throws -> T {
            precondition(U >= 0)
            let m = min(n1, n2)
            let n = max(n1, n2)
            let maxU = m * n
            precondition(U <= maxU)
            let cdf = try self.cdf(n1: m, n2: n)
            let u = U
            switch alternative {
            case .less:
                return cdf[u]
            case .greater:
                if u == 0 {
                    return .one
                }
                else {
                    return max(.zero, .one - cdf[u - 1])
                }
            case .twoSided:
                let left = cdf[u]
                let right: T
                if u == 0 {
                    right = .one
                }
                else {
                    right = T.maximum(.zero, .one - cdf[u - 1])
                }
                let p = .two * T.minimum(left, right)
                return T.minimum(.one, p)
            }
        }
    }
    
    /// Wilcoxon signed-rank test for paired samples.
    ///
    /// - Parameters:
    ///   - tier1: First paired sample.
    ///   - tier2: Second paired sample (same length as `tier1`).
    ///   - zeroHandling: Zero-difference strategy (`discardZeroDifferences` for classic Wilcoxon or `pratt` to keep zeros with rank zero and tie correction).
    ///   - cc: Apply a ±0.5 continuity correction to the normal approximation.
    ///   - exact: Compute the exact two-sided p-value when eligible.
    ///   - exactMaxN: Maximum non-zero paired count allowed for exact enumeration.
    /// - Returns: Test result (z, approximate and optional exact p-values, effect size, rank counts) or `nil` when all differences are zero.
    public static func wilcoxonSignedRankTest(tier1: SSExamine<T, T>, tier2: SSExamine<T, T>, zeroHandling: WilcoxonZeroHandling = .discardZeroDifferences, continuityCorrection cc: Bool = true, exact: Bool = true, exactMaxN: Int = 30) throws -> WilcoxonMatchedPairsTestResult<T>? {
        precondition(tier1.sampleSize >= 2 && tier1.sampleSize == tier2.sampleSize && tier1.isNumeric && tier2.isNumeric, "Both datasets must have the same size and must contain numeric values.")
        var i: Int = 0
        var temp: T
        var diff: [T] = []
        var signs: [T] = []
        var zeroDiffCount: Int = 0
        let a1: [T] = tier1.itemsAsArray(sorted: .original)
        let a2: [T] = tier2.itemsAsArray(sorted: .original)
        while i < tier1.sampleSize {
            temp = a1[i] - a2[i]
            if temp.isZero {
                zeroDiffCount += 1
                if zeroHandling == .discardZeroDifferences {
                    i += 1
                    continue
                }
            }
            diff.append(temp)
            i += 1
        }
        if diff.count == 0 {
            return nil
        }
        let diffSorted: [T] = diff.sorted(by: { $0.magnitude < $1.magnitude })
        var absDiffSorted: [T] = []
        for i in 0..<diffSorted.count {
            if diffSorted[i] > T.zero {
                signs.append(T.one)
            }
            else if diffSorted[i] < T.zero {
                signs.append(-T.one)
            }
            else {
                signs.append(.zero)
            }
            absDiffSorted.append(diffSorted[i].magnitude)
        }
        let rankInfo: Rank<T, String, T> = .init(sortedData: absDiffSorted)
        let ranks = rankInfo.ranks
        var usedRanks = ranks
        if zeroHandling == .pratt {
            for idx in 0..<signs.count where signs[idx].isZero {
                usedRanks[idx] = .zero
            }
        }
        let tieBlocks = rankInfo.tieBlocks
        let n = absDiffSorted.count
        let nf: T = T(n)
        var nP: Int = 0, nN: Int = 0
        var sumPosRanks: T = 0, sumNegRanks: T = 0
        var meanPosRank: T = 0, meanNegRank: T = 0
        for i in 0..<n {
            let rankValue = usedRanks[i]
            if signs[i] == T.one {
                nP += 1
                sumPosRanks += rankValue
            }
            else if signs[i] == -T.one {
                nN += 1
                sumNegRanks += rankValue
            }
        }
        guard (nP + nN) > 0 else {
            return nil
        }
        let nnp1: T = nf * (nf + T.one)
        let rankSumCheck: T = usedRanks.reduce(.zero, +)
        precondition(rankSumCheck.isApproximatelyEqual(to: sumNegRanks + sumPosRanks, absoluteTolerance: .ulpOfOne * T.maximum(.one, rankSumCheck)), "Something went wrong with the calculation of the expected number of ties.")
        meanNegRank = nN > 0 ? (sumNegRanks / T(nN)) : 0
        meanPosRank = nP > 0 ? (sumPosRanks / T(nP)) : 0
        var z: T
        var ts: T = .zero
        for block in tieBlocks {
            ts += block.correction / T(48)
        }
        let n1n21n1: T = nnp1 * (T.two * nf + T.one)
        let sigma: T = T.sqrt(n1n21n1 / T(24) - ts)
        let Tmin: T = T.minimum(sumNegRanks, sumPosRanks)
        let z0: T = Tmin - nnp1 / T(4)
        let correction: T = z0 < .zero ? -.half : .half
        z = cc ? (z0 - correction) / sigma : z0 / sigma
        let dist: Distribution.Normal<T> = try Distribution.Normal<T>()
        let pp: T = try dist.cdf(z)
        let pUpper: T = .one - pp
        let p2Normal: T = .two * .minimum(pp, pUpper)
        let cohenD: T = z.magnitude / T.sqrt(.two * T(nP + nN))
        var result: WilcoxonMatchedPairsTestResult<T> = WilcoxonMatchedPairsTestResult<T>()
        result.p2Value = p2Normal
        result.sampleSize = tier1.sampleSize
        result.nPosRanks = nP
        result.nNegRanks = nN
        result.nTies = rankInfo.numberOfTies
        result.nZeroDiff = zeroDiffCount
        result.sumNegRanks = sumNegRanks
        result.sumPosRanks = sumPosRanks
        result.meanNegRank = meanNegRank
        result.meanPosRank = meanPosRank
        result.cohenD = cohenD
        result.zStat = z
        result.zeroHandling = zeroHandling.description
        if exact && n <= exactMaxN {
            var ranksScaled = [Int]()
            ranksScaled.reserveCapacity(n)
            for r in usedRanks {
                let scaled = .two * r
                let intScaled: Int = Helpers.integerValue(scaled.rounded())
                ranksScaled.append(intScaled)
            }
            var TplusScaled = 0
            for i in 0..<n {
                if signs[i] == .one {
                    TplusScaled += ranksScaled[i]
                }
            }
            let pmfScaled = Self.wilcoxonExactPMFScaled(ranksScaled: &ranksScaled)
            result.p2ValueExact = Self.wilcoxonTwoSidedExactScaled(Tobs: TplusScaled, pmf: pmfScaled)
            result.Tplus = T(TplusScaled) * .half
        }
        return result
    }
    
    /// Convenience overload accepting raw arrays for the Wilcoxon signed-rank test.
    ///
    /// - Parameters:
    ///   - tier1: First paired sample.
    ///   - tier2: Second paired sample.
    ///   - zeroHandling: Zero-difference strategy (`discardZeroDifferences` or Pratt).
    ///   - cc: Apply a ±0.5 continuity correction to the normal approximation.
    ///   - exact: Compute the exact two-sided p-value when eligible.
    ///   - exactMaxN: Maximum non-zero paired count to allow exact enumeration.
    /// - Returns: Test result or `nil` when all differences are zero.
    public static func wilcoxonSignedRankTest(tier1: [T], tier2: [T], zeroHandling: WilcoxonZeroHandling = .discardZeroDifferences, continuityCorrection cc: Bool = true, exact: Bool = true, exactMaxN: Int = 30) throws -> WilcoxonMatchedPairsTestResult<T>? {
        let ex1: SSExamine<T, T> = .init(usingArray: tier1, name: "tier1", characterSet: nil)
        let ex2: SSExamine<T, T> = .init(usingArray: tier2, name: "tier2", characterSet: nil)
        return try wilcoxonSignedRankTest(tier1: ex1, tier2: ex2, zeroHandling: zeroHandling, continuityCorrection: cc, exact: exact, exactMaxN: exactMaxN)
    }
    private static func wilcoxonExactPMFScaled(ranksScaled: inout [Int]) -> [T] {
        let m = ranksScaled.count
        if m == 0 {
            return []
        }
        let maxSumScaled: Int = ranksScaled.reduce(0, +)
        var ways = [Int](repeating: 0, count: maxSumScaled + 1)
        ways[0] = 1
        var nonZeroRanksCount: Int = 0
        for r in ranksScaled {
            guard r > 0 else { continue }
            nonZeroRanksCount += 1
            if maxSumScaled >= r {
                for s in stride(from: maxSumScaled - r, through: 0, by: -1) {
                    let c = ways[s]
                    if c != 0 {
                        ways[s + r] &+= c
                    }
                }
            }
        }
        
        let totalPatterns = T(1 << nonZeroRanksCount)
        var pmf = [T](repeating: .zero, count: maxSumScaled + 1)
        for s in 0...maxSumScaled {
            if ways[s] > 0 {
                pmf[s] = T(ways[s]) / totalPatterns
            }
        }
        return pmf
    }
    
    private static func wilcoxonOneSidedExactScaled(Tobs: Int, pmf: [T], direction: Alternative) -> T {
        let maxS = pmf.count - 1
        var p = T.zero
        switch direction {
        case .greater:
            if Tobs <= maxS {
                for s in Tobs...maxS {
                    p += pmf[s]
                }
            }
        case .less:
            let upper = min(Tobs, maxS)
            if upper >= 0 {
                for s in 0...upper {
                    p += pmf[s]
                }
            }
        case .twoSided:
            return .nan
        }
        return p
    }
    
    private static func wilcoxonTwoSidedExactScaled(Tobs: Int, pmf: [T]) -> T {
        let maxS = pmf.count - 1
        let center = T.half * T(maxS)
        let distObs = (T(Tobs) - center).magnitude
        var p = T.zero
        for s in 0...maxS {
            let dist = (T(s) - center).magnitude
            if dist >= distObs - (.two * T.ulpOfOne) {
                p += pmf[s]
            }
        }
        return T.minimum(p, .one)
    }
    
/// Holds the results of the Wilcoxon test for matched pairs
    public struct WilcoxonMatchedPairsTestResult<FP: RealLike>: CustomStringConvertible, Codable {
        /// two sided p-value
        public var p2Value: FP?
        /// exact two sided p-value
        public var p2ValueExact: FP?
        /// sample size
        public var sampleSize: Int?
        /// number of ranks > 0
        public var nPosRanks: Int?
        /// number of ranks < 0
        public var nNegRanks: Int?
        /// number of ties
        public var nTies: Int?
        /// number of zero valued differences
        public var nZeroDiff: Int?
        /// sum of negative ranks
        public var sumNegRanks: FP?
        /// sum of positive ranks
        public var sumPosRanks: FP?
        /// mean of positive ranks
        public var meanPosRank: FP?
        /// mean of negative ranks
        public var meanNegRank: FP?
        /// z statistic
        public var zStat: FP?
        /// Cohen's d
        public var cohenD: FP?
        /// Sum of positive ranks (`T+`) after applying the selected zero-difference handling.
        public var Tplus: FP?
        /// Zero-difference handling strategy used for the test (serialised form).
        public var zeroHandling: String?
        /// Returns a description
        public var description: String {
            get {
                var descr = String()
                if let p2 = self.p2Value, let p2e = self.p2ValueExact, let n = self.sampleSize, let np = self.nPosRanks, let nn = self.nNegRanks, let nt = self.nTies, let nz = self.nZeroDiff, let snr = self.sumNegRanks, let spr = self.sumPosRanks, let mpr = self.meanPosRank, let mnr = self.meanNegRank, let z = self.zStat, let dc = self.cohenD, let tp = self.Tplus, let zh = self.zeroHandling {
                    descr.append("WILCOXON TEST FOR MATCHED PAIRS\n")
                    descr.append("*******************************\n")
                    descr.append("exact two sided p-Value: \(niceNumber(p2e))\n")
                    descr.append("approx. two sided p-value: \(niceNumber(p2))\n")
                    descr.append("z: \(niceNumber(z))\n")
                    descr.append("Cohen's d: \(niceNumber(dc))\n")
                    descr.append("sample size: \(n)\n")
                    descr.append("mean of ranks > 0: \(niceNumber(mpr))\n")
                    descr.append("mean of ranks < 0: \(niceNumber(mnr))\n")
                    descr.append("sum of ranks > 0: \(niceNumber(spr))\n")
                    descr.append("sum of ranks < 0: \(niceNumber(snr))\n")
                    descr.append("count of ranks > 0: \(np)\n")
                    descr.append("count of ranks < 0: \(nn)\n")
                    descr.append("count of tied ranks: \(nt)\n")
                    descr.append("count of zeros: \(nz)\n")
                    descr.append("sum of positive ranks after zero-handling: \(niceNumber(tp))\n")
                    descr.append("zero handling: \(zh)\n")
                }
                return descr
            }
        }
    }
    
    /// Holds the results of the Mann-Whitney U test
    public struct MannWhitneyUTestResult<FP: RealLike>: CustomStringConvertible, Codable {
        /// sum of ranks in set 1
        public var sumRanks1: FP?
        /// sum of ranks in set 2
        public var sumRanks2: FP?
        /// mean of ranks in set 1
        public var meanRank1: FP?
        /// mean of ranks in set 2
        public var meanRank2: FP?
        /// sum of tied ranks
        public var sumTiedRanks: FP?
        /// number of ties
        public var nTies: Int?
        /// U
        public var U: FP?
        /// Wilcoxon W
        public var W: FP?
        /// Z
        public var zStat: FP?
        /// two sided approximated p-value
        public var pValueAsymTwoSided: FP?
        /// two sided exact p-value
        public var pValueExactTwoSided: FP?
        /// one sided approximated p-value
        public var pValueAsym: FP?
        /// one sided exact p-value
        public var pValueExact: FP?
        /// effect size
        public var effectSize: FP?
        /// Whether the asymptotic z statistic applied a continuity correction (typically `0.5`).
        public var continuityCorrected: Bool?
        /// Returns a description
        public var description: String {
            get {
                var descr = String()
                if let sr1 = self.sumRanks1, let sr2 = self.sumRanks2, let mr1 = self.meanRank1, let mr2 = self.meanRank2, let str = self.sumTiedRanks, let nt = self.nTies, let U = self.U, let W = self.W, let z = self.zStat, let p1a = self.pValueAsym, let p2a = self.pValueAsymTwoSided, let p1e = self.pValueExact, let p2e = self.pValueExactTwoSided, let es = self.effectSize, let cr = self.continuityCorrected {
                    descr.append("MANN-WHITNEY U TEST\n")
                    descr.append("*******************\n")
                    descr.append("Mann-Whitney U: \(niceNumber(U))\n")
                    descr.append("Wilcoxon W: \(niceNumber(W))\n")
                    descr.append("z: \(niceNumber(z))\n")
                    descr.append("one sided exact p-value: \(niceNumber(p1e))\n")
                    descr.append("two sided exact p-value: \(niceNumber(p2e))\n")
                    descr.append("one sided asymp. p-value: \(niceNumber(p1a))\n")
                    descr.append("two sided asymp. p-value: \(niceNumber(p2a))\n")
                    descr.append("sum ranks group 1: \(niceNumber(sr1))\n")
                    descr.append("sum ranks group 2: \(niceNumber(sr2))\n")
                    descr.append("mean rank group 1: \(niceNumber(mr1))\n")
                    descr.append("mean rank group 2: \(niceNumber(mr2))\n")
                    descr.append("sum of tied ranks: \(niceNumber(str))\n")
                    descr.append("count of ties: \(nt)\n")
                    descr.append("effect size: \(niceNumber(es))\n")
                    descr.append("Continuity corrected: \(cr)\n")
                }
                return descr
            }
        }
    }
    
    /// Zero-difference handling for the Wilcoxon signed-rank test.
    ///
    /// - ``discardZeroDifferences`` matches the classical Wilcoxon approach by removing zero differences before ranking.
    /// - ``pratt`` keeps zero differences, assigns them rank zero, and includes them in tie corrections (aligns with R's Pratt option).
    public enum WilcoxonZeroHandling: String, CustomStringConvertible, Sendable, Codable {
        case discardZeroDifferences
        case pratt
        /// Human-readable label for the handling strategy.
        public var description: String {
            switch self {
            case .discardZeroDifferences: return "discardZeroDifferences (Wilcoxon)"
            case .pratt: return "Pratt"
            }
        }
    }
}
