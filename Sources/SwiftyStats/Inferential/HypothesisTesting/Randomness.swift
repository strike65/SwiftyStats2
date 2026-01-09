//
//  Created by VT on 01.12.25.
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

/// Randomness and runs test helpers for binary and continuous data streams.

import SwiftyStatsPrelude

extension Inferential.HypothesisTesting.RandomnessTests {
    /// Performs the Wald–Wolfowitz runs test to assess whether an observed
    /// sequence follows a random ordering around a cutting point.
    ///
    /// Values equal to the cutting point are ignored, the remaining samples are
    /// converted into sign changes, and the observed run count is compared
    /// against the normal approximation and exact distribution.
    /// - Parameters:
    ///   - data: Sample container that supplies the ordered observations and
    ///     summary statistics for the cutting point.
    ///   - cp: Strategy for choosing the split point (median, mean,
    ///     mode, or user-provided value).
    ///   - alpha: Two-sided significance level for the approximate normal test.
    ///   - coco: Whether to apply a 0.5 continuity
    ///     correction to the z statistic.
    /// - Returns: A ``RunsTestResult`` describing the observed runs, exact and
    ///   approximate p-values, and randomness decision, or `nil` if there are
    ///   fewer than two usable observations after tie handling.
    /// - Throws: ``SSError`` or distribution errors if the normal CDF or
    ///   quantile cannot be evaluated.
    public static func runsTest(data: SSExamine<T, T>, cuttingPoint cp: RunsTestCuttingPoint<T>, alpha: T, withContinuityCorrection coco: Bool = false) throws -> RunsTestResult<T>? {
        guard data.sampleSize >= 2 else {
            return nil
        }
        let items = data.itemsAsArray(sorted: .original)
        var cuttingPoint: T = .nan
        var signs: [Int] = []
        switch cp {
        case .median:
            if let cp = data.median {
                cuttingPoint = cp
            }
            else {
                SSLog.statisticsError("Median in runstest not available")
                return nil
            }
        case .mean:
            if let cp = data.arithmeticMean {
                cuttingPoint = cp
            }
            else {
                SSLog.statisticsError("Mean in runstest not available")
                return nil
            }
        case .mode:
            if let a = data.mode {
                if let cp = a.sorted(by: { $0 < $1 }).first {
                    cuttingPoint = cp
                }
                else {
                    SSLog.statisticsError("Mode in runstest not available")
                    return nil
                }
            }
            else {
                SSLog.statisticsError("Mode in runstest not available")
                return nil
            }
        case .userDefined(cuttingPoint: let cp):
            cuttingPoint = cp
        }
        for e in items {
            if e < cuttingPoint {
                signs.append(-1)
            }
            else if e > cuttingPoint {
                signs.append(1)
            }
        }
        let N = signs.count
        if N < 2 {
            SSLog.statisticsError("All values are equal or too few after tie handling.")
            return nil
        }
        let n1 = signs.filter({$0 == 1}).count
        let n2 = signs.filter({$0 == -1}).count
        var Runs = 1
        for i in 1...N-1 {
            if signs[i] != signs[i - 1] {
                Runs += 1
            }
        }
        let n1f: T = T(n1)
        let n2f: T = T(n2)
        let Nf: T = T(N)
        let Runsf: T = T(Runs)
        let meanR = (T.two * n1f * n2f) / Nf + T.one
        let varR:T = (T.two * n1f * n2f * (T.two * n1f * n2f - n1f - n2f)) / ((Nf * Nf) * (Nf - T.one))
        if varR.isZero {
            SSLog.statisticsError("Var(R) is zero")
            return nil
        }
        let sdR = T.sqrt(varR)
        let diff = Runsf - meanR
        var z: T
        if coco {
            z = (diff.magnitude - .half) / sdR
            z = diff.sign == .minus ? -z : z
        }
        else {
            z = diff / meanR
        }
        let pTwoApprox: T = try { () -> T in
            let p1 = try Distribution.Normal<T>().cdf(z.magnitude)
            return T.two * (T.one - p1)
        }()
        let cv: T = try Distribution.Normal<T>().quantile(T.one - alpha / T.two)
        let pExact = try runsPExact(r: Runs, n1: n1, n2: n2)
        let res: RunsTestResult<T> = .init(continuityCorrection: coco, alpha: alpha, cp: cuttingPoint, nAboveCP: n1, nBelowcp: n2, nRuns: Runs, meanRuns: meanR, varianceRuns: varR, zStatApprox: z, criticalValue: cv, pValueApprox: pTwoApprox, pTwoSidedExact: pExact.pTwoSidedExact, pTwoSidedExactSymmetric: pExact.pTwoSidedSym, pUpperExact: pExact.pUpper, pLowerExact: pExact.pLower, signs: signs, randomness: pTwoApprox >= alpha)
        return res
    }
    
    /// Convenience overload of ``runsTest(data:cuttingPoint:alpha:withContinuityCorrection:)`` for raw arrays.
    ///
    /// This initialises an ``SSExamine`` view of `data` and forwards to the main implementation.
    ///
    /// - Parameters:
    ///   - data: Sample to test for randomness in ordering.
    ///   - cp: Cutting point definition used to transform values into two symbols.
    ///   - alpha: Significance level in `(0, 1)` used for the normal critical value.
    ///   - cc: Whether to apply a continuity correction in the asymptotic normal approximation.
    /// - Returns: A ``RunsTestResult`` or `nil` when the sample cannot be split into two non-empty parts.
    /// - Throws: ``SSError`` if exact probabilities cannot be computed.
    public static func runsTest(data: Array<T>, cuttingPoint cp: RunsTestCuttingPoint<T>, alpha: T, withContinuityCorrection cc: Bool = false) throws -> RunsTestResult<T>? {
        let examine: SSExamine<T, T> = .init(usingArray: data, name: nil, characterSet: nil)
        return try runsTest(data: examine, cuttingPoint: cp, alpha: alpha, withContinuityCorrection: cc)
    }
    
    /// Probability mass function for the run count given the number of
    /// observations above and below the cutting point.
    /// - Parameters:
    ///   - r: Observed number of runs.
    ///   - n1: Count of values above the cutting point.
    ///   - n2: Count of values below the cutting point.
    /// - Returns: The probability of observing `r` runs for the provided split.
    /// - Throws: ``SSError`` if binomial coefficients cannot be computed.
    internal static func runsPMF(r: Int, n1: Int, n2: Int) throws -> T {
        if n1 <= 0 || n2 <= 0 {
            return T.zero
        }
        let nMin = min(n1, n2)
        let rMin = 2
        let rMax = 2 * nMin + (n1 == n2 ? 0 : 1)
        if r < rMin || r > rMax {
            return T.zero
        }
        let denom = T(try SwiftyBoost.SpecialFunctions.binomial_coeff(UInt32(n1 + n2), UInt32(n1)))
        if denom.isZero { return T.zero }
        if r % 2 == 0 {
            let k = r / 2
            let f1 = try SpecialFunctions.binomial_coeff(UInt32(n1 - 1), UInt32(k - 1))
            let f2 = try SpecialFunctions.binomial_coeff(UInt32(n2 - 1), UInt32(k - 1))
            let num = T.two * T(f1) * T(f2)
            return num / denom
        }
        else {
            let k = (r - 1) / 2
            let f1 = T(try SpecialFunctions.binomial_coeff(UInt32(n1 - 1), UInt32(k)))
            let f2 = T(try SpecialFunctions.binomial_coeff(UInt32(n2 - 1), UInt32(k - 1)))
            let f3 = T(try SpecialFunctions.binomial_coeff(UInt32(n2 - 1), UInt32(k)))
            let f4 = T(try SpecialFunctions.binomial_coeff(UInt32(n1 - 1), UInt32(k - 1)))
            let num = f1 * f2 + f3 * f4
            return num / denom
        }
    }
    
    /// Computes exact lower, upper, symmetric, and two-sided p-values for an
    /// observed run count.
    /// - Parameters:
    ///   - r: Observed number of runs.
    ///   - n1: Count of values above the cutting point.
    ///   - n2: Count of values below the cutting point.
    /// - Returns: Tail probabilities covering lower, upper, symmetric, and
    ///   two-sided exact conventions.
    /// - Throws: ``SSError`` if binomial coefficients cannot be computed.
    internal static func runsPExact(r: Int, n1: Int, n2: Int) throws -> (pLower: T, pUpper: T, pTwoSidedSym: T, pTwoSidedExact: T) {
        let nMin = min(n1, n2)
        let rMin: Int = 2
        let rMax: Int = 2 * nMin + (n1 == n2 ? 0 : 1)
        var probs:[Int:T] = [:]
        for i in rMin...rMax {
            let p = try runsPMF(r: i, n1: n1, n2: n2)
            probs[Int(i)] = p
        }
        let pLower = (rMin...r).reduce(T.zero) { acc, r in
            acc + (probs[Int(r)] ?? T.zero)
        }
        let pUpper = (r...rMax).reduce(T.zero) { acc, r in
            acc + (probs[Int(r)] ?? T.zero)
        }
        
        let twoSidedSym = min(T.one, T.two * min(pLower, pUpper))
        let pObs = probs[Int(r)] ?? T.zero
        let pTwoSidedExact = probs.values.reduce(T.zero) { acc, p in
            acc + (p <= pObs + T.ulpOfOne ? p : T.zero)
        }
        return (pLower, pUpper, twoSidedSym, pTwoSidedExact)
    }
    
    /// Holds the results of the runs test
    public struct RunsTestResult<FP: RealLike>: CustomStringConvertible, Codable {
        /// Whether a continuity correction was applied in the asymptotic calculation.
        public var continuityCorrection: Bool?
        /// Significance level used by the test.
        public var alpha: FP?
        // cutting point used
        /// Cutting point used to split the series into binary signs.
        public var cp: FP?
        /// Number of items > cutting point
        public var nAboveCP: Int?
        /// Number of items < cutting point
        public var nBelowcp: Int?
        /// Number of runs
        public var nRuns: Int?
        /// Expected number of runs under the null.
        public var meanRuns: FP?
        /// Variance of the number of runs under the null.
        public var varianceRuns: FP?
        /// z value
        public var zStatApprox: FP?
        /// critical value
        public var criticalValue: FP?
        /// two-sided p-value asymptotic
        public var pValueApprox: FP?
        /// exact two-sided p-value
        public var pTwoSidedExact: FP?
        /// symmetric two-sided p-value
        public var pTwoSidedExactSymmetric: FP?
        /// exact upper p-value
        public var pUpperExact: FP?
        /// exact lower p-value
        public var pLowerExact: FP?
        /// Array of signs
        public var signs: Array<Int>?
        /// Randomness?
        public var randomness: Bool?
        /// Returns a description
        public var description: String {
            get {
                var descr = String()
                if  let cc = self.continuityCorrection,
                    let a = self.alpha,
                    let ccp = self.cp,
                    let ng = self.nAboveCP,
                    let nl = self.nBelowcp,
                    let runs = self.nRuns,
                    let mR = self.meanRuns,
                    let vR = self.varianceRuns,
                    let z = self.zStatApprox,
                    let cv = self.criticalValue,
                    let pa = self.pValueApprox,
                    let pe2 = self.pTwoSidedExact,
                    let pe2s = self.pTwoSidedExactSymmetric,
                    let peu = self.pUpperExact,
                    let pel = self.pLowerExact,
                    let s = self.signs,
                    let rnd = self.randomness {
                    descr.append("*******************\n")
                    descr.append("RUNS TEST (WALD-WOLFOWITZ)\n")
                    descr.append("*******************\n")
                    descr.append("Use continuity correction: \(cc)\n")
                    descr.append("alpha: \(niceNumber(a))\n")
                    descr.append("cutting point: \(niceNumber(ccp))\n")
                    descr.append("number of values > cutting point: \(ng)\n")
                    descr.append("number of values < cutting point: \(nl)\n")
                    descr.append("number of runs: \(runs)\n")
                    descr.append("mean of runs: \(niceNumber(mR))\n")
                    descr.append("variance of runs: \(niceNumber(vR))\n")
                    descr.append("z: \(niceNumber(z))\n")
                    descr.append("critical value: \(niceNumber(cv))\n")
                    descr.append("asymp. p-value (two-sided): \(niceNumber(pa))\n")
                    descr.append("exact p-value (two-sided): \(niceNumber(pe2))\n")
                    descr.append("exact p-value (two-sided symmetric): \(niceNumber(pe2s))\n")
                    descr.append("exact lower p-value: \(niceNumber(pel))\n")
                    descr.append("exact upper p-value: \(niceNumber(peu))\n")
                    descr.append("signs: \(s)\n")
                    descr.append("randomness?: \(rnd)\n")
                }
                return descr
            }
        }
    }
    
    /// Cutting-point strategies for deriving binary sequences in the runs test.
    public enum RunsTestCuttingPoint<FP: RealLike>: CustomStringConvertible, Codable {
        /// Human-readable label for the selected cutting point.
        public var description: String {
            switch self {
            case .median: return "median"
            case .mean: return "mean"
            case .mode: return "mode"
            case .userDefined(let cuttingPoint):return "\(niceNumber(cuttingPoint))"
            }
        }
        case median
        case mean
        case mode
        case userDefined(cuttingPoint: FP)
    }
}
