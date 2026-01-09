//
//  Created by VT on 09.12.25.
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
    
    /// Performs the Wald-Wolfowitz two-sample runs test.
    ///
    /// The two samples are pooled and sorted by value; the sequence of sample labels is then
    /// analysed for runs (maximal contiguous blocks of identical labels). Under the null
    /// hypothesis of equal distributions, the run count admits an exact distribution
    /// conditioned on `(n1, n2)` and a normal approximation for larger samples.
    ///
    /// This implementation reports both the exact two-sided p-value and the asymptotic
    /// z-based p-value, and can apply a `0.5` continuity correction to the z statistic.
    ///
    /// - Parameters:
    ///   - set1: First sample.
    ///   - set2: Second sample.
    ///   - alpha: Significance level in `(0, 1)` used to compute the two-sided critical value.
    ///   - coco: Apply a continuity correction to the asymptotic z statistic.
    /// - Returns: A ``WaldWolfowitzTwoSampleTestResult``.
    /// - Throws: ``SSError`` if inputs are invalid or probabilities cannot be computed.
    public static func waldWolfowitzTwoSampleTest(set1: SSExamine<T,T>, set2: SSExamine<T,T>, alpha: T, withContinuityCorrection coco: Bool = true) throws -> WaldWolfowitzTwoSampleTestResult<T>? {
        guard alpha > .zero, alpha < .one else {
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Alpha must be in (0,1).")
        }
        guard set1.sampleSize > 0, set2.sampleSize > 0 else {
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Both samples must contain at least one observation.")
        }
        let groups: [Int] = .init(repeating: 1, count: set1.sampleSize) + .init(repeating: 2, count: set2.sampleSize)
        let data = set1.itemsAsArray(sorted: .original) + set2.itemsAsArray(sorted: .original)
        let sorter = GroupedData<T, Int>(data: data, groups: groups)
        let (sortedData, sortedGroups) = sorter.sorted()
        let totalN = sortedGroups.count
        guard totalN >= 2 else {
            throw SSError.init(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Total sample must be >= 2.")
        }
        var runs = 1
        for i in 1..<sortedGroups.count {
            if sortedGroups[i] != sortedGroups[i - 1] {
                runs += 1
            }
        }
        
        let n1Int = sortedGroups.reduce(0) { $0 + ($1 == 1 ? 1 : 0)}
        let n2Int = totalN - n1Int
        if n1Int == 0 || n2Int == 0 {
            throw SSError.init(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Both samples must contain at least one element")
        }
        var nTiesIntergroup: Int = 0
        var nTiedCases: Int = 0
        if !sortedData.isEmpty {
            var currentValue = sortedData[0]
            var currentGroups = Set<Int>([sortedGroups[0]])
            var currentCount = 1
            for i in 1..<sortedData.count {
                if sortedData[i] == currentValue {
                    currentGroups.insert(sortedGroups[i])
                    currentCount += 1
                }
                else {
                    if currentGroups.count > 1 {
                        nTiesIntergroup += 1
                        nTiedCases += currentCount
                    }
                    currentValue = sortedData[i]
                    currentGroups = Set<Int>([sortedGroups[i]])
                    currentCount = 1
                }
            }
            if currentGroups.count > 1 {
                nTiesIntergroup += 1
                nTiedCases += currentCount
            }
        }
        let n1 = T(n1Int)
        let n2 = T(n2Int)
        let N = n1 + n2
        let numVar: T = .two * n1 * n2 * (.two * n1 * n2 - n1 - n2)
        let denVar: T = N * N * (N - .one)
        if denVar.isZero {
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Variance can not be computed for this sample.")
        }
        let varR = numVar / denVar
        if varR.isZero || !varR.isFinite {
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Variance can not be computed for this sample (variance is zero or infinite).")
        }
        let sigma: T = varR.squareRoot()
        let meanR: T = (.two * n1 * n2) / N + .one
        
        let runsF: T = T(runs)
        let diff = runsF - meanR
        let zNumerator: T = coco ? (diff.sign == .minus ? -(diff.magnitude - .half) : (diff.magnitude - .half)) : diff
        let z: T = zNumerator / sigma
        let zAbs: T = z.magnitude
        let stdN: Distribution.Normal<T> = try .init()
        let criticalValue = try stdN.quantile(.one - alpha / .two)
        let phi: T = try stdN.cdf(zAbs)
        let pAsymp: T = .two * (.one - phi)
        let rejectNull = zAbs >= criticalValue
        let exact = try Randomness<T>.runsPExact(r: runs, n1: n1Int, n2: n2Int)
        let result = WaldWolfowitzTwoSampleTestResult<T>(nRuns: runs, zStat: z, pValueExact: exact.pTwoSidedExact, pValueAsymp: pAsymp, continuityCorrection: coco, alpha: alpha, criticalValue: criticalValue, rejectNull: rejectNull, mean: meanR, variance: varR, nTiesIntergroup: nTiesIntergroup, nTiedCases: nTiedCases, sampleSize1: n1Int, sampleSize2: n2Int)
        return result
    }
    
    
    /// Convenience overload of ``waldWolfowitzTwoSampleTest(set1:set2:alpha:withContinuityCorrection:)`` for raw arrays.
    ///
    /// - Parameters:
    ///   - set1: First sample.
    ///   - set2: Second sample.
    ///   - alpha: Significance level in `(0, 1)` used to compute the two-sided critical value.
    ///   - coco: Apply a continuity correction to the asymptotic z statistic.
    /// - Returns: A ``WaldWolfowitzTwoSampleTestResult``.
    /// - Throws: ``SSError`` if inputs are invalid or probabilities cannot be computed.
    public static func waldWolfowitzTwoSampleTest(set1: [T], set2: [T], alpha: T, withContinuityCorrection coco: Bool = true) throws -> WaldWolfowitzTwoSampleTestResult<T>? {
        let ex1 = SSExamine<T,T>.init(usingArray: set1, name:  nil, characterSet: nil)
        let ex2 = SSExamine<T,T>.init(usingArray: set2, name:  nil, characterSet: nil)
        return try waldWolfowitzTwoSampleTest(set1: ex1, set2: ex2, alpha: alpha, withContinuityCorrection: coco)
    }
    
    
    /// Holds the results of the two sample runs test
    public struct WaldWolfowitzTwoSampleTestResult<FP: RealLike>: CustomStringConvertible, Codable {
        /// Number of runs
        public var nRuns: Int?
        /// z value
        public var zStat: FP?
        // critical value
        //    public var criticalValue: FPT?
        /// p-value
        public var pValueExact: FP?
        /// p-value asymptotic
        public var pValueAsymp: FP?
        /// Whether a continuity correction was applied to the asymptotic z statistic.
        public var continuityCorrection: Bool?
        /// Significance level used for the asymptotic decision.
        public var alpha: FP?
        /// Two-sided critical value from the standard normal distribution.
        public var criticalValue: FP?
        /// Whether the null hypothesis (equal distributions) is rejected using the asymptotic p-value.
        public var rejectNull: Bool?
        /// mean
        public var mean: FP?
        /// variance
        public var variance: FP?
        /// number of intergroup ties
        public var nTiesIntergroup: Int?
        /// number of inner group ties
        public var nTiedCases: Int?
        /// size of sample 1
        public var sampleSize1: Int?
        /// size of sample 2
        public var sampleSize2: Int?
        /// Returns a description
        public var description: String {
            get {
                var descr = String()
                if let nr = self.nRuns,
                   let z = self.zStat,
                   let pe = self.pValueExact,
                   let pa = self.pValueAsymp,
                   let cc = self.continuityCorrection,
                   let a = self.alpha,
                   let cv = self.criticalValue,
                   let rn = self.rejectNull,
                   let m = self.mean,
                   let v = self.variance,
                   let n1 = self.sampleSize1,
                   let n2 = self.sampleSize2,
                   let nt = self.nTiesIntergroup,
                   let tc = self.nTiedCases {
                    descr.append("WALD-WOLFOWITZ TWO SAMPLE TEST\n")
                    descr.append("******************************\n")
                    descr.append("p-value exact: \(niceNumber(pe))\n")
                    descr.append("p-value asymp: \(niceNumber(pa))\n")
                    descr.append("use continuity correction: \(cc)\n")
                    descr.append("alpha: \(niceNumber(a))\n")
                    descr.append("critical value (|z|): \(niceNumber(cv))\n")
                    descr.append("reject null?: \(rn)\n")
                    descr.append("z: \(niceNumber(z))\n")
                    descr.append("n1: \(n1)\n")
                    descr.append("n1: \(n2)\n")
                    descr.append("mean: \(niceNumber(m))\n")
                    descr.append("variance: \(niceNumber(v))\n")
                    descr.append("count of intergroup ties: \(nt)\n")
                    descr.append("count of tied cases: \(tc)\n")
                    descr.append("count runs: \(nr)\n")
                }
                return descr
            }
        }
        
    }
}
