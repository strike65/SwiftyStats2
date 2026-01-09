//
//  Created by VT on 07.12.25.
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

/// Kolmogorov-Smirnov two-sample test against paired empirical distributions.
extension Inferential.HypothesisTesting.NonParametricTests {
    
    /// Runs a two-sample Kolmogorov-Smirnov test between two empirical datasets.
    ///
    /// ### Note ###
    /// H<sub>0</sub>: The two samples are from populations with same distribution function (F<sub>1</sub>(x) = F<sub>2</sub>(x))
    ///
    /// H<sub>a1</sub>:(F<sub>1</sub>(x) > F<sub>2</sub>(x))
    ///
    /// H<sub>a2</sub>:(F<sub>1</sub>(x) < F<sub>2</sub>(x))
    ///
    /// The routine walks both empirical CDFs, tracks the largest positive and negative
    /// deviations, and derives the two-sided p-value from the Kolmogorov distribution with
    /// an effective sample size of `n1 * n2 / (n1 + n2)`.
    ///
    /// - Parameters:
    ///   - set1: First empirical dataset with precomputed empirical CDF support.
    ///   - set2: Second empirical dataset to compare against `set1`.
    /// - Returns: KS summary containing `D+`, `D-`, absolute `D`, normalised `z`, and the two-sided p-value.
    /// - Throws: `SSError.invalidArgument` when either dataset is empty, smaller than three
    ///   observations, or so uniform that the KS statistic cannot be computed.
    public static func ksTest(set1: SSExamine<T, T>, set2: SSExamine<T, T>) throws -> KSTwoSampleTestResult<T>? {
        guard !set1.isEmpty, !set2.isEmpty else {
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "One or both sample sets are empty.")
        }
        guard set1.sampleSize > 2, set2.sampleSize > 2 else {
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Sample sizes must be greater than 2.")
        }
        
        let data = set1.itemsAsArray(sorted: .ascending) + set2.itemsAsArray(sorted: .ascending)
        let n1: T = T(set1.sampleSize)
        let n2: T = T(set2.sampleSize)
        var dcdf: T
        var maxNegDev: T = 0
        var maxPosDev: T = 0
        if n1 > n2 {
            for element in data {
                if let cdf1 = set1.empiricalCDF(item: element), let cdf2 = set2.empiricalCDF(item: element) {
                    dcdf = cdf1 - cdf2
                    if dcdf < 0 {
                        maxNegDev = dcdf < maxNegDev ? dcdf : maxNegDev
                    }
                    else {
                        maxPosDev = dcdf > maxPosDev ? dcdf : maxPosDev
                    }
                }
            }
        }
        else {
            for element in data {
                if let cdf1 = set1.empiricalCDF(item: element), let cdf2 = set2.empiricalCDF(item: element) {
                    dcdf = cdf2 - cdf1
                    if dcdf < 0 {
                        maxNegDev = dcdf < maxNegDev ? dcdf : maxNegDev
                    }
                    else {
                        maxPosDev = dcdf > maxPosDev ? dcdf : maxPosDev
                    }
                }
            }
        }
        let maxD: T = maxNegDev.magnitude > maxPosDev.magnitude ? maxNegDev.magnitude : maxPosDev.magnitude
        var z: T
        guard !maxD.isNaN else {
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Data too uniform to calculate KS test")
        }
        z = maxD * T.sqrt((n1 * n2) / (n1 + n2))
        // Approximated p-Values Smirnov 1948
        // archived for historical reasons
        //        var p: T = .zero, q: T = .zero, term1: T
        //        if ((z >= .zero) && (z < T(0.27))) {
        //            p = .one
        //        }
        //        else if ((z >= T(0.27 )) && (z < .one)) {
        //            term1 = T(-1.233701) * T.pow(z, 2)
        //            q = T.exp(term1)
        //            term1 = (T.sqrt(T.twoPi) * (q + T.pow(q, 25) + T.pow(q, 9))) / z
        //            p = .one - term1
        //        }
        //        else if ((z >= 1) && z < T(3.1)) {
        //            q = T.exp(-.two * T.pow(z, 2))
        //            term1 = T.pow(q, 4) + T.pow(q, 9) - T.pow(q, 16)
        //            p = .two * (q - term1)
        //        }
        //        else if z >= T(3.1) {
        //            p = .zero
        //        }
        let Neff = (n1 * n2) / (n1 + n2)
        let ksDist: Distribution.KolmogorovSmirnov<T> = try .init(numberOfObservations: Neff)
        var p = try ksDist.cdf(maxD)
        p = T.one - p
        return KSTwoSampleTestResult(dMaxPos: maxPosDev, dMaxNeg: maxNegDev, dMaxAbs: maxD, zStatistic: z, p2Value: p, sampleSize1: set1.sampleSize, sampleSize2: set2.sampleSize)
    }
    
    /// Convenience overload that performs the KS test on raw array inputs.
    ///
    /// The arrays are wrapped into `SSExamine` instances before delegating to the primary
    /// KS implementation so both overloads share validation and statistic computation.
    ///
    /// - Parameters:
    ///   - set1: First collection of observations.
    ///   - set2: Second collection of observations.
    /// - Returns: KS summary containing `D+`, `D-`, absolute `D`, normalised `z`, and the two-sided p-value.
    /// - Throws: `SSError.invalidArgument` when either dataset is empty, smaller than three
    ///   observations, or so uniform that the KS statistic cannot be computed.
    public static func ksTest(set1: [T], set2: [T]) throws -> KSTwoSampleTestResult<T>? {
        let ex1: SSExamine<T,T> = .init(usingArray: set1, name: nil, characterSet: nil)
        let ex2: SSExamine<T,T> = .init(usingArray: set2, name: nil, characterSet: nil)
        return try ksTest(set1: ex1, set2: ex2)
    }
    
    /// Results of the KS-2-Sample test
    public struct KSTwoSampleTestResult<FP: RealLike>: CustomStringConvertible, Codable {
        /// max pos diff
        var dMaxPos: FP?
        /// max neg diff
        var dMaxNeg: FP?
        /// max abs diff
        var dMaxAbs: FP?
        /// z value
        var zStatistic: FP?
        /// p valie
        var p2Value: FP?
        /// size of sample 1
        var sampleSize1: Int?
        /// size of sample 2
        var sampleSize2: Int?
        /// Returns a description
        public var description: String {
            get {
                var descr = String()
                if let dmp = self.dMaxPos, let dmn = self.dMaxNeg, let dma = self.dMaxAbs, let z = self.zStatistic, let p2 = self.p2Value, let n1 = self.sampleSize1, let n2 = self.sampleSize2 {
                    descr.append("KOLMOGOROV-SMIRNOV TWO SAMPLE TEST\n")
                    descr.append("**********************************\n")
                    descr.append("two sided p-value: \(niceNumber(p2))\n")
                    descr.append("z: \(niceNumber(z))\n")
                    descr.append("n1: \(n1)\n")
                    descr.append("n1: \(n2)\n")
                    descr.append("D(-)max: \(niceNumber(dmn))\n")
                    descr.append("D(+)max: \(niceNumber(dmp))\n")
                    descr.append("|D|max: \(niceNumber(dma))\n")
                }
                return descr
            }
        }
    }
}
