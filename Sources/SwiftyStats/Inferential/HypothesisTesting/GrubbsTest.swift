//
//  Created by VT on 04.11.25.
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
// Grubbs' test utilities for detecting single outliers in normal samples.

import SwiftyStatsPrelude


extension Inferential.HypothesisTesting.OutliersTests {
    
    /// Performs Grubbs' test for a single outlier in a univariate sample.
    ///
    /// Grubbs' test evaluates whether the most extreme observation in the sample
    /// (i.e., the one with the largest absolute deviation from the mean) is an outlier
    /// under the assumption of normality. The test statistic is:
    ///
    /// G = max(|xᵢ − x̄|) / s
    ///
    /// where x̄ is the sample mean and s is the (sample) standard deviation.
    /// The critical value is computed using the Student’s t distribution with
    /// n − 2 degrees of freedom at tail probability α / (2n).
    ///
    /// - Parameters:
    ///   - v: The sample values to be tested. Must contain at least 3 observations.
    ///   - mean: The mean of the sample `v` (x̄). Must correspond to `v`.
    ///   - sd: The (sample) standard deviation of the sample `v` (s). Must correspond to `v` and be positive.
    ///   - a: The significance level α. Must lie in the open interval (0, 1).
    ///
    /// - Returns: A `SSGrubbsTestResult` describing the test, or `nil` if the computation of the required
    ///   t-quantile fails. The `hasOutliers` flag indicates whether the most extreme observation is detected
    ///   as an outlier (i.e., `G > criticalValue`).
    ///
    /// - Throws: `SSError` of type `.missingData` if there are fewer than 3 observations.
    ///           `SSError` of type `.invalidArgument` if `alpha` is not in (0, 1) or if `standardDeviation` is not positive/finite.
    ///   Errors are also logged via `SSLog.statisticsError`.
    ///
    /// - Important: `mean` and `standardDeviation` must be computed from exactly the same `values` supplied to this function.
    ///   Supplying inconsistent summary statistics will invalidate the test.
    ///
    /// - Precondition: The input sample is approximately normally distributed; Grubbs' test assumes normality.
    ///
    /// - Note: If multiple outliers are suspected, remove one detected outlier and re-run iteratively (with caution).
    ///
    /// - SeeAlso: `SSGrubbsTestResult`, `SwiftyBoost.Distribution.StudentT`
    ///
    /// References:
    /// - Grubbs, F. E. (1950). Sample criteria for testing outlying observations. The Annals of Mathematical Statistics, 21(1), 27–58.
    /// - ASTM E178 – Standard Practice for Dealing With Outlying Observations.
    /// - NIST/SEMATECH e-Handbook of Statistical Methods: 1.3.5.17. Grubbs' Test for Outliers.
    public static func grubbsTest(values v: [T], mean: T, standardDeviation sd: T, alpha a: T) throws -> GrubbsTestResult<T>? {
        // Need at least three observations
        guard v.count >= 3 else {
            let msg = "Not enough data to perform Grubbs' test - at least 3 elements are needed"
            SSLog.statisticsError(msg)
            throw SSError(type: .missingData, file: #fileID, line: #line, function: #function, reason: msg)
        }
        // Validate alpha in (0,1)
        guard a.liesInOpenRange(from: 0, to: 1) else {
            let msg = "alpha must be in the open interval (0, 1)"
            SSLog.statisticsError(msg)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: msg)
        }
        // Validate standard deviation
        guard sd.isFinite && sd > 0 else {
            let msg = "standardDeviation must be positive and finite"
            SSLog.statisticsError(msg)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: msg)
        }
        
        // Compute G statistic
        let n: T = T(v.count)
        // min/max are defined for non-empty arrays; v.count >= 3 here
        guard let mi:T = v.min(), let ma:T = v.max() else {
            // Defensive fallback; should not happen with n >= 3
            return nil
        }
        let maxDiff: T = max((ma - mean).magnitude, (mi - mean).magnitude)
        let G: T = maxDiff / sd
        
        // Get t-quantile for df = n-2 at upper tail α/(2n)
        let tQuantile: T
        do {
            // Upper-tail critical value t_{ν, 1 - α/(2n)}
            let upperTailP: T = 1 - a / (T(2) * n)
            tQuantile = T(try SwiftyBoost.Distribution.StudentT(degreesOfFreedom: n - 2).quantile(upperTailP))
        } catch {
            return nil
        }
        // Critical value for Grubbs' test
        let t2: T = T.pow(tQuantile, 2)
        let t3: T = (n - 1) / n.squareRoot()
        let t4: T = n - 2 + T(t2)
        let t5: T = T(t2) / t4
        let tCrit: T = t3 * t5.squareRoot()
        
        var res: GrubbsTestResult<T> = GrubbsTestResult<T>()
        res.sampleSize = Int(n)
        res.maxDiff = maxDiff
        res.largest = ma
        res.smallest = mi
        res.criticalValue = tCrit
        res.mean = mean
        res.G = G
        res.stdDev = sd
        res.hasOutliers = G > tCrit
        return res
    }
    
    
    /// Result container for the single-outlier Grubbs test.
    public struct GrubbsTestResult<F: RealLike>: CustomStringConvertible {
        /// Critical value G_crit at the requested alpha level.
        public var criticalValue: F?
        /// Largest observation in the sample.
        public var largest: F?
        /// Smallest observation in the sample.
        public var smallest: F?
        /// Number of observations used in the test.
        public var sampleSize: Int?
        /// Maximum absolute deviation |x_i − mean|.
        public var maxDiff: F?
        /// Sample mean of the data.
        public var mean: F?
        /// Observed Grubbs statistic G = maxDiff / stdDev.
        public var G: F?
        /// Sample standard deviation.
        public var stdDev: F?
        /// Indicates whether the most extreme point is flagged as an outlier (G > G_crit).
        public var hasOutliers: Bool?
        /// Pretty-printed multiline description of the result.
        public var description: String {
            get {
                var descr = String()
                descr.append("GRUBBS TEST FOR OUTLIERS\n")
                descr.append("***********************\n")
                if let cv = self.criticalValue, let l = self.largest, let s = self.smallest, let n = self.sampleSize, let md = self.maxDiff, let m = self.mean, let g = self.G, let sd = self.stdDev, let ho = self.hasOutliers {
                    descr.append("has outliers: \(ho)\n")
                    descr.append("Grubbs G: \(niceNumber(g))\n")
                    descr.append("max diff: \(niceNumber(md))\n")
                    descr.append("critical value: \(niceNumber(cv))\n")
                    descr.append("mean: \(niceNumber(m))\n")
                    descr.append("sd: \(niceNumber(sd))\n")
                    descr.append("largest value: \(niceNumber(l))\n")
                    descr.append("smallest value: \(niceNumber(s))\n")
                    descr.append("sample size: \(n)\n")
                }
                return descr
            }
        }
    }
}
