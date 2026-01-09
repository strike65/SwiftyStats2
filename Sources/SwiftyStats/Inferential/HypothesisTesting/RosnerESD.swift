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
//
//  Rosner's ESD (Extreme Studentized Deviate) test for outliers
//

/// Rosner's generalised ESD (Extreme Studentized Deviate) outlier test.
///
/// Implements the iterative procedure from Rosner (1983) to detect up to `k`
/// outliers in approximately normally distributed data by comparing `R_i`
/// statistics against critical values `lambda_i`.

//  Overview
//  - Implements Rosner's generalized ESD procedure to detect up to k outliers in
//    approximately normally distributed data.
//  - Supports both two-sided and one-sided variants via ESDTestType.
//  - Uses Rosner's backward decision rule: after iteratively removing the most
//    extreme candidate(s), select the largest i such that R_i > λ_i.
//
//  Assumptions
//  - The underlying distribution is approximately normal after removing up to
//    the true outliers.
//  - Observations are independent.
//  - The test is sensitive to non-normality; consider transformations when needed.
//
//  Key Notation
//  - n: original sample size (constant across runs).
//  - i: run index (1-based in Rosner's description).
//  - R_i: test statistic at run i, defined as max deviation from the current mean
//         divided by the current sample standard deviation.
//  - λ_i: critical value computed from a t-quantile with adjusted p and degrees
//         of freedom df = n - i - 1.
//  - Tail handling:
//      * bothTails: uses a two-sided p-adjustment (alpha is halved).
//      * upperTail: only positive deviations (x - mean).
//      * lowerTail: only negative deviations (mean - x).
//
//  References
//  - B. Rosner (1983). Percentage Points for a Generalized ESD Many-Outlier
//    Procedure. Technometrics, 25(2), 165–172.
//  - F. E. Grubbs (1950, 1969) for single-outlier tests (related but not identical).
//
//  Notes on implementation
//  - The original sample size n is used in all λ_i computations (as in Rosner).
//  - The current mean and standard deviation are recomputed after each removal.
//  - Degrees of freedom are validated to avoid invalid t-quantile requests.
//  - For one-sided tests, if all deviations are non-positive w.r.t. the chosen
//    tail, the procedure stops early.
//
//  Error handling
//  - Invalid alpha bounds or sample sizes throw SSError(.invalidArgument).
//  - When degrees of freedom become non-positive for required i, an error is
//    thrown with guidance on maxOutliers ≤ n − 2.
//
//  Performance
//  - Sorting once: O(n log n).
//  - Up to k iterations of scanning and recomputation. Each recomputation uses
//    SSExamine to obtain mean and sample SD over the current dataset.
//  - Typical overall complexity: O(n log n + k n).
//
//  Example
//  let data: [Double] = [1, 2, 2, 3, 100, 3, 2, 1]
//  let result = try Outliers<Double>.rosnerESD(
//      data: data,
//      alpha: 0.05,
//      maxOutliers: 2,
//      testType: .bothTails
//  )
//  // Inspect result?.outliers, result?.countOfOutliers, etc.
//
/// Rosner's generalized ESD multi-outlier detection routines.

import SwiftyStatsPrelude
import SwiftyBoost

extension Inferential.HypothesisTesting.OutliersTests {
    /// Computes Rosner's adjusted p for the t-quantile at run i.
    ///
    /// This utility implements the p-adjustment specified by Rosner for the
    /// generalized ESD procedure. The adjustment depends on the tail selection:
    /// - bothTails: p = 1 − α / [2 (n − i + 1)]
    /// - one-sided (lowerTail or upperTail): p = 1 − α / (n − i + 1)
    ///
    /// - Parameters:
    ///   - alpha: Significance level in (0, 1).
    ///   - sampleSize: Original sample size n (held constant across runs).
    ///   - i: Run index (1-based as per Rosner). Valid range is 1...k.
    ///   - testType: Tail configuration. For `.bothTails`, alpha is halved.
    /// - Returns: The adjusted p in (0, 1) for use with a Student-t quantile.
    fileprivate static func rosnerP(alpha: T, sampleSize: Int, run i: Int, testType: ESDTestType) -> T {
        let n = T(sampleSize)
        let ii = T(i)
        // ex1 = n - i
        let ex1 = n - ii
        // ex2 = ex1 + 1 = n - i + 1
        let ex2 = ex1 + T.one
        switch testType {
        case .bothTails:
            // p = 1 - alpha / [2 (n - i + 1)]
            let denom = T(2) * ex2
            return T.one - alpha / denom
        case .lowerTail, .upperTail:
            // p = 1 - alpha / (n - i + 1)
            return T.one - alpha / ex2
        }
    }
    
    /// Computes Rosner's critical value λ_i for run i.
    ///
    /// λ_i is derived from a t-quantile with:
    /// - Adjusted p: see `rosnerP(alpha:sampleSize:run:testType:)`.
    /// - Degrees of freedom: df = n − i − 1 (must be > 0).
    ///
    /// Formula:
    /// - Let t = t_{df}^{−1}(p). Then
    ///   λ_i = (n − i) * t / sqrt( (n − i − 1 + t^2) * (n − i + 1) )
    ///
    /// Preconditions:
    /// - df = n − i − 1 > 0. If violated, an error is thrown and logged.
    ///
    /// - Parameters:
    ///   - alpha: Significance level in (0, 1).
    ///   - sampleSize: Original sample size n (constant across runs).
    ///   - i: Run index (1-based).
    ///   - testType: Tail configuration.
    /// - Returns: The critical value λ_i.
    /// - Throws: Rethrows any error from the Student-t quantile computation,
    ///           or `SSError(.invalidArgument)` if df <= 0.
    fileprivate static func rosnerLambdaRun(alpha: T, sampleSize: Int, run i: Int, testType: ESDTestType) throws -> T {
        let p: T = rosnerP(alpha: alpha, sampleSize: sampleSize, run: i, testType: testType)
        // df = n - i - 1 must be > 0
        let df = T(sampleSize - i - 1)
        guard df > 0 else {
            let msg = "Invalid degrees of freedom for Rosner ESD: df = n - i - 1 = \(sampleSize - i - 1). Ensure maxOutliers <= n - 2."
            SSLog.statisticsError(msg)
            throw SSError.init(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: msg)
        }
        let tQuantile = try SwiftyBoost.Distribution.StudentT(degreesOfFreedom: df).quantile(p)
        // λ_i = (n - i) * t / sqrt((n - i - 1 + t^2) * (n - i + 1))
        let ni = T(sampleSize - i)
        let num = ni * tQuantile
        let ex1 = ni - T.one                  // n - i - 1
        let ex2 = ex1 + T.pow(tQuantile, 2)  // n - i - 1 + t^2
        let ex3 = df + T(2)                  // n - i + 1
        let denom = (ex2 * ex3).squareRoot()
        return num / denom
    }
    
    /// Rosner's ESD test for up to `maxOutliers` potential outliers.
    ///
    /// This function implements the generalized ESD algorithm with support for both
    /// two-sided and one-sided alternatives. The procedure:
    /// 1. Sorts the data and initializes an internal SSExamine structure.
    /// 2. For k = 1...maxOutliers:
    ///    - Recompute the current mean and sample SD from the remaining data.
    ///    - Identify the most extreme candidate according to `testType`:
    ///      • bothTails: maximize |x − mean|
    ///      • upperTail: maximize (x − mean)
    ///      • lowerTail: maximize (mean − x)
    ///    - Compute the test statistic R_k = maxDeviation / currentSD (0 if SD == 0).
    ///    - Compute λ_k via `rosnerLambdaRun(alpha:sampleSize:run:testType:)` using
    ///      the original n (Rosner's prescription).
    ///    - Remove the identified item and repeat.
    /// 3. Backward decision rule:
    ///    - Let i* be the largest index such that R_i > λ_i. Then i* is the number
    ///      of outliers detected. The outliers are the first i* removed items.
    ///
    /// Preconditions:
    /// - data.count ≥ 3
    /// - alpha ∈ (0, 1)
    /// - 1 ≤ maxOutliers ≤ n − 2 (to ensure df ≥ 1 for the final required run)
    ///
    /// Tail behavior and early termination:
    /// - For one-sided tests, if all deviations are ≤ 0 with respect to the chosen
    ///   tail at some iteration, the procedure stops early (no further candidates).
    ///
    /// Return semantics:
    /// - Returns `nil` only if `data` is empty.
    /// - Otherwise returns an `ESDTestResult` populated with per-run statistics:
    ///   • alpha, testType, maxOutliers
    ///   • arrays of means, standard deviations, test statistics (R_k), and λ_k
    ///   • removed items in order and the final set of detected outliers
    ///   • countOfOutliers (the selected i*)
    ///
    /// Errors:
    /// - Throws `SSError(.invalidArgument)` for invalid alpha, insufficient n,
    ///   or invalid maxOutliers bounds.
    /// - Rethrows any error that occurs while evaluating t-quantiles for λ_k.
    ///
    /// Numerical considerations:
    /// - If the current sample SD becomes 0.0 at some run, R_k is set to 0.0,
    ///   which typically prevents further outlier detection.
    ///
    /// - Parameters:
    ///   - data: Sample values. At least 3 observations are required.
    ///   - alpha: Significance level in (0, 1).
    ///   - maxOutliers: Upper bound on outliers to test (recommended ≤ n − 2).
    ///   - testType: Tail configuration (.bothTails, .upperTail, .lowerTail).
    /// - Returns: `ESDTestResult<T>` describing the run, or `nil` for empty input.
    /// - Throws: `SSError` for invalid inputs, or distribution errors from quantile evaluation.
    public static func rosnerESD(data: [T], alpha: T, maxOutliers: Int, testType: ESDTestType) throws -> ESDTestResult<T>? {
        guard !data.isEmpty else { return nil }
        // alpha bounds
        guard alpha > 0 && alpha < 1 else {
            let msg = "Alpha must be in (0, 1). Provided alpha = \(alpha)."
            SSLog.statisticsError(msg)
            throw SSError.init(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: msg)
        }
        // Ensure at least 3 observations and maxOutliers ≤ n - 2 so that df = n - i - 1 ≥ 1
        guard data.count >= 3 else {
            let msg = "Rosner ESD requires at least 3 observations. Provided n = \(data.count)."
            SSLog.statisticsError(msg)
            throw SSError.init(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: msg)
        }
        guard maxOutliers <= data.count - 2 && maxOutliers >= 1 else {
            let msg = "maxOutliers must satisfy 1 <= maxOutliers <= n - 2. Provided n = \(data.count), maxOutliers = \(maxOutliers)."
            SSLog.statisticsError(msg)
            throw SSError.init(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: msg)
        }
        
        var sorted = data.sorted(by: <)
        let examine: SSExamine<T,T>
        do {
            examine = try SSExamine<T, T>.init(using: sorted, levelOfMeasurement: .ratio, name: nil, characterSet: nil)
        }
        catch let error {
            throw error
        }
        // Store original n; Rosner uses original n in λ_i at all runs
        let originalN = examine.sampleSize
        
        let sd = examine.sampleStandardDeviation!
        let orgMean = examine.arithmeticMean!
        var currentSd: T = sd
        var currentMean: T = orgMean
        
        var itemsRemoved = Array<T>()
        var means = Array<T>()
        var stdDevs = Array<T>()
        var testStats = Array<T>()
        var lambdas = Array<T>()
        var outliers = Array<T>()
        
        let currentData: SSExamine<T, T> = SSExamine<T, T>()
        let _ = currentData.append(contentOf: sorted)
        
        var k = 1
        while k <= maxOutliers {
            // Reset per-run sentinels
            var maxDev: T = -T.infinity
            var currentIndex: Int = -1
            
            // Compute mean from current data
            guard let runMean = currentData.arithmeticMean else { break }
            currentMean = runMean
            
            // Scan for most extreme candidate according to the selected tail
            for (idx, t) in sorted.enumerated() {
                let deviation: T
                switch testType {
                case .bothTails:
                    deviation = abs(t - currentMean)
                case .upperTail:
                    deviation = t - currentMean
                case .lowerTail:
                    deviation = currentMean - t
                }
                if deviation > maxDev {
                    maxDev = deviation
                    currentIndex = idx
                }
            }
            
            // If no valid candidate found (e.g., one-sided and all deviations ≤ 0), stop early
            if currentIndex < 0 || maxDev <= 0 {
                // Append nothing further; break the outer loop
                break
            }
            
            // Compute current SD and R_k
            guard let runSd = currentData.sampleStandardDeviation else { break }
            currentSd = runSd
            let currentTestStat: T = currentSd == 0 ? T(0) : (maxDev / currentSd)
            
            // Compute λ_k with correct tail p
            let currentLambda = try rosnerLambdaRun(alpha: alpha, sampleSize: originalN, run: k, testType: testType)
            
            // Record run statistics before removal
            testStats.append(currentTestStat)
            lambdas.append(currentLambda)
            means.append(currentMean)
            stdDevs.append(currentSd)
            
            // Remove the identified element and update current dataset
            let itemToRemove = sorted[currentIndex]
            itemsRemoved.append(itemToRemove)
            sorted.remove(at: currentIndex)
            currentData.removeAll()
            let _ = currentData.append(contentOf: sorted)
            
            k += 1
        }
        
        // Backward decision: find the largest i such that R_i > λ_i
        var countOfOL = 0
        if !testStats.isEmpty {
            for i in stride(from: testStats.count - 1, through: 0, by: -1) {
                if testStats[i] > lambdas[i] {
                    countOfOL = i + 1
                    break
                }
            }
        }
        
        // Outliers are the first countOfOL removed items
        if countOfOL > 0 {
            outliers.append(contentsOf: itemsRemoved.prefix(countOfOL))
        }
        
        var res: ESDTestResult<T> = ESDTestResult<T>()
        res.alpha = alpha
        res.countOfOutliers = countOfOL
        res.itemsRemoved = itemsRemoved
        res.lambdas = lambdas
        res.maxOutliers = maxOutliers
        res.means = means
        res.outliers = outliers
        res.stdDeviations = stdDevs
        res.testStatistics = testStats
        res.testType = testType
        return res
    }
    
    /// Result container for the Rosner ESD iterative outlier test.
    public struct ESDTestResult<FP: RealLike>: CustomStringConvertible  {
        /// Standard deviation snapshots for each iteration.
        public var stdDeviations: Array<FP>?
        /// Items removed at each iteration in the order they were flagged.
        public var itemsRemoved: Array<FP>?
        /// Test statistics R_i computed per iteration.
        public var testStatistics: Array<FP>?
        /// Critical values λ_i matched to each statistic.
        public var lambdas: Array<FP>?
        /// Final count of detected outliers.
        public var countOfOutliers: Int?
        /// Values identified as outliers after applying Rosner’s decision rule.
        public var outliers: Array<FP>?
        /// Significance level α used for the run.
        public var alpha: FP?
        /// Maximum number of outliers k that were tested.
        public var maxOutliers: Int?
        /// Tail configuration (two-sided, upper, or lower).
        public var testType: ESDTestType?
        /// Mean snapshots for each iteration.
        public var means: Array<FP>?
        /// Pretty-printed multiline description of the iterative run.
        public var description: String {
            get {
                var descr = String()
                descr.append("ESD TEST FOR OUTLIERS\n")
                descr.append("*********************\n")
                if let sd = self.stdDeviations, let ir = self.itemsRemoved, let s = self.testStatistics, let l = self.lambdas, let no = self.countOfOutliers, let o = self.outliers, let a = self.alpha, let tt = self.testType, let me = self.means {
                    descr.append("test type: \(tt)\n")
                    descr.append("alpha: \(niceNumber(a))\n")
                    descr.append("number of detected outliers: \(no)\n")
                    descr.append("outliers: [\(o.map { niceNumber($0) }.joined(separator: ", "))]\n")
                    descr.append("statistics: [\(s.map { niceNumber($0) }.joined(separator: ", "))]\n")
                    descr.append("removed items: [\(ir.map { niceNumber($0) }.joined(separator: ", "))]\n")
                    descr.append("lambdas: [\(l.map { niceNumber($0) }.joined(separator: ", "))]\n")
                    descr.append("means: [\(me.map { niceNumber($0) }.joined(separator: ", "))]\n")
                    descr.append("sd: [\(sd.map { niceNumber($0) }.joined(separator: ", "))]\n")
                }
                return descr
            }
        }
    }
}
