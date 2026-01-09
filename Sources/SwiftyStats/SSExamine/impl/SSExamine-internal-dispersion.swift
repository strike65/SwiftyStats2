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

/// Internal dispersion, variance, and moment helpers for SSExamine.

import SwiftyStatsPrelude
import SwiftyBoost

// MARK: Dispersion

/// Internal, lock-requiring statistical routines for SSExamine.
///
/// Quantile operators are defined against the empirical CDF `F_n` with configurable
/// Hyndman–Fan interpolation, while deviation measures compute norms of centered
/// residuals such as `median(|x - median(x)|)`. This extension provides:
/// - Median/quantile-related operations (including Hyndman–Fan Types 1–9),
/// - Absolute deviation measures (MAD and mean absolute deviation),
/// - Two-sided confidence intervals for the mean under Normal and Student’s t assumptions.
///
/// Locking:
/// - All functions assume the caller already holds the instance lock. Public entry points
///   acquire the lock and delegate to these routines.
///
/// Return conventions:
/// - Unless specified otherwise, functions return `nil` for non-numeric datasets or empty collections.
///
/// Errors:
/// - Functions that throw use `SSError(type: .invalidArgument, ...)` for out-of-range arguments.
///
/// Complexity notes:
/// - Complexity annotations refer to the work after any cached values are available.
/// - Quantile computations require sorting a snapshot: O(n log n).
///
/// Quantiles (Hyndman–Fan Types 1–9):
/// - We support the nine standard types via an a–b–c–d parameterization (Hyndman & Fan, 1996).
/// - Unified formula on sorted values x_(1) ≤ … ≤ x_(n) (1-based order statistics):
///     r = a + (n + b) · p
///     f = fractional part of r
///     i = floor(r) clamped to [1, n], j = ceil(r) clamped to [1, n]
///     t = c + d · f
///     Q = x_(i) + (x_(j) − x_(i)) · t
/// - Type 2 (averaged ECDF) is implemented via a dedicated branch for exact midpoint behavior
///   at integer positions; other types use the unified parameterization.
extension SSExamine {
    
    /// Interquantile range |Q(lower) − Q(upper)| under lock using a selected Hyndman–Fan type.
    ///
    /// Parameters:
    /// - lower: Lower quantile probability in [0, 1].
    /// - upper: Upper quantile probability in [0, 1].
    /// - quantileType: Hyndman–Fan type to use (default `.averagedECDF_HF2`).
    ///
    /// Returns:
    /// - The absolute difference |Q(lower) − Q(upper)|, or `nil` if the dataset is empty or non-numeric.
    ///
    /// Throws:
    /// - `SSError(.invalidArgument, ...)` if `lower` or `upper` are outside [0, 1].
    ///
    /// Thread-safety:
    /// - The caller must hold the lock.
    ///
    /// Notes:
    /// - If `lower == upper`, returns `0`.
    /// - Quantiles within this routine are computed on a sorted snapshot with the specified type.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting for the underlying quantiles.
    @inline(__always)
    internal func interquantileRangeLocked(lower: FP, upper: FP, quantileType: HyndmanFanQuantileType = .type2) throws -> FP? {
        guard upper.liesInClosedRange(from: 0, to: 1), lower.liesInClosedRange(from: 0, to: 1) else {
            let message = "Quantile must be in the range [0, 1]"
            SSLog.statisticsError(message)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: message)
        }
        // Domain checks
        if !self.isNumeric || self.isEmpty {
            SSLog.statisticsError("Quantile is not defined for dataset")
            return nil
        }
        if lower == upper {
            return FP.zero
        }
        if let q1 = try self.quantileLocked(q: lower, quantileType: .type7), let q2 = try self.quantileLocked(q: upper, quantileType: .type7) {
            return (q1 - q2).magnitude
        }
        else {
            return nil
        }
    }
    
    /// Median Absolute Deviation (MAD) around a specified center under lock.
    ///
    /// Parameters:
    /// - p: The center value around which absolute deviations are computed (often the sample median).
    ///
    /// Returns:
    /// - The median of |x − p| over all numeric items, or `nil` if not applicable.
    ///
    /// Quantile semantics:
    /// - Uses Hyndman–Fan Type 2 (averaged ECDF) median on the deviation sample.
    ///
    /// Thread-safety:
    /// - The caller must hold the lock.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting of the absolute deviations.
    @inline(__always)
    internal func medianAbsoluteDeviationLocked(center p: FP) -> FP? {
        guard !self.isEmpty, self.isNumeric else {
            return self.logNil("medianAbsoluteDeviationLocked requires a non-empty numeric dataset")
        }
        guard let values = self.itemsAsNumericArrayLocked() else {
            return self.logNil("medianAbsoluteDeviationLocked could not project items to numeric values")
        }
        let deviations = values.map { ($0 - p).magnitude }.sorted()
        return SSExamine<FP, FP>.quantileHF2(sorted: deviations, p: FP(0.5))
    }
    /// Quantile under lock using a specified Hyndman–Fan convention.
    ///
    /// Parameters:
    /// - q: Probability in [0, 1]. Values outside this range throw `SSError(.invalidArgument, ...)`.
    /// - quantileType: The Hyndman–Fan type to use (HF1–HF9). Default is `.averagedECDF_HF2`.
    ///
    /// Returns:
    /// - `nil` for non-numeric or empty datasets.
    /// - Otherwise, the requested quantile as `FP`.
    ///
    /// Thread-safety:
    /// - The caller must hold the lock.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting; constant-time operations after sorting.
    @inline(__always)
    internal func quantileLocked(q: FP, quantileType: HyndmanFanQuantileType = .type2) throws -> FP? {
        guard q.liesInClosedRange(from: 0, to: 1) else {
            let message = "Quantile must be in the range [0, 1]"
            SSLog.statisticsError(message)
            throw SSError(type: .invalidArgument, file: #file, line: #line, function: #function, reason: message)
        }
        if !self.isNumeric || self.isEmpty {
            SSLog.statisticsError("Quantile is not defined for dataset")
            return nil
        }
        guard let values = self.itemsAsNumericArrayLocked() else {
            return self.logNil("quantileLocked could not project items to numeric values")
        }
        let sorted = values.sorted()
        
        switch quantileType {
        case .type1:
            return SSExamine<SSElement, FP>.quantileHF1(sorted: sorted, p: q)
        case .type2:
            return SSExamine<SSElement, FP>.quantileHF2(sorted: sorted, p: q)
        case .type3:
            return SSExamine<SSElement, FP>.quantileHF3(sorted: sorted, p: q)
        case .type4, .type5, .type6, .type7, .type8, .type9:
            let params = try self.quantileParams(quantileType: quantileType)
            return SSExamine<SSElement, FP>.quantileWithParams(p: q, sorted: sorted, params: params)
        }
    }
    
    /// Hyndman–Fan parameter preset for quantile computation.
    ///
    /// - Returns: The `(a, b)` pair used in the unified interpolation formula.
    @inline(__always)
    fileprivate func quantileParams(quantileType: HyndmanFanQuantileType) throws -> (a: FP, b: FP) {
        switch quantileType {
        case .type4:
            return (FP.zero, FP.zero)
        case .type5:
            return (FP.one / FP(2), FP.zero)
        case .type6:
            return (FP.zero, FP.one)
        case .type7:
            return (FP.one, -FP.one)
        case .type8:
            return (FP.one / FP(3), FP.one / FP(3))
        case .type9:
            return (FP(3) / FP(8), FP.one / FP(4))
        default:
            let message = "quantileParams lacks coefficients for type \(quantileType)"
            SSLog.statisticsError(message)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: message)
        }
    }
    
    /// Hyndman–Fan Type 1 (inverse empirical CDF).
    @inline(__always)
    internal static func quantileHF1(sorted: [FP], p: FP) -> FP? {
        guard !sorted.isEmpty else {
            SSLog.statisticsError("quantileHF1 requires a non-empty sorted sample")
            return nil
        }
        if p <= 0 { return sorted.first }
        if p >= 1 { return sorted.last }
        let n = FP(sorted.count)
        var index = Int(ceil(p * n))
        index = max(1, min(index, sorted.count))
        return sorted[index - 1]
    }
    
    /// Hyndman–Fan Type 2 (averaged ECDF).
    @inline(__always)
    internal static func quantileHF2(sorted: [FP], p: FP) -> FP? {
        guard !sorted.isEmpty else {
            SSLog.statisticsError("quantileHF2 requires a non-empty sorted sample")
            return nil
        }
        if p <= 0 { return sorted.first }
        if p >= 1 { return sorted.last }
        let n = FP(sorted.count)
        let h = p * n
        let floorH = FP(Int(h.rounded(.down)))
        if h == floorH {
            let idx = Int(floorH)
            if idx <= 0 { return sorted.first }
            if idx >= sorted.count { return sorted.last }
            return (sorted[idx - 1] + sorted[idx]) / FP(2)
        } else {
            var index = Int(ceil(h))
            index = max(1, min(index, sorted.count))
            return sorted[index - 1]
        }
    }
    
    /// Hyndman–Fan Type 3 (nearest order statistic).
    @inline(__always)
    internal static func quantileHF3(sorted: [FP], p: FP) -> FP? {
        guard !sorted.isEmpty else {
            SSLog.statisticsError("quantileHF3 requires a non-empty sorted sample")
            return nil
        }
        if p <= 0 { return sorted.first }
        if p >= 1 { return sorted.last }
        let n = FP(sorted.count)
        var index = Int((p * n).rounded())
        index = max(1, min(index, sorted.count))
        return sorted[index - 1]
    }
    
    /// Hyndman–Fan Types 4–9 via unified interpolation parameters.
    @inline(__always)
    internal static func quantileWithParams(
        p: FP,
        sorted: [FP],
        params: (a: FP, b: FP)
    ) -> FP? {
        guard !sorted.isEmpty else {
            SSLog.statisticsError("quantileWithParams requires a non-empty sorted sample")
            return nil
        }
        if p <= 0 { return sorted.first }
        if p >= 1 { return sorted.last }
        let n = FP(sorted.count)
        let (a, b) = params
        let h = a + (n + b) * p
        let j = Int(h.rounded(.down))
        let gamma = h - FP(j)
        if j < 1 {
            return sorted.first
        }
        if j >= sorted.count {
            return sorted.last
        }
        let lower = sorted[j - 1]
        let upper = sorted[j]
        return lower + gamma * (upper - lower)
    }
    
    
    /// Standardized Median Absolute Deviation (MAD) around a specified center under lock.
    ///
    /// Scales the MAD by a consistency factor for the Normal distribution so that,
    /// for normal data, the standardized MAD is a consistent estimator of the standard deviation.
    ///
    /// Parameters:
    /// - p: The center value around which absolute deviations are computed.
    ///
    /// Returns:
    /// - The standardized MAD, or `nil` if not applicable.
    ///
    /// Thread-safety:
    /// - The caller must hold the lock.
    ///
    /// Notes:
    /// - Uses the constant 1.482602218505601860547 (≈ 1 / Φ⁻¹(3/4)).
    @inline(__always)
    internal func standardizedMedianAbsoluteDeviationLocked(center p: FP) -> FP? {
        guard !self.isEmpty, self.isNumeric else {
            return self.logNil("standardizedMedianAbsoluteDeviationLocked requires a non-empty numeric dataset")
        }
        if let MAD = self.medianAbsoluteDeviationLocked(center: p) {
            return MAD * FP(1.482602218505601860547)
        }
        else {
            return self.logNil("standardizedMedianAbsoluteDeviationLocked could not compute the base MAD")
        }
    }
    
    /// Mean Absolute Deviation (about a specified center) under lock.
    ///
    /// Parameters:
    /// - p: The center value around which absolute deviations are computed.
    ///
    /// Returns:
    /// - The mean of |x − p| over all numeric items, or `nil` if not applicable.
    ///
    /// Thread-safety:
    /// - The caller must hold the lock.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting of the absolute deviations; O(n) to accumulate the mean.
    @inline(__always)
    internal func meanAbsoluteDeviationLocked(center p: FP) -> FP? {
        guard !self.isEmpty, self.isNumeric else {
            return self.logNil("meanAbsoluteDeviationLocked requires a non-empty numeric dataset")
        }
        if let values = self.itemsAsNumericArrayLocked() {
            var diffArray: [FP] = values.map { ($0 - p).magnitude }
            let a = diffArray.sorted { $0 < $1 }
            let n = a.count
            if n == 0 {
                return self.logNil("meanAbsoluteDeviationLocked found no numeric values")
            }
            if n == 1 { return a[0] }
            let sum = Helpers.sum(&diffArray)
            return sum / FP(n)
        }
        else {
            return self.logNil("meanAbsoluteDeviationLocked could not project items to numeric values")
        }
    }
    
    /// Two-sided confidence interval for the mean using the Normal distribution (locked variant).
    ///
    /// Computes a (1 − α) confidence interval for the population mean under the assumption
    /// that either:
    /// - the population standard deviation `sd` is known (provided via `populationSD`), or
    /// - the sample size is sufficiently large for the sample standard deviation to be used
    ///   as a plug-in estimate (CLT/normal approximation).
    ///
    /// Formula:
    /// - CI = mean ± z_{1−α/2} × (σ / √n)
    ///
    /// Parameters:
    /// - a: Significance level α in [0, 1]. The coverage is (1 − α).
    /// - sd: Optional known population standard deviation σ. If `nil`, uses
    ///       `populationStandardDeviation` from the dataset as an estimate.
    ///
    /// Returns:
    /// - A `ConfidenceInterval<FP>` on success; `nil` if the dataset is empty or non-numeric.
    ///
    /// Throws:
    /// - `SSError(.invalidArgument, ...)` if `alpha` is outside [0, 1].
    ///
    /// Thread-safety:
    /// - The caller must hold the lock.
    ///
    /// Complexity:
    /// - O(1) after mean and standard deviation are available; those may be cached.
    @inline(__always)
    internal func normalCILocked(alpha a: FP, populationSD sd: FP? = nil) throws -> ConfidenceInterval<FP>? {
        guard !self.isEmpty, self.isNumeric else {
            return self.logNil("normalCILocked requires a non-empty numeric dataset")
        }
        guard a.liesInClosedRange(from: 0, to: 1) else {
            let message = "alpha must be in [0,1]"
            SSLog.statisticsError(message)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: message)
        }
        var upper, lower, t1: FP
        do {
            if let m = self.arithmeticMeanLocked() {
                let u: FP = try SwiftyBoost.Distribution.Normal().quantile(1 - a / 2)
                let n = FP(self.count)
                if let s = sd {
                    t1 = s / n.squareRoot()
                }
                else {
                    if let sd = self.standardDeviationLocked(type: .population) {
                        t1 = sd / n.squareRoot()
                    }
                    else {
                        return self.logNil("normalCILocked requires a population standard deviation estimate")
                    }
                }
                upper = m + u * t1
                lower = m - u * t1
                let result: ConfidenceInterval<FP> = ConfidenceInterval(lowerBound: lower, upperBound: upper)
                return result
            }
            else {
                return self.logNil("normalCILocked could not compute the arithmetic mean")
            }
        }
        catch {
            return self.logNil("normalCILocked failed with error: \(error)")
        }
    }
    
    /// Arithmetic mean (average) under lock.
    ///
    /// Returns:
    /// - `nil` when the dataset is non-numeric or empty; otherwise the cached or freshly
    ///   computed total divided by the sample size.
    ///
    /// Thread-safety:
    /// - The caller must hold the lock.
    @inline(__always)
    internal func arithmeticMeanLocked() -> FP? {
        guard self.isNumeric, !self.isEmpty else {
            return self.logNil("arithmeticMeanLocked requires a non-empty numeric dataset")
        }
        updateTotalLocked()
        guard let t = self.totalStorage else {
            return self.logNil("arithmeticMeanLocked could not access totalStorage cache")
        }
        return t / FP(self.count)
    }
    
    /// Two-sided confidence interval for the mean using Student’s t distribution (locked variant).
    ///
    /// Computes a (1 − α) confidence interval for the population mean when the population
    /// standard deviation is unknown and the sample is assumed to be drawn from a normal
    /// distribution (or the sample size is moderate and t-approximation is acceptable).
    ///
    /// Formula:
    /// - CI = mean ± t_{ν, 1−α/2} × (s / √n), with ν = n − 1 degrees of freedom
    ///
    /// Parameters:
    /// - a: Significance level α in [0, 1]. The coverage is (1 − α).
    ///
    /// Returns:
    /// - A `ConfidenceInterval<FP>` on success; `nil` if the dataset is empty,
    ///   non-numeric, or the variance is undefined (e.g., n < 2).
    ///
    /// Throws:
    /// - `SSError(.invalidArgument, ...)` if `alpha` is outside [0, 1].
    ///
    /// Thread-safety:
    /// - The caller must hold the lock.
    ///
    /// Complexity:
    /// - O(1) after mean and sample standard deviation are available; those may be cached.
    @inline(__always)
    internal func studentTCILocked(alpha a: FP) throws -> ConfidenceInterval<FP>? {
        guard !self.isEmpty, self.isNumeric else {
            return self.logNil("studentTCILocked requires a non-empty numeric dataset")
        }
        guard a.liesInClosedRange(from: 0, to: 1) else {
            let message = "alpha must be in [0,1]"
            SSLog.statisticsError(message)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: message)
        }
        var upper, lower: FP
        do {
            if let m = self.arithmeticMeanLocked(), let s = self.standardDeviationLocked(type: .sample) {
                let t: FP = try SwiftyBoost.Distribution.StudentT(degreesOfFreedom: FP(self.count - 1)).quantile(1 - a / 2)
                let nsqrt = FP(self.count).squareRoot()
                lower = m - t * s / nsqrt
                upper = m + t * s / nsqrt
                let result: ConfidenceInterval<FP> = ConfidenceInterval(lowerBound: lower, upperBound: upper)
                return result
            }
            else {
                return self.logNil("studentTCILocked requires both mean and sample standard deviation")
            }
        }
        catch {
            return self.logNil("studentTCILocked failed with error: \(error)")
        }
    }
    
    /// Box–Whisker summary under lock.
    ///
    /// - Returns: A box-and-whisker descriptor with hinges, IQR, whiskers, notches, fences, and outliers,
    ///            or `nil` if the dataset is empty or non-numeric.
    /// - Throws: `SSError(.invalidArgument, ...)` on insufficient data for summary construction.
    /// - Complexity: O(n log n) for sorting; linear scans thereafter.
    internal func boxWhiskerLocked() throws -> BoxWhisker<FP>? {
        guard !self.isEmpty, self.isNumeric else {
            return nil
        }
        var res: BoxWhisker<FP> = BoxWhisker<FP>(median: FP.nan, q25: FP.nan, q75: FP.nan, iqr: FP.nan, lWhiskerExtreme: FP.nan, uWhiskerExtreme: FP.nan, fences: nil, outliers: nil, uNotch: FP.nan, lNotch: FP.nan)
        do {
            if let q = try self.quartile(), let aa = self.itemsAsNumericArrayLocked() {
                let q25 = q.lower
                let median = q.median
                let q75 = q.upper
                let iqr = q75 - q25
                let N = FP(self.sampleSize)
                let iqr3h = FP(1.5) * iqr
                let iqr3t = FP.two * iqr3h
                let notchCoeff = FP(1.57) * iqr / FP.sqrt(N)
                let upperLimitOutliers = q75 + iqr3t
                let upperLimitExtreme = q75 + iqr3h
                let lowerLimitOutliers = q25 - iqr3t
                let lowerLimitExtreme = q25 - iqr3h
                let a = aa.filter({!$0.isNaN }).sorted(by: >)
                if a.isEmpty {
                    throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Box-whisker calculation failed due to insufficient data")
                }
                res.median = median
                res.q25 = q25
                res.q75 = q75
                res.iqr = iqr
                res.uNotch = median + notchCoeff
                res.lNotch = median - notchCoeff
                res.fences = Array<FP>()
                res.outliers = Array<FP>()
                res.lWhiskerExtreme = q25
                res.uWhiskerExtreme = q75
                for value: FP in a {
                    if value > upperLimitOutliers {
                        res.outliers?.append(value)
                    }
                    else if value > upperLimitExtreme {
                        res.fences?.append(value)
                    }
                    else {
                        res.uWhiskerExtreme = value
                        break
                    }
                }
                for i in stride(from: a.count - 1, through: 0, by: -1) {
                    let value:FP = a[i]
                    if value < lowerLimitOutliers {
                        res.outliers?.append(value)
                    }
                    else if value < lowerLimitExtreme {
                        res.fences?.append(value)
                    }
                    else {
                        res.lWhiskerExtreme = value
                        break
                    }
                }
                return res
            }
            else {
                return nil
            }
        }
        catch {
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Box-whisker calculation failed due to insufficient data")
        }
    }
    
    /// Returns the median of a sorted numeric array.
    ///
    /// - Parameter sorted: A non-empty array sorted in ascending order.
    /// - Returns: The median value. For even length, the average of the two middle values (stable formula).
    /// - Precondition: `sorted.count > 0`.
    /// - Complexity: O(1).
    @inline(__always)
    fileprivate func medianOfSorted(_ sorted: [FP]) -> FP {
        let n = sorted.count
        precondition(n > 0, "medianOfSorted requires a non-empty array")
        
        let mid = n / 2
        if n % 2 == 1 {
            return sorted[mid]
        } else {
            let a = sorted[mid - 1]
            let b = sorted[mid]
            // more stable than (a + b) / 2
            return a + (b - a) / FP.two
        }
    }
    
    /// Base Sn computation (Rousseeuw & Croux): median_i median_j |x_i − x_j|.
    ///
    /// - Returns: The unscaled Sn base (no finite-sample correction, no asymptotic factor),
    ///   or `nil` if data are non-numeric or size < 2.
    /// - Complexity: O(n^2 log n) due to per-row median over n distances and final median.
    fileprivate func snBase() -> FP? {
        guard self.isNumeric, self.sampleSize >= 2, let a = self.itemsAsNumericArrayLocked() else {
            return nil
        }
        let n = self.sampleSize
        let sortedData = a.sorted(by: <)
        var medians = [FP]()
        medians.reserveCapacity(n)
        
        for i in 0..<n {
            let xi = sortedData[i]
            var row = [FP]()
            row.reserveCapacity(n)
            
            for j in 0..<n {
                let xj = sortedData[j]
                row.append((xi - xj).magnitude)
            }
            
            row.sort()
            let mi = medianOfSorted(row)
            medians.append(mi)
        }
        medians.sort()
        let baseSn = medianOfSorted(medians)
        return baseSn
    }
    
    /// Sn (Rousseeuw & Croux) with finite-sample correction.
    ///
    /// Definition:
    /// - Sn = 1.1926 · c_n · median_i median_j |x_i − x_j|
    ///
    /// Notes:
    /// - 1.1926 is the asymptotic consistency constant for normal data.
    /// - c_n is a finite-sample correction; for small n it follows tabulated values,
    ///   for larger n it approaches 1 (odd n uses n/(n−0.9)).
    ///
    /// Behavior:
    /// - Uses the exact O(n^2 log n) base for n ≤ 1000, otherwise a faster approximation.
    ///
    /// - Returns: Robust scale estimate Sn, or `nil` if data are invalid.
    internal func SnLocked() -> FP? {
        guard self.isNumeric, self.sampleSize >= 2 else { return nil }
        
        func sn(_ sn0: FP) -> FP? {
            let n = self.sampleSize
            let scale = FP(1.1926)   // asymptotic consistency
            
            // finite-sample correction c_n
            var cn = FP.one
            
            if n <= 9 {
                switch n {
                case 2: cn = FP(0.743)
                case 3: cn = FP(1.851)
                case 4: cn = FP(0.954)
                case 5: cn = FP(1.351)
                case 6: cn = FP(0.993)
                case 7: cn = FP(1.198)
                case 8: cn = FP(1.005)
                case 9: cn = FP(1.131)
                default: break
                }
            } else if n % 2 == 1 {
                // n odd, >= 11  --> cn = n / (n - 0.9)
                cn = FP(n) / (FP(n) - FP(0.9))
            }
            // for n even >= 10, cn = 1
            
            return scale * cn * sn0
        }
        
        if self.sampleSize <= 1000 {
            guard let sn0 = snBase() else { return nil }
            return sn(sn0)
        }
        else {
            guard let sn0 = snBaseFast() else { return nil }
            return sn(sn0)
        }
    }
    

    /// Base Qn computation (Rousseeuw & Croux): k-th order statistic of pairwise distances.
    ///
    /// Definition:
    /// - Qn0 is the k-th smallest value among all |x_i − x_j|, i < j, with
    ///   k = h(h−1)/2 and h = ⌊n/2⌋ + 1.
    ///
    /// - Returns: The unscaled Qn base (no finite-sample correction, no asymptotic factor),
    ///   or `nil` if data are non-numeric or size < 2.
    /// - Complexity: O(n^2 log n) to build and sort all pairwise distances.
    fileprivate func qnBase() -> FP? {
        guard self.isNumeric, self.sampleSize >= 2, let a = self.itemsAsNumericArrayLocked() else {
            return nil
        }
        
        let n = sampleSize
        let sortedData = a.sorted(by: <)
        
        let pairCount = n * (n - 1) / 2
        var distances = [FP]()
        distances.reserveCapacity(pairCount)
        
        for i in 0..<(n - 1) {
            let xi = sortedData[i]
            for j in (i + 1)..<n {
                let xj = sortedData[j]
                distances.append(abs(xi - xj))
            }
        }
        
        distances.sort()
        
        let h = n / 2 + 1
        let k = h * (h - 1) / 2   // 1-based
        let index = max(0, min(distances.count - 1, k - 1))
        let qn0 = distances[index]
        
        return qn0
    }

    /// Qn (Rousseeuw & Croux) with finite-sample correction.
    ///
    /// Definition:
    /// - Qn = 2.21914 · d_n · Qn0
    ///
    /// Notes:
    /// - 2.21914 is the asymptotic consistency constant for normal data.
    /// - d_n is a finite-sample correction depending on n (tabulated for n ≤ 12;
    ///   rational approximations for larger n; tends to 1 as n → ∞).
    ///
    /// Behavior:
    /// - Uses the exact O(n^2 log n) base for n < 500, otherwise a faster O(n log n) search.
    ///
    /// - Returns: Robust scale estimate Qn, or `nil` if data are invalid.
    internal func QnLocked() -> FP? {
        func qn(_ qn0: FP) -> FP? {
            let n = sampleSize
            let nn = FP(n)
            
            // Asymptotic consistency for sigma
            let scale = FP(2.21914)
            
            // Finite-sample correction d_n
            var dn = FP.one
            
            if n <= 12 {
                switch n {
                case 2:  dn = FP(0.399356)
                case 3:  dn = FP(0.99365)
                case 4:  dn = FP(0.51321)
                case 5:  dn = FP(0.84401)
                case 6:  dn = FP(0.61220)
                case 7:  dn = FP(0.85877)
                case 8:  dn = FP(0.66993)
                case 9:  dn = FP(0.87344)
                case 10: dn = FP(0.72014)
                case 11: dn = FP(0.88906)
                case 12: dn = FP(0.75743)
                default: break
                }
            } else {
                if n % 2 == 1 {
                    // odd n >= 13
                    let ex1 = FP(5.172) / nn
                    let ex2 = FP(-2.1284) - ex1
                    let ex3 = ex2 / nn
                    dn = FP(1.60188) + ex3
                } else {
                    // even n >= 14
                    let ex1 = FP(77.0) / nn
                    let ex2 = FP(6.987) - ex1
                    let ex3 = FP(1.9654) + ex2 / nn
                    dn = FP(3.67561) + ex3 / nn
                }
                // final transformation
                dn = FP.one / (dn / nn + FP.one)
            }
            
            return scale * dn * qn0
        }
        if self.sampleSize < 500 {
            guard let qn0 = qnBase() else { return nil }
            return qn(qn0)
        }
        else {
            guard let qn0 = qnBaseFast() else { return nil }
            return qn(qn0)
        }

    }
    /// Fast Qn base using a two-pointer pair-count and binary search on distances.
    ///
    /// Definition:
    /// - Finds the smallest t such that at least k pairwise distances |x_j − x_i| ≤ t,
    ///   where k = h(h−1)/2 and h = ⌊n/2⌋ + 1, on sorted data.
    ///
    /// - Returns: The approximate/unscaled Qn base, or `nil` if data are invalid.
    /// - Complexity: O(n log R) where R is the distance range (binary search iterations fixed ~60).
    fileprivate func qnBaseFast() -> FP? {
        guard self.isNumeric, sampleSize >= 2, let a = self.itemsAsNumericArrayLocked() else {
            return nil
        }
        
        let n = sampleSize
        // Sorted numeric data
        let sorted = a.sorted(by: <)
        
        // h = floor(n/2) + 1, k = h(h-1)/2
        let h = n / 2 + 1
        let k = h * (h - 1) / 2
        
        // Count number of pairs with distance <= t in O(n)
        @inline(__always)
        func countPairs(within t: FP) -> Int {
            var count = 0
            var j = 1
            
            for i in 0..<n {
                if j < i + 1 { j = i + 1 }
                while j < n && sorted[j] - sorted[i] <= t {
                    j += 1
                }
                count += max(0, j - i - 1)
            }
            return count
        }
        
        // Binary search bounds
        var low: FP = FP.zero
        var high: FP = sorted.last! - sorted.first!
        
        // If all values equal: Qn0 = 0
        if high <= 0 {
            return FP.zero
        }
        
        // Binary search for smallest t with countPairs(t) >= k
        // 60 iterations is enough to exhaust Double mantissa.
        let maxIter = 60
        for _ in 0..<maxIter {
            let mid = low + (high - low) / FP.two
            let cnt = countPairs(within: mid)
            if cnt >= k {
                high = mid
            } else {
                low = mid
            }
            
            // Stop if interval is on the order of machine epsilon
            let tol = FP.ulpOfOne * FP.maximum(FP.one, high.magnitude)
            if (high - low) <= tol {
                break
            }
        }
        
        // high is the smallest value for which countPairs(high) >= k (up to float precision)
        return high
    }

    /// Fast Sn base using per-row median of distances via Quickselect.
    ///
    /// - Returns: The unscaled Sn base computed with average O(n^2) time and reduced constants,
    ///   or `nil` if data are invalid.
    /// - Complexity: O(n^2) average using Quickselect for row medians, then O(n log n) for final median.
    fileprivate func snBaseFast() -> FP? {
        guard self.isNumeric, sampleSize >= 2, let a = self.itemsAsNumericArrayLocked() else {
            return nil
        }
        
        let n = sampleSize
        let sorted = a.sorted(by: <)
        
        var rowMedians = [FP]()
        rowMedians.reserveCapacity(n)
        
        for i in 0..<n {
            let xi = sorted[i]
            
            // reuse a single buffer
            var distances = [FP]()
            distances.reserveCapacity(n)
            for j in 0..<n {
                distances.append(abs(sorted[j] - xi))
            }
            
            // use Quickselect to get median in O(n) average:
            let mi = quickSelectMedian(&distances)
            rowMedians.append(mi)
        }
        
        // median of rowMedians (size n) – full sort is fine (O(n log n)):
        rowMedians.sort()
        let baseSn = medianOfSorted(rowMedians)
        
        return baseSn
    }

    /// Returns the k-th smallest element (0-based) of `array` in-place.
    ///
    /// - Parameters:
    ///   - array: Mutable array whose contents may be permuted.
    ///   - k: Zero-based order statistic index to select.
    /// - Returns: The k-th smallest value.
    /// - Preconditions: `array` must be non-empty and `0 ≤ k < array.count`.
    /// - Complexity: Average O(n); worst-case O(n^2).
    @inline(__always)
    fileprivate func quickSelectNth(
        _ array: inout [FP],
        k: Int
    ) -> FP {
        precondition(!array.isEmpty, "quickSelectNth requires a non-empty array")
        precondition(k >= 0 && k < array.count, "k out of range in quickSelectNth")

        var left = 0
        var right = array.count - 1

        while true {
            if left == right {
                return array[left]
            }

            // Choose pivot index (here: middle element)
            var pivotIndex = (left + right) / 2

            // Partition around pivot and get its final index
            pivotIndex = partition(&array, left: left, right: right, pivotIndex: pivotIndex)

            if k == pivotIndex {
                return array[k]
            } else if k < pivotIndex {
                right = pivotIndex - 1
            } else {
                left = pivotIndex + 1
            }
        }
    }

    /// Partitions `array[left...right]` around element at `pivotIndex` (Lomuto-style).
    ///
    /// - Parameters:
    ///   - array: Mutable array to partition (permuted in-place).
    ///   - left: Left boundary index (inclusive).
    ///   - right: Right boundary index (inclusive).
    ///   - pivotIndex: Index of the pivot element to partition around.
    /// - Returns: The final index of the pivot element after partitioning.
    /// - Complexity: O(n) for the subarray length.
    @inline(__always)
    fileprivate func partition(
        _ array: inout [FP],
        left: Int,
        right: Int,
        pivotIndex: Int
    ) -> Int {
        let pivotValue = array[pivotIndex]
        // Move pivot to the end
        array.swapAt(pivotIndex, right)
        var storeIndex = left

        for i in left..<right {
            if array[i] < pivotValue {
                array.swapAt(storeIndex, i)
                storeIndex += 1
            }
        }

        // Move pivot to its final place
        array.swapAt(storeIndex, right)
        return storeIndex
    }

    /// Computes the median using Quickselect (in-place).
    ///
    /// - Parameter array: Mutable array whose contents may be permuted.
    /// - Returns: The median value. For even length, the average of the two central order statistics.
    /// - Preconditions: `array` must be non-empty.
    /// - Complexity: Average O(n).
    fileprivate func quickSelectMedian(
        _ array: inout [FP]
    ) -> FP {
        precondition(!array.isEmpty, "quickSelectMedian requires a non-empty array")

        let n = array.count
        let mid = n / 2

        if n % 2 == 1 {
            // Odd length: the median is the element with index mid
            return quickSelectNth(&array, k: mid)
        } else {
            // Even length: median is average of order statistics at indices mid-1 and mid
            let lower = quickSelectNth(&array, k: mid - 1)
            let upper = quickSelectNth(&array, k: mid)
            // Numerically a bit safer than (lower + upper) / 2
            return lower + (upper - lower) / FP(2)
        }
    }

}
