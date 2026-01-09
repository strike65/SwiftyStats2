//
//  Created by VT on 28.10.25.
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

/// Scatter and dispersion measures for SSExamine collections.

import SwiftyStatsPrelude
import SwiftyBoost

extension SSExamine {
    /// Largest element in the dataset (by natural order of `SSElement`).
    ///
    /// Returns:
    /// - `nil` if the dataset is empty; otherwise the maximum element.
    ///
    /// Thread-safety:
    /// - Obtains a sorted snapshot under `withLock`.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting.
    public var largest: SSElement? {
        return withLock { () -> SSElement? in
            guard !self.isEmpty else {
                return self.logNil("largest requires a non-empty dataset")
            }
            return self.itemsAsArray(sorted: .descending).first
        }
    }
    

    /// Smallest element in the dataset (by natural order of `SSElement`).
    ///
    /// Returns:
    /// - `nil` if the dataset is empty; otherwise the minimum element.
    ///
    /// Thread-safety:
    /// - Obtains a sorted snapshot under `withLock`.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting.
    public var smallest: SSElement? {
        return withLock { () -> SSElement? in
            guard !self.isEmpty else {
                return self.logNil("smallest requires a non-empty dataset")
            }
            return self.itemsAsArray(sorted: .ascending).first
        }
    }
    

    /// Smallest and largest elements in one pass.
    ///
    /// Returns:
    /// - `(smallest, largest)` tuple when non-empty; otherwise `nil`.
    ///
    /// Thread-safety:
    /// - Obtains a sorted ascending snapshot under `withLock`.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting.
    public var smallestAndLargest: (smallest: SSElement, largest: SSElement)? {
        return withLock { () -> (smallest: SSElement, largest: SSElement)? in
            guard !self.isEmpty else {
                return self.logNil("smallestAndLargest requires a non-empty dataset")
            }
            let sorted = self.itemsAsArray(sorted: .ascending)
            if let first = sorted.first, let last = sorted.last {
                return (smallest: first, largest: last)
            }
            else {
                return self.logNil("smallestAndLargest could not determine extremes")
            }
        }
    }
    

    /// Range of the dataset (max − min) for numeric values.
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets.
    /// - Otherwise, `max - min` using the numeric representation of elements.
    ///
    /// Thread-safety:
    /// - Operates under `withLock`.
    ///
    /// Complexity:
    /// - O(n) via a single pass to find min and max.
    public var range: FP? {
        return withLock { () -> FP? in
            guard !self.isEmpty else {
                return self.logNil("range requires a non-empty dataset")
            }
            guard let values: [FP] = self.itemsAsNumericArrayLocked(), !values.isEmpty else {
                return self.logNil("range requires numeric data")
            }
            if let minV = values.min(), let maxV = values.max() {
                return maxV - minV
            }
            else {
                return self.logNil("range could not determine min/max")
            }
        }
    }
    

    /// Interquantile range |Q(lower) − Q(upper)| for arbitrary quantiles.
    ///
    /// - Parameters:
    ///   - lower: Lower quantile in [0, 1].
    ///   - upper: Upper quantile in [0, 1]; must satisfy `upper ≥ lower`.
    /// - Returns: The absolute difference |Q(lower) − Q(upper)|, or `nil` if not applicable or quantiles are undefined.
    ///
    /// Notes:
    /// - Uses Hyndman–Fan Type 2 quantiles internally via `interquantileRangeLocked`.
    /// - For discrete or small samples, quantile definitions may produce ties or stepwise results.
    ///
    /// Throws:
    /// - `SSError(.invalidArgument, ...)` when `lower` or `upper` fall outside [0, 1] or when type checks fail upstream.
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and delegates to the locked variant.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting for quantile computation.
    public func interquantileRange(lower: FP, upper: FP) throws -> FP? {
        return try withLock { () -> FP? in
            try self.interquantileRangeLocked(lower: lower, upper: upper)
        }
    }
    

    /// Interquartile range (IQR) = Q3 − Q1 using Hyndman–Fan Type 2 quantiles.
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets; otherwise the IQR.
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and computes via `interquantileRangeLocked(lower: 0.25, upper: 0.75)`.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting for quantile computation.
    public var interquartileRange: FP? {
        return withLock { () -> FP? in
            return try? self.interquantileRangeLocked(lower: 0.25, upper: 0.75)
        }
    }
    

    /// Mid-range = (min + max) / 2 for numeric datasets.
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets, or if either bound cannot be converted to `FP`.
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and uses `smallest` and `largest` under consistent snapshot semantics.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting for min/max retrieval.
    public var midRange: FP? {
        return withLock { () -> FP? in
            do {
                guard !self.isEmpty, self.isNumeric else {
                    return self.logNil("midRange requires a non-empty numeric dataset")
                }
                guard self.count > 1 else {
                    return self.logNil("midRange requires at least two observations")
                }
                if let min = self.smallest, let max = self.largest {
                    let minVal: FP = try RealConverter.to(min)
                    let maxVal: FP = try RealConverter.to(max)
                    return (minVal + maxVal) / 2
                }
                else {
                    return self.logNil("midRange could not determine smallest/largest values")
                }
            }
            catch _ {
                return self.logNil("midRange failed to convert bounds to numeric values")
            }
        }
    }
    

    /// Quartile deviation (semi-interquartile range) = IQR / 2.
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets or when IQR is unavailable.
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and derives from `interquartileRange`.
    ///
    /// Complexity:
    /// - O(n log n) dominated by IQR computation.
    public var quartileDeviation: FP? {
        get {
            return withLock { () -> FP? in
                guard !self.isEmpty, self.isNumeric else {
                    return self.logNil("quartileDeviation requires a non-empty numeric dataset")
                }
                guard self.count > 1 else {
                    return self.logNil("quartileDeviation requires at least two observations")
                }
                if let iqr = self.interquartileRange {
                    return iqr / 2
                }
                else {
                    return self.logNil("quartileDeviation could not compute interquartile range")
                }
            }
        }
    }
    

    /// Relative quartile distance = (Q3 − Q1) / median.
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets, or if the median is zero/undefined.
    ///
    /// Notes:
    /// - Uses Hyndman–Fan Type 2 quantiles and median.
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and computes from a consistent snapshot.
    ///
    /// Complexity:
    /// - O(n log n) due to quantile computations.
    public var relativeQuartileDistance: FP? {
        get {
            return withLock { () -> FP? in
                guard !self.isEmpty, self.isNumeric else {
                    return self.logNil("relativeQuartileDistance requires a non-empty numeric dataset")
                }
                guard self.count > 1 else {
                    return self.logNil("relativeQuartileDistance requires at least two observations")
                }
                if let quart = try? self.quartile() {
                    if quart.median.isZero {
                        return self.logNil("relativeQuartileDistance undefined when the median is zero")
                    }
                    return (quart.upper - quart.lower) / quart.median
                }
                else {
                    return self.logNil("relativeQuartileDistance could not compute quartiles")
                }
            }
        }
    }
    
    /// Sample variance s² (unbiased estimator).
    ///
    /// Definition:
    /// - Uses the unbiased denominator (n − 1).
    ///
    /// Returns:
    /// - `nil` for non-numeric datasets or when `count < 2`.
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and delegates to `varianceLocked(type: .sample)`.
    ///
    /// Complexity:
    /// - O(n log n) dominated by moment/mean computations if not cached.
    public var sampleVariance: FP? {
        return withLock { () -> FP? in
            if self.isEmpty || !self.isNumeric {
                return self.logNil("sampleVariance requires a non-empty numeric dataset")
            }
            return varianceLocked(type: .sample)
        }
    }
    

    /// Population variance σ² (biased estimator).
    ///
    /// Definition:
    /// - Uses the biased denominator n.
    ///
    /// Returns:
    /// - `nil` for non-numeric datasets or when `count < 2` (undefined variance).
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and delegates to `varianceLocked(type: .population)`.
    ///
    /// Complexity:
    /// - O(n log n) dominated by moment/mean computations if not cached.
    public var populationVariance: FP? {
        return withLock { () -> FP? in
            if self.isEmpty || !self.isNumeric {
                return self.logNil("populationVariance requires a non-empty numeric dataset")
            }
            return varianceLocked(type: .population)
        }
    }
    

    /// Alias for `populationVariance` (biased estimator).
    public var biasedVariance: FP? {
        return populationVariance
    }
    

    /// Alias for `sampleVariance` (unbiased estimator).
    public var unbiasedVariance: FP? {
        return sampleVariance
    }
    

    /// Population standard deviation (sqrt of population variance).
    ///
    /// Returns:
    /// - `nil` for non-numeric datasets or when variance is undefined.
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and derives from `varianceLocked(type: .population)`.
    public var populationStandardDeviation: FP? {
        return withLock { () -> FP? in
            return standardDeviationLocked(type: .population)
        }
    }
    

    /// Sample standard deviation (sqrt of sample variance).
    ///
    /// Returns:
    /// - `nil` for non-numeric datasets or when variance is undefined (e.g., `count < 2`).
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and derives from `varianceLocked(type: .sample)`.
    public var sampleStandardDeviation: FP? {
        return withLock { () -> FP? in
            return standardDeviationLocked(type: .sample)
        }
    }
    
    /// Alias for ``sampleStandardDeviation`` (uses Bessel's correction; divides by `n - 1`).
    public var standardDeviationUnbiased: FP? {
        return withLock { () -> FP? in
            return sampleStandardDeviation
        }
    }
    
    /// Alias for ``populationStandardDeviation`` (divides by `n`).
    public var standardDeviationBiased: FP? {
        return withLock { () -> FP? in
            return populationStandardDeviation
        }
    }

    /// Upper semi-variance (upside): mean of squared deviations for observations strictly above the mean.
    ///
    /// Returns:
    /// - `nil` for non-numeric datasets or when fewer than two observations are present.
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and delegates to `semiVarianceLocked(type: .upper)`.
    public var semiVarianceUpper: FP? {
        return withLock { () -> FP? in
            if self.isEmpty || !self.isNumeric {
                return self.logNil("semiVarianceUpper requires a non-empty numeric dataset")
            }
            guard let varSam = self.semiVarianceLocked(type: .upper) else {
                return self.logNil("semiVarianceUpper requires at least two observations above the mean")
            }
            return varSam
        }
    }
    

    /// Lower semi-variance (downside relative to mean): mean of squared deviations for observations strictly below the mean.
    ///
    /// Returns:
    /// - `nil` for non-numeric datasets or when fewer than two observations are present.
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and delegates to `semiVarianceLocked(type: .lower)`.
    public var semiVarianceLower: FP? {
        return withLock { () -> FP? in
            if self.isEmpty || !self.isNumeric {
                return self.logNil("semiVarianceLower requires a non-empty numeric dataset")
            }
            guard let varSam = self.semiVarianceLocked(type: .lower) else {
                return self.logNil("semiVarianceLower requires at least two observations below the mean")
            }
            return varSam
        }
    }
    

    /// Downside semi-variance: mean of squared negative departures from the mean (positive deviations treated as zero).
    ///
    /// Returns:
    /// - `nil` for non-numeric datasets or when fewer than two observations are present.
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and delegates to `semiVarianceLocked(type: .downside)`.
    public var semiVarianceDownside: FP? {
        return withLock { () -> FP? in
            if self.isEmpty || !self.isNumeric {
                return self.logNil("semiVarianceDownside requires a non-empty numeric dataset")
            }
            guard let varSam = self.semiVarianceLocked(type: .downside) else {
                return self.logNil("semiVarianceDownside requires at least two observations")
            }
            return varSam
        }
    }
    

    /// Upper mean absolute deviation about a reference point (locked).
    ///
    /// Computes the mean of absolute deviations |x − from| over the subset of
    /// observations strictly greater than `from`.
    ///
    /// - Parameters:
    ///   - from: The reference value about which absolute deviations are measured.
    ///
    /// - Returns:
    ///   - The average absolute deviation of all values `x` with `x > from`,
    ///     or `nil` if the dataset is empty, non-numeric, or numeric projection fails.
    ///     If no observations satisfy `x > from`, this method returns `0` (the mean
    ///     over an empty subset is treated as zero in this implementation).
    ///
    /// - Thread safety:
    ///   - The caller must hold the `SSExamine` lock. Public wrappers acquire the lock
    ///     and delegate to this method.
    ///
    /// - Complexity: O(n), dominated by a linear scan and a stable summation.
    public func upperMeanAbsoluteDeviation(from: FP) -> FP? {
        return withLock { () -> FP? in
            return self.upperMeanAbsoluteDeviationLocked(from: from)
        }
    }

    /// Lower mean absolute deviation about a reference point (locked).
    ///
    /// Computes the mean of absolute deviations |x − from| over the subset of
    /// observations strictly less than `from`.
    ///
    /// - Parameters:
    ///   - from: The reference value about which absolute deviations are measured.
    ///
    /// - Returns:
    ///   - The average absolute deviation of all values `x` with `x < from`,
    ///     or `nil` if the dataset is empty, non-numeric, or numeric projection fails.
    ///     If no observations satisfy `x < from`, this method returns `0` (the mean
    ///     over an empty subset is treated as zero in this implementation).
    ///
    /// - Thread safety:
    ///   - The caller must hold the `SSExamine` lock. Public wrappers acquire the lock
    ///     and delegate to this method.
    ///
    /// - Complexity: O(n), dominated by a linear scan and a stable summation.
    public func lowerMeanAbsoluteDeviation(from: FP) -> FP? {
        return withLock { () -> FP? in
            return self.lowerMeanAbsoluteDeviationLocked(from: from)
        }
    }
    

    /// Standard error of the mean (SEM) = s / √n using sample standard deviation.
    ///
    /// Returns:
    /// - `nil` for non-numeric datasets or when `sampleStandardDeviation` is unavailable.
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and reads a consistent `count` and `sampleStandardDeviation`.
    public var standardError: FP? {
        return withLock { () -> FP? in
            if self.isEmpty || !self.isNumeric {
                return self.logNil("standardError requires a non-empty numeric dataset")
            }
            guard let sampleSd = self.standardDeviationLocked(type: .sample) else {
                return self.logNil("standardError requires an available sample standard deviation")
            }
            let n = FP(self.count)
            return sampleSd / n.squareRoot()
        }
    }
        

    /// Median Absolute Deviation (MAD) about a specified center.
    ///
    /// Definition:
    /// - MAD(p) = median(|x_i − p|), where p is a supplied center (often the sample median).
    /// - Uses the same quantile semantics as `quantile(q: 0.5, quantileType: .averagedECDF_HF2)` for the deviation sample.
    ///
    /// Returns:
    /// - `nil` if the dataset is empty or non-numeric.
    /// - Otherwise, the 50th percentile of absolute deviations from `p`.
    ///
    /// Thread-safety:
    /// - Acquires a consistent snapshot via `withLock`.
    ///
    /// Notes:
    /// - For a robust estimate of scale, see `standardizedMedianAbsoluteDeviation(center:)`.
    public func medianAbsoluteDeviation(center p: FP) -> FP? {
        return withLock { () -> FP? in
            medianAbsoluteDeviationLocked(center: p)
        }
    }
    

    /// Standardized Median Absolute Deviation (MAD) about a specified center.
    ///
    /// Definition:
    /// - sMAD(p) = 1.482602218505601860547 × MAD(p), a consistency scaling factor for Normal data,
    ///   making sMAD a consistent estimator of the standard deviation under normality.
    ///
    /// Returns:
    /// - `nil` if the dataset is empty or non-numeric.
    /// - Otherwise, the standardized MAD as `FP`.
    ///
    /// Thread-safety:
    /// - Acquires a consistent snapshot via `withLock`.
    ///
    /// See also:
    /// - `medianAbsoluteDeviation(center:)`
    public func standardizedMedianAbsoluteDeviation(center p: FP) -> FP? {
        return withLock { () -> FP? in
            standardizedMedianAbsoluteDeviationLocked(center: p)
        }
    }
    

    /// Mean absolute deviation about a supplied center.
    ///
    /// Definition:
    /// - MAD₁(p) = mean(|x_i − p|) computed on the numeric projection of the dataset.
    ///
    /// Returns:
    /// - `nil` when the dataset is empty or non-numeric.
    /// - Otherwise, the arithmetic mean of absolute deviations from `p`.
    ///
    /// Thread-safety:
    /// - Acquires the container lock and defers to `meanAbsoluteDeviationLocked(center:)`.
    public func meanAbsoluteDeviation(center p: FP) -> FP? {
        return withLock { () -> FP? in
            meanAbsoluteDeviationLocked(center: p)
        }
    }

    

    /// Coefficient of variation (CV) = s / mean using sample standard deviation.
    ///
    /// Returns:
    /// - `nil` for non-numeric datasets or when either the sample standard deviation
    ///   or the arithmetic mean is unavailable (or mean is zero).
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and reads a consistent snapshot.
    ///
    /// Note:
    /// - See `cv` as a preferred alias.
    public var coefficientOfVariation: FP? {
        return withLock { () -> FP? in
            if self.isEmpty || !self.isNumeric {
                return self.logNil("coefficientOfVariation requires a non-empty numeric dataset")
            }
            guard let sampleSd = self.standardDeviationLocked(type: .sample),
                  let m = self.arithmeticMeanLocked() else {
                return self.logNil("coefficientOfVariation requires both sample SD and mean")
            }
            return sampleSd / m
        }
    }
    

    /// Preferred alias for the coefficient of variation (CV) = s / mean.
    ///
    /// Returns:
    /// - Same as `coefficientOfVariation`.
    ///
    /// Thread-safety:
    /// - Delegates to `coefficientOfVariation` under lock.
    public var cv: FP? {
        return withLock { () -> FP? in
            if self.isEmpty || !self.isNumeric {
                return self.logNil("cv requires a non-empty numeric dataset")
            }
            guard let sampleSd = self.standardDeviationLocked(type: .sample),
                  let m = self.arithmeticMeanLocked() else {
                return self.logNil("cv requires both sample SD and mean")
            }
            return sampleSd / m
        }
    }
    

    /// Detects outliers using Rosner’s generalized ESD procedure.
    ///
    /// This is a convenience wrapper that:
    /// - extracts the numeric values from this `SSExamine` instance, and
    /// - calls `Inferential.HypothesisTesting.Outliers.rosnerESD(data:alpha:maxOutliers:testType:)`.
    ///
    /// - Parameters:
    ///   - alpha: Significance level in (0, 1).
    ///   - maxOutliers: Upper bound on the number of outliers to test (1 ≤ k ≤ n − 2).
    ///   - tails: Tail configuration (`.bothTails`, `.upperTail`, `.lowerTail`).
    /// - Returns: Array of detected outliers in the order they were removed, or `nil` if unavailable.
    ///
    /// Assumptions:
    /// - The underlying distribution is approximately normal after removing up to the true outliers.
    ///
    /// Notes:
    /// - If you need single-outlier testing, see `Inferential.HypothesisTesting.Outliers.grubbsTest(...)`.
    /// - Rosner’s ESD mitigates masking when multiple outliers are present.
    ///
    /// Thread-safety:
    /// - Operates under `withLock` and works on a snapshot of the current data.
    public func outliers(alpha: FP, maxOutliers: Int, tails: ESDTestType) -> [FP]? {
        return withLock { () -> [FP]? in
            if let data = self.itemsAsNumericArrayLocked() {
                do {
                    if let rt: Outliers<FP>.ESDTestResult<FP>  = try Outliers<FP>.rosnerESD(data: data, alpha: alpha, maxOutliers: maxOutliers, testType: tails) {
                        if let outliers = rt.outliers {
                            return outliers
                        }
                        else {
                            return self.logNil("outliers(alpha:maxOutliers:tails:) returned no outlier list from Rosner ESD")
                        }
                    }
                    else {
                        return self.logNil("outliers(alpha:maxOutliers:tails:) received nil result from Rosner ESD")
                    }
                }
                catch {
                    return self.logNil("outliers(alpha:maxOutliers:tails:) failed due to Rosner ESD error: \(error)")
                }
            }
            else {
                return self.logNil("outliers(alpha:maxOutliers:tails:) requires numeric data")
            }
        }
    }
    

    /// Runs Rosner’s ESD test and returns the full result record.
    ///
    /// This convenience method:
    /// - extracts a numeric snapshot from this `SSExamine` instance, and
    /// - invokes `Inferential.HypothesisTesting.Outliers.rosnerESD(...)` using the instance’s `alpha`
    ///   (if available) alongside the provided `maxOutliers` and `tails`.
    ///
    /// - Parameters:
    ///   - maxOutliers: Upper bound on outliers to search for (1 ≤ k ≤ n − 2).
    ///   - tails: Tail configuration (`.bothTails`, `.upperTail`, `.lowerTail`).
    /// - Returns: `ESDTestResult<FP>` describing the sequential statistics and detected outliers, or `nil` if unavailable.
    ///
    /// See also:
    /// - `Inferential.HypothesisTesting.Outliers.rosnerESD(data:alpha:maxOutliers:testType:)`
    /// - `Inferential.HypothesisTesting.Outliers.grubbsTest(values:mean:standardDeviation:alpha:)`
    ///
    /// Thread-safety:
    /// - Operates under `withLock` and uses a consistent snapshot of the data.
    public func rosnerESDTest(maxOutliers: Int, tails: ESDTestType) -> Outliers<FP>.ESDTestResult<FP>? {
        return withLock { () -> Outliers<FP>.ESDTestResult<FP>? in
            do {
                if let data = self.itemsAsNumericArrayLocked() {
                    if let esdResult: Outliers<FP>.ESDTestResult<FP> = try Outliers.rosnerESD(data: data, alpha: alpha, maxOutliers: maxOutliers, testType: tails) {
                        return esdResult
                    }
                    else {
                        return self.logNil("rosnerESDTest(maxOutliers:tails:) received nil result from Rosner ESD")
                    }
                }
                else {
                    return self.logNil("rosnerESDTest(maxOutliers:tails:) requires numeric data")
                }
                
            }
            catch {
                return self.logNil("rosnerESDTest(maxOutliers:tails:) failed due to Rosner ESD error: \(error)")
            }
        }
    }
    

    /// Box–Whisker summary under lock.
    ///
    /// - Returns: A box-and-whisker descriptor with hinges, IQR, whiskers, notches, fences, and outliers,
    ///            or `nil` if the dataset is empty or non-numeric.
    /// - Throws: `SSError(.invalidArgument, ...)` on insufficient data for summary construction.
    /// - Complexity: O(n log n) for sorting; linear scans thereafter.
    public func boxWhisker() throws ->  BoxWhisker<FP>? {
        return try withLock { () -> BoxWhisker<FP>? in
            do {
                if let res = try boxWhiskerLocked() {
                    return res
                }
                else {
                    return nil
                }
            }
            catch let e {
                throw e
            }
        }
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
    public var Sn: FP? {
        return withLock { () -> FP? in
            SnLocked()
        }
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
    public var Qn: FP? {
        return withLock { () -> FP? in
            QnLocked()
        }
    }
}

