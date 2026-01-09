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

/// Location estimators and quantile helpers for SSExamine data.

import SwiftyStatsPrelude
import SwiftyBoost

extension SSExamine {
    /// Arithmetic mean (average) of all numeric observations.
    ///
    /// Definition:
    /// - mean = (Σ value × frequency) / (Σ frequency) = `total / sampleSize`.
    ///
    /// Preconditions:
    /// - The dataset must be numeric (`isNumeric == true`) and non-empty.
    ///
    /// Return value:
    /// - `nil` for non-numeric datasets or when the dataset is empty.
    /// - Otherwise returns the mean as `FP`. If the cached `total` is `±∞` or `NaN`,
    ///   those special values propagate according to floating-point division rules.
    ///
    /// Caching and numerical stability:
    /// - The numerator `total` is maintained as a cached value that is recomputed on demand
    ///   using a numerically stable, magnitude-sorted compensated summation
    ///   (see `updateTotalLocked()` and `sum(_:)`).
    /// - Mutations invalidate the cache; reads re-use it when valid.
    ///
    /// Thread-safety:
    /// - Evaluated under `withLock` to read a consistent `count` and `total`.
    ///
    /// Complexity:
    /// - O(1) when `total` is already cached.
    /// - O(k log k) in the worst case if `total` must be recomputed for k distinct items
    ///   (due to stable summation by increasing magnitude).
    ///
    /// See also:
    /// - `total`, `squareTotal`, `poweredTotal(power:)`
    /// - `geometricMean`, `harmonicMean`, `powerMean(_:)`, `lehmerMean(_:)`
    public var arithmeticMean: FP? {
        get {
            return withLock { arithmeticMeanLocked() }
        }
    }
    

    /// Returns all items that occur with the highest frequency (the mode).
    ///
    /// Returns:
    /// - `nil` for empty datasets.
    /// - Otherwise, an array of items sharing the maximum absolute frequency.
    ///
    /// Notes:
    /// - Uses `.frequencyDescending` to find the maximum frequency, then collects ties.
    /// - The order of results follows the order provided by the frequency table.
    /// - If all items are equally frequent, all unique items are returned.
    ///
    /// Complexity:
    /// - O(k log k) for k distinct items, due to frequency table sorting.
    public var mode: Array<SSElement>? {
        get {
            return withLock { () -> [SSElement]? in
                if self.count == 0 { return nil }
                var result: [SSElement] = []
                let ft = self.frequencyTableLocked(sorted: .frequencyDescending)
                guard let top = ft.first else { return nil }
                let maxFrequency = top.frequency
                for e in ft where e.frequency == maxFrequency {
                    result.append(e.item)
                }
                return result
            }
        }
    }
    

    /// Alias for `mode` (the most common items).
    public var commonest: Array<SSElement>? {
        get {
            return mode
        }
    }
    

    /// Returns all items that occur with the lowest frequency (the scarcest).
    ///
    /// Returns:
    /// - `nil` for empty datasets.
    /// - Otherwise, an array of items sharing the minimum absolute frequency.
    ///
    /// Notes:
    /// - Uses `.frequencyAscending` to find the minimum frequency, then collects ties.
    /// - The order of results follows the order provided by the frequency table.
    ///
    /// Complexity:
    /// - O(k log k) for k distinct items, due to frequency table sorting.
    public var scarcest: Array<SSElement>? {
        get {
            return withLock { () -> [SSElement]? in
                if self.count == 0 { return nil }
                var result: [SSElement] = []
                let ft = self.frequencyTableLocked(sorted: .frequencyAscending)
                guard let low = ft.first else { return nil }
                let lowestFrequency = low.frequency
                for e in ft where e.frequency == lowestFrequency {
                    result.append(e.item)
                }
                return result
            }
        }
    }

    /// Quantile at probability q using a specified Hyndman–Fan convention.
    ///
    /// Parameters:
    /// - q: Probability in [0, 1]. Values outside this range throw `SSError(.invalidArgument, ...)`.
    /// - quantileType: The Hyndman–Fan type to use (HF1–HF9). Default is `.averagedECDF_HF2`.
    ///
    /// Returns:
    /// - `nil` for non-numeric or empty datasets.
    /// - Otherwise, the requested quantile as `FP`. Boundary behavior depends on the type:
    ///   - For HF2 (default): at q = 0 returns min; at q = 1 returns max; at integer n·q in 1...(n−1) returns midpoint of central two order stats; otherwise returns floor(n·q)-th order statistic (0-based).
    ///   - For HF1 and HF3–HF9: computed via the a–b–c–d parameterization on order statistics (see internal docs).
    ///
    /// Throws:
    /// - `SSError(.invalidArgument, ...)` when `q` is outside [0, 1] or the quantile type is unsupported.
    ///
    /// Thread-safety:
    /// - Computed under `withLock` from a sorted snapshot to avoid mutation races.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting; constant-time operations after sorting.
    public func quantile(q: FP, quantileType: HyndmanFanQuantileType = .type7) throws -> FP? {
        return try withLock { () -> FP? in
            try self.quantileLocked(q: q, quantileType: quantileType)
        }
    }
    

    /// Convenience accessor for the 25th, 50th (median), and 75th percentiles.
    ///
    /// Parameters:
    /// - quantileType: Hyndman–Fan type for all three quantiles. Default is HF2.
    ///
    /// Returns:
    /// - `Quartile<FP>` when the dataset is numeric and non-empty; otherwise `nil`.
    /// - Each quantile is evaluated on the same snapshot semantics as `quantile(q:quantileType:)`.
    ///
    /// Throws:
    /// - Propagates `SSError(.invalidArgument, ...)` from the underlying quantile computations.
    ///
    /// Thread-safety:
    /// - Each quantile is obtained under `withLock` using the same stable snapshot semantics as `quantile(_:)`.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting (work may be repeated per call; consider caching if needed).
    public func quartile(quantileType: HyndmanFanQuantileType = .type7) throws -> Quartile<FP>? {
        return try withLock { () -> Quartile<FP>? in
            if let q25 = try self.quantileLocked(q: 0.25, quantileType: quantileType),
               let q50 = try self.quantileLocked(q: 0.5, quantileType: quantileType),
               let q75 = try self.quantileLocked(q: 0.75, quantileType: quantileType) {
                return Quartile(q25: q25, q50: q50, q75: q75)
            }
            else {
                return nil
            }
        }
    }

    

    /// The median (50th percentile).
    ///
    /// Behavior:
    /// - Uses Hyndman–Fan Type 2 (averaged ECDF) semantics, matching `quantile(q: 0.5, quantileType: .averagedECDF_HF2)`.
    ///
    /// Returns:
    /// - `nil` for non-numeric or empty datasets.
    ///
    /// Thread-safety:
    /// - Computed under `withLock` via `quantileLocked(0.5)`.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting.
    public var median: FP? {
        withLock {
            do {
                return try self.quantileLocked(q: 0.5, quantileType: .type2)
            }
            catch _ {
                return nil
            }
        }
    }
    

    /// Geometric mean of all numeric values.
    ///
    /// Definition:
    /// - GM = exp((1/n) × Σ log(x_i)).
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets.
    /// - `0` if any value is exactly zero (since log(0) = −∞ ⇒ exp(−∞) = 0),
    ///   subject to the behavior of `logSumLocked()`.
    ///
    /// Thread-safety:
    /// - Uses `withLock` to read a consistent `sampleSize` and `logSum`.
    ///
    /// Complexity:
    /// - O(n log n) (log-sum uses stable summation).
    public var geometricMean: FP? {
        return withLock { () -> FP? in
            guard !self.isEmpty else { return nil }
            if let lp: FP = self.logSumLocked() {
                let c: FP = FP(self.sampleSize)
                return FP.exp(lp / c)
            }
            else {
                return nil
            }
        }
    }
    

    /// Harmonic mean (classical definition).
    ///
    /// Definition:
    /// - HM = n / (Σ 1/x_i), where the sum runs over all observations (with multiplicities).
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets.
    /// - May return ±∞ if any value is zero (since 1/0 is ±∞ and included in the sum).
    /// - May return `NaN` if reciprocals produce mixed infinities or `NaN`.
    ///
    /// Thread-safety:
    /// - Uses `withLock` to evaluate `inverseTotalLocked()` (Σ 1/x_i) and `sampleSize` atomically.
    ///
    /// Complexity:
    /// - O(k log k) for k distinct items (stable summation of reciprocals).
    public var harmonicMean: FP? {
        return withLock { () -> FP? in
            guard !self.isEmpty else { return nil }
            if let it: FP = self.inverseTotalLocked() {
                let c: FP = FP(self.sampleSize)
                return  c / it
            }
            else {
                return nil
            }
        }
    }
    

    /// Contra-harmonic mean of all numeric values.
    ///
    /// Definition:
    /// - C = (Σ x_i^2) / (Σ x_i)
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets.
    ///
    /// - May be ±∞ or `NaN` when `Σ x_i` is zero or when special values propagate.
    ///
    /// Thread-safety:
    /// - Uses `withLock` to read a consistent `squareTotal` and `total`.
    ///
    /// Complexity:
    /// - O(k log k) for k distinct items.
    public var contraharmonicMean: FP? {
        return withLock { () -> FP? in
            guard !self.isEmpty else { return nil }
            if let st = self.squareTotal, let total = self.total {
                return st / total
            }
            else {
                return nil
            }
        }
    }
    

    /// Lehmer mean of order `power` (L_p).
    ///
    /// Definition:
    /// - L_p = (Σ x_i^p) / (Σ x_i^(p-1))
    ///
    /// Special cases:
    /// - p = 1 → arithmetic mean.
    /// - p = 2 → contra-harmonic mean.
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets.
    /// - `nil` if the denominator Σ x_i^(p-1) is zero (undefined).
    /// - Otherwise, a finite value, ±∞, or `NaN` depending on inputs.
    ///
    /// Thread-safety:
    /// - Uses `withLock` to compute both sums from a consistent snapshot.
    ///
    /// Complexity:
    /// - O(k log k) for k distinct items.
    public func lehmerMean(_ power: FP) -> FP? {
        return withLock { () -> FP? in
            guard !self.isEmpty else { return nil }
            // Fast paths
            if power == 1 {
                return self.arithmeticMean
            } else if power == 2 {
                return self.contraharmonicMean
            }
            if let pSum = self.powerSumLocked(power: power),
               let pM1Sum = self.powerSumLocked(power: power - FP.one) {
                // Undefined when denominator is zero.
                if pM1Sum == 0 { return nil }
                return pSum / pM1Sum
            }
            else {
                return nil
            }
        }
    }
    

    /// Power mean (generalized mean) of order `power` (M_p).
    ///
    /// Definition (for p ≠ 0):
    /// - M_p = ( (1/n) × Σ x_i^p )^(1/p)
    ///
    /// Special cases:
    /// - p = 1 → arithmetic mean.
    /// - p = 2 → quadratic (RMS) mean.
    /// - p → 0 (limit) → geometric mean. We return `geometricMean` when `power == 0`.
    ///
    /// Domain notes:
    /// - For non-integer powers, negative x_i produce `NaN` via `pow`.
    /// - For negative powers, zeros produce ±∞ in the sum.
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets.
    /// - Otherwise, a finite value, ±∞, or `NaN` depending on inputs.
    ///
    /// Thread-safety:
    /// - Uses `withLock` to ensure a consistent snapshot of inputs.
    ///
    /// Complexity:
    /// - O(k log k) for k distinct items.
    public func powerMean(_ power: FP) -> FP? {
        return withLock { () -> FP? in
            guard !self.isEmpty else { return nil }
            // Special cases
            if power == 0 {
                return self.geometricMean
            } else if power == 1 {
                return self.arithmeticMean
            } else if power == 2 {
                // RMS: sqrt( (1/n) Σ x^2 )
                if let st = self.squareTotal {
                    let n: FP = FP(self.sampleSize)
                    return FP.sqrt(st / n)
                } else {
                    return nil
                }
            }
            if let pSum = self.powerSumLocked(power: power) {
                let n: FP = FP(self.sampleSize)
                return FP.pow(pSum / n, 1 / power)
            }
            else {
                return nil
            }
        }
    }
    

    /// Alpha-trimmed mean (two-sided), removing the lowest and highest `alpha` fraction.
    ///
    /// Definition:
    /// - Sort values ascending, remove `k = floor(alpha * n)` items from each tail,
    ///   and average the remaining middle slice.
    ///
    /// Parameters:
    /// - a: Trimming proportion in [0, 0.5].
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets.
    /// - `nil` if trimming removes all elements.
    ///
    /// Thread-safety:
    /// - Computes from a sorted snapshot under `withLock`.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting; O(m log m) for the stable sum over the middle m items.
    public func trimmedMean(alpha a: FP) -> FP? {
        return withLock { () -> FP? in
            guard !self.isEmpty else { return nil }
            // Validate alpha without throwing from this non-throwing API.
            guard a.liesInClosedRange(from: 0, to: FP(0.5)) else {
                SSLog.statisticsError("alpha must be in [0, 0.5]")
                return nil
            }
            if let v: [FP] = self.itemsAsNumericArrayLocked() {
                let values = v.sorted(by: <)
                let n = values.count
                if n == 0 { return nil }
                // k = floor(n * alpha)
                let v: FP = floor(FP(n) * a)
                let k: Int = Helpers.integerValue(v)
                let start = k
                let end = n - k - 1
                // If trimming removes all elements, undefined -> nil.
                if end < start { return nil }
                var terms: [FP] = Array<FP>(values[start...end]) // Convert slice to Array
                if terms.isEmpty { return nil }
                // Use the numerically stable sum helper.
                let s = Helpers.sum(&terms)
                return s / FP(terms.count)
            } else {
                return nil
            }
        }
    }
    

    /// Winsorized mean with two-sided winsorization at proportion `alpha`.
    ///
    /// Definition:
    /// - Let n = sample size, k = floor(alpha * n), capped at n/2.
    /// - Replace the smallest k values by the value at index k,
    ///   and the largest k values by the value at index n - k - 1 (both 0-based, ascending order).
    /// - Compute the arithmetic mean of the modified data (sample size remains n).
    ///
    /// Efficient computation:
    /// - Let s = sum(values[k ... n-k-1]) (the trimmed middle slice, possibly empty).
    /// - Then winsorized sum = s + k * (lowReplacement + highReplacement).
    ///
    /// Parameters:
    /// - a: Winsorization proportion in [0, 0.5].
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets.
    /// - Otherwise, a finite value, ±∞, or `NaN` depending on inputs.
    ///
    /// Thread-safety:
    /// - Operates on a sorted snapshot under `withLock`.
    ///
    /// Complexity:
    /// - O(n log n) due to sorting; summations are O(1) or O(m log m) for the middle slice of length m.
    public func winsorizedMean(alpha a: FP) -> FP? {
        return withLock { () -> FP? in
            guard !self.isEmpty else { return nil }
            guard a.liesInClosedRange(from: 0, to: FP(0.5)) else {
                SSLog.statisticsError("alpha must be in [0, 0.5]")
                return nil
            }
            guard let v: [FP] = self.itemsAsNumericArrayLocked() else {
                return nil
            }
            let values = v.sorted(by: <)
            let n = values.count
            if n == 0 { return nil }
            
            // k = floor(alpha * n), capped at n/2
            var k: Int = Helpers.integerValue(floor(a * FP(n)))
            if k > n / 2 { k = n / 2 }
            
            if k == 0 {
                // No winsorization; arithmetic mean via stable sum.
                var all = values
                let s = Helpers.sum(&all)
                return s / FP(n)
            }
            
            // Replacement values from order statistics
            let lowReplacement = values[k]
            let highReplacement = values[n - k - 1]
            
            // Sum the middle slice (if any) stably
            let start = k
            let end = n - k - 1
            var middleSum = FP.zero
            if start <= end {
                var middle = Array(values[start...end])
                middleSum = Helpers.sum(&middle)
            }
            
            // Winsorized sum = middleSum + k*(lowReplacement + highReplacement)
            var terms: [FP] = []
            terms.reserveCapacity(3)
            terms.append(middleSum)
            terms.append(FP(k) * lowReplacement)
            terms.append(FP(k) * highReplacement)
            let wSum = Helpers.sum(&terms)
            return wSum / FP(n)
        }
    }
    

    /// British spelling alias for `winsorizedMean(alpha:)`.
    ///
    /// - Parameter a: Winsorization proportion in [0, 0.5].
    /// - Returns: Same as `winsorizedMean(alpha:)`.
    public func winsorisedMean(alpha a: FP) -> FP? {
        return winsorizedMean(alpha: a)
    }
    

    /// Gastwirth robust location estimator.
    ///
    /// Definition (weights 1/3, 2/5, 1/3 on quantiles and median):
    /// - G = (1/3)·Q(1/3) + (2/5)·median + (1/3)·Q(2/3)
    ///
    /// Notes:
    /// - Quantiles and the median here use Hyndman–Fan Type 2 by default.
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets.
    ///
    /// Thread-safety:
    /// - Computes quantiles and median under `withLock` using stable snapshot semantics.
    ///
    /// Complexity:
    /// - O(n log n) due to quantile computations.
    public func gastwirthEstimate() -> FP? {
        return withLock { () -> FP? in
            guard !self.isEmpty else { return nil }
            let twothird = FP(2) / FP(3)
            let third = FP.one / FP(3)
            if let q_33 = try? self.quantile(q: third), let q_66 = try? self.quantile(q: twothird), let m = self.median {
                return third * q_33 + FP(2)/FP(5) * m + third * q_66
            }
            else {
                return nil
            }
        }
    }
    

    /// Two-sided confidence interval for the mean using the Normal distribution.
    ///
    /// Wrapper around `normalCILocked(alpha:populationSD:)` that acquires the instance lock,
    /// executes the locked variant, captures any thrown error, and rethrows it after releasing the lock.
    ///
    /// Assumptions:
    /// - Either a known population standard deviation `sd` is provided, or a normal approximation
    ///   is acceptable with the sample standard deviation as a plug-in estimate.
    ///
    /// Parameters:
    /// - a: Significance level α in [0, 1]. Coverage is (1 − α).
    /// - sd: Optional known population standard deviation σ. If `nil`, uses
    ///       `populationStandardDeviation` from the dataset as an estimate.
    ///
    /// Returns:
    /// - A `ConfidenceInterval<FP>` on success; `nil` if the dataset is empty or non-numeric.
    ///
    /// Throws:
    /// - `SSError(.invalidArgument, ...)` if `alpha` is outside [0, 1] (propagated from the locked variant).
    ///
    /// Thread-safety:
    /// - Acquires the instance lock via `withLock` to operate on a consistent snapshot.
    public func normalCI(alpha a: FP, sd: FP? = nil) throws -> ConfidenceInterval<FP>? {
        var capturedError: Error?
        let result: ConfidenceInterval<FP>? = withLock { () -> ConfidenceInterval<FP>? in
            do {
                return try normalCILocked(alpha: a, populationSD: sd)
            } catch {
                capturedError = error
                return nil
            }
        }
        if let e = capturedError {
            throw e
        }
        return result
    }
    

    /// Two-sided confidence interval for the mean using Student’s t distribution.
    ///
    /// Wrapper around `studentTCILocked(alpha:)` that acquires the instance lock,
    /// executes the locked variant, captures any thrown error, and rethrows it after releasing the lock.
    ///
    /// Assumptions:
    /// - Population standard deviation is unknown; data are approximately normal
    ///   (or t-approximation is acceptable). Degrees of freedom ν = n − 1.
    ///
    /// Parameters:
    /// - a: Significance level α in [0, 1]. Coverage is (1 − α).
    ///
    /// Returns:
    /// - A `ConfidenceInterval<FP>` on success; `nil` if the dataset is empty, non-numeric,
    ///   or variance is undefined (e.g., n < 2).
    ///
    /// Throws:
    /// - `SSError(.invalidArgument, ...)` if `alpha` is outside [0, 1] (propagated from the locked variant).
    ///
    /// Thread-safety:
    /// - Acquires the instance lock via `withLock` to operate on a consistent snapshot.
    public func studentTCI(alpha a: FP) throws -> ConfidenceInterval<FP>? {
        var capturedError: Error?
        let result: ConfidenceInterval<FP>? = withLock { () -> ConfidenceInterval<FP>? in
            do {
                return try studentTCILocked(alpha: a)
            } catch {
                capturedError = error
                return nil
            }
        }
        if let e = capturedError {
            throw e
        }
        return result
    }

    /// Convenience 95% two-sided confidence interval for the mean using Student’s t distribution.
    ///
    /// Behavior:
    /// - Computes a (1 − α) CI with α = 0.05 by delegating to `studentTCILocked(alpha:)`
    ///   under the instance lock.
    ///
    /// Returns:
    /// - A `ConfidenceInterval<FP>` on success.
    /// - `nil` if the dataset is empty, non-numeric, or the variance is undefined (e.g., sample size < 2).
    ///
    /// Notes:
    /// - This is a read-only convenience accessor equivalent to `try? studentTCILocked(alpha: 0.05)`.
    /// - It does not consult or modify any stored alpha; the level is fixed at 0.05 here.
    ///
    /// Thread-safety:
    /// - Acquires the instance lock via `withLock` to operate on a consistent snapshot.
    public var meanCI: ConfidenceInterval<FP>? {
        return withLock { () -> ConfidenceInterval<FP>? in
            do {
                return try studentTCILocked(alpha: 0.05)
            }
            catch _ {
                return nil
            }
        }
    }

}
