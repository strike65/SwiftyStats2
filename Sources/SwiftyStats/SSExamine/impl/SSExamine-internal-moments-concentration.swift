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

/// Internal moment and concentration calculations for SSExamine.

import SwiftyStatsPrelude
import SwiftyBoost

// MARK: moments, concentration
/// Internal kernels for dispersion metrics computed under the SSExamine lock.
///
/// These helpers evaluate moments `Σ (x-μ)^p`, semivariance partitions, and quantile-based
/// robust statistics while keeping numerical stability via integer exponentiation and
/// compensated products. They provide the algebraic backbone for public variance,
/// kurtosis, and outlier detectors without re-locking the container.
extension SSExamine {
    /// Computes the raw (unnormalized) moment Σ (v − x)^p under lock.
    ///
    /// Optimizations:
    /// - If `x == 0`, delegates to `powerSumLocked(power:)` for a faster path.
    /// - If `p == 0`, returns the count as `FP`.
    /// - Specialized loops for `p == 1` and `p == 2` avoid generic exponentiation.
    ///
    /// - Parameters:
    ///   - p: Order of the moment (non-negative integer as `UInt`).
    ///   - x: Reference point (e.g., 0 for raw moments, mean for central moments).
    /// - Returns: The raw (unnormalized) moment or `nil` if no values exist. Returns `FP.nan` if any term is `NaN`.
    /// - Thread safety: The caller must hold the lock.
    /// - Complexity: O(n).
    @inline(__always)
    internal func rawMomentLocked(order p: UInt, about x: FP) -> FP? {
        guard let values: [FP] = self.itemsAsNumericArrayLocked() else {
            return self.logNil("rawMomentLocked requires numeric data but itemsAsNumericArrayLocked returned nil")
        }
        guard !values.isEmpty else {
            return self.logNil("rawMomentLocked requires at least one observation")
        }
        if x.isZero {
            return self.powerSumLocked(power: FP(p))
        }
        if p == 0 {
            return FP(values.count)
        }
        /// Fast integer exponentiation for `FP` base with a non-negative integer exponent.
        @inline(__always)
        func integerPower(_ value: FP, exponent: UInt) -> FP {
            if exponent == 0 { return FP.one }
            var result: FP = FP.one
            var base = value
            var exp = exponent
            while exp > 0 {
                if exp & 1 == 1 {
                    result *= base
                }
                exp >>= 1
                if exp > 0 {
                    base *= base
                }
            }
            return result
        }
        var terms: [FP] = []
        terms.reserveCapacity(values.count)
        if p == 2 {
            for v in values {
                let diff = v - x
                let term = diff * diff
                if term.isNaN { return FP.nan }
                terms.append(term)
            }
        } else if p == 1 {
            for v in values {
                let term = v - x
                if term.isNaN { return FP.nan }
                terms.append(term)
            }
        } else {
            for v in values {
                let diff = v - x
                let term = integerPower(diff, exponent: p)
                if term.isNaN { return FP.nan }
                terms.append(term)
            }
        }
        return Helpers.sum(&terms)
    }
    
    /// Computes variance (population or sample) under lock.
    ///
    /// Process:
    /// 1. Computes the arithmetic mean.
    /// 2. Computes the second raw moment about the mean (Σ (v − μ)^2).
    /// 3. Divides by `n` for `.population` or by `n − 1` for `.sample`.
    ///
    /// - Parameter type: `.population` (σ²) or `.sample` (s²).
    /// - Returns: The variance or `nil` if undefined (e.g., fewer than 2 items).
    /// - Thread safety: The caller must hold the lock.
    /// - Complexity: O(n).
    @inline(__always)
    internal func varianceLocked(type: VarianceType) -> FP? {
        guard self.count > 1 else {
            return self.logNil("varianceLocked requires at least two observations")
        }
        guard let m = self.arithmeticMeanLocked() else {
            return self.logNil("varianceLocked cannot compute mean for the current data")
        }
        guard let M2 = self.rawMomentLocked(order: 2, about: m) else {
            return self.logNil("varianceLocked failed to compute the second raw moment about the mean")
        }
        if M2.isZero { return FP.zero }
        switch type {
        case .population:
            return M2 / FP(self.count)
        case .sample:
            return M2 / (FP(self.count) - FP.one)
        }
    }
    
    /// Computes standard deviation as the square root of the specified variance under lock.
    ///
    /// - Parameter type: `.population` or `.sample`.
    /// - Returns: The standard deviation or `nil` if variance is undefined.
    /// - Thread safety: The caller must hold the lock.
    /// - Complexity: O(n).
    @inline(__always)
    internal func standardDeviationLocked(type: VarianceType) -> FP? {
        guard let v = varianceLocked(type: type) else {
            return self.logNil("standardDeviationLocked unavailable because varianceLocked returned nil")
        }
        return v.squareRoot()
    }
    
    /// Computes semi-variance under lock.
    ///
    /// Definitions:
    /// - `.lower`: Mean of squared deviations for observations strictly below the mean.
    /// - `.upper`: Mean of squared deviations for observations strictly above the mean.
    /// - `.downside`: Mean of squared negative departures from the mean (positive deviations treated as zero).
    ///
    /// Edge cases:
    /// - For `.lower` and `.upper`, if no observations fall on the selected side, returns `0`.
    ///
    /// - Parameter type: The semi-variance definition to apply.
    /// - Returns: The semi-variance or `nil` if undefined (e.g., fewer than 2 items).
    /// - Thread safety: The caller must hold the lock.
    /// - Complexity: O(n).
    @inline(__always)
    internal func semiVarianceLocked(type: SemiVarianceType) -> FP? {
        guard self.count > 1 else {
            return self.logNil("semiVarianceLocked requires at least two observations")
        }
        guard let values: [FP] = self.itemsAsNumericArrayLocked() else {
            return self.logNil("semiVarianceLocked requires numeric data but itemsAsNumericArrayLocked returned nil")
        }
        guard let m = self.arithmeticMeanLocked() else {
            return self.logNil("semiVarianceLocked cannot compute the arithmetic mean")
        }
        var terms: [FP]
        var n: FP
        switch type {
        case .lower:
            terms = (values.filter { $0 < m }).map { ($0 - m ) }.map { $0 * $0 }
            if terms.count == 0 { return FP.zero }
        case .upper:
            terms = (values.filter { $0 > m }).map { ($0 - m ) }.map { $0 * $0 }
            if terms.count == 0 { return FP.zero }
        case .downside:
            terms = values.map { FP.minimum(0, $0 - m) }.map { $0 * $0 }
        }
        let sum = Helpers.sum(&terms)
        n = FP(terms.count)
        return sum / n
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
    internal func lowerMeanAbsoluteDeviationLocked(from: FP) -> FP? {
        guard !self.isEmpty, self.isNumeric else {
            return self.logNil("lowerMeanAbsoluteDeviationLocked requires a non-empty numeric dataset")
        }
        if let values = self.itemsAsNumericArrayLocked() {
            let lowerPart = values.filter({ $0 < from })
            var terms: [FP] = []
            terms.reserveCapacity(lowerPart.count)
            for x in lowerPart {
                terms.append((x - from).magnitude)
            }
            return Helpers.sum(&terms) / FP(terms.count)
        }
        else {
            return self.logNil("lowerMeanAbsoluteDeviationLocked requires a non-empty numeric dataset")
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
    internal func upperMeanAbsoluteDeviationLocked(from: FP) -> FP? {
        guard !self.isEmpty, self.isNumeric else {
            return self.logNil("upperMeanAbsoluteDeviationLocked requires a non-empty numeric dataset")
        }
        if let values = self.itemsAsNumericArrayLocked() {
            let upperPart = values.filter({ $0 > from })
            var terms: [FP] = []
            terms.reserveCapacity(upperPart.count)
            for x in upperPart {
                terms.append((x - from).magnitude)
            }
            return Helpers.sum(&terms) / FP(terms.count)
        }
        else {
            return self.logNil("lowerMeanAbsoluteDeviationLocked requires a non-empty numeric dataset")
        }
    }

    /// Computes an empirical moment of order `p` under lock.
    ///
    /// Types:
    /// - `.central`: Σ (xᵢ − μ)^p / n
    /// - `.origin`: Σ xᵢ^p / n
    /// - `.standardized`: Σ (xᵢ − μ)^p / (n · s^p), where `s` is the sample standard deviation.
    ///
    /// Notes:
    /// - Returns `nil` for `.standardized` when variance is zero (division by zero).
    ///
    /// - Parameters:
    ///   - p: Order of the moment (non-negative integer).
    ///   - kind: The empirical moment type.
    /// - Returns: The empirical moment or `nil` on failure.
    /// - Thread safety: The caller must hold the lock.
    /// - Complexity: O(n).
    @inline(__always)
    internal func empiricalMomentLocked(order p: UInt, kind t: EmpiricalMomentType) -> FP? {
        guard self.count > 0 else {
            return self.logNil("empiricalMomentLocked requires at least one observation")
        }
        let n = FP(self.count)
        switch t {
        case .central:
            guard let m = self.arithmeticMeanLocked() else {
                return self.logNil("empiricalMomentLocked(.central) cannot compute mean")
            }
            guard let mp = self.rawMomentLocked(order: p, about: m) else {
                return self.logNil("empiricalMomentLocked(.central) could not compute the raw moment about the mean")
            }
            return mp / n
        case .origin:
            guard let mp = self.rawMomentLocked(order: p, about: 0) else {
                return self.logNil("empiricalMomentLocked(.origin) could not compute the raw moment about zero")
            }
            return mp / n
        case .standardized:
            guard let m = self.arithmeticMeanLocked() else {
                return self.logNil("empiricalMomentLocked(.standardized) cannot compute mean")
            }
            guard let s = self.varianceLocked(type: .sample) else {
                return self.logNil("empiricalMomentLocked(.standardized) cannot compute sample variance")
            }
            guard let mp = self.rawMomentLocked(order: p, about: m) else {
                return self.logNil("empiricalMomentLocked(.standardized) could not compute the raw moment about the mean")
            }
            if s.isZero {
                return self.logNil("empiricalMomentLocked(.standardized) undefined because sample variance is zero")
            }
            let sd = s.squareRoot()
            return (mp / FP.pow(sd, FP(p))) / n
        }
    }
    
    /// Computes Moors' kurtosis (quantile-based kurtosis) under lock.
    ///
    /// Formula:
    /// - (Q7 − Q5 + Q3 − Q1) / (Q6 − Q2), using octiles at probabilities
    ///   0.125, 0.25, 0.375, 0.625, 0.75, 0.875.
    ///
    /// Requirements:
    /// - `levelOfMeasurement` must be `.ordinal`.
    ///
    /// - Returns: Moors' kurtosis or `nil` if quantiles are unavailable or the level is unsuitable.
    /// - Thread safety: The caller must hold the lock.
    /// - Complexity: Dominated by quantile computations (typically O(n log n)).
    @inline(__always)
    internal func moorsKurtosisLocked() -> FP? {
        guard self.levelOfMeasurement == .ordinal else {
            return self.logNil("moorsKurtosisLocked requires ordinal level of measurement")
        }
        do {
            if let q1 = try self.quantileLocked(q: 0.125),
               let q2 = try self.quantileLocked(q: 0.25),
               let q3 = try self.quantileLocked(q: 0.375),
               let q5 = try self.quantileLocked(q: 0.625),
               let q6 = try self.quantileLocked(q: 0.75),
               let q7 = try self.quantileLocked(q: 0.875) {
                let n = q7 - q5 + q3 - q1
                let d = q6 - q2
                return n / d
            }
            else {
                return self.logNil("moorsKurtosisLocked could not obtain the required quantiles")
            }
        }
        catch let error {
            return self.logNil("Quantile computation failed for moorsKurtosisLocked: \(error)")
        }
    }
    
    /// Computes Bowley's skewness (quantile-based skewness) under lock.
    ///
    /// Formula:
    /// - (Q3 + Q1 − 2·Q2) / (Q3 − Q1), using quartiles at probabilities
    ///   0.125, 0.25, 0.375 (consistent with the quantile accessor used here).
    ///
    /// Requirements:
    /// - `levelOfMeasurement` must be `.ordinal`.
    ///
    /// - Returns: Bowley's skewness or `nil` if quantiles are unavailable or the level is unsuitable.
    /// - Thread safety: The caller must hold the lock.
    /// - Complexity: Dominated by quantile computations (typically O(n log n)).
    @inline(__always)
    internal func bowleySkewnessLocked() -> FP? {
        guard self.levelOfMeasurement == .ordinal else {
            return self.logNil("bowleySkewnessLocked requires ordinal level of measurement")
        }
        do {
            if let q1 = try self.quantileLocked(q: 0.125),
               let q2 = try self.quantileLocked(q: 0.25),
               let q3 = try self.quantileLocked(q: 0.375) {
                let n = q3 + q1 - 2 * q2
                let d = q3 - q1
                return n / d
            }
            else {
                return self.logNil("bowleySkewnessLocked could not obtain the required quantiles")
            }
        }
        catch let error {
            return self.logNil("Quantile computation failed for bowleySkewnessLocked: \(error)")
        }
    }
    
    /// Gini coefficient computed on a locked snapshot of the numeric data.
    ///
    /// Formula:
    /// - G = (1 / (n² · mean)) Σ_{i=1..n} (2i − n − 1) x_{(i)}, with `x_{(i)}` ascending order statistics.
    ///
    /// Preconditions:
    /// - Requires at least two observations.
    /// - Rejects negative or non-finite values because the Lorenz-based derivation assumes non-negative mass.
    ///
    /// Returns:
    /// - `nil` when preconditions fail; otherwise the Gini coefficient in [0, 1].
    internal func giniCoefficientLocked() -> FP? {
        guard self.count >= 2 else {
            return self.logNil("Gini coefficient requires at least two observations")
        }
        if let values = self.itemsAsNumericArrayLocked() {
            let sorted = values.sorted(by: <)
            // Reject negatives early
            if sorted.contains(where: { $0 < FP.zero }) {
                return self.logNil("Negative values not allowed for Gini coefficient")
            }
            // All zeros → Gini = 0
            if sorted.allSatisfy({ $0.isZero }) { return FP.zero }
            // Guard mean
            guard let m = self.arithmeticMeanLocked() else {
                return self.logNil("Mean not available for Gini coefficient calculation")
            }
            if m.isZero {
                return self.logNil("Mean is zero for Gini coefficient calculation")
            }
            // Optional: reject non-finite values to avoid undefined results
            if sorted.contains(where: { !$0.isFinite }) {
                return self.logNil("Non-finite values not allowed for Gini coefficient")
            }
            let N: FP = FP(self.count)
            var terms: [FP] = []
            terms.reserveCapacity(self.count)
            for i in 1...self.count {
                let x = sorted[i - 1]
                // coeff = (2i − N − 1)
                let coeff = (FP.two * FP(i)) - N - FP.one
                terms.append(coeff * x)
            }
            let s = Helpers.sum(&terms)
            let p = FP.pow(N, FP.two) * m
            return s / p
        }
        else {
            return self.logNil("No numeric data available for Gini coefficient calculation")
        }
    }
    
    /// Normalized Gini coefficient `G_n = G · n / (n − 1)` to bound results in [0, 1].
    ///
    /// - Returns: Rescales the raw Gini so a perfectly unequal sample (all mass on one item)
    ///   maps to 1. Returns `nil` when the base Gini cannot be computed.
    internal func giniNormLocked() -> FP? {
        if let g = self.giniCoefficientLocked() {
            let N = FP(self.count)
            return g * N / (N - FP.one)
        }
        else {
            return self.logNil("giniNormLocked unavailable because the Gini coefficient could not be computed")
        }
    }
    /// Gini's mean difference 2 · mean · G · n / (n - 1)
    ///
    /// - Returns: Gini's mean difference or `nil` if G or mean not available
    internal func giniMeanDifferenceLocked() -> FP? {
        if let g = self.giniCoefficientLocked(), let m = self.arithmeticMeanLocked() {
            let N = FP(self.count)
            return g * m * FP.two * N / (N - FP.one)
        }
        else {
            return self.logNil("giniMeanDifferenceLocked unavailable because the Gini coefficient or arithmetic mean could not be computed")
        }
    }
    
    /// Concentration ratio CR\_k: sum of the k largest values.
    ///
    /// - Parameter k: Number of top observations to include (1 ≤ k ≤ n).
    /// - Returns: Σ_{i=1..k} x\_{(i)} where x\_{(i)} are order statistics sorted descending,
    ///   or `nil` when inputs are invalid or non-numeric.
    internal func crLocked(k: Int) throws -> FP? {
        guard k > 0, k <= self.count else {
            let message = "k must be a positive integer no larger than the number of observations"
            SSLog.statisticsError(message)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: message)
        }
        if let a = self.itemsAsNumericArrayLocked() {
            let values = a.sorted(by: >)
            var terms = Array<FP>(values[0...k])
            if terms.allSatisfy( { !$0.isNaN } ) {
                let total = Helpers.sum(&terms)
                return total
            }
            else {
               return self.logNil("Insuffiently many non-NaN values")
            }
        }
        else {
            return self.logNil("itemsAsNumericArrayLocked failed")
        }
    }
    
    /// Herfindahl–Hirschman index (HHI) computed on the sample.
    ///
    /// - Returns: Σ (xᵢ / Σ xⱼ)² for positive numeric values, or `nil` if values are invalid
    ///   or totals cannot be established.
    internal func hhiLocked() -> FP? {
        if let a = self.itemsAsNumericArrayLocked() {
            if let N = self.total {
                let n = FP(N)
                if a.allSatisfy( { !$0.isNaN && $0.isFinite && $0 > FP.zero } ) {
                    var terms: [FP] = []
                    terms.reserveCapacity(a.count)
                    for item in a {
                        terms.append(FP.pow(item / n, FP.two))
                    }
                    let s = Helpers.sum(&terms)
                    return s
                }
                else {
                    return self.logNil("HHI cannot be calculated when there are NaN, negative or infinite values")
                }
            }
            else {
                return self.logNil( "total is nil")
            }
        }
        else {
            return self.logNil( "HHI can only be calculated for numeric data.")
        }
    }
    
    /// Computes the Lorenz curve points (cumulative unit versus cumulative value shares).
    /// - Returns: A sequence of `(u, v)` pairs where `u` is the population share and
    ///   `v` is the cumulative share of the observed values.
    /// - Note: Requires a non-empty, numeric, non-negative dataset.
    internal var lorenzCurveLocked: [(u: FP, v: FP)]? {
        guard self.count > 0 else {
            return self.logNil("lorenzCurveLocked requires at least one observation")
        }
        guard let rawValues = self.itemsAsNumericArrayLocked() else {
            return self.logNil("Lorenz curve can only be computed for numeric values")
        }
        
        var values: [FP] = []
        values.reserveCapacity(rawValues.count)
        for value in rawValues {
            if value.isNaN {
                return self.logNil("Lorenz curve cannot be computed when values contain NaN")
            }
            if !value.isFinite {
                return self.logNil("Lorenz curve cannot be computed when values contain infinite entries")
            }
            if value < FP.zero {
                return self.logNil("Lorenz curve requires non-negative values")
            }
            values.append(value)
        }
        guard !values.isEmpty else {
            return self.logNil("Lorenz curve requires at least one numeric value")
        }
        
        var sumTerms = values
        let total = Helpers.sum(&sumTerms)
        if total.isNaN {
            return self.logNil("Lorenz curve cannot be computed because the total is undefined")
        }
        
        let sampleSize = FP(values.count)
        var points: [(u: FP, v: FP)] = []
        points.reserveCapacity(values.count + 1)
        points.append((u: FP.zero, v: FP.zero))
        
        if total.isZero {
            for index in 1...values.count {
                let share = FP(index) / sampleSize
                points.append((u: share, v: share))
            }
            return points
        }
        
        let sortedValues = values.sorted(by: <)
        var cumulative: FP = .zero
        let invTotal = FP.one / total
        for (index, element) in sortedValues.enumerated() {
            cumulative += element
            let u = FP(index + 1) / sampleSize
            let scaled = cumulative * invTotal
            let v = FP.minimum(FP.one, FP.maximum(FP.zero, scaled))
            points.append((u: u, v: v))
        }
        // Ensure the Lorenz curve terminates exactly at (1, 1)
        if !points.isEmpty {
            points[points.count - 1] = (u: FP.one, v: FP.one)
        }
        return points
    }
}
