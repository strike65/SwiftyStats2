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

/// Moment calculations (raw, central, standardized) over SSExamine observations.

import SwiftyStatsPrelude
import SwiftyBoost

extension SSExamine {
    /// The sum of all numeric items, weighted by their absolute frequencies.
    ///
    /// Definition:
    /// - For numeric datasets, computes Σ (value × frequency) over distinct items.
    ///
    /// Returns:
    /// - For numeric datasets: a finite value, ±∞, or `NaN` depending on the inputs.
    /// - For non-numeric datasets: `nil`.
    ///
    /// Caching and stability:
    /// - Maintains an internal cache (`totalStorage`) that is updated under lock on access.
    /// - Uses a magnitude-sorted, compensated summation strategy to reduce round-off error.
    ///
    /// Thread-safety:
    /// - Access is synchronized via `withLock`.
    ///
    /// Complexity:
    /// - Amortized O(k log k) for k distinct items when recomputed; cached thereafter until data mutates.
    public var total: FP? {
       get {
            return withLock { () -> FP? in
                // Keep cached total consistent with current items.
                updateTotalLocked()
                return self.totalStorage
            }
        }
    }
    

    /// The sum of squares of all numeric items, weighted by their absolute frequencies.
    ///
    /// Equivalent to `poweredTotal(power: 2)` but implemented with a fast path that
    /// avoids `pow(_, 2)` and uses direct squaring for performance and precision.
    ///
    /// Returns:
    /// - `nil` for non-numeric datasets.
    /// - Otherwise, a finite value, ±∞, or `NaN` depending on the inputs.
    ///
    /// Thread-safety:
    /// - Uses `withLock` to ensure a consistent snapshot of `items`.
    ///
    /// Complexity:
    /// - O(k log k) for k distinct items.
    public var squareTotal: FP? {
        get {
            withLock {
                // Delegate to the locked powered sum with explicit power = 2.
                powerSumLocked(power: 2)
            }
        }
    }
    

    /// The sum of reciprocals of all numeric items, weighted by their absolute frequencies.
    ///
    /// Computes Σ (frequency / value) over distinct items.
    ///
    /// Edge cases:
    /// - If any value is zero, the corresponding term is ±∞ and is included in the sum.
    /// - If both +∞ and −∞ terms appear, the sum is undefined and yields `NaN`.
    ///
    /// Thread-safety:
    /// - Uses `withLock` to ensure a consistent snapshot of `items`.
    ///
    /// Complexity:
    /// - O(k log k) for k distinct items (due to stable summation).
    public var inverseTotal: FP? {
        get {
            withLock {
                inverseTotalLocked()
            }
        }
    }
    

    /// The sum of powers of all numeric items, weighted by their absolute frequencies.
    ///
    /// Computes Σ (value^power × frequency) over distinct items.
    ///
    /// Special cases:
    /// - `power == 0`: Returns the sample size (Σ 1 × frequency).
    /// - `power == 1`: Returns `total`.
    /// - `power == 2`: Uses a fast path (direct squaring) for performance and precision.
    ///
    /// Edge cases:
    /// - Non-integer powers of negative values produce `NaN` via `pow` and short-circuit the sum to `NaN`.
    /// - Zero values with negative powers produce ±∞; these are included in the sum.
    ///
    /// Returns:
    /// - `nil` for non-numeric datasets.
    /// - Otherwise, a finite value, ±∞, or `NaN` depending on the inputs.
    ///
    /// Thread-safety:
    /// - Uses `withLock` to ensure a consistent snapshot of `items`.
    ///
    /// Complexity:
    /// - O(k log k) for k distinct items.
    public func poweredTotal(power : FP) -> FP? {
        withLock {
            powerSumLocked(power: power)
        }
    }
    

    /// Sum of natural logarithms of all numeric values: Σ log(x_i).
    ///
    /// Returns:
    /// - `nil` for non-numeric datasets or if any value is negative/`NaN`.
    /// - `-∞` if any value is exactly zero.
    /// - For empty datasets, returns 0.
    ///
    /// Notes:
    /// - When finite, `exp(logSum)` equals the product of all values. This is useful for
    ///   geometric means and logarithmic product accumulation.
    ///
    /// Thread-safety:
    /// - Uses `withLock` to ensure a consistent snapshot of `items`.
    ///
    /// Complexity:
    /// - O(n log n) for n observations (via stable summation).
    public var logSum: FP? {
        withLock {
            return self.logSumLocked()
        }
    }

    /// Computes the product of all numeric values using the specified method.
    ///
    /// - Parameter method: Strategy to compute the product (`.naive`, `.logarithmic`, `.scaled`,
    ///   `.exponentMantissa`, `.compensated`).
    /// - Returns: `nil` for non-numeric datasets; otherwise a `ProductResult` with the product
    ///   and optional diagnostics (e.g., `logValue`, `overflowRisk`, `underflowRisk`, `exponent`, `mantissa`).
    ///
    /// Special values and edge cases:
    /// - Detects signed zero, `NaN`, and ±∞ and propagates them sensibly per method.
    /// - `.logarithmic` and `.scaled` reduce overflow/underflow risk; `.compensated` aims to reduce rounding error.
    ///
    /// Thread-safety:
    /// - Computed under `withLock` to avoid races over the underlying items.
    ///
    /// Complexity:
    /// - O(n) for n observations for all strategies (some do additional constant-time diagnostics).
    public func product(method: ProductMethod) -> ProductResult<FP>? {
        return withLock{ () -> ProductResult<FP>? in
            productLocked(method: method)
        }
    }

    /// Unnormalized moment Σ (v − x)^p for numeric datasets.
    ///
    /// - Parameters:
    ///   - p: Order of the moment (non-negative integer).
    ///   - x: Reference point (0 for raw moments, mean for central moments, etc.).
    /// - Returns: The unnormalized moment, or `nil` if not applicable.
    ///
    /// Notes:
    /// - For p = 2 and x = mean, this equals n × variance (population or sample depending on normalization).
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and delegates to `rawMomentLocked`.
    ///
    /// Complexity:
    /// - O(n log n) in worst case due to stable summation; conversions are O(n).
    public func unnormalizedMoment(order p: UInt, about x: FP) -> FP? {
        return withLock { () -> FP? in
            self.rawMomentLocked(order: p, about: x)
        }
    }
    

    /// Central moment of order `p`, computed as Σ((xᵢ − μ)^p) / n.
    ///
    /// - Parameter p: Moment order (non-negative integer).
    /// - Returns: The central moment, or `nil` if the dataset is empty or non-numeric.
    public func centralMoment(order p: UInt) -> FP? {
        return withLock { () -> FP? in
            return empiricalMomentLocked(order: p, kind: .central)
        }
    }
    

    /// Raw moment about the origin, Σ(xᵢ^p) / n.
    ///
    /// - Parameter p: Moment order (non-negative integer).
    /// - Returns: The raw moment, or `nil` if the dataset is empty or non-numeric.
    public func originMoment(order p: UInt) -> FP? {
        return withLock { () -> FP? in
            return empiricalMomentLocked(order: p, kind: .origin)
        }
    }
    

    /// Standardized moment Σ((xᵢ − μ)^p) / (n·s^p) where `s` is the sample standard deviation.
    ///
    /// - Parameter p: Moment order (non-negative integer).
    /// - Returns: The standardized moment, or `nil` if the dataset is empty, non-numeric,
    ///   or if variance is zero (which would make the denominator undefined).
    public func standardizedMoment(order p: UInt) -> FP? {
        return withLock { () -> FP? in
            return empiricalMomentLocked(order: p, kind: .standardized)
        }
    }
    

    /// Raw moment Σ((xᵢ − x)^p) without normalization.
    ///
    /// - Parameters:
    ///   - p: Moment order (non-negative integer).
    ///   - x: Reference point (e.g., 0 for origin moments, μ for central moments).
    /// - Returns: The unnormalized sum, or `nil` if the dataset is empty or non-numeric.
    public func rawMoment(order p: UInt, about x: FP) -> FP? {
        return withLock { () -> FP? in
            return rawMomentLocked(order: p, about: x)
        }
    }
    
}
