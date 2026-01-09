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

/// Internal accumulators and tally helpers backing SSExamine statistics.

import SwiftyStatsPrelude
import SwiftyBoost

// MARK: Accumulators

/// Internal additive/multiplicative accumulators that run under the SSExamine lock.
///
/// These helpers implement stable summation (sorted Kahan), compensated products, and
/// logarithmic accumulation to preserve sign, exponent, and mantissa of products across
/// extreme magnitudes. They form the numerical core for totals, powered totals, and
/// product diagnostics without reacquiring locks.
extension SSExamine {
    /// Core product implementation under lock.
    ///
    /// - Applies method-specific strategies and returns `ProductResult` with diagnostics.
    /// - Handles special values (NaN, ±∞, ±0) consistently across strategies.
    /// - Thread safety: The caller must hold the lock.
    @inline(__always)
    internal func productLocked(method: ProductMethod) -> ProductResult<FP>? {
        guard self.isNumeric else { return nil } // Important branch: only defined for numeric datasets.
        if self.isEmpty { return ProductResult(value: 1) } // Important branch: empty product is 1.
        guard let values: [FP] = self.itemsAsNumericArrayLocked() else { return nil}
        
        // Count sign (including signed zeros) once for all methods.
        let negativeSignCount = values.reduce(into: 0) { acc, v in
            if v.sign == .minus { acc &+= 1 }
        }
        let sign: FP = negativeSignCount.isMultiple(of: 2) ? 1.0 : -1.0
        
        // Short-circuit zeros with correct signed zero.
        if values.contains(where: { $0 == FP.zero }) {
            // Important branch: any zero in product yields signed zero (respecting overall sign).
            let zero = FP.zero
            return ProductResult(value: sign == 1 ? zero : -zero, exponent: 0, mantissa: 0)
        }
        
        switch method {
        case .naive:
            // Simple left-to-right multiplication; fastest but least robust to overflow/underflow.
            var product = FP.one
            for value in values {
                product *= value
            }
            let result = ProductResult(value: product, logValue: FP.log(product), overflowRisk: false, underflowRisk: false, exponent: nil, mantissa: nil)
            return result
            
        case .logarithmic:
            // Sum logs of magnitudes and apply sign separately.
            guard let logSumMag = self.logSumOfMagnitudes(values) else { return nil }
            
            // Risk assessment and clamping using representable bounds.
            let logMax = FP.log(FP.greatestFiniteMagnitude)
            let logLeastNormal = FP.log(FP.leastNormalMagnitude)
            let logLeastNonzero = FP.log(FP.leastNonzeroMagnitude)
            
            if logSumMag > logMax {
                // Important branch: definitely overflows to ±∞.
                return ProductResult(value: sign * FP.infinity, logValue: logSumMag, overflowRisk: true, underflowRisk: false, exponent: nil, mantissa: nil)
            }
            if logSumMag < logLeastNonzero {
                // Definitely underflows to (signed) zero.
                let zero = FP.zero
                return ProductResult(value: sign == 1 ? zero : -zero, logValue: logSumMag, overflowRisk: false, underflowRisk: true, exponent: nil, mantissa: nil)
            }
            // May be subnormal if below leastNormalMagnitude.
            let ofRisk = logSumMag > logMax
            let ufRisk = logSumMag < logLeastNormal // may be subnormal
            let value = sign * FP.exp(logSumMag)
            return ProductResult(value: value, logValue: logSumMag, overflowRisk: ofRisk, underflowRisk: ufRisk, exponent: nil, mantissa: nil)
            
        case .scaled, .exponentMantissa:
            // Use exponent/mantissa accumulation to mitigate overflow/underflow.
            return scaledProduct(values)
            
        case .compensated:
            // Pairwise error tracking in multiplicative domain.
            let corrected = compensatedProduct(values)
            return ProductResult(value: corrected)
        }
    }
    
    /// Scaled product via exponent/mantissa accumulation (frexp/ldexp).
    ///
    /// - Keeps the mantissa normalized during accumulation to reduce risk of overflow/underflow.
    /// - Tracks exponent sum separately and reconstructs at the end.
    /// - Returns diagnostic info in `ProductResult`.
    /// - Thread safety: Pure function over the provided values.
    @inline(__always)
    internal func scaledProduct(_ values: [FP]) -> ProductResult<FP>? {
        guard !values.isEmpty else { return ProductResult(value: 1.0, logValue: 0.0, overflowRisk: false, underflowRisk: false, exponent: 0, mantissa: 1.0) }
        var mantissa: FP = 1.0
        var exponentSum: Int = 0
        // Compute overall sign once.
        let negativeSignCount = values.reduce(into: 0) { acc, v in
            if v.sign == .minus { acc &+= 1 }
        }
        let sign: FP = negativeSignCount.isMultiple(of: 2) ? 1.0 : -1.0
        
        for value in values {
            let absValue = value.magnitude
            let frex = Helpers.frexp(absValue)
            let m = frex.mantissa
            let exp = frex.exponent
            mantissa *= m
            exponentSum &+= exp
            
            // Keep mantissa normalized in [0.5, 1) each step.
            let norm = Helpers.frexp(mantissa.magnitude)
            if norm.mantissa != 0 {
                mantissa = norm.mantissa
                exponentSum &+= norm.exponent
            }
        }
        
        // Bound checks (approximate): given mantissa ∈ [0.5, 1), compare by exponent.
        let maxExp = Helpers.frexp(FP.greatestFiniteMagnitude).exponent
        let minNormalExp = Helpers.frexp(FP.leastNormalMagnitude).exponent
        let minNonZeroExp = Helpers.frexp(FP.leastNonzeroMagnitude).exponent
        
        if exponentSum > maxExp {
            // Important branch: exponent beyond max -> overflow to ±∞.
            return ProductResult(value: sign * FP.infinity, logValue: nil, overflowRisk: true, underflowRisk: false, exponent: exponentSum, mantissa: mantissa)
        }
        if exponentSum < minNonZeroExp {
            // Definitely underflows to signed zero.
            let zero = FP.zero
            return ProductResult(value: sign == 1 ? zero : -zero, logValue: nil, overflowRisk: false, underflowRisk: true, exponent: exponentSum, mantissa: mantissa)
        }
        
        // Risk flags (subnormal if below normal threshold).
        let ofRisk = exponentSum > maxExp
        let ufRisk = exponentSum < minNormalExp // may be subnormal
        
        let result = sign * mantissa * FP.pow(2, FP(exponentSum))
        return ProductResult(value: result, logValue: nil, overflowRisk: ofRisk, underflowRisk: ufRisk, exponent: exponentSum, mantissa: mantissa)
    }
    
    /// Dekker/Veltkamp split into high/low parts for error-free transformations.
    ///
    /// - Returns: (hi, lo) such that hi+lo = x and hi carries the leading bits.
    /// - Thread safety: Pure function.
    @inline(__always)
    internal func split<T: RealLike>(_ x: T) -> (hi: T, lo: T) {
        let p = T.significandBitCount
        let s = (p + 1) / 2
        // c is a power-of-two-ish scaling factor to extract high bits.
        let c = T(sign: .plus, exponent: T.Exponent(s), significand: 1) + 1
        let t = c * x
        let hi = t - (t - x)
        let lo = x - hi
        return (hi, lo)
    }
    
    /// Error-free product of two floating-point numbers: returns (product, rounding error).
    ///
    /// - Parameters:
    ///   - a: First factor.
    ///   - b: Second factor.
    /// - Returns: The rounded product and the associated rounding error.
    /// - Thread safety: Pure function.
    @inline(__always)
    internal func twoProd<T: RealLike>(_ a: T, _ b: T) -> (product: T, error: T) {
        let p = a * b
        let (ah, al) = split(a)
        let (bh, bl) = split(b)
        let err = ((ah * bh - p) + ah * bl + al * bh) + al * bl
        return (p, err)
    }

    /// Fast two-sum variant assuming |a| >= |b|: returns (sum, rounding error).
    ///
    /// - Parameters:
    ///   - a: First addend (with magnitude ≥ |b|).
    ///   - b: Second addend.
    /// - Returns: The rounded sum and the associated rounding error.
    /// - Thread safety: Pure function.
    @inline(__always)
    internal func fastSum<T: RealLike>(_ a: T, _ b: T) -> (sum: T, error: T) {
        // precond: abs(a) >= abs(b)
        let s = a + b
        let z = s - a
        let err = b - z
        return (s, err)
    }
    
    /// Compensated product accumulation with normalization and special value handling.
    ///
    /// - Handles NaN, ±∞, and signed zeros.
    /// - Tracks sign separately and keeps a (hi, lo) pair to reduce multiplicative error.
    /// - Thread safety: Pure function over the provided values.
    @inline(__always)
    internal func compensatedProduct<T: RealLike>(_ xs: [T]) -> T {
        if xs.isEmpty { return 1 } // Important branch: empty product is 1.
        var signMinus = false
        var expSum: Int = 0
        var hi = T(1)
        var lo = T(0)
        var sawNaN = false
        var sawPosInf = false
        var sawNegInf = false
        var sawZero = false

        for x in xs {
            if x.isNaN { sawNaN = true; continue } // Important branch: defer NaN decision until after loop.
            if x == 0 {
                sawZero = true
                signMinus.toggleIf(x.sign == .minus)
                continue
            }
            if x.isInfinite {
                if x.sign == .minus { sawNegInf = true } else { sawPosInf = true }
                signMinus.toggleIf(x.sign == .minus)
                continue
            }
            if x.sign == .minus { signMinus.toggle() }

            let xa = x.magnitude
            let (m, e) = Helpers.frexp(xa)
            expSum &+= e

            let (p, pe) = twoProd(hi, m)
            let q = lo * m
            let (s, e1) = fastSum(p, pe + q)
            hi = s
            lo = e1

            // Normalize mantissa to keep magnitude in a stable range.
            let (mn, eOff) = Helpers.frexp(hi.magnitude)
            if mn != 0 {
                let scale = Helpers.ldexp(T(1), -eOff)
                hi = hi * scale
                lo = lo * scale
                expSum &+= eOff
            }
        }

        // Final special value resolution:
        if sawNaN { return .nan }
        if sawPosInf || sawNegInf {
            if sawZero { return .nan } // Important branch: 0 * inf is undefined -> NaN.
            let inf = T.infinity
            return signMinus ? -inf : inf
        }
        if sawZero {
            let zero = T(0)
            return signMinus ? -zero : zero
        }

        // Recompose from compensated mantissa and exponent.
        let mant = hi + lo
        var result = Helpers.ldexp(mant, expSum)
        if signMinus { result = -result }
        return result
    }
}

// Helpers (fileprivate)
fileprivate extension Bool {
    /// Conditionally toggles the Boolean value.
    ///
    /// - Parameter condition: If `true`, `self` is toggled; otherwise, no effect.
    /// - Example:
    ///   - `var b = false; b.toggleIf(true)  // b == true`
    ///   - `var b = false; b.toggleIf(false) // b == false`
    mutating func toggleIf(_ condition: Bool) {
        if condition { self.toggle() }
    }
}
