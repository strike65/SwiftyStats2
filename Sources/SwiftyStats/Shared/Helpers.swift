//
//  Created by VT on 29.10.25.
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
/// Helpers and runtime utilities shared across SwiftyStats:
/// - A tolerant, locale-aware converter (RealConverter) for turning “any” values
///   into floating-point numbers.
/// - Lightweight runtime numeric checking.
/// - Convenience helpers for parity and fractional/integer-part operations.
///
/// Notes:
/// - The converter is designed to be strict about non-finite values (NaN/Inf) and
///   will throw if such values are encountered in textual or numeric input.
/// - For string inputs, common locale formats and separators (including thin/NNBSP spaces
///   and quote-like group separators) are recognized.
/// - On Intel architectures, Float80 is handled with care to preserve extended precision
///   where possible.
import SwiftyStatsPrelude
import Foundation


/// Format a floating-point value to a fixed number of significant figures for user-facing output.
@inlinable
internal func niceNumber<T: RealLike>(_ value: T, sigFigs: Int = 6) -> String {
    let clampedFigs = Swift.max(1, Swift.min(sigFigs, 15))
    guard value.isFinite else { return String(describing: value) }
    if value == .zero { return "0" }
    #if DEBUG
    return String(format: "%.*g", 15, Double(value))
    #else
    return String(format: "%.*g", clampedFigs, Double(value))
    #endif
}

/// Format an optional floating-point value to a fixed number of significant figures for user-facing output.
@inlinable
internal func niceNumber<T: RealLike>(_ value: T?, sigFigs: Int = 6) -> String {
    guard let unwrapped = value else { return "nil" }
    return niceNumber(unwrapped, sigFigs: sigFigs)
}

/// Namespace that hosts runtime helper utilities shared across SwiftyStats.
internal enum Helpers {
    /// Unwraps an Optional value using reflection, returning the wrapped value or the original.
    static func unwrapOptional(_ value: Any) -> Any {
        let m = Mirror(reflecting: value)
        if m.displayStyle == .optional, let wrapped = m.children.first {
            return wrapped.value
        }
        return value
    }
    
    /// Returns true if the value (or wrapped Optional) conforms to `RuntimeNumeric`.
    static func _isRuntimeNumeric(_ value: Any) -> Bool {
        let v = unwrapOptional(value)
        return v is any RuntimeNumeric
    }
    
    /// Public-facing runtime numeric check.
    static func isNumeric(_ value: Any) -> Bool {
        return _isRuntimeNumeric(value)
    }
    
    /// Returns true if the floating-point value is an integer (fractional part is zero).
    static func isInteger<FP: Real & BinaryFloatingPoint>(_ value: FP) -> Bool {
        return value.truncatingRemainder(dividingBy: 1).isZero
    }
    
    /// Clamps and converts a floating-point value to a fixed-width integer `I`.
    static func integerValue<FP: BinaryFloatingPoint, I: FixedWidthInteger>(_ x: FP) -> I {
        guard x.isFinite else { return .zero }
        let dx = Double(x)
        let minD = Double(I.min)
        let maxD = Double(I.max)
        if dx <= minD { return I.min }
        if dx >= maxD { return I.max }
        return I(dx)
    }
    
    /// Returns true if the value is an odd integer (non-integers return false).
    static func isOdd<FP: Real & BinaryFloatingPoint>(_ value: FP) -> Bool {
        let (ip, frac) = modf(value)
        if !frac.isZero {
            return false
        }
        return !ip.truncatingRemainder(dividingBy: 2).isZero
    }
    
    /// Returns true if the value is an even integer (non-integers return false).
    static func isEven<FP: Real & BinaryFloatingPoint>(_ value: FP) -> Bool {
        return !isOdd(value)
    }
    
    /// Returns the integer part of a floating-point value (truncation toward zero).
    static func integerPart<FP: Real & BinaryFloatingPoint>(_ value: FP) -> FP {
        let (ip, _) = modf(value)
        return ip
    }
    
    /// Returns the fractional part of a floating-point value.
    static func fractionalPart<FP: Real & BinaryFloatingPoint>(_ value: FP) -> FP {
        let (_, frac) = modf(value)
        return frac
    }
    
    /// Converts an arbitrary value to a `RealLike` floating-point type, throwing on failure.
    static func toReal<FP: RealLike>(_ value: Any, locale: Locale? = nil) throws -> FP {
        return try RealConverter.to(value, locale: locale)
    }
    
    
    /// Decomposes `x` into mantissa and exponent such that `x = mantissa * 2^exponent`,
    /// with `mantissa` in `[0.5, 1)` (or `(-1, -0.5]` for negatives).
    static func frexp<T: RealLike>(_ x: T) -> (mantissa: T, exponent: Int) {
        if x == 0 || !x.isFinite {
            return (x, 0)
        }
        
        let signMinus = x.sign == .minus
        var y = signMinus ? -x : x
        
        var expAdjust = 0
        
        if y.isSubnormal {
            let k = T.significandBitCount
            let twoToK = T(sign: .plus, exponent: T.Exponent(k), significand: 1)
            repeat {
                y *= twoToK
                expAdjust -= k
            } while y.isSubnormal
        }
        
        let e = y.exponent
        let twoToE = T(sign: .plus, exponent: e, significand: 1)
        var s = y / twoToE
        let outExp = Int(e) + 1 + expAdjust
        s *= 0.5
        
        if signMinus { s = -s }
        
        return (s, outExp)
    }
    
    /// Computes `m * 2^e` for `RealLike` floating-point types.
    static func ldexp<T: RealLike>(_ m: T, _ e: Int) -> T {
        m * T(sign: .plus, exponent: T.Exponent(e), significand: 1)
    }
    
    /// Numerically stable summation of a list of BinaryFloatingPoint terms.
    ///
    /// Strategy:
    /// - For 0 terms, returns 0.
    /// - For 1–2 terms, returns the direct sum (micro-optimization).
    /// - For 3+ terms, sorts by increasing magnitude and applies a Neumaier-style
    ///   compensated summation to reduce round-off error.
    ///
    /// - Parameter terms: The terms to sum. Callers should ensure no `NaN` values are present.
    ///                    Finite and infinite values are allowed.
    /// - Returns: The sum as `FP`. Returns `FP.nan` if the accumulation yields `NaN`
    ///            (e.g., mixing +∞ and −∞).
    ///
    /// Implementation note:
    /// - Uses `BinaryFloatingPoint.magnitude` instead of `abs(_)` for clarity and to avoid generic overloads.
    /// - Sorting is O(n log n). If you prefer speed over maximal precision for large n,
    ///   consider a single-pass Neumaier/Kahan without sorting.
    static func sum<T: RealLike>(_ terms: inout [T]) -> T {
        // Defensive: if empty, return additive identity (callers normally avoid empty)
        if terms.isEmpty { return T.zero }
        
        // Tiny-n fast path: avoid sort/compensation overhead for 1–2 terms
        if terms.count <= 2 {
            var s = T.zero
            for t in terms {
                s += t
            }
            // Important branch: if sum is NaN, propagate as FP.nan (e.g., +∞ + -∞).
            return s.isNaN ? T.nan : s
        }
        
        // General case: sort by increasing magnitude and use compensated summation
        terms.sort { $0.magnitude < $1.magnitude }
        var sum = T.zero
        var comp = T.zero
        for term in terms {
            // Allow finite or infinite terms; only NaN is rejected (already checked upstream)
            let t = sum + term
            // Neumaier compensation: accumulate the low-order error into `comp`.
            if sum.magnitude >= term.magnitude {
                comp += (sum - t) + term
            } else {
                comp += (term - t) + sum
            }
            sum = t
        }
        // Important branch: if either sum or compensation is NaN (e.g., undefined ∞ arithmetic), return NaN.
        if sum.isNaN || comp.isNaN {
            return T.nan
        } else {
            return T(sum + comp)
        }
    }
    
    static func softplus<T: RealLike>( _ x: T) -> T {
        let y:T = x > T.zero ? x : T.zero
        let t = -x.magnitude
        return y + T.log(onePlus: T.exp(t))
    }
    
    static func invSoftplus<T: RealLike>(_ theta: T) -> T {
        if theta > 20 {
            return theta + T.log(onePlus: -T.exp(-theta))
        } else {
            return T.log(T.expMinusOne(theta))
        }
    }
    
    static func logit<T: RealLike>( _ p: T) -> T {
        precondition(p.liesInOpenRange(from: 0, to: 1), "logit domain (0,1)")
        return T.log(p / (T.one - p ))
    }
    
    
    static func logistic<T: RealLike>( _ t: T) -> T {
        if t >= 0 {
            let e = T.exp(-t)
            return T.one / (T.one + e)
        }
        else {
            let e = T.exp(t)
            return e / (T.one + e)
        }
    }

}

