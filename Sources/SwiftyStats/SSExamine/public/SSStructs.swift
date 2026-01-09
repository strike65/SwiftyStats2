//
//  Created by VT on 27.10.25.
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

/// Public support structs such as frequency tables, quartiles, and product results.

import SwiftyStatsPrelude

/// A single entry of a (non-cumulative) frequency table.
///
/// Represents the absolute and relative frequency for a particular value
/// observed in a sample.
/// It corresponds to the empirical probability mass function value
/// `(x, f(x), f(x)/n)` used when summing expectations over discrete support.
///
/// - Generics:
///   - SSElement: The value type being counted. Must be `Hashable` and `Comparable`
///     to support stable ordering and lookups, and `Sendable` for concurrency safety.
///   - FP: The floating-point type used for relative frequencies
///     (e.g. `Double`, `Float`). Must conform to `Real`, `BinaryFloatingPoint`,
///     and `Sendable`.
public struct FrequencyTableItem<SSElement, FP>: CustomStringConvertible
where SSElement: Hashable & Comparable & Sendable,
      FP: RealLike
{
    /// Backing storage for the counted value.
    private var itemStorage: SSElement
    
    /// The value this entry summarizes.
    public var item: SSElement { itemStorage }
    
    /// The relative frequency of `item` in the sample, typically in the range [0, 1].
    ///
    /// Implementations that construct this type should ensure the value is derived as
    /// `FP(frequency) / FP(sampleSize)` when `sampleSize > 0`, and commonly `0` when
    /// `sampleSize == 0`.
    public var relativeFrequency: FP
    
    /// The absolute number of occurrences (a non-negative integer) of `item` in the sample.
    public var frequency: Int
    
    /// Creates a frequency table entry.
    ///
    /// - Parameters:
    ///   - item: The value being summarized.
    ///   - relativeFrequency: The relative frequency (proportion) of the item.
    ///   - absoluteFrequency: The absolute count of the item.
    public init(item: SSElement, relativeFrequency: FP, absoluteFrequency: Int) {
        self.itemStorage = item
        self.relativeFrequency = relativeFrequency
        self.frequency = absoluteFrequency
    }
    
    /// A human-readable description containing the value, absolute count, and relative frequency.
    public var description: String {
        "Value: \(String(describing: item)), abs: \(frequency), rel: \(niceNumber(relativeFrequency))"
    }
}

/// A single entry of a cumulative frequency table.
///
/// Represents cumulative absolute and relative frequencies up to and including
/// a particular value in a sorted order.
/// Each entry is a point on the empirical CDF `(x, F_n(x))`, enabling
/// piecewise-constant integration and quantile lookups.
///
/// - Generics:
///   - SSElement: The value type being summarized. Must be `Hashable`, `Comparable`,
///     and `Sendable`.
///   - FP: The floating-point type used for cumulative relative frequencies
///     (e.g. `Double`, `Float`). Must conform to `Real`, `BinaryFloatingPoint`,
///     and `Sendable`.
public struct CumulativeFrequencyTableItem<SSElement, FP>: CustomStringConvertible
where SSElement: Hashable & Comparable & Sendable,
      FP: RealLike
{
    /// Backing storage for the cumulative boundary value.
    private var itemStorage: SSElement
    /// Backing storage for the cumulative relative count.
    private var cumRelCountStorage: FP
    /// Backing storage for the cumulative absolute count.
    private var cumCountStorage: Int

    /// The item value this entry summarizes (the upper bound of the cumulative count).
    public var item: SSElement {
        get { itemStorage }
        set { itemStorage = newValue }
    }

    /// The cumulative absolute count up to and including `item`.
    ///
    /// This is the sum of absolute frequencies for all values less than or equal to `item`
    /// under the chosen ordering.
    public var cumulativeCount: Int {
        get { cumCountStorage }
        set { cumCountStorage = newValue }
    }

    /// The cumulative relative frequency up to and including `item`.
    ///
    /// Typically computed as `FP(cumulativeCount) / FP(sampleSize)` when `sampleSize > 0`.
    public var cumulativeRelativeCount: FP {
        get { cumRelCountStorage }
        set { cumRelCountStorage = newValue }
    }

    /// Creates a cumulative frequency table entry.
    ///
    /// - Parameters:
    ///   - item: The value defining the cumulative boundary.
    ///   - cumulativeCount: The cumulative absolute count up to and including `item`.
    ///   - cumulativeRelativeCount: The cumulative relative frequency up to and including `item`.
    public init(item: SSElement, cumulativeCount: Int, cumulativeRelativeCount: FP) {
        self.itemStorage = item
        self.cumCountStorage = cumulativeCount
        self.cumRelCountStorage = cumulativeRelativeCount
    }

    /// A human-readable description containing the value, cumulative absolute count,
    /// and cumulative relative frequency.
    public var description: String {
        "Value: \(String(describing: item)), cumAbs: \(cumulativeCount), cumRel: \(niceNumber(cumulativeRelativeCount))"
    }
}

/// Diagnostic result for product computations.
///
/// Encapsulates the computed product and optional diagnostics that depend on the
/// chosen multiplication strategy (e.g., logarithmic or exponent/mantissa scaling).
///
/// - Generics:
///   - T: The floating-point type of the result. Must conform to `RealLike`,
///        `BinaryFloatingPoint`, and `Sendable`.
///
/// Fields:
/// - `value`: The final product (may be finite, ±infinity, or NaN).
/// - `logValue`: Optional natural logarithm of the magnitude of the product, when available
///   (e.g., for logarithmic methods).
/// - `overflowRisk`: Heuristic flag indicating a risk of overflow for certain strategies.
/// - `underflowRisk`: Heuristic flag indicating a risk of underflow or subnormal results.
/// - `exponent`: Accumulated exponent used by scaled strategies (if applicable).
/// - `mantissa`: Accumulated mantissa used by scaled strategies (if applicable).
public struct ProductResult<T: RealLike> {
    /// The final product.
    public let value: T
    /// Optional log(product magnitude), if provided by the computation method.
    public let logValue: T?
    /// Heuristic overflow risk indicator (method-dependent).
    public let overflowRisk: Bool?
    /// Heuristic underflow risk indicator (method-dependent).
    public let underflowRisk: Bool?
    /// Accumulated exponent for scaled strategies (if provided).
    public let exponent: Int?
    /// Accumulated mantissa for scaled strategies (if provided).
    public let mantissa: T?

    /// Creates a `ProductResult`.
    ///
    /// - Parameters:
    ///   - value: The product value.
    ///   - logValue: Optional log(product magnitude).
    ///   - overflowRisk: Optional overflow risk flag (defaults to `false` when present).
    ///   - underflowRisk: Optional underflow risk flag (defaults to `false` when present).
    ///   - exponent: Optional accumulated exponent (scaled strategies).
    ///   - mantissa: Optional accumulated mantissa (scaled strategies).
    public init (value: T, logValue: T? = nil, overflowRisk: Bool = false, underflowRisk: Bool = false, exponent: Int? = nil, mantissa: T? = nil) {
        self.value = value
        self.logValue = logValue
        self.overflowRisk = overflowRisk
        self.underflowRisk = underflowRisk
        self.exponent = exponent
        self.mantissa = mantissa
    }
}

/// Quartiles (25th, 50th/median, and 75th percentile).
///
/// Useful for describing the spread of a distribution. This type is agnostic to
/// the specific quartile definition used (inclusive/exclusive, interpolation rules);
/// it simply holds values provided by the caller.
///
/// - Generics:
///   - T: The floating-point type (e.g., `Double`, `Float`). Must conform to `RealLike`,
///        `BinaryFloatingPoint`, and `Sendable`.
public struct Quartile<T: RealLike>: CustomStringConvertible {
    /// 25th percentile (first quartile).
    public let lower: T
    /// 50th percentile (median).
    public let median: T
    /// 75th percentile (third quartile).
    public let upper: T
    
    /// Creates a quartile container.
    ///
    /// - Parameters:
    ///   - q25: The 25th percentile (first quartile).
    ///   - q50: The 50th percentile (median).
    ///   - q75: The 75th percentile (third quartile).
    public init(q25: T, q50: T, q75: T) {
        self.lower = q25
        self.median = q50
        self.upper = q75
    }
    
    /// A multi-line human-readable description of quartile values.
    public var description: String {
        var descr = String()
        descr.append("QUARTILES\n")
        descr.append("*********\n")
        descr.append("q25: \(niceNumber(self.lower))\nq50: \(niceNumber(self.median))\nq75: \(niceNumber(self.upper))\n")
        return descr
    }
}

/// A two-sided confidence interval.
///
/// Represents an interval estimate for a parameter (e.g., a population mean),
/// typically constructed from sample statistics and an assumed distribution.
/// This type is agnostic to how the interval was computed (normal approximation,
/// Student’s t, bootstrap, etc.); it simply holds the bounds and derived width.
///
/// - Generics:
///   - T: The floating-point type (e.g., `Double`, `Float`). Must conform to `RealLike`,
///        `BinaryFloatingPoint`, and `Sendable`.
///
/// Semantics and invariants:
/// - No ordering is enforced between `lowerBound` and `upperBound`. If a producer
///   supplies `lowerBound > upperBound`, the `width` will be negative. Callers that
///   require non-negative widths should normalize or check ordering upstream.
/// - `width` is computed as `upperBound - lowerBound` at initialization time and
///   stored for convenience.
///
/// Example:
/// ```swift
/// let ci: ConfidenceInterval<Double> = try examine.studentTCI(alpha: 0.05)!
/// print(ci.lowerBound, ci.upperBound, ci.width)
/// ```
///
/// Related:
/// - See `SSExamine.normalCI(alpha:populationSD:)` and `SSExamine.studentTCI(alpha:)`
///   for interval construction helpers.
public struct ConfidenceInterval <T: RealLike>: CustomStringConvertible, Codable {
    /// Lower interval bound.
    public var lowerBound: T
    /// Upper interval bound.
    public var upperBound: T
    /// Interval width, computed as `upperBound - lowerBound`.
    ///
    /// Note: This value may be negative if `lowerBound > upperBound`.
    public var width: T
    
    /// Creates a confidence interval from lower and upper bounds.
    ///
    /// - Parameters:
    ///   - lowerBound: The lower bound of the interval.
    ///   - upperBound: The upper bound of the interval.
    ///
    /// Notes:
    /// - This initializer does not reorder the bounds; it records them as given.
    /// - If you need a normalized interval, reorder externally before initializing.
    public init(lowerBound: T, upperBound: T) {
        self.lowerBound = lowerBound
        self.upperBound = upperBound
        self.width = upperBound - lowerBound
    }
    
    /// A multi-line human-readable description of the interval.
    public var description: String {
        var descr = String()
        descr.append("\nCONFIDENCE INTERVAL\n")
        descr.append("*******************\n")
        descr.append("lower bound: \(niceNumber(self.lowerBound))\n")
        descr.append("upper bound: \(niceNumber(self.upperBound))\n")
        descr.append("width: \(niceNumber(self.width))\n")
        return descr
    }
}

/// Box-and-whisker summary used to render boxplots.
///
/// Values are computed by SSExamine.boxWhiskerLocked as follows:
/// - Hinges: q25 and q75 are the first and third quartiles; median is q50. Quartile
///   definitions depend on the underlying quantile method used by the dataset (default tests expect Hyndman–Fan Type 7).
/// - IQR: iqr = q75 − q25.
/// - Whiskers (Tukey): lWhiskerExtreme is the smallest observation ≥ q25 − 1.5·iqr;
///   uWhiskerExtreme is the largest observation ≤ q75 + 1.5·iqr. If no interior value exists on a side, the hinge itself remains the whisker extreme.
/// - Extremes (mild outliers): values strictly between the 1.5·iqr and 3·iqr fences
///   on either side, i.e. q75 + 1.5·iqr < x ≤ q75 + 3·iqr or q25 − 3·iqr ≤ x < q25 − 1.5·iqr.
/// - Outliers (far outliers): values beyond the 3·iqr fences, i.e. x > q75 + 3·iqr or x < q25 − 3·iqr.
/// - Notches: uNotch = median + 1.57 · iqr / √n and lNotch = median − 1.57 · iqr / √n,
///   where n is the sample size used in the computation.
///
/// Notes:
/// - NaN values are excluded before classification; infinite values are not explicitly filtered by the current implementation.
/// - extrema and outliers contain raw values (not indices) in unspecified order.
/// - All properties are optional to allow partial construction or nil-returning callers.
///
/// - Generics:
///   - T: The floating-point type (e.g., `Double`, `Float`). Must conform to `RealLike`,
///        `BinaryFloatingPoint`, and `Sendable`.
public struct BoxWhisker<T: RealLike>: CustomStringConvertible, Codable {
    /// Median (q50).
    public var median: T?
    /// Lower hinge (q25).
    public var q25: T?
    /// Upper hinge (q75).
    public var q75: T?
    /// Interquartile range (q75 − q25).
    public var iqr: T?
    /// Smallest observation ≥ q25 − 1.5·IQR (Tukey lower whisker).
    public var lWhiskerExtreme: T?
    /// Largest observation ≤ q75 + 1.5·IQR (Tukey upper whisker).
    public var uWhiskerExtreme: T?
    /// Mild outliers: values between the 1.5·IQR and 3·IQR fences on either side
    /// (q25 − 3·IQR ≤ x < q25 − 1.5·IQR or q75 + 1.5·IQR < x ≤ q75 + 3·IQR).
    public var fences: Array<T>?
    /// Far outliers: values beyond the 3·IQR fences
    /// (x < q25 − 3·IQR or x > q75 + 3·IQR).
    public var outliers: Array<T>?
    /// Upper notch: median + 1.57 · IQR / √n.
    public var uNotch: T?
    /// Lower notch: median − 1.57 · IQR / √n.
    public var lNotch: T?
    
    /// Returns a description
    public var description: String {
        get {
            var descr = String()
            if let m = self.median, let q25 = self.q25, let q75 = self.q75, let iqr = self.iqr, let e = self.fences, let o = self.outliers, let lw = self.lWhiskerExtreme, let uw = self.uWhiskerExtreme, let ln = self.lNotch, let un = self.uNotch {
                let fencesStr: String = e.map { niceNumber($0) }.joined(separator: ", ")
                let outliersStr: String = o.map { niceNumber($0) }.joined(separator: ", ")
                descr.append("extreme of lower whisker: \(niceNumber(lw))\n")
                descr.append("q25: \(niceNumber(q25))\n")
                descr.append("lower notch: \(niceNumber(ln))\n")
                descr.append("median: \(niceNumber(m))\n")
                descr.append("upper notch: \(niceNumber(un))\n")
                descr.append("q75: \(niceNumber(q75))\n")
                descr.append("extreme of upper whisker: \(niceNumber(uw))\n")
                descr.append("iqr: \(niceNumber(iqr))\n")
                descr.append("fences (mild outliers): [\(fencesStr)]\n")
                descr.append("outliers (far outliers): [\(outliersStr)]\n")
            }
            return descr
        }
    }
}

internal struct ParseError: LocalizedError {
    let token: String
    let index: Int
    let path: String
    var errorDescription: String? {
        "Unable to parse token '\(token)' at index \(index) in file \(path)"
    }
}
