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

/// Shared enums for SSExamine options such as sort orders and binning strategies.

import SwiftyStatsPrelude

/// Strategies for computing the product of a collection while preserving
/// numerical stability and performance.
public enum ProductMethod {
    case naive
    case logarithmic
    case scaled
    case exponentMantissa
    case compensated
}

/// Sort order directives applied to statistic data arrays.
///
/// Conforms to:
/// - `Codable` for persistence and interop
/// - `CustomStringConvertible` for readable logging and UI display
public enum SortOrderDataArray: Int, Codable, CustomStringConvertible {
    /// Ascending order.
    ///
    /// Example: `[1, 2, 3, 4]`
    case ascending
    /// Descending order.
    ///
    /// Example: `[4, 3, 2, 1]`
    case descending
    /// Original (raw) order as provided by the source.
    case original
    /// Undefined or not specified.
    case none

    /// Human-readable string describing the associated sort order.
    public var description: String {
        switch self {
        case .ascending:
            return "ascending"
        case .descending:
            return "descending"
        case .original:
            return "original"
        case .none:
            return "none"
        }
    }
}

/// Sort order applied to unique item sets.
///
/// This enum is tailored for contexts where only unique values are considered
/// (e.g., distinct categories or deduplicated keys).
///
/// Raw values:
/// - `ascending = 1`
/// - `descending = 2`
/// - `none = 0xff` (undefined)
public enum SortOrderUniqeItems: Int, CustomStringConvertible, Codable {
    /// Ascending order of unique items.
    case ascending = 1
    /// Descending order of unique items.
    case descending = 2
    /// Undefined/not determined.
    case none = 0xff
    
    /// Human-readable string describing the associated sort order.
    public var description: String {
        switch self {
        case .ascending:
            return "ascending"
        case .descending:
            return "descending"
        case .none:
            return "none"
        }
    }
}

/// Sort orders for frequency tables.
///
/// Determines whether a frequency table is sorted by frequency or value,
/// and in which direction. Useful for presenting histograms, categorical
/// counts, and similar summaries.
public enum SortOrderFrequencyTable: Int, Codable, CustomStringConvertible {
    /// Sorts by frequency ascending.
    ///
    /// Example: rarest items first.
    case frequencyAscending = 1
    /// Sorts by frequency descending.
    ///
    /// Example: most frequent items first.
    case frequencyDescending = 2
    /// Sorts by value ascending.
    ///
    /// Example: alphabetical or numeric ascending by the item key/value.
    case valueAscending = 3
    /// Sorts by value descending.
    ///
    /// Example: reverse alphabetical or numeric descending by the item key/value.
    case valueDescending = 4
    /// Undefined or not specified.
    case none = 0

    /// Human-readable string describing the associated sort order.
    public var description: String {
        switch self {
        case .frequencyAscending:
            return "ascending by frequency"
        case .frequencyDescending:
            return "descending by frequency"
        case .valueAscending:
            return "ascending by value"
        case .valueDescending:
            return "descending by value"
        case .none:
            return "none"
        }
    }
}

/// Output format options for cumulative frequency tables.
///
/// Controls whether each unique item appears once or is expanded to its
/// observed multiplicity in the cumulative output.
public enum CumulativeFrequencyTableFormat: Int, Codable {
    /// Each unique item appears as many times as observed.
    ///
    /// Example: for data `[a, a, b]` you may see `a, a, b` in sequence.
    case eachUniqueItem = 1
    /// Each unique item appears exactly once.
    ///
    /// Example: for data `[a, a, b]` you may see `a, b` in sequence.
    case eachItem = 2
}

/// Statistical measurement scales recognised by SwiftyStats analytics.
///
/// Use these to declare the nature of variables so that appropriate
/// statistical procedures and visualisations can be selected.
///
/// - nominal: Unordered categorical labels (e.g., colors, names).
/// - ordinal: Ordered categories without equal intervals (e.g., rankings).
/// - interval: Ordered with equal intervals but no true zero (e.g., Celsius).
/// - ratio: Interval scale with an absolute zero (e.g., length, mass).
public enum LevelOfMeasurement: String, Codable {
    /// Unordered categorical scale with mutually exclusive labels.
    case nominal = "nominal"
    /// Ordered categories without fixed spacing between ranks.
    case ordinal = "ordinal"
    /// Ordered scale with consistent intervals but no true zero.
    case interval = "interval"
    /// Ordered scale with absolute zero enabling ratio comparisons.
    case ratio = "ratio"
}

/// Variance estimators supported by SwiftyStats.
public enum VarianceType: String, Codable {
    case population = "population"
    case sample = "sample"
}

/// Semi-variance perspectives for risk/dispersion analyses.
public enum SemiVarianceType: String, Codable {
    case upper = "upper"
    case lower = "lower"
    case downside = "downside"
}

/// Empirical moment definitions relative to the distribution.
public enum EmpiricalMomentType: String, Codable {
    case central = "central"
    case origin = "origin"
    case standardized = "standardized"
}

/// Kurtosis classifications derived from excess kurtosis.
public enum KurtosisType: String, Codable {
    case platykurtic = "platykurtic"
    case mesokurtic = "mesokurtic"
    case leptokurtic = "leptokurtic"
}

/// Skewness direction inferred from third-order moment analysis.
public enum SkewnessType: String, Codable {
    case leftSkewed = "leftSkewed"
    case rightSkewed = "rightSkewed"
    case symmetric = "symmetric"
}

/// Hyndman–Fan quantile definitions (Types 1…9).
///
/// These correspond to the plotting-position formulas described in
/// Hyndman & Fan (1996), “Sample Quantiles in Statistical Packages”.
/// They are commonly exposed by statistical packages such as R (`type = 1...9`).
public enum HyndmanFanQuantileType: CustomStringConvertible {
    case type1      // {{0,0},{1,0}}  (Type 1)
    case type2      // Type 2
    case type3      // {{1/2,0},{0,0}} (Type 3)
    case type4      // {{0,0},{0,1}}  (Type 4)
    case type5      // {{1/2,0},{0,1}} (Type 5)
    case type6      // {{0,1},{0,1}}  (Type 6)
    case type7      // {{1,-1},{0,1}} (Type 7)
    case type8      // {{1/3,1/3},{0,1}} (Type 8)
    case type9      // {{3/8,1/4},{0,1}} (Type 9)
    
    /// Describes the quantile algorithm as a human-readable label.
    public var description: String {
        switch self {
        case .type1:
            return "Type 1 (inverse empirical CDF)"
        case .type2:
            return "Type 2 (median-unbiased ECDF)"
        case .type3:
            return "Type 3 (nearest order statistic)"
        case .type4:
            return "Type 4 (linear interpolation of the empirical CDF)"
        case .type5:
            return "Type 5 (Hazen)"
        case .type6:
            return "Type 6 (Weibull)"
        case .type7:
            return "Type 7 (R default)"
        case .type8:
            return "Type 8 (median-unbiased Normal)"
        case .type9:
            return "Type 9 (approximately unbiased Normal)"
        }
    }
}

/// Type of the Rosner test for outliers (ESD test)
public enum ESDTestType: Int, Codable, CustomStringConvertible {
    /// consider lower tail only
    case lowerTail
    /// consider upper tail only
    case upperTail
    /// consider both tails
    case bothTails
    /// Human-readable description of which distribution tail(s) will be probed.
    public var description: String {
        switch self {
        case .lowerTail:
            return "lower tail"
        case .upperTail:
            return "upper tail"
        case .bothTails:
            return "both tails"
        }
    }
}
// MARK: - Binning

/// Strategies for picking the number of bins used by histograms.
public enum BinningStrategy {
    case auto
    case count(Int)
    case sturges
    case scott
    case freedmanDiaconis
    case sqrt
    case rice
}

/// Ways to evaluate Shannon entropy given numeric or categorical data.
public enum ShannonMethod<FP: RealLike & BinaryFloatingPoint & Sendable & RuntimeNumeric> {
    case histogram
    case frequency
    case kde(bandwidth: FP? = nil)
    case millerMadow
    case discrete  // For SSElement directly
}
