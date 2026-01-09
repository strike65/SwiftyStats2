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
/// Namespace enums and options for inferential statistics and hypothesis testing.

import SwiftyStatsPrelude

/// Root namespace for inferential statistics (hypothesis testing, goodness-of-fit, etc.).
///
/// Generic over `T: RealLike` so the same API surface can be specialized to `Double`, `Float`,
/// or platform-supported extended precision types. This is a static namespace and contains
/// sub-namespaces only.
public enum Inferential<T: RealLike> {
    /// Sub-namespace grouping hypothesis-testing helpers.
    ///
    /// Organizes routines by method family (parametric, non-parametric, variance tests, outliers,
    /// goodness-of-fit, and randomness). Each method typically returns a strongly-typed result
    /// structure that includes the test statistic, p-value, and method metadata.
    public enum HypothesisTesting {
        /// Marker namespace for outlier-detection routines (e.g., Grubbs, Rosner ESD).
        public enum OutliersTests {}
        /// Marker namespace for variance-comparison tests (e.g., Bartlett, Levene).
        public enum VarianceTests {}
        /// Marker namespace for goodness-of-fit tests (e.g., KS).
        public enum GoodnessOfFitTests{}
        /// Marker namespace for randomness and runs tests.
        public enum RandomnessTests{}
        /// Marker namespace for parametric hypothesis tests (e.g., t-tests, ANOVA, binomial).
        public enum ParametricTests{}
        /// Marker namespace for non-parametric tests (e.g., sign, Wilcoxon).
        public enum NonParametricTests{}

        public enum PostHocTests{}
    }
}

/// Shorthand for the parametric testing namespace specialized to floating-point type `T`.
///
/// Use, for example, `Parametric<Double>` to access Double-precision routines.
public typealias Parametric<T: RealLike> = Inferential<T>.HypothesisTesting.ParametricTests

/// Shorthand for the non-parametric testing namespace specialized to `T`.
public typealias NonParametric<T: RealLike> = Inferential<T>.HypothesisTesting.NonParametricTests

/// Shorthand for goodness-of-fit testing namespace specialized to `T`.
public typealias GoodnessOfFit<T: RealLike> = Inferential<T>.HypothesisTesting.GoodnessOfFitTests

/// Shorthand for variance-comparison testing namespace specialized to `T`.
public typealias Variance<T: RealLike> = Inferential<T>.HypothesisTesting.VarianceTests

/// Shorthand for outlier-detection namespace specialized to `T`.
public typealias Outliers<T: RealLike> = Inferential<T>.HypothesisTesting.OutliersTests

/// Shorthand for randomness and runs testing namespace specialized to `T`.
public typealias Randomness<T: RealLike> = Inferential<T>.HypothesisTesting.RandomnessTests

public typealias PostHoc<T: RealLike> = Inferential<T>.HypothesisTesting.PostHocTests

/// Supported target distributions for goodness-of-fit tests.
public enum GOFTestTarget: CustomStringConvertible {
    case arcsine
    case bernoulli
    case beta
    case binomial
    case cauchy
    case chisquared
    case exponential
    case gumbel
    case fisherF
    case gamma
    case geometric
    case holtsmark
    case inverseChiSquared
    case inverseGamma
    case inverseNormal
    case landau
    case laplace
    case logistic
    case lognormal
    case mapAiry
    case negativeBinomial
    case noncentralBeta
    case noncentralChiSquared
    case noncentralStudentT
    case noncentralFisherF
    case normalSTD
    case pareto
    case poisson
    case rayleigh
    case sasPoint5
    case skewNormal
    case StudentT
    case triangular
    case uniform
    case weibull
    case normal
//    case truncated(baseDistribution: any DistributionProtocol, lower: T, upper: T)
    
    /// User-facing display label for the distribution.
    public var description: String {
        switch self {
        case .arcsine:
            return "Arcsine"
        case .bernoulli:
            return "Bernoulli"
        case .beta:
            return "Beta"
        case .binomial:
            return "Binomial"
        case .cauchy:
            return "Cauchy"
        case .chisquared:
            return "Chi-squared"
        case .exponential:
            return "Exponential"
        case .gumbel:
            return "Gumbel (Extreme Value)"
        case .fisherF:
            return "Fisher F"
        case .gamma:
            return "Gamma"
        case .geometric:
            return "Geometric"
        case .holtsmark:
            return "Holtsmark"
        case .inverseChiSquared:
            return "Inverse Chi-squared"
        case .inverseGamma:
            return "Inverse Gamma"
        case .inverseNormal:
            return "Inverse Normal (Wald)"
        case .landau:
            return "Landau"
        case .laplace:
            return "Laplace"
        case .logistic:
            return "Logistic"
        case .lognormal:
            return "Log-normal"
        case .mapAiry:
            return "Airy"
        case .negativeBinomial:
            return "Negative Binomial"
        case .noncentralBeta:
            return "Noncentral Beta"
        case .noncentralChiSquared:
            return "Noncentral Chi-squared"
        case .noncentralStudentT:
            return "Noncentral Student's t"
        case .noncentralFisherF:
            return "Noncentral Fisher F"
        case .normalSTD:
            return "Std. Gaussian N(0,1)"
        case .normal:
            return "Normal N(mu, sigma)"
        case .pareto:
            return "Pareto"
        case .poisson:
            return "Poisson"
        case .rayleigh:
            return "Rayleigh"
        case .sasPoint5:
            return "SAS 0.5"
        case .skewNormal:
            return "Skew-Normal"
        case .StudentT:
            return "Student's t"
        case .triangular:
            return "Triangular"
        case .uniform:
            return "Uniform"
        case .weibull:
            return "Weibull"
//       case .truncated:
//            return "Truncated"
        }
    }
}

/// Alternative hypothesis for a test.
///
/// - twoSided: H₁: parameter ≠ reference
/// - greater: H₁: parameter > reference (upper-tailed)
/// - less:    H₁: parameter < reference (lower-tailed)
public enum Alternative: CustomStringConvertible, Codable {
    case twoSided
    case greater
    case less
    
    /// Human-readable label for the alternative hypothesis.
    public var description: String {
        switch self {
        case .twoSided: return "two-sided"
        case .greater: return "greater"
        case .less: return "less"
        }
    }
}

/// Tail selector used by some distribution helpers (lower/upper).
public enum Tail: Int, Codable, CustomStringConvertible {
    case lower, upper
    /// Human-readable label for the selected tail.
    public var description: String {
        switch self {
        case .lower: return "lower"
        case .upper: return "upper"
        }
    }
}
