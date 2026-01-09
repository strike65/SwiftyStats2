//
//  Created by VT on 09.11.25.
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
/// Problem definitions, convergence reasons, and results for MLE runs.

import SwiftyStatsPrelude

/// MLE problem definition: sample data, per-observation log-pdf, and parameter specs.
///
/// This struct encapsulates everything needed by the solver to maximize
/// the log-likelihood:
/// - data: The sample values.
/// - logpdf: A function that returns log f(x | θ) for a single observation x.
/// - gradlogpdf: Optional gradient of logpdf with respect to θ for a single observation.
/// - paramSpecs: Parameter constraints and initial values.
///
/// - Note: The solver operates in an internal unconstrained u-space and maps
///   to θ-space using smooth transforms derived from `paramSpecs`.
internal struct MLEProblem<T: RealLike>: Sendable {
    /// Sample data (must be non-empty).
    public let data: [T]
    /// logpdf(x, theta) -> log f(x | θ); θ is a vector [T].
    public let logpdf: @Sendable (_ x: T, _ theta: [T]) -> T
    /// Optional gradient of logpdf wrt θ for a single observation.
    /// If provided, the solver will use it to build the total gradient of the log-likelihood.
    public let gradlogpdf: (@Sendable (_ x: T, _ theta: [T]) -> [T])?
    /// Parameter constraints, initial values and steps.
    public let paramSpecs: [ParamSpec<T>]

    /// Create a new MLE problem instance.
    ///
    /// - Parameters:
    ///   - data: Sample data (must be non-empty).
    ///   - logpdf: Per-observation log density `log f(x | θ)`.
    ///   - gradlogpdf: Optional per-observation gradient `∂/∂θ log f(x | θ)`.
    ///   - paramSpecs: Parameter constraints and initial values.
    public init(
        data: [T],
        logpdf: @escaping @Sendable (_ x: T, _ theta: [T]) -> T,
        gradlogpdf: (@Sendable (_ x: T, _ theta: [T]) -> [T])? = nil,
        paramSpecs: [ParamSpec<T>]
    ) {
        precondition(!data.isEmpty, "data must not be empty")
        precondition(!paramSpecs.isEmpty, "need at least one parameter")
        self.data = data
        self.logpdf = logpdf
        self.gradlogpdf = gradlogpdf
        self.paramSpecs = paramSpecs
    }
}

/// Convergence reason reported by the local optimizer.
public enum ConvergenceReason: String, Codable, Sendable {
    case gradient
    case step
    case function
    case stall
    case maxIter
}

/// Result of an MLE fit via transformed coordinates.
///
/// Contains the final estimate, objective value, convergence information,
/// and optional covariance and diagnostics.
public struct MLEResult<T: RealLike>: Sendable {
    /// Estimated parameters in θ-space.
    public let thetaHat: [T]
    /// Maximized log-likelihood at `thetaHat`.
    public let logLik: T
    /// Number of iterations performed by the local optimizer for the selected solution.
    public let iterations: Int
    /// Whether the optimizer reported convergence under the configured criteria.
    public let converged: Bool
    /// Total number of objective evaluations across all starts.
    public let nEval: Int
    /// Optional covariance matrix (θ-space), obtained via the observed Hessian of the NLL.
    public let cov: [[T]]?
    /// Optional robust “sandwich” covariance (θ-space), when gradlogpdf is available:
    /// Cov(θ) ≈ H⁻¹ J H⁻¹ mapped to θ with Jacobian, where J is the sum of score outer products.
    public let robustCov: [[T]]?

    // Diagnostics (optional; populated when enabled in options)
    /// All terminal solutions across multi-starts, as `(θ̂, logLik)`.
    public let allSolutions: [([T], T)]?
    /// Count of distinct terminal solutions (in u-space) up to tolerance.
    public let uniqueSolutionCount: Int?
    /// Count of distinct terminal solutions (in θ-space) using range-aware scaling.
    public let uniqueSolutionCountTheta: Int?
    /// Norm of the gradient of NLL in u-space at the selected optimum.
    public let gradientNormAtOpt: T?
    /// Whether the observed Hessian (in u) was positive definite at the optimum.
    public let hessianPositiveDefinite: Bool?
    /// Estimated condition number κ(H_u) if available.
    public let conditionNumberEstimateHu: T?
    /// Source label for condition estimate (e.g., "H_u-Gershgorin" or "H_u-Power").
    public let conditionSource: String?
    /// Convergence reason for the selected local optimizer run, if available.
    public let convergenceReason: ConvergenceReason?

    /// Creates an `MLEResult` value from raw solver outputs.
    ///
    /// - Parameters:
    ///   - thetaHat: Final θ estimates.
    ///   - logLik: Maximized log-likelihood at `thetaHat`.
    ///   - iterations: Iteration count for the chosen solution.
    ///   - converged: Whether convergence tolerances were satisfied.
    ///   - nEval: Total objective evaluations across all starts.
    ///   - cov: Observed Hessian covariance estimate (optional).
    ///   - robustCov: Sandwich covariance estimate (optional).
    ///   - allSolutions: Terminal θ/log-likelihood pairs.
    ///   - uniqueSolutionCount: Unique terminal solutions in u-space.
    ///   - uniqueSolutionCountTheta: Unique terminal solutions in θ-space.
    ///   - gradientNormAtOpt: ||∇_u NLL|| at the optimizer.
    ///   - hessianPositiveDefinite: Whether Hessian was SPD.
    ///   - conditionNumberEstimateHu: Estimated κ(H_u).
    ///   - conditionSource: Source label for the condition estimate.
    ///   - convergenceReason: Reason code emitted by the optimizer.
    public init(thetaHat: [T],
                logLik: T,
                iterations: Int,
                converged: Bool,
                nEval: Int,
                cov: [[T]]?,
                robustCov: [[T]]? = nil,
                allSolutions: [([T], T)]? = nil,
                uniqueSolutionCount: Int? = nil,
                uniqueSolutionCountTheta: Int? = nil,
                gradientNormAtOpt: T? = nil,
                hessianPositiveDefinite: Bool? = nil,
                conditionNumberEstimateHu: T? = nil,
                conditionSource: String? = nil,
                convergenceReason: ConvergenceReason? = nil) {
        self.thetaHat = thetaHat
        self.logLik = logLik
        self.iterations = iterations
        self.converged = converged
        self.nEval = nEval
        self.cov = cov
        self.robustCov = robustCov
        self.allSolutions = allSolutions
        self.uniqueSolutionCount = uniqueSolutionCount
        self.uniqueSolutionCountTheta = uniqueSolutionCountTheta
        self.gradientNormAtOpt = gradientNormAtOpt
        self.hessianPositiveDefinite = hessianPositiveDefinite
        self.conditionNumberEstimateHu = conditionNumberEstimateHu
        self.conditionSource = conditionSource
        self.convergenceReason = convergenceReason
    }
}

/// Optional context for special cases and warm-start helpers (e.g., fixed dfs for noncentral distributions).
internal struct MLEStartContext<T: RealLike> {
    /// For GPD: exceedance data y = x - u (y ≥ 0). If nil, x are treated as exceedances directly.
    public var gpdExcessData: [T]? = nil

    /// For Binomial: known number of trials n (if available).
    public var binomialN: Int? = nil

    /// For Hypergeometric: known population size N (required for good starts).
    public var hyperN: Int? = nil
    /// For Hypergeometric: known sample size n (number drawn per observation).
    public var hyperSampleSize: Int? = nil

    /// For noncentral F: known degrees of freedom (df1, df2). We estimate only λ.
    public var ncF_df1: Int? = nil
    /// Optional denominator degrees of freedom for noncentral F.
    public var ncF_df2: Int? = nil

    /// For noncentral χ²: known degrees of freedom k. We estimate only λ.
    public var ncChi_k: Int? = nil

    /// For noncentral Student’s t: known ν. We estimate only δ.
    public var ncT_nu: Int? = nil

    /// For Log-Normal starts: optional overrides for μ (mean of log X) and σ (sd of log X).
    public var ln_mu: T? = nil
    /// Optional override for the log-normal log standard deviation.
    public var ln_sigma: T? = nil
    
    /// Cached minimum observation to avoid recomputation in start heuristics.
    public var minX: T? = nil
    /// Cached maximum observation to avoid recomputation in start heuristics.
    public var maxX: T? = nil

    /// Creates an empty context with no overrides.
    public init() {}
}
