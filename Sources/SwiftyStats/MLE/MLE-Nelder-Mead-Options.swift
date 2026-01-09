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

import SwiftyStatsPrelude

/// Design for multi-start candidate generation in θ-space (feasible parameter space).
public enum MultiStartDesign: Sendable {
    case lhs           // Latin Hypercube in per-parameter ranges
    case sobol         // Sobol-like low-discrepancy (implemented as Halton fallback)
    case random        // IID uniform in ranges
    case grid          // Coarse grid (small dimensions only)
}

/// Strategy for setting initial simplex steps per dimension in u-space.
public enum InitialSimplexStepStrategy<T: RealLike>: Sendable {
    /// Use the per-dimension step stored in ParamSpec.step (default).
    case spec
    /// Use the same absolute u-step for all dimensions.
    case absoluteU(T)
    /// Use a relative step in θ: Δθ = r * max(1, |θ₀|) and map to u via Δu ≈ Δθ / (dθ/du)|_{u₀}.
    case relativeTheta(T)
}

/// Options for Nelder–Mead and gradient-based optimization used in the MLE solver.
///
/// This type centralizes all user-tunable knobs for the solver, including:
/// - Selecting the optimization algorithm (Nelder–Mead, BFGS, or L-BFGS).
/// - Convergence tolerances for the chosen algorithm.
/// - Line search and finite-difference controls for gradient-based methods.
/// - Numerical Hessian controls used to approximate the covariance of the estimates.
/// - Multi-start and warm-start orchestration, and basic global candidate generation.
///
/// All tolerances and steps are expressed in the unconstrained u-coordinate system
/// (i.e., the internal optimization space before parameter constraints are applied).
/// Defaults are chosen to be robust for a broad range of problems; for ill-scaled
/// or very noisy likelihoods you may need to relax tolerances or increase `maxIter`.
public struct MLEOptimizationOpts<T: RealLike>: Sendable {
    // MARK: Common

    /// Primary solver to employ (e.g., Nelder–Mead, BFGS, or grid search).
    public var optimizer: OptimizerKind = .nelderMead
    /// Maximum objective evaluations/iterations allowed before aborting.
    public var maxIter: Int = 10_000

    /// Relative tolerance on log-likelihood change (applied to NLL internally).
    public var relTolLogLik: T = T(1e-8)

    // MARK: Nelder–Mead

    /// Simplex diameter tolerance in u-space; below this the simplex is considered converged.
    public var tolSimplexNM: T = T(1e-8)
    /// Absolute f-spread tolerance (legacy).
    public var tolFNM: T = T(1e-10)
    /// Relative f-spread tolerance: (fmax - fmin)/max(1, |fmin|) < relTolFSpread.
    public var relTolFSpreadNM: T = T(1e-10)
    /// Stall stopping: if relative improvement of best f is below `stallRelTol` for this many iterations, stop.
    public var stallIterationsNM: Int = 200
    /// Relative improvement threshold used together with `stallIterationsNM` to detect stagnation.
    public var stallRelTolNM: T = T(1e-10)

    /// Reflection coefficient α used to mirror the worst vertex.
    public var alpha: T = T.one
    /// Expansion coefficient γ controlling how aggressively the simplex expands along improving directions.
    public var gamma: T = T.two
    /// Contraction coefficient ρ used when the reflected point fails.
    public var rho: T = T.half
    /// Shrink coefficient σ that collapses the simplex towards the best vertex.
    public var sigma: T = T.half

    /// Strategy for per-dimension initial simplex step selection.
    public var initialStepStrategy: InitialSimplexStepStrategy<T> = .spec

    // MARK: Gradient-based (BFGS / L-BFGS / Newton)

    /// Infinity-norm tolerance on the gradient in unconstrained space.
    public var tolGrad: T = T(1e-6)
    /// Relative step-size tolerance for consecutive iterates in u-space.
    public var tolStep: T = T(1e-10)

    // Strong Wolfe line search parameters
    /// Armijo parameter c1 in (0, 1e-1].
    public var lineSearchC1: T = T(1e-4)
    /// Curvature parameter c2 in (c1, 1). Typical 0.9 for quasi-Newton.
    public var lineSearchC2: T = T(0.9)
    /// Backoff factor used only when bracketing fails to find sufficient decrease quickly.
    public var lineSearchTau: T = T(0.5)
    /// Maximum total iterations (function/gradient eval points) in line search including zoom.
    public var lineSearchMaxIter: Int = 50
    /// Maximum iterations for the zoom phase.
    public var lineSearchZoomMaxIter: Int = 50
    /// Lower/upper bounds for step length during line search.
    public var lineSearchStepMin: T = T(1e-20)
    /// Upper bound on line-search step length to guard against runaway extrapolation.
    public var lineSearchStepMax: T = T(1e20)

    /// Finite-difference perturbation used to approximate gradients when an analytic gradient is unavailable.
    public var gradStep: T = T(1e-6)
    /// Number of correction pairs retained in the limited-memory BFGS history.
    public var lbfgsMemory: Int = 10
    /// Enables a final damped Newton polish after convergence of the quasi-Newton phase.
    public var enableNewtonRefinement: Bool = true
    /// Initial damping factor for Newton steps to avoid overshooting in poorly conditioned regions.
    public var newtonInitialDamping: T = T(1e-3)

    // MARK: Soft trust-region/backoff controls

    /// When line search or Newton encounters invalid evaluations, multiply step by this factor and retry.
    /// Kept separate from lineSearchTau to allow independent tuning for invalid-only backoffs.
    public var invalidBackoffTau: T = T(0.5)
    /// Maximum number of consecutive invalid backoffs before giving up in line search.
    public var invalidBackoffMax: Int = 20
    /// Maximum damping to try in Newton refinement when encountering invalid evaluations.
    public var newtonMaxDamping: T = T(1e12)

    // MARK: Trust-region refinement (new)

    /// Initial LM-style damping λ for the local quadratic model.
    public var trInitialLambda: T = T(1e-2)
    /// Minimum λ clamp.
    public var trLambdaMin: T = T(1e-12)
    /// Maximum λ clamp.
    public var trLambdaMax: T = T(1e12)
    /// Initial trust-region radius Δ.
    public var trInitialRadius: T = T(1.0)
    /// Acceptance threshold: if ρ ≥ η2, accept and expand radius.
    public var trEta2: T = T(0.75)
    /// Acceptance threshold: if ρ ≥ η1, accept (keep radius).
    public var trEta1: T = T(0.25)
    /// Radius expansion factor when model is good.
    public var trGammaDec: T = T(2.0)
    /// Radius shrink factor when model is poor.
    public var trGammaInc: T = T(0.5)
    /// Max CG iterations factor: max iters = min(k, trCGMaxIterFactor * k)
    public var trCGMaxIterFactor: Int = 5
    /// Finite-difference step for Hessian–vector product (central difference on grad).
    public var trHvEps: T = T(1e-4)
    /// Max number of refinement iterations.
    public var trMaxRefineIters: Int = 5

    // MARK: Covariance

    /// When `true`, approximate the parameter covariance via a numerical Hessian around θ̂.
    public var computeCovariance: Bool = false
    /// Step size used for finite-difference Hessian construction.
    public var hessianStep: T = T(1e-4)
    /// Initial diagonal regularization added to Hessian before Cholesky (λ I). If 0, try without reg first.
    public var hessianRegInitial: T = T(0)
    /// Factor by which λ is multiplied on each failed Cholesky attempt.
    public var hessianRegFactor: T = T(10)
    /// Maximum λ for regularization attempts.
    public var hessianRegMax: T = T(1e12)

    // MARK: Multi-start / Warm-start / Global candidates

    /// Number of randomly generated starting points to attempt.
    public var multiStartCount: Int = 1
    /// Sampling design used when drawing starting points.
    public var multiStartDesign: MultiStartDesign = .lhs
    /// Optional user-specified θ guesses evaluated before other candidates.
    public var warmStartTheta: [T]? = nil
    /// Number of additional starts spawned mid-optimization when stagnation occurs.
    public var randomRestartCount: Int = 0
    /// Radius in u-space around the incumbent solution for jittering new starts.
    public var multiStartURadius: T = T(2)
    /// Log-scale jitter applied when sampling new starts to better cover scale differences.
    public var multiStartLogScale: T = T(2)
    /// Number of grid samples per dimension when `multiStartDesign == .grid`.
    public var gridPointsPerDim: Int = 3
    /// Enable Swift Concurrency powered parallel evaluation of multi-start candidates on supported platforms.
    public var enableParallelStarts: Bool = true

    /// Optional RNG seed used for candidate generation and jitter.
    public var rngSeed: UInt64? = nil

#if swift(>=5.5)
    /// Task priority used for `parallelLocalRuns` when spawning detached tasks.
    /// Defaults to `.userInitiated` so solver work inherits the caller’s QoS.
    public var parallelTaskPriority: TaskPriority? = .userInitiated
#endif

    /// Maximum resampling attempts per start when initial candidate maps to invalid or excessively penalized NLL.
    /// Set to 0 to disable resampling.
    public var maxStartResampleAttempts: Int = 10

    /// Optional NLL cutoff used to prune “bad” starts: if NLL(u) > cutoff, discard and resample.
    /// If nil, only +∞ NLL is considered invalid.
    public var startNLLCutoff: T? = nil

    // MARK: Diagnostics / penalties

    /// Emits solver progress and summary diagnostics when `true`.
    public var diagnosticsEnabled: Bool = true

    /// θ-distance tolerance for considering two terminal solutions “the same” (legacy).
    public var uniquenessTolTheta: T = T(1e-6)
    /// Relative u-space uniqueness tolerance: ||Δu|| / max(1, ||u||, ||v||) < uniquenessRelTolU.
    public var uniquenessRelTolU: T = T(1e-6)

    /// Optional finite penalty added to NLL when logpdf returns invalid (NaN/±Inf).
    /// If nil, invalid evaluations yield +∞ (default behavior).
    public var invalidLogPdfPenalty: T? = nil

    /// Enable a smooth, scale-aware boundary proximity penalty in u-space to gently push back from θ-constraints.
    /// Set to 0 to disable. Typical small values: 1e-4 … 1e-1 depending on scaling.
    public var boundaryPenaltyWeight: T = T(0)
    /// Small epsilon to avoid singularities in boundary penalties.
    public var boundaryPenaltyEps: T = T(1e-12)
    /// Shape control for boundary penalties (power > 0). 1 = linear-like, 2 = quadratic-like.
    public var boundaryPenaltyPower: T = T(1)

    // MARK: Quasi-Newton safeguards

    /// Curvature threshold ε for accepting BFGS/L-BFGS updates: require sᵀy > ε ||s|| ||y||.
    public var curvatureEps: T = T(1e-10)
    /// Powell damping parameter θ in (0,1]; when curvature is weak or negative, use ȳ = θ y + (1−θ) H s.
    public var powellTheta: T = T(0.2)
    /// Floor/Ceiling for initial inverse-Hessian scalar γ used in L-BFGS two-loop recursion.
    public var h0ScaleFloor: T = T(1e-6)
    /// Upper cap on γ to avoid over-scaling the identity preconditioner.
    public var h0ScaleCeil: T = T(1e6)

    // MARK: New θ-scaled tolerances and windowed stall

    /// Scale-aware gradient infinity-norm threshold in θ scaling:
    /// max_i |g_u,i| / max(1, |θ_i|) < tolGradInfTheta
    public var tolGradInfTheta: T = T(1e-6)

    /// θ-space step tolerance (relative): ||Δθ||2 / max(1, ||θ||2) < tolStepTheta
    public var tolStepTheta: T = T(1e-8)

    /// Window size for recent relative NLL improvements to detect stalls (median-based).
    public var stallWindow: Int = 50

    /// Median threshold over the window to declare a stall.
    public var stallMedianTol: T = T(1e-10)

    /// Creates an option set populated with the recommended defaults.
    public init() {}
}
