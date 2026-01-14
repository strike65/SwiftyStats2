# Maximum Likelihood Estimation with the Transformed-Parameter Solver

This guide explains how to define and solve MLE problems with the transformed-parameter solver used inside SwiftyStats, how to choose optimizers and configure multi-start/warm-start, how to provide gradients, and how to interpret results. The solver types (`MLEProblem`, `MLESolver`, `ParamSpec`, `ParamConstraint`) are internal to SwiftyStats and are not exported. External users should use `MLEFitter` with `MLEOptimizationOpts` and `OptimizerKind` instead. Sections that use `MLEProblem` or `ParamSpec` only build inside the SwiftyStats package (or a fork), because those types and `SwiftyStatsPrelude` are internal.

Public API quick start:

```swift
import SwiftyStats

let data = [1.2, 0.9, 1.1, 1.3, 0.95]
let result = MLEFitter<Double>.fitGamma(data: data)
print(result.thetaHat, result.logLik)
```

---

## 1. Overview

- What this is
  - A generic, transform-based maximum likelihood solver that operates in an unconstrained parameter space u and maps smoothly to your constrained model parameters θ. This lets you fit models with positivity, interval, and other constraints without manual reparameterization.

- Core pieces used inside the SwiftyStats module (internal types)
  - MLEProblem<T>: Defines your data, a per-observation log-pdf logpdf(x, θ), an optional per-observation gradient gradlogpdf(x, θ), and parameter specifications (constraints + initial values).
  - MLEOptimizationOpts<T>: All solver options (optimizer choice, tolerances, line search, Hessian/covariance, multi-start/warm-start/diagnostics).
  - MLESolver.fit(_:options:): The entry point that returns MLEResult<T> with θ̂, log-likelihood, diagnostics, and (optionally) a covariance matrix.

- Supported numeric types
  - The solver is generic over T: RealLike (Float, Double, and Float80 on Intel). Prefer Double unless you have a reason to use others.

---

## 2. Concepts

- Transformed coordinates
  - The solver optimizes in u ∈ R^k. Each parameter θ_i is constrained via a smooth, invertible transform θ_i = f_i(u_i):
    - real: identity
    - positive: softplus(u)
    - unitInterval: logistic(u)
    - lowerBound(a): a + softplus(u)
    - upperBound(b): b − softplus(u)
    - interval(a, b): a + (b − a) * logistic(u)
  - The Jacobian diagonal dθ/du is available and used for proper gradient transformation and covariance mapping.

- Objective
  - The solver minimizes the negative log-likelihood (NLL). Your logpdf returns log f(x | θ); the solver sums over data and negates internally.

- Optimizers
  - Nelder–Mead: Derivative-free, robust, slower; good fallback.
  - BFGS: Quasi-Newton with dense inverse Hessian; fast when gradients are smooth. Can use analytic gradient if provided, else finite differences.
  - L-BFGS: Memory-limited BFGS; better for higher-dimensional problems.
  - Optional damped Newton refinement (using numerical Hessian) can polish a solution when gradients are small.

---

## 3. Defining an MLE Problem

- MLEProblem<T> requires:
  - data: [T] — your observations
  - logpdf: (x: T, theta: [T]) -> T — per-observation log density/mass
  - gradlogpdf (optional): (x: T, theta: [T]) -> [T] — per-observation gradient with respect to θ
  - paramSpecs: [ParamSpec<T>] — one per parameter, with constraint and initial value (and an initial step in u)

Note: The examples in this section require working inside the SwiftyStats package so that `SwiftyStatsPrelude` is available.

Example: A custom two-parameter model

```swift
import SwiftyStatsPrelude

// Suppose θ = [mu, sigma] with mu ∈ R, sigma > 0.
// We'll write a simple logpdf for a Normal-like model (for demo; use your own distribution code).

let data: [Double] = /* your observations */

func logpdfNormal(_ x: Double, _ theta: [Double]) -> Double {
    let mu = theta[0]
    let sigma = theta[1]
    // Guard against invalid σ (the transform should ensure positivity, but keep it safe)
    if sigma <= 0 || !sigma.isFinite { return -.infinity }
    let z = (x - mu) / sigma
    return -0.5 * (log(2 * .pi) + 2 * log(sigma) + z*z)
}

// Optional analytic gradient wrt θ = [mu, sigma]
func gradlogpdfNormal(_ x: Double, _ theta: [Double]) -> [Double] {
    let mu = theta[0]
    let sigma = theta[1]
    let z = (x - mu) / sigma
    // ∂/∂mu: (x - mu) / sigma^2
    let dmu = (x - mu) / (sigma * sigma)
    // ∂/∂sigma: -1/sigma + (x - mu)^2 / sigma^3
    let dsigma = -1.0 / sigma + (z*z) / sigma
    return [dmu, dsigma]
}

// Parameter specs: mu ∈ R, sigma > 0
let specs = [
    ParamSpec<Double>(.real,     initial: 0.0),   // mu
    ParamSpec<Double>(.positive, initial: 1.0)    // sigma
]

// Assemble the problem
let problem = MLEProblem<Double>(
    data: data,
    logpdf: logpdfNormal,
    gradlogpdf: gradlogpdfNormal, // or nil to use finite differences
    paramSpecs: specs
)
```
## 4. Using Distribution-Based Parameter Starts (Factory)

- For many common distributions, you can generate robust parameter specs from data:
   - MLEProblem.makeParamSpecs(for: SwiftyBoostContDist, data: [T], ctx: MLEStartContext<T> = .init()) -> [ParamSpec<T>]

- Supported families include beta, gamma, weibull, cauchy, logistic, extreme_value (Gumbel), students_t, fisher_f, inverse_gamma, skew_normal, skew_t, gev, gpd, log_logistic, generalized_gamma, generalized_normal, burr_xii, johnson_su, noncentral variants (nc_*), and discrete families (bernoulli, binomial, negative_binomial, geometric, hypergeometric, poisson, holtsmark).

- Some families require context (ctx):
   - Binomial: ctx.binomialN (known n)
   - Hypergeometric: ctx.hyperN (population size), ctx.hyperSampleSize (draws per observation)
   - GPD: ctx.gpdExcessData if you already computed exceedances y = x − u

Example: Poisson(λ)
```swift
let x: [Double] = [0,1,1,2,3,0,4,2,1,0]
let blank = MLEStartContext<Double>()
let poissonSpecs = problem.makeParamSpecs(for: .poisson, data: x, ctx: blank)
// poissonSpecs contains [.positive(initial: sampleMeanClamped)]
```

---

## 5. Running the Solver

- Call MLESolver.fit(problem, options:):
   - Returns MLEResult<T> with:
      - thetaHat: [T] — estimated parameters
      - logLik: T — maximized log-likelihood
      - iterations, converged, nEval — run diagnostics
      - cov: [[T]]? — optional covariance (see options)
      - allSolutions, uniqueSolutionCount (u-space), uniqueSolutionCountTheta (θ-space) — optional multi-start diagnostics
      - gradientNormAtOpt, hessianPositiveDefinite, conditionNumberEstimateHu (+ conditionSource), convergenceReason — optional diagnostics

Example: Basic Nelder–Mead fit
```swift
var opt = MLEOptimizationOpts<Double>()
opt.optimizer = .nelderMead
opt.maxIter = 10_000
opt.diagnosticsEnabled = true
opt.computeCovariance = true // if you want a covariance estimate (numerical Hessian)

let result = MLESolver<Double>.fit(problem, options: opt)
print("θ̂ =", result.thetaHat, "logLik =", result.logLik)
print("Converged:", result.converged, "iters:", result.iterations, "evals:", result.nEval)
if let cov = result.cov { print("Covariance:", cov) }
```

### Solver pipeline (under the hood)

- **Transforms & scaling:** Parameters live in constrained θ-space; the solver maps them to unconstrained u, optimizes there, then pushes gradients/steps back via `dθ/du`. Step tolerances and gradient infinity norms are evaluated after scaling so stopping criteria remain meaningful even when θ spans multiple orders of magnitude.
- **Multi-start orchestration:** Candidate θ₀ vectors come from LHS/Sobol/random/grid designs. Starts that immediately map to invalid or high-penalty regions are resampled (up to the configured cap) before entering the local optimizer. On Apple platforms with Swift Concurrency available, the solver fans out these starts via TaskGroups to reduce wall-clock time.
- **Local solvers:** Nelder–Mead uses θ-relative simplex steps plus stall detection in both f and θ. BFGS/L-BFGS rely on strong Wolfe line search with safeguarded interpolation, curvature checks, and Powell damping. Optional trust-region refinement (truncated CG on (H+λI)p = −g) polishes stationary points when the gradient norm is already tiny.
- **Diagnostics:** When requested, the solver estimates Gershgorin bounds for κ(H_u), checks Hessian definiteness, counts unique solutions in both u and θ (range-normalized), and records why the selected run converged (gradient, step, function, stall, or maxIter).

---

## 6. Choosing an Optimizer

- Nelder–Mead (default)
   - No gradients needed. Good when logpdf is rough or non-smooth. Slower in high dimensions.

- BFGS
   - Fast and accurate when gradients are smooth. Provide gradlogpdf to avoid finite difference noise and extra function evaluations.

- L-BFGS
   - Use for higher-dimensional parameter vectors to limit memory footprint. Set options.lbfgsMemory (e.g., 5–20).

Example: BFGS with line search
```swift
var opt = MLEOptimizationOpts<Double>()
opt.optimizer = .bfgs
opt.tolGrad = 1e-6
opt.relTolLogLik = 1e-8
opt.lineSearchC1 = 1e-4
opt.lineSearchTau = 0.5

let res = MLESolver<Double>.fit(problem, options: opt)
```

---

## 7. Multi-start, Warm-start, and Random Restarts

- Multi-start
   - options.multiStartCount: number of starting points (including warm-start if provided).
   - options.multiStartDesign: .lhs, .sobol, .random, or .grid (grid is useful only for small k).
   - options.multiStartLogScale: controls span for positive/bounded parameters (log-symmetric sampling around base).
   - options.gridPointsPerDim: for grid design.

- Warm-start
   - options.warmStartTheta: provide a θ vector that satisfies constraints. Used as the first candidate.

- Random restarts
   - options.randomRestartCount: add purely random starts (in addition to multiStartCount).

- Diagnostics
   - Enable options.diagnosticsEnabled to collect all terminal solutions and count uniqueness.

Example: Multi-start L-BFGS with warm-start
```swift
var opt = MLEOptimizationOpts<Double>()
opt.optimizer = .lbfgs
opt.multiStartCount = 8
opt.multiStartDesign = .lhs
opt.randomRestartCount = 2
opt.warmStartTheta = [0.1, 2.0] // must satisfy constraints
opt.diagnosticsEnabled = true

let res = MLESolver<Double>.fit(problem, options: opt)
print("Best θ̂:", res.thetaHat)
if let sols = res.allSolutions {
    print("All terminal solutions (θ̂, logLik):")
    for (th, ll) in sols { print(th, ll) }
    print("Unique solutions (tol):", res.uniqueSolutionCount ?? 0)
    print("Unique θ solutions (scaled):", res.uniqueSolutionCountTheta ?? 0)
}
if let kappa = res.conditionNumberEstimateHu {
    print("Condition estimate [\(res.conditionSource ?? "H_u")]:", kappa)
}
if let reason = res.convergenceReason {
    print("Convergence reason:", reason)
}
```
---

## 8. Covariance and Hessian Diagnostics

- Covariance
   - Set options.computeCovariance = true.
   - The solver computes a numerical Hessian of the NLL in u, inverts it to get Cov(u), then maps it to Cov(θ) via diag(dθ/du).

- Hessian PD check
   - When diagnostics are enabled, the solver checks whether the observed Hessian in u is positive definite at the selected optimum (a weak proxy for a local maximum).

- Step size
   - options.hessianStep controls finite difference spacing. If your likelihood is noisy, consider increasing this (e.g., 1e-3 to 1e-2).

---

## 9. Providing Gradients

- If you provide gradlogpdf(x, θ), the solver:
   - Accumulates gradients over data to build ∇θ logLik.
   - Converts to ∇u NLL using the chain rule: ∇u NLL = − diag(dθ/du) ⋅ ∇θ logLik.

- If you do not provide gradients:
   - The solver can still run BFGS/L-BFGS using finite differences in u (options.gradStep).
   - For rough/noisy objectives, consider Nelder–Mead.

---

## 10. Tuning and Stability Tips

- Start values
   - Good initial values are crucial. Use the factory method makeParamSpecs(for:data:ctx:) where possible.

- Scaling and constraints
   - Choose constraints that match the true domain to help the transform and optimizer (e.g., positive, interval).

- Tolerances
   - For ill-scaled models or noisy data, relax tolerances:
      - Increase options.tolGrad, options.relTolLogLik, and options.tolStep.
      - Increase options.maxIter.

- Line search
   - If line search often fails, reduce initial step (implicitly via BFGS direction quality), adjust line search constants (c1 or tau), or try Nelder–Mead.

- Numerical Hessian
   - If covariance fails (non-invertible Hessian), try:
      - Increasing hessianStep.
      - Switching optimizer to reach a cleaner local maximum.
      - Using more data or better parameterization.

---

## 11. Examples

### 11.1 Poisson(λ)


```swift
let x = [0.0,1,1,2,3,0,4,2,1,0]

// log pmf: log f(x|λ) = x*log λ - λ - log(x!)
func logpmfPoisson(_ x: Double, _ theta: [Double]) -> Double {
    let lambda = theta[0]
    if lambda <= 0 { return -.infinity }
    // crude log-factorial via Stirling for demo; you can use a better implementation
    let lfac = x <= 1 ? 0.0 : (x*log(x) - x + 0.5*log(2 * .pi * x))
    return x * log(lambda) - lambda - lfac
}

let specs = [ParamSpec<Double>(.positive, initial: max(x.reduce(0,+)/Double(x.count), 1e-9))]
let prob = MLEProblem<Double>(data: x, logpdf: logpmfPoisson, paramSpecs: specs)

var opt = MLEOptimizationOpts<Double>()
opt.optimizer = .nelderMead
opt.computeCovariance = true

let res = MLESolver<Double>.fit(prob, options: opt)
print("λ̂ =", res.thetaHat[0], "logLik =", res.logLik)```

### 11.2 Beta(α, β) via factory starts

```swift
let x = (0..<200).map { _ in Double.random(in: 0.05...0.95) } // fake beta-like data

// Use the factory to get robust start values and constraints
let dummy = MLEProblem<Double>(data: x, logpdf: { _,_ in 0 }, paramSpecs: [])
let specs = dummy.makeParamSpecs(for: .beta, data: x)

// Define the logpdf (per observation)
func logpdfBeta(_ x: Double, _ theta: [Double]) -> Double {
    let a = theta[0], b = theta[1]
    if a <= 0 || b <= 0 || x <= 0 || x >= 1 { return -.infinity }
    // log Beta density = (a-1)log x + (b-1)log(1-x) - log B(a,b)
    // B(a,b) = Γ(a)Γ(b)/Γ(a+b)
    let lgB = try! SwiftyBoost.SpecialFunctions.lnBeta(a, b) // or your own implementation
    return (a-1)*log(x) + (b-1)*log(1-x) - lgB
}

let problem = MLEProblem<Double>(data: x, logpdf: logpdfBeta, paramSpecs: specs)

var opt = MLEOptimizationOpts<Double>()
opt.optimizer = .lbfgs
opt.multiStartCount = 6
opt.diagnosticsEnabled = true
opt.computeCovariance = true

let res = MLESolver<Double>.fit(problem, options: opt)
print("α̂, β̂ =", res.thetaHat, "logLik =", res.logLik)
```
### 11.3 Binomial(n, p) with known n

```swift
let counts: [Double] = [3,2,4,1,0,2,5,3,4,1]
let n = 6 // known number of trials

// log pmf per observation: log C(n, x) + x log p + (n-x) log(1-p)
func logpmfBinomial(_ x: Double, _ theta: [Double]) -> Double {
    let p = theta[0]
    if p <= 0 || p >= 1 { return -.infinity }
    let xi = Int(x.rounded())
    // log C(n, x)
    let lC = logGamma(Double(n+1)) - logGamma(Double(xi+1)) - logGamma(Double(n-xi+1))
    return lC + Double(xi)*log(p) + Double(n-xi)*log(1-p)
}

// Factory spec uses ctx.binomialN to produce a good start for p
let ctx = MLEStartContext<Double>(gpdExcessData: nil, binomialN: n, hyperN: nil, hyperSampleSize: nil)
let dummy = MLEProblem<Double>(data: counts, logpdf: {_,_ in 0}, paramSpecs: [])
let specs = dummy.makeParamSpecs(for: .binomial, data: counts, ctx: ctx)

let prob = MLEProblem<Double>(data: counts, logpdf: logpmfBinomial, paramSpecs: specs)

var opt = MLEOptimizationOpts<Double>()
opt.optimizer = .bfgs
opt.computeCovariance = true

let res = MLESolver<Double>.fit(prob, options: opt)
print("p̂ =", res.thetaHat[0], "logLik =", res.logLik)
```
### 11.4 GPD(σ, ξ) with exceedances

```swift
let exceedances: [Double] = /* y = x - u, y >= 0 */

func logpdfGPD(_ y: Double, _ theta: [Double]) -> Double {
    let sigma = theta[0], xi = theta[1]
    if sigma <= 0 || y < 0 { return -.infinity }
    let t = 1 + xi * y / sigma
    if t <= 0 { return -.infinity }
    // log density: -log σ - (1/ξ + 1) * log(1 + ξ y / σ) for ξ != 0
    if abs(xi) > 1e-12 {
        return -log(sigma) - (1/xi + 1) * log(t)
    } else {
        // limit xi -> 0 gives exponential
        return -log(sigma) - y / sigma
    }
}

let ctx = MLEStartContext<Double>(gpdExcessData: exceedances)
let dummy = MLEProblem<Double>(data: exceedances, logpdf: {_,_ in 0}, paramSpecs: [])
let specs = dummy.makeParamSpecs(for: .gpd, data: exceedances, ctx: ctx)

let prob = MLEProblem<Double>(data: exceedances, logpdf: logpdfGPD, paramSpecs: specs)

var opt = MLEOptimizationOpts<Double>()
opt.optimizer = .lbfgs
opt.multiStartCount = 8
opt.computeCovariance = true

let res = MLESolver<Double>.fit(prob, options: opt)
print("σ̂, ξ̂ =", res.thetaHat, "logLik =", res.logLik)
```

---

### 11.5 A small regression-style custom model (per-observation log-likelihood)

```swift
// Logistic regression-like log-likelihood for binary y ∈ {0,1} with scalar feature x
struct Obs { let x: Double; let y: Double }
let data: [Obs] = /* ... */

// θ = [β0, β1] both real
let specs = [ParamSpec<Double>(.real, initial: 0.0),
             ParamSpec<Double>(.real, initial: 0.0)]

// We’ll encode the observation as a single Double by packing index or use a closure capturing data.
// Here, we ignore the input x passed to logpdf and use the problem.data index trick:
let problem = MLEProblem<Double>(
    data: (0..<data.count).map { Double($0) },
    logpdf: { idxD, theta in
        let i = Int(idxD)
        let obs = data[i]
        let z = theta[0] + theta[1] * obs.x
        // logistic
        let p = 1.0 / (1.0 + exp(-z))
        // log Bernoulli
        if p <= 0 || p >= 1 { return -.infinity }
        return obs.y * log(p) + (1 - obs.y) * log(1 - p)
    },
    // Optional gradient in θ:
    gradlogpdf: { idxD, theta in
        let i = Int(idxD)
        let obs = data[i]
        let z = theta[0] + theta[1] * obs.x
        let p = 1.0 / (1.0 + exp(-z))
        let diff = obs.y - p
        // ∂/∂β0 = y - p; ∂/∂β1 = (y - p) * x
        return [diff, diff * obs.x]
    },
    paramSpecs: specs
)

var opt = MLEOptimizationOpts<Double>()
opt.optimizer = .bfgs
opt.computeCovariance = true

let res = MLESolver<Double>.fit(problem, options: opt)
print("β̂ =", res.thetaHat, "logLik =", res.logLik)
```

---

## 12. Interpreting Results

- θ̂: The MLE parameter vector in model space.
- logLik: The maximized log-likelihood (higher is better).
- iterations, nEval: How long the local solver worked; nEval counts objective calls.
- converged: Whether the stopping criteria were met (tolerances on function change and/or gradient norm).
- cov: Optional covariance matrix of θ̂ (approximate; based on numerical Hessian in u and delta-method mapping).
- allSolutions plus `uniqueSolutionCount`/`uniqueSolutionCountTheta`: If diagnostics are enabled you receive every terminal θ̂ and counts for distinct solutions in u-space and range-aware θ-space.
- gradientNormAtOpt: Norm of ∇u NLL at the chosen optimum; small values indicate stationarity.
- hessianPositiveDefinite: Whether the observed Hessian in u at the optimum appears PD (Cholesky success). For a maximum, the Hessian of NLL should be PD, but this is a local numerical check.
- conditionNumberEstimateHu & conditionSource: Gershgorin-based estimate of κ(H_u) (or its inverse) plus the matrix used, useful for diagnosing ill-conditioning.
- convergenceReason: Why the local solver stopped (`gradient`, `step`, `function`, `stall`, or `maxIter`), improving reproducibility when comparing runs.

### 12.1 Validating fitted models with bootstrap KS tests

After fitting a distribution, you can validate the fit (heuristically) by bootstrapping the one-sample Kolmogorov-Smirnov statistic. SwiftyStats provides:

- `GoodnessOfFit<T>.KSTest.ksTestOneSample(data:testTarget:)` for built-in targets (fits internally and returns a bootstrap p-value).
- `GoodnessOfFit<T>.KSTest.bootstrapOneSampleKS(...)` for custom fitters and custom CDFs.

Example: use an MLE-derived Normal model inside a deterministic bootstrap harness.

```swift
import SwiftyStats
import SwiftyBoost

let x: [Double] = [-1.2, -0.8, -0.5, -0.2, 0.0, 0.1, 0.4, 0.7, 1.1, 1.4]

// Solve a two-parameter Normal model with the transformed solver.
func fitNormal(_ data: [Double]) -> [Double] {
    let specs = [
        ParamSpec<Double>(.real, initial: 0.0),
        ParamSpec<Double>(.positive, initial: 1.0)
    ]
    let problem = MLEProblem<Double>(
        data: data,
        logpdf: { x, theta in
            let mu = theta[0]
            let sigma = theta[1]
            let z = (x - mu) / sigma
            return -0.5 * (log(2 * .pi) + 2 * log(sigma) + z * z)
        },
        paramSpecs: specs
    )
    var opt = MLEOptimizationOpts<Double>()
    opt.optimizer = .bfgs
    let result = MLESolver<Double>.fit(problem, options: opt)
    return result.thetaHat
}

let ks = try GoodnessOfFit<Double>.KSTest.bootstrapOneSampleKS(
    data: x,
    targetDistribution: .normal,
    fit: fitNormal,
    cdfFrom: { theta in
        let dist = try? SwiftyBoost.Distribution.Normal<Double>(mean: theta[0], sd: theta[1])
        return { value in
            guard let d = dist else { return .nan }
            return (try? d.cdf(value)) ?? .nan
        }
    },
    sampler: { n, theta, seed in
        RNGSampler<Double>.randoms(
            n: n,
            dist: .gaussian(mean: theta[0], standardDeviation: theta[1]),
            seed: seed
        )
    },
    samplerSeed: { 0xD00D_BEEF },
    reestimateEachReplicate: false,
    B: 500
)

print("D =", ks.d, "p_bootstrap =", ks.pBootstrap)
```

---

## 13. Platform and Performance Notes

- Accelerate usage
   - On Apple platforms with Accelerate, BLAS/LAPACK is used for dot products, norms, GEMM/GEMV, outer products, LU inversion, and Cholesky checks (Double/Float).

- Portable fallback
   - On platforms without Accelerate, portable code paths are used (including a naive symmetric inversion).

- Type choice
   - Prefer Double for general use. Float may be faster but less stable; Float80 (Intel only) can be more precise but slower.

---

## 14. Common Pitfalls

- Domain violations
   - Your logpdf must return a finite log density/mass for valid inputs. If it returns NaN/±∞, the solver penalizes that region and may backtrack frequently.

- Inconsistent gradients
   - If you provide gradlogpdf, ensure it matches logpdf exactly. Mismatches can cause line search failures or divergence.

- Boundary issues
   - For unitInterval parameters, do not return logpdf at exactly 0 or 1; use open intervals and stable transforms (already handled by ParamSpec and transforms).

- Covariance non-invertibility
   - If the Hessian is singular or near-singular, covariance may be nil. Try better starting values, more data, weaker correlations, or adjust hessianStep.

---

## 15. Testing Examples (Swift Testing)

```swift
import Testing

@Suite("MLE Poisson MLE basic sanity")
struct PoissonMLETests {

    @Test
    func estimateLambda() throws {
        let x: [Double] = [0,1,1,2,3,0,4,2,1,0]

        func logpmfPoisson(_ x: Double, _ theta: [Double]) -> Double {
            let lambda = theta[0]
            if lambda <= 0 { return -.infinity }
            let lfac = x <= 1 ? 0.0 : (x*log(x) - x + 0.5*log(2 * .pi * x))
            return x * log(lambda) - lambda - lfac
        }

        let lambda0 = max(x.reduce(0,+)/Double(x.count), 1e-9)
        let specs = [ParamSpec<Double>(.positive, initial: lambda0)]

        let prob = MLEProblem<Double>(data: x, logpdf: logpmfPoisson, paramSpecs: specs)

        var opt = MLEOptimizationOpts<Double>()
        opt.optimizer = .nelderMead

        let res = MLESolver<Double>.fit(prob, options: opt)
        #expect(res.converged)
        #expect(res.thetaHat[0] > 0)
    }
}
```

---

## 16. API Summary (What You Will Typically Touch)

### External (public)

- MLEFitter<T>
   - fit... methods such as fitGamma, fitNormal, fitPoisson, fitBinomial

- MLEOptimizationOpts<T>
   - optimizer: .nelderMead | .bfgs | .lbfgs
   - maxIter, relTolLogLik, tolGrad, tolStep, lineSearchC1, lineSearchTau, gradStep, lbfgsMemory
   - enableNewtonRefinement, newtonInitialDamping
   - computeCovariance, hessianStep
   - multiStartCount, multiStartDesign, warmStartTheta, randomRestartCount, multiStartURadius, multiStartLogScale, gridPointsPerDim
   - enableParallelStarts (default true) to fan out multi-start solves with Swift Concurrency on supported Apple platforms
   - diagnosticsEnabled, uniquenessTolTheta

- OptimizerKind
   - .nelderMead | .bfgs | .lbfgs

- MLEResult<T>
   - thetaHat, logLik, iterations, converged, nEval, cov
   - allSolutions, uniqueSolutionCount, uniqueSolutionCountTheta, gradientNormAtOpt, hessianPositiveDefinite, conditionNumberEstimateHu/conditionSource, convergenceReason

### Internal (within SwiftyStats)

- MLEProblem<T>
   - init(data: [T], logpdf: (T, [T]) -> T, gradlogpdf: ((T, [T]) -> [T])? = nil, paramSpecs: [ParamSpec<T>])
   - makeParamSpecs(for: SwiftyBoostContDist, data: [T], ctx: MLEStartContext<T> = .init()) -> [ParamSpec<T>]

- ParamSpec<T>
   - init(_ constraint: ParamConstraint<T>, initial: T, step: T = 0.2)
   - Constraints: .real, .positive, .unitInterval, .lowerBound(a), .upperBound(b), .interval(a,b)

- MLESolver<T>
   - static func fit(_ problem: MLEProblem<T>, options: MLEOptimizationOpts<T> = .init()) -> MLEResult<T>

---

## 17. Final Recommendations

- Start simple with Nelder–Mead and factory starts to ensure a robust baseline.
- If performance is an issue and your logpdf is smooth, switch to BFGS or L-BFGS and provide gradlogpdf.
- Use multi-start for multimodal or tricky problems; enable diagnostics to understand solution diversity.
- When you need uncertainty estimates, enable computeCovariance and check the Hessian PD flag.
