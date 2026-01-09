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
/// Parameter specifications, constraints, and sample summaries for supported distributions.

import SwiftyStatsPrelude


/// A parameter-domain constraint that informs how unconstrained optimizer
/// coordinates (u-space) are mapped to valid model parameters (theta-space)
/// via smooth, invertible transforms.
///
/// The solver works in u-space. Each `ParamConstraint` selects an appropriate
/// transform (e.g., identity, softplus, logit, or affine-logit) so that
/// the mapped theta respects the distribution’s domain.
///
/// - Note: See MLE-transform.swift for the concrete transforms used.
///
/// - real: No constraint on theta (−inf, +inf).
/// - positive: theta > 0, typically using softplus.
/// - unitInterval: 0 < theta < 1, typically using logistic/logit.
/// - lowerBound(a): theta > a.
/// - upperBound(b): theta < b.
/// - interval(a,b): a < theta < b (requires a < b).
internal enum ParamConstraint<T: RealLike>: Sendable {
    case real                     // (-inf, inf)
    case positive                 // (0, inf)
    case unitInterval             // (0, 1)
    case lowerBound(T)            // (a, inf)
    case upperBound(T)            // (-inf, b)
    case interval(T, T)           // (a, b) with a < b
}

/// A parameter specification consisting of a domain constraint, an initial
/// value (in theta-space), and an initial step size used to seed the optimizer
/// in u-space.
///
/// The `step` is interpreted as an initial characteristic scale in u-space
/// for building starting simplices or line-search steps. `makeParamSpecs`
/// computes robust starts and calls `safeStep(_:)` to derive a reasonable
/// magnitude.
///
/// - Parameters:
///   - constraint: The domain constraint for theta.
///   - initial: Initial guess in theta-space.
///   - step: Initial step size in transformed coordinates (u). Defaults to 0.2.
internal struct ParamSpec<T: RealLike>: Sendable {
    /// Domain restriction applied to this parameter.
    public let constraint: ParamConstraint<T>
    /// Initial theta-space guess supplied to the solver.
    public let initial: T
    /// Initial step size in u-coordinates used when seeding simplices.
    public let step: T              // Initial step size in u-coordinates

    /// Creates a parameter specification with the provided domain, guess, and step.
    public init(_ c: ParamConstraint<T>, initial: T, step: T = T(0.2)) {
        self.constraint = c
        self.initial = initial
        self.step = step
    }
}


/// Supported distributions for which parameter specifications (constraints and
/// robust initial values) can be generated.
///
/// This enum is used by `MLEProblem.makeParamSpecs(for:data:ctx:)` to select
/// a factory tailored to the family, using simple method-of-moments or robust
/// estimators (median/MAD) and guarding edge cases.
///
/// The list covers continuous, discrete, and several noncentral families.
internal enum SwiftyBoostContDist {
    // Common families
    case beta                   // alpha>0, beta>0
    case gamma                  // k>0, theta>0
    case weibull                // k>0, lambda>0
    case cauchy                 // mu in R, gamma>0
    case logistic               // muinR, s>0
    case extreme_value          // (Gumbel) muinR, beta>0
    case students_t             // muinR, sigma>0, nu>0
    case fisher_f               // d1>0, d2>0  (real)
    case chiSquared
    case inverse_gamma          // alpha>0, beta>0
    case inverse_chi_squared
    case lognormal
    case inverseNormal

    // Skew / heavy tails
    case skew_normal            // xi in R, omega>0, alpha in R
    case skew_t                 // muinR, sigma>0, nu>0, alphainR

    // Extreme value / tail models
    case gev                    // muinR, sigma>0, xiinR
    case gpd                    // sigma>0, xiinR   (for exceedances y>=0)

    // Additional commonly used families
    case log_logistic           // alpha>0 (scale), beta>0 (shape)
    case generalized_gamma      // k>0, theta>0, p>0
    case generalized_normal     // (Exp-Power) muinR, alpha>0, beta>0
    case burr_xii               // c>0, k>0, s>0
    case johnson_su             // gamma in R, delta>0, xi in R, lambda>0

    // Noncentral families (always numerical)
    case nc_chi_squared         // k>0, lambda>0  (input data y>0)
    case nc_students_t          // mu in R, sigma>0, nu>0, delta in R
    case nc_fisher_f            // d1>0, d2>0, lambda>0
    case nc_beta                // alpha>0, beta>0, lambda>0
    
    case holtsmark
    
    case bernoulli
    case binomial
    case negative_binomial
    case geometric
    case hypergeometric
    case poisson
    case pareto
    case rayleigh
    case triangular
    case saspoint5
    case landau
    case mapairy

}


// ============================================================
// MARK: - Factory
// ============================================================
extension MLEProblem {
    /// Build parameter specifications (constraints, initial values, and u-steps)
    /// for a given distribution and dataset.
    ///
    /// This function selects a robust, distribution-appropriate initializer:
    /// - Filters input to the support (e.g., x>0 where required).
    /// - Uses robust or method-of-moments estimates (median/MAD, trimmed stats, or MoM).
    /// - Guards degeneracies (e.g., zero variance) via `safePositive` and clamps.
    ///
    /// - Parameters:
    ///   - dist: The target distribution family.
    ///   - x: The sample data. Must be non-empty.
    ///   - ctx: Optional start-context with hints (e.g., known dfs, bounds, or binomial N).
    /// - Returns: An ordered array of `ParamSpec` describing theta-parameters for `dist`.
    internal func makeParamSpecs(
        for dist: SwiftyBoostContDist,
        data x: [T],
        ctx: MLEStartContext<T> = .init()
    ) -> [ParamSpec<T>] {
        precondition(!x.isEmpty, "data must not be empty")
        
        switch dist {
        case .landau: return specs_Landau(x, ctx: ctx)
        case .saspoint5: return specs_SASPoint5(x, ctx: ctx)
        case .triangular: return specs_Triangular(x, ctx: ctx)
        case .inverseNormal: return specs_inorm(x)
        case .pareto: return specs_pareto(x, ctx: ctx)
            
        case .lognormal: return specs_LogNormal(x, ctx: ctx)
            // ---------- Beta(alpha,beta) ----------
        case .beta: return specs_beta(x)
            
            // ---------- Gamma(k,theta) ----------
        case .gamma: return specs_gamma(x)
            
            // ---------- Weibull(k,lambda) ----------
        case .weibull: return specs_weibull(x)
            
            // ---------- Cauchy(mu,γ) ----------
        case .cauchy: return specs_cauchy(x)
            
            // ---------- Logistic(mu,s) ----------
        case .logistic: return specs_logistic(x)
            
            // ---------- Extreme Value / Gumbel(mu,beta) ----------
        case .extreme_value: return specs_gumbel(x)
            
            // ---------- Student-t(mu,sigma,nu) ----------
        case .students_t: return specs_studentsT(x)
            
            // ---------- F(d1,d2) ----------
        case .fisher_f: return specs_f(x)
        case .chiSquared: return specs_ChSq(x, ctx: ctx)
            
            // ---------- Inverse-Gamma(alpha,beta) ----------
        case .inverse_gamma: return specs_invGamma(x)
            
            // ---------- Inverse-Chi-Squared(nu, τ²) ----------
        case .inverse_chi_squared: return specs_inverseChiSquared(x)
            
            // ---------- Skew-Normal(xi,ω,alpha) ----------
        case .skew_normal: return specs_skewNormal(x)
            
            // ---------- Skew-t(mu,sigma,nu,alpha) ----------
        case .skew_t: return specs_skewT(x)
            
            // ---------- GEV(mu,sigma,xi) ----------
        case .gev: return specs_gev(x)
            
            // ---------- GPD(sigma,xi), exceedances y = x-u (y>=0) ----------
        case .gpd:
            let y = ctx.gpdExcessData ?? x
            return specs_gpd(y)
            
            // ---------- Log-logistic(alpha,beta) ----------
        case .log_logistic: return specs_logLogistic(x)
            
            // ---------- Generalized Gamma(k,theta,p) ----------
        case .generalized_gamma: return specs_genGamma(x)
            
            // ---------- Generalized Normal / Exp-Power(mu,alpha,beta) ----------
        case .generalized_normal: return specs_genNormal(x)
            
            // ---------- Burr XII(c,k,s) ----------
        case .burr_xii: return specs_burrXII(x)
            
            // ---------- Johnson SU(γ,δ,xi,lambda) ----------
        case .johnson_su: return specs_johnsonSU(x)
            
            // ---------- Noncentral families ----------
        case .nc_chi_squared: return specs_ncChiSq(x, ctx: ctx)
        case .nc_students_t:  return specs_ncT(x, ctx: ctx)
        case .nc_fisher_f:    return specs_ncF(x, ctx: ctx)
        case .nc_beta:        return specs_ncBeta(x)

            // ---------- Discrete families ----------
        case .bernoulli:
            return specs_bernoulli(x)
        case .binomial:
            return specs_binomial(x, ctx: ctx)
        case .negative_binomial:
            return specs_negativeBinomial(x)
        case .geometric:
            return specs_geometric(x)
        case .hypergeometric:
            return specs_hypergeometric(x, ctx: ctx)
        case .poisson:
            return specs_poisson(x)
        case .holtsmark:
            return specs_holtsmark(x)
        case .rayleigh: return specs_rayleigh(x, ctx: ctx)
        case .mapairy: return specs_MapAiry(x, ctx: ctx)
        }
    }
}

// MARK: - Spec factories & helpers

private extension MLEProblem {
    
    /// Map-Airy location-scale start: robust center via median and scale via MAD.
    /// - Returns: [mu (real), sigma (positive)]
    func specs_MapAiry(_ x: [T], ctx: MLEStartContext<T>) -> [ParamSpec<T>] {
        let mu0 = sampleMedian(x, fallback: T.zero)
        let ex: SSExamine<T,T> = .init(usingArray: x, name: nil, characterSet: nil)
        let s0  = T.maximum(ex.medianAbsoluteDeviation(center: mu0)!, T(1e-3))
        return [
            ParamSpec(.real, initial: mu0, step: safeStep(mu0)),
            ParamSpec(.positive, initial: s0,  step: safeStep(s0))
        ]
    }
    
    /// SAS 0.5-point (Laplace-like) start: mean and standard deviation.
    /// - Returns: [mu (real), sigma (positive)]
    func specs_SASPoint5(_ x: [T], ctx: MLEStartContext<T>) -> [ParamSpec<T>] {
        let m = sampleMean(x, fallback: T.zero)
        let s = T.maximum(T.sqrt(sampleVariance(x, mean: m)), T(1e-6))
        return [
            ParamSpec(.real,      initial: m,  step: safeStep(m)),
            ParamSpec(.positive, initial: s,  step: safeStep(s))
        ]
    }
    
    /// Landau start: location via median, scale via MAD (mean/variance undefined).
    /// - Returns: [mu (real), sigma (positive)]
    func specs_Landau(_ x: [T], ctx: MLEStartContext<T>) -> [ParamSpec<T>] {
        // Landau has undefined mean/variance; use median-like start and a scale from MAD
        let mu0 = sampleMedian(x, fallback: T.zero)
        let ex: SSExamine<T,T> = .init(usingArray: x, name: nil, characterSet: nil)
        let s0  = T.maximum(ex.medianAbsoluteDeviation(center: mu0)!, T(1e-3))
        return [
            ParamSpec(.real,      initial: mu0, step: safeStep(mu0)),
            ParamSpec(.positive, initial: s0,  step: safeStep(s0))
        ]
    }
    
    /// Triangular(a,b,c) with known a=min(x), b=max(x), estimate mode c in (a,b).
    /// - Returns: [c (interval(a,b))]
    func specs_Triangular(_ x: [T], ctx: MLEStartContext<T>) -> ([ParamSpec<T>]) {
        precondition(!x.isEmpty)
        let minX = ctx.minX ?? (x.min() ?? T.zero)
        let maxX = ctx.maxX ?? (x.max() ?? T.one)
        let mid = sampleMean(x, fallback: (minX + maxX) / T(2))
        let mode0 = T.minimum(maxX - T(1e-6), T.maximum(minX + T(1e-6), mid))

        let lowerStep = safeStep(minX == T.zero ? T(1) : minX)
        let upperStep = safeStep(maxX == T.zero ? T(1) : maxX)

        return [
            ParamSpec(.real, initial: mode0, step: safeStep(mode0)),
            ParamSpec(.real, initial: minX, step: lowerStep),
            ParamSpec(.real, initial: maxX, step: upperStep)
        ]
    }

    
    // Inverse Gaussian (Wald) ParamSpecs
    // Parameters: mu > 0, lambda > 0
    // Starts:
    //   mu0     = mean(x)
    //   lambda0 = n / sum( (x_i - mu0)^2 / (mu0^2 * x_i) )  (MLE plug-in)
    // Fallbacks guard positivity and degeneracy.
    ///
    /// - Returns: [mu (positive), lambda (positive)]
    func specs_inorm(_ x: [T]) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        precondition(!data.isEmpty, "inverse normal requires x > 0")

        // mu start: sample mean on positive support
        let muRaw = sampleMean(data, fallback: T(1))
        let mu0 = safePositive(muRaw, fallback: T(1))

        // lambda start: plug-in MLE denominator
        var denom: T = .zero
        let mu2 = mu0 * mu0
        for v in data {
            let d = v - mu0
            denom += (d * d) / (mu2 * v)
        }
        let n = T(data.count)
        let lamRaw = denom > .zero ? n / denom : T(1)
        let lambda0 = safePositive(lamRaw, fallback: T(1))

        return [
            ParamSpec(.positive, initial: mu0,     step: safeStep(mu0)),
            ParamSpec(.positive, initial: lambda0, step: safeStep(lambda0))
        ]
    }

    /// Pareto-I with xm fixed to min(x): estimate alpha from log-ratios.
    /// - Returns: [alpha (positive)]
    func specs_pareto(_ x: [T], ctx: MLEStartContext<T>) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        precondition(!data.isEmpty)
        let xmHat = data.min()!
        // MLE: alpha_hat = n / sum(log(x/xmHat))
        var s = T.zero
        for v in data { s += T.log(v / xmHat) }
        let alphaHat = T(data.count) / T.maximum(s, T(1e-8))
        return [
            // xm fixed to min(data) in Pareto-I MLE; if you want xm free, reparam with logXm and box xmax.
            ParamSpec(.positive, initial: alphaHat, step: safeStep(alphaHat))
        ]
    }
    
    /// Beta(alpha,beta) start via method-of-moments using mean and variance on (0,1).
    /// - Returns: [alpha (positive), beta (positive)]
    func specs_beta(_ x: [T]) -> [ParamSpec<T>] {
        let data = clamp01Samples(x)
        let mean = sampleMean(data, fallback: T(0.5))
        let variance = max(sampleVariance(data, mean: mean), T(1e-3))
        let common = max(mean * (1 - mean) / variance - 1, T(0.5))
        let alpha = safePositive(mean * common, fallback: T(1))
        let beta = safePositive((1 - mean) * common, fallback: T(1))
        return [
            ParamSpec(.positive, initial: alpha, step: safeStep(alpha)),
            ParamSpec(.positive, initial: beta, step: safeStep(beta))
        ]
    }

    /// Gamma(k,theta) start via mean/variance MoM.
    /// - Returns: [k (positive), theta (positive)]
    func specs_gamma(_ x: [T]) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let mean = sampleMean(data, fallback: T(1))
        let variance = max(sampleVariance(data, mean: mean), mean * T(1e-3) + T(1e-3))
        let shape = safePositive(mean * mean / variance, fallback: T(1))
        let scale = safePositive(variance / max(mean, T(1e-3)), fallback: T(1))
        return [
            ParamSpec(.positive, initial: shape, step: safeStep(shape)),
            ParamSpec(.positive, initial: scale, step: safeStep(scale))
        ]
    }

    /// Weibull(k,lambda) start: scale from mean, shape set to a moderate value.
    /// - Returns: [k (positive), lambda (positive)]
    func specs_weibull(_ x: [T]) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let scale = safePositive(sampleMean(data, fallback: T(1)), fallback: T(1))
        let shape = safePositive(T(1.5), fallback: T(1))
        return [
            ParamSpec(.positive, initial: shape, step: T(0.5)),
            ParamSpec(.positive, initial: scale, step: safeStep(scale))
        ]
    }

    /// Cauchy(mu,γ) start: median and MAD for robust center/scale.
    /// - Returns: [mu (real), γ (positive)]
    func specs_cauchy(_ x: [T]) -> [ParamSpec<T>] {
        let clean = finiteSamples(x)
        let location = sampleMedian(clean, fallback: sampleMean(clean, fallback: T.zero))
        let mad = medianAbsoluteDeviation(clean, median: location)
        let scale = safePositive(mad, fallback: safePositive(sampleStd(clean, mean: location), fallback: T(1)))
        return [
            ParamSpec(.real, initial: location, step: safeStep(location == 0 ? T(1) : location)),
            ParamSpec(.positive, initial: scale, step: safeStep(scale))
        ]
    }

    /// Logistic location-scale start: mean and std translated to scale.
    /// - Returns: [mu (real), s (positive)]
    func specs_logistic(_ x: [T]) -> [ParamSpec<T>] {
        let clean = finiteSamples(x)
        let location = sampleMean(clean, fallback: T.zero)
        let std = sampleStd(clean, mean: location)
        let scale = safePositive(std * T.sqrt(T(3)) / T.pi, fallback: T(1))
        return [
            ParamSpec(.real, initial: location, step: safeStep(location == 0 ? T(1) : location)),
            ParamSpec(.positive, initial: scale, step: safeStep(scale))
        ]
    }

    // LogNormal(meanlog = mu, sdlog = sigma) ParamSpecs
    // mu is unbounded (any), sigma is positive.
    // Initials from log-data: mu0 = mean(log x), sigma0 = sqrt(var(log x)).
    // Optional overrides via ctx.ln_mu / ctx.ln_sigma if your context provides them.
    ///
    /// - Returns: [mu (real), sigma (positive)]
    func specs_LogNormal(_ x: [T], ctx: MLEStartContext<T>) -> [ParamSpec<T>] {
        // Support filter: x > 0
        let data = positiveSamples(x)

        // Compute log-data statistics with safe fallbacks
        let logs: [T] = data.map { T.log($0) }
        let mu0_raw = sampleMean(logs, fallback: T.zero)
        let v0_raw: T = {
            // If you have sampleVariance(_:fallback:), use it. Otherwise, simple two-pass:
            let m = mu0_raw
            let v = logs.reduce(T.zero) { $0 + ($1 - m) * ($1 - m) } / T(max(logs.count, 1))
            return v
        }()
        let sigma0_raw = T.sqrt(T.maximum(v0_raw, T(1e-8)))

        // Allow optional overrides from context if present
        let mu0 = ctx.ln_mu.map { $0 }.map { $0 } ?? mu0_raw
        let sigma0 = safePositive(ctx.ln_sigma.map { $0 }.map { $0 } ?? sigma0_raw, fallback: T(0.5))

        return [
            ParamSpec(.real,      initial: mu0,    step: safeStep(mu0)),
            ParamSpec(.positive, initial: sigma0, step: safeStep(sigma0))
        ]
    }

    /// Gumbel(mu,beta) start: translate mean/std to location/scale using Euler’s γ.
    /// - Returns: [mu (real), beta (positive)]
    func specs_gumbel(_ x: [T]) -> [ParamSpec<T>] {
        let clean = finiteSamples(x)
        let mean = sampleMean(clean, fallback: T.zero)
        let std = sampleStd(clean, mean: mean)
        let scale = safePositive(std * T.sqrt(T(6)) / T.pi, fallback: T(1))
        let eulerGamma = T(0.5772156649015329)
        let location = mean - eulerGamma * scale
        return [
            ParamSpec(.real, initial: location, step: safeStep(location == 0 ? T(1) : location)),
            ParamSpec(.positive, initial: scale, step: safeStep(scale))
        ]
    }

    /// Student’s t start: estimate nu from excess kurtosis if available; otherwise use a moderate nu.
    /// - Returns: [nu (positive)]
    func specs_studentsT(_ x: [T]) -> [ParamSpec<T>] {
        let clean = finiteSamples(x)
        let kurt = sampleExcessKurtosis(clean)
        var nuGuess = T(5)
        if kurt.isFinite && kurt > T.zero {
            // excess kurtosis of Student's t: 6/(nu−4) for nu>4
            let candidate = (T(6) / kurt) + T(4)
            if candidate.isFinite && candidate > T(0.5) {
                nuGuess = candidate
            }
        }
        nuGuess = safePositive(nuGuess, fallback: T(5))
        return [ParamSpec(.positive, initial: nuGuess, step: T(1))]
    }

    /// Fisher F(d1,d2) start via mean and variance relations with safeguards.
    /// - Returns: [d1 (positive), d2 (positive)]
    func specs_f(_ x: [T]) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let mean = sampleMean(data, fallback: T(2))
        let variance = sampleVariance(data, mean: mean)

        var d2Guess = safePositive(T(6), fallback: T(6))
        if mean > T.one + T(1e-3) {
            let candidate = (T(2) * mean) / (mean - T.one)
            if candidate.isFinite && candidate > T(4.5) {
                d2Guess = candidate
            }
        }

        var d1Guess = safePositive(T(5), fallback: T(5))
        if variance > T(1e-6) && d2Guess > T(4.5) {
            let numerator = T(2) * d2Guess * d2Guess * (d2Guess - T(2))
            let a = (d2Guess - T(2)) * (d2Guess - T(2)) * (d2Guess - T(4))
            let denominator = variance * a - T(2) * d2Guess * d2Guess
            if denominator > T(1e-6) {
                let candidate = numerator / denominator
                if candidate.isFinite && candidate > T(0.5) {
                    d1Guess = candidate
                }
            }
        }

        d1Guess = safePositive(d1Guess, fallback: T(5))
        d2Guess = safePositive(d2Guess, fallback: T(7))

        return [
            ParamSpec(.positive, initial: d1Guess, step: safeStep(d1Guess)),
            ParamSpec(.positive, initial: d2Guess, step: safeStep(d2Guess))
        ]
    }

    /// Inverse-Gamma(alpha,beta) start via MoM with clamps for alpha and beta.
    /// - Returns: [alpha (positive), beta (positive)]
    func specs_invGamma(_ x: [T]) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let mean = sampleMean(data, fallback: T(1))
        let variance = max(sampleVariance(data, mean: mean), T(1e-6))

        // Method-of-moments: alphâ = 2 + m² / s², betâ = m (alphâ - 1)
        var shape = T.two + (mean * mean) / variance
        if !shape.isFinite || shape <= T(1.1) {
            shape = safePositive(T(2.5), fallback: T(2.5))
        }
        var scale = mean * (shape - T.one)
        if !scale.isFinite || scale <= T(1e-6) {
            scale = safePositive(T.one, fallback: T.one)
        }

        return [
            ParamSpec(.positive, initial: shape, step: safeStep(shape)),
            ParamSpec(.positive, initial: scale, step: safeStep(scale))
        ]
    }

    /// Inverse-χ²(nu, τ²) start: nu from sample size and τ² from mean.
    /// - Returns: [nu (positive), τ² (positive)]
    func specs_inverseChiSquared(_ x: [T]) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let scale = safePositive(sampleMean(data, fallback: T(1)), fallback: T(1))
        let nu = safePositive(T(max(5, data.count)), fallback: T(5))
        return [
            ParamSpec(.positive, initial: nu, step: T(1)),
            ParamSpec(.positive, initial: scale, step: safeStep(scale))
        ]
    }

    /// Skew-Normal(xi,ω,alpha) start: mean/std for xi,ω and alpha=0.
    /// - Returns: [xi (real), ω (positive), alpha (real)]
    func specs_skewNormal(_ x: [T]) -> [ParamSpec<T>] {
        let clean = finiteSamples(x)
        let location = sampleMean(clean, fallback: T.zero)
        let scale = safePositive(sampleStd(clean, mean: location), fallback: T(1))
        return [
            ParamSpec(.real, initial: location, step: safeStep(location == 0 ? T(1) : location)),
            ParamSpec(.positive, initial: scale, step: safeStep(scale)),
            ParamSpec(.real, initial: T.zero, step: T(1))
        ]
    }

    /// Skew-t(mu,sigma,nu,alpha) start: mean/std for mu,sigma; nu moderate; alpha=0.
    /// - Returns: [mu (real), sigma (positive), nu (positive), alpha (real)]
    func specs_skewT(_ x: [T]) -> [ParamSpec<T>] {
        let clean = finiteSamples(x)
        let location = sampleMean(clean, fallback: T.zero)
        let scale = safePositive(sampleStd(clean, mean: location), fallback: T(1))
        return [
            ParamSpec(.real, initial: location, step: safeStep(location == 0 ? T(1) : location)),
            ParamSpec(.positive, initial: scale, step: safeStep(scale)),
            ParamSpec(.positive, initial: T(5), step: T(1)),
            ParamSpec(.real, initial: T.zero, step: T(1))
        ]
    }

    /// GEV(mu,sigma,xi) start: mean/std for mu,sigma; xi=0 (Gumbel) as a neutral start.
    /// - Returns: [mu (real), sigma (positive), xi (real)]
    func specs_gev(_ x: [T]) -> [ParamSpec<T>] {
        let clean = finiteSamples(x)
        let location = sampleMean(clean, fallback: T.zero)
        let scale = safePositive(sampleStd(clean, mean: location), fallback: T(1))
        let shape = T.zero
        return [
            ParamSpec(.real, initial: location, step: safeStep(location == 0 ? T(1) : location)),
            ParamSpec(.positive, initial: scale, step: safeStep(scale)),
            ParamSpec(.real, initial: shape, step: T(0.5))
        ]
    }

    /// GPD(sigma,xi) start on exceedances y>=0: sigma from mean(y), xi=0.
    /// - Returns: [sigma (positive), xi (real)]
    func specs_gpd(_ x: [T]) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let sigma = safePositive(sampleMean(data, fallback: T(1)), fallback: T(1))
        let xi = T.zero
        return [
            ParamSpec(.positive, initial: sigma, step: safeStep(sigma)),
            ParamSpec(.real, initial: xi, step: T(0.5))
        ]
    }

    /// Log-Logistic(alpha,beta) start: alpha from median/mean, beta moderate.
    /// - Returns: [alpha (positive), beta (positive)]
    func specs_logLogistic(_ x: [T]) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let scale = safePositive(sampleMedian(data, fallback: sampleMean(data, fallback: T(1))), fallback: T(1))
        let shape = safePositive(T(2), fallback: T(1))
        return [
            ParamSpec(.positive, initial: scale, step: safeStep(scale)),
            ParamSpec(.positive, initial: shape, step: T(0.5))
        ]
    }

    /// Generalized Gamma(k,theta,p) start: theta from mean, k≈1.5, p≈1.
    /// - Returns: [k (positive), theta (positive), p (positive)]
    func specs_genGamma(_ x: [T]) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let theta = safePositive(sampleMean(data, fallback: T(1)), fallback: T(1))
        return [
            ParamSpec(.positive, initial: T(1.5), step: T(0.5)),
            ParamSpec(.positive, initial: theta, step: safeStep(theta)),
            ParamSpec(.positive, initial: T(1), step: T(0.5))
        ]
    }

    /// Generalized Normal / Exponential Power(mu,alpha,beta) start: mean/std for mu,alpha; beta≈2.
    /// - Returns: [mu (real), alpha (positive), beta (positive)]
    func specs_genNormal(_ x: [T]) -> [ParamSpec<T>] {
        let clean = finiteSamples(x)
        let mu = sampleMean(clean, fallback: T.zero)
        let alpha = safePositive(sampleStd(clean, mean: mu), fallback: T(1))
        let beta = safePositive(T(2), fallback: T(1))
        return [
            ParamSpec(.real, initial: mu, step: safeStep(mu == 0 ? T(1) : mu)),
            ParamSpec(.positive, initial: alpha, step: safeStep(alpha)),
            ParamSpec(.positive, initial: beta, step: T(0.5))
        ]
    }

    /// Burr XII(c,k,s) start: moderate shapes and scale from mean.
    /// - Returns: [c (positive), k (positive), s (positive)]
    func specs_burrXII(_ x: [T]) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let s = safePositive(sampleMean(data, fallback: T(1)), fallback: T(1))
        return [
            ParamSpec(.positive, initial: T(2), step: T(0.5)),
            ParamSpec(.positive, initial: T(2), step: T(0.5)),
            ParamSpec(.positive, initial: s, step: safeStep(s))
        ]
    }

    /// Johnson SU(γ,δ,xi,lambda) start: xi≈mean, lambda≈std, γ≈0, δ≈1.
    /// - Returns: [γ (real), δ (positive), xi (real), lambda (positive)]
    func specs_johnsonSU(_ x: [T]) -> [ParamSpec<T>] {
        let clean = finiteSamples(x)
        let xi = sampleMean(clean, fallback: T.zero)
        let lambda = safePositive(sampleStd(clean, mean: xi), fallback: T(1))
        return [
            ParamSpec(.real, initial: T.zero, step: T(0.5)),
            ParamSpec(.positive, initial: T(1), step: T(0.5)),
            ParamSpec(.real, initial: xi, step: safeStep(xi == 0 ? T(1) : xi)),
            ParamSpec(.positive, initial: lambda, step: safeStep(lambda))
        ]
    }

    /// Noncentral χ²: estimate lambda with k possibly provided in `ctx`; data must be x>0.
    /// - Returns: [lambda (positive)]
    func specs_ncChiSq(_ x: [T], ctx: MLEStartContext<T>) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let mean = sampleMean(data, fallback: T(2))
        let k = safePositive(ctx.ncChi_k.map { T($0) } ?? mean, fallback: T(2))
        let lambda = safePositive(max(mean - k, T(1e-3)), fallback: T(1))
        return [ParamSpec(.positive, initial: lambda, step: safeStep(lambda))]
    }

    /// Noncentral t: estimate δ from sample mean (mu-like), with step scaled by magnitude.
    /// - Returns: [δ (real)]
    func specs_ncT(_ x: [T], ctx: MLEStartContext<T>) -> [ParamSpec<T>] {
        let clean = finiteSamples(x)
        let deltaGuess = sampleMean(clean, fallback: T.zero)
        let stepBase = T.maximum(abs(deltaGuess), T.one)
        return [ParamSpec(.real, initial: deltaGuess, step: safeStep(stepBase))]
    }
    
    /// Central χ²(k) start: k≈mean(x) or overridden via context.
    /// - Returns: [k (positive)]
    func specs_ChSq(_ x: [T], ctx: MLEStartContext<T>) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let mean = sampleMean(data, fallback: T(2))
        let k0 = ctx.ncChi_k
            .map { T($0) }                // optional integer df from context
            .map { safePositive($0, fallback: T(2)) }
            ?? safePositive(mean, fallback: T(2))
        return [ParamSpec(.positive, initial: k0, step: safeStep(k0))]
    }


    /// Noncentral F: estimate lambda from mean with df possibly provided in `ctx`.
    /// - Returns: [lambda (positive)]
    func specs_ncF(_ x: [T], ctx: MLEStartContext<T>) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let mean = sampleMean(data, fallback: T(2))
        let d1 = safePositive(ctx.ncF_df1.map { T($0) } ?? (mean + T(2)), fallback: T(5))
        let d2 = safePositive(ctx.ncF_df2.map { T($0) } ?? (mean + T(4)), fallback: T(7))
        var lambda = T(1)
        let denom = d1 * (d2 - T(2))
        if denom > T.zero {
            let candidate = (mean * denom / d2) - d1
            if candidate.isFinite {
                lambda = candidate
            }
        }
        lambda = safePositive(lambda, fallback: T(1))
        return [ParamSpec(.positive, initial: lambda, step: safeStep(lambda))]
    }

    /// Noncentral Beta: extend Beta(alpha,beta) with lambda≈0.5 as a neutral start.
    /// - Returns: [alpha (positive), beta (positive), lambda (positive)]
    func specs_ncBeta(_ x: [T]) -> [ParamSpec<T>] {
        let base = specs_beta(x)
        let lambda = safePositive(T(0.5), fallback: T(0.5))
        return base + [ParamSpec(.positive, initial: lambda, step: safeStep(lambda))]
    }

    /// Bernoulli(p) start: p≈mean(x) clamped to (ε, 1−ε).
    /// - Returns: [p (unit interval)]
    func specs_bernoulli(_ x: [T]) -> [ParamSpec<T>] {
        let clean = clamp01Samples(x)
        let mean = sampleMean(clean, fallback: T(0.5))
        let p = clampProbability(mean)
        return [ParamSpec(.unitInterval, initial: p, step: T(0.1))]
    }

    /// Binomial(n,p) start: n from context or max(x), p≈mean(x/n) clamped.
    /// - Returns: [p (unit interval)]
    func specs_binomial(_ x: [T], ctx: MLEStartContext<T>) -> [ParamSpec<T>] {
        let n = max(ctx.binomialN ?? Int(x.max() ?? 0), 1)
        let total = sampleMean(clamp01Samples(x.map { max($0, T.zero) / T(n) }), fallback: T(0.5))
        let p = clampProbability(total)
        return [ParamSpec(.unitInterval, initial: p, step: T(0.1))]
    }

    /// Negative Binomial(r,p) start: r≈mean, p≈r/(r+mean) clamped.
    /// - Returns: [r (positive), p (unit interval)]
    func specs_negativeBinomial(_ x: [T]) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let mean = sampleMean(data, fallback: T(1))
        let r = safePositive(mean, fallback: T(1))
        let p = clampProbability(r / (r + mean))
        return [
            ParamSpec(.positive, initial: r, step: T(1)),
            ParamSpec(.unitInterval, initial: p, step: T(0.1))
        ]
    }

    /// Geometric(p) start: p≈1/(mean+1) clamped.
    /// - Returns: [p (unit interval)]
    func specs_geometric(_ x: [T]) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let mean = sampleMean(data, fallback: T(1))
        let p = clampProbability(T(1) / (mean + T(1)))
        return [ParamSpec(.unitInterval, initial: p, step: T(0.1))]
    }

    /// Hypergeometric(N,K,n) start: derive N,K,n from sample mean or context hints.
    /// - Returns: [N (positive), K (positive), n (positive)]
    func specs_hypergeometric(_ x: [T], ctx: MLEStartContext<T>) -> [ParamSpec<T>] {
        let mean = sampleMean(finiteSamples(x), fallback: T(1))
        let population = max(ctx.hyperN ?? Int(mean * T(4) + T(10)), 10)
        let draws = max(ctx.hyperSampleSize ?? Int(mean + T(1)), 1)
        let successes = min(population - 1, max(Int(mean), 1))
        return [
            ParamSpec(.positive, initial: T(population), step: T(1)),
            ParamSpec(.positive, initial: T(successes), step: T(1)),
            ParamSpec(.positive, initial: T(draws), step: T(1))
        ]
    }

    /// Poisson(lambda) start: lambda≈mean(x) with positivity guard.
    /// - Returns: [lambda (positive)]
    func specs_poisson(_ x: [T]) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        let lambda = safePositive(sampleMean(data, fallback: T(1)), fallback: T(1))
        return [ParamSpec(.positive, initial: lambda, step: safeStep(lambda))]
    }

    /// Holtsmark location-scale start: mean/std with sigma>0 clamp.
    /// - Returns: [mu (real), sigma (positive)]
    func specs_holtsmark(_ x: [T]) -> [ParamSpec<T>] {
        let clean = finiteSamples(x)
        let location = sampleMean(clean, fallback: T.zero)
        let scale = safePositive(sampleStd(clean, mean: location), fallback: T(1))
        return [
            ParamSpec(.real, initial: location, step: safeStep(location == 0 ? T(1) : location)),
            ParamSpec(.positive, initial: scale, step: safeStep(scale))
        ]
    }

    /// Rayleigh(sigma) start from second moment: sigma = sqrt(mean(x²)/2).
    /// - Returns: [sigma (positive)]
    func specs_rayleigh(_ x: [T], ctx: MLEStartContext<T>) -> [ParamSpec<T>] {
        let data = positiveSamples(x)
        precondition(!data.isEmpty)
        // MLE: sigma^2 = mean(x^2)/2
        var s2 = T.zero
        for v in data { s2 += v*v }
        let sigma0 = T.sqrt(T.maximum(s2 / (T(2) * T(data.count)), T(1e-8)))
        return [ParamSpec(.positive, initial: sigma0, step: safeStep(sigma0))]
    }
    
    // MARK: - Helpers

    /// Filter finite values; if empty, return [0] as a neutral fallback.
    func finiteSamples(_ x: [T]) -> [T] {
        let filtered = x.filter { $0.isFinite }
        return filtered.isEmpty ? [T.zero] : filtered
    }

    /// Filter strictly positive finite values; if empty, return [1] as fallback.
    func positiveSamples(_ x: [T]) -> [T] {
        let filtered = x.filter { $0.isFinite && $0 > 0 }
        return filtered.isEmpty ? [T(1)] : filtered
    }

    /// Keep values strictly inside (0,1); if empty, return [0.5] as fallback.
    func clamp01Samples(_ x: [T]) -> [T] {
        let filtered = x.filter { $0.isFinite && $0 > 0 && $0 < 1 }
        return filtered.isEmpty ? [T(0.5)] : filtered
    }

    /// Arithmetic mean with fallback for empty input.
    func sampleMean(_ x: [T], fallback: T) -> T {
        guard !x.isEmpty else { return fallback }
        let sum = x.reduce(T.zero, +)
        return sum / T(x.count)
    }

    /// Population variance with fallback for n<=1.
    func sampleVariance(_ x: [T], mean: T) -> T {
        guard x.count > 1 else { return T.zero }
        let sum = x.reduce(T.zero) { $0 + ( $1 - mean) * ( $1 - mean) }
        return sum / T(x.count)
    }

    /// Standard deviation from `sampleVariance(_:mean:)`.
    func sampleStd(_ x: [T], mean: T) -> T {
        T.sqrt(max(sampleVariance(x, mean: mean), T.zero))
    }

    /// Convenience overload: std computed around the sample mean.
    func sampleStd(_ x: [T]) -> T {
        let mean = sampleMean(x, fallback: T.zero)
        return sampleStd(x, mean: mean)
    }

    /// Median of finite values; falls back to provided value if no finite values exist.
    func sampleMedian(_ x: [T], fallback: T) -> T {
        var sorted = x.filter { $0.isFinite }
        guard !sorted.isEmpty else { return fallback }
        sorted.sort()
        let mid = sorted.count / 2
        if sorted.count % 2 == 0 {
            return (sorted[mid - 1] + sorted[mid]) / T(2)
        } else {
            return sorted[mid]
        }
    }

    /// Median absolute deviation around a supplied median, with fallback to std.
    func medianAbsoluteDeviation(_ x: [T], median: T) -> T {
        let deviations = x.map { abs($0 - median) }
        return sampleMedian(deviations, fallback: sampleStd(x, mean: median))
    }

    /// Compute mean/variance after trimming a fraction from both tails.
    func trimmedStatistics(_ x: [T], trimFraction: T) -> (mean: T, variance: T) {
        let trimmed = trimmedData(x, trimFraction: trimFraction)
        let mean = sampleMean(trimmed, fallback: T.zero)
        let variance = sampleVariance(trimmed, mean: mean)
        return (mean, variance)
    }

    /// Trim tails and return the median of the trimmed sample.
    func trimmedMedian(_ x: [T], trimFraction: T = T(0.05)) -> T {
        let trimmed = trimmedData(x, trimFraction: trimFraction)
        return sampleMedian(trimmed, fallback: sampleMean(trimmed, fallback: T.zero))
    }

    /// Return the central slice after removing a fraction from both tails.
    /// The trim fraction is clamped to [0, 0.45] to ensure at least one element remains.
    func trimmedData(_ x: [T], trimFraction: T) -> [T] {
        guard !x.isEmpty else { return [T.zero] }
        var sorted = x.filter { $0.isFinite }
        if sorted.isEmpty { return [T.zero] }
        sorted.sort()
        let count = sorted.count
        let frac = max(0.0, min(0.45, Double(trimFraction)))
        let trim = Int((frac * Double(count)).rounded(.towardZero))
        let start = min(max(trim, 0), count - 1)
        let end = max(count - trim, start + 1)
        return Array(sorted[start..<end])
    }

    /// Excess kurtosis (population definition) with safe fallbacks.
    func sampleExcessKurtosis(_ x: [T]) -> T {
        guard x.count > 3 else { return T.zero }
        let mean = sampleMean(x, fallback: T.zero)
        var m2: T = 0
        var m4: T = 0
        for value in x {
            let d = value - mean
            let d2 = d * d
            m2 += d2
            m4 += d2 * d2
        }
        let n = T(x.count)
        if m2 == 0 { return T.zero }
        let var2 = (m2 / n) * (m2 / n)
        let kurt = (m4 / n) / var2
        return kurt - T(3)
    }

    /// Clamp to a strictly positive value, preserving finiteness and providing a minimum floor.
    func safePositive(_ value: T, fallback: T) -> T {
        let v = value.isFinite ? value : fallback
        return max(v, fallback * T(1e-3) + T(1e-3))
    }

    /// Derive a reasonable initial step size in u-space from a theta-space magnitude.
    /// Ensures a minimum step to avoid degenerate simplices or zero-length directions.
    func safeStep(_ value: T) -> T {
        let magnitude = max(abs(value), T(1))
        return max(magnitude * T(0.1), T(0.2))
    }

    /// Clamp a probability to (ε, 1−ε) to avoid boundary issues in logit/log transforms.
    func clampProbability(_ value: T) -> T {
        let eps = T(1e-3)
        return min(max(value, eps), T.one - eps)
    }
}
