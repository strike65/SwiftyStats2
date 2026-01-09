//
//  Created by VT on 08.11.25.
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

/// Random sampling utilities and generators used across SwiftyStats distributions.

import SwiftyStatsPrelude

/// Random number generation helpers for various distributions used in tests or sampling.
///
/// Overview:
/// - All functions are generic over `T: RealLike`.
/// - Reproducibility is supported via an internal deterministic LCG (`SeededGenerator`) when a seed is provided.
/// - Unless otherwise noted, parameters use the common statistics conventions for each distribution:
///   - Gamma uses shape k and scale θ.
///   - Normal uses mean μ and standard deviation σ.
///   - Student’s t uses location μ, scale σ, and degrees of freedom ν.
///   - Logistic uses location μ and scale s.
///   - Cauchy uses location μ and scale γ.
///   - GEV uses location μ, scale σ, and shape ξ.
///   - GPD uses scale σ and shape ξ.
/// - Each function returns a single variate per call; `randoms(n:dist:seed:)` returns an array.
///
/// Algorithms used (concise):
/// - Normal: Marsaglia polar (Box–Muller polar form) for standard normal, then affine transform.
/// - Gamma: Marsaglia–Tsang; for k < 1 uses Johnk’s trick to reduce to k ≥ 1 case.
/// - Beta: Ratio of two Gamma variables.
/// - Chi-square: Gamma with (df/2, 2).
/// - Fisher F: Ratio of scaled chi-squares.
/// - Exponential: Inverse transform.
/// - Weibull: Inverse transform.
/// - Arcsine: Scaled Beta(1/2, 1/2).
/// - Student’s t: Normal divided by sqrt(χ²/ν) with affine transform (μ, σ).
/// - Cauchy: Inverse CDF (quantile).
/// - Logistic: Inverse CDF (quantile).
/// - Gumbel: Inverse CDF (quantile).
/// - GEV / GPD: Inverse CDF with ξ→0 limit handled by Gumbel/exponential forms.
/// - Skew-Normal: Azzalini’s construction.
/// - Noncentral χ²: Poisson-mixture; large-λ branch uses normal approximation for Poisson count.
/// - Noncentral t: Shifted normal over scaled χ².
/// - Inverse Gaussian (Wald): Michael–Schucany–Haas acceptance rule.
/// - Log-normal: Box–Muller (classical form) then exponentiate.
/// - Discrete (Bernoulli, Binomial, Geometric, Negative-Binomial, Poisson): standard inverse/compound constructions.
/// - Rayleigh: Inverse transform.
/// - Landau, Map Airy, SAS 0.5, Holtsmark: inversion via SwiftyBoost quantiles or Chambers–Mallows–Stuck-style formulas where applicable.
public enum RNGSampler<T: RealLike> {
    
    /// Enumerates supported distributions for `RNGSampler.randoms(n:dist:seed:)`.
    ///
    /// Cases group parameters in conventional order. See each case for the parameterization.
    public enum RandomNumbers {
        case uniform01
        case bernoulli(p: T)
        case binomial(nTrials:Int, p: T)
        case negativeBinomial(r: T, p: T)
        case negativeBinomialIntR(r: Int, p: T)
        case gaussian(mean: T, standardDeviation: T)
        case gamma(shape: T, scale: T)
        case beta(a: T, b: T)
        case chiSquare(df: T)
        case f(d1: T, d2: T)
        case arcsine(a: T, b: T)
        case exponential(rate: T)
        case weibull(k: T, lambda: T)
        case studentT(mu: T, sigma: T, nu: T)
        case cauchy(mu: T, gamma: T)
        case logistic(mu: T, s: T)
        case gumbel(mu: T, beta: T)
        case gev(mu: T, sigma: T, xi: T)
        case gpd(sigma: T, xi: T)
        case skewNormal(xi: T, omega: T, alpha: T)
        case noncentralChiSquare(k: T, lambda: T)
        case noncentralT(mu: T, sigma: T, nu: T, delta: T)
        case sasPoint5(location: T, sigma: T)
        case holtsmark(mu: T, c: T)
        case lognormal(lmean: T, lsigma: T)
        case pareto(scale: T, shape: T)
        case poisson(lambda: T)
        case rayleigh(scale: T)
        case geometric(p: T)
        case inverseNormal(mu: T, sigma: T)
        case triangular(lower: T, upper: T, mode: T)
        case landau(location: T, scale: T)
        case mapAiry(location: T, scale: T)
    }
    
    /// Deterministic 64-bit LCG for reproducible streams when a seed is provided.
    ///
    /// Notes:
    /// - Fast, simple, not cryptographically secure.
    /// - Used internally by `randoms(n:dist:seed:)` when `seed != nil`.
    private struct SeededGenerator: RandomNumberGenerator {
        private var state: UInt64
        init(seed: UInt64) {
            self.state = seed != 0 ? seed : 0x9E3779B185EBCA87
        }
        mutating func next() -> UInt64 {
            state &*= 2862933555777941757
            state &+= 3037000493
            return state
        }
    }
    
    /// Generate `n` random variates from the specified distribution.
    ///
    /// - Parameters:
    ///   - n: Number of variates to generate (`n > 0`).
    ///   - dist: Distribution and parameters (see `RandomNumbers`).
    ///   - seed: Optional seed for reproducibility. When provided, a deterministic LCG is used;
    ///           otherwise `SystemRandomNumberGenerator` is used.
    /// - Returns: An array of `n` variates of type `T`.
    ///
    /// Distribution methods (selection):
    /// - `.gaussian`: Marsaglia polar (Box–Muller polar) + affine transform.
    /// - `.gamma`: Marsaglia–Tsang; Johnk’s reduction when shape < 1.
    /// - `.beta`: Ratio of Gammas.
    /// - `.chiSquare`: Gamma(df/2, 2).
    /// - `.f`: Ratio of scaled chi-squares.
    /// - `.exponential`, `.weibull`, `.rayleigh`, `.cauchy`, `.logistic`, `.gumbel`, `.gev`, `.gpd`: inverse CDF.
    /// - `.studentT`: Normal over sqrt(χ²/ν) with location/scale.
    /// - `.skewNormal`: Azzalini.
    /// - `.noncentralChiSquare`: Poisson mixture with large-λ acceleration.
    /// - `.noncentralT`: Shifted normal over scaled χ².
    /// - `.inverseNormal`: Michael–Schucany–Haas.
    /// - `.lognormal`: Box–Muller + exp.
    /// - Discrete families use standard constructions (Knuth for small-λ Poisson, inverse for geometric, IID Bernoulli sum for binomial, sum of geometrics for negative binomial).
    public static func randoms(n: Int, dist: RandomNumbers, seed: UInt64? = nil) -> [T] {
        precondition(n > 0)
        
        if let s = seed {
            var rng = SeededGenerator(seed: s)
            switch dist {
            case .mapAiry(location: let mu, scale: let sigma):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, mu: T, scale: T) -> T in
                                        rmapAiry(mu: mu, sigma: scale, &rng)
                                    },
                                rng: &rng,
                                mu,
                                sigma)
            case .landau(location: let mu, scale: let scale):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, mu: T, scale: T) -> T in
                                        rlandau(mu: mu, c: scale, &rng)
                                    },
                                rng: &rng,
                                mu,
                                scale)
            case .triangular(lower: let a, upper: let b, mode: let c):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, a: T, b: T, c: T) -> T in
                                    rtriangular(lower: a, upper: b, mode: c, &rng)
                                    },
                                rng: &rng,
                                a,
                                b,
                                c)
            case .negativeBinomial(r: let r, p: let p):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, r: T, p: T) -> T in
                                    rnegativeBinomial(r: r, p: p, &rng)
                                    },
                                rng: &rng,
                                r,
                                p)
            case .negativeBinomialIntR(r: let r, p: let p):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, nTrials: Int, p: T) -> T in
                                    rnegativeBinomial(r: r, p: p, &rng)
                                    },
                                rng: &rng,
                                r,
                                p)

            case .inverseNormal(mu: let mu, sigma: let sigma):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, mu: T, sigma: T) -> T in
                                    rinorm(mu, lambda: sigma, &rng)
                                    },
                                rng: &rng,
                                mu,
                                sigma)
            case .geometric(p: let p):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, p: T) -> T in
                    rgeometric(p: p, &rng)
                },
                                rng: &rng,
                                p)
            case .rayleigh(scale: let sigma):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, sigma: T) -> T in
                    rrayleigh(sigma: sigma, &rng)
                },
                                rng: &rng,
                                sigma)
            case .poisson(lambda: let l):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, l: T) -> T in
                    rpoisson(lambda: l, &rng)
                },
                                rng: &rng,
                                l)
            case .pareto(scale: let a, shape: let b):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, a: T, b: T) -> T in
                    rpareto(scale: a, shape: b, &rng)
                },
                                rng: &rng,
                                a,
                                b)
            case .lognormal(lmean: let lmean, lsigma: let lsigma):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, lmean: T, lsigma: T) -> T in
                    rlognormal(meanlog: lmean, sdlog: lsigma, &rng)
                },
                                rng: &rng,
                                lmean,
                                lsigma)
            case .uniform01:
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator) -> T in
                    u01(&rng)
                },
                                rng: &rng)
            case .bernoulli(p: let p):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, p: T) -> T in
                    rbernoulli(p, &rng)
                },
                                rng: &rng,
                                p)
            case .binomial(nTrials: let nTrials, p: let p):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, nTrials: Int, p: T) -> T in
                    rbinomial(numberOfTrials: nTrials, probabilityOfSuccess: p, &rng)
                },
                                rng: &rng,
                                nTrials,
                                p)
            case .gaussian(mean: let mu, standardDeviation: let sigma):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, mu: T, s: T) -> T in
                    rnorm(mu, s, &rng)
                },
                                rng: &rng,
                                mu,
                                sigma)
            case .gamma(shape: let k, scale: let theta):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, k: T, theta: T) -> T in
                    rgamma(k, theta, &rng)
                },
                                rng: &rng,
                                k,
                                theta)
            case .beta(a: let a, b: let b):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, a: T, b: T) -> T in
                    rbeta(a, b, &rng)
                },
                                rng: &rng,
                                a,
                                b)
            case .chiSquare(df: let df):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, df: T) -> T in
                    rchisq(df, &rng)
                },
                                rng: &rng,
                                df)
            case .f(d1: let d1, d2: let d2):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, d1: T, d2: T) -> T in
                    rf(d1, d2, &rng)
                },
                                rng: &rng,
                                d1,
                                d2)
            case .arcsine(a: let a, b: let b):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, a: T, b: T) -> T in
                    rarcsine(a, b, &rng)
                },
                                rng: &rng,
                                a,
                                b)
            case .exponential(rate: let rate):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, rate: T) -> T in
                    rexp(rate, &rng)
                },
                                rng: &rng,
                                rate)
            case .weibull(k: let k, lambda: let lambda):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, k: T, lambda: T) -> T in
                    rweibull(k, lambda, &rng)
                },
                                rng: &rng,
                                k,
                                lambda)
            case .studentT(mu: let mu, sigma: let sigma, nu: let nu):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, mu: T, sigma: T, nu: T) -> T in
                    rt(mu, sigma, nu, &rng)
                },
                                rng: &rng,
                                mu,
                                sigma,
                                nu)
            case .cauchy(mu: let mu, gamma: let gamma):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, mu: T, gamma: T) -> T in
                    rcauchy(mu, gamma, &rng)
                },
                                rng: &rng,
                                mu,
                                gamma)
            case .logistic(mu: let mu, s: let sParam):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, mu: T, s: T) -> T in
                    rlogistic(mu, s, &rng)
                },
                                rng: &rng,
                                mu,
                                sParam)
            case .gumbel(mu: let mu, beta: let beta):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, mu: T, beta: T) -> T in
                    rgumbel(mu, beta, &rng)
                },
                                rng: &rng,
                                mu,
                                beta)
            case .gev(mu: let mu, sigma: let sigma, xi: let xi):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, mu: T, sigma: T, xi: T) -> T in
                    rgev(mu, sigma, xi, &rng)
                },
                                rng: &rng,
                                mu,
                                sigma,
                                xi)
            case .gpd(sigma: let sigma, xi: let xi):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, sigma: T, xi: T) -> T in
                    rgpd(sigma, xi, &rng)
                },
                                rng: &rng,
                                sigma,
                                xi)
            case .skewNormal(xi: let xi, omega: let omega, alpha: let alpha):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, xi: T, omega: T, alpha: T) -> T in
                    rskewNormal(xi, omega, alpha, &rng)
                },
                                rng: &rng,
                                xi,
                                omega,
                                alpha)
            case .noncentralChiSquare(k: let k, lambda: let lambda):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, k: T, lambda: T) -> T in
                    rncchisq(k, lambda, &rng)
                },
                                rng: &rng,
                                k,
                                lambda)
            case .noncentralT(mu: let mu, sigma: let sigma, nu: let nu, delta: let delta):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, mu: T, sigma: T, nu: T, delta: T) -> T in
                    rnct(mu, sigma, nu, delta, &rng)
                },
                                rng: &rng,
                                mu,
                                sigma,
                                nu,
                                delta)
            case .sasPoint5(location: let location, sigma: let sigma):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, location: T, sigma: T) -> T in
                                    rsasPoint5(mu: location, sigma: sigma, &rng)
                                    },
                                rng: &rng,
                                location,
                                sigma)
            case .holtsmark(mu: let mu, c: let c):
                return rsampler(n,
                                using:
                                    { (rng: inout SeededGenerator, mu: T, c: T) -> T in
                    rholtsmark(mu, c, &rng)
                },
                                rng: &rng,
                                mu,
                                c)
            }
        } else {
            var rng = SystemRandomNumberGenerator()
            switch dist {
            case .mapAiry(location: let mu, scale: let sigma):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, mu: T, scale: T) -> T in
                                        rmapAiry(mu: mu, sigma: scale, &rng)
                                    },
                                rng: &rng,
                                mu,
                                sigma)

            case .landau(location: let mu, scale: let scale):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, mu: T, scale: T) -> T in
                                        rlandau(mu: mu, c: scale, &rng)
                                    },
                                rng: &rng,
                                mu,
                                scale)
            case .triangular(lower: let a, upper: let b, mode: let c):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, a: T, b: T, c: T) -> T in
                                    rtriangular(lower: a, upper: b, mode: c, &rng)
                                    },
                                rng: &rng,
                                a,
                                b,
                                c)

            case .negativeBinomial(r: let r, p: let p):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, r: T, p: T) -> T in
                                    rnegativeBinomial(r: r, p: p, &rng)
                                    },
                                rng: &rng,
                                r,
                                p)
            case .negativeBinomialIntR(r: let r, p: let p):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, r: Int, p: T) -> T in
                                    rnegativeBinomial(r: r, p: p, &rng)
                                    },
                                rng: &rng,
                                r,
                                p)

            case .inverseNormal(mu: let mu, sigma: let sigma):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, mu: T, sigma: T) -> T in
                                    rinorm(mu, lambda: sigma, &rng)
                                    },
                                rng: &rng,
                                mu,
                                sigma)
            case .geometric(p: let p):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, p: T) -> T in
                    rgeometric(p: p, &rng)
                },
                                rng: &rng,
                                p)
                
            case .rayleigh(scale: let sigma):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, sigma: T) -> T in
                    rrayleigh(sigma: sigma, &rng)
                },
                                rng: &rng,
                                sigma)
            case .poisson(lambda: let l):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, l: T) -> T in
                    rpoisson(lambda: l, &rng)
                },
                                rng: &rng,
                                l)
            case .pareto(scale: let a, shape: let b):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, a: T, b: T) -> T in
                    rpareto(scale: a, shape: b, &rng)
                },
                                rng: &rng,
                                a,
                                b)
            case .lognormal(lmean: let mu, lsigma: let s):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, mu: T, s: T) -> T in
                    rlognormal(meanlog: mu, sdlog: s, &rng)
                },
                                rng: &rng,
                                mu,
                                s)
            case .uniform01:
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator) -> T in
                    u01(&rng)
                },
                                rng: &rng)
            case .bernoulli(p: let p):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, p: T) -> T in
                    rbernoulli(p, &rng)
                },
                                rng: &rng,
                                p)
            case .binomial(nTrials: let nTrials, p: let p):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator,nTrials: Int, p: T) -> T in
                    rbinomial(numberOfTrials: nTrials, probabilityOfSuccess: p, &rng)
                },
                                rng: &rng,
                                nTrials,
                                p)
            case .gaussian(mean: let mu, standardDeviation: let sigma):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, mu: T, s: T) -> T in
                    rnorm(mu, s, &rng)
                },
                                rng: &rng,
                                mu,
                                sigma)
            case .gamma(shape: let k, scale: let theta):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, k: T, theta: T) -> T in
                    rgamma(k, theta, &rng)
                },
                                rng: &rng,
                                k,
                                theta)
            case .beta(a: let a, b: let b):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, a: T, b: T) -> T in
                    rbeta(a, b, &rng)
                },
                                rng: &rng,
                                a,
                                b)
            case .chiSquare(df: let df):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, df: T) -> T in
                    rchisq(df, &rng)
                },
                                rng: &rng,
                                df)
            case .f(d1: let d1, d2: let d2):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, d1: T, d2: T) -> T in
                    rf(d1, d2, &rng)
                },
                                rng: &rng,
                                d1,
                                d2)
            case .arcsine(a: let a, b: let b):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, a: T, b: T) -> T in
                    rarcsine(a, b, &rng)
                },
                                rng: &rng,
                                a,
                                b)
            case .exponential(rate: let rate):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, rate: T) -> T in
                    rexp(rate, &rng)
                },
                                rng: &rng,
                                rate)
            case .weibull(k: let k, lambda: let lambda):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, k: T, lambda: T) -> T in
                    rweibull(k, lambda, &rng)
                },
                                rng: &rng,
                                k,
                                lambda)
            case .studentT(mu: let mu, sigma: let sigma, nu: let nu):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, mu: T, sigma: T, nu: T) -> T in
                    rt(mu, sigma, nu, &rng)
                },
                                rng: &rng,
                                mu,
                                sigma,
                                nu)
            case .cauchy(mu: let mu, gamma: let gamma):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, mu: T, gamma: T) -> T in
                    rcauchy(mu, gamma, &rng)
                },
                                rng: &rng,
                                mu,
                                gamma)
            case .logistic(mu: let mu, s: let sParam):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, mu: T, s: T) -> T in
                    rlogistic(mu, s, &rng)
                },
                                rng: &rng,
                                mu,
                                sParam)
            case .gumbel(mu: let mu, beta: let beta):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, mu: T, beta: T) -> T in
                    rgumbel(mu, beta, &rng)
                },
                                rng: &rng,
                                mu,
                                beta)
            case .gev(mu: let mu, sigma: let sigma, xi: let xi):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, mu: T, sigma: T, xi: T) -> T in
                    rgev(mu, sigma, xi, &rng)
                },
                                rng: &rng,
                                mu,
                                sigma,
                                xi)
            case .gpd(sigma: let sigma, xi: let xi):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, sigma: T, xi: T) -> T in
                    rgpd(sigma, xi, &rng)
                },
                                rng: &rng,
                                sigma,
                                xi)
            case .skewNormal(xi: let xi, omega: let omega, alpha: let alpha):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, xi: T, omega: T, alpha: T) -> T in
                    rskewNormal(xi, omega, alpha, &rng)
                },
                                rng: &rng,
                                xi,
                                omega,
                                alpha)
            case .noncentralChiSquare(k: let k, lambda: let lambda):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, k: T, lambda: T) -> T in
                    rncchisq(k, lambda, &rng)
                },
                                rng: &rng,
                                k,
                                lambda)
            case .noncentralT(mu: let mu, sigma: let sigma, nu: let nu, delta: let delta):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, mu: T, sigma: T, nu: T, delta: T) -> T in
                    rnct(mu, sigma, nu, delta, &rng)
                },
                                rng: &rng,
                                mu,
                                sigma,
                                nu,
                                delta)
            case .sasPoint5(location: let location, sigma: let sigma):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator,location: T, sigma: T) -> T in
                                rsasPoint5(mu: location, sigma: sigma, &rng)
                    },
                                rng: &rng,
                                location,
                                sigma)
            case .holtsmark(mu: let mu, c: let c):
                return rsampler(n,
                                using:
                                    { (rng: inout SystemRandomNumberGenerator, mu: T, c: T) -> T in
                    rholtsmark(mu, c, &rng)
                },
                                rng: &rng,
                                mu,
                                c)
            }
        }
    }
    
    /// Internal helper to loop `n` times calling a distribution-specific generator.
    ///
    /// - Parameters:
    ///   - n: Count of samples.
    ///   - f: Sampling closure accepting the RNG and distribution parameters, returning one `T`.
    ///   - rng: Mutable random generator.
    ///   - params: Distribution parameters to pass through.
    /// - Returns: Array of `n` draws.
    private static func rsampler<RNG: RandomNumberGenerator, each P>(_ n: Int,
                                                                     using f: (inout RNG, repeat each P) -> T,
                                                                     rng: inout RNG,
                                                                     _ params: repeat each P) -> [T] {
        precondition(n > 0)
        var out:[T] = []
        out.reserveCapacity(n)
        for _ in 0..<n {
            out.append(f(&rng, repeat each params))
        }
        return out
    }
    
    // MARK: - Specialized samplers (one variate each)
    
    /// Landau variate via inverse-CDF sampling using SwiftyBoost quantiles.
    ///
    /// - Parameters:
    ///   - mu: Location.
    ///   - c: Scale (> 0).
    ///   - rng: Random generator.
    /// - Returns: One Landau draw.
    public static func rlandau (
        mu: T, c: T, _ rng: inout some RandomNumberGenerator
    ) -> T {
        let u = T.random(in: 0..<1, using: &rng)
        return try! SwiftyBoost.Distribution.Landau<T>(location: mu, scale: c).quantile(u)
    }
    
    /// Map Airy variate via the distribution’s quantile function (inverse CDF).
    public static func rmapAiry(
        mu: T, sigma: T, _ rng: inout some RandomNumberGenerator
    ) -> T {
        let u = T.random(in: 0..<1, using: &rng)
        return try! SwiftyBoost.Distribution.MapAiry(location: mu, scale: sigma).quantile(u)
    }
    
    /// Triangular distribution sampler with parameters (lower, upper, mode) using inverse CDF.
    ///
    /// - Precondition: `lower < mode < upper`.
    public static func rtriangular(
        lower: T, upper: T, mode: T, _ rng: inout some RandomNumberGenerator
    ) -> T {
        precondition(lower < mode && mode < upper)
        let u = T.random(in: 0..<1, using: &rng)
        let Fc = (mode - lower) / (upper - lower)
        if u < Fc {
            return lower + T.sqrt(u * (upper - lower) * (mode - lower))
        } else {
            return upper - T.sqrt((T.one - u) * (upper - lower) * (upper - mode))
        }
    }
    
    /// Uniform(0,1) variate.
    ///
    /// - Returns: A single `U(0,1)` draw.
    private static func u01(_ rng: inout some RandomNumberGenerator) -> T {
        return T.random(in: 0..<1, using: &rng)
    }
    
    /// Bernoulli(p) variate (0/1) via inverse transform.
    ///
    /// - Precondition: `0 ≤ p ≤ 1`.
    private static func rbernoulli(_ p: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(p >= .zero && p <= T.one, "Bernoulli probability must lie in [0,1]")
        return u01(&rng) < p ? T.one : .zero
    }
    
    /// Standard normal N(0,1) via Marsaglia polar (Box–Muller polar form).
    private static func rnorm01(_ rng: inout some RandomNumberGenerator) -> T {
        var u: T = .zero
        var v: T = .zero
        var s: T = .zero
        repeat {
            let u01a = u01(&rng)
            let u01b = u01(&rng)
            u = T(2) * u01a - T.one
            v = T(2) * u01b - T.one
            s = u*u + v*v
        } while s >= T.one || s == .zero
        let factor = T.sqrt(-T(2) * T.log(s) / s)
        return u * factor
    }
    
    // ---------- Gamma(shape k, scale θ) : Marsaglia–Tsang ----------
    
    /// Gamma(k, θ) variate using the Marsaglia–Tsang method with Johnk’s trick for k < 1.
    ///
    /// - Precondition: `k > 0`, `theta > 0`.
    private static func rgamma(_ k: T, _ theta: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(k > 0 && theta > 0)
        if k < 1 {
            // Johnk’s trick: Gamma(k,θ) = Gamma(k+1,θ) * U^(1/k)
            let x = rgamma(k + 1, theta, &rng)
            let u = T.maximum(u01(&rng), .leastNonzeroMagnitude)
            return x * T.pow(u, 1.0 / k)
        }
        // k >= 1
        let d = k - T.one/T(3)
        let c = T.one / T.sqrt(T(9) * d)
        while true {
            let z = rnorm01(&rng)
            let v = T.pow(1.0 + c * z, 3)
            if v <= 0 { continue }
            let u = u01(&rng)
            let z2 = z * z
            if u < 1.0 - 0.0331 * z2 * z2 { return theta * d * v }
            if T.log(u) < T.half * z2 + d * (T.one - v + T.log(v)) { return theta * d * v }
        }
    }
    
    /// Chi-square variate via Gamma(df/2, 2).
    private static func rchisq(_ df: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(df > 0)
        return rgamma(df / T.two, T.two, &rng)
    }
    
    /// Fisher’s F variate via ratio of scaled chi-squares.
    private static func rf(_ d1: T, _ d2: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(d1 > 0 && d2 > 0)
        let x1 = rchisq(d1, &rng) / d1
        let x2 = rchisq(d2, &rng) / d2
        return x1 / x2
    }
    
    /// Beta(α, β) variate via two Gammas.
    private static func rbeta(_ a: T, _ b: T, _ rng: inout some RandomNumberGenerator) -> T {
        let g1 = rgamma(a, T.one, &rng)
        let g2 = rgamma(b, T.one, &rng)
        return g1 / (g1 + g2)
    }
    
    /// Arcsine distribution on (a, b) via Beta(1/2, 1/2) scaling.
    private static func rarcsine(_ a: T, _ b: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(b > a, "Arcsine support requires b > a")
        let y = rbeta(T.half, T.half, &rng)
        return a + (b - a) * y
    }
    
    /// Exponential(rate) variate using inverse transform.
    private static func rexp(_ rate: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(rate > 0)
        let u = T.maximum(u01(&rng), .leastNonzeroMagnitude)
        return -T.log(u) / rate
    }
    
    /// Weibull(k, λ) variate via inversion.
    private static func rweibull(_ k: T, _ lambda: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(k > 0 && lambda > 0)
        let u = T.maximum(T.one - u01(&rng), .leastNonzeroMagnitude)
        return lambda * T.pow(-T.log(u), T.one / k)
    }
    
    /// Normal(μ, σ) variate using a standard normal and affine transform.
    private static func rnorm(_ mu: T, _ sigma: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(sigma >= 0)
        return mu + sigma * rnorm01(&rng)
    }
    
    /// Student’s t variate with location, scale, and degrees of freedom.
    private static func rt(_ mu: T, _ sigma: T, _ nu: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(sigma >= 0 && nu > 0)
        let z = rnorm01(&rng)
        let v = rchisq(nu, &rng)
        return mu + sigma * z / T.sqrt(v / nu)
    }
    
    /// Cauchy(μ, γ) variate via inverse CDF.
    private static func rcauchy(_ mu: T, _ gamma: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(gamma > 0)
        let eps = T.leastNonzeroMagnitude
        let uRaw = u01(&rng)
        let u = min(max(uRaw, eps), T.one - eps)
        let centered = u - T.half
        return mu + gamma * T.tan(T.pi * centered)
    }
    
    /// Logistic(μ, s) variate via inverse CDF.
    private static func rlogistic(_ mu: T, _ s: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(s > 0)
        let eps = T.leastNonzeroMagnitude
        let uRaw = u01(&rng)
        let u = min(max(uRaw, eps), T.one - eps)
        return mu + s * T.log(u / (T.one - u))
    }
    
    /// Gumbel (extreme value type I) variate via inverse CDF.
    private static func rgumbel(_ mu: T, _ beta: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(beta > 0)
        let u = T.maximum(u01(&rng), .leastNonzeroMagnitude)
        return mu - beta * T.log(-T.log(u))
    }
    
    /// Generalized Extreme Value (GEV) variate via inverse CDF, with ξ ≈ 0 handled by Gumbel limit.
    private static func rgev(_ mu: T, _ sigma: T, _ xi: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(sigma > 0)
        let u = T.maximum(u01(&rng), .leastNonzeroMagnitude)
        if xi.magnitude < 1e-10 {
            return mu - sigma * T.log(-T.log(u))
        } else {
            return mu + (sigma/xi) * (T.pow(-T.log(u), -xi) - T.one)
        }
    }
    
    /// Generalized Pareto Distribution (GPD) variate via inverse CDF, with ξ ≈ 0 handled by exponential limit.
    private static func rgpd(_ sigma: T, _ xi: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(sigma > 0)
        let u = T.maximum(u01(&rng), .leastNonzeroMagnitude)
        if xi.magnitude < 1e-10 {
            return -sigma * T.log(1.0 - u)
        } else {
            return sigma * (T.pow(T.one - u, -xi) - T.one) / xi
        }
    }
    
    /// Skew-Normal(ξ, ω, α) variate using Azzalini’s construction.
    private static func rskewNormal(_ xi: T, _ omega: T, _ alpha: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(omega > 0)
        let delta = alpha / T.sqrt(T.one + alpha*alpha)
        let z0 = rnorm01(&rng)
        let z1 = rnorm01(&rng)
        let y  = delta * z0.magnitude + T.sqrt(T.one - delta * delta) * z1
        return xi + omega * y
    }
    
    /// Noncentral chi-square variate via Poisson mixture with a large-λ acceleration.
    ///
    /// - For small λ/2, samples N ~ Poisson(λ/2) exactly (Knuth) then returns χ²(k + 2N).
    /// - For large λ/2, uses Normal approximation for the Poisson count to avoid long loops.
    private static func rncchisq(_ k: T, _ lambda: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(k > 0 && lambda >= 0)
        let poisMean = lambda / T.two
        
        if poisMean > T(50) {
            let z = rnorm01(&rng)
            var nReal = poisMean + T.sqrt(poisMean) * z
            if nReal < 0 { nReal = 0 }
            let nInt = Int(nReal.rounded())
            let df = k + T.two * T(nInt)
            return rchisq(df, &rng)
        } else {
            let L = T.exp(-poisMean)
            var p = T.one
            var n = 0
            repeat {
                n += 1
                p *= u01(&rng)
            } while p > L
            let df = k + T.two * T(n - 1)
            return rchisq(df, &rng)
        }
    }
    
    /// Noncentral Student’s t variate via shifted normal and scaled chi-square.
    private static func rnct(_ mu: T, _ sigma: T, _ nu: T, _ delta: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(sigma >= 0 && nu > 0)
        let z = rnorm01(&rng) + delta
        let v = rchisq(nu, &rng)
        let t = z / T.sqrt(v / nu)
        return mu + sigma * t
    }
    
    /// SAS 0.5 (symmetric α-stable with α=0.5) variate using the inverse CDF via SwiftyBoost quantiles.
    public static func rsasPoint5(
        mu: T, sigma: T, _ rng: inout some RandomNumberGenerator
    ) -> T {
        let u = T.random(in: 0..<1, using: &rng)
        return try! SwiftyBoost.Distribution.SASPoint5(location: mu, scale: sigma).quantile(u)
    }
    
    /// Holtsmark variate using a Chambers–Mallows–Stuck–style construction specialized for α = 3/2.
    ///
    /// - Parameters:
    ///   - mu: Location.
    ///   - c: Scale (> 0).
    /// - Returns: One Holtsmark draw.
    private static func rholtsmark(_ mu: T, _ c: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(c > T.zero, "scale must be > 0")
        let alpha: T = 1.5
        let invAlpha: T = T.one / alpha
        let exp2 = (T.one - alpha) * invAlpha        // = -1/3
        
        while true {
            let V = T.pi * (u01(&rng) - 0.5)
            let cosV = T.cos(V)
            if cosV.magnitude < 1e-12 { continue }
            
            // W ~ Exp(1)
            let U = T.maximum(u01(&rng), .leastNonzeroMagnitude)
            let W = -T.log(U)
            
            let s1 = T.sin(alpha * V) / T.pow(cosV, invAlpha)
            let c2 = T.maximum(T.cos((T.one - alpha) * V).magnitude, 1e-300) * (T.cos((T.one - alpha) * V) >= T.zero ? T.one : -T.one)
            let s2 = T.pow(c2 / W, exp2)
            
            let S = s1 * s2
            return mu + c * S
        }
    }
    
    /// Binomial variate drawn as a sum of IID Bernoulli trials.
    ///
    /// - Parameters:
    ///   - n: Number of trials (≥ 0).
    ///   - p: Success probability in [0, 1].
    ///   - rng: A generator conforming to ``RandomNumberGenerator``
    /// - Returns: Count of successes as `T`.
    public static func rbinomial(
        numberOfTrials n: Int,
        probabilityOfSuccess p: T,
        _ rng: inout some RandomNumberGenerator
    ) -> T {
        precondition(n >= 0, "numberOfTrials must be >= 0")
        precondition(p >= .zero && p <= T.one, "probabilityOfSuccess must lie in [0,1]")
        
        if n == 0 { return .zero }
        if p == .zero { return .zero }
        if p == T.one { return T(n) }
        
        var count: Int = 0
        for _ in 0..<n {
            if u01(&rng) < p { count &+= 1 }
        }
        return T(count)
    }
    
    /// Log-normal variate parameterized by log-mean (`meanlog`) and log-sd (`sdlog`),
    /// using the classical Box–Muller transform.
    ///
    /// - Precondition: `sdlog > 0`.
    public static func rlognormal(
        meanlog mu: T,
        sdlog sigma: T,
        _ rng: inout some RandomNumberGenerator
    ) -> T {
        precondition(sigma > .zero && sigma.isFinite, "sdlog must be > 0 and finite")
        
        // Draw a standard normal using Box-Muller (classical form)
        var r = rng
        var u1: T = 0.0
        var u2: T = 0.0
        repeat {
            u1 = T.random(in: T.zero..<T.one, using: &r)
        } while u1 <= 0.0
        u2 = T.random(in: T.zero..<T.one, using: &r)
        
        let z = T.sqrt(-T.two * T.log(u1)) * T.cos(T.two * T.pi * u2) // N(0,1)
        
        rng = r
        let val = T.exp(mu + sigma * z)
        return val
    }
    
    /// Pareto(scale, shape) variate via inverse transform sampling.
    public static func rpareto(
        scale: T, shape: T, _ rng: inout some RandomNumberGenerator
    ) -> T {
        precondition(scale > .zero && shape > .zero)
        let u = T.random(in: 0..<1, using: &rng)
        // inverse CDF: x = xm * (1 - u)^(-1/alpha)
        return scale * T.pow(T.one - u, -T.one/shape)
    }
    
    /// Poisson sampler using Knuth’s multiplication algorithm (accurate for small λ).
    ///
    /// - Precondition: `lambda ≥ 0`.
    /// - Note: For large λ, consider switching to a more efficient sampler if needed.
    public static func rpoisson(
        lambda: T, _ rng: inout some RandomNumberGenerator
    ) -> T {
        precondition(lambda >= .zero && lambda.isFinite)
        if lambda == .zero { return .zero }
        let L = T.exp(-lambda)
        var p = T.one
        var k: Int = 0
        var r = rng
        while p > L {
            k &+= 1
            p *= T.random(in: 0..<1, using: &r)
        }
        rng = r
        return T(k - 1)
    }
    
    /// Rayleigh sampler: X = σ sqrt(-2 log(1 - U)).
    public static func rrayleigh(
        sigma: T, _ rng: inout some RandomNumberGenerator
    ) -> T {
        precondition(sigma > .zero)
        let u = T.random(in: 0..<1, using: &rng)
        return sigma * T.sqrt(-T(2) * T.log(T.one - u))
    }
    
    /// Geometric(p) sampler returning the number of failures before the first success (support {0,1,...}).
    ///
    /// - Precondition: `0 < p < 1`.
    /// - Method: Inverse CDF.
    public static func rgeometric(
        p: T, _ rng: inout some RandomNumberGenerator
    ) -> T {
        precondition(p > .zero && p < T.one)
        // Inverse CDF: k = floor( log(1-U)/log(1-p) )
        let u = T.random(in: 0..<1, using: &rng)
        return T(floor(T.log(T.one - u) / T.log(T.one - p)))
    }
    
    /// Inverse Gaussian (Wald) sampler using the Michael–Schucany–Haas acceptance rule.
    ///
    /// - Parameters:
    ///   - mu: Mean (> 0).
    ///   - lambda: Shape (> 0).
    ///   - rng: A generator conforming to ``RandomNumberGenerator``
    /// - Returns: One Wald draw.
    public static func rinorm(_ mu: T, lambda: T, _ rng: inout some RandomNumberGenerator) -> T {
        precondition(mu > .zero && mu.isFinite, "mu must be > 0 and finite")
        precondition(lambda > .zero && lambda.isFinite, "lambda must be > 0 and finite")

        // 1) Draw V ~ N(0,1) (Marsaglia polar method to avoid trig)
        let V = rnorm01(&rng)
        let Y = V*V

        // 2) Candidate root
        //    X = mu + (mu^2 * Y)/(2*lambda) - (mu/(2*lambda)) * sqrt(4*mu*lambda*Y + mu^2 * Y^2)
        let mu2 = mu*mu
        let lam2 = T(2) * lambda
        let term = T.sqrt(T(4) * mu * lambda * Y + mu2 * Y * Y)
        let X = mu + (mu2 * Y)/lam2 - (mu/lam2) * term

        // 3) Accept with probability mu/(mu+X), else return mu^2 / X
        let U = u01(&rng)
        if U <= mu / (mu + X) {
            return X
        } else {
            return mu2 / X
        }
    }
    
    /// Negative-binomial variate sampled as the sum of `r` IID geometric draws.
    ///
    /// - Parameters:
    ///   - r: Number of successes (> 0).
    ///   - p: Success probability in (0, 1).
    ///   - rng: A generator conforming to ``RandomNumberGenerator``
    /// - Returns: Count of failures before the r-th success.
    public static func rnegativeBinomial(
        r: Int, p: T, _ rng: inout some RandomNumberGenerator
    ) -> T {
        precondition(r > 0 && p > .zero && p < T.one)
        var total: T = .zero
        var numbers: [T] = []
        numbers.reserveCapacity(r)
        for _ in 0..<r {
            numbers.append(rgeometric(p: p, &rng))
        }
        total = Helpers.sum(&numbers)
        return total
    }
    public static func rnegativeBinomial(
        r: T, p: T, _ rng: inout some RandomNumberGenerator
    ) -> T {
        precondition(r > .zero)
        precondition(p > .zero && p < T.one)
        let scale = (T.one - p) / p
        let lambda = rgamma(r, scale, &rng)
        return rpoisson(lambda: lambda, &rng)
    }
}

