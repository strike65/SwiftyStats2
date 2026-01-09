//
//  Created by VT on 08.11.25.
//  © 2025 Volker Thieme. All rights reserved.
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

/// Maximum likelihood fitter for the NormalSkew distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {
    /// Fit a Skew-Normal(μ, ω, α) distribution using numerical MLE with gradient support and robust warm-starts.
    ///
    /// Parameterization: location μ ∈ ℝ, scale ω > 0, shape α ∈ ℝ.
    ///
    /// - Parameters:
    ///   - data: Sample values.
    ///   - optimizer: Optimizer to use (default: `.lbfgs`).
    ///   - options: Optional optimizer configuration. If `warmStartTheta` is not set, a method-of-moments warm start is used.
    /// - Returns: `MLEResult` with `thetaHat = [μ̂, ω̂, α̂]`.
    /// - Discussion:
    ///   - Uses Azzalini identities to construct initial guesses from sample mean, variance, and skewness.
    ///   - Provides analytic score for improved convergence when the optimizer supports gradients.
    ///   - Uses a small grid of candidates to robustify warm-start selection.
    public static func fitSkewNormal(
        _ data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)

        // Parameter specs for Skew-Normal.
        let dummy = MLEProblem<T>(data: data, logpdf: { _, _ in 0 }, paramSpecs: [.init(.real, initial: 0)])
        let specs = dummy.makeParamSpecs(for: .skew_normal, data: data)

        // Method-of-moments warm start from sample mean, variance, skewness.
        let n = data.count
        let nT = T(n)

        // Sample mean.
        let mean: T = data.reduce(0, +) / nT

        // Population variance (divide by n) for method-of-moments.
        var s2acc: T = 0
        for x in data {
            let d = x - mean
            s2acc += d * d
        }
        let varPop = s2acc / nT
        let vSample = max(varPop, T(1e-12))

        // Population skewness (divide by n).
        var m3acc: T = 0
        for x in data {
            let d = x - mean
            m3acc += d * d * d
        }
        let skewSample: T = {
            let s = T.sqrt(varPop)
            if s == 0 { return 0 }
            return (m3acc / nT) / (s * s * s)
        }()

        // Azzalini’s identities relating δ and α, and moments.
        func skewnessForDelta(_ delta: T) -> T {
            let t = delta * T.sqrt(T.two / T.pi)
            let num = (T(4) - T.pi) / T.two * (t * t * t)
            let den = T.pow(T.one - t * t, T(3) / T.two)
            return num / den
        }

        // Invert γ1(δ) ≈ target via bisection on |δ| with sign restoration.
        func deltaFromSkewness(_ g1: T) -> T {
            if g1 == .zero { return .zero }
            let sgn: T = g1 > 0 ? T.one : -T.one
            let target = abs(g1)
            var lo: T = 0
            var hi: T = T(0.999) // avoid singularity near |delta| = 1
            for _ in 0..<60 {
                let mid = (lo + hi) / T.two
                let val = abs(skewnessForDelta(mid))
                if val >= target {
                    hi = mid
                } else {
                    lo = mid
                }
            }
            return sgn * (lo + hi) / T.two
        }

        // Construct μ0, ω0, α0 from moments.
        let delta0 = deltaFromSkewness(skewSample)
        let denomVar = max(T.one - (T.two * delta0 * delta0 / T.pi), T(1e-6))
        let omega0 = T.sqrt(vSample / denomVar)
        let mu0 = mean - omega0 * delta0 * T.sqrt(T.two / T.pi)
        let alpha0: T = {
            let d = max(min(delta0, T(0.999999)), T(-0.999999))
            return d / T.sqrt(T.one - d * d)
        }()

        // Log-pdf and gradient (score) for Skew-Normal.
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let mu    = theta[0]
            let omega = theta[1]
            let alpha = theta[2]
            if !(mu.isFinite && omega.isFinite && alpha.isFinite && omega > 0) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.SkewNormal(location: mu, scale: omega, shape: alpha).logPdf(x)
            } catch {
                return T.infinity
            }
        }

        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let mu    = theta[0]
            let sigma = theta[1]
            let alpha = theta[2]
            guard mu.isFinite, sigma.isFinite, alpha.isFinite, sigma > 0 else {
                return [T.nan, T.nan, T.nan]
            }
            do {
                let sn: SwiftyBoost.Distribution.SkewNormal<T> = try SwiftyBoost.Distribution.SkewNormal<T>(location: mu, scale: sigma, shape: alpha)
                let g = sn.score(x: x, mu: mu, sigma: sigma, alpha: alpha)
                return [g.dmu, g.dsigma, g.dalpha]
            }
            catch _ {
                return [T.nan, T.nan, T.nan]
            }
        }

        // Assemble problem.
        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)

        // Build a small grid of warm-start candidates for robustness across shapes.
        func skewNormalWarmStarts(_ base: [T], grid: [T]) -> [[T]] {
            var starts: [[T]] = []
            let mu0b = base[0], om0b = base[1], a0b = base[2]
            for g in grid {
                starts.append([mu0b, om0b, a0b + g])
            }
            return starts
        }

        // Pick the best warm start by scanning candidates and selecting the highest feasible log-likelihood.
        func bestWarmStart() -> [T] {
            let base = specs.map { $0.initial }
            var candidates: [[T]] = []
            let grid: [T] = [-3, -1.5, 0, 1.5, 3]
            candidates.append(contentsOf: skewNormalWarmStarts(base, grid: grid))

            // Moment-based start and nearby shape perturbations.
            let mom = [mu0, max(omega0, T(1e-6)), alpha0]
            candidates.append(mom)
            for g in grid {
                candidates.append([mu0, max(omega0, T(1e-6)), alpha0 + g * T(0.25)])
            }

            // Evaluate log-likelihood across candidates; keep the best feasible.
            var best: [T] = candidates.first!
            var bestLL: T = -.infinity
            for th in candidates {
                var s: T = 0
                var valid = true
                for x in data {
                    let lp = logPDF(x, th)
                    if !lp.isFinite { valid = false; break }
                    s += lp
                }
                if valid && s > bestLL {
                    bestLL = s
                    best = th
                }
            }
            return best
        }

        // Optimizer defaults tuned for a moderately challenging, 3-parameter likelihood.
        let opts: MLEOptimizationOpts<T> = {
            if var o = options {
                o.optimizer = optimizer
                o.computeCovariance = o.computeCovariance
                if o.warmStartTheta == nil {
                    o.warmStartTheta = bestWarmStart()
                }
                if o.gradStep < T(1e-5) { o.gradStep = T(1e-5) }
                o.diagnosticsEnabled = true
                return o
            } else {
                var o = MLEOptimizationOpts<T>(); o.computeCovariance = true
                o.optimizer = optimizer
                o.multiStartCount = 12
                o.multiStartDesign = .lhs
                o.randomRestartCount = 3
                o.initialStepStrategy = .relativeTheta(T(0.35))
                o.gradStep = T(1e-5)
                o.hessianStep = T(1e-3)
                o.relTolLogLik = T(1e-7)
                o.diagnosticsEnabled = true
                o.warmStartTheta = bestWarmStart()
                return o
            }
        }()
        return MLESolver.fit(problem, options: opts)
    }
}
