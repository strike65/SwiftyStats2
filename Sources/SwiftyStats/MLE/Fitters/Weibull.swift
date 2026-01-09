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

/// Maximum likelihood fitter for the Weibull distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {

    /// Fit a Weibull(k, λ) distribution using a hybrid approach:
    /// Newton iteration for shape k followed by numerical optimization if needed.
    ///
    /// Parameterization: shape k > 0, scale λ > 0.
    ///
    /// - Parameters:
    ///   - data: Sample values. Must be strictly positive.
    ///   - optimizer: Optimizer to use (default: `.lbfgs`).
    ///   - options: Optional optimizer configuration. If not provided, a robust default with warm-start and diagnostics is used.
    /// - Returns: `MLEResult` with `thetaHat = [k̂, λ̂]`.
    /// - Discussion:
    ///   - Uses Newton updates for `k` based on log-transformed data and moment identities.
    ///   - If Newton converges, returns the analytic log-likelihood at `[k̂, λ̂]`.
    ///   - Otherwise, falls back to numerical MLE with warm-start.
    public static func fitWeibull(
        _ data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)
        precondition(data.allSatisfy { $0 > .zero }, "Weibull MLE requires strictly positive data")

        // Precompute log(x) and its mean for Newton updates on k.
        let n = data.count
        let nT = T(n)
        let lnX: [T] = data.map { T.log($0) }
        let meanLnX: T = lnX.reduce(0, +) / nT

        // Initial shape from parameter specs; clamp to a reasonable range.
        let dummy = MLEProblem<T>(data: data, logpdf: {_,_ in 0}, paramSpecs: [.init(.real, initial: 0)])
        let baseSpecs = dummy.makeParamSpecs(for: .weibull, data: data)
        var k = max(baseSpecs[0].initial, T(1e-3))
        k = min(k, T(100))

        // Newton refinement for shape k.
        // f(k) = 1/k + E[log X] - E_w[log X], where weights are proportional to X^k.
        var it = 0
        var convergedNewton = false
        let maxIt = 60
        while it < maxIt {
            it += 1

            // Weighted moments with weights w_i ∝ X_i^k.
            var S0: T = 0, S1: T = 0, S2: T = 0
            for i in 0..<n {
                let t = T.exp(k * lnX[i])
                S0 += t
                let lx = lnX[i]
                S1 += t * lx
                S2 += t * lx * lx
            }
            if !(S0.isFinite && S1.isFinite && S2.isFinite) { break }

            // Mean and variance of log X under weights.
            let m2 = S1 / S0
            let varW = max(S2 / S0 - m2 * m2, T.zero)

            // Newton step: k_{new} = k - f / f'.
            let f = T.one / k + meanLnX - m2
            let fp = -T.one / (k * k) - varW
            if !fp.isFinite || fp == .zero { break }

            let step = f / fp
            var kNew = k - step
            if !(kNew.isFinite) { break }

            // Guard against collapse or explosion.
            if kNew <= T(1e-6) { kNew = T(1e-6) }
            if kNew > T(1e6) { kNew = T(1e6) }

            // Relative convergence check on k.
            let rel = abs(kNew - k) / max(T.one, abs(k))
            k = kNew
            if rel < T(1e-10) {
                convergedNewton = true
                break
            }
        }

        // Scale estimate λ̂ from k̂ via mean of X^k: λ̂ = (E[X^k])^{1/k}.
        var sumXk: T = 0
        for i in 0..<n {
            sumXk += T.exp(k * lnX[i])
        }
        let lambdaHat = T.exp((T.log(sumXk / nT)) / k)

        // Log-likelihood for Weibull(k, λ) for short-circuit return.
        func logLikWeibull(_ k: T, _ lam: T) -> T {
            let nT = T(n)
            let sumLnX = lnX.reduce(0, +)
            let lnLam = T.log(lam)
            var sPow: T = 0
            for lx in lnX {
                sPow += T.exp(k * (lx - lnLam))
            }
            return nT * T.log(k) - nT * k * lnLam + (k - T.one) * sumLnX - sPow
        }

        // Warm-start for the numerical fallback.
        let warmTheta = [k, lambdaHat]

        // Parameter specs and distribution closures.
        let specs = baseSpecs
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let kk = theta[0]
            let lam = theta[1]
            if !(kk.isFinite && lam.isFinite && kk > 0 && lam > 0) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.Weibull(shape: kk, scale: lam).logPdf(x)
            } catch {
                return T.infinity
            }
        }
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let kk = theta[0]
            let lam = theta[1]
            guard kk.isFinite, lam.isFinite, kk > 0, lam > 0 else { return [T.nan, T.nan] }
            do {
                let dist: SwiftyBoost.Distribution.Weibull<T> = try SwiftyBoost.Distribution.Weibull(shape: kk, scale: lam)
                let g = dist.score(x: x, k: kk, lambda: lam)
                return [g.dk, g.dlambda]
            } catch {
                return [T.nan, T.nan]
            }
        }
        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)

        // If Newton converged sufficiently, short-circuit and return.
        if convergedNewton {
            let ll = logLikWeibull(k, lambdaHat)
            return MLEResult(thetaHat: [k, lambdaHat],
                             logLik: ll,
                             iterations: it,
                             converged: true,
                             nEval: 0,
                             cov: nil)
        }

        // Otherwise, fall back to numerical MLE with warm-start and reasonable defaults.
        let opts: MLEOptimizationOpts<T> = {
            if var o = options {
                if o.warmStartTheta == nil { o.warmStartTheta = warmTheta }
                o.computeCovariance = o.computeCovariance
                o.optimizer = optimizer
                return o
            } else {
                var o = MLEOptimizationOpts<T>()
                o.computeCovariance = true
                o.optimizer = optimizer
                o.multiStartCount = 6
                o.multiStartDesign = .lhs
                o.randomRestartCount = 2
                o.initialStepStrategy = .relativeTheta(T(0.35))
                o.gradStep = T(5e-6)
                o.hessianStep = T(1e-3)
                o.relTolLogLik = T(1e-7)
                o.diagnosticsEnabled = true
                o.warmStartTheta = warmTheta
                return o
            }
        }()
        return MLESolver.fit(problem, options: opts)
    }
}
