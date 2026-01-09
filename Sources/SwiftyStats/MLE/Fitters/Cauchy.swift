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

/// Maximum likelihood fitter for the Cauchy distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {

    /// Fit a Cauchy(μ, γ) distribution using numerical MLE.
    ///
    /// Parameterization: location μ ∈ ℝ, scale γ > 0.
    ///
    /// - Parameters:
    ///   - data: Sample values.
    ///   - optimizer: Optimizer to use (default: `.nelderMead`).
    ///   - options: Optional optimizer configuration.
    /// - Returns: `MLEResult` with `thetaHat = [μ̂, γ̂]`.
    public static func fitCauchy(
        _ data: [T],
        optimizer: OptimizerKind = .nelderMead,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)

        // Parameter specs for Cauchy.
        let dummy = MLEProblem<T>(data: data, logpdf: {_,_ in 0}, paramSpecs: [.init(.real, initial: 0)])
        let specs = dummy.makeParamSpecs(for: .cauchy, data: data)

        // Log-pdf with domain checks.
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let mu = theta[0]
            let gamma = theta[1]
            if !(mu.isFinite && gamma.isFinite && gamma > 0) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.Cauchy(location: mu, scale: gamma).logPdf(x)
            } catch {
                return T.infinity
            }
        }

        // Analytic score for (μ, γ).
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let mu = theta[0]
            let s  = theta[1]
            guard mu.isFinite, s.isFinite, s > 0 else { return [T.nan, T.nan] }
            do {
                let dist: SwiftyBoost.Distribution.Cauchy<T> = try SwiftyBoost.Distribution.Cauchy(location: mu, scale: s)
                let g = dist.score(x: x, mu: mu, sigma: s)
                return [g.dmu, g.dsigma]
            } catch {
                return [T.nan, T.nan]
            }
        }

        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)

        // Default optimizer setup.
        let opts: MLEOptimizationOpts<T> = {
            if let o = options {
                return o
            } else {
                var o = MLEOptimizationOpts<T>(); o.computeCovariance = true
                o.optimizer = optimizer
                return o
            }
        }()
        return MLESolver.fit(problem, options: opts)
    }
}
