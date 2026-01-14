//
//  Created by VT on 16.11.25.
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

/// Maximum likelihood fitter for the SASPoint5 distribution in SwiftyStats.

import SwiftyStatsPrelude

extension MLEFitter {
    /// Fit a Rayleigh distribution via its analytic MLE.
    ///
    /// Parameterization: variance-like parameter v = σ² (not the scale σ).
    /// The MLE is v̂ = (1/(2n)) * Σ x_i².
    ///
    /// - Parameters:
    ///   - data: Sample values. Must be non-empty and nonnegative.
    ///   - options: Optional optimizer options (unused for analytic fits).
    ///   - optimizer: ``OptimizerKind`` setting (defaults to ``OptimizerKind.lbfgs``).
    /// - Returns: `MLEResult` with `thetaHat = [v̂]` where `v̂ = (1/(2n)) * Σ x_i²` (the variance-like parameterization).
    /// - Precondition: `!data.isEmpty && data.allSatisfy { $0 >= 0 }`.
    /// - Note: If any `x_i == 0` or `Σ x_i² == 0`, the likelihood can be unbounded; returns `logLik = -∞` accordingly.
    public static func fitSASPoint5(
        data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
        
    ) -> MLEResult<T> {
        precondition(!data.isEmpty, "data must not be empty")

        let dummy = MLEProblem<T>(data: data, logpdf: { _, _ in 0 }, paramSpecs: [.init(.real, initial: 0)])
        let specs = dummy.makeParamSpecs(for: .saspoint5, data: data)

        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            do {
                let mu = theta[0]
                let scale = theta[1]
                if !(scale > .zero) { return T.infinity }
                let dist = try SwiftyBoost.Distribution.SASPoint5<T>(location: mu, scale: scale)
                let res = try dist.logPdf(x)
                return res
            }
            catch _ {
                return T.infinity
            }
        }

        let gradLogPdf: @Sendable (T, [T]) -> [T] = { x, theta in
            do {
                let mu = theta[0], scale = theta[1]
                if !(scale > T.zero) { return [T.nan, .nan] }
                let dist = try SwiftyBoost.Distribution.SASPoint5<T>(location: mu, scale: scale)
                let res = dist.score(x: x, mu: mu, sigma: scale)
                return [res.dmu, res.dsigma]
            }
            catch _ {
                return [T.nan, .nan]
            }
        }
        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPdf, paramSpecs: specs)
        var opts = options ?? MLEOptimizationOpts<T>()
        opts.computeCovariance = true
        opts.optimizer = optimizer
        opts.diagnosticsEnabled = true
        return MLESolver.fit(problem, options: opts)
    }
}
