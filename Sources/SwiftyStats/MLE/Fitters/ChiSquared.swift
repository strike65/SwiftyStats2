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
/// Maximum likelihood fitter for the ChiSquared distribution in SwiftyStats.

import SwiftyStatsPrelude

extension MLEFitter {
    /// Fit a central Chi-Squared(k) distribution using numerical MLE.
    ///
    /// Parameterization: degrees of freedom k > 0 (unknown).
    ///
    /// - Parameters:
    ///   - data: Sample values (support x >= 0).
    ///   - optimizer: Optimizer to use (default: `.lbfgs`).
    ///   - options: Optional optimizer configuration.
    /// - Returns: `MLEResult` with `thetaHat = [kHat]`.
    public static func fitChiSquared(
        _ data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty, "data must be non-empty")

        // Optional start-context hook for your spec builder.
        let ctx = MLEStartContext<T>()

        // Parameter spec via factory: one positive parameter (k) with robust start.
        // Create a dummy problem just to access the spec factory extension.
        let dummy = MLEProblem<T>(data: data, logpdf: {_,_ in 0}, paramSpecs: [
            .init(.positive, initial: T.one)
        ])
        let specs = dummy.makeParamSpecs(for: .chiSquared, data: data, ctx: ctx)

        // Log-pdf with domain checks for support and parameter.
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let k = theta[0]
            // For invalid θ or x, return -∞ logpdf so NLL becomes +∞ (rejected by solver).
            if !(k.isFinite && k > 0) { return -T.infinity }
            if !(x.isFinite && x >= .zero) { return -T.infinity }
            do {
                return try SwiftyBoost.Distribution.ChiSquared(degreesOfFreedom: k).logPdf(x)
            } catch {
                // Treat construction errors as invalid evaluation.
                return -T.infinity
            }
        }

        // Analytic gradient using the implemented score d/dk log f(x; k).
        // On invalid inputs, return a finite vector (zeros) to avoid NaNs.
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let k = theta[0]
            guard k.isFinite, k > 0, x.isFinite, x >= .zero else { return [T.zero] }
            do {
                let dist: SwiftyBoost.Distribution.ChiSquared<T> =
                    try SwiftyBoost.Distribution.ChiSquared(degreesOfFreedom: k)
                // Score API: score(x:k) -> T
                let g = dist.score(x: x, k: k)
                return [g]
            } catch {
                return [T.zero]
            }
        }

        let problem = MLEProblem<T>(
            data: data,
            logpdf: logPDF,
            gradlogpdf: gradLogPDF,
            paramSpecs: specs
        )

        // Robust defaults for a single-parameter likelihood.
        let opts: MLEOptimizationOpts<T> = {
            if let o = options { return o }
            var o = MLEOptimizationOpts<T>()
            o.computeCovariance = true
            o.optimizer = optimizer
            o.multiStartCount = 6
            o.multiStartDesign = .lhs
            o.randomRestartCount = 2
            o.initialStepStrategy = .relativeTheta(T(0.25))
            o.gradStep = T(1e-5)
            o.hessianStep = T(5e-4)
            o.relTolLogLik = T(1e-7)
            o.diagnosticsEnabled = true
            return o
        }()

        return MLESolver.fit(problem, options: opts)
    }
}

