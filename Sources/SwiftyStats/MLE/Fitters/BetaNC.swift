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

/// Maximum likelihood fitter for the BetaNC distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {
    /// Fit a Noncentral Beta(α, β, λ) distribution using numerical MLE.
    ///
    /// Parameterization: α > 0, β > 0, λ ≥ 0.
    ///
    /// - Parameters:
    ///   - data: Sample values.
    ///   - optimizer: Optimizer to use (default: `.lbfgs`).
    ///   - options: Optional optimizer configuration.
    /// - Returns: `MLEResult` with `thetaHat = [α̂, β̂, λ̂]`.
    /// - Implementation notes:
    ///   - This is a 3-parameter noncentral likelihood; defaults use multi-start to improve robustness.
    public static func fitNCBeta(
        _ data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)

        // Parameter specs for Noncentral Beta.
        let dummy = MLEProblem<T>(data: data, logpdf: { _, _ in 0 }, paramSpecs: [.init(.positive, initial: 1)])
        let specs = dummy.makeParamSpecs(for: .nc_beta, data: data)

        // Log-pdf with domain checks.
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let a   = theta[0]
            let b   = theta[1]
            let lam = theta[2]
            if !(a.isFinite && b.isFinite && lam.isFinite && a > 0 && b > 0 && lam >= 0) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.NonCentralBeta(alpha: a, beta: b, lambda: lam).logPdf(x)
            } catch {
                return T.infinity
            }
        }

        // Function-only optimization (no analytic gradient provided here).
        let problem = MLEProblem<T>(data: data, logpdf: logPDF, paramSpecs: specs)

        // Robust optimizer defaults for a 3-parameter noncentral likelihood.
        let opts: MLEOptimizationOpts<T> = {
            if let o = options {
                return o
            } else {
                var o = MLEOptimizationOpts<T>(); o.computeCovariance = true
                o.optimizer = optimizer
                o.multiStartCount = 10
                o.multiStartDesign = .lhs
                o.randomRestartCount = 2
                o.initialStepStrategy = .relativeTheta(T(0.30))
                o.gradStep = T(2e-5)
                o.hessianStep = T(5e-4)
                o.relTolLogLik = T(1e-7)
                o.diagnosticsEnabled = true
                return o
            }
        }()
        return MLESolver.fit(problem, options: opts)
    }
}
