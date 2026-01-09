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

/// Maximum likelihood fitter for the FisherFNC distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {
    /// Fit a Noncentral F(d1, d2, λ) distribution with fixed degrees of freedom using numerical MLE.
    ///
    /// Parameterization: fixed df1 > 0, df2 > 0, noncentrality λ ≥ 0.
    ///
    /// - Parameters:
    ///   - data: Sample values.
    ///   - df1: Fixed numerator degrees of freedom.
    ///   - df2: Fixed denominator degrees of freedom.
    ///   - optimizer: Optimizer to use (default: `.lbfgs`).
    ///   - options: Optional optimizer configuration.
    /// - Returns: `MLEResult` with `thetaHat = [λ̂]`.
    public static func fitNCFisherF(
        _ data: [T],
        df1: T,
        df2: T,
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)

        // Provide integer dfs to the start context if needed by the spec builder.
        var ctx = MLEStartContext<T>()
        ctx.ncF_df1 = Helpers.integerValue(df1.rounded(.toNearestOrAwayFromZero))
        ctx.ncF_df2 = Helpers.integerValue(df2.rounded(.toNearestOrAwayFromZero))

        // Parameter specs for Noncentral F with fixed dfs.
        let dummy = MLEProblem<T>(data: data, logpdf: { _, _ in 0 }, paramSpecs: [.init(.positive, initial: 1)])
        let specs = dummy.makeParamSpecs(for: .nc_fisher_f, data: data, ctx: ctx)

        // Log-pdf with domain checks for λ.
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let lam = theta[0]
            if !(lam.isFinite && lam >= 0) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.NonCentralF(degreesOfFreedom1: df1, degreesOfFreedom2: df2, nonCentrality: lam).logPdf(x)
            } catch {
                return T.infinity
            }
        }

        // Noncentral F score not provided; use function-only optimization.
        let problem = MLEProblem<T>(data: data, logpdf: logPDF, paramSpecs: specs)

        // Robust optimizer defaults for a 1-parameter noncentral likelihood.
        let opts: MLEOptimizationOpts<T> = {
            if let o = options {
                return o
            } else {
                var o = MLEOptimizationOpts<T>()
                o.computeCovariance = true
                o.optimizer = optimizer
                o.multiStartCount = 8
                o.multiStartDesign = .lhs
                o.randomRestartCount = 2
                o.initialStepStrategy = .relativeTheta(T(0.25))
                o.gradStep = T(1e-5)
                o.hessianStep = T(5e-4)
                o.relTolLogLik = T(1e-7)
                o.diagnosticsEnabled = true
                return o
            }
        }()
        return MLESolver.fit(problem, options: opts)
    }
}
