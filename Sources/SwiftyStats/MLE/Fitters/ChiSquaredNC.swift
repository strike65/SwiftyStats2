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

/// Maximum likelihood fitter for the ChiSquaredNC distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {
    /// Fit a Noncentral Chi-Squared(k, λ) distribution with fixed degrees of freedom using numerical MLE.
    ///
    /// Parameterization: fixed k > 0 (degrees of freedom), noncentrality λ ≥ 0.
    ///
    /// - Parameters:
    ///   - data: Sample values.
    ///   - k: Fixed degrees of freedom.
    ///   - optimizer: Optimizer to use (default: `.lbfgs`).
    ///   - options: Optional optimizer configuration.
    /// - Returns: `MLEResult` with `thetaHat = [λ̂]`.
    public static func fitNCChiSquared(
        _ data: [T],
        degreesOfFreedom k: T,
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)

        // Provide integer df to the start context if needed by the spec builder.
        var ctx = MLEStartContext<T>()
        ctx.ncChi_k = Helpers.integerValue(k.rounded(.toNearestOrAwayFromZero))

        // Parameter specs for Noncentral Chi-Squared with fixed df.
        let dummy = MLEProblem<T>(data: data, logpdf: { _, _ in 0 }, paramSpecs: [.init(.positive, initial: 1)])
        let specs = dummy.makeParamSpecs(for: .nc_chi_squared, data: data, ctx: ctx)

        // Log-pdf with domain checks for λ.
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let lam = theta[0]
            if !(lam.isFinite && lam >= 0) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.NonCentralChiSquared(degreesOfFreedom: k, nonCentrality: lam).logPdf(x)
            } catch {
                return T.infinity
            }
        }
        
        // Analytic gradient using the implemented score
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let lam = theta[0]
            guard lam.isFinite, lam > 0 else { return [T.nan, T.nan] }
            do {
                let dist: SwiftyBoost.Distribution.NonCentralChiSquared<T> = try SwiftyBoost.Distribution.NonCentralChiSquared(degreesOfFreedom: k, nonCentrality: lam)
                let g = dist.score(x: x, df: k, lambda: lam)
                return [k, g.dlambda]
            } catch {
                return [T.nan, T.nan]
            }
        }


        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)

        // Robust optimizer defaults for a 1-parameter noncentral likelihood.
        let opts: MLEOptimizationOpts<T> = {
            if let o = options { return o }
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
        }()
        return MLESolver.fit(problem, options: opts)
    }
}
