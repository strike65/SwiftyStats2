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

/// Maximum likelihood fitter for the Triangular distribution in SwiftyStats.

import SwiftyStatsPrelude

extension MLEFitter {

    /// Fit a Triangular(lower, mode, upper) distribution via numerical MLE.
    ///
    /// Parameterization: lower bound a, upper bound b, and mode c with `a < c < b`.
    ///
    /// - Parameters:
    ///   - data: Sample values used for estimation.
    ///   - optimizer: Optimizer passed to the solver (default `.lbfgs`).
    ///   - options: Optional solver configuration (covariance and diagnostics enabled by default).
    /// - Returns: `MLEResult` with `thetaHat = [ĉ, â, b̂]` (matching the solver’s parameter ordering).
    public static func fitTriangular(
        _ data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)
        // Provide integer df to the start context if needed by the spec builder.
        var ctx = MLEStartContext<T>()
        ctx.minX = data.min()!
        ctx.maxX = data.max()!
        // Parameter specs for the Triangular distribution leverage the observed bounds.
        let dummy = MLEProblem<T>(data: data, logpdf: { _, _ in 0 }, paramSpecs: [.init(.positive, initial: 1)])
        let specs = dummy.makeParamSpecs(for: .triangular, data: data, ctx: ctx)

        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let a = theta[1], b = theta[2], c = theta[0]
            if !(a < c && c < b && x >= a && x <= b) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.Triangular(lower: a, mode: c, upper: b).logPdf(x)
            } catch {
                return T.infinity
            }
        }

        // Gradient (score) for ν.
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let a = theta[1], b = theta[2], c = theta[0]
            if !(a < c && c < b && x >= a && x <= b) { return [T.nan, T.nan, T.nan] }
            do {
                let dist: SwiftyBoost.Distribution.Triangular<T> = try SwiftyBoost.Distribution.Triangular(lower: a, mode: c, upper: b)
                let g = dist.score(x: x, lower: a, upper: b, mode: c)
                return [g.da, g.db, g.dc]
            } catch {
                return [T.nan,T.nan,T.nan]
            }
        }


        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)
        var opts = options ?? MLEOptimizationOpts<T>()
        opts.computeCovariance = true
        opts.optimizer = optimizer
        opts.diagnosticsEnabled = true
        return MLESolver.fit(problem, options: opts)
    }
}
