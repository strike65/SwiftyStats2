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
/// Maximum likelihood fitter for the Landau distribution in SwiftyStats.

import SwiftyStatsPrelude

extension MLEFitter {

    /// Fit a Landau(location, scale) distribution via numerical MLE.
    ///
    /// Parameterization: location μ ∈ ℝ, scale c > 0.
    ///
    /// - Parameters:
    ///   - data: Sample values (any real numbers).
    ///   - optimizer: Optimizer used for the two-parameter solve. Defaults to `.lbfgs`
    ///     to exploit the analytic score.
    ///   - options: Optional solver configuration. When `nil`, covariance estimation
    ///     and diagnostics are enabled automatically.
    /// - Returns: `MLEResult` whose `thetaHat = [μ̂, ĉ]`.
    public static func fitLandau(
        _ data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)
        let dummy = MLEProblem<T>(data: data, logpdf: {_,_ in 0}, paramSpecs: [.init(.real, initial: 0)])
        let specs = dummy.makeParamSpecs(for: .landau, data: data)

        let logPDF: @Sendable (T, [T]) -> T = { x, th in
            if th[1] <= T.zero { return T.infinity }
            do { return try SwiftyBoost.Distribution.Landau(location: th[0], scale: th[1]).logPdf(x) }
            catch { return T.infinity }
        }
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, th in
            if th[1] <= T.zero { return [.nan, .nan] }
            do {
                let dist = try SwiftyBoost.Distribution.Landau(location: th[0], scale: th[1])
                let res = dist.score(x: x, mu: th[0], scale: th[1])
                return [res.dmu, res.dc]
            }
            catch { return [.nan, .nan] }
        }
        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)
        var opts = options ?? MLEOptimizationOpts<T>()
        opts.computeCovariance = true
        opts.optimizer = optimizer
        opts.diagnosticsEnabled = true
        return MLESolver.fit(problem, options: opts)
    }

}
