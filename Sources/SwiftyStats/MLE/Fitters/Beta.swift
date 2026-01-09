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

/// Maximum likelihood fitter for the Beta distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {

    /// Fit a Beta(α, β) distribution using numerical MLE.
    ///
    /// - Parameters:
    ///   - data: Sample values (typically in (0, 1); the distribution’s logPdf governs validity).
    ///   - optimizer: Optimizer to use (default: `.nelderMead`).
    ///   - options: Optional optimizer configuration. If `nil`, a sensible default is used and covariance is computed.
    /// - Returns: `MLEResult` with `thetaHat = [α̂, β̂]`.
    /// - Implementation notes:
    ///   - Builds an `MLEProblem` with the Beta log-pdf from `SwiftyBoost.Distribution`.
    ///   - Provides an analytic score function for gradient-based optimizers.
    public static func fitBeta(
        _ data: [T],
        optimizer: OptimizerKind = .nelderMead,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        // Use a "dummy" problem to obtain consistent parameter specs (bounds, initial values).
        let dummy = MLEProblem<T>(data: data, logpdf: {_,_ in 0}, paramSpecs: [.init(.positive, initial: 0)])
        let specs = dummy.makeParamSpecs(for: .beta, data: data)

        // Log-pdf closure with parameter domain checks (return +∞ for infeasible).
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let alpha = theta[0]
            let beta = theta[1]
            if !(alpha.isFinite && beta.isFinite && alpha > 0 && beta > 0) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.Beta(alpha: alpha, beta: beta).logPdf(x)
            }
            catch _ {
                return T.infinity
            }
        }

        // Analytic gradient (score) to accelerate convergence.
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let a = theta[0]
            let b = theta[1]
            guard a.isFinite, b.isFinite, a > 0, b > 0 else { return [T.nan, T.nan] }
            do {
                let dist: SwiftyBoost.Distribution.Beta<T> = try SwiftyBoost.Distribution.Beta(alpha: a, beta: b)
                let g = dist.score(x: x, a: a, b: b)
                return [g.da, g.db]
            } catch {
                return [T.nan, T.nan]
            }
        }

        // Assemble the MLE problem.
        let problem = MLEProblem(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)

        // Default optimizer configuration: compute covariance and use requested optimizer.
        let opts: MLEOptimizationOpts<T> = {
            if let o = options {
                return o
            } else {
                var o = MLEOptimizationOpts<T>()
                o.computeCovariance = true
                o.optimizer = optimizer
                return o
            }
        }()
        return MLESolver<T>.fit(problem, options: opts)
    }
}
