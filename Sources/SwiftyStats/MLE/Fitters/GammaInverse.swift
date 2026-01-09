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

/// Maximum likelihood fitter for the GammaInverse distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {
    /// Fit an Inverse-Gamma(α, β) distribution using numerical MLE.
    ///
    /// Parameterization: shape α > 0, scale β > 0.
    ///
    /// - Parameters:
    ///   - data: Sample values (typically positive).
    ///   - optimizer: Optimizer to use (default: `.lbfgs`).
    ///   - options: Optional optimizer configuration (defaults emphasize robustness).
    /// - Returns: `MLEResult` with `thetaHat = [α̂, β̂]`.
    /// - Implementation notes:
    ///   - Uses multi-start and random restarts to mitigate local optima / poor conditioning.
    public static func fitInverseGamma(
        _ data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)

        // Parameter specs for Inverse-Gamma.
        let dummy = MLEProblem<T>(data: data, logpdf: { _, _ in 0 }, paramSpecs: [.init(.real, initial: 0)])
        let specs = dummy.makeParamSpecs(for: .inverse_gamma, data: data)

        // Log-pdf with domain checks.
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let alpha = theta[0]
            let beta  = theta[1]
            if !(alpha.isFinite && beta.isFinite && alpha > 0 && beta > 0) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.InverseGamma(shape: alpha, scale: beta).logPdf(x)
            } catch {
                return T.infinity
            }
        }

        // Gradient (score) for (α, β).
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let alpha = theta[0]
            let beta = theta[1]
            guard alpha > 0, beta > 0, alpha.isFinite, beta.isFinite else { return [T.nan, T.nan] }
            do {
                let dist: SwiftyBoost.Distribution.InverseGamma<T> = try SwiftyBoost.Distribution.InverseGamma(shape: alpha, scale: beta)
                let g = dist.score(x: x, alpha: alpha, beta: beta)
                return [g.dalpha, g.dbeta]
            } catch {
                return [T.nan, T.nan]
            }
        }

        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)

        // Robust optimizer defaults for potentially tricky likelihoods.
        // Use user options when available; otherwise rely on the solver defaults.
        return MLESolver.fit(problem, options: options ?? MLEOptimizationOpts<T>())
    }
}
