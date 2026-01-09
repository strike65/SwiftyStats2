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
/// Maximum likelihood fitter for the Pareto distribution in SwiftyStats.

import SwiftyStatsPrelude

extension MLEFitter {
    /// Fits a Pareto-I distribution with the scale fixed to the sample minimum (MLE property) and estimates the shape α.
    ///
    /// Parameterization and support
    /// - Pareto-I with scale `x_m > 0` and shape `α > 0`.
    /// - Support: `x ≥ x_m`.
    ///
    /// Density
    /// - f(x | x_m, α) = α x_m^α / x^{α+1}, for `x ≥ x_m`; 0 otherwise.
    ///
    /// Estimation strategy
    /// - For Pareto-I, the MLE of the scale is `x̂_m = min(data)`. This routine fixes `x_m` to the
    ///   sample minimum and maximizes the likelihood over `α > 0` only.
    /// - Given observations `x₁,…,x_n ≥ x_m`, the log-likelihood in α is
    ///   ℓ(α) = n · log α + n · α · log x_m − (α + 1) · Σ log x_i.
    /// - Closed-form MLE (with fixed `x_m`): `α̂ = n / Σ log(x_i / x_m)`. The numerical optimizer
    ///   should recover this under regular conditions.
    ///
    /// Implementation details
    /// - Data are filtered to positive finite values; `x_m` is set to their minimum.
    /// - Parameter constraints (`α > 0`) are enforced via the solver’s positive transform using
    ///   `makeParamSpecs(for: .pareto, data:)`.
    /// - Per-observation log-pdf and analytic score are delegated to
    ///   `SwiftyBoost.Distribution.Pareto(scale: x_m, shape: α)`.
    /// - This wrapper enables covariance (θ-space) and diagnostics by default.
    ///
    /// - Parameters:
    ///   - data: Sample values; must be non-empty. Only positive, finite observations are considered.
    ///           The scale is fixed to their minimum.
    ///   - optimizer: Local optimizer (`.nelderMead`, `.bfgs`, `.lbfgs`). Default is `.lbfgs`.
    ///   - options: Optional solver options. If `nil`, sensible defaults are used and adjusted to compute covariance and enable diagnostics.
    ///
    /// - Returns: An `MLEResult<T>` with:
    ///   - `thetaHat = [α̂]` (shape estimate; the fitted scale is `x_m = min(filtered data)`),
    ///   - `logLik` (maximized log-likelihood),
    ///   - convergence flags and evaluation counts,
    ///   - optional covariance and diagnostics when enabled.
    ///
    /// - Precondition: `data` must be non-empty and contain at least one positive finite value.
    ///
    /// - Note: This routine fits Pareto-I with fixed scale. Estimating `x_m` jointly with `α` is not well-posed
    ///   in pure MLE (it collapses to the sample minimum); use a reparameterization and constraints if you need a free scale.
    public static func fitPareto(
        _ data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)
        let x = data.filter { $0.isFinite && $0 > 0 }
        precondition(!x.isEmpty)
        let xm = x.min()!
        // Parameter specs for Pareto-I (estimate α > 0; x_m fixed to xm).
        let dummy = MLEProblem<T>(data: data, logpdf: { _, _ in 0 }, paramSpecs: [.init(.real, initial: 0)])
        let specs = dummy.makeParamSpecs(for: .pareto, data: data)

        let logPDF: @Sendable (T, [T]) -> T = { v, theta in
            do {
                let alpha = theta[0]
                if !(v >= xm && alpha > .zero) { return T.infinity }
                let dist = try SwiftyBoost.Distribution.Pareto<T>(scale: xm, shape: alpha)
                let res = try dist.logPdf(v)
                return res
            }
            catch _ {
                return T.infinity
            }
        }

        // Analytic score with respect to α (scale fixed).
        let gradLogPdf: @Sendable (T, [T]) -> [T] = { v, theta in
            do {
                let alpha = theta[0]
                if !(v >= xm && alpha > .zero) { return [T.nan, .nan] }
                let dist = try SwiftyBoost.Distribution.Pareto<T>(scale: xm, shape: alpha)
                let res = dist.score(x: v, scale: xm, shape: alpha)
                return [res.dalpha]
            }
            catch _ {
                return [T.nan, .nan]
            }
        }

        let problem = MLEProblem<T>(data: x, logpdf: logPDF, gradlogpdf: gradLogPdf, paramSpecs: specs)
        var opts = options ?? MLEOptimizationOpts<T>()
        opts.computeCovariance = true
        opts.optimizer = optimizer
        opts.diagnosticsEnabled = true
        return MLESolver.fit(problem, options: opts)
    }
}
