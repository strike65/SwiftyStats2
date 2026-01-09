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

/// Maximum likelihood fitter for the ChiSquaredInverse distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {
    /// Fit a Scaled Inverse Chi-Squared(ν, τ²) distribution using numerical MLE.
    ///
    /// Parameterization: degrees of freedom ν > 0, scale τ² > 0.
    ///
    /// - Parameters:
    ///   - data: Sample values (strictly positive).
    ///   - optimizer: Optimizer to use (default: `.lbfgs` to leverage the analytic score).
    ///   - options: Optional optimizer configuration (sensible robust defaults when `nil`).
    /// - Returns: `MLEResult` with `thetaHat = [ν̂, τ̂²]`.
    /// - Implementation notes:
    ///   - Uses factory starts from `.inverse_chi_squared`.
    ///   - Uses the analytic score from Score-functions.swift: `score(x:nu:s2:) -> (dNu, dS2)`.
    public static func fitInverseChiSquared(
        _ data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty, "data must not be empty")

        // Parameter specs for Scaled Inverse Chi-Squared.
        let dummy = MLEProblem<T>(data: data, logpdf: { _, _ in 0 }, paramSpecs: [.init(.positive, initial: 1)])
        let specs = dummy.makeParamSpecs(for: .inverse_chi_squared, data: data)

        // Log-pdf with domain checks (ν > 0, τ² > 0).
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let nu = theta[0]
            let s2 = theta[1]
            if !(nu.isFinite && s2.isFinite && nu > 0 && s2 > 0) { return T.infinity }
            do {
                // Assumed API: degreesOfFreedom, scale
                return try SwiftyBoost.Distribution.InverseChiSquared(degreesOfFreedom: nu, scale: s2).logPdf(x)
            } catch {
                return T.infinity
            }
        }

        // Analytic gradient using the implemented score(x:nu:s2:)
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let nu = theta[0]
            let s2 = theta[1]
            guard nu.isFinite, s2.isFinite, nu > 0, s2 > 0 else { return [T.nan, T.nan] }
            do {
                let dist: SwiftyBoost.Distribution.InverseChiSquared<T> = try SwiftyBoost.Distribution.InverseChiSquared(degreesOfFreedom: nu, scale: s2)
                let g = dist.score(x: x, nu: nu, s2: s2) // (dNu, dS2)
                return [g.dNu, g.dS2]
            } catch {
                return [T.nan, T.nan]
            }
        }

        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)

        // Robust optimizer defaults (multi-start, covariance, diagnostics).
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
                o.initialStepStrategy = .relativeTheta(T(0.35))
                o.gradStep = T(5e-6)
                o.hessianStep = T(1e-3)
                o.relTolLogLik = T(1e-7)
                o.diagnosticsEnabled = true
                return o
            }
        }()
        return MLESolver<T>.fit(problem, options: opts)
    }
}
