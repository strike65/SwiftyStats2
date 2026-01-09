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
/// Maximum likelihood fitter for the NormalInverse distribution in SwiftyStats.

import SwiftyStatsPrelude

extension MLEFitter {

    /// Fit an Inverse Normal (Wald) distribution via numerical MLE.
    ///
    /// Parameterization: mean μ > 0, shape λ > 0.
    ///
    /// - Parameters:
    ///   - data: Sample values (strictly positive).
    ///   - optimizer: Optimizer powering the solver (defaults to `.lbfgs` to use the
    ///     analytic score).
    ///   - options: Optional solver settings. When `nil`, the solver falls back to the
    ///     default Nelder–Mead options.
    /// - Returns: `MLEResult` with `thetaHat = [μ̂, λ̂]`.
    public static func fitInverseNormal(
        _ data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)
        
        // Parameter specs for Inverse Normal (Wald): μ > 0, λ > 0.
        let dummy = MLEProblem<T>(data: data, logpdf: { _, _ in 0 }, paramSpecs: [.init(.real, initial: 0)])
        let specs = dummy.makeParamSpecs(for: .inverseNormal, data: data)
        
        // Log-pdf with domain checks.
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let mu = theta[0]
            let lambda  = theta[1]
            if !(mu.isFinite && lambda.isFinite && mu > 0 && lambda > 0) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.InverseNormal(mean: mu, shape: lambda).logPdf(x)
            } catch {
                return T.infinity
            }
        }
        
        // Analytic score ∂ℓ/∂(μ, λ) per observation.
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let mu = theta[0]
            let lambda = theta[1]
            guard mu > 0, lambda > 0, mu.isFinite, lambda.isFinite else { return [T.nan, T.nan] }
            do {
                let dist: SwiftyBoost.Distribution.InverseNormal<T> = try SwiftyBoost.Distribution.InverseNormal(mean: mu, shape: lambda)
                let g = dist.score(x: x, mu: mu, lambda: lambda)
                return [g.dmu, g.dlambda]
            } catch {
                return [T.nan, T.nan]
            }
        }
        
        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)
        
        // Use provided options or the solver defaults.
        return MLESolver.fit(problem, options: options ?? MLEOptimizationOpts<T>())
    }
}
