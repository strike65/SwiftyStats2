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

/// Maximum likelihood fitter for the Rayleigh distribution in SwiftyStats.

import SwiftyStatsPrelude

extension MLEFitter {

    /// Fit a Rayleigh distribution via its analytic MLE.
    ///
    /// Parameterization: variance-like parameter `v = σ²` (not the scale `σ`).
    /// The MLE is `v̂ = (1/(2n)) * Σ x_i²`.
    ///
    /// - Parameters:
    ///   - data: Sample values. Must be non-empty and nonnegative.
    ///   - optimizer: ``OptimzerKind`` setting (defaults to ``OptimizerKind.lbfgs``
    ///   - options: An optional ``NelderMeadOptions`` struct.
    /// - Returns: `MLEResult` with `thetaHat = [v̂]`.
    /// - Precondition: `!data.isEmpty && data.allSatisfy { $0 >= 0 }`.
    /// - Note: Degenerate samples with `Σ x_i² == 0` return `logLik = -∞` and `converged = false`.
    public static func fitRayleigh(
        data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        _ = optimizer
        _ = options
        precondition(!data.isEmpty, "data must not be empty")
        precondition(data.allSatisfy { $0 >= .zero }, "Rayleigh MLE requires nonnegative data")

        let n = T(data.count)
        let sumSquares = data.reduce(T.zero) { $0 + $1 * $1 }

        // Degenerate sample: no positive mass.
        guard sumSquares > .zero else {
            return MLEResult(thetaHat: [.zero], logLik: -.infinity, iterations: 0, converged: false, nEval: 0, cov: nil)
        }

        let sigmaSquared = sumSquares / (T.two * n)

        // Log-likelihood at the optimum: Σ log x - n log σ² - Σ x² / (2 σ²)
        var logLikelihood: T = .zero
        var hasFiniteLog = true
        for x in data {
            if x <= .zero {
                hasFiniteLog = false
                break
            }
            logLikelihood += T.log(x) - T.log(sigmaSquared) - (x * x) / (T.two * sigmaSquared)
        }
        if !hasFiniteLog {
            logLikelihood = -.infinity
        }

        return MLEResult(
            thetaHat: [T.sqrt(sigmaSquared)],
            logLik: logLikelihood,
            iterations: 0,
            converged: true,
            nEval: 0,
            cov: nil
        )
    }
}
