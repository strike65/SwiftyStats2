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

/// Maximum likelihood fitter for the Exponential distribution in SwiftyStats.

import SwiftyStatsPrelude

extension MLEFitter {
    /// Fit an Exponential(rate λ) distribution using its analytic MLE.
    ///
    /// The MLE for the rate parameter is λ̂ = n / Σ x_i.
    ///
    /// - Parameters:
    ///   - data: Sample values. Must be non-empty and nonnegative.
    ///   - options: Optional optimizer options (unused for analytic fits; accepted for API symmetry).
    /// - Returns: `MLEResult` with `thetaHat = [λ̂]` where `λ̂ = n / sum(data)`. If `sum(data) == 0`, returns `λ̂ = +∞` and `logLik = +∞`.
    /// - Precondition: `!data.isEmpty && data.allSatisfy { $0 >= 0 }`.
    /// - Note: This is a closed-form solution; no optimizer iterations are performed.
    public static func fitExponential(
        data: [T],
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty, "data must not be empty")
        precondition(data.allSatisfy { $0 >= .zero }, "Exponential MLE requires nonnegative data")

        // Sum and count in native T for numerical stability/consistency.
        var d = data
        let s = Helpers.sum(&d)
        let nT = T(data.count)

        // Degenerate case: all zeros ⇒ λ̂ = +∞ and log-likelihood = +∞ under the usual convention.
        if s == .zero {
            return MLEResult(thetaHat: [.infinity],
                             logLik: .infinity,
                             iterations: 0, converged: true, nEval: 0, cov: nil)
        }

        // Closed-form MLE and corresponding log-likelihood.
        let lambdaHat: T = nT / s
        let logLik = nT * T.log(lambdaHat) - lambdaHat * s
        return MLEResult(thetaHat: [lambdaHat], logLik: logLik, iterations: 0, converged: true, nEval: 0, cov: nil)
    }
}
