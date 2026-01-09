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

/// Maximum likelihood fitter for the Uniform distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {

    /// Fit a Uniform(a, b) distribution via order statistics.
    ///
    /// MLEs are the minimum and maximum of the sample.
    ///
    /// - Parameters:
    ///   - data: Sample values. Must be non-empty.
    ///   - options: Optional optimizer options (unused for analytic fits).
    /// - Returns: `MLEResult` with `thetaHat = [â, b̂]` where `â = min(data)` and `b̂ = max(data)`.
    /// - Note: If `â == b̂`, the log-likelihood is treated as `+∞` (degenerate interval).
    public static func fitUniform(
        data: [T],
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty, "data must not be empty")

        // Order statistics give the MLEs directly.
        let a: T = data.min()!
        let b: T = data.max()!
        let width = b - a
        let nT = T(data.count)

        // Log-likelihood for Uniform(a, b): -n log(b - a), assuming all points in [a, b].
        let logLik: T = (width == .zero) ? .infinity : -nT * T.log(width)
        return MLEResult(thetaHat: [a, b], logLik: logLik, iterations: 0, converged: true, nEval: 0, cov: nil)
    }
}
