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

/// Maximum likelihood fitter for the Arcsine distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {

    /// Fit an Arcsine distribution using its analytic MLE for the support [a, b].
    ///
    /// The MLE of the support endpoints is given by the sample min and max.
    ///
    /// - Parameters:
    ///   - data: Sample values. Must be non-empty.
    ///   - options: Optional optimizer options (unused for analytic fits).
    /// - Returns: `MLEResult` with `thetaHat = [â, b̂]` where `â = min(data)` and `b̂ = max(data)`.
    /// - Note:
    ///   - If `â == b̂` or data include exact boundary hits, the log-likelihood is treated as `+∞`.
    ///   - The log-likelihood is Σ[-log(π) - 1/2 log(x-a) - 1/2 log(b-x)].
    public static func fitArcsine(
        data: [T],
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty, "data must not be empty")

        // Endpoints from order statistics.
        let a = data.min()!
        let b = data.max()!
        let width = b - a

        // Degenerate or boundary-hitting cases are treated as +∞ (unbounded likelihood).
        if width == .zero {
            return MLEResult(thetaHat: [a, b], logLik: .infinity, iterations: 0, converged: true, nEval: 0, cov: nil)
        }
        if data.contains(where: { $0 == a || $0 == b }) {
            return MLEResult(thetaHat: [a, b], logLik: .infinity, iterations: 0, converged: true, nEval: 0, cov: nil)
        }

        // Sum log-likelihood terms: -log(π) - 0.5 log(x-a) - 0.5 log(b-x)
        var logLik: T = 0
        let c0 = -T.log(T.pi)
        for x in data {
            let t1 = x - a
            let t2 = b - x
            logLik += c0 - T.half * T.log(t1) - T.half * T.log(t2)
        }
        return MLEResult(thetaHat: [a, b], logLik: logLik, iterations: 0, converged: true, nEval: 0, cov: nil)
    }
}
