//
//  Created by VT on 09.11.25.
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

/// Maximum likelihood fitter for the Bernoulli distribution in SwiftyStats.

import SwiftyStatsPrelude

extension MLEFitter {
    /// Fit the success probability of a Bernoulli(p) distribution.
    ///
    /// - Parameters:
    ///   - data: Array of 0/1 outcomes.
    ///   - optimizer: Local optimizer to use for the single-parameter fit.
    ///   - options: Optional Nelder–Mead configuration (covariance enabled automatically).
    /// - Returns: `MLEResult` containing the estimated probability `p`.
    public static func fitBernoulli(
        data: [T],
        optimizer: OptimizerKind = .nelderMead,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        _ = optimizer
        _ = options
        precondition(!data.isEmpty, "data must not be empty")
        precondition(
            data.allSatisfy { $0 == .zero || $0 == T.one },
            "Bernoulli data must consist of 0/1 outcomes"
        )

        let successes = data.reduce(T.zero, +)
        let count = T(data.count)
        let pHat = successes / count

        var logLik: T = .zero
        if pHat <= .zero || pHat >= T.one {
            logLik = .infinity
        } else {
            for x in data {
                if x == .one {
                    logLik += T.log(pHat)
                } else {
                    logLik += T.log(T.one - pHat)
                }
            }
        }

        return MLEResult(
            thetaHat: [pHat],
            logLik: logLik,
            iterations: 0,
            converged: true,
            nEval: 0,
            cov: nil
        )
    }
}
