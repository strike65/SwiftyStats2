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

/// Maximum likelihood fitter for the Laplace distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {
    /// Fit a Laplace (double exponential) distribution using robust statistics (median and MAD).
    ///
    /// Parameterization: Laplace(m, b), where m is the location (median)
    /// and b is the scale (mean absolute deviation about the median).
    ///
    /// - Parameters:
    ///   - data: Sample values. Must be non-empty.
    ///   - options: Optional optimizer options (unused for analytic fits).
    /// - Returns: `MLEResult` with `thetaHat = [m̂, b̂]` where `m̂ = median(data)` and `b̂ = mean absolute deviation` around the median.
    /// - Note:
    ///   - This is the standard closed-form MLE for Laplace.
    ///   - If `b̂ == 0`, the log-likelihood is treated as `+∞` (degenerate).
    public static func fitLaplace(
        data: [T],
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty, "data must not be empty")

        // Use robust descriptive stats for the closed-form MLE.
        let ex: SSExamine<T,T> = try! .init(using: data, levelOfMeasurement: .ratio, name: nil, characterSet: nil)
        let m = ex.median!
        let bHat: T = ex.meanAbsoluteDeviation(center: m)!
        let nT = T(ex.sampleSize)

        // Laplace log-likelihood at the MLE: -n log(2b) - n.
        let logLik: T
        if bHat == .zero {
            logLik = .infinity
        } else {
            logLik = -nT * T.log(T.two * bHat) - nT
        }
        return MLEResult(thetaHat: [m, bHat], logLik: logLik, iterations: 0, converged: true, nEval: 0, cov: nil)
    }
}
