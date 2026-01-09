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

/// Maximum likelihood fitter for the Wald distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {

    /// Fit a Wald (Inverse Gaussian) distribution via analytic formulas.
    ///
    /// Parameterization: Wald(μ, λ) where μ is the mean and λ is the shape.
    ///
    /// - Parameters:
    ///   - data: Sample values. Must be strictly positive.
    ///   - options: Optional optimizer options (unused for analytic fits).
    /// - Returns: `MLEResult` with `thetaHat = [μ̂, λ̂]` where `μ̂ = mean(data)` and `λ̂` from the standard closed-form MLE.
    /// - Precondition: `data.allSatisfy { $0 > 0 } && mean(data) > 0`.
    public static func fitWald(
        data: [T],
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty, "data must not be empty")
        precondition(data.allSatisfy { $0 > .zero }, "Wald (Inverse Gaussian) MLE requires strictly positive data")

        // Compute mean and required sums for the Wald closed-form MLE.
        let ex: SSExamine<T,T> = try! .init(using: data, levelOfMeasurement: .ratio, name: nil, characterSet: nil)
        let m = ex.arithmeticMean!
        precondition(m > .zero, "Wald (Inverse Gaussian) MLE requires positive mean")

        let nT = T(ex.sampleSize)
        let m2 = m * m

        // s = Σ (x_i - μ)² / (μ² x_i)
        var s: T = 0
        var sumLogX: T = 0
        for x in data {
            let d = x - m
            s += (d * d) / (m2 * x)
            sumLogX += T.log(x)
        }

        // λ̂ = n / s (with the usual ∞ convention if s == 0).
        let lambdaHat: T = (s == .zero) ? .infinity : (nT / s)

        // Log-likelihood at the MLE (up to constants) for reporting.
        let logLik: T
        if s == .zero {
            logLik = .infinity
        } else {
            let nOver2 = nT / T.two
            let constTerm = nOver2 * (T.log(lambdaHat) - T.log(T.twoPi))
            let logXTerm = (T(3) / T.two) * sumLogX
            let quadTerm = (lambdaHat / T.two) * s
            logLik = constTerm - logXTerm - quadTerm
        }
        return MLEResult(thetaHat: [m, lambdaHat], logLik: logLik, iterations: 0, converged: true, nEval: 0, cov: nil)
    }
}
