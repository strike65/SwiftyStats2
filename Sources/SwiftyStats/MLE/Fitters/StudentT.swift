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

/// Maximum likelihood fitter for the StudentT distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {
    /// Fit a Student’s t(ν) distribution using numerical MLE.
    ///
    /// Parameterization: degrees of freedom ν > 0.
    ///
    /// - Parameters:
    ///   - data: Sample values.
    ///   - optimizer: Optimizer to use (default: `.nelderMead`).
    ///   - options: Optional optimizer configuration.
    /// - Returns: `MLEResult` with `thetaHat = [ν̂]`.
    public static func fitStudentsT(
        _ data: [T],
        optimizer: OptimizerKind = .nelderMead,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)

        // Parameter specs for Student's t.
        let dummy = MLEProblem<T>(data: data, logpdf: {_,_ in 0}, paramSpecs: [.init(.real, initial: 0)])
        let specs = dummy.makeParamSpecs(for: .students_t, data: data)

        // Log-pdf with domain checks.
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let nu = theta[0]
            if !(nu.isFinite && nu > 0) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.StudentT(degreesOfFreedom: nu).logPdf(x)
            } catch {
                return T.infinity
            }
        }

        // Gradient (score) for ν.
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let nu = theta[0]
            guard nu.isFinite else { return [T.nan] }
            do {
                let dist: SwiftyBoost.Distribution.StudentT<T> = try SwiftyBoost.Distribution.StudentT(degreesOfFreedom: nu)
                let g = dist.score(x: x, nu: nu)
                return [g]
            } catch {
                return [T.nan]
            }
        }

        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)

        // Default optimizer setup.
        let opts: MLEOptimizationOpts<T> = {
            if let o = options {
                return o
            } else {
                var o = MLEOptimizationOpts<T>(); o.computeCovariance = true
                o.optimizer = optimizer
                return o
            }
        }()
        return MLESolver.fit(problem, options: opts)
    }
}
