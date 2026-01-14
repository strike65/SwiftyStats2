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

/// Maximum likelihood fitter for the Binomial distribution in SwiftyStats.

import SwiftyStatsPrelude

extension MLEFitter {

    /// Fits a Binomial(n, p) distribution to count data using maximum likelihood.
    ///
    /// This routine estimates the success probability `p` given a known number of trials `n`
    /// for each observation in `data`. Each element of `data` is treated as a count of
    /// successes in `n` Bernoulli trials.
    ///
    /// - Parameters:
    ///   - data: The observed counts of successes. Each `x` in `data` must satisfy `0 ≤ x ≤ n`.
    ///   - nInt: The fixed number of trials per observation (non-negative integer).
    ///   - optimizer: Optimizer to use. Defaults to `.nelderMead`. See `MLEOptimizationOpts`.
    ///   - options: Optional optimizer configuration. When provided, `computeCovariance` is
    ///              forced to `true` so that standard errors can be derived from the Hessian.
    ///
    /// - Returns: An `MLEResult<T>` with:
    ///   - `thetaHat = [p̂]` (the MLE of the success probability),
    ///   - `logLik` (maximized log-likelihood),
    ///   - solver diagnostics (iterations, convergence flags, optional covariance).
    ///
    /// - Preconditions:
    ///   - `n >= 0`
    ///   - Every `x` in `data` satisfies `0 ≤ x ≤ n`.
    ///
    /// - Discussion:
    ///   The parameterization is the canonical Binomial with a single free parameter `p ∈ (0, 1)`.
    ///   The function builds a one-parameter MLE problem with domain checks and an analytic score.
    ///
    /// - Important:
    ///   Calls must use the labeled parameter `n:`:
    ///   ```
    ///   let fit = MLEFitter<Double>.fitBinomial(myCounts, n: 12)
    ///   ```
    ///   Omitting the label (e.g., `fitBinomial(myCounts, 12)`) will result in
    ///   “Missing argument label 'n:' in call”.
    ///
    /// - Example:
    ///   ```
    ///   let counts: [Double] = [5, 7, 6, 8, 4, 9]   // successes out of n=10 trials
    ///   let result = MLEFitter<Double>.fitBinomial(counts, n: 10)
    ///   let pHat = result.thetaHat[0]               // estimated success probability
    ///   ```
    public static func fitBinomial(
        _ data: [T],
        n nInt: Int,
        optimizer: OptimizerKind = .nelderMead,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(nInt >= 0, "n must be >= 0")
        let n = T(nInt)
        
        let dummy = MLEProblem<T>(data: data, logpdf: {_,_ in 0}, paramSpecs: [.init(.unitInterval, initial: T.half)])
        // Provide context with known n so the factory can seed p well.
        var ctx = MLEStartContext<T>()
        ctx.binomialN = nInt
        let specs = dummy.makeParamSpecs(for: .binomial, data: data, ctx: ctx)
        
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            // p ∈ (0,1)
            let p = theta[0]
            guard p.isFinite, p > 0, p < 1, n.isFinite, n >= 0 else { return T.infinity }
            guard x >= T.zero, x <= n else { return T.infinity }
            do {
                let dist = try SwiftyBoost.Distribution.Binomial<T>(numberOfTrials: n, probabibilityOfSuccess: p)
                return try dist.logPdf(x)
            } catch {
                return T.infinity
            }
        }
        
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let p = theta[0]
            guard p.isFinite, p > 0, p < 1, n.isFinite, n >= 0 else { return [T.nan] }
            do {
                let dist = try SwiftyBoost.Distribution.Binomial<T>(numberOfTrials: n, probabibilityOfSuccess: p)
                let g = dist.score(x: x, n: n, p: p) // (dp)
                return [g]
            } catch {
                return [T.nan]
            }
        }
        
        let problem = MLEProblem(
            data: data,
            logpdf: logPDF,
            gradlogpdf: gradLogPDF,
            paramSpecs: specs
        )
        
        let opts: MLEOptimizationOpts<T> = {
            if var o = options {
                o.computeCovariance = true
                o.optimizer = optimizer
                return o
            } else {
                var o = MLEOptimizationOpts<T>()
                o.computeCovariance = true
                o.optimizer = optimizer
                return o
            }
        }()
        
        return MLESolver<T>.fit(problem, options: opts)
    }
}
