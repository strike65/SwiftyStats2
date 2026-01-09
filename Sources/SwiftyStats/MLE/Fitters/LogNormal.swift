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
/// Maximum likelihood fitter for the LogNormal distribution in SwiftyStats.

import SwiftyStatsPrelude

extension MLEFitter {
    /// Fit a LogNormal(meanlog = mu, sdlog = sigma) via numerical MLE with analytic gradient.
    ///
    /// - Parameters:
    ///   - data: Sample values (support x > 0).
    ///   - optimizer: Optimizer to use (default: `.lbfgs`).
    ///   - options: Optional optimizer configuration.
    /// - Returns: `MLEResult` with `thetaHat = [mu_hat, sigma_hat]`.
    public static func fitLogNormal(
        _ data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty, "data must be non-empty")
        
        // Filter valid support (x > 0). If your pipeline guarantees this already, you can skip.
        let xPos = data.filter { $0 > .zero && $0.isFinite }
        precondition(!xPos.isEmpty, "lognormal MLE requires x > 0")
        // Parameter specs for Logistic.
        let dummy = MLEProblem<T>(data: data, logpdf: {_,_ in 0}, paramSpecs: [.init(.real, initial: 0)])
        let specs = dummy.makeParamSpecs(for: .lognormal, data: data)
        
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let mu = theta[0]
            let s  = theta[1]
            if !(mu.isFinite && s.isFinite && s > 0) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.LogNormal(location: mu, scale: s).logPdf(x)
            } catch {
                return T.infinity
            }
        }
        
        // Analytic score for (μ, s).
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let mu = theta[0]
            let s  = theta[1]
            guard mu.isFinite, s.isFinite, s > 0 else { return [T.nan, T.nan] }
            do {
                let dist: SwiftyBoost.Distribution.LogNormal<T> = try SwiftyBoost.Distribution.LogNormal(location: mu, scale: s)
                let g = dist.score(x: x, mu: mu, sigma: s)
                return [g.dmu, g.dsigma]
            } catch {
                return [T.nan, T.nan]
            }
        }
        
        let problem = MLEProblem<T>(
            data: xPos,
            logpdf: logPDF,
            gradlogpdf: gradLogPDF,
            paramSpecs: specs
        )
        
        // Reasonable defaults; adjust to your global conventions if needed
        let opts: MLEOptimizationOpts<T> = {
            if let o = options { return o }
            var o = MLEOptimizationOpts<T>()
            o.computeCovariance = true
            o.optimizer = optimizer          // .lbfgs recommended
            o.multiStartCount = 4
            o.multiStartDesign = .lhs
            o.randomRestartCount = 1
            o.initialStepStrategy = .relativeTheta(T(0.25))
            o.gradStep = T(1e-5)
            o.hessianStep = T(5e-4)
            o.relTolLogLik = T(1e-7)
            o.diagnosticsEnabled = true
            return o
        }()
        
        return MLESolver.fit(problem, options: opts)
    }
}
