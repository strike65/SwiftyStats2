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

/// Maximum likelihood fitter for the Holtsmark distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {
    /// Fit a Holtsmark(μ, c) distribution using numerical MLE with robust warm-starts.
    ///
    /// Parameterization: location μ ∈ ℝ, scale c > 0.
    ///
    /// - Parameters:
    ///   - data: Sample values.
    ///   - optimizer: Optimizer to use (default: `.lbfgs`).
    ///   - options: Optional optimizer configuration. If `warmStartTheta` is not set, a robust start from median and MAD/IQR is used.
    /// - Returns: `MLEResult` with `thetaHat = [μ̂, ĉ]`.
    /// - Implementation notes:
    ///   - Uses robust descriptive statistics for warm starts (median, MAD, IQR).
    ///   - Provides analytic score for gradient-based optimization.
    public static func fitHoltsmark(
        _ data: [T],
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)

        // Parameter specs for Holtsmark.
        let dummy = MLEProblem<T>(data: data, logpdf: { _, _ in 0 }, paramSpecs: [.init(.positive, initial: 1)])
        let specs = dummy.makeParamSpecs(for: .holtsmark, data: data)

        // Robust warm start: median for location, MAD/IQR for scale.
        let ex: SSExamine<T, T> = try! .init(using: data, levelOfMeasurement: .ratio, name: nil, characterSet: nil)
        let muWarm = ex.median ?? (ex.arithmeticMean ?? T.zero)
        let mad = ex.medianAbsoluteDeviation(center: muWarm) ?? T.zero
        let quart = try? ex.quartile()
        let iqr = (quart?.upper ?? muWarm) - (quart?.lower ?? muWarm)
        let cWarm = max(mad, iqr * T(0.5), T(1e-6))

        // Log-pdf with domain checks.
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let mu0  = theta[0]
            let c0   = theta[1]
            if !(mu0.isFinite && c0.isFinite && c0 > 0) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.Holtsmark(loc: mu0, scale: c0).logPdf(x)
            } catch {
                return T.infinity
            }
        }

        // Gradient (score) for (μ, c).
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let mu    = theta[0]
            let sigma = theta[1]
            guard mu.isFinite, sigma.isFinite, sigma > 0 else {
                return [T.nan, T.nan]
            }
            do {
                let hm: SwiftyBoost.Distribution.Holtsmark<T> = try SwiftyBoost.Distribution.Holtsmark<T>(loc: mu, scale: sigma)
                let g = hm.score(x: x, mu: mu, sigma: sigma)
                return [g.dmu, g.dsigma]
            }
            catch _ {
                return [T.nan, T.nan]
            }
        }

        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)

        // Optimizer defaults and warm-start injection.
        let opts: MLEOptimizationOpts<T> = {
            if var o = options {
                o.optimizer = optimizer
                o.computeCovariance = o.computeCovariance
                if o.warmStartTheta == nil {
                    o.warmStartTheta = [muWarm, cWarm]
                }
                if o.gradStep < T(1e-5) { o.gradStep = T(1e-5) }
                o.diagnosticsEnabled = true
                return o
            } else {
                var o = MLEOptimizationOpts<T>(); o.computeCovariance = true
                o.optimizer = optimizer
                o.multiStartCount = 12
                o.multiStartDesign = .lhs
                o.randomRestartCount = 3
                o.initialStepStrategy = .relativeTheta(T(0.35))
                o.gradStep = T(1e-5)
                o.hessianStep = T(1e-3)
                o.relTolLogLik = T(1e-7)
                o.diagnosticsEnabled = true
                o.warmStartTheta = [muWarm, cWarm]
                return o
            }
        }()
        return MLESolver.fit(problem, options: opts)
    }
}
