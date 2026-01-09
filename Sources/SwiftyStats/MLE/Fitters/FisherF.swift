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

/// Maximum likelihood fitter for the FisherF distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {

    private static func fisherFWarmStart(data: [T]) -> [T]? {
        let filtered = data.filter { $0.isFinite && $0 > .zero }
        guard filtered.count >= 4 else { return nil }

        let n = T(filtered.count)
        let mean = filtered.reduce(T.zero, +) / n
        guard mean > T.one + T(1e-3) else { return nil }

        let variance = filtered.reduce(T.zero) { acc, value in
            let delta = value - mean
            return acc + delta * delta
        } / n
        guard variance > T(1e-8) else { return nil }

        let d2Candidate = (T.two * mean) / (mean - T.one)
        guard d2Candidate.isFinite && d2Candidate > T(4.5) else { return nil }

        let numerator = T.two * d2Candidate * d2Candidate * (d2Candidate - T.two)
        let denomHelper = (d2Candidate - T.two) * (d2Candidate - T.two) * (d2Candidate - T(4))
        let denominator = variance * denomHelper - T.two * d2Candidate * d2Candidate
        guard denominator > T(1e-6) else { return nil }

        let d1Candidate = numerator / denominator
        guard d1Candidate.isFinite && d1Candidate > T(0.5) else { return nil }

        return [d1Candidate, d2Candidate]
    }

    /// Fit a Fisher F(d1, d2) distribution using numerical MLE.
    ///
    /// Parameterization: numerator df d1 > 0, denominator df d2 > 0.
    ///
    /// - Parameters:
    ///   - data: Sample values.
    ///   - optimizer: Optimizer to use (default: `.nelderMead`).
    ///   - options: Optional optimizer configuration.
    /// - Returns: `MLEResult` with `thetaHat = [d1̂, d2̂]`.
    public static func fitFisherF(
        _ data: [T],
        optimizer: OptimizerKind = .nelderMead,
        options: MLEOptimizationOpts<T>? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty)

        // Parameter specs for F distribution.
        let dummy = MLEProblem<T>(data: data, logpdf: {_,_ in 0}, paramSpecs: [.init(.real, initial: 0)])
        let specs = dummy.makeParamSpecs(for: .fisher_f, data: data)

        // Log-pdf with domain checks.
        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            let d1 = theta[0]
            let d2 = theta[1]
            if !(d1.isFinite && d2.isFinite && d1 > 0 && d2 > 0) { return T.infinity }
            do {
                return try SwiftyBoost.Distribution.FisherF(degreesOfFreedom1: d1, degreesOfFreedom2: d2).logPdf(x)
            } catch {
                return T.infinity
            }
        }

        // Gradient (score) for (d1, d2).
        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            let df1 = theta[0]
            let df2 = theta[1]
            guard df1 > 0, df2 > 0, df1.isFinite, df2.isFinite else { return [T.nan, T.nan] }
            do {
                let dist: SwiftyBoost.Distribution.FisherF<T> = try SwiftyBoost.Distribution.FisherF(degreesOfFreedom1: df1, degreesOfFreedom2: df2)
                let g = dist.score(x: x, d1: df1, d2: df2)
                return [g.dd1, g.dd2]
            } catch {
                return [T.nan, T.nan]
            }
        }

        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)

        var resolvedOptions: MLEOptimizationOpts<T>
        if let custom = options {
            resolvedOptions = custom
        } else {
            var o = MLEOptimizationOpts<T>()
            o.optimizer = optimizer
            o.computeCovariance = true
            o.multiStartCount = 8
            o.randomRestartCount = 3
            o.multiStartDesign = .sobol
            o.initialStepStrategy = .relativeTheta(T(0.3))
            o.gradStep = T(5e-6)
            o.hessianStep = T(1e-3)
            o.relTolLogLik = T(1e-7)
            o.diagnosticsEnabled = true
            resolvedOptions = o
        }

        if resolvedOptions.warmStartTheta == nil,
           let warmStart = fisherFWarmStart(data: data) {
            resolvedOptions.warmStartTheta = warmStart
        }

        return MLESolver.fit(problem, options: resolvedOptions)
    }
}
