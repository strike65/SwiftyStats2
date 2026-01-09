//
//  Created by VT on 19.11.25.
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

/// Score helpers for the InverseNormal distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the InverseNormal distribution used by MLEFitter.
extension SwiftyBoost.Distribution.InverseNormal {
    /// Score for Inverse-Gaussian (Wald) distribution X ~ IG(μ, λ), with μ > 0, λ > 0, x > 0.
    ///
    /// - Density: f(x) = sqrt(λ/(2π x³)) exp(−λ (x − μ)² / (2 μ² x)).
    /// - Log-likelihood (up to x-only terms): ℓ = 0.5 log λ − 0.5 log(2π) − 1.5 log x − λ (x − μ)² / (2 μ² x).
    /// - Score:
    ///   - ∂ℓ/∂μ = λ (x − μ)/μ³
    ///   - ∂ℓ/∂λ = 0.5/λ − (x − μ)²/(2 μ² x)
    public func score(x: T, mu: T, lambda: T) -> (dmu: T, dlambda: T) {
        precondition(x > 0 && mu > 0 && lambda > 0, "x>0, mu>0, lambda>0")
        let dmu: T = lambda * (x - mu) / (mu * mu * mu)// d log-likelihood / d mu
        let ex1 = ((x - mu) * (x - mu))
        let ex2 = (T.two * mu * mu * x)
        let dlambda: T  = T.half / lambda - ex1 / ex2 // d log-likelihood / d lambda
        return (dmu, dlambda)
    }

    /// Sum of Inverse-Gaussian scores across data.
    public func totalScore(data: [T], mu: T, lambda: T) -> (dmu: T, dlambda: T) {
        var sMu = T.zero, sLa = T.zero
        for x in data {
            let (a,b) = score(x: x, mu: mu, lambda: lambda)
            sMu += a; sLa += b
        }
        return (sMu, sLa)
    }

    /// Score with log-μ reparameterization (μ = exp(logMu)).
    ///
    /// - Returns: `(∂ℓ/∂logMu, ∂ℓ/∂lambda)` for one observation.
    @inline(__always)
    public func scoreWithLogMu(x: T, logMu: T, lambda: T) -> (dlogMu: T, dlambda: T) {
        let mu = T.exp(logMu)
        let (dmu, dl) = score(x: x, mu: mu, lambda: lambda)
        return (mu * dmu, dl)                        // chain rule: d/d log mu = mu * d/d mu
    }

    /// Sum of log-μ scores across all observations.
    @inline(__always)
    public func totalScoreWithLogMu(data: [T], logMu: T, lambda: T) -> (dlogMu: T, dlambda: T) {
        let mu = T.exp(logMu)
        var sLm = T.zero, sLa = T.zero
        for x in data {
            let (dmu, dl) = score(x: x, mu: mu, lambda: lambda)
            sLm += mu * dmu
            sLa += dl
        }
        return (sLm, sLa)
    }

    /// Score with log-λ reparameterization (λ = exp(logLambda)).
    ///
    /// - Returns: `(∂ℓ/∂mu, ∂ℓ/∂logLambda)` for one observation.
    @inline(__always)
    public func scoreWithLogLambda(x: T, mu: T, logLambda: T) -> (dmu: T, dlogLambda: T) {
        let lambda = T.exp(logLambda)
        let (dmu, dl) = score(x: x, mu: mu, lambda: lambda)
        return (dmu, lambda * dl)                    // d/d log lambda = lambda * d/d lambda
    }

    /// Sum of log-λ scores across all observations.
    @inline(__always)
    public func totalScoreWithLogLambda(data: [T], mu: T, logLambda: T) -> (dmu: T, dlogLambda: T) {
        let lambda = T.exp(logLambda)
        var sMu = T.zero, sLl = T.zero
        for x in data {
            let (dmu, dl) = score(x: x, mu: mu, lambda: lambda)
            sMu += dmu
            sLl += lambda * dl
        }
        return (sMu, sLl)
    }

    /// Score with log-μ and log-λ reparameterization.
    ///
    /// - Returns: `(∂ℓ/∂logMu, ∂ℓ/∂logLambda)` for one observation.
    @inline(__always)
    public func scoreWithLogParams(x: T, logMu: T, logLambda: T) -> (dlogMu: T, dlogLambda: T) {
        let mu = T.exp(logMu), lambda = T.exp(logLambda)
        let (dmu, dl) = score(x: x, mu: mu, lambda: lambda)
        return (mu * dmu, lambda * dl)
    }

    /// Sum of log-μ/log-λ scores across all observations.
    @inline(__always)
    public func totalScoreWithLogParams(data: [T], logMu: T, logLambda: T) -> (dlogMu: T, dlogLambda: T) {
        let mu = T.exp(logMu), lambda = T.exp(logLambda)
        var sLm = T.zero, sLl = T.zero
        for x in data {
            let (dmu, dl) = score(x: x, mu: mu, lambda: lambda)
            sLm += mu * dmu
            sLl += lambda * dl
        }
        return (sLm, sLl)
    }
}
