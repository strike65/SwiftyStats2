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

/// Score helpers for the Normal distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Normal distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Normal {
    /// Score for Normal(x | μ, σ), σ > 0.
    ///
    /// - Log-likelihood: ℓ = −log(σ√(2π)) − (x − μ)²/(2σ²).
    /// - Score:
    ///   - ∂ℓ/∂μ = (x − μ)/σ²
    ///   - ∂ℓ/∂σ = −1/σ + (x − μ)²/σ³
    public func score(x: T, mu: T, sigma: T) -> (dmu: T, dsigma: T) {
        precondition(sigma > 0, "sigma must be > 0")
        let z = x - mu
        let s2 = sigma * sigma
        let dmu =  z / s2
        let ds  = -T.one / sigma + (z*z) / (sigma*s2)
        return (dmu, ds)
    }

    /// Sum of Normal scores across data.
    public func totalScore(data: [T], mu: T, sigma: T) -> (dmu: T, dsigma: T) {
        var sMu = T.zero, sS = T.zero
        for x in data {
            let (a,b) = score(x: x, mu: mu, sigma: sigma)
            sMu += a; sS += b
        }
        return (sMu, sS)
    }

    /// Score with log-scale σ = exp(logSigma).
    ///
    /// - Chain rule: ∂/∂log σ = σ ∂/∂σ.
    public func scoreWithLogSigma(x: T, mu: T, logSigma: T) -> (dmu: T, dlogSigma: T) {
        let sigma = T.exp(logSigma)
        let (dmu, ds) = score(x: x, mu: mu, sigma: sigma)
        return (dmu, sigma * ds)
    }

    /// Total score with log-scale across data.
    public func totalScoreWithLogSigma(data: [T], mu: T, logSigma: T) -> (dmu: T, dlogSigma: T) {
        var sMu = T.zero, sL = T.zero
        let sigma = T.exp(logSigma)
        for x in data {
            let (a,b) = score(x: x, mu: mu, sigma: sigma)
            sMu += a; sL += sigma * b
        }
        return (sMu, sL)
    }

}
