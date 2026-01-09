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

/// Score helpers for the Rayleigh distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Rayleigh distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Rayleigh {
    /// Score for Rayleigh(x | σ), x ≥ 0, σ > 0.
    ///
    /// - Density: f(x) = x/σ² exp(−x²/(2σ²)).
    /// - Log-likelihood: ℓ = log x − 2 log σ − x²/(2σ²).
    /// - Score: ∂ℓ/∂σ = −2/σ + x²/σ³.
    public func score(x: T, sigma: T) -> T {
        precondition(sigma > 0 && x >= 0, "sigma>0, x>=0")
        return (-2 / sigma + (x*x) / (sigma*sigma*sigma))
    }

    /// Sum of Rayleigh scores across data.
    public func totalScore(data: [T], sigma: T) -> T {
        var s = T.zero
        for x in data { s += score(x: x, sigma: sigma) }
        return (s)
    }

    /// Score with log-scale σ = exp(logSigma).
    ///
    /// - Chain rule: ∂/∂log σ = σ ∂/∂σ ⇒ −2 + (x/σ)².
    public func scoreWithLogSigma(x: T, logSigma: T) -> T {
        let sigma = T.exp(logSigma)
        return sigma * score(x: x, sigma: sigma)   // = -2 + (x/sigma)^2
    }

    /// Total score with log-scale across data.
    public func totalScoreWithLogSigma(data: [T], logSigma: T) -> T {
        let sigma = T.exp(logSigma)
        var s = T.zero
        for x in data { s += sigma * score(x: x, sigma: sigma) }
        return (s)
    }
}
