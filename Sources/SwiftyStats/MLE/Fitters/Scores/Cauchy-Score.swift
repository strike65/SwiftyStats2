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

/// Score helpers for the Cauchy distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude


/// Score helpers for the Cauchy distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Cauchy {
    /// Score for Cauchy(x | μ, σ), σ > 0.
    ///
    /// - Density: f(x) = 1/(πσ (1 + z²)), z = (x − μ)/σ.
    /// - Log-likelihood: ℓ = −log(πσ) − log(1 + z²).
    /// - Score:
    ///   - ∂ℓ/∂μ = 2z / (σ (1 + z²))
    ///   - ∂ℓ/∂σ = (z² − 1) / (σ (1 + z²))
    /// - Constraints: σ > 0.
    public func score(x: T, mu: T, sigma: T) -> (dmu: T, dsigma: T) {
        precondition(sigma > 0, "sigma must be > 0")
        let z: T = (x - mu) / sigma
        let denom: T = T.one + z * z
        let dmu: T = (T.two * z) / (sigma * denom)
        let ds  = (z * z - T.one) / (sigma * denom)
        return (dmu, ds)
    }

    /// Sum of Cauchy scores across the dataset.
    public func totalScore(data: [T], mu: T, sigma: T) -> (dmu: T, dsigma: T) {
        precondition(sigma > 0, "sigma must be > 0")
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

    /// Sum of scores with log-scale across data.
    public func totalScoreWithLogSigma(data: [T], mu: T, logSigma: T) -> (dmu: T, dlogSigma: T) {
        let sigma = T.exp(logSigma)
        var sMu = T.zero, sL = T.zero
        for x in data {
            let (a,b) = score(x: x, mu: mu, sigma: sigma)
            sMu += a; sL += sigma * b
        }
        return (sMu, sL)
    }
}
