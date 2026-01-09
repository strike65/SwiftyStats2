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

/// Score helpers for the SkewNormal distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the SkewNormal distribution used by MLEFitter.
extension SwiftyBoost.Distribution.SkewNormal {
    /// Standard Normal PDF φ(t).
    ///
    /// - Parameter t: Evaluation point.
    /// - Returns: φ(t), or `.nan` if the internal normal construction fails.
    @inline(__always)
    private func stdNormalPDF(_ t: T) -> T {
        do {
            let sn: SwiftyBoost.Distribution.Normal<T> = try SwiftyBoost.Distribution.Normal()
            return try sn.pdf(t)
        }
        catch _ {
            return .nan
        }
    }

    /// Standard Normal CDF Φ(t).
    ///
    /// - Parameter t: Evaluation point.
    /// - Returns: Φ(t), or `.nan` if the internal normal construction fails.
    @inline(__always)
    private func stdNormalCDF(_ t: T) -> T {
        do {
            let sn: SwiftyBoost.Distribution.Normal<T> = try SwiftyBoost.Distribution.Normal()
            return try sn.cdf(t)
        }
        catch _ {
            return .nan
        }
    }

    /// Numerically stable Mills ratio λ(t) = φ(t)/Φ(t).
    ///
    /// - Parameter t: Evaluation point.
    /// - Returns: λ(t). For very negative `t`, uses an asymptotic expansion to avoid underflow.
    /// - Numerical stability: avoids catastrophic cancellation when Φ(t) ≈ 0.
    @inline(__always)
    private func millsRatio(_ t: T) -> T {
        if t < T(-10) {
            // Asymptotic expansion for large negative t to avoid catastrophic cancellation.
            let inv: T = T.one / t
            let inv2: T = inv * inv
            let inv3: T = inv2 * inv
            let inv5: T = inv3 * inv2
            let inv7: T = inv5 * inv2
            let ex1: T = T(2) * inv3
            let ex2: T = T(5) * inv5
            let ex3: T = T(14) * inv7
            return (-t) - inv + ex1 - ex2 + ex3
        } else {
            let phi: T = stdNormalPDF(t)
            let Phi: T = stdNormalCDF(t)
            let denom: T = T.maximum(Phi, T.leastNonzeroMagnitude)
            return phi / denom
        }
    }

    /// Pointwise score (gradient of log-likelihood) for Skew-Normal(μ, σ, α) at observation `x`.
    ///
    /// - Distribution and parameterization:
    ///   - X has density f(x) = (2/σ) φ(z) Φ(α z), with z = (x − μ)/σ.
    ///   - Parameters: μ ∈ ℝ, σ > 0, α ∈ ℝ.
    /// - Score (closed form):
    ///   - ∂ℓ/∂μ = (z − α λ)/σ
    ///   - ∂ℓ/∂σ = (−1 + z² − α z λ)/σ
    ///   - ∂ℓ/∂α = z λ
    ///   where λ = φ(α z)/Φ(α z) is the Mills ratio.
    /// - Constraints: σ > 0.
    /// - Numerical stability: uses stable Mills ratio evaluation.
    public func score(x: T, mu: T, sigma: T, alpha: T)
    -> (dmu: T, dsigma: T, dalpha: T)
    {
        precondition(sigma > 0, "sigma must be > 0")
        let z = (x - mu) / sigma
        let t = alpha * z
        let lambda = self.millsRatio(t)
        // d/dμ, d/dσ, d/dα from standard skew-normal score identities.
        let dmu    = ( z - alpha * lambda ) / sigma
        let dsigma = ( -T.one + z*z - alpha * z * lambda ) / sigma
        let dalpha = z * lambda
        return (dmu, dsigma, dalpha)
    }

    /// Sample total score for Skew-Normal(μ, σ, α) across data points.
    public func totalScore(data: [T], mu: T, sigma: T, alpha: T)
    -> (dmu: T, dsigma: T, dalpha: T)
    {
        var sMu: T = T.zero, sSigma: T = T.zero, sAlpha: T = T.zero
        for x in data {
            let (dmu, dsigma, dalpha) = score(x: x, mu: mu, sigma: sigma, alpha: alpha)
            sMu += dmu; sSigma += dsigma; sAlpha += dalpha
        }
        return (sMu, sSigma, sAlpha)
    }

    /// Score with log-scale parameterization (log σ) for improved numerical stability.
    ///
    /// - Reparameterization: σ = exp(logSigma).
    /// - Chain rule: ∂/∂log σ = σ ∂/∂σ.
    public func scoreWithLogSigma(x: T, mu: T, logSigma: T, alpha: T)
    -> (dmu: T, dlogSigma: T, dalpha: T)
    {
        let sigma = T.exp(logSigma)
        let (dmu, dsigma, dalpha) = score(x: x, mu: mu, sigma: sigma, alpha: alpha)
        // Chain rule: ∂/∂ log σ = σ * ∂/∂σ
        return (dmu, sigma * dsigma, dalpha)
    }
    
    /// Total score with log-scale parameterization across data.
    public func totalScoreWithLogSigma(data: [T], mu: T, logSigma: T, alpha: T)
        -> (dmu: T, dlogSigma: T, dalpha: T)
        {
            var sMu: T = T.zero, sLogS: T = T.zero, sAlpha: T = T.zero
            for x: T in data {
                let (dmu,dlogSigma,dalpha) = scoreWithLogSigma(x: x, mu: mu, logSigma: logSigma, alpha: alpha)
                sMu += dmu; sLogS += dlogSigma; sAlpha += dalpha
            }
            return (sMu, sLogS, sAlpha)
        }
}

