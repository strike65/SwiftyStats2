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

/// Score helpers for the Holtsmark distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Holtsmark distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Holtsmark {

    /// Heuristic central-difference step for x-derivatives, scaled by magnitude and sigma.
    ///
    /// - Stability: uses √eps scaling to balance truncation and rounding errors.
    private func diffStep(_ x: T, scale: T) -> T {
        // central-difference step ~ sqrt(eps) * (|x| + scale)
        let eps = T.ulpOfOne
        return T.sqrt(eps) * (x.magnitude + T.maximum(scale, T.one))
    }
    
    /// Finite-difference approximation of ∂/∂x log f(x | μ, σ) using central differences with
    /// one step of Richardson refinement for improved accuracy.
    ///
    /// - Numerical stability:
    ///   - Step size chosen via √eps scaling.
    ///   - Richardson refinement reduces O(h²) error.
    private func dLogPdf_dx(
        _ x: T, mu: T, sigma: T
    ) -> T {
        var h = diffStep(x, scale: sigma)
        do {
            // central difference (coarse)
            let fP: T = try self.logPdf(x + h)
            let fM: T = try self.logPdf(x - h)
            var deriv = (fP - fM) / (T.two * h)
            // Optional: one step Richardson refinement for more accuracy
            h *= T.half
            let fP2: T = try self.logPdf(x + h)
            let fM2: T = try self.logPdf(x - h)
            let derivFine: T = (fP2 - fM2) / (T.two * h)
            // Combine O(h^2) terms: d ≈ d_fine + (d_fine - d_coarse)/3
            deriv = derivFine + (derivFine - deriv) / 3.0
            return deriv
        }
        catch _ {
            return T.nan
        }
    }
    
    /// Score (∂ℓ/∂μ, ∂ℓ/∂σ) for Holtsmark(μ, σ) using scale-location derivative identities:
    ///
    /// - For location-scale families: if ψ(z) = σ ∂/∂x log f(x | μ, σ) at x, with z = (x − μ)/σ,
    ///   then:
    ///   - ∂ℓ/∂μ = −ψ/σ
    ///   - ∂ℓ/∂σ = (−1 − z ψ)/σ
    /// - This method estimates ψ via finite differences of logPdf.
    /// - Constraints: σ > 0.
    public func score(
        x: T, mu: T, sigma: T
    ) -> (dmu: T, dsigma: T) {
        precondition(sigma > 0, "sigma must be > 0")
        let z: T = (x - mu) / sigma

        // ψ(z) = σ * (∂/∂x) log f(x | μ, σ) at x
        let dlogf_dx: T = dLogPdf_dx(x, mu: mu, sigma: sigma)
        let psi: T = sigma * dlogf_dx

        let dmu: T = -(psi) / sigma
        let dsigma: T = (-T.one - z * psi) / sigma
        return (dmu, dsigma)
    }
    
    /// Sum of Holtsmark scores across data points.
    public func totalScore(
        data: [T], mu: T, sigma: T
    ) -> (dmu: T, dsigma: T) {
        var sMu: T = T.zero, sSigma: T = T.zero
        for x in data {
            let (dmu, dsigma) = self.score(x: x, mu: mu, sigma: sigma)
            sMu += dmu; sSigma += dsigma
        }
        return (sMu, sSigma)
    }

    /// Score with log-scale σ = exp(logSigma).
    ///
    /// - Chain rule: ∂/∂log σ = σ ∂/∂σ.
    public func scoreWithLogSigma(
        x: T, mu: T, logSigma: T
    ) -> (dmu: T, dlogSigma: T) {
        let sigma: T = T.exp(logSigma)
        let (dmu, dsigma) = score(x: x, mu: mu, sigma: sigma)
        return (dmu, sigma * dsigma)
    }

    /// Total score with log-scale across data.
    public func totalScoreWithLogSigma(
        data: [T], mu: T, logSigma: T
    ) -> (dmu: T, dlogSigma: T) {
        var sMu: T = T.zero, sLogS: T = T.zero
        let sigma: T = T.exp(logSigma)
        for x in data {
            let (dmu, dsigma) = self.score(x: x, mu: mu, sigma: sigma)
            sMu += dmu; sLogS += sigma * dsigma
        }
        return (sMu, sLogS)
    }
}
