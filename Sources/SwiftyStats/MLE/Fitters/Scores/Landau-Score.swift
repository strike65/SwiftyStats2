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

/// Score helpers for the Landau distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Landau distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Landau {
    // Typical 2-parameter: location mu, scale c>0
    /// Central-difference estimate of the Landau log-likelihood gradient at a single sample.
    public func score(x: T, mu: T, scale: T) -> (dmu: T, dc: T) {
        precondition(scale > 0)
        let f: (T, T) -> T = { m, s in
            do {
                return try self.logPdf(x)
            }
            catch {
                return T.nan
            }
        }
        let hm = diffStep(mu), hc = diffStep(scale)
        let (dm, dc) = centralDiff2(f, mu, scale, hm, hc)
        return (dm, dc)
    }
    /// Aggregate score for Landau parameters obtained by summing per-sample gradients.
    public func totalScore(data: [T], mu: T, scale: T) -> (dmu: T, dc: T) {
        var sm=T.zero, sc=T.zero
        for v in data {
            let g = self.score(x: v, mu: mu, scale: scale)
            sm += g.dmu; sc += g.dc
        }
        return (sm, sc)
    }

    /// Score with log-scale reparameterization: `scale = exp(logScale)`.
    ///
    /// - Chain rule: ∂ℓ/∂logScale = scale · ∂ℓ/∂scale.
    public func scoreWithLog(x: T, mu: T, logScale: T) -> (dmu: T, dlogScale: T) {
        let scale = T.exp(logScale)
        let (dmu, dscale) = score(x: x, mu: mu, scale: scale)
        return (dmu, scale * dscale)
    }

    /// Sum of log-scale scores across the dataset.
    public func totalScoreWithLog(data: [T], mu: T, logScale: T) -> (dmu: T, dlogScale: T) {
        let scale = T.exp(logScale)
        var sMu = T.zero, sLog = T.zero
        for v in data {
            let (dmu, dscale) = score(x: v, mu: mu, scale: scale)
            sMu += dmu
            sLog += scale * dscale
        }
        return (sMu, sLog)
    }
    @inline(__always)
    private func diffStep(_ scale: T) -> T {
        // central-difference step ~ sqrt(ulp) * max(1, |scale|)
        let eps = T.ulpOfOne
        return T.sqrt(eps) * T.maximum(T.one, scale.magnitude)
    }

    // Two-parameter central difference for log-likelihood gradients
    @inline(__always)
    private func centralDiff2(
        _ f: (T, T) -> T, _ a: T, _ b: T, _ ha: T, _ hb: T
    ) -> (da: T, db: T) {
        let fa_p = f(a + ha, b)
        let fa_m = f(a - ha, b)
        let da = (fa_p - fa_m) / (T(2) * ha)

        let fb_p = f(a, b + hb)
        let fb_m = f(a, b - hb)
        let db = (fb_p - fb_m) / (T(2) * hb)

        return (da, db)
    }
}
