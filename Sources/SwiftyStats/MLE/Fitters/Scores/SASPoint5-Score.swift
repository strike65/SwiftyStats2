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

/// Score helpers for the SASPoint5 distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

// ============================================================
// SASPoint5 (two-parameter variant): score functions
// Parameters: location mu (R), scale sigma > 0
// Implementation: numerical score via central differences on log-pdf.
// No fallback to Double; uses only generic T operations.
// Chain-rule variants for log-scale are provided.
// ============================================================
/// Score helpers for the SASPoint5 distribution used by MLEFitter.
extension SwiftyBoost.Distribution.SASPoint5 {
    // Score w.r.t. mu and sigma (two-parameter interface)
    // Requires: self.logPdf(x, mu:sigma:) to be available.
    /// Central-difference gradient of the SAS 0.5 log-density at a single observation.
    public func score(
        x: T, mu: T, sigma: T
    ) -> (dmu: T, dsigma: T) {
        precondition(sigma > 0, "sigma>0")

        let f: (T, T) -> T = { m, s in
            do { return try self.logPdf(x) }
            catch { return T.nan }
        }
        let hm = _sasp5_diffStep(mu)
        let hs = _sasp5_diffStep(sigma)
        let (dm, ds) = _sasp5_cdiff2(f, mu, sigma, hm, hs)
        return (dm, ds)
    }

    // Sum of scores across data
    /// Accumulates SAS 0.5 scores across a dataset.
    public func totalScore(
        data: [T], mu: T, sigma: T
    ) -> (dmu: T, dsigma: T) {
        precondition(sigma > 0, "sigma>0")
        var sMu = T.zero, sSi = T.zero
        for xi in data {
            let g = self.score(x: xi, mu: mu, sigma: sigma)
            sMu += g.dmu
            sSi += g.dsigma
        }
        return (sMu, sSi)
    }

    // Log-scale score for sigma: let logSigma be unconstrained, sigma = exp(logSigma)
    // Chain rule: d/d logSigma = sigma * d/d sigma
    /// Convenience variant where σ is parameterized via an unconstrained log-scale.
    public func scoreWithLogSigma(
        x: T, mu: T, logSigma: T
    ) -> (dmu: T, dlogSigma: T) {
        let sigma = T.exp(logSigma)
        let (dm, ds) = score(x: x, mu: mu, sigma: sigma)
        return (dm, sigma * ds)
    }

    /// Sum of log-scale SAS 0.5 scores across an entire dataset.
    public func totalScoreWithLogSigma(
        data: [T], mu: T, logSigma: T
    ) -> (dmu: T, dlogSigma: T) {
        let sigma = T.exp(logSigma)
        var sMu = T.zero, sLs = T.zero
        for xi in data {
            let (dm, ds) = score(x: xi, mu: mu, sigma: sigma)
            sMu += dm
            sLs += sigma * ds
        }
        return (sMu, sLs)
    }
    @inline(__always)
    private func _sasp5_diffStep(_ scale: T) -> T {
        // central-difference step ~ sqrt(ulp) * max(1, |scale|)
        let eps = T.ulpOfOne
        return T.sqrt(eps) * T.maximum(T.one, scale.magnitude)
    }

    @inline(__always)
    private func _sasp5_cdiff2(
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
