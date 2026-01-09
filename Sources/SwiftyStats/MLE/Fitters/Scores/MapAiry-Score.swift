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

/// Score helpers for the MapAiry distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the MapAiry distribution used by MLEFitter.
extension SwiftyBoost.Distribution.MapAiry {
    // Example params: location mu, scale sigma>0 (adjust to your API)
    /// Central-difference gradient for the Map Airy density at a single sample.
    public func score(x: T, mu: T, sigma: T) -> (dmu: T, dsigma: T) {
        precondition(sigma > 0)
        let f: (T, T) -> T = { m, s in
            do { return try self.logPdf(x) } catch { return T.nan }
        }
        let hm = diffStep(mu), hs = diffStep(sigma)
        let (dm, ds) = centralDiff2(f, mu, sigma, hm, hs)
        return (dm, ds)
    }
    /// Dataset-level Map Airy score obtained by summing per-point gradients.
    public func totalScore(data: [T], mu: T, sigma: T) -> (dmu: T, dsigma: T) {
        var sm=T.zero, ss=T.zero
        for v in data {
            let g = self.score(x: v, mu: mu, sigma: sigma)
            sm += g.dmu; ss += g.dsigma
        }
        return (sm, ss)
    }

    /// Score with log-scale reparameterization `sigma = exp(logSigma)`.
    ///
    /// - Chain rule: ∂ℓ/∂logSigma = sigma · ∂ℓ/∂sigma.
    public func scoreWithLog(x: T, mu: T, logSigma: T) -> (dmu: T, dlogSigma: T) {
        let sigma = T.exp(logSigma)
        let (dmu, dsigma) = score(x: x, mu: mu, sigma: sigma)
        return (dmu, sigma * dsigma)
    }

    /// Dataset-level score with the log-scale parameterization.
    public func totalScoreWithLog(data: [T], mu: T, logSigma: T) -> (dmu: T, dlogSigma: T) {
        let sigma = T.exp(logSigma)
        var sm = T.zero, sl = T.zero
        for v in data {
            let (dmu, dsigma) = score(x: v, mu: mu, sigma: sigma)
            sm += dmu
            sl += sigma * dsigma
        }
        return (sm, sl)
    }
    
    @inline(__always)
    private func diffStep(_ scale: T) -> T {
        // central-difference step ~ sqrt(ulp) * max(1, |scale|)
        let eps = T.ulpOfOne
        return T.sqrt(eps) * T.maximum(T.one, scale.magnitude)
    }
    
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

internal extension SwiftyBoost.Distribution {
    @inline(__always)
    static func dg<T: Real & BinaryFloatingPoint & Sendable>(_ x: T) -> T {
        do {
            return try SwiftyBoost.SpecialFunctions.digamma(x)
        }
        catch _ {
            return T.nan
        }
    }
}
































