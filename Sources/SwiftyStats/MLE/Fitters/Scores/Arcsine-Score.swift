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
/// Score helpers for the Arcsine distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// MARK: Extend SwiftyBoost.Distribution to provide Score-functions if needed

/// Score helpers for the Arcsine distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Arcsine {
    /// Score for the Arcsine distribution on (a, b) with density f(x) = 1/(π√((x-a)(b-x))) for a < x < b.
    ///
    /// - Distribution and parameterization:
    ///   - Support: x ∈ (a, b), with a < b.
    ///   - Parameters: lower bound a, upper bound b.
    /// - Log-likelihood (single x): ℓ = -log π - 0.5 log(x - a) - 0.5 log(b - x).
    /// - Score definition:
    ///   - ∂ℓ/∂a =  1/(2(x - a))
    ///   - ∂ℓ/∂b = -1/(2(b - x))
    /// - Parameters/data constraints:
    ///   - a < b, a < x < b.
    public func score(x: T, a: T, b: T) -> (da: T, db: T) {
        precondition(a < b && x > a && x < b, "a<b and a<x<b")
        let da =  T.one / (2 * (x - a))
        let db = -T.one / (2 * (b - x))
        return (da, db)
    }

    /// Sum of scores for the Arcsine distribution across a dataset.
    ///
    /// - Parameters: see `score(x:a:b:)`.
    /// - Returns: Sum over observations of (∂ℓ/∂a, ∂ℓ/∂b).
    public func totalScore(data: [T], a: T, b: T) -> (da: T, db: T) {
        var A = T.zero, B = T.zero
        for x in data {
            let (u,v) = score(x: x, a: a, b: b)
            A += u; B += v
        }
        return (A,B)
    }

    /// Score under reparameterization by center and log-width:
    /// a = m − w/2, b = m + w/2, w = exp(logW).
    ///
    /// - Reparameterization and chain rule:
    ///   - ∂ℓ/∂m = ∂ℓ/∂a * ∂a/∂m + ∂ℓ/∂b * ∂b/∂m = da + db
    ///   - ∂ℓ/∂w = ∂ℓ/∂a * (−1/2) + ∂ℓ/∂b * (1/2)
    ///   - ∂ℓ/∂logW = w * ∂ℓ/∂w
    /// - Numerical stability: requires w = exp(logW) > 0 and finite.
    public func scoreWithCenterLogWidth(x: T, m: T, logW: T) -> (dm: T, dlogW: T) {
        let w = T.exp(logW)
        precondition(w > 0 && w.isFinite, "exp(logW) must be finite and > 0")
        let a = m - w / 2
        let b = m + w / 2
        let (da, db) = score(x: x, a: a, b: b)
        // a = m - w/2 → ∂a/∂m = 1, ∂a/∂w = -1/2
        // b = m + w/2 → ∂b/∂m = 1, ∂b/∂w =  1/2
        let dm = da + db
        let dw = (-da + db) / 2
        let dlogW = w * dw
        return (dm, dlogW)
    }

    /// Total score under center/log-width reparameterization.
    ///
    /// - See `scoreWithCenterLogWidth(x:m:logW:)`.
    public func totalScoreWithCenterLogWidth(data: [T], m: T, logW: T) -> (dm: T, dlogW: T) {
        var M = T.zero, L = T.zero
        let w = T.exp(logW)
        precondition(w > 0 && w.isFinite, "exp(logW) must be finite and > 0")
        let a = m - w / 2
        let b = m + w / 2
        for x in data {
            let (da, db) = score(x: x, a: a, b: b)
            M += (da + db)
            L += w * ((-da + db) / 2)
        }
        return (M, L)
    }
}

