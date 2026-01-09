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

/// Score helpers for the Triangular distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Triangular distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Triangular {
    /// Numerically stable sigmoid for logit parameterization of τ ∈ [0,1].
    ///
    /// - Stability: avoids overflow for large |t|.
    @inline(__always) private func sigmoidStable(_ t: T) -> T {
        if t >= T.zero {
            let e = T.exp(-t)
            return T.one / (T.one + e)
        } else {
            let e = T.exp(t)
            return e / (T.one + e)
        }
    }

    /// Score for Triangular(x | a, b, c), with a < b, a ≤ c ≤ b, and x ∈ [a, b].
    ///
    /// - Density: piecewise linear; mode at c.
    /// - Log-likelihood derivatives (closed form), split by x ≤ c and x > c.
    /// - Constraints: a < b, a ≤ c ≤ b, a ≤ x ≤ b.
    @inline(__always)
    public func score(x: T, lower: T, upper: T, mode: T) -> (da: T, db: T, dc: T) {
        precondition(lower < upper && lower <= mode && mode <= upper && x >= lower && x <= upper, "a<b, a<=c<=b, a<=x<=b")
        if x <= mode {
            let da:T = -T.one/(x - lower) + T.one/(upper - lower) + T.one/(mode - lower)
            let db = -T.one/(upper - lower)
            let dc = -T.one/(mode - lower)
            return (da, db, dc)
        } else {
            let da =  T.one/(upper - lower)
            let db =  T.one/(upper - x) - T.one/(upper - lower) - T.one/(upper - mode)
            let dc =  T.one/(upper - mode)
            return (da, db, dc)
        }
    }
    
    /// Sum of Triangular scores across data.
    @inline(__always)
    public func totalScore(data: [T], a: T, b: T, c: T) -> (da: T, db: T, dc: T) {
        var A = T.zero, B = T.zero, C = T.zero
        for x in data {
            let (u,v,w) = score(x: x, lower: a, upper: b, mode: c)
            A += u; B += v; C += w
        }
        return (A,B,C)
    }
    
    /// Score with reparameterization via m = a, w = b − a > 0, τ = (c − a)/(b − a) ∈ [0,1].
    ///
    /// - Optimization parameters: (m, logW, logitTau) with w = exp(logW), τ = sigmoid(logitTau).
    /// - Chain rule:
    ///   - a = m
    ///   - b = m + w        ⇒ ∂b/∂m = 1, ∂b/∂w = 1
    ///   - c = m + w τ      ⇒ ∂c/∂m = 1, ∂c/∂w = τ, ∂c/∂τ = w
    ///   - ∂/∂logW = w ∂/∂w, ∂/∂logitTau = τ(1 − τ) ∂/∂τ
    /// - Numerical stability: stable sigmoid for τ.
    @inline(__always)
    public func scoreWithReparam(x: T, m: T, logW: T, logitTau: T) -> (dm: T, dlogW: T, dlogitTau: T) {
        let w = T.exp(logW)
        let tau = sigmoidStable(logitTau)
        let a = m
        let b = m + w
        let c = m + w * tau
        let (da, db, dc) = score(x: x, lower: a, upper: b, mode: c)
        // a=m → ∂a/∂m=1
        // b=m+w → ∂b/∂m=1, ∂b/∂w=1
        // c=m+w*tau → ∂c/∂m=1, ∂c/∂w=tau, ∂c/∂tau=w
        let dm = da + db + dc
        let dw = db + dc * tau
        let dtau = dc * w
        let dlogW = w * dw
        let dlogitTau = tau * (T.one - tau) * dtau
        return (dm, dlogW, dlogitTau)
    }
    
    /// Total score with (m, logW, logitTau) across data.
    @inline(__always)
    public func totalScoreWithReparam(data: [T], m: T, logW: T, logitTau: T) -> (dm: T, dlogW: T, dlogitTau: T) {
        var M = T.zero, L = T.zero, R = T.zero
        let w = T.exp(logW)
        let tau = sigmoidStable(logitTau)
        let a = m
        let b = m + w
        let c = m + w * tau
        for x in data {
            let (da, db, dc) = score(x: x, lower: a, upper: b, mode: c)
            M += (da + db + dc)
            L += w * (db + dc * tau)
            R += tau * (T.one - tau) * (dc * w)
        }
        return (M, L, R)
    }
}
