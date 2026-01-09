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

/// Score helpers for the NegativeBinomial distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the NegativeBinomial distribution used by MLEFitter.
extension SwiftyBoost.Distribution.NegativeBinomial {
    /// Score contribution for one observation parameterized by shape `r` and success probability `p`.
    public func score(x: T, r: T, p: T) -> (dr: T, dp: T) {
        precondition(r > 0 && p > 0 && p < 1 && x >= 0)
        let psi = SwiftyBoost.Distribution.dg(r)
        let psi2 = SwiftyBoost.Distribution.dg(x + r)
        let dr = psi2 - psi + T.log(p)
        let dp = r/p - x/(T.one - p)
        return (dr, dp)
    }
    /// Batch score for the negative binomial obtained by summing single-point scores.
    public func totalScore(data: [T], r: T, p: T) -> (dr: T, dp: T) {
        let logp = T.log(p)
        var sr = T.zero, sp = T.zero
        for v in data {
            let psi2 = SwiftyBoost.Distribution.dg(v + r)
            let psi = SwiftyBoost.Distribution.dg(r)
            sr += psi2 - psi + logp
            sp += r/p - v/(T.one - p)
        }
        return (sr, sp)
    }

    /// Score with log-parameterization `r = exp(logR)`.
    ///
    /// - Chain rule: ∂ℓ/∂logR = r · ∂ℓ/∂r.
    public func scoreWithLog(x: T, logR: T, p: T) -> (dlogR: T, dp: T) {
        let r = T.exp(logR)
        let (dr, dp) = score(x: x, r: r, p: p)
        return (r * dr, dp)
    }

    /// Batch score using the log-parameterization for `r`.
    public func totalScoreWithLog(data: [T], logR: T, p: T) -> (dlogR: T, dp: T) {
        let r = T.exp(logR)
        var sr = T.zero, sp = T.zero
        for v in data {
            let (dr, dp) = score(x: v, r: r, p: p)
            sr += r * dr
            sp += dp
        }
        return (sr, sp)
    }
}
