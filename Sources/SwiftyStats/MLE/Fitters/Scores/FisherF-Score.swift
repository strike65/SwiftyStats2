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

/// Score helpers for the FisherF distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude


/// Score helpers for the FisherF distribution used by MLEFitter.
extension SwiftyBoost.Distribution.FisherF {

    /// Score for Fisher-Snedecor F distribution F(d1, d2) with x > 0 and d1, d2 > 0.
    ///
    /// - Density: f(x) ∝ x^{d1/2 − 1} (1 + (d1/d2)x)^{−(d1+d2)/2}.
    /// - Log-likelihood (up to additive constants in d1,d2):
    ///   ℓ = (d1/2 − 1) log x − (d1 + d2)/2 log(1 + (d1/d2) x) + const(d1,d2).
    /// - Score: closed forms implemented in numerically stable grouped terms for d1 and d2.
    /// - Constraints: x > 0, d1 > 0, d2 > 0.
    public func score(x: T, d1: T, d2: T) -> (dd1: T, dd2: T) {
        precondition(x > 0 && d1 > 0 && d2 > 0, "x>0, d1>0, d2>0")

        let h = T.half
        let a = d1 * h                 // d1/2
        let b = d2 * h                 // d2/2
        let r = d1 / d2                // d1/d2

        let log1prx = T.log(1 + r * x)
        let psiA  = Distribution.dg(a)
        let psiB  = Distribution.dg(b)
        let psiAB = Distribution.dg(a + b)

        // --- dd1 parts (split to keep ≤3 terms each) ---
        let t1 = h * (T.log(d1) - T.log(d2))
        let t2 = a / d1
        let t3 = h * T.log(x)

        let t4 = -h * (psiA - psiAB)
        let t5 = -h * log1prx

        let denom = d2 * (1 + r * x)
        let frac  = x / denom
        let t6 = -(a + b) * frac

        let s1 = t1 + t2
        let s2 = t3 + t4
        let s3 = t5 + t6
        let dd1 = s1 + s2 + s3

        // --- dd2 parts ---
        let u1 = -(a / d2)
        let u2 = -h * (psiB - psiAB)
        let u3 = -h * log1prx
        let u4 = (a + b) * ( (r * x) / denom )

        let v1 = u1 + u2
        let v2 = u3 + u4
        let dd2 = v1 + v2

        return (dd1, dd2)
    }

    /// Sum of Fisher F scores across data (d1, d2).
    public func totalScore(data: [T], d1: T, d2: T) -> (dd1: T, dd2: T) {
        precondition(d1 > 0 && d2 > 0, "d1>0, d2>0")

        let h = T.half
        let a = d1 * h
        let b = d2 * h
        let r = d1 / d2

        let logD1 = T.log(d1)
        let logD2 = T.log(d2)
        let psiA  = Distribution.dg(a)
        let psiB  = Distribution.dg(b)
        let psiAB = Distribution.dg(a + b)

        var S1 = T.zero, S2 = T.zero
        for x in data {
            precondition(x > 0, "x>0")

            let logx = T.log(x)
            let log1prx = T.log(1 + r * x)
            let denom = d2 * (1 + r * x)
            let frac  = x / denom

            // dd1 per x
            let t1 = h * (logD1 - logD2)
            let t2 = a / d1
            let t3 = h * logx
            let t4 = -h * (psiA - psiAB)
            let t5 = -h * log1prx
            let t6 = -(a + b) * frac
            let s1 = t1 + t2
            let s2 = t3 + t4
            let s3 = t5 + t6
            let d1_i = s1 + s2 + s3

            // dd2 per x
            let u1 = -(a / d2)
            let u2 = -h * (psiB - psiAB)
            let u3 = -h * log1prx
            let u4 = (a + b) * ((r * x) / denom)
            let v1 = u1 + u2
            let v2 = u3 + u4
            let d2_i = v1 + v2

            S1 += d1_i
            S2 += d2_i
        }
        return (S1, S2)
    }

    /// Score with log-parameters d1 = exp(logD1), d2 = exp(logD2).
    ///
    /// - Chain rule: ∂/∂log d1 = d1 ∂/∂d1, ∂/∂log d2 = d2 ∂/∂d2.
    public func scoreWithLogParams(x: T, logD1: T, logD2: T) -> (dlogD1: T, dlogD2: T) {
        let d1 = T.exp(logD1)
        let d2 = T.exp(logD2)
        let (u, v) = score(x: x, d1: d1, d2: d2)
        let w1 = d1 * u
        let w2 = d2 * v
        return (w1, w2)
    }

    /// Total score with log-parameters across data.
    public func totalScoreWithLogParams(data: [T], logD1: T, logD2: T) -> (dlogD1: T, dlogD2: T) {
        let d1 = T.exp(logD1)
        let d2 = T.exp(logD2)
        let (U, V) = totalScore(data: data, d1: d1, d2: d2)
        let W1 = d1 * U
        let W2 = d2 * V
        return (W1, W2)
    }
}

