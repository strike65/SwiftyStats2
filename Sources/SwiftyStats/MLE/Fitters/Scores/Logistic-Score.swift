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

/// Score helpers for the Logistic distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Logistic distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Logistic {

    /// Score for Logistic(x | μ, s), s > 0.
    ///
    /// - Density: f(x) = exp(−t) / (s (1 + exp(−t))²), t = (x − μ)/s.
    /// - Log-likelihood: ℓ = −log s − 2 log(1 + e^{−t}) − t.
    /// - Using tanh(t/2) = (e^t − 1)/(e^t + 1) for stable expressions.
    /// - Score:
    ///   - ∂ℓ/∂μ = tanh(t/2)/s
    ///   - ∂ℓ/∂s = (−1 + t tanh(t/2))/s
    /// - Constraints: s > 0.
    public func score(x: T, mu: T, s: T) -> (dmu: T, ds: T) {
        precondition(s > 0, "s must be > 0")
        let t = (x - mu) / s
        // tanh(t/2) = (e^t - 1)/(e^t + 1)
        let e = T.exp(t)
        let tanhHalf = (e - T.one) / (e + T.one)
        let dmu =  tanhHalf / s
        let ds  = (-T.one + t * tanhHalf) / s
        return (dmu, ds)
    }

    /// Sum of Logistic scores across data.
    public func totalScore(data: [T], mu: T, s: T) -> (dmu: T, ds: T) {
        var a = T.zero, b = T.zero
        for x in data {
            let (u,v) = score(x: x, mu: mu, s: s)
            a += u; b += v
        }
        return (a,b)
    }

    /// Score with log-scale s = exp(logS).
    ///
    /// - Chain rule: ∂/∂log s = s ∂/∂s.
    public func scoreWithLogS(x: T, mu: T, logS: T) -> (dmu: T, dlogS: T) {
        let s = T.exp(logS)
        let (dmu, ds) = score(x: x, mu: mu, s: s)
        return (dmu, s * ds)
    }

    /// Total score with log-scale across data.
    public func totalScoreWithLogS(data: [T], mu: T, logS: T) -> (dmu: T, dlogS: T) {
        var a = T.zero, b = T.zero
        let s = T.exp(logS)
        for x in data {
            let (u,v) = score(x: x, mu: mu, s: s)
            a += u; b += s * v
        }
        return (a,b)
    }
}
