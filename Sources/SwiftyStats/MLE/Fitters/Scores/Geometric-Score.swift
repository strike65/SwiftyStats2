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

/// Score helpers for the Geometric distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude


/// Score helpers for the Geometric distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Geometric {
    /// Score for Geometric(k | p) with support k ∈ {0,1,2,...}, success prob p ∈ (0,1).
    ///
    /// - PMF (failures before first success): P(K=k) = (1 − p)^k p.
    /// - Log-likelihood: ℓ = k log(1 − p) + log p.
    /// - Score: ∂ℓ/∂p = 1/p − k/(1 − p).
    public func score(k: T, p: T) -> T {
        precondition(k >= 0 && p > 0 && p < T.one, "k>=0, p in (0,1)")
        return ( T.one / p - k / (T.one - p) )
    }

    /// Sum of Geometric scores across the dataset.
    public func totalScore(data: [T], p: T) -> T {
        var s = T.zero
        for k in data { s += score(k: k, p: p) }
        return (s)
    }

    /// Numerically stable sigmoid for logit parameterization.
    @inline(__always) private func sigmoidStable(_ t: T) -> T {
        if t >= T.zero {
            let e = T.exp(-t)
            return T.one / (T.one + e)
        } else {
            let e = T.exp(t)
            return e / (T.one + e)
        }
    }

    /// Score with logit reparameterization p = sigmoid(η).
    ///
    /// - Chain rule: ∂/∂η = p(1 − p) ∂/∂p = 1 − p (k + 1).
    public func scoreWithLogitP(k: T, logitP: T) -> T {
        let p = sigmoidStable(logitP)
        return ( T.one - p * (k + T.one) )
    }

    /// Total score with logit parameterization across data.
    public func totalScoreWithLogitP(data: [T], logitP: T) -> T {
        let p = sigmoidStable(logitP)
        var s = T.zero
        for k in data { s += (T.one - p * (k + T.one)) }
        return (s)
    }
}
