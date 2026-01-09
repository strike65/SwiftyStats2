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
/// Score helpers for the Bernoulli distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude


/// Score helpers for the Bernoulli distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Bernoulli {
    /// Score for Bernoulli(x | p), x ∈ {0,1}, p ∈ (0,1).
    ///
    /// - Log-likelihood: ℓ = x log p + (1 − x) log(1 − p).
    /// - Score: ∂ℓ/∂p = x/p − (1 − x)/(1 − p).
    /// - Constraints: p ∈ (0,1), x ∈ {0,1}.
    public func score(x: T, p: T) -> T {
        precondition(p > 0 && p < T.one && (x == 0 || x == T.one), "valid p and x in {0,1}")
        return ( x / p - (T.one - x) / (T.one - p) )
    }

    /// Sum of Bernoulli scores across the dataset.
    public func totalScore(data: [T], p: T) -> T {
        var s = T.zero
        for x in data { s += score(x: x, p: p) }
        return (s)
    }

    /// Numerically stable sigmoid to map logit to probability without overflow.
    ///
    /// - Stability: avoids overflow for large |t| by switching the evaluation form.
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
    /// - Chain rule: ∂/∂η = p(1 − p) ∂/∂p; with Bernoulli score this simplifies to x − p.
    /// - Stability: uses `sigmoidStable` to avoid overflow in exp.
    public func scoreWithLogitP(x: T, logitP: T) -> T {
        let p = sigmoidStable(logitP)
        return (x - p)
    }

    /// Total score with logit parameterization across data.
    public func totalScoreWithLogitP(data: [T], logitP: T) -> T {
        let p = sigmoidStable(logitP)
        var s = T.zero
        for x in data { s += (x - p) }
        return (s)
    }
}
