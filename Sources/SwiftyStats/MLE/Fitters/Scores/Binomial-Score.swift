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

/// Score helpers for the Binomial distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Binomial distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Binomial {
    /// Score for Binomial(x | n, p), with fixed n and 0 < p < 1.
    ///
    /// - PMF:     f(x) = C(n, x) p^x (1-p)^(n-x),  x = 0,...,n
    /// - Log-lik: l = x log p + (n - x) log(1 - p) + const
    /// - Score:
    ///   - dl/dp = x/p - (n - x)/(1 - p) = (x - n p)/(p (1 - p))
    /// - Constraints: n >= 0 (integer), 0 < p < 1
    public func score(x: T, n: T, p: T) -> T {
        precondition(n >= 0, "n must be >= 0")
        precondition(p > 0 && p < 1, "p must be in (0, 1)")
        let dp: T = (x - n * p) / (p * (T.one - p))
        return (dp)
    }
    
    /// Sum of Binomial scores across the dataset (same n and p for all x).
    public func totalScore(data: [T], n: T, p: T) -> T {
        precondition(n >= 0, "n must be >= 0")
        precondition(p > 0 && p < 1, "p must be in (0, 1)")
        var s: T = .zero
        let denom: T = p * (T.one - p)
        for x in data {
            s += (x - n * p) / denom
        }
        return (s)
    }
    
    /// Score with logit parameter eta = log(p/(1-p)).
    ///
    /// - Reparameterization: p = sigmoid(eta) = 1 / (1 + exp(-eta))
    /// - Chain rule: dl/deta = (dl/dp) * (dp/deta) = (x - n p)
    ///   (since dp/deta = p(1-p))
    public func scoreWithLogitP(x: T, n: T, logitP: T) -> T {
        precondition(n >= 0, "n must be >= 0")
        let p = T.one / (T.one + T.exp(-logitP))
        let dEta: T = x - n * p
        return (dEta)
    }
    
    /// Sum of scores with logit parameter across data.
    public func totalScoreWithLogitP(data: [T], n: T, logitP: T) -> T {
        precondition(n >= 0, "n must be >= 0")
        let p = T.one / (T.one + T.exp(-logitP))
        var s: T = .zero
        for x in data {
            s += x - n * p
        }
        return (s)
    }
    
    /// Convenience overload forwarding integer counts to the generic binomial score.
    public func score(x: Int, n: Int, p: T) -> T {
        return score(x: T(x), n: T(n), p: p)
    }
    /// Sum of integer-valued scores mapped into the floating-point generic implementation.
    public func totalScore(data: [Int], n: Int, p: T) -> T {
        return totalScore(data: data.map(T.init), n: T(n), p: p)
    }
    /// Integer overload for the logit-parameterized binomial score.
    public func scoreWithLogitP(x: Int, n: Int, logitP: T) -> T {
        return scoreWithLogitP(x: T(x), n: T(n), logitP: logitP)
    }
    /// Integer overload for the aggregated logit-parameterized binomial score.
    public func totalScoreWithLogitP(data: [Int], n: Int, logitP: T) -> T {
        return totalScoreWithLogitP(data: data.map(T.init), n: T(n), logitP: logitP)
    }
}
