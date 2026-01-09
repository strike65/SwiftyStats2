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

/// Score helpers for the ChiSquared distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude


/// Score helpers for the ChiSquared distribution used by MLEFitter.
extension SwiftyBoost.Distribution.ChiSquared {

    /// Score for Chi-square(x | k), k > 0, x > 0, w.r.t. k.
    ///
    /// - Log-likelihood: ℓ = (k/2 − 1) log x − x/2 − (k/2) log 2 − log Γ(k/2).
    /// - Score:
    ///   - ∂ℓ/∂k = 0.5 log x − 0.5 log 2 − 0.5 ψ(k/2)
    /// - Constraints: x > 0, k > 0.
    public func score(x: T, k: T) -> T {
        precondition(x > 0 && k > 0, "x>0, k>0")
        let h = T.half
        let a = k * h                  // k/2
        let term1 = h * T.log(x)
        let term2 = -h * T.lnTwo
        let term3 = -h * Distribution.dg(a)
        let dk = (term1 + term2) + term3
        return dk
    }

    /// Sum of Chi-square scores across the dataset (w.r.t. k).
    public func totalScore(data: [T], k: T) -> T {
        precondition(k > 0, "k>0")
        let h = T.half
        let a = k * h
        let c = -h * ( T.log(2) + Distribution.dg(a) )   // constant across x
        var s = T.zero
        for x in data {
            precondition(x > 0, "x>0")
            let t = h * T.log(x)
            s += t + c
        }
        return s
    }

    /// Score with log-parameterization k = exp(logK).
    ///
    /// - Chain rule: ∂/∂log k = k ∂/∂k.
    public func scoreWithLogK(x: T, logK: T) -> T {
        let k = T.exp(logK)
        let (dk) = score(x: x, k: k)
        let dlogK = k * dk
        return dlogK
    }

    /// Total score with log-parameterization across data.
    public func totalScoreWithLogK(data: [T], logK: T) -> T {
        let k = T.exp(logK)
        let (s) = totalScore(data: data, k: k)
        let dlogK = k * s
        return dlogK
    }
}
