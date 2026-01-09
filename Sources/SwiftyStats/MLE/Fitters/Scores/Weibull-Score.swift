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

/// Score helpers for the Weibull distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Weibull distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Weibull {
    /// Score for Weibull(x | k, λ), k > 0 (shape), λ > 0 (scale), x ≥ 0.
    ///
    /// - Density: f(x) = (k/λ) (x/λ)^{k−1} exp(−(x/λ)^k).
    /// - Log-likelihood: ℓ = log k − log λ + (k − 1) log x − k log λ − (x/λ)^k.
    /// - Score:
    ///   - ∂ℓ/∂k = 1/k + log(x/λ) − (x/λ)^k log(x/λ)
    ///   - ∂ℓ/∂λ = (k/λ)((x/λ)^k − 1)
    public func score(x: T, k: T, lambda: T) -> (dk: T, dlambda: T) {
        precondition(x >= 0 && k > 0 && lambda > 0, "x>=0, k>0, lambda>0")
        let logx = T.log(x)
        let logl = T.log(lambda)
        let logxOverL = logx - logl                      // log(x/lambda); -inf wenn x=0
        let t = T.exp(k * logxOverL)                     // (x/lambda)^k  (0^k = 0 bei x=0)
        let dk = T.one / k + logxOverL - t * logxOverL
        let dl = (k / lambda) * (t - T.one)
        return (dk, dl)
    }
    
    /// Sum of Weibull scores across data.
    public func totalScore(data: [T], k: T, lambda: T) -> (dk: T, dlambda: T) {
        var a = T.zero, b = T.zero
        for x in data {
            let (u,v) = score(x: x, k: k, lambda: lambda)
            a += u; b += v
        }
        return (a,b)
    }
    
    /// Score with log-parameters: k = exp(logK), λ = exp(logLambda).
    ///
    /// - Chain rule: ∂/∂log k = k ∂/∂k, ∂/∂log λ = λ ∂/∂λ.
    public func scoreWithLogParams(x: T, logK: T, logLambda: T) -> (dlogK: T, dlogLambda: T) {
        let k = T.exp(logK), lambda = T.exp(logLambda)
        let (dk, dl) = score(x: x, k: k, lambda: lambda)
        return (k * dk, lambda * dl)
    }
    
    /// Total score with log-parameters across data.
    public func totalScoreWithLogParams(data: [T], logK: T, logLambda: T) -> (dlogK: T, dlogLambda: T) {
        let k = T.exp(logK), lambda = T.exp(logLambda)
        var a = T.zero, b = T.zero
        for x in data {
            let (u,v) = score(x: x, k: k, lambda: lambda)
            a += k * u; b += lambda * v
        }
        return (a,b)
    }
}
