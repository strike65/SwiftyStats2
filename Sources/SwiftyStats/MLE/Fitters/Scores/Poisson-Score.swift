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

/// Score helpers for the Poisson distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Poisson distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Poisson {
    /// Score for Poisson(x | λ), x ∈ {0,1,2,...}, λ > 0.
    ///
    /// - Log-likelihood: ℓ = x log λ − λ − log x!.
    /// - Score: ∂ℓ/∂λ = x/λ − 1.
    public func score(x: T, lambda: T) -> T {
        precondition(lambda > 0, "lambda>0")
        return ( x / lambda - T.one )
    }

    /// Sum of Poisson scores across data.
    public func totalScore(data: [T], lambda: T) -> T {
        var s = T.zero
        for x in data { s += score(x: x, lambda: lambda) }
        return (s)
    }

    /// Score with log-parameterization λ = exp(logLambda).
    ///
    /// - Chain rule: ∂/∂log λ = λ ∂/∂λ ⇒ x − λ.
    public func scoreWithLogLambda(x: T, logLambda: T) -> T {
        let lambda = T.exp(logLambda)
        return ( x - lambda )
    }

    /// Total score with log-parameterization across data.
    public func totalScoreWithLogLambda(data: [T], logLambda: T) -> T {
        let lambda = T.exp(logLambda)
        var s = T.zero
        for x in data { s += (x - lambda) }
        return (s)
    }
}
