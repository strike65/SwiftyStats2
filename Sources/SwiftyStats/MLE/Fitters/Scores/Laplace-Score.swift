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

/// Score helpers for the Laplace distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Laplace distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Laplace {
    @inline(__always) private func sgn(_ z: T) -> T { z > 0 ? T.one : (z < 0 ? -T.one : 0) }

    /// Score for Laplace(x | μ, b), b > 0.
    ///
    /// - Density: f(x) = (1/(2b)) exp(−|x − μ|/b).
    /// - Log-likelihood: ℓ = −log(2b) − |x − μ|/b.
    /// - Score:
    ///   - ∂ℓ/∂μ = sign(x − μ)/b
    ///   - ∂ℓ/∂b = −1/b + |x − μ|/b²
    public func score(x: T, mu: T, b: T) -> (dmu: T, db: T) {
        precondition(b > 0, "b must be > 0")
        let dmu: T = sgn(x - mu) / b
        let db:T  = -T.one / b + (x - mu).magnitude / (b * b)
        return (dmu, db)
    }

    /// Sum of Laplace scores across data.
    public func totalScore(data: [T], mu: T, b: T) -> (dmu: T, db: T) {
        var a = T.zero, c = T.zero
        for x in data {
            let (u,v) = score(x: x, mu: mu, b: b)
            a += u; c += v
        }
        return (a,c)
    }

    /// Score with log-scale b = exp(logB).
    ///
    /// - Chain rule: ∂/∂log b = b ∂/∂b.
    public func scoreWithLogB(x: T, mu: T, logB: T) -> (dmu: T, dlogB: T) {
        let b = T.exp(logB)
        let (u,v) = score(x: x, mu: mu, b: b)
        return (u, b * v)
    }

    /// Total score with log-scale across data.
    public func totalScoreWithLogB(data: [T], mu: T, logB: T) -> (dmu: T, dlogB: T) {
        var a = T.zero, c = T.zero
        let b = T.exp(logB)
        for x in data {
            let (u,v) = score(x: x, mu: mu, b: b)
            a += u; c += b * v
        }
        return (a,c)
    }
}

