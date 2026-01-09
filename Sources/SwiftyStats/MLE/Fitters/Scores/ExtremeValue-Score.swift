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

/// Score helpers for the ExtremeValue distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude


/// Score helpers for the ExtremeValueGumbel distribution used by MLEFitter.
extension SwiftyBoost.Distribution.ExtremeValueGumbel {
    /// Score for Gumbel/Extreme-Value Type I: X ~ Gumbel(μ, β), β > 0.
    ///
    /// - Density: f(x) = (1/β) exp(−z − e^{−z}), z = (x − μ)/β.
    /// - Log-likelihood: ℓ = −log β − z − e^{−z}.
    /// - Score:
    ///   - ∂ℓ/∂μ = (1 − e^{−z}) / β
    ///   - ∂ℓ/∂β = (−1 + z − z e^{−z}) / β
    /// - Constraints: β > 0.
    public func score(x: T, mu: T, beta: T) -> (dmu: T, dbeta: T) {
        precondition(beta > 0, "beta>0")
        let z = (x - mu) / beta
        let ez = T.exp(-z)
        let dmu   = (T.one - ez) / beta
        let dbeta = (-T.one + z - z * ez) / beta
        return (dmu, dbeta)
    }

    /// Sum of Gumbel scores across data.
    public func totalScore(data: [T], mu: T, beta: T) -> (dmu: T, dbeta: T) {
        var a = T.zero, b = T.zero
        for x in data {
            let (u,v) = score(x: x, mu: mu, beta: beta)
            a += u; b += v
        }
        return (a,b)
    }

    /// Score with log-scale β = exp(logBeta).
    ///
    /// - Chain rule: ∂/∂log β = β ∂/∂β.
    public func scoreWithLogBeta(x: T, mu: T, logBeta: T) -> (dmu: T, dlogBeta: T) {
        let beta = T.exp(logBeta)
        let (u,v) = score(x: x, mu: mu, beta: beta)
        return (u, beta * v)
    }

    /// Total score with log-scale across data.
    public func totalScoreWithLogBeta(data: [T], mu: T, logBeta: T) -> (dmu: T, dlogBeta: T) {
        let beta = T.exp(logBeta)
        var a = T.zero, b = T.zero
        for x in data {
            let (u,v) = score(x: x, mu: mu, beta: beta)
            a += u; b += beta * v
        }
        return (a,b)
    }
}
