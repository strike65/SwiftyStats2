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

/// Score helpers for the GammaInverse distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the InverseGamma distribution used by MLEFitter.
extension SwiftyBoost.Distribution.InverseGamma {

    /// Score for Inverse-Gamma(x | α, β), α > 0, β > 0, x > 0.
    ///
    /// - Density: f(x) ∝ β^α x^{−α−1} exp(−β/x).
    /// - Log-likelihood: ℓ = α log β − log Γ(α) − (α + 1) log x − β/x.
    /// - Score:
    ///   - ∂ℓ/∂α = log β − ψ(α) − log x
    ///   - ∂ℓ/∂β = α/β − 1/x
    public func score(x: T, alpha: T, beta: T) -> (dalpha: T, dbeta: T) {
        precondition(x > 0 && alpha > 0 && beta > 0, "x>0, alpha>0, beta>0")

        // dalpha parts
        let p1 = T.log(beta)
        let p2 = -Distribution.dg(alpha)
        let p3 = -T.log(x)
        let dalpha = (p1 + p2) + p3

        // dbeta parts
        let q1 = alpha / beta
        let q2 = -T.one / x
        let dbeta = q1 + q2

        return (dalpha, dbeta)
    }

    /// Sum of Inverse-Gamma scores across data.
    public func totalScore(data: [T], alpha: T, beta: T) -> (dalpha: T, dbeta: T) {
        precondition(alpha > 0 && beta > 0, "alpha>0, beta>0")

        let logBeta = T.log(beta)
        let psiA = Distribution.dg(alpha)

        var SA = T.zero, SB = T.zero
        for x in data {
            precondition(x > 0, "x>0")

            // dalpha per x
            let p1 = logBeta
            let p2 = -psiA
            let p3 = -T.log(x)
            let da = (p1 + p2) + p3

            // dbeta per x
            let q1 = alpha / beta
            let q2 = -T.one / x
            let db = q1 + q2

            SA += da
            SB += db
        }
        return (SA, SB)
    }

    // ---------- Log-reparameterizations ----------
    /// Score with log-shape α = exp(logAlpha).
    ///
    /// - Chain rule: ∂/∂log α = α ∂/∂α.
    public func scoreWithLogAlpha(x: T, logAlpha: T, beta: T) -> (dlogAlpha: T, dbeta: T) {
        let alpha = T.exp(logAlpha)
        let (da, db) = score(x: x, alpha: alpha, beta: beta)
        let r1 = alpha * da
        return (r1, db)
    }

    /// Total score with log-shape across data.
    public func totalScoreWithLogAlpha(data: [T], logAlpha: T, beta: T) -> (dlogAlpha: T, dbeta: T) {
        let alpha = T.exp(logAlpha)
        precondition(beta > 0 && alpha > 0, "alpha>0, beta>0")
        do {
            let logBeta = T.log(beta)
            let psiA = try SwiftyBoost.SpecialFunctions.digamma(alpha)
            var SA = T.zero, SB = T.zero
            for x in data {
                precondition(x > 0, "x>0")
                // d/d log α = α * (∂ℓ/∂α)
                SA += alpha * (logBeta - psiA - T.log(x))
                // ∂ℓ/∂β (unchanged)
                SB += alpha / beta - T.one / x
            }
            return (SA, SB)
        }
        catch _ {
            return (T.nan, T.nan)
        }
    }

    /// Score with log-scale β = exp(logBeta).
    ///
    /// - Chain rule: ∂/∂log β = β ∂/∂β.
    public func scoreWithLogBeta(x: T, alpha: T, logBeta: T) -> (dalpha: T, dlogBeta: T) {
        let beta = T.exp(logBeta)
        let (da, db) = score(x: x, alpha: alpha, beta: beta)
        let r2 = beta * db
        return (da, r2)
    }

    /// Total score with log-scale across data.
    public func totalScoreWithLogBeta(data: [T], alpha: T, logBeta: T) -> (dalpha: T, dlogBeta: T) {
        let beta = T.exp(logBeta)
        precondition(alpha > 0 && beta > 0, "alpha>0, beta>0")
        do {
            let logBeta = T.log(beta) // equals input logBeta; kept for symmetry/readability
            let psiA = try SwiftyBoost.SpecialFunctions.digamma(alpha)
            var SA = T.zero, SB = T.zero
            for x in data {
                precondition(x > 0, "x>0")
                // ∂ℓ/∂α
                SA += logBeta - psiA - T.log(x)
                // d/d log β = β * (∂ℓ/∂β) = α − β/x
                SB += alpha - beta / x
            }
            return (SA, SB)
        }
        catch _ {
            return (T.nan, T.nan)
        }
    }

    /// Score with both logs: α = exp(logAlpha), β = exp(logBeta).
    ///
    /// - Chain rule: ∂/∂log α = α ∂/∂α, ∂/∂log β = β ∂/∂β.
    public func scoreWithLogParams(x: T, logAlpha: T, logBeta: T) -> (dlogAlpha: T, dlogBeta: T) {
        let alpha = T.exp(logAlpha)
        let beta  = T.exp(logBeta)
        let (da, db) = score(x: x, alpha: alpha, beta: beta)
        let r1 = alpha * da
        let r2 = beta * db
        return (r1, r2)
    }

    /// Total score with both logs across data.
    public func totalScoreWithLogParams(data: [T], logAlpha: T, logBeta: T) -> (dlogAlpha: T, dlogBeta: T) {
        let alpha = T.exp(logAlpha)
        let beta  = T.exp(logBeta)
        precondition(alpha > 0 && beta > 0, "alpha>0, beta>0")
        do {
            let logBeta = T.log(beta) // equals input logBeta; kept for symmetry/readability
            let psiA = try SwiftyBoost.SpecialFunctions.digamma(alpha)
            var SA = T.zero, SB = T.zero
            for x in data {
                precondition(x > 0, "x>0")
                // d/d log α
                SA += alpha * (logBeta - psiA - T.log(x))
                // d/d log β
                SB += alpha - beta / x
            }
            return (SA, SB)
        }
        catch _ {
            return (T.nan, T.nan)
        }
    }
}
