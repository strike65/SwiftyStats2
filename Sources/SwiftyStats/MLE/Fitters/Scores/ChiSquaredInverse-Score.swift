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

/// Score helpers for the ChiSquaredInverse distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the InverseChiSquared distribution used by MLEFitter.
extension SwiftyBoost.Distribution.InverseChiSquared {
    // The scaled inverse-chi-squared with df nu > 0 and scale s2 > 0:
    // X ~ InvChiSq(nu, s2)  ⇔  X ~ InvGamma(alpha = nu/2, beta = nu*s2/2)
    //
    // PDF:
    //   f(x) = [beta^alpha / Gamma(alpha)] * x^(-alpha-1) * exp(-beta / x),  x > 0
    // Log-likelihood for one x:
    //   ell = alpha*log(beta) - lgamma(alpha) - (alpha+1)*log(x) - beta/x
    //
    // Scores (derivatives of ell):
    //   alpha = nu/2, beta = nu*s2/2
    //   d ell / d alpha = log(beta) - psi(alpha) - log(x)
    //   d ell / d beta  = alpha/beta - 1/x  = 1/s2 - 1/x
    //   d ell / d nu    = (1/2)*[log(beta) - psi(alpha) - log(x)] + (s2/2)*(1/s2 - 1/x)
    //                    = (1/2)*[log(beta) - psi(alpha) - log(x) + 1 - s2/x]
    //   d ell / d s2    = (nu/2)*(1/s2 - 1/x)
    //
    /// Score for a single observation.
    public func score(x: T, nu: T, s2: T) -> (dNu: T, dS2: T) {
        precondition(x > .zero && x.isFinite, "x must be > 0 and finite")
        precondition(nu > .zero && nu.isFinite, "nu must be > 0 and finite")
        precondition(s2 > .zero && s2.isFinite, "s2 must be > 0 and finite")

        let alpha: T = nu * T.half
        let beta:  T = nu * s2 * T.half

        // core pieces
        let logBeta: T = T.log(beta)
        let psiAlpha: T = Distribution.dg(alpha) // digamma(alpha)
        let logx: T = T.log(x)

        // d ell / d nu
        let dNu: T = T.half * (logBeta - psiAlpha - logx + T.one - (s2 / x))

        // d ell / d s2
        let dS2: T = (nu * T.half) * (T.one / s2 - T.one / x)

        return (dNu, dS2)
    }

    /// Sum of scores across the dataset.
    public func totalScore(data: [T], nu: T, s2: T) -> (dNu: T, dS2: T) {
        precondition(!data.isEmpty, "data must be non-empty")
        precondition(nu > .zero && s2 > .zero, "nu, s2 must be > 0")

        var sNu: T = .zero
        var sS2: T = .zero
        for xi in data {
            let (a, b) = score(x: xi, nu: nu, s2: s2)
            sNu += a
            sS2 += b
        }
        return (sNu, sS2)
    }

    // MARK: - Log-parameterization (nu = exp(logNu), s2 = exp(logS2))

    /// Score with log-parameters to remove positivity constraints.
    ///
    /// Chain rule:
    ///   d ell / d logNu = nu * (d ell / d nu)
    ///   d ell / d logS2 = s2 * (d ell / d s2) = (nu/2) * (1 - s2/x)
    public func scoreWithLogParams(x: T, logNu: T, logS2: T) -> (dLogNu: T, dLogS2: T) {
        let nu = T.exp(logNu)
        let s2 = T.exp(logS2)
        let (dNu, _) = score(x: x, nu: nu, s2: s2)
        let dLogNu: T = nu * dNu
        // Use simplified closed form for d/d logS2
        let dLogS2: T = (nu * T.half) * (T.one - s2 / x)
        return (dLogNu, dLogS2)
    }

    /// Sum of scores with log-parameters across the dataset.
    public func totalScoreWithLogParams(data: [T], logNu: T, logS2: T) -> (dLogNu: T, dLogS2: T) {
        let nu = T.exp(logNu)
        let s2 = T.exp(logS2)

        var sLnNu: T = .zero
        var sLnS2: T = .zero
        for xi in data {
            // d/d logNu = nu * d/d nu; d/d logS2 = (nu/2)*(1 - s2/x)
            let (dNu, _) = score(x: xi, nu: nu, s2: s2)
            sLnNu += nu * dNu
            sLnS2 += (nu * T.half) * (T.one - s2 / xi)
        }
        return (sLnNu, sLnS2)
    }
}
