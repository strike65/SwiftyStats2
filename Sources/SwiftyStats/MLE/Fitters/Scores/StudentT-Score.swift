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

/// Score helpers for the StudentT distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the StudentT distribution used by MLEFitter.
extension SwiftyBoost.Distribution.StudentT {
    /// Score for standardized Student-t(x | ν) with μ = 0, σ = 1, ν > 0, w.r.t. ν.
    ///
    /// - Log-likelihood:
    ///   ℓ = log Γ((ν+1)/2) − log Γ(ν/2) − 0.5 log(π ν) − ((ν+1)/2) log(1 + x²/ν).
    /// - Score:
    ///   ∂ℓ/∂ν = 0.5(ψ((ν+1)/2) − ψ(ν/2)) − 0.5/ν − 0.5 log(1 + x²/ν) + ((ν+1)/2) x² / (ν(ν + x²))
    /// - Constraints: ν > 0.
    @inline(__always)
    public func score(x: T, nu: T) -> T {
        precondition(nu > 0, "nu must be > 0")
        let a = (nu + 1) * T.half        // (nu+1)/2
        let b = nu * T.half            // nu/2
        let psiDiff = Distribution.dg(a) - Distribution.dg(b)
        let z2 = x * x
        // d/dν log(1 + x^2/ν) = - x^2 / (ν(ν + x^2))
        let logTerm = T.log(1 + z2/nu)
        let frac = z2 / (nu * (nu + z2))
        // ∂ℓ/∂ν
        let dnu = T.half * psiDiff - T.half / nu - T.half * logTerm + a * frac
        return dnu
    }
    
    /// Sum of standardized Student-t scores across data (w.r.t. ν).
    @inline(__always)
    public func totalScore(data: [T], nu: T) -> T {
        precondition(nu > 0, "nu must be > 0")
        let a = (nu + 1) * T.half
        let b = nu * T.half
        let psiDiff = Distribution.dg(a) - Distribution.dg(b)
        var s = T.zero
        for x in data {
            let z2 = x * x
            let logTerm = T.log(1 + z2/nu)
            let frac = z2 / (nu * (nu + z2))
            s += T.half * psiDiff - T.half / nu - T.half * logTerm + a * frac
        }
        return s
    }
    
    /// Score with log-parameterization ν = exp(logNu).
    ///
    /// - Chain rule: ∂/∂log ν = ν ∂/∂ν.
    @inline(__always)
    public func scoreWithLogNu(x: T, logNu: T) -> T {
        let nu = T.exp(logNu)
        let (dnu) = score(x: x, nu: nu)
        return nu * dnu                 // chain rule
    }
    
    /// Total score with log-parameterization across data.
    @inline(__always)
    public func totalScoreWithLogNu(data: [T], logNu: T) -> T {
        let nu = T.exp(logNu)
        let (s) = totalScore(data: data, nu: nu)
        return nu * s
    }
}
