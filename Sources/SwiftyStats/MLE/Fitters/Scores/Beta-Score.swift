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
/// Score helpers for the Beta distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Beta distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Beta {
    /// Score for Beta(x | a, b), x ∈ (0,1), a > 0, b > 0.
    ///
    /// - Log-likelihood: ℓ = (a − 1) log x + (b − 1) log(1 − x) − log B(a, b).
    /// - Score:
    ///   - ∂ℓ/∂a = log x − (ψ(a) − ψ(a + b))
    ///   - ∂ℓ/∂b = log(1 − x) − (ψ(b) − ψ(a + b))
    /// - Constraints: x ∈ (0,1), a > 0, b > 0.
    /// - Numerical functions: uses digamma ψ.
    public func score(x: T, a: T, b: T) -> (da: T, db: T) {
        precondition(x > 0 && x < T.one && a > 0 && b > 0, "x in (0,1), a>0, b>0")
        do {
            let psiA = try SwiftyBoost.SpecialFunctions.digamma(a)
            let psiB = try SwiftyBoost.SpecialFunctions.digamma(b)
            let psiAB = try SwiftyBoost.SpecialFunctions.digamma(a + b)
            let da = T.log(x) - (psiA - psiAB)
            let db = T.log(T.one - x) - (psiB - psiAB)
            return (da, db)
        }
        catch _ {
            return (T.nan, T.nan)
        }
    }
    
    /// Sum of Beta scores across a dataset.
    public func totalScore(data: [T], a: T, b: T) -> (da: T, db: T) {
        do {
            var A = T.zero, B = T.zero
            let psiA = try SwiftyBoost.SpecialFunctions.digamma(a)
            let psiB = try SwiftyBoost.SpecialFunctions.digamma(b)
            let psiAB = try SwiftyBoost.SpecialFunctions.digamma(a + b)
            for x in data {
                precondition(x > 0 && x < T.one, "x in (0,1)")
                A += T.log(x) - (psiA - psiAB)
                B += T.log(T.one - x) - (psiB - psiAB)
            }
            return (A,B)
        }
        catch _ {
            return (T.nan, T.nan)
        }
    }
    
    /// Score with log-parameters a = exp(logA), b = exp(logB).
    ///
    /// - Chain rule: ∂/∂log a = a ∂/∂a, ∂/∂log b = b ∂/∂b.
    public func scoreWithLogParams(x: T, logA: T, logB: T) -> (dlogA: T, dlogB: T) {
        let a = T.exp(logA), b = T.exp(logB)
        let (da, db) = score(x: x, a: a, b: b)
        return (a * da, b * db)
    }
    
    /// Sum of scores with log-parameters across the dataset.
    public func totalScoreWithLogParams(data: [T], logA: T, logB: T) -> (dlogA: T, dlogB: T) {
        do {
            let a = T.exp(logA), b = T.exp(logB)
            var A = T.zero, B = T.zero
            let psiA = try SwiftyBoost.SpecialFunctions.digamma(a)
            let psiB = try SwiftyBoost.SpecialFunctions.digamma(b)
            let psiAB = try SwiftyBoost.SpecialFunctions.digamma(a + b)
            for x in data {
                precondition(x > 0 && x < T.one, "x in (0,1)")
                A += a * (T.log(x) - (psiA - psiAB))
                B += b * (T.log(T.one - x) - (psiB - psiAB))
            }
            return (A,B)
        }
        catch _ {
            return (T.nan, T.nan)
        }
    }
}
