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

/// Score helpers for the Uniform distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Uniform distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Uniform {
    /// Score for Uniform(x | a, b) with a < x < b (open interval for differentiability).
    ///
    /// - Density: f(x) = 1/(b − a) on [a, b].
    /// - Log-likelihood: ℓ = −log(b − a), for x ∈ [a, b].
    /// - Score (interior, a < x < b):
    ///   - ∂ℓ/∂a =  1/(b − a)
    ///   - ∂ℓ/∂b = −1/(b − a)
    /// - Constraints: a < b, and we require a < x < b to avoid boundary issues.
    @inline(__always)
    public func score(x: T, a: T, b: T) -> (da: T, db: T) {
        precondition(a < b && a < x && x < b, "a<b and inside open interval")
        let invW = T.one / (b - a)
        return ( invW, -invW )
    }

    /// Sum of Uniform scores across data (interior points).
    @inline(__always)
    public func totalScore(data: [T], a: T, b: T) -> (da: T, db: T) {
        var A = T.zero, B = T.zero
        for x in data {
            let (u,v) = score(x: x, a: a, b: b)
            A += u; B += v
        }
        return (A,B)
    }
}
