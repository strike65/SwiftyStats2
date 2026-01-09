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

/// Score helpers for the Pareto distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Pareto distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Pareto {
    /// Score for Pareto(x | scale, shape), x ≥ scale > 0, shape > 0.
    ///
    /// - Density: f(x) = shape scale^shape / x^{shape+1}.
    /// - Log-likelihood: ℓ = log shape + shape log scale − (α + 1) log x.
    /// - Score:
    ///   - ∂ℓ/∂scale = shape / scale
    ///   - ∂ℓ/∂shape = 1/shape + log scale − log x
    public func score(x: T, scale: T, shape: T) -> (dxm: T, dalpha: T) {
        precondition(scale > 0 && shape > 0 && x >= scale, "xm>0, alpha>0, x>=xm")
        let dxm = shape / scale
        let dal = T.one / shape + (T.log(scale) - T.log(x))
        return (dxm, dal)
    }
    
    /// Sum of Pareto scores across data.
    public func totalScore(data: [T], scale: T, shape: T) -> (dxm: T, dalpha: T) {
        var a = T.zero, b = T.zero
        for x in data {
            let (u,v) = score(x: x, scale: scale, shape: shape)
            a += u; b += v
        }
        return (a,b)
    }
    
    /// Score with log-parameters x_m = exp(logXm), α = exp(logAlpha).
    ///
    /// - Chain rule: ∂/∂log x_m = x_m ∂/∂x_m, ∂/∂log α = α ∂/∂α.
    public func scoreWithLogParams(x: T, logScale: T, logShape: T) -> (dlogXm: T, dlogAlpha: T) {
        let xm = T.exp(logScale), alpha = T.exp(logShape)
        let (u,v) = score(x: x, scale: xm, shape: alpha)
        return (xm * u, alpha * v)
    }
    
    /// Total score with log-parameters across data.
    public func totalScoreWithLogParams(data: [T], logScale: T, logShape: T) -> (dlogXm: T, dlogAlpha: T) {
        let xm = T.exp(logScale), alpha = T.exp(logShape)
        var a = T.zero, b = T.zero
        for x in data {
            let (u,v) = score(x: x, scale: xm, shape: alpha)
            a += xm * u; b += alpha * v
        }
        return (a,b)
    }
}
