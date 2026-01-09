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

/// Score helpers for the Gamma distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Gamma distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Gamma {
    /// Score for Gamma(x | k, θ), k > 0 (shape), θ > 0 (scale).
    ///
    /// - Density: f(x) ∝ x^{k−1} e^{−x/θ} / (θ^k Γ(k)), x > 0.
    /// - Log-likelihood: ℓ = (k − 1) log x − x/θ − k log θ − log Γ(k).
    /// - Score:
    ///   - ∂ℓ/∂k = log x − log θ − ψ(k)
    ///   - ∂ℓ/∂θ = x/θ² − k/θ
    /// - Constraints: x > 0, k > 0, θ > 0.
    public func score(x: T, k: T, theta: T) -> (dk: T, dtheta: T) {
        precondition(x > 0 && k > 0 && theta > 0, "x>0, k>0, theta>0")
        do {
            let dg = try SwiftyBoost.SpecialFunctions.digamma(k)
            let dk = T.log(x) - T.log(theta) - dg
            let dtheta  = x/(theta*theta) - k/theta
            return (dk, dtheta)
        }
        catch _ {
            return (T.nan, T.nan)
        }
    }
    
    /// Sum of Gamma scores across data.
    public func totalScore(data: [T], k: T, theta: T) -> (dk: T, dtheta: T) {
        precondition(k > 0 && theta > 0, "k>0, theta>0")
        do {
            let psiK = try SwiftyBoost.SpecialFunctions.digamma(k)
            let logTheta = T.log(theta)
            var sK = T.zero, sTh = T.zero
            for x in data {
                precondition(x > 0, "x>0")
                sK  += T.log(x) - logTheta - psiK
                sTh += x/(theta*theta) - k/theta
            }
            return (sK, sTh)
        }
        catch _ {
            return (T.nan, T.nan)
        }
    }
    
    /// Score with log-scale: θ = exp(logTheta).
    ///
    /// - Chain rule: ∂/∂log θ = θ ∂/∂θ.
    public func scoreWithLogTheta(x: T, k: T, logTheta: T) -> (dk: T, dlogTheta: T) {
        let theta = T.exp(logTheta)
        let (dk, dth) = score(x: x, k: k, theta: theta)
        return (dk, theta * dth) // chain rule
    }
    
    /// Total score with log-scale across data.
    public func totalScoreWithLogTheta(data: [T], k: T, logTheta: T) -> (dk: T, dlogTheta: T) {
        do {
            let theta = T.exp(logTheta)
            let psiK = try SwiftyBoost.SpecialFunctions.digamma(k)
            let logThetaVal = T.log(theta) // equals logTheta, kept for symmetry
            var sK = T.zero, sL = T.zero
            for x in data {
                precondition(x > 0, "x>0")
                sK += T.log(x) - logThetaVal - psiK
                sL += theta * (x/(theta*theta) - k/theta) // = x/theta - k
            }
            return (sK, sL)
        }
        catch _ {
            return (T.nan, T.nan)
        }
    }
    
    /// Score with log-shape: k = exp(logK).
    ///
    /// - Chain rule: ∂/∂log k = k ∂/∂k.
    public func scoreWithLogK(x: T, logK: T, theta: T) -> (dlogK: T, dtheta: T) {
        let k = T.exp(logK)
        let (dk, dth) = score(x: x, k: k, theta: theta)
        return (k * dk, dth)
    }
    
    /// Total score with log-shape across data.
    public func totalScoreWithLogK(data: [T], logK: T, theta: T) -> (dlogK: T, dtheta: T) {
        do {
            let k = T.exp(logK)
            let psiK = try SwiftyBoost.SpecialFunctions.digamma(k)
            let logTheta = T.log(theta)
            var sLK = T.zero, sTh = T.zero
            for x in data {
                precondition(x > 0, "x>0")
                sLK += k * (T.log(x) - logTheta - psiK)
                sTh += x/(theta*theta) - k/theta
            }
            return (sLK, sTh)
        }
        catch _ {
            return (T.nan, T.nan)
        }
    }
    
    /// Score with both logs: k = exp(logK), θ = exp(logTheta).
    ///
    /// - Chain rule: ∂/∂log k = k ∂/∂k, ∂/∂log θ = θ ∂/∂θ.
    public func scoreWithLogParams(x: T, logK: T, logTheta: T) -> (dlogK: T, dlogTheta: T) {
        let k = T.exp(logK), theta = T.exp(logTheta)
        let (dk, dth) = score(x: x, k: k, theta: theta)
        return (k * dk, theta * dth)
    }
    
    /// Total score with both logs across data.
    public func totalScoreWithLogParams(data: [T], logK: T, logTheta: T) -> (dlogK: T, dlogTheta: T) {
        do {
            let k = T.exp(logK), theta = T.exp(logTheta)
            let psiK = try SwiftyBoost.SpecialFunctions.digamma(k)
            let logThetaVal = T.log(theta)
            var sLK = T.zero, sLTh = T.zero
            for x in data {
                precondition(x > 0, "x>0")
                sLK  += k * (T.log(x) - logThetaVal - psiK)
                sLTh += theta * (x/(theta*theta) - k/theta) // = x/theta - k
            }
            return (sLK, sLTh)
        }
        catch _ {
            return (T.nan, T.nan)
        }
    }
}
