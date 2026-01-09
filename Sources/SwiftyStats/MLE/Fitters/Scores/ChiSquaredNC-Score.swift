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

/// Score helpers for the ChiSquaredNC distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

// ====================================
// NONCENTRAL CHI-SQUARED (df>0, lambda>0), x>0
//
// pdf: f(x) = 0.5 * exp(-(x + lambda)/2)
//             * (x/lambda)^(df/4 - 1/2)
//             * I_{nu}(sqrt(lambda * x)),
// where nu = df/2 - 1.
//
// log f(x) = -log 2
//            - (x + lambda)/2
//            + (df/4 - 1/2) * (log x - log lambda)
//            + log I_{nu}(z), z = sqrt(lambda * x)
//
// Score w.r.t. (df, lambda):
//
// nu = df/2 - 1
// a  = df/4 - 1/2
// z  = sqrt(lambda * x)
//
// ∂ℓ/∂df
//   = 0.25 * (log x - log lambda)
//     + 0.5 * ∂/∂nu [ log I_nu(z) ]
//
// ∂ℓ/∂lambda
//   = -1/2
//     - a / lambda
//     + 0.5 * (z / lambda) * ∂/∂z [ log I_nu(z) ]
// ====================================
/// Score helpers for the NonCentralChiSquared distribution used by MLEFitter.
extension SwiftyBoost.Distribution.NonCentralChiSquared {
    
    @inline(__always)
    /// Single-observation non-central χ² score with respect to the degrees of freedom and non-centrality.
    public func score(x: T, df: T, lambda: T) -> (ddf: T, dlambda: T) {
        precondition(x > .zero && df > .zero && lambda > .zero, "x>0, df>0, lambda>0")
        
        let halfDf = df * T.half         // df/2
        let nu = halfDf - T.one       // df/2 - 1
        let a = df * T.quarter - T.half   // df/4 - 1/2
        
        let logX = T.log(x)
        let logLambda = T.log(lambda)
        let logRatio = logX - logLambda
        
        let z = T.sqrt(lambda * x)
        
        // Placeholders for special-function derivatives:
        // dLogI_dNu = ∂/∂nu log I_nu(z)
        // dLogI_dZ  = ∂/∂z  log I_nu(z)
        let dLogI_dNu = besselILogDerivativeWrtOrder(nu, z)
        let dLogI_dZ  = besselILogDerivativeWrtArgument(nu, z)
        
        // ∂ℓ/∂df
        let ddf = T.quarter * logRatio + T.half * dLogI_dNu
        
        // ∂ℓ/∂lambda
        let base = -T.half - a / lambda
        let dlambda = base + T.half * (z / lambda) * dLogI_dZ
        
        return (ddf, dlambda)
    }
    
    @inline(__always)
    /// Accumulates the non-central χ² scores over an entire dataset.
    public func totalScore(data: [T], df: T, lambda: T) -> (ddf: T, dlambda: T) {
        precondition(df > .zero && lambda > .zero, "df>0, lambda>0")
        
        let halfDf = df * T.half          // df/2
        let nu = halfDf - 1        // df/2 - 1
        let a = df * T.quarter - T.half   // df/4 - 1/2
        let logLambda = T.log(lambda)
        
        var sumDf: T = .zero
        var sumLambda: T = .zero
        
        for x in data {
            precondition(x > 0, "x>0")
            
            let logX = T.log(x)
            let logRatio = logX - logLambda
            
            let z = T.sqrt(lambda * x)
            
            let dLogI_dNu = besselILogDerivativeWrtOrder(nu, z)
            let dLogI_dZ  = besselILogDerivativeWrtArgument(nu, z)
            
            let ddf_i = T.quarter * logRatio + T.half * dLogI_dNu
            
            let base = -T.half - a / lambda
            let dlambda_i = base + T.half * (z / lambda) * dLogI_dZ
            
            sumDf += ddf_i
            sumLambda += dlambda_i
        }
        
        return (sumDf, sumLambda)
    }
    
    /// Score under the log-parameterization where df = exp(logDf) and λ = exp(logLambda).
    @inline(__always)
    public func scoreWithLogParams(
        x: T,
        logDf: T,
        logLambda: T
    ) -> (dlogDf: T, dlogLambda: T) {
        let df     = T.exp(logDf)
        let lambda = T.exp(logLambda)
        let (u, v) = score(x: x, df: df, lambda: lambda)
        return (df * u, lambda * v) // chain rule
    }
    
    /// Total score under the log-parameterization across the dataset.
    @inline(__always)
    public func totalScoreWithLogParams(
        data: [T],
        logDf: T,
        logLambda: T
    ) -> (dlogDf: T, dlogLambda: T) {
        let df = T.exp(logDf)
        let lambda = T.exp(logLambda)
        let (U, V) = totalScore(data: data, df: df, lambda: lambda)
        return (df * U, lambda * V)
    }
    
    /// Computes ∂/∂z log I_ν(z) using the ratio of the derivative and the value of
    /// the modified Bessel function of the first kind.
    ///
    /// The helper wraps SwiftyBoost’s Bessel implementations and returns `NaN` if the
    /// backend reports an error so callers can fall back to numerical guards.
    @inline(__always)
    private func besselILogDerivativeWrtArgument(_ nu: T, _ z: T) -> T {
        do {
            let I: T = try SwiftyBoost.SpecialFunctions.modifiedBesselI(v: nu, x: z)
            let dI_dZ: T = try SwiftyBoost.SpecialFunctions.modifiedBesselIPrime(v: nu, x: z)
            return dI_dZ / I
        }
        catch _ {
            return .nan
        }
    }
    
    /// Finite-difference step for derivatives w.r.t. the Bessel order.
    ///
    /// Uses √ε scaling with the current magnitude to balance truncation and rounding errors.
    @inline(__always)
    private func besselOrderDiffStep(_ order: T) -> T {
        let eps = T.ulpOfOne
        return T.sqrt(eps) * (order.magnitude + T.one)
    }
    
    /// Approximates ∂/∂ν I_ν(z) via a central finite difference in the order parameter.
    ///
    /// - Parameters:
    ///   - order: Bessel order ν.
    ///   - z: Argument of the function.
    /// - Returns: Central difference estimate or `NaN` if the backend fails.
    @inline(__always)
    private func besselIOrderDerivative(_ order: T, _ z: T) -> T {
        let h = besselOrderDiffStep(order)
        // Simple guard to avoid degenerate h; usually not hit in practice.
        precondition(h > 0, "step size for Bessel I order derivative must be positive")
        
        let orderPlus = order + h
        let orderMinus = order - h
        do {
            // Evaluate I at order +/- h
            let IPlus  = try SwiftyBoost.SpecialFunctions.modifiedBesselI(v: orderPlus, x: z)
            let IMinus = try SwiftyBoost.SpecialFunctions.modifiedBesselI(v: orderMinus, x: z)
            // Central difference approximation:
            // dI/dorder(order, z) ≈ (I(order+h, z) - I(order-h, z)) / (2h)
            return (IPlus - IMinus) / (2 * h)
        }
        catch _ {
            return .nan
        }
    }
    
    /// Computes ∂/∂ν log I_ν(z) by dividing the order derivative by I_ν(z).
    ///
    /// Returns `NaN` when the underlying Bessel evaluation fails so upstream callers
    /// can drop the contribution gracefully.
    func besselILogDerivativeWrtOrder(_ order: T, _ z: T) -> T {
        do {
            let I0: T = try SwiftyBoost.SpecialFunctions.modifiedBesselI(v: order, x: z)
            precondition(I0 > 0, "Bessel I must be positive for log derivative")
            let dI_dOrder = besselIOrderDerivative(order, z)
            return dI_dOrder / I0
        }
        catch _ {
            return .nan
        }
    }
}
