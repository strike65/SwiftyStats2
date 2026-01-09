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

/// Score helpers for the StudentTNC distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the NonCentralStudentT distribution used by MLEFitter.
extension SwiftyBoost.Distribution.NonCentralStudentT {

    // ---------- finite-difference utilities ----------
    /// Heuristic central-difference step for parameter derivatives, scaled by magnitude.
    ///
    /// - Stability: uses √eps scaling to balance truncation and rounding errors.
    @inline(__always)
    private func diffStepParam(_ p: T, scale: T) -> T {
        // step ~ sqrt(eps) * (|p| + scale)
        let eps = T.ulpOfOne
        let mag = p.magnitude
        let sc  = T.maximum(scale, T.one)
        let s   = T.sqrt(eps) * (mag + sc)
        return s
    }

    /// Clamp ν to a small positive floor to keep degrees of freedom valid during perturbation.
    @inline(__always)
    private func clampNuPositive(_ nu: T) -> T {
        // Use a small, safe positive floor to avoid invalid ν during perturbation
        let floor: T = T(1e-8)
        return nu > floor ? nu : floor
    }

    /// Helper to instantiate a NonCentralStudentT safely, returning nil upon construction failure.
    @inline(__always)
    private func makeDist(nu: T, delta: T) -> SwiftyBoost.Distribution.NonCentralStudentT<T>? {
        do {
            return try SwiftyBoost.Distribution.NonCentralStudentT<T>(degreesOfFreedom: nu, nonCentrality: delta)
        } catch {
            return nil
        }
    }

    /// Finite-difference approximation of ∂/∂ν log f(x | ν, δ) with Richardson refinement.
    ///
    /// - Numerical stability:
    ///   - Uses √eps step scaling.
    ///   - Clamps ν to remain positive.
    ///   - Reuses perturbed distributions.
    @inline(__always)
    private func dLogPdf_dnu(_ x: T, nu: T, delta: T) -> T {
        // Central difference in ν with Richardson refinement; instantiate perturbed distributions.
        var h = diffStepParam(nu, scale: T.one)
        let nuP1 = clampNuPositive(nu + h)
        let nuM1 = clampNuPositive(nu - h)
        guard let dP1 = makeDist(nu: nuP1, delta: delta),
              let dM1 = makeDist(nu: nuM1, delta: delta) else { return T.nan }
        do {
            let fP1: T = try dP1.logPdf(x)
            let fM1: T = try dM1.logPdf(x)
            let d1 = (fP1 - fM1) / (T.two * h)

            h *= T.half
            let nuP2 = clampNuPositive(nu + h)
            let nuM2 = clampNuPositive(nu - h)
            guard let dP2 = makeDist(nu: nuP2, delta: delta),
                  let dM2 = makeDist(nu: nuM2, delta: delta) else { return T.nan }
            let fP2: T = try dP2.logPdf(x)
            let fM2: T = try dM2.logPdf(x)
            let d2 = (fP2 - fM2) / (T.two * h)

            // Richardson: d ≈ d2 + (d2 - d1)/3
            return d2 + (d2 - d1) / T(3)
        } catch {
            return T.nan
        }
    }

    /// Finite-difference approximation of ∂/∂δ log f(x | ν, δ) with Richardson refinement.
    ///
    /// - Numerical stability:
    ///   - Uses √eps step scaling and reused perturbed distributions.
    @inline(__always)
    private func dLogPdf_ddelta(_ x: T, nu: T, delta: T) -> T {
        // Central difference in δ with Richardson refinement; instantiate perturbed distributions.
        var h = diffStepParam(delta, scale: T.one)
        let dP1v = delta + h
        let dM1v = delta - h
        guard let dP1 = makeDist(nu: clampNuPositive(nu), delta: dP1v),
              let dM1 = makeDist(nu: clampNuPositive(nu), delta: dM1v) else { return T.nan }
        do {
            let fP1: T = try dP1.logPdf(x)
            let fM1: T = try dM1.logPdf(x)
            let d1 = (fP1 - fM1) / (T.two * h)

            h *= T.half
            let dP2v = delta + h
            let dM2v = delta - h
            guard let dP2 = makeDist(nu: clampNuPositive(nu), delta: dP2v),
                  let dM2 = makeDist(nu: clampNuPositive(nu), delta: dM2v) else { return T.nan }
            let fP2: T = try dP2.logPdf(x)
            let fM2: T = try dM2.logPdf(x)
            let d2 = (fP2 - fM2) / (T.two * h)

            // Richardson: d ≈ d2 + (d2 - d1)/3
            return d2 + (d2 - d1) / T(3)
        } catch {
            return T.nan
        }
    }

    // ---------- scores (per-observation) ----------
    /// Score w.r.t. degrees of freedom ν and noncentrality δ using finite differences.
    ///
    /// - Constraints: ν > 0.
    /// - Numerical stability: central differences with Richardson refinement.
    @inline(__always)
    public func score(x: T, nu: T, delta: T) -> (dnu: T, ddelta: T) {
        precondition(nu > 0, "nu must be > 0")
        let s1 = dLogPdf_dnu(x, nu: nu, delta: delta)
        let s2 = dLogPdf_ddelta(x, nu: nu, delta: delta)
        let dnu     = s1
        let ddelta  = s2
        return (dnu, ddelta)
    }

    /// Sum of scores (ν, δ) across data using reused finite-difference evaluations.
    ///
    /// - Numerical stability: reuse perturbed distributions, apply Richardson refinement on sums.
    @inline(__always)
    public func totalScore(data: [T], nu: T, delta: T) -> (dnu: T, ddelta: T) {
        precondition(nu > 0, "nu must be > 0")

        // Reuse perturbed distributions across the sample (coarse and fine) for both parameters.

        // ν-derivative setup
        var hNu = diffStepParam(nu, scale: T.one)
        let nuP1 = clampNuPositive(nu + hNu)
        let nuM1 = clampNuPositive(nu - hNu)
        let nuP2 = clampNuPositive(nu + hNu * T.half)
        let nuM2 = clampNuPositive(nu - hNu * T.half)

        // δ-derivative setup
        var hDe = diffStepParam(delta, scale: T.one)
        let deP1 = delta + hDe
        let deM1 = delta - hDe
        let deP2 = delta + hDe * T.half
        let deM2 = delta - hDe * T.half

        guard let dNuP1 = makeDist(nu: nuP1, delta: delta),
              let dNuM1 = makeDist(nu: nuM1, delta: delta),
              let dNuP2 = makeDist(nu: nuP2, delta: delta),
              let dNuM2 = makeDist(nu: nuM2, delta: delta),
              let dDeP1 = makeDist(nu: clampNuPositive(nu), delta: deP1),
              let dDeM1 = makeDist(nu: clampNuPositive(nu), delta: deM1),
              let dDeP2 = makeDist(nu: clampNuPositive(nu), delta: deP2),
              let dDeM2 = makeDist(nu: clampNuPositive(nu), delta: deM2)
        else { return (T.nan, T.nan) }

        var sNu1 = T.zero, sNu2 = T.zero
        var sDe1 = T.zero, sDe2 = T.zero

        do {
            for x in data {
                // ν coarse/fine
                let fNuP1: T = try dNuP1.logPdf(x)
                let fNuM1: T = try dNuM1.logPdf(x)
                sNu1 += (fNuP1 - fNuM1)

                let fNuP2: T = try dNuP2.logPdf(x)
                let fNuM2: T = try dNuM2.logPdf(x)
                sNu2 += (fNuP2 - fNuM2)

                // δ coarse/fine
                let fDeP1: T = try dDeP1.logPdf(x)
                let fDeM1: T = try dDeM1.logPdf(x)
                sDe1 += (fDeP1 - fDeM1)

                let fDeP2: T = try dDeP2.logPdf(x)
                let fDeM2: T = try dDeM2.logPdf(x)
                sDe2 += (fDeP2 - fDeM2)
            }
        } catch {
            return (T.nan, T.nan)
        }

        // Central differences
        let dNu1 = sNu1 / (T.two * hNu)
        hNu *= T.half
        let dNu2 = sNu2 / (T.two * hNu)

        let dDe1 = sDe1 / (T.two * hDe)
        hDe *= T.half
        let dDe2 = sDe2 / (T.two * hDe)

        // Richardson refinement (sum-wise)
        let dNu = dNu2 + (dNu2 - dNu1) / T(3)
        let dDe = dDe2 + (dDe2 - dDe1) / T(3)

        return (dNu, dDe)
    }

    // ---------- log-reparameterizations ----------
    /// Score with ν = exp(logNu).
    ///
    /// - Chain rule: ∂/∂log ν = ν ∂/∂ν; δ unchanged.
    @inline(__always)
    public func scoreWithLogNu(x: T, logNu: T, delta: T) -> (dlogNu: T, ddelta: T) {
        let nu = T.exp(logNu)
        let s1 = dLogPdf_dnu(x, nu: nu, delta: delta)
        let s2 = dLogPdf_ddelta(x, nu: nu, delta: delta)
        let dlogNu = nu * s1
        let ddel   = s2
        return (dlogNu, ddel)
    }

    /// Total score with ν = exp(logNu) across data.
    @inline(__always)
    public func totalScoreWithLogNu(data: [T], logNu: T, delta: T) -> (dlogNu: T, ddelta: T) {
        let nu = T.exp(logNu)

        // Reuse perturbed distributions across sample, then map ν-derivative by ν.
        var hNu = diffStepParam(nu, scale: T.one)
        let nuP1 = clampNuPositive(nu + hNu)
        let nuM1 = clampNuPositive(nu - hNu)
        let nuP2 = clampNuPositive(nu + hNu * T.half)
        let nuM2 = clampNuPositive(nu - hNu * T.half)

        var hDe = diffStepParam(delta, scale: T.one)
        let deP1 = delta + hDe
        let deM1 = delta - hDe
        let deP2 = delta + hDe * T.half
        let deM2 = delta - hDe * T.half

        guard let dNuP1 = makeDist(nu: nuP1, delta: delta),
              let dNuM1 = makeDist(nu: nuM1, delta: delta),
              let dNuP2 = makeDist(nu: nuP2, delta: delta),
              let dNuM2 = makeDist(nu: nuM2, delta: delta),
              let dDeP1 = makeDist(nu: clampNuPositive(nu), delta: deP1),
              let dDeM1 = makeDist(nu: clampNuPositive(nu), delta: deM1),
              let dDeP2 = makeDist(nu: clampNuPositive(nu), delta: deP2),
              let dDeM2 = makeDist(nu: clampNuPositive(nu), delta: deM2)
        else { return (T.nan, T.nan) }

        var sNu1 = T.zero, sNu2 = T.zero
        var sDe1 = T.zero, sDe2 = T.zero

        do {
            for x in data {
                let fNuP1: T = try dNuP1.logPdf(x)
                let fNuM1: T = try dNuM1.logPdf(x)
                sNu1 += (fNuP1 - fNuM1)

                let fNuP2: T = try dNuP2.logPdf(x)
                let fNuM2: T = try dNuM2.logPdf(x)
                sNu2 += (fNuP2 - fNuM2)

                let fDeP1: T = try dDeP1.logPdf(x)
                let fDeM1: T = try dDeM1.logPdf(x)
                sDe1 += (fDeP1 - fDeM1)

                let fDeP2: T = try dDeP2.logPdf(x)
                let fDeM2: T = try dDeM2.logPdf(x)
                sDe2 += (fDeP2 - fDeM2)
            }
        } catch {
            return (T.nan, T.nan)
        }

        let dNu1 = sNu1 / (T.two * hNu)
        hNu *= T.half
        let dNu2 = sNu2 / (T.two * hNu)
        let dNu = dNu2 + (dNu2 - dNu1) / T(3)

        let dDe1 = sDe1 / (T.two * hDe)
        hDe *= T.half
        let dDe2 = sDe2 / (T.two * hDe)
        let dDe = dDe2 + (dDe2 - dDe1) / T(3)

        return (nu * dNu, dDe)
    }

    /// Score with δ = signDelta * exp(logAbsDelta), allowing unconstrained signed δ.
    ///
    /// - Chain rule: ∂/∂log|δ| = |δ| sign(δ) ∂/∂δ = exp(logAbsDelta) * signDelta * ∂/∂δ.
    @inline(__always)
    public func scoreWithSignedLogAbsDelta(
        x: T, nu: T, signDelta: T, logAbsDelta: T
    ) -> (dnu: T, dlogAbsDelta: T) {
        precondition(nu > 0, "nu must be > 0")
        let absDelta = T.exp(logAbsDelta)
        let delta    = signDelta * absDelta

        let g1 = dLogPdf_dnu(x, nu: nu, delta: delta)
        let g2 = dLogPdf_ddelta(x, nu: nu, delta: delta)

        let dnu        = g1
        let dlogAbsDel = absDelta * (g2 * signDelta)   // d delta / d log|delta| = sign * |delta|
        return (dnu, dlogAbsDel)
    }

    /// Total score with signed log-absolute-δ across data.
    @inline(__always)
    public func totalScoreWithSignedLogAbsDelta(
        data: [T], nu: T, signDelta: T, logAbsDelta: T
    ) -> (dnu: T, dlogAbsDelta: T) {
        precondition(nu > 0, "nu must be > 0")
        let absDelta = T.exp(logAbsDelta)
        let delta    = signDelta * absDelta

        // Reuse perturbed distributions for ν and δ
        var hNu = diffStepParam(nu, scale: T.one)
        let nuP1 = clampNuPositive(nu + hNu)
        let nuM1 = clampNuPositive(nu - hNu)
        let nuP2 = clampNuPositive(nu + hNu * T.half)
        let nuM2 = clampNuPositive(nu - hNu * T.half)

        var hDe = diffStepParam(delta, scale: T.one)
        let deP1 = delta + hDe
        let deM1 = delta - hDe
        let deP2 = delta + hDe * T.half
        let deM2 = delta - hDe * T.half

        guard let dNuP1 = makeDist(nu: nuP1, delta: delta),
              let dNuM1 = makeDist(nu: nuM1, delta: delta),
              let dNuP2 = makeDist(nu: nuP2, delta: delta),
              let dNuM2 = makeDist(nu: nuM2, delta: delta),
              let dDeP1 = makeDist(nu: clampNuPositive(nu), delta: deP1),
              let dDeM1 = makeDist(nu: clampNuPositive(nu), delta: deM1),
              let dDeP2 = makeDist(nu: clampNuPositive(nu), delta: deP2),
              let dDeM2 = makeDist(nu: clampNuPositive(nu), delta: deM2)
        else { return (T.nan, T.nan) }

        var sNu1 = T.zero, sNu2 = T.zero
        var sDe1 = T.zero, sDe2 = T.zero

        do {
            for x in data {
                let fNuP1: T = try dNuP1.logPdf(x)
                let fNuM1: T = try dNuM1.logPdf(x)
                sNu1 += (fNuP1 - fNuM1)

                let fNuP2: T = try dNuP2.logPdf(x)
                let fNuM2: T = try dNuM2.logPdf(x)
                sNu2 += (fNuP2 - fNuM2)

                let fDeP1: T = try dDeP1.logPdf(x)
                let fDeM1: T = try dDeM1.logPdf(x)
                sDe1 += (fDeP1 - fDeM1)

                let fDeP2: T = try dDeP2.logPdf(x)
                let fDeM2: T = try dDeM2.logPdf(x)
                sDe2 += (fDeP2 - fDeM2)
            }
        } catch {
            return (T.nan, T.nan)
        }

        let dNu1 = sNu1 / (T.two * hNu)
        hNu *= T.half
        let dNu2 = sNu2 / (T.two * hNu)
        let dNu = dNu2 + (dNu2 - dNu1) / T(3)

        let dDe1 = sDe1 / (T.two * hDe)
        hDe *= T.half
        let dDe2 = sDe2 / (T.two * hDe)
        let dDe = dDe2 + (dDe2 - dDe1) / T(3)

        // Chain for log|δ|: ∂/∂log|δ| = |δ| * sign * ∂/∂δ
        let dlogAbsDel = absDelta * (dDe * signDelta)

        return (dNu, dlogAbsDel)
    }
}
