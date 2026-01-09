//
//  Created by VT on 15.12.25.
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

/// Tukey's studentised range distribution helpers (`ptukey`).
///
/// Provides a numerically stable approximation of the lower-tail CDF for
/// Tukey's studentised range statistic, following the structure of R's
/// `stats::ptukey` (Gauss-Legendre quadrature for both the normal-range CDF
/// and the finite-`df` mixture).

//  Mathematical summary
//  --------------------
//  Let Q denote Tukey's studentized range statistic for c means and df error
//  degrees of freedom. The lower-tail CDF P(Q <= q) can be expressed via a
//  mixture over a scale parameter arising from the t/chi-square components.
//
//  For df -> infinity (known variance), the distribution reduces to the
//  range of i.i.d. standard normals. Using Hartley's form, the CDF for a
//  single range of width w is
//
//     W(w; c) = [ 2 Phi(w/2) - 1 ]^c
//               + c * ∫ phi(x) * [ Phi(x) - Phi(x - w) ]^(c-1) dx,  x ∈ [w/2, ∞),
//
//  where Phi and phi denote the standard normal CDF and PDF. For rr
//  independent ranges, the probability becomes W(w; c)^rr. We evaluate W via
//  Gauss–Legendre quadrature on sub-intervals of [w/2, 8], plus the closed
//  form term.
//
//  For finite df, we integrate W over a mixing density that depends on df.
//  This file computes
//
//     P(Q <= q) ≈ ∫ W( q * sqrt(u/2) ; c )^rr * g_df(u) du
//
//  using 16-point Gauss–Legendre quadrature in the outer integral with a
//  stable log-weight formulation. The implementation follows the same
//  structure used in common statistical libraries but is written afresh for
//  Swift with careful cutoffs for numerical stability (underflow/overflow).
//
//  Numerical details
//  -----------------
//  - Inner integral (W): 12-point Gauss–Legendre on 2–3 sub-intervals.
//  - Outer integral (mixture over df): 16-point Gauss–Legendre; sub-interval
//    size chosen by df; large-df shortcut uses the normal-range approximation.
//  - Thresholds on exponential terms and small contributions are applied to
//    avoid denormals while keeping relative error small across the domain.
//  - Tail/log-p handling is delegated to r_derived_macros for consistency.
//
//  References
//  ----------
//  - Hartley, H. O. (1950s/60s work on the range distribution)
//  - Copenhaver, M. D., & Holland, B. S. (1988). Multiple comparisons of
//    simple effects in two-way ANOVA with fixed effects.
//  - Stroud, A. H., & Secrest, D. (1966). Gaussian Quadrature Formulas.
//  - Algorithmic structure inspired by Copenhaver & Holland (1988) and
//    the reference implementation in R's stats::ptukey.
//  - Gauss–Legendre abscissae/weights taken from Stroud & Secrest (1966).
//

import SwiftyStatsPrelude


extension Inferential.HypothesisTesting.ParametricTests {
    
    internal static func ptukey(q: T, nranges: T, numberOfMeans: T, df: T, tail: Tail, returnLogP: Bool) throws -> T {
        return try ss_ptukey(q: q, rr: nranges, cc: numberOfMeans, df: df, tail: tail, logP: returnLogP)
    }
    
    // MARK: - Core numerics
    
    /// Compute W(w; c)^rr where W(w; c) is the lower-tail CDF of the normal-range
    /// distribution (df = ∞ case) for range width w and number of means c.
    ///
    /// Implements Hartley’s decomposition:
    ///   W(w; c) = [2 Phi(w/2) - 1]^c + c ∫ phi(x) [Phi(x) - Phi(x-w)]^(c-1) dx,
    /// with the integral evaluated by Gauss–Legendre over x ∈ [w/2, 8].
    /// Returns W(w; c)^rr with early cutoffs for very small/large probabilities.
    fileprivate static func ss_wprob(_ w: T, rr: T, cc: T) throws -> T {
        // Guard trivial/limit cases
        let qsqz = 0.5 * w
        if qsqz <= 0 { return 0 }
        if qsqz >= 8.0 { return 1.0 }
        
        // Gauss–Legendre (order 12) nodes/weights on (0, 1)
        let nleg = 12
        let ihalf = nleg / 2
        let xleg: [T] = [
            0.981560634246719250690549090149,
            0.904117256370474856678465866119,
            0.769902674194304687036893833213,
            0.587317954286617447296702418941,
            0.367831498998180193752691536644,
            0.125233408511468915472441369464
        ]
        let aleg: [T] = [
            0.047175336386511827194615961485,
            0.106939325995318430960254718194,
            0.160078328543346226334652529543,
            0.203167426723065921749064455810,
            0.233492536538354808760849898925,
            0.249147045813402785000562436043
        ]
        
        // First Hartley term: [2 Phi(w/2) - 1]^c
        let gaussian: Distribution.Normal<T> = try .init(mean: .zero, sd: .one)
        var pr_w: T = T.two * (try gaussian.cdf(qsqz)) - T.one
        if pr_w > .zero {
            pr_w = T.pow(pr_w, cc)
        } else {
            pr_w = .zero
        }
        
        // Choose number of subintervals for the second term
        let wlar: T = 3.0
        let wincr: T = (w > wlar) ? T.two : T(3.0)
        
        // Integrate f(x) = c * phi(x) * [Phi(x) - Phi(x - w)]^(c-1) over x in [w/2, 8]
        var blb: T = qsqz
        var bub = T.minimum(8.0, blb + wincr)
        var einsum = T.zero
        let cc1 = cc - T.one
        while blb < 8.0 - T.ulpOfOne {
            // Midpoint/half-length for affine map from u in [-1,1]
            let a = T.half * (bub + blb)
            let b = T.half * (bub - blb)
            var elsum = T.zero
            var jj = 1
            while jj <= nleg {
                let j = (jj > ihalf) ? (nleg - jj + 1) : jj
                let xx: T = (jj > ihalf) ? xleg[j - 1] : -xleg[j - 1]
                let ac: T = a + b * xx
                // If exp(-ac^2/2) is too small, break
                let qexpo: T = ac * ac
                if qexpo > 90.0 { break }
                // Compute rinsum = Phi(ac) - Phi(ac - w)
                let pplus: T = try gaussian.cdf(ac)
                let pminus: T = try gaussian.cdf(ac - w)
                let rinsum: T = pplus - pminus
                // Ignore small contributions (more permissive)
                if rinsum > .zero {
                    // Add Legendre contribution
                    let term: T = aleg[j - 1] * T.exp(-T.half * qexpo) * T.pow(rinsum, cc1)
                    elsum += term
                }
                jj += 1
            }
            // Scale by mapping Jacobian and constants
            elsum *= (2.0 * b) * cc * T.sqrt2PiInv
            einsum += elsum
            blb = bub
            bub = T.minimum(8.0, bub + wincr)
        }
        pr_w += einsum
        // Compute W^rr in log space to avoid underflow
        if pr_w <= .zero { return .zero }
        let logW = T.log(pr_w)
        let logWrr = rr * logW
        // Guard against overflow/underflow in exp
        if logWrr <= -T(745) { return .zero } // ~ Double underflow threshold
        let wrr = T.exp(logWrr)
        return wrr >= .one ? .one : wrr
    }
    
    /// Lower/upper tail CDF for Tukey's studentized range with parameters:
    ///  - q: statistic value (q > 0)
    ///  - rr: number of independent ranges (nranges ≥ 1, internally floored)
    ///  - cc: number of means per range (c ≥ 2, internally floored)
    ///  - df: degrees of freedom for the error term (df ≥ 2)
    ///  - tail: lower vs upper tail selection
    ///  - logP: whether to return log-probability
    ///
    /// For df > 25000 we use the normal-range approximation W, otherwise a
    /// Gauss–Legendre mixture integral in a stable log-weight form.
    fileprivate static func ss_ptukey(q: T, rr: T, cc: T, df: T, tail: Tail, logP: Bool) throws -> T {
        // Guard invalid input
        if q <= .zero {
            switch tail {
            case .lower: return logP ?  -T.infinity : T.zero
            case .upper: return logP ? T.zero : T.one
            }
        }
        if df < 2 || rr < 1 || cc < 2 { return T.nan }
        if !q.isFinite {
            switch tail {
            case .lower: return logP ? T.zero : T.one
            case .upper: return logP ? -T.infinity : T.zero
            }
        }
        
        // Large df: approximate with range distribution of standard normal
        if df > 25000 {
            let w = try ss_wprob(q, rr: rr, cc: cc)
            switch tail {
            case .lower: return logP ? T.log(T.maximum(w, .zero)) : w
            case .upper: return logP ? T.log(T.maximum(.one - w, .zero)) : .one - w
            }
        }
        
        // Gauss–Legendre (order 16) for outer integral
        let nlegq = 16
        let ihalfq = 8
        let xlegq: [T] = [
            0.989400934991649932596154173450,
            0.944575023073232576077988415535,
            0.865631202387831743880467897712,
            0.755404408355003033895101194847,
            0.617876244402643748446671764049,
            0.458016777657227386342419442984,
            0.281603550779258913230460501460,
            0.0950125098376374401853193354250
        ]
        let alegq: [T] = [
            0.0271524594117540948517805724560,
            0.0622535239386478928628438369944,
            0.0951585116824927848099251076022,
            0.124628971255533872052476282192,
            0.149595988816576732081501730547,
            0.169156519395002538189312079030,
            0.182603415044923588866763667969,
            0.189450610455068496285396723208
        ]
        
        // Integration parameters per df (more permissive thresholds)
        let eps1: T = -80.0
        let eps2: T = T(1e-16)
        let dhaf: T = 100.0
        let dquar: T = 800.0
        let deigh: T = 5000.0
        let ulen1: T = .one, ulen2: T = T.half, ulen3: T = .quarter, ulen4: T = T(0.125)
        
        let f2: T = .half * df
        let f2lf: T = (f2 * T.log(df)) - (df * T.lnTwo) - (try SpecialFunctions.logGamma(f2))
        let f21: T = f2 - .one
        
        var ulen: T
        if df <= dhaf { ulen = ulen1 }
        else if df <= dquar { ulen = ulen2 }
        else if df <= deigh { ulen = ulen3 }
        else { ulen = ulen4 }
        
        func integrateOnce(ulen: T) throws -> T {
            var ansLoc: T = T.zero
            var i = 1
            while i <= 100 {
                var otsum: T = .zero
                let twa1: T = (.two * T(i) - .one) * ulen
                var jj = 1
                while jj <= nlegq {
                    let j: Int
                    let signPlus: Bool = (jj > ihalfq)
                    if signPlus {
                        j = jj - ihalfq - 1
                    } else {
                        j = jj - 1
                    }
                    let x: T = xlegq[j]
                    let xl: T = x * ulen
                    let expoTerm: T
                    let arg: T
                    if signPlus {
                        expoTerm = f2lf + f21 * T.log(twa1 + xl) - (xl + twa1) * (T.quarter * df)
                        arg = (x * ulen + twa1) * T.half
                    } else {
                        expoTerm = f2lf + f21 * T.log(twa1 - xl) + (xl - twa1) * (T.quarter * df)
                        arg = (-x * ulen + twa1) * T.half
                    }
                    if expoTerm >= eps1 {
                        let a: T = T.maximum(.zero, arg)
                        let qsqz: T = q * T.sqrt(a)
                        let wprb: T = try ss_wprob(qsqz, rr: rr, cc: cc)
                        let weight = alegq[j] * T.exp(expoTerm)
                        if wprb.isFinite && weight.isFinite {
                            let contrib: T = wprb * weight
                            otsum += contrib
                        }
                    }
                    jj += 1
                }
                if T(i) * ulen > 1.0 && otsum <= eps2 {
                    // perform one more iteration to reduce bias before breaking
                    ansLoc += (ulen * otsum)
                    break
                }
                // Gauss–Legendre on [-1,1] scaled to segment length; multiply by ulen (Jacobian)
                ansLoc += (ulen * otsum)
                i += 1
            }
            return T.minimum(1.0, ansLoc)
        }
        var ans: T = try integrateOnce(ulen: ulen)
        if ans == 0.0 {
            let finer: T = ulen * T.half
            ans = try integrateOnce(ulen: finer)
        }
        if !(ans.isFinite) { ans = .nan }
        if ans < .zero { ans = .zero }
        if ans > .one { ans = .one }
        switch tail {
        case .lower: return logP ? T.log(T.maximum(ans, .zero)) : ans
        case .upper: let comp = T.one - ans; return logP ? T.log(T.maximum(comp, .zero)) : comp
        }
    }
    /// Quantile of Tukey's studentized range distribution.
    ///
    /// Solves for q such that P(Q <= q) = p (for the requested tail) using a
    /// robust bracketing+bisection scheme. The CDF P(Q <= q) is evaluated by
    /// the clean implementation `ptukey_new` described in SSTukeyNew.swift.
    ///
    /// Mathematical notes
    /// ------------------
    /// - The CDF in q is continuous and strictly increasing on q > 0, so
    ///   bracketing by doubling the upper endpoint followed by bisection is
    ///   guaranteed to converge.
    /// - Input p is converted to the lower-tail probability via r_derived
    ///   helpers, supporting both tail selection and log-prob input.
    /// - Boundary handling: p <= 0 → q = 0; p >= 1 → q = +∞.
    internal static func qtukey(p: T, nranges: T, numberOfMeans: T, df: T, tail: Tail, log_p: Bool) throws -> T {
        // NaN handling
        if p.isNaN || nranges.isNaN || numberOfMeans.isNaN || df.isNaN { return .nan }
        if df < 2 || nranges < 1 || numberOfMeans < 2 { return .nan }

        // Convert input probability to lower-tail probability in [0, 1]
        let pLinear: T = {
            if log_p {
                return T.exp(p)
            }
            else {
                return p
            }
        }()
        let lowerTailP: T = {
            switch tail {
            case .lower: return pLinear
            case .upper: return T.one - pLinear
            }
        }()
        let pl: T = T.minimum(.one, T.maximum(.zero, lowerTailP))
        if pl <= 0 { return 0 }
        if pl >= 1 || !pl.isFinite { return T.infinity }

        // Normalize discrete params to valid domain
        let rr = T.maximum(.one, floor(nranges))
        let cc = T.maximum(.two, floor(numberOfMeans))

        // Bracket [lo, hi] with F(lo) <= pl <= F(hi)
        var lo: T = 0.0
        var hi: T = 1.0
        func cdf(_ q: T) throws -> T { try ptukey(q: q, nranges: rr, numberOfMeans: cc, df: df, tail: .lower, returnLogP: false) }

        var fhi = try cdf(hi)
        var iter = 0
        while fhi < pl && iter < 64 {
            hi *= 2
            fhi = try cdf(hi)
            iter += 1
        }

        // Bisection solve (monotone in q)
        let tolProb: T = 1e-10
        let tolQ: T = 1e-10
        var mid: T = hi
        while iter < 256 {
            mid = 0.5 * (lo + hi)
            let fmid: T = try cdf(mid)
            if (fmid - pl).magnitude <= tolProb || (hi - lo) <= tolQ { break }
            if fmid < pl {
                lo = mid
            } else {
                hi = mid
                fhi = fmid
            }
            iter += 1
        }
        return mid
    }
}
