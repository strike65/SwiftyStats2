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

/// Score helpers for the Hyperxponential distribution used by the maximum likelihood fitters.

import SwiftyStatsPrelude

/// Score helpers for the Hyperexponential distribution used by MLEFitter.
extension SwiftyBoost.Distribution.Hyperexponential {

    // MARK: - Helpers

    @inline(__always)
    private func softmax(logits: [T]) -> [T] {
        let m = logits.max() ?? .zero
        var e = [T](repeating: .zero, count: logits.count)
        var s: T = .zero
        for i in 0..<logits.count { e[i] = T.exp(logits[i] - m); s += e[i] }
        var p = e
        for i in 0..<p.count { p[i] = e[i] / s }
        return p
    }

    @inline(__always)
    private func responsibilities(x: T, probs: [T], rates: [T]) -> (w: [T], f: T, g: [T]) {
        // g_k = rate_k * exp(-rate_k * x), f = sum probs_k * g_k, w_k = probs_k * g_k / f
        let K = probs.count
        var g = [T](repeating: .zero, count: K)
        var f: T = .zero
        for k in 0..<K {
            let t = rates[k] * T.exp(-rates[k] * x)
            g[k] = t
            f += probs[k] * t
        }
        var w = g
        if f > .zero && f.isFinite {
            for k in 0..<K { w[k] = (probs[k] * g[k]) / f }
        } else {
            for k in 0..<K { w[k] = .nan }
        }
        return (w, f, g)
    }

    // MARK: - Constrained parameters: probs on simplex, rates > 0

    /// Score for a single observation x >= 0.
    ///
    /// Density: f(x) = sum_k probs_k * rates_k * exp(-rates_k * x)
    ///
    /// Score:
    ///  - d/d(probs_k): g_k / f(x), where g_k = rates_k * exp(-rates_k * x)
    ///  - d/d(rates_k): w_k * (1/rates_k - x), where w_k = posterior responsibility
    ///
    /// Note: The gradient w.r.t. probs lives on the simplex (sum_k probs_k * d/d(probs_k) = 1).
    public func score(
        x: T,
        probabilities probs: [T],
        rates: [T]
    ) -> (dProbabilities: [T], dRates: [T]) {
        precondition(x >= .zero, "x must be >= 0")
        precondition(!probs.isEmpty && probs.count == rates.count, "size mismatch")
        for r in rates { precondition(r > .zero && r.isFinite, "rates must be > 0 and finite") }
        for p in probs { precondition(p >= .zero && p.isFinite, "probabilities must be >= 0 and finite") }

        let (w, f, g) = responsibilities(x: x, probs: probs, rates: rates)
        precondition(f > .zero && f.isFinite, "mixture density must be positive and finite")

        // d/d(probs_k) = g_k / f
        var dP = g
        for k in 0..<g.count { dP[k] = g[k] / f }

        // d/d(rates_k) = w_k * (1/rates_k - x)
        var dR = w
        for k in 0..<w.count { dR[k] = w[k] * (T.one / rates[k] - x) }

        return (dP, dR)
    }

    /// Sum of scores over data (constrained params).
    public func totalScore(
        data: [T],
        probabilities probs: [T],
        rates: [T]
    ) -> (dProbabilities: [T], dRates: [T]) {
        precondition(!data.isEmpty, "data must be non-empty")
        let K = probs.count
        var sP = [T](repeating: .zero, count: K)
        var sR = [T](repeating: .zero, count: K)
        for x in data {
            let (a, b) = score(x: x, probabilities: probs, rates: rates)
            for k in 0..<K { sP[k] += a[k]; sR[k] += b[k] }
        }
        return (sP, sR)
    }

    // MARK: - Unconstrained reparameterization: logits for probs, logRates for rates

    /// Score with logits (probs = softmax(logits)) and logRates (rates = exp(logRates)).
    ///
    /// Gradients:
    ///  - d/d(logits_k)  = w_k - probs_k
    ///  - d/d(logRates_k)= w_k * (1 - rates_k * x)
    public func scoreWithLogitsAndLogRates(
        x: T,
        logits: [T],
        logRates: [T]
    ) -> (dLogits: [T], dLogRates: [T]) {
        precondition(x >= .zero, "x must be >= 0")
        precondition(!logits.isEmpty && logits.count == logRates.count, "size mismatch")

        // probs and rates
        let probs = softmax(logits: logits)
        var rates = [T](repeating: .zero, count: logits.count)
        for k in 0..<rates.count { rates[k] = T.exp(logRates[k]) }

        // responsibilities
        let (w, f, _) = responsibilities(x: x, probs: probs, rates: rates)
        guard f > .zero && f.isFinite else {
            let nanVec = [T](repeating: .nan, count: logits.count)
            return (nanVec, nanVec)
        }

        // gradients
        var dLogits = w
        for k in 0..<dLogits.count { dLogits[k] = w[k] - probs[k] }

        var dLogRates = w
        for k in 0..<dLogRates.count { dLogRates[k] = w[k] * (T.one - rates[k] * x) }

        return (dLogits, dLogRates)
    }

    /// Sum of scores over data (unconstrained params).
    public func totalScoreWithLogitsAndLogRates(
        data: [T],
        logits: [T],
        logRates: [T]
    ) -> (dLogits: [T], dLogRates: [T]) {
        precondition(!data.isEmpty, "data must be non-empty")
        let K = logits.count
        var sA = [T](repeating: .zero, count: K)
        var sR = [T](repeating: .zero, count: K)
        for x in data {
            let (a, b) = scoreWithLogitsAndLogRates(x: x, logits: logits, logRates: logRates)
            for k in 0..<K { sA[k] += a[k]; sR[k] += b[k] }
        }
        return (sA, sR)
    }
}
