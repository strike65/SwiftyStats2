//
//  Created by VT on 08.11.25.
//  © 2025 Volker Thieme. All rights reserved.
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

/// Maximum likelihood fitter for the Hyperexponential distribution in SwiftyStats.

import SwiftyStatsPrelude
extension MLEFitter {
    /// Fit a Hyperexponential mixture (finite mixture of exponentials) using numerical MLE with logits/log-rate reparameterization.
    ///
    /// Parameterization: probabilities `p_k` on the simplex (Σ p_k = 1) and strictly positive rates `λ_k`.
    ///
    /// - Parameters:
    ///   - data: Sample values. Must be non-empty, finite, and ≥ 0.
    ///   - phaseCount: Number of exponential phases in the mixture (≥ 1).
    ///   - optimizer: Optimizer to use (default: `.lbfgs` to leverage the analytic score).
    ///   - options: Optional optimizer configuration; sensible defaults with multi-starts are used when `nil`.
    ///   - initialProbabilities: Optional warm-start probabilities; will be normalised and clipped to the simplex.
    ///   - initialRates: Optional warm-start rates per phase.
    /// - Returns: `MLEResult` with `thetaHat = [p̂₀,…,p̂_{K-1}, λ̂₀,…,λ̂_{K-1}]`.
    /// - Implementation notes:
    ///   - Internal optimization happens in unconstrained logits/log-rates; outputs are mapped back to simplex/rates.
    ///   - Uses the analytic Hyperexponential score from `Score-functions.swift` for faster convergence.
    public static func fitHyperexponential(
        _ data: [T],
        phaseCount: Int,
        optimizer: OptimizerKind = .lbfgs,
        options: MLEOptimizationOpts<T>? = nil,
        initialProbabilities: [T]? = nil,
        initialRates: [T]? = nil
    ) -> MLEResult<T> {
        precondition(!data.isEmpty, "data must not be empty")
        precondition(phaseCount >= 1, "phaseCount must be >= 1")
        precondition(data.allSatisfy { $0.isFinite && $0 >= .zero }, "Hyperexponential fits require finite, nonnegative data")

        let phases = phaseCount
        let totalParams = phases * 2
        let probFloor: T = T(1e-6)
        let rateFloor: T = T(1e-6)

        @inline(__always)
        func percentile(_ sorted: [T], _ p: T) -> T {
            guard !sorted.isEmpty else { return T.zero }
            if sorted.count == 1 { return sorted[0] }
            let clamped = max(T.zero, min(T.one, p))
            let r = T(sorted.count - 1) * clamped
            var lo = Int(r.rounded(.down))
            var hi = Int(r.rounded(.up))
            lo = max(0, min(lo, sorted.count - 1))
            hi = max(0, min(hi, sorted.count - 1))
            let lower = sorted[lo]
            if lo == hi { return lower }
            let upper = sorted[hi]
            let frac = r - T(lo)
            return lower + (upper - lower) * frac
        }

        var sumDataInput = data
        let sumData = Helpers.sum(&sumDataInput)
        let meanData = max(sumData / T(data.count), rateFloor)
        var sortedData = data
        sortedData.sort()

        let baseProbabilities: [T] = {
            if let provided = initialProbabilities {
                precondition(provided.count == phases, "initialProbabilities.count must match phaseCount")
                var sanitized = [T](repeating: T.zero, count: phases)
                var total: T = .zero
                for i in 0..<phases {
                    let value = provided[i]
                    precondition(value.isFinite, "initialProbabilities must be finite")
                    let clipped = max(value, probFloor)
                    sanitized[i] = clipped
                    total += clipped
                }
                if total <= .zero {
                    return Array(repeating: T.one / T(phases), count: phases)
                }
                for i in 0..<phases { sanitized[i] /= total }
                return sanitized
            } else {
                return Array(repeating: T.one / T(phases), count: phases)
            }
        }()

        let baseRates: [T] = {
            if let provided = initialRates {
                precondition(provided.count == phases, "initialRates.count must match phaseCount")
                var sanitized = [T](repeating: rateFloor, count: phases)
                for i in 0..<phases {
                    let value = provided[i]
                    precondition(value.isFinite, "initialRates must be finite")
                    sanitized[i] = max(value, rateFloor)
                }
                return sanitized
            } else {
                var seeds: [T] = []
                seeds.reserveCapacity(phases)
                for k in 0..<phases {
                    let p = (T(k) + T.one) / T(phases + 1)
                    var q = percentile(sortedData, p)
                    if !q.isFinite || q <= rateFloor {
                        q = meanData * (T.one + T(k) / T(max(phases, 1)))
                    }
                    let denom = max(q, meanData * T(0.5))
                    let candidate = denom > rateFloor ? max(T.one / denom, rateFloor) : max(T.one / meanData, rateFloor)
                    seeds.append(candidate)
                }
                return seeds
            }
        }()

        var probSeed = baseProbabilities.map { max($0, probFloor) }
        var probTotal: T = .zero
        for p in probSeed { probTotal += p }
        if probTotal <= .zero {
            probSeed = Array(repeating: T.one / T(phases), count: phases)
            probTotal = T.one
        }
        for i in 0..<phases { probSeed[i] /= probTotal }
        let logits0 = probSeed.map { T.log($0) }

        var rateSeed = baseRates.map { max($0, rateFloor) }
        for i in 0..<phases {
            if !rateSeed[i].isFinite || rateSeed[i] <= rateFloor {
                rateSeed[i] = rateFloor * T(10)
            }
        }
        let logRates0 = rateSeed.map { T.log($0) }

        var specs: [ParamSpec<T>] = []
        specs.reserveCapacity(totalParams)
        for value in logits0 {
            specs.append(.init(.real, initial: value, step: T(0.5)))
        }
        for value in logRates0 {
            specs.append(.init(.real, initial: value, step: T(0.5)))
        }

        typealias DecodedTheta = (logits: [T], logRates: [T], probs: [T], rates: [T], logitNormalizer: T)

        @Sendable @inline(__always)
        func decodeTheta(_ theta: [T]) -> DecodedTheta? {
            guard theta.count == totalParams else { return nil }
            let logits = Array(theta[0..<phases])
            let logRates = Array(theta[phases..<totalParams])
            guard logits.allSatisfy({ $0.isFinite }) else { return nil }
            var rates = [T](repeating: .zero, count: phases)
            for i in 0..<phases {
                let lr = logRates[i]
                guard lr.isFinite else { return nil }
                let rate = T.exp(lr)
                guard rate > .zero && rate.isFinite else { return nil }
                rates[i] = rate
            }
            let maxLogit = logits.max() ?? T.zero
            guard maxLogit.isFinite else { return nil }
            var shifted = [T](repeating: .zero, count: phases)
            var sumExp: T = .zero
            for i in 0..<phases {
                let delta = logits[i] - maxLogit
                guard delta.isFinite else { return nil }
                let e = T.exp(delta)
                shifted[i] = e
                sumExp += e
            }
            guard sumExp > .zero && sumExp.isFinite else { return nil }
            var probs = [T](repeating: .zero, count: phases)
            for i in 0..<phases {
                let p = shifted[i] / sumExp
                guard p > .zero && p.isFinite else { return nil }
                probs[i] = p
            }
            let logNorm = maxLogit + T.log(sumExp)
            return (logits, logRates, probs, rates, logNorm)
        }

        let logPDF: @Sendable (T, [T]) -> T = { x, theta in
            guard x.isFinite && x >= .zero else { return T.infinity }
            guard let decoded = decodeTheta(theta) else { return T.infinity }
            var terms = [T](repeating: .zero, count: phases)
            var maxTerm: T = -.infinity
            for idx in 0..<phases {
                let logProb = decoded.logits[idx] - decoded.logitNormalizer
                let logRate = decoded.logRates[idx]
                let term = logProb + logRate - decoded.rates[idx] * x
                terms[idx] = term
                if term > maxTerm { maxTerm = term }
            }
            guard maxTerm.isFinite else { return T.infinity }
            var sumExp: T = .zero
            for term in terms {
                let delta = term - maxTerm
                guard delta.isFinite else { return T.infinity }
                sumExp += T.exp(delta)
            }
            guard sumExp > .zero && sumExp.isFinite else { return T.infinity }
            return maxTerm + T.log(sumExp)
        }

        let gradLogPDF: @Sendable (T, [T]) -> [T] = { x, theta in
            guard x.isFinite && x >= .zero else { return Array(repeating: T.nan, count: totalParams) }
            guard let decoded = decodeTheta(theta) else { return Array(repeating: T.nan, count: totalParams) }
            do {
                let dist = try SwiftyBoost.Distribution.Hyperexponential(probabilities: decoded.probs, rates: decoded.rates)
                let g = dist.scoreWithLogitsAndLogRates(x: x, logits: decoded.logits, logRates: decoded.logRates)
                var grad = g.dLogits
                grad.append(contentsOf: g.dLogRates)
                return grad
            } catch {
                return Array(repeating: T.nan, count: totalParams)
            }
        }

        let problem = MLEProblem<T>(data: data, logpdf: logPDF, gradlogpdf: gradLogPDF, paramSpecs: specs)

        // -------------------------
        // Internal EM bootstrap (logits/log-rate space)
        // -------------------------
        @inline(__always)
        func emWarmStart(_ theta0: [T]) -> [T] {
            // Decode; if fails use theta0 directly.
            guard var dec = decodeTheta(theta0) else { return theta0 }

            let n = data.count
            let nT = T(n)
            // EM hyperparameters
            let maxIter = 200
            let tolRelLL: T = T(1e-7)
            let tolParam: T = T(1e-6)

            // Working arrays
            var logits = dec.logits
            var logRates = dec.logRates
            var probs = dec.probs
            var rates = dec.rates

            // Precompute nothing; iterate
            var prevAvgLL: T = -.infinity

            // Utility: compute average log-likelihood under current params
            func avgLogLik() -> T {
                var s: T = 0
                for x in data {
                    var maxTerm: T = -.infinity
                    var terms = [T](repeating: .zero, count: phases)
                    for k in 0..<phases {
                        let term = logits[k] - dec.logitNormalizer + logRates[k] - rates[k] * x
                        terms[k] = term
                        if term > maxTerm { maxTerm = term }
                    }
                    var sumExp: T = 0
                    for t in terms { sumExp += T.exp(t - maxTerm) }
                    s += maxTerm + T.log(sumExp)
                }
                return s / nT
            }

            // Keep logits centered to avoid drift: subtract log-sum-exp of logits each M-step
            @inline(__always)
            func recenterLogits(_ v: inout [T]) -> T {
                var m = v.max() ?? T.zero
                if !m.isFinite { m = 0 }
                var sumExp: T = 0
                for i in 0..<phases { sumExp += T.exp(v[i] - m) }
                let lse = m + T.log(sumExp)
                for i in 0..<phases { v[i] -= lse }
                return lse
            }

            // Start with consistent normalizer
            dec.logitNormalizer = recenterLogits(&logits)
            // Rates already physical; ensure floors
            for k in 0..<phases { rates[k] = max(rates[k], rateFloor) }
            for k in 0..<phases { logRates[k] = T.log(rates[k]) }
            // Probabilities from centered logits
            var shifted = [T](repeating: .zero, count: phases)
            var sumExpLogits: T = 0
            for k in 0..<phases { shifted[k] = T.exp(logits[k]); sumExpLogits += shifted[k] }
            for k in 0..<phases { probs[k] = shifted[k] / sumExpLogits }

            var iter = 0
            while iter < maxIter {
                iter += 1

                // E-step: responsibilities r_{i,k} proportional to p_k * λ_k * exp(-λ_k x_i)
                // We do it stably in log space using logits/logRates.
                var Nk = [T](repeating: 0, count: phases)
                var Sxk = [T](repeating: 0, count: phases)

                for x in data {
                    // Compute component log-terms and stabilize
                    var maxTerm: T = -.infinity
                    var terms = [T](repeating: .zero, count: phases)
                    for k in 0..<phases {
                        let term = logits[k] + logRates[k] - rates[k] * x // since logits already centered
                        terms[k] = term
                        if term > maxTerm { maxTerm = term }
                    }
                    var denom: T = 0
                    for k in 0..<phases { denom += T.exp(terms[k] - maxTerm) }

                    // Accumulate responsibilities and sufficient statistics
                    for k in 0..<phases {
                        let rik = T.exp(terms[k] - maxTerm) / denom
                        Nk[k] += rik
                        Sxk[k] += rik * x
                    }
                }

                // M-step:
                // p_k = Nk / n, λ_k = Nk / Sxk
                var maxLogitRaw: T = -.infinity
                var newLogits = [T](repeating: .zero, count: phases)
                var newLogRates = [T](repeating: .zero, count: phases)
                var newRates = [T](repeating: .zero, count: phases)
                for k in 0..<phases {
                    var pk = Nk[k] / nT
                    pk = max(pk, probFloor)
                    // we store logits unnormalized; recenter later
                    let logitRaw = T.log(pk)
                    newLogits[k] = logitRaw
                    if logitRaw > maxLogitRaw { maxLogitRaw = logitRaw }

                    let denom = max(Sxk[k], rateFloor)
                    var lam = Nk[k] / denom
                    lam = max(lam, rateFloor)
                    newRates[k] = lam
                    newLogRates[k] = T.log(lam)
                }

                // Recenter logits to keep numerical stability
                let lse = recenterLogits(&newLogits)

                // Convergence checks: parameters and average log-likelihood
                var maxDelta: T = 0
                for k in 0..<phases {
                    maxDelta = max(maxDelta, abs(newLogits[k] - logits[k]))
                    maxDelta = max(maxDelta, abs(newLogRates[k] - logRates[k]))
                }

                logits = newLogits
                logRates = newLogRates
                rates = newRates
                probs.removeAll(keepingCapacity: true)
                var sumExpL: T = 0
                for k in 0..<phases { sumExpL += T.exp(logits[k]) }
                for k in 0..<phases { probs.append(T.exp(logits[k]) / sumExpL) }
                dec.logitNormalizer = lse

                let avgLL = avgLogLik()
                let rel = abs(avgLL - prevAvgLL) / max(T.one, abs(avgLL))
                prevAvgLL = avgLL

                if rel < tolRelLL || maxDelta < tolParam {
                    break
                }
            }

            // Return unconstrained θ = [logits || logRates]
            var theta = logits
            theta.append(contentsOf: logRates)
            return theta
        }

        // Seed θ0 in unconstrained space from specs.initial
        var theta0: [T] = []
        theta0.reserveCapacity(totalParams)
        for i in 0..<phases { theta0.append(specs[i].initial) }              // logits
        for i in 0..<phases { theta0.append(specs[phases + i].initial) }     // logRates

        // Run EM bootstrap unless caller supplied a warm start
        let thetaEM = emWarmStart(theta0)

        // Build options, injecting EM warm-start when appropriate
        let opts: MLEOptimizationOpts<T> = {
            if var o = options {
                o.optimizer = optimizer
                o.computeCovariance = true
                // Only set warmStart if not provided by caller
                if o.warmStartTheta == nil {
                    o.warmStartTheta = thetaEM
                }
                o.multiStartCount = max(o.multiStartCount, max(14, phases * 3))
                o.randomRestartCount = max(o.randomRestartCount, max(4, phases))
                if o.gradStep < T(1e-6) { o.gradStep = T(1e-6) }
                if o.hessianStep < T(5e-4) { o.hessianStep = T(5e-4) }
                o.diagnosticsEnabled = true
                if o.invalidLogPdfPenalty == nil { o.invalidLogPdfPenalty = T(1e6) }
                return o
            } else {
                var o = MLEOptimizationOpts<T>()
                o.optimizer = optimizer
                o.computeCovariance = true
                // Strengthened defaults
                o.relTolLogLik = T(1e-9)
                o.tolGrad = T(1e-7)
                o.enableNewtonRefinement = true
                o.multiStartCount = max(14, phases * 3)
                o.randomRestartCount = max(4, phases)
                o.multiStartDesign = .lhs
                o.initialStepStrategy = .relativeTheta(T(0.35))
                o.gradStep = T(1e-5)
                o.hessianStep = T(5e-4)
                o.diagnosticsEnabled = true
                o.invalidLogPdfPenalty = T(1e6)
                o.warmStartTheta = thetaEM
                return o
            }
        }()

        let rawResult = MLESolver<T>.fit(problem, options: opts)
        guard let decodedHat = decodeTheta(rawResult.thetaHat) else {
            return rawResult
        }

        @inline(__always)
        func packPhysical(_ decoded: DecodedTheta) -> [T] {
            var vector = decoded.probs
            vector.append(contentsOf: decoded.rates)
            return vector
        }

        func makeJacobian(for decoded: DecodedTheta) -> [[T]] {
            var jac = Array(repeating: Array(repeating: T.zero, count: totalParams), count: totalParams)
            for i in 0..<phases {
                for j in 0..<phases {
                    let delta = (i == j) ? T.one : T.zero
                    jac[i][j] = decoded.probs[i] * (delta - decoded.probs[j])
                }
            }
            for i in 0..<phases {
                let row = phases + i
                for j in 0..<phases {
                    let col = phases + j
                    jac[row][col] = (i == j) ? decoded.rates[i] : T.zero
                }
            }
            return jac
        }

        func transpose(_ matrix: [[T]]) -> [[T]] {
            guard let first = matrix.first else { return [] }
            var result = Array(repeating: Array(repeating: T.zero, count: matrix.count), count: first.count)
            for i in 0..<matrix.count {
                for j in 0..<first.count {
                    result[j][i] = matrix[i][j]
                }
            }
            return result
        }

        func multiply(_ A: [[T]], _ B: [[T]]) -> [[T]]? {
            guard !A.isEmpty, !B.isEmpty else { return nil }
            let m = A.count
            let k = A[0].count
            guard B.count == k else { return nil }
            let n = B[0].count
            guard A.allSatisfy({ $0.count == k }), B.allSatisfy({ $0.count == n }) else { return nil }
            var result = Array(repeating: Array(repeating: T.zero, count: n), count: m)
            for i in 0..<m {
                for j in 0..<n {
                    var sum: T = .zero
                    for l in 0..<k {
                        sum += A[i][l] * B[l][j]
                    }
                    result[i][j] = sum
                }
            }
            return result
        }

        func transformCovariance(_ cov: [[T]]?, using decoded: DecodedTheta) -> [[T]]? {
            guard let cov = cov, cov.count == totalParams else { return nil }
            guard cov.allSatisfy({ $0.count == totalParams }) else { return nil }
            let jac = makeJacobian(for: decoded)
            guard let temp = multiply(jac, cov) else { return nil }
            return multiply(temp, transpose(jac))
        }

        let convertedSolutions = rawResult.allSolutions?.map { pair -> ([T], T) in
            if let decoded = decodeTheta(pair.0) {
                return (packPhysical(decoded), pair.1)
            } else {
                return pair
            }
        }

        return MLEResult(
            thetaHat: packPhysical(decodedHat),
            logLik: rawResult.logLik,
            iterations: rawResult.iterations,
            converged: rawResult.converged,
            nEval: rawResult.nEval,
            cov: transformCovariance(rawResult.cov, using: decodedHat),
            robustCov: transformCovariance(rawResult.robustCov, using: decodedHat),
            allSolutions: convertedSolutions,
            uniqueSolutionCount: rawResult.uniqueSolutionCount,
            uniqueSolutionCountTheta: rawResult.uniqueSolutionCountTheta,
            gradientNormAtOpt: rawResult.gradientNormAtOpt,
            hessianPositiveDefinite: rawResult.hessianPositiveDefinite,
            conditionNumberEstimateHu: rawResult.conditionNumberEstimateHu,
            conditionSource: rawResult.conditionSource,
            convergenceReason: rawResult.convergenceReason
        )
    }
}
