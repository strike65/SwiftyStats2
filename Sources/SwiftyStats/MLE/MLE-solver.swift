//
//  Created by VT on 07.11.25.
//  © 2025 Volker Thieme. All rights reserved.
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to do so, subject to the following conditions:
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
/// Optimization engines (Nelder–Mead, LBFGS-like, trust region) powering MLE fits.

import SwiftyStatsPrelude
#if canImport(Accelerate)
import Accelerate
#endif
import Dispatch

/// Choice of optimizer.
public enum OptimizerKind: Sendable {
    case nelderMead
    case bfgs
    case lbfgs
}

// MARK: - MLE Solver

/// Nelder–Mead or (L-)BFGS based MLE solver operating in unconstrained u-coordinates with
/// smooth transforms to constrained parameter theta-space.
internal struct MLESolver<T: RealLike> {

    /// Fit an MLE problem by maximizing the joint log-likelihood in θ-space.
    ///
    /// The solver works in unconstrained u-coordinates, applies smooth transforms to
    /// enforce bounds, and orchestrates multi-start refinement so the numerically
    /// largest log-likelihood survives. Diagnostics such as Hessian conditioning
    /// and uniqueness counts are produced when the corresponding options are enabled.
    public static func fit(
        _ problem: MLEProblem<T>,
        options: MLEOptimizationOpts<T> = .init()
    ) -> MLEResult<T> {

        let k = problem.paramSpecs.count
        var cachedNumericH_u: [[T]]? = nil
        // Build transforms and initial u0 from theta0
        let transforms = problem.paramSpecs.map { Transform<T>.make($0.constraint) }
        let theta0 = problem.paramSpecs.map { $0.initial }
        var u0 = [T](repeating: 0, count: k)
        for i in 0..<k { u0[i] = transforms[i].toU(theta0[i]) }

        // Generate candidate starts in theta and map to u
        var thetaCandidates: [[T]] = []

        // Always include base theta0 to preserve historical behavior
        thetaCandidates.append(theta0)

        // Include warm-start if provided and distinct from theta0
        if let warm = options.warmStartTheta, warm.count == k {
            let sameAsBase = zip(warm, theta0).allSatisfy { $0 == $1 }
            if !sameAsBase {
                thetaCandidates.append(warm)
            }
        }

        // Global candidates according to design
        let needCount = max(options.multiStartCount, 1)
        if needCount > thetaCandidates.count {
            let extra = needCount - thetaCandidates.count
            let generated = generateCandidates(problem.paramSpecs,
                                               baseTheta: theta0,
                                               count: extra,
                                               design: options.multiStartDesign,
                                               logScale: options.multiStartLogScale,
                                               gridPerDim: options.gridPointsPerDim,
                                               rngSeed: options.rngSeed,
                                               problem: problem,
                                               transforms: transforms,
                                               options: options)
            thetaCandidates.append(contentsOf: generated)
        }

        // Inject random restarts if requested
        if options.randomRestartCount > 0 {
            let rnd = generateCandidates(problem.paramSpecs,
                                         baseTheta: theta0,
                                         count: options.randomRestartCount,
                                         design: .random,
                                         logScale: options.multiStartLogScale,
                                         gridPerDim: options.gridPointsPerDim,
                                         rngSeed: options.rngSeed,
                                         problem: problem,
                                         transforms: transforms,
                                         options: options)
            thetaCandidates.append(contentsOf: rnd)
        }

        // Map to u-space
        var uStarts: [[T]] = thetaCandidates.map { th in
            zip(th, transforms).map { (theta, tf) in tf.toU(theta) }
        }

        // Also consider jitter around u0 for unbounded params to avoid lock-in
        if needCount == 1 && options.randomRestartCount == 0 {
            // nothing
        } else {
            let jittered = jitterAroundU(u0, problem.paramSpecs, transforms, count: 0, radius: options.multiStartURadius, rngSeed: options.rngSeed)
            uStarts.append(contentsOf: jittered)
        }

        // Local refinement for each start; keep best (highest log-likelihood)
        var allSolutions: [([T], T)] = []
        var uSolutions: [[T]] = []
        var bestTheta: [T] = theta0
        var bestLogLik: T = -.infinity
        var bestIter = 0
        var bestConverged = false
        var totalEval = 0
        var bestU: [T] = u0
        var bestGradNorm: T? = nil
        var bestInvH_u: [[T]]? = nil
        var bestReason: ConvergenceReason? = nil

        /// Reduce a single local solve into the global incumbent.
        func recordResult(_ r: LocalResult) {
            totalEval += r.nEval
            let thetaHat = zip(r.uHat, transforms).map { (ui, tf) in tf.toTheta(ui) }
            let logLik = -r.fHat
            allSolutions.append((thetaHat, logLik))
            uSolutions.append(r.uHat)
            if logLik > bestLogLik {
                bestLogLik = logLik
                bestTheta = thetaHat
                bestIter = r.iter
                bestConverged = r.converged
                bestU = r.uHat
                bestGradNorm = r.gradNorm
                bestInvH_u = r.invH_u
                bestReason = r.reason
            }
        }

        var ranParallel = false
        #if swift(>=5.5)
        if options.enableParallelStarts, uStarts.count > 1 {
            if #available(macOS 12.0, iOS 15.0, tvOS 15.0, watchOS 8.0, *) {
                let results = parallelLocalRuns(problem: problem, transforms: transforms, uStarts: uStarts, options: options)
                for r in results { recordResult(r) }
                ranParallel = true
            }
        }
        #endif
        if !ranParallel {
            for uStart in uStarts {
                let local = runLocal(problem, transforms, uStart, options)
                recordResult(local)
            }
        }

        // Diagnostics: count unique terminal solutions in u-space (relative) and θ-space (range-scaled)
        var uniqueCount: Int? = nil
        var uniqueCountTheta: Int? = nil
        if options.diagnosticsEnabled {
            uniqueCount = countUniqueUSolutions(uSolutions, tolRel: options.uniquenessRelTolU)
            let thetaSolutions = allSolutions.map { $0.0 }
            uniqueCountTheta = countUniqueThetaSolutions(thetaSolutions,
                                                        specs: problem.paramSpecs,
                                                        tol: options.uniquenessTolTheta)
        }

        // Optional Hessian PD check at optimum (symmetrized)
        var hessianPD: Bool? = nil
        var conditionEstimate: T? = nil
        var conditionSource: String? = nil
        if options.diagnosticsEnabled {
            if let invH = bestInvH_u {
                let invSym = symmetrize(invH)
                hessianPD = isSymPosDef(invSym)
                if let estimate = estimateConditionNumberGershgorin(invSym) {
                    conditionEstimate = estimate
                    conditionSource = "H_u^{-1}-Gershgorin"
                }
            } else {
                if cachedNumericH_u == nil {
                    var H_u = numericHessian(problem, transforms, bestU, step: options.hessianStep, options: options)
                    H_u = symmetrize(H_u)
                    cachedNumericH_u = H_u
                }
                if let Hu = cachedNumericH_u {
                    hessianPD = isSymPosDef(Hu)
                    if let estimate = estimateConditionNumberGershgorin(Hu) {
                        conditionEstimate = estimate
                        conditionSource = "H_u-Gershgorin"
                    }
                }
            }
        }

        // Optional covariance
        var cov: [[T]]? = nil
        var robustCov: [[T]]? = nil
        if options.computeCovariance {
            // Prefer inverse Hessian from optimizer if available and PD
            var Cov_u: [[T]]? = nil
            if let invH = bestInvH_u {
                let Hsym = symmetrize(invH)
                if isSymPosDef(Hsym) {
                    Cov_u = Hsym
                }
            }
            // Fallback to numeric Hessian inversion if unavailable or not PD
            if Cov_u == nil {
                if cachedNumericH_u == nil {
                    var H_u = numericHessian(problem, transforms, bestU, step: options.hessianStep, options: options)
                    H_u = symmetrize(H_u)
                    cachedNumericH_u = H_u
                }
                if let Hu = cachedNumericH_u,
                   let Cu = invertSymmetricCholesky(Hu,
                                                    regInit: options.hessianRegInitial,
                                                    regFactor: options.hessianRegFactor,
                                                    regMax: options.hessianRegMax) ?? invertSymmetric(Hu) {
                    Cov_u = Cu
                }
            } else {
                // If we already have Cov_u from optimizer, we still need H_u for robust covariance scaling if we want to verify PD etc.
                // Not strictly necessary; robust covariance uses H^{-1} directly. We can skip H_u here.
            }

            // Map Cov_u to θ-space with Jacobian
            if let Cov_u {
                let jacDiag: [T] = zip(bestU, transforms).map { (u, tr) in tr.dToTheta(u) }
                var Cov_theta = Cov_u
                let k = jacDiag.count
                for i in 0..<k {
                    let di = jacDiag[i]
                    for j in 0..<k { Cov_theta[i][j] *= di }
                }
                for j in 0..<k {
                    let dj = jacDiag[j]
                    for i in 0..<k { Cov_theta[i][j] *= dj }
                }
                cov = Cov_theta
            }

            // Robust “sandwich” covariance if per-observation gradient is available
            if problem.gradlogpdf != nil {
                // Need H_u^{-1} (already Cov_u) and the score outer products in u-space
                if let Hinv_u = Cov_u {
                    // Compute J_u = sum_x s_u(x) s_u(x)^T, where s_u(x) = -diag(dθ/du) * gradlogpdf(x, θ̂)
                    let thetaHat = zip(bestU, transforms).map { (ui, tf) in tf.toTheta(ui) }
                    let jacDiagU: [T] = zip(bestU, transforms).map { (u, tr) in tr.dToTheta(u) }
                    var J_u = Array(repeating: Array(repeating: T.zero, count: k), count: k)
                    for x in problem.data {
                        if let gtheta = problem.gradlogpdf?(x, thetaHat), gtheta.count == k,
                           gtheta.allSatisfy({ $0.isFinite }) {
                            var su = [T](repeating: 0, count: k)
                            for i in 0..<k {
                                su[i] = -(jacDiagU[i] * gtheta[i])
                            }
                            // Outer product su su^T and accumulate
                            for i in 0..<k {
                                let sui = su[i]
                                for j in 0..<k {
                                    J_u[i][j] += sui * su[j]
                                }
                            }
                        } else {
                            // If any gradient is invalid, skip robust covariance
                            J_u = []
                            break
                        }
                    }
                    if !J_u.isEmpty {
                        // Cov_u,robust = Hinv_u * J_u * Hinv_u
                        let temp = matmul(Hinv_u, J_u)
                        let Cov_u_robust = matmul(temp, Hinv_u)
                        // Map to θ-space: Cov_θ,robust = G * Cov_u_robust * G, where G = diag(dθ/du at û)
                        let Gdiag = jacDiagU
                        var Cov_theta_robust = Cov_u_robust
                        for i in 0..<k {
                            let di = Gdiag[i]
                            for j in 0..<k { Cov_theta_robust[i][j] *= di }
                        }
                        for j in 0..<k {
                            let dj = Gdiag[j]
                            for i in 0..<k { Cov_theta_robust[i][j] *= dj }
                        }
                        robustCov = Cov_theta_robust
                    }
                } else {
                    // If we couldn’t form H^{-1}, robust covariance cannot be formed reliably.
                    robustCov = nil
                }
            }
        }

        return MLEResult(thetaHat: bestTheta,
                         logLik: bestLogLik,
                         iterations: bestIter,
                         converged: bestConverged,
                         nEval: totalEval,
                         cov: cov,
                         robustCov: robustCov,
                         allSolutions: options.diagnosticsEnabled ? allSolutions : nil,
                         uniqueSolutionCount: uniqueCount,
                         uniqueSolutionCountTheta: uniqueCountTheta,
                         gradientNormAtOpt: bestGradNorm,
                         hessianPositiveDefinite: hessianPD,
                         conditionNumberEstimateHu: conditionEstimate,
                         conditionSource: conditionSource,
                         convergenceReason: bestReason)
    }

    /// Bundle of statistics returned by each local optimizer call.
    private typealias LocalResult = (uHat: [T], fHat: T, iter: Int, converged: Bool, nEval: Int, gradNorm: T?, invH_u: [[T]]?, reason: ConvergenceReason?)

    // MARK: - Local refinement wrapper

    private static func runLocal(
        _ problem: MLEProblem<T>,
        _ transforms: [Transform<T>],
        _ u0: [T],
        _ options: MLEOptimizationOpts<T>
    ) -> LocalResult {

        // Choose optimizer
        var core: (uHat: [T], fHat: T, iter: Int, converged: Bool, nEval: Int, invH_u: [[T]]?, reason: ConvergenceReason?)
        switch options.optimizer {
        case .nelderMead:
            let nm = nelderMead(problem, transforms, u0, options)
            core = (nm.uHat, nm.fHat, nm.iter, nm.converged, nm.nEval, nil, nm.reason)
        case .bfgs:
            if let r = bfgsLike(problem, transforms, u0, options, useLBFGS: false) {
                core = r
            } else {
                let nm = nelderMead(problem, transforms, u0, options)
                core = (nm.uHat, nm.fHat, nm.iter, nm.converged, nm.nEval, nil, nm.reason)
            }
        case .lbfgs:
            if let r = bfgsLike(problem, transforms, u0, options, useLBFGS: true) {
                core = r
            } else {
                let nm = nelderMead(problem, transforms, u0, options)
                core = (nm.uHat, nm.fHat, nm.iter, nm.converged, nm.nEval, nil, nm.reason)
            }
        }

        // Optional trust-region refinement if gradient small
        var gradNormOut: T? = nil
        if options.enableNewtonRefinement {
            if let (g, gEvals) = gradNLL_u(problem, transforms, core.uHat, options) {
                let gNorm = norm2(g)
                gradNormOut = gNorm
                // If gradient small, try trust-region refinement
                if gNorm < max(options.tolGrad, T(10) * .leastNonzeroMagnitude) {
                    if let tr = trustRegionRefinement(problem, transforms, core.uHat, core.fHat, g, options) {
                        // Count the gradient evals we already did to compute g at core.uHat
                        let adjustedEval = core.nEval + gEvals + tr.nEval
                        // Accept if non-worse
                        if tr.fHat <= core.fHat {
                            core = (tr.uHat, tr.fHat, core.iter, true, adjustedEval, core.invH_u, core.reason)
                        } else {
                            core = (core.uHat, core.fHat, core.iter, core.converged, adjustedEval, core.invH_u, core.reason)
                        }
                    } else {
                        core.nEval += gEvals
                    }
                } else {
                    core.nEval += gEvals
                }
            }
        }

        return (core.uHat, core.fHat, core.iter, core.converged, core.nEval, gradNormOut, core.invH_u, core.reason)
    }

#if swift(>=5.5)
    /// Evaluate multiple starts concurrently (Apple platforms) so the wall-clock
    /// cost of multi-start search scales sub-linearly with the number of starts.
    @available(macOS 12.0, iOS 15.0, tvOS 15.0, watchOS 8.0, *)
    private static func parallelLocalRuns(
        problem: MLEProblem<T>,
        transforms: [Transform<T>],
        uStarts: [[T]],
        options: MLEOptimizationOpts<T>
    ) -> [LocalResult] {
#if swift(>=5.5)
        runBlocking(priority: options.parallelTaskPriority) {
            await parallelLocalRunsAsync(problem: problem,
                                         transforms: transforms,
                                         uStarts: uStarts,
                                         options: options)
        }
#else
        runBlocking {
            await parallelLocalRunsAsync(problem: problem,
                                         transforms: transforms,
                                         uStarts: uStarts,
                                         options: options)
        }
#endif
    }

    @available(macOS 12.0, iOS 15.0, tvOS 15.0, watchOS 8.0, *)
    private static func parallelLocalRunsAsync(
        problem: MLEProblem<T>,
        transforms: [Transform<T>],
        uStarts: [[T]],
        options: MLEOptimizationOpts<T>
    ) async -> [LocalResult] {
        await withTaskGroup(of: (Int, LocalResult).self) { group in
            for (idx, start) in uStarts.enumerated() {
                group.addTask {
                    let result = runLocal(problem, transforms, start, options)
                    return (idx, result)
                }
            }
            var tmp = Array(
                repeating: LocalResult(uHat: [], fHat: T.zero, iter: 0, converged: false, nEval: 0, gradNorm: nil, invH_u: nil, reason: nil),
                count: uStarts.count
            )
            for await (idx, value) in group {
                tmp[idx] = value
            }
            return tmp
        }
    }
#endif

#if swift(>=5.5)
    private static func runBlocking<R>(
        priority: TaskPriority? = nil,
        _ operation: @escaping @Sendable () async -> R
    ) -> R {
        let semaphore = DispatchSemaphore(value: 0)
        let box = BlockingResultBox<R>()
        let taskBody = {
            let value = await operation()
            box.value = value
            semaphore.signal()
        }
        if let priority = priority {
            Task.detached(priority: priority, operation: taskBody)
        } else {
            Task.detached(operation: taskBody)
        }
        semaphore.wait()
        return box.value!
    }
#else
    private static func runBlocking<R>(
        _ operation: @escaping () -> R
    ) -> R {
        return operation()
    }
#endif

    // MARK: - Nelder–Mead

    /// Nelder–Mead with θ-aware step control and additional stall criteria.
    private static func nelderMead(
        _ problem: MLEProblem<T>,
        _ transforms: [Transform<T>],
        _ u0: [T],
        _ options: MLEOptimizationOpts<T>
    ) -> (uHat: [T], fHat: T, iter: Int, converged: Bool, nEval: Int, reason: ConvergenceReason?) {
        let k = u0.count

        // Initial simplex in u-space with per-dimension step selection
        var simplex = [[T]]()
        simplex.append(u0)

        // Build theta0 and dtheta/du at u0 for relative step strategy
        let theta0 = zip(u0, transforms).map { (ui, tf) in tf.toTheta(ui) }
        let dtheta_du0 = zip(u0, transforms).map { (ui, tf) in tf.dToTheta(ui) }

        for j in 0..<k {
            var uj = u0
            let du: T
            switch options.initialStepStrategy {
            case .spec:
                du = problem.paramSpecs[j].step
            case .absoluteU(let s):
                du = s
            case .relativeTheta(let r):
                let scaleTheta = r * T.maximum(T.one, abs(theta0[j]))
                let deriv = dtheta_du0[j]
                let safeDeriv = (abs(deriv) > T(1e-12)) ? deriv : (deriv.sign == .minus ? -T(1e-12) : T(1e-12))
                du = scaleTheta / safeDeriv
            }
            uj[j] += du
            simplex.append(uj)
        }

        var fvals = simplex.map { nll(problem, transforms, $0, options) }
        var nEval = fvals.count

        func order() {
            let zipped = zip(fvals, simplex).sorted { $0.0 < $1.0 }
            for i in 0..<zipped.count {
                fvals[i] = zipped[i].0
                simplex[i] = zipped[i].1
            }
        }

        order()
        var iter = 0
        var converged = false
        var prevBest = fvals[0]
        var stallCount = 0
        var reason: ConvergenceReason? = nil

        // New: track best vertex movement in θ-space and a window of relative improvements
        var bestU_prev = simplex[0]
        var bestTheta_prev = zip(bestU_prev, transforms).map { (u, tr) in tr.toTheta(u) }
        var stallWindowBuf: [T] = []
        stallWindowBuf.reserveCapacity(max(1, options.stallWindow))

        while iter < options.maxIter {
            iter += 1

            // Stopping criteria
            let fmin = fvals.min()!
            let fmax = fvals.max()!
            let fSpreadAbs = fmax - fmin
            let fSpreadRel = fSpreadAbs / T.maximum(T.one, abs(fmin))
            let size: T = simplexDiameter(simplex)
            let relChange = relChangeVal(prevBest, fvals[0])

            // Maintain windowed stall buffer on best f improvements
            if options.stallWindow > 0 {
                if stallWindowBuf.count == options.stallWindow {
                    stallWindowBuf.removeFirst()
                }
                stallWindowBuf.append(relChange)
            }

            // Stall counter on best-f relative improvement
            if relChange < options.stallRelTolNM {
                stallCount += 1
            } else {
                stallCount = 0
            }

            // θ-step tolerance based on movement of best vertex
            let bestU_curr = simplex[0]
            let jacPrev: [T] = zip(bestU_prev, transforms).map { (u, tr) in tr.dToTheta(u) }
            let dUbest = zip(bestU_curr, bestU_prev).map { $0 - $1 }
            let dThetaApprox = zip(jacPrev, dUbest).map { $0 * $1 }
            let thetaNorm = T.sqrt(bestTheta_prev.map { $0 * $0 }.reduce(0, +))
            let dThetaNorm = T.sqrt(dThetaApprox.map { $0 * $0 }.reduce(0, +))
            let thetaStepRel = dThetaNorm / T.maximum(T.one, thetaNorm)

            // Windowed median stall decision
            var stallMedianStop = false
            if options.stallWindow > 0, stallWindowBuf.count == options.stallWindow {
                var tmp = stallWindowBuf
                tmp.sort()
                let mid = options.stallWindow / 2
                let median: T = (options.stallWindow % 2 == 1) ? tmp[mid] : (tmp[mid-1] + tmp[mid]) / 2
                stallMedianStop = median < options.stallMedianTol
            }

            let simplexCollapse = ((fSpreadAbs < options.tolFNM || fSpreadRel < options.relTolFSpreadNM) && size < options.tolSimplexNM)
            let relStop = relChange < options.relTolLogLik
            let stallIterStop = stallCount >= options.stallIterationsNM
            let stepStop = thetaStepRel < options.tolStepTheta
            if simplexCollapse || relStop || stallIterStop || stepStop || stallMedianStop {
                converged = true
                if simplexCollapse || relStop {
                    reason = .function
                } else if stepStop {
                    reason = .step
                } else if stallIterStop || stallMedianStop {
                    reason = .stall
                }
                break
            }

            // Centroid of the best n points (excluding the worst)
            let n = k
            var centroid = [T](repeating: 0, count: k)
            for i in 0..<n {
                for j in 0..<k { centroid[j] += simplex[i][j] }
            }
            for j in 0..<k { centroid[j] /= T(n) }

            let worst = k
            let bestVal = fvals[0]
            let worstVal = fvals[worst]

            // Reflection
            var xr = [T](repeating: 0, count: k)
            for j in 0..<k { xr[j] = centroid[j] + options.alpha * (centroid[j] - simplex[worst][j]) }
            let fr = nll(problem, transforms, xr, options); nEval += 1

            if fr < fvals[0] {
                // Expansion
                var xe = [T](repeating: 0, count: k)
                for j in 0..<k { xe[j] = centroid[j] + options.gamma * (xr[j] - centroid[j]) }
                let fe = nll(problem, transforms, xe, options); nEval += 1
                if fe < fr {
                    simplex[worst] = xe; fvals[worst] = fe
                } else {
                    simplex[worst] = xr; fvals[worst] = fr
                }
            } else if fr < fvals[k-1] {
                // Accept reflection
                simplex[worst] = xr; fvals[worst] = fr
            } else {
                // Contraction
                let xc: [T]
                if fr < worstVal {
                    // Outside contraction
                    xc = zip(centroid, xr).map { (c, r) in c + options.rho * (r - c) }
                } else {
                    // Inside contraction
                    xc = zip(centroid, simplex[worst]).map { (c, w) in c + options.rho * (w - c) }
                }
                let fc = nll(problem, transforms, xc, options); nEval += 1
                if fc < min(fr, worstVal) {
                    simplex[worst] = xc; fvals[worst] = fc
                } else {
                    // Shrink
                    for i in 1...k {
                        for j in 0..<k {
                            simplex[i][j] = simplex[0][j] + options.sigma * (simplex[i][j] - simplex[0][j])
                        }
                        fvals[i] = nll(problem, transforms, simplex[i], options); nEval += 1
                    }
                }
            }

            order()
            prevBest = bestVal

            // Update best-vertex θ tracking for θ-step tolerance
            if simplex[0] != bestU_prev {
                bestU_prev = simplex[0]
                bestTheta_prev = zip(bestU_prev, transforms).map { (u, tr) in tr.toTheta(u) }
            }
        }

        if !converged {
            reason = .maxIter
        }
        return (simplex[0], fvals[0], iter, converged, nEval, reason)
    }

    /// Smooth, scale-aware boundary penalty in u-space.
    /// Uses the constraint shape and the local derivative dθ/du at u to build a gentle repulsion near boundaries.
    /// Returns a nonnegative penalty to be added to NLL.
    @inline(__always)
    private static func boundaryPenalty(
        _ transforms: [Transform<T>],
        _ specs: [ParamSpec<T>],
        _ u: [T],
        weight: T,
        eps: T,
        power: T
    ) -> T {
        if weight == 0 { return 0 }
        var pen: T = 0
        for i in 0..<u.count {
            let tr = transforms[i]
            let ui = u[i]
            let d = tr.dToTheta(ui)
            let di = abs(d) + eps // avoid zero
            let theta = tr.toTheta(ui)

            switch specs[i].constraint {
            case .real:
                // no boundaries, no penalty
                break
            case .positive:
                // distance to boundary θ=0 measured via soft barrier on small θ
                // u-scale: divide by |dθ/du| to normalize
                let t = max(theta, eps)
                pen += weight * T.pow((T.one / t) / di, power)
            case .lowerBound(let a):
                let dist = max(theta - a, eps)
                pen += weight * T.pow((T.one / dist) / di, power)
            case .upperBound(let b):
                let dist = max(b - theta, eps)
                pen += weight * T.pow((T.one / dist) / di, power)
            case .unitInterval:
                let t = min(max(theta, eps), 1 - eps)
                let distLo = max(t - 0, eps)
                let distHi = max(1 - t, eps)
                pen += weight * (T.pow((T.one / distLo) / di, power) + T.pow((T.one / distHi) / di, power))
            case .interval(let a, let b):
                let t = min(max(theta, a + eps), b - eps)
                let distLo = max(t - a, eps)
                let distHi = max(b - t, eps)
                pen += weight * (T.pow((T.one / distLo) / di, power) + T.pow((T.one / distHi) / di, power))
            }
        }
        return pen
    }

    /// Negative log-likelihood in u-coordinates (with stabilization/penalty).
    private static func nll(
        _ problem: MLEProblem<T>,
        _ transforms: [Transform<T>],
        _ u: [T],
        _ options: MLEOptimizationOpts<T>
    ) -> T {
        // θ(u)
        let theta = zip(u, transforms).map { (ui, tf) in tf.toTheta(ui) }

        // Accumulate −sum logpdf(x, θ) using compensated summation
        var terms = [T]()
        terms.reserveCapacity(problem.data.count)
        for x in problem.data {
            let lp = problem.logpdf(x, theta)
            if !lp.isFinite || lp.isNaN || lp == T.infinity || lp == -T.infinity {
                if let pen = options.invalidLogPdfPenalty {
                    // finite penalty for invalid evaluations; treat as adding a single penalty term
                    terms.append(pen)
                    continue
                } else {
                    // strong penalty for invalid evaluations
                    return T.infinity
                }
            }
            terms.append(-lp)
        }
        var tmp = terms
        var acc = Helpers.sum(&tmp)
        if !acc.isFinite { return T.infinity }

        // Add smooth, scale-aware boundary penalty in u-space (optional)
        if options.boundaryPenaltyWeight > 0 {
            let pen = boundaryPenalty(transforms, problem.paramSpecs, u,
                                      weight: options.boundaryPenaltyWeight,
                                      eps: options.boundaryPenaltyEps,
                                      power: options.boundaryPenaltyPower)
            acc += pen
        }
        return acc
    }

    // MARK: - Gradient-based (BFGS / L-BFGS)

    /// BFGS / L-BFGS with Powell damping, curvature tests, and strong Wolfe line search.
    private static func bfgsLike(
        _ problem: MLEProblem<T>,
        _ transforms: [Transform<T>],
        _ u0: [T],
        _ options: MLEOptimizationOpts<T>,
        useLBFGS: Bool
    ) -> (uHat: [T], fHat: T, iter: Int, converged: Bool, nEval: Int, invH_u: [[T]]?, reason: ConvergenceReason?)? {
        let k = u0.count
        var u = u0
        var f = nll(problem, transforms, u, options)
        var nEval = 1

        // Initialize gradient
        guard let g0 = gradNLL_u(problem, transforms, u, options) else { return nil }
        nEval += g0.1
        var grad = g0.0
        var gradNorm = norm2(grad)
        if !gradNorm.isFinite { return nil }

        // BFGS state
        var H = eye(k) // inverse Hessian approximation
        // L-BFGS history
        var sList = [[T]]()
        var yList = [[T]]()
        var rhoList = [T]()

        var iter = 0
        var converged = false
        var prevF = f
        var reason: ConvergenceReason? = nil

        // New: windowed stall buffer for relative NLL improvements
        var stallWindowBuf: [T] = []
        stallWindowBuf.reserveCapacity(max(1, options.stallWindow))

        while iter < options.maxIter {
            iter += 1

            // θ-scaled gradient infinity norm
            let thetaHere = zip(u, transforms).map { (ui, tf) in tf.toTheta(ui) }
            var gScaledInf: T = 0
            for i in 0..<k {
                let denom = T.maximum(T.one, abs(thetaHere[i]))
                gScaledInf = T.maximum(gScaledInf, abs(grad[i]) / denom)
            }
            if gScaledInf < options.tolGradInfTheta || gradNorm < options.tolGrad {
                converged = true
                reason = .gradient
                break
            }

            if iter > 1 {
                // Relative log-likelihood change via NLL
                let rel = relChangeVal(prevF, f)

                // Update windowed stall buffer and check median-based stall
                if options.stallWindow > 0 {
                    if stallWindowBuf.count == options.stallWindow { stallWindowBuf.removeFirst() }
                    stallWindowBuf.append(rel)
                    if stallWindowBuf.count == options.stallWindow {
                        var tmp = stallWindowBuf
                        tmp.sort()
                        let mid = options.stallWindow / 2
                        let median: T = (options.stallWindow % 2 == 1) ? tmp[mid] : (tmp[mid-1] + tmp[mid]) / 2
                        if median < options.stallMedianTol {
                            converged = true
                            reason = .stall
                            break
                        }
                    }
                }

                if rel < options.relTolLogLik {
                    converged = true
                    reason = .function
                    break
                }
            }

            // Compute search direction p = -H g (BFGS) or two-loop recursion (L-BFGS)
            let p: [T]
            if useLBFGS && !sList.isEmpty {
                p = lbfgsDirectionWithGammaClamped(grad, sList, yList, rhoList, floor: options.h0ScaleFloor, ceil: options.h0ScaleCeil)
            } else {
                p = matvec(H, grad).map { -$0 }
            }

            // Ensure descent; if not, reset
            let descent = dot(p, grad)
            var pDir = p
            if !(descent < 0 && descent.isFinite) {
                // reset
                H = eye(k)
                pDir = grad.map { -$0 }
            }

            // Strong Wolfe line search on NLL with soft backoff on invalid evaluations
            let c1 = options.lineSearchC1
            let c2 = options.lineSearchC2
            var (step, fNew, gNew, lsEvals, ok) = strongWolfeLineSearch(problem, transforms, u, f, grad, pDir, c1: c1, c2: c2, options: options)
            nEval += lsEvals

            // If failed, try a soft invalid backoff loop (shrink step until valid)
            if !ok {
                var alpha = T.one
                var attempts = 0
                let tau = T.maximum(T(1e-6), options.invalidBackoffTau)
                var bestF = f
                var bestG = grad
                var bestU = u
                while attempts < options.invalidBackoffMax {
                    attempts += 1
                    alpha *= tau
                    let uTry = zip(u, pDir).map { $0 + alpha * $1 }
                    let fTry = nll(problem, transforms, uTry, options); nEval += 1
                    if fTry.isFinite, let gPair = gradNLL_u(problem, transforms, uTry, options) {
                        let gTry = gPair.0
                        step = alpha
                        fNew = fTry
                        gNew = gTry
                        ok = true
                        bestF = fTry; bestG = gTry; bestU = uTry
                        break
                    }
                }
                if ok == false {
                    // If still not ok, give up this local run
                    return nil
                } else {
                    // Apply the backoff result
                    // θ-step tolerance check using Jacobian at current u and step s
                    let sVec = zip(bestU, u).map { $0 - $1 }
                    let jacHere: [T] = zip(u, transforms).map { (ui, tr) in tr.dToTheta(ui) }
                    let dTheta = zip(jacHere, sVec).map { $0 * $1 }
                    let thetaNorm = T.sqrt(thetaHere.map { $0 * $0 }.reduce(0, +))
                    let dThetaNorm = T.sqrt(dTheta.map { $0 * $0 }.reduce(0, +))
                    let thetaStepRel = dThetaNorm / T.maximum(T.one, thetaNorm)
                    if thetaStepRel < options.tolStepTheta {
                        converged = true
                        u = bestU
                        f = bestF
                        grad = bestG
                        gradNorm = norm2(grad)
                        break
                    }

                    u = bestU
                    prevF = f
                    f = bestF
                    grad = bestG
                    gradNorm = norm2(grad)
                    // No BFGS update possible without s,y; continue to next iter
                    continue
                }
            }

            // Convergence on step in θ-space using Jacobian at current u
            let uNew = zip(u, pDir).map { (ui, pi) in ui + step * pi }
            let s = zip(uNew, u).map { (a, b) in a - b }
            let jacHere: [T] = zip(u, transforms).map { (ui, tr) in tr.dToTheta(ui) }
            let dTheta = zip(jacHere, s).map { $0 * $1 }
            let thetaNorm = T.sqrt(thetaHere.map { $0 * $0 }.reduce(0, +))
            let dThetaNorm = T.sqrt(dTheta.map { $0 * $0 }.reduce(0, +))
            let thetaStepRel = dThetaNorm / T.maximum(T.one, thetaNorm)
            if thetaStepRel < options.tolStepTheta {
                prevF = f
                u = uNew
                f = fNew
                grad = gNew
                gradNorm = norm2(grad)
                converged = true
                break
            }

            // BFGS/L-BFGS update with curvature check and Powell damping
            let sVec = s   // s_k = u_{k+1} - u_k
            var y = zip(gNew, grad).map { (a, b) in a - b } // y_k = g_{k+1} - g_k
            var ys = dot(y, sVec)
            let sn = norm2(sVec)
            let yn = norm2(y)
            let curvThresh = options.curvatureEps * sn * yn

            // If curvature too small or negative, apply cautious update:
            if !(ys > curvThresh && ys.isFinite) {
                // Attempt Powell damping
                let thetaPow = T.minimum(T.one, T.maximum(T.zero, options.powellTheta))
                let Hy = useLBFGS ? sVec : matvec(H, sVec)
                let yBar = zip(y, Hy).map { thetaPow * $0 + (T.one - thetaPow) * $1 }
                let ysBar = dot(yBar, sVec)
                let ynBar = norm2(yBar)

                if ysBar > curvThresh && ysBar.isFinite && ynBar.isFinite {
                    // Accept damped pair
                    y = yBar
                    ys = ysBar
                } else {
                    if useLBFGS {
                        // Skip update
                    } else {
                        // Reset H to identity to recover
                        H = eye(k)
                    }
                    // Move state and continue
                    prevF = f
                    u = uNew
                    f = fNew
                    grad = gNew
                    gradNorm = norm2(grad)
                    let stepNorm = norm2(sVec)
                    if stepNorm < options.tolStep {
                        converged = true
                        reason = .step
                        break
                    }
                    continue
                }
            }

            // Safe update now that curvature is adequate
            if useLBFGS {
                let rho = T.one / ys
                sList.append(sVec)
                yList.append(y)
                rhoList.append(rho)
                if sList.count > options.lbfgsMemory {
                    sList.removeFirst()
                    yList.removeFirst()
                    rhoList.removeFirst()
                }
            } else {
                // Full BFGS update: H_{k+1} = (I - ρ s y^T) H (I - ρ y s^T) + ρ s s^T
                let rho = T.one / ys
                let I = eye(k)
                let syT = outer(sVec, y)
                let ysT = outer(y, sVec)
                let A = add(I, scale(syT, -rho))
                let B = add(I, scale(ysT, -rho))
                let AHB = matmul(A, matmul(H, B))
                let ssT = outer(sVec, sVec)
                H = add(AHB, scale(ssT, rho))
            }

            // Move
            prevF = f
            u = uNew
            f = fNew
            grad = gNew
            gradNorm = norm2(grad)

            // Step tolerance (legacy u-based as a secondary guard)
            let stepNorm = norm2(sVec)
            if stepNorm < options.tolStep {
                converged = true
                reason = .step
                break
            }
        }

        // Prepare inverse Hessian approximation to return (if possible)
        var invH_u: [[T]]? = nil
        if useLBFGS {
            invH_u = reconstructDenseLBFGSInverseH(sList: sList, yList: yList, rhoList: rhoList, dim: k, floor: options.h0ScaleFloor, ceil: options.h0ScaleCeil)
        } else {
            // Symmetrize H for safety
            invH_u = symmetrize(H)
        }

        if !converged {
            reason = .maxIter
        }
        return (u, f, iter, converged, nEval, invH_u, reason)
    }

    /// Strong Wolfe line search (NLL) with safeguarded interpolation and Armijo fallback.
    /// 1) f(u+αp) ≤ f(u) + c1 α gᵀp (sufficient decrease)
    /// 2) |∇f(u+αp)ᵀp| ≤ c2 |gᵀp| (curvature)
    private static func strongWolfeLineSearch(
        _ problem: MLEProblem<T>,
        _ transforms: [Transform<T>],
        _ u: [T],
        _ f0: T,
        _ g0: [T],
        _ p: [T],
        c1: T,
        c2: T,
        options: MLEOptimizationOpts<T>
    ) -> (alpha: T, fNew: T, gNew: [T], evals: Int, ok: Bool) {
        // Directional derivative at 0
        let phi0 = f0
        let dphi0 = dot(g0, p) // should be negative for descent
        var evals = 0

        // Guard: if not descent, fail quickly
        if !(dphi0 < 0 && dphi0.isFinite) {
            return (T.zero, f0, g0, evals, false)
        }

        // Bracketing phase
        var alphaPrev: T = 0
        var phiPrev: T = phi0

        // Minimal expansion to avoid stalling
        let minExpand: T = T(1.2)
        var alpha: T = T.one
        // Allow expansion up to stepMax; if f is ill-scaled, back off
        let stepMax = options.lineSearchStepMax
        let baseStepMin = options.lineSearchStepMin
        let stepMinScaled = max(baseStepMin, baseStepMin * norm2(p)) // scale by ||p||

        // Helper to evaluate f and g at u + a p, with invalid backoff
        func evalAt(_ a: T) -> (f: T, g: [T])? {
            var aTry = a
            var attempts = 0
            let tauInv = T.maximum(T(1e-6), options.invalidBackoffTau)
            while attempts <= options.invalidBackoffMax {
                let uA = zip(u, p).map { (ui, pi) in ui + aTry * pi }
                let fA = nll(problem, transforms, uA, options); evals += 1
                if fA.isFinite, let gPair = gradNLL_u(problem, transforms, uA, options) {
                    return (fA, gPair.0)
                }
                // shrink and retry
                aTry *= tauInv
                attempts += 1
                if aTry < stepMinScaled { break }
            }
            return nil
        }

        // Try a few bracket iterations
        var phiAlpha: T = phi0
        var gAlpha: [T] = g0
        for _ in 0..<options.lineSearchMaxIter {
            // Ensure alpha within bounds
            if alpha < stepMinScaled { alpha = stepMinScaled }
            if alpha > stepMax { alpha = stepMax }

            guard let (fA, gA) = evalAt(alpha) else {
                // If invalid after backoffs, reduce alpha and continue
                alpha *= max(T(0.1), options.invalidBackoffTau)
                continue
            }
            phiAlpha = fA
            gAlpha = gA

            // Armijo check
            if phiAlpha > phi0 + c1 * alpha * dphi0 || (phiAlpha >= phiPrev && alpha > T.zero) {
                // Zoom between [alphaPrev, alpha]
                let zoom = wolfeZoom(problem, transforms, u, p, phi0, dphi0, alphaPrev, alpha, phiPrev, phiAlpha, c1, c2, options, stepMinScaled: stepMinScaled)
                evals += zoom.evals
                if let z = zoom.res {
                    return (z.alpha, z.fNew, z.gNew, evals, true)
                } else {
                    // Final Armijo backtracking if zoom failed
                    let back = finalArmijoBacktrack(problem, transforms, u, p, phi0, dphi0, c1, options, stepMinScaled: stepMinScaled)
                    evals += back.evals
                    if let r = back.res { return (r.alpha, r.fNew, r.gNew, evals, true) }
                    return (alpha, phiAlpha, gAlpha, evals, false)
                }
            }

            // Curvature condition (with safeguard: accept if derivative is sufficiently small relative to |dphi0|)
            let dphiAlpha = dot(gAlpha, p)
            if abs(dphiAlpha) <= c2 * abs(dphi0) {
                return (alpha, phiAlpha, gAlpha, evals, true)
            }

            // If derivative positive or noisy sign flip, minimum lies between alphaPrev and alpha: zoom
            if dphiAlpha >= 0 || !dphiAlpha.isFinite {
                let zoom = wolfeZoom(problem, transforms, u, p, phi0, dphi0, alpha, alphaPrev, phiAlpha, phiPrev, c1, c2, options, stepMinScaled: stepMinScaled)
                evals += zoom.evals
                if let z = zoom.res {
                    return (z.alpha, z.fNew, z.gNew, evals, true)
                } else {
                    // Final Armijo backtracking if zoom failed
                    let back = finalArmijoBacktrack(problem, transforms, u, p, phi0, dphi0, c1, options, stepMinScaled: stepMinScaled)
                    evals += back.evals
                    if let r = back.res { return (r.alpha, r.fNew, r.gNew, evals, true) }
                    return (alpha, phiAlpha, gAlpha, evals, false)
                }
            }

            // Otherwise, increase alpha (expand) with minimal factor
            alphaPrev = alpha
            phiPrev = phiAlpha
            alpha = min(stepMax, max(alpha * T(2), alpha * minExpand))
        }

        // Fallback: cautious backtracking if bracketing failed
        let back = finalArmijoBacktrack(problem, transforms, u, p, phi0, dphi0, c1, options, stepMinScaled: stepMinScaled, extraIters: 10)
        evals += back.evals
        if let r = back.res { return (r.alpha, r.fNew, r.gNew, evals, true) }
        return (alpha, phiAlpha, gAlpha, evals, false)
    }

    // Safeguarded interpolation helpers (cubic/quadratic with bounds)
    @inline(__always) private static func safeguardedInterp(
        alo: T, phiLo: T, dphiLo: T,
        ahi: T, phiHi: T, dphiHi: T
    ) -> T? {
        // Try cubic interpolation first if derivatives on both ends are available and finite
        if dphiLo.isFinite && dphiHi.isFinite {
            let d1 = dphiLo + dphiHi - T(3) * (phiLo - phiHi) / (alo - ahi)
            let rad = d1 * d1 - dphiLo * dphiHi
            if rad.isFinite, rad >= 0 {
                let denom = dphiHi - dphiLo + T(2) * d1
                if denom != 0 {
                    let z = (dphiHi + d1 - T.sqrt(rad)) / denom
                    // Candidate in [0,1] for alpha = alo + z(ahi - alo)
                    if z.isFinite {
                        let a = alo + z * (ahi - alo)
                        if a > min(alo, ahi) && a < max(alo, ahi) {
                            return a
                        }
                    }
                }
            }
        }
        // Quadratic interpolation using phi and derivative at lower end
        if dphiLo.isFinite {
            let a = alo
            let b = ahi
            let denom = (phiHi - phiLo - dphiLo * (b - a))
            if denom != 0 {
                let z = -dphiLo * (b - a) * (b - a) / (T(2) * denom)
                let cand = a + z
                if cand.isFinite, cand > min(a, b) && cand < max(a, b) {
                    return cand
                }
            }
        }
        return nil
    }

    private static func wolfeZoom(
        _ problem: MLEProblem<T>,
        _ transforms: [Transform<T>],
        _ u: [T],
        _ p: [T],
        _ phi0: T,
        _ dphi0: T,
        _ aloIn: T,
        _ ahiIn: T,
        _ phiLoIn: T,
        _ phiHiIn: T,
        _ c1: T,
        _ c2: T,
        _ options: MLEOptimizationOpts<T>,
        stepMinScaled: T
    ) -> (res: (alpha: T, fNew: T, gNew: [T])?, evals: Int) {
        var evals = 0
        var alo = aloIn
        var ahi = ahiIn
        var phiLo = phiLoIn
        var phiHi = phiHiIn

        // Curvature safeguard parameters
        let curvatureRelax: T = T(0.5) // allow some noise when |dphi| shrinks slowly

        func evalAt(_ a: T) -> (f: T, g: [T])? {
            var aTry = a
            var attempts = 0
            let tauInv = T.maximum(T(1e-6), options.invalidBackoffTau)
            while attempts <= options.invalidBackoffMax {
                let uA = zip(u, p).map { (ui, pi) in ui + aTry * pi }
                let fA = nll(problem, transforms, uA, options); evals += 1
                if fA.isFinite, let gPair = gradNLL_u(problem, transforms, uA, options) {
                    return (fA, gPair.0)
                }
                aTry *= tauInv
                attempts += 1
                if aTry < stepMinScaled { break }
            }
            return nil
        }

        for _ in 0..<options.lineSearchZoomMaxIter {
            // Safeguarded interpolation attempt
            let dphiLo = { () -> T in
                // Approximate directional derivative at alo using gradient if available
                if let (g, _) = gradNLL_u(problem, transforms, zip(u, p).map { $0 + alo * $1 }, options) {
                    return dot(g, p)
                }
                return dphi0 // fallback to initial slope
            }()
            let dphiHi = { () -> T in
                if let (g, _) = gradNLL_u(problem, transforms, zip(u, p).map { $0 + ahi * $1 }, options) {
                    return dot(g, p)
                }
                return dphi0
            }()

            var alpha: T? = safeguardedInterp(alo: alo, phiLo: phiLo, dphiLo: dphiLo,
                                              ahi: ahi, phiHi: phiHi, dphiHi: dphiHi)

            // If interpolation failed or produced an out-of-bounds value, use bisection
            if alpha == nil {
                alpha = (alo + ahi) / T.two
            }

            guard let aEval = alpha, let (phi, g) = evalAt(aEval) else {
                // If invalid, shrink interval (bisection)
                let mid = (alo + ahi) / T.two
                ahi = mid
                continue
            }

            // Standard Wolfe zoom logic with curvature safeguard
            if phi > phi0 + c1 * aEval * dphi0 || phi >= phiLo {
                ahi = aEval
                phiHi = phi
            } else {
                let dphi = dot(g, p)
                // Curvature condition (relaxed slightly)
                if abs(dphi) <= c2 * abs(dphi0) || abs(dphi) <= curvatureRelax * abs(dphi0) * T(1e-3) {
                    return ((aEval, phi, g), evals)
                }
                if dphi * (ahi - alo) >= 0 {
                    ahi = alo
                    phiHi = phiLo
                }
                alo = aEval
                phiLo = phi
            }

            // If interval collapsed too much, accept bisection point if Armijo holds
            if abs(ahi - alo) <= stepMinScaled {
                let amid = (alo + ahi) / T.two
                if let (phiM, gM) = evalAt(amid) {
                    let _ = dot(gM, p)
                    if phiM <= phi0 + c1 * amid * dphi0 {
                        return ((amid, phiM, gM), evals)
                    }
                }
                break
            }
        }
        return (nil, evals)
    }

    // Final Armijo-only backtracking with capped iterations, using scaled stepMin
    private static func finalArmijoBacktrack(
        _ problem: MLEProblem<T>,
        _ transforms: [Transform<T>],
        _ u: [T],
        _ p: [T],
        _ phi0: T,
        _ dphi0: T,
        _ c1: T,
        _ options: MLEOptimizationOpts<T>,
        stepMinScaled: T,
        extraIters: Int = 20
    ) -> (res: (alpha: T, fNew: T, gNew: [T])?, evals: Int) {
        var evals = 0
        var a = T.one
        let tau = options.lineSearchTau
        let maxIter = max(20, extraIters) // widen a bit beyond previous fallback

        for _ in 0..<maxIter {
            if a < stepMinScaled { break }
            let uA = zip(u, p).map { $0 + a * $1 }
            let fA = nll(problem, transforms, uA, options); evals += 1
            if fA.isFinite, let gPair = gradNLL_u(problem, transforms, uA, options) {
                let gA = gPair.0
                if fA <= phi0 + c1 * a * dphi0 {
                    return ((a, fA, gA), evals)
                }
            }
            a *= tau
        }
        return (nil, evals)
    }

    // Compute gradient of NLL in u (analytic via theta if available, else finite differences)
    private static func gradNLL_u(
        _ problem: MLEProblem<T>,
        _ transforms: [Transform<T>],
        _ u: [T],
        _ options: MLEOptimizationOpts<T>
    ) -> ([T], Int)? {
        let k = u.count
        func gTheta(_ theta: [T]) -> [T]? {
            if let glp = problem.gradlogpdf {
                var g = [T](repeating: 0, count: k)
                for x in problem.data {
                    let gi = glp(x, theta)
                    if gi.count != k || gi.contains(where: { !$0.isFinite }) {
                        return nil
                    }
                    for i in 0..<k { g[i] += gi[i] }
                }
                return g
            }
            return nil
        }
        // Try analytic gradient via theta first
        let theta = zip(u, transforms).map { (ui, tf) in tf.toTheta(ui) }
        if let g_theta = gTheta(theta) {
            // grad_u NLL = - diag(dtheta/du) * g_theta
            var gu = [T](repeating: 0, count: k)
            for i in 0..<k {
                let dth = transforms[i].dToTheta(u[i])
                gu[i] = -(dth * g_theta[i])
            }
            if gu.contains(where: { !$0.isFinite }) { return nil }
            return (gu, 0)
        }
        // Fallback: adaptive finite-difference gradient in u with forward/central switching
        let hBase = options.gradStep
        var g = [T](repeating: 0, count: k)
        var evals = 0

        // f(u) for forward differences
        let f0 = nll(problem, transforms, u, options); evals += 1
        if !f0.isFinite { return nil }

        for i in 0..<k {
            // Adaptive step based on |u_i| and |dθ/du|
            let ui = u[i]
            let dth_du = transforms[i].dToTheta(ui)
            let scaleU = T.maximum(T.one, abs(ui))
            let scaleJac = T.maximum(T.one, abs(dth_du))
            var hi = hBase * (scaleU / scaleJac)

            // Try to get valid evaluations; shrink if invalid
            var tried = 0
            let maxTries = max(5, options.invalidBackoffMax / 2)
            let shrink = T.maximum(T(1e-3), options.invalidBackoffTau)

            var fPlus: T = .infinity
            var fMinus: T = .infinity
            var plusValid = false
            var minusValid = false

            while tried <= maxTries {
                // u + h e_i
                var up = u; up[i] = ui + hi
                fPlus = nll(problem, transforms, up, options); evals += 1
                plusValid = fPlus.isFinite

                // u - h e_i
                var um = u; um[i] = ui - hi
                fMinus = nll(problem, transforms, um, options); evals += 1
                minusValid = fMinus.isFinite

                if plusValid && minusValid { break }
                if plusValid && !minusValid {
                    // forward will be possible; still attempt a smaller h to see if central becomes valid
                    hi *= shrink
                    tried += 1
                    continue
                }
                if !plusValid {
                    // shrink and retry
                    hi *= shrink
                    tried += 1
                    continue
                }
            }

            if plusValid && minusValid {
                // Central difference
                g[i] = (fPlus - fMinus) / (T.two * hi)
            } else if plusValid {
                // Forward difference
                g[i] = (fPlus - f0) / hi
            } else {
                // Cannot produce a valid difference for this coordinate
                return nil
            }
        }
        if g.contains(where: { !$0.isFinite }) { return nil }
        return (g, evals)
    }

    // MARK: - Trust-region refinement (replacing dampedNewtonStep)

    /// Perform up to trMaxRefineIters trust-region steps using CG with finite-diff Hessian–vector products.
    /// Returns improved (u, f) if successful with proper nEval accounting.
    private static func trustRegionRefinement(
        _ problem: MLEProblem<T>,
        _ transforms: [Transform<T>],
        _ u: [T],
        _ f: T,
        _ grad: [T],
        _ options: MLEOptimizationOpts<T>
    ) -> (uHat: [T], fHat: T, nEval: Int)? {
        let k = u.count
        var nEval = 0

        var uCurr = u
        var fCurr = f
        var gCurr = grad

        var lambda = min(max(options.trInitialLambda, options.trLambdaMin), options.trLambdaMax)
        var delta = max(options.trInitialRadius, T(1e-12))

        // Helper for model prediction m(p) = f + g^T p + 0.5 p^T (H + λ I) p (approx via Hv)
        func modelReduction(_ p: [T], _ Hv: [T]) -> T {
            let gTp = dot(gCurr, p)
            let pTHp = dot(p, Hv) + lambda * dot(p, p)
            return -(gTp + T(0.5) * pTHp) // predicted reduction (positive if good)
        }

        let maxCGIter = min(k, options.trCGMaxIterFactor * k)

        for _ in 0..<options.trMaxRefineIters {
            // Compute a step p by truncated CG on (H + λ I) p = -g with trust region radius delta.
            let cg = truncatedCG(
                problem, transforms,
                uCurr, gCurr,
                lambda: lambda,
                delta: delta,
                hvEps: options.trHvEps,
                maxIter: maxCGIter,
                options: options
            )
            guard let p = cg.p else {
                nEval += cg.nEval
                // If CG failed to produce a direction, try increasing lambda and shrinking delta
                lambda = min(lambda * T(10), options.trLambdaMax)
                delta *= options.trGammaInc
                continue
            }
            nEval += cg.nEval // gradient evals used by Hv products

            // Evaluate candidate
            let uNew = zip(uCurr, p).map { $0 + $1 }
            let fNew = nll(problem, transforms, uNew, options); nEval += 1

            // If invalid, shrink delta and increase lambda; retry next iteration
            if !fNew.isFinite {
                delta *= options.trGammaInc
                lambda = min(lambda * T(10), options.trLambdaMax)
                continue
            }

            // Compute actual reduction and predicted reduction
            let ared = fCurr - fNew
            // For mred, we need H p; recompute Hv at p direction for prediction (counts extra evals)
            let (_, hvEvalsOk, okHv) = hvProduct(problem, transforms, uCurr, p, options.trHvEps, options)
            nEval += hvEvalsOk
            let Hv = okHv ?? [T](repeating: 0, count: k) // if failed, fall back to 0
            let mred = modelReduction(p, Hv)

            // Ratio ρ = ared / mred
            let rho = (mred > 0 && mred.isFinite) ? (ared / mred) : -T.infinity

            // Accept/reject and update radius
            if rho >= options.trEta2 {
                // Good model: accept and expand radius
                uCurr = uNew
                fCurr = fNew
                // update gradient at new point
                if let (gNew, gE) = gradNLL_u(problem, transforms, uCurr, options) {
                    gCurr = gNew; nEval += gE
                } else {
                    // if gradient fails, stop refinement and return current
                    return (uCurr, fCurr, nEval)
                }
                delta = max(delta, norm2(p)) * options.trGammaDec
                // Optionally decrease lambda
                lambda = max(lambda / T(2), options.trLambdaMin)
            } else if rho >= options.trEta1 {
                // Accept but keep radius
                uCurr = uNew
                fCurr = fNew
                if let (gNew, gE) = gradNLL_u(problem, transforms, uCurr, options) {
                    gCurr = gNew; nEval += gE
                } else {
                    return (uCurr, fCurr, nEval)
                }
                // Slightly adjust lambda toward smaller
                lambda = max(lambda / T(2), options.trLambdaMin)
            } else {
                // Reject: shrink radius and increase lambda
                delta *= options.trGammaInc
                lambda = min(lambda * T(2), options.trLambdaMax)
            }

            // Terminate if step is tiny or gradient small
            if norm2(p) < options.tolStep || norm2(gCurr) < options.tolGrad {
                break
            }
        }

        return (uCurr, fCurr, nEval)
    }

    /// Hessian–vector product via central finite differences of gradNLL_u:
    /// H·v ≈ (g(u + ε v) − g(u − ε v)) / (2 ε)
    /// Returns (Hv, nEval, HvOrNilIfFailed)
    private static func hvProduct(
        _ problem: MLEProblem<T>,
        _ transforms: [Transform<T>],
        _ u: [T],
        _ v: [T],
        _ eps: T,
        _ options: MLEOptimizationOpts<T>
    ) -> ([T], Int, [T]?) {
        let k = u.count
        var nEval = 0

        // Build u± = u ± ε v
        let uPlus = zip(u, v).map { $0 + eps * $1 }
        let uMinus = zip(u, v).map { $0 - eps * $1 }

        guard let (gPlus, e1) = gradNLL_u(problem, transforms, uPlus, options) else {
            return ([T](repeating: 0, count: k), nEval, nil)
        }
        nEval += e1

        guard let (gMinus, e2) = gradNLL_u(problem, transforms, uMinus, options) else {
            return ([T](repeating: 0, count: k), nEval, nil)
        }
        nEval += e2

        var Hv = [T](repeating: 0, count: k)
        let denom = T.two * eps
        for i in 0..<k { Hv[i] = (gPlus[i] - gMinus[i]) / denom }
        return (Hv, nEval, Hv)
    }

    /// Truncated Conjugate Gradient to approximately solve (H + λ I) p = −g subject to ||p|| ≤ Δ.
    /// Uses Hessian–vector products computed by hvProduct.
    private static func truncatedCG(
        _ problem: MLEProblem<T>,
        _ transforms: [Transform<T>],
        _ u: [T],
        _ g: [T],
        lambda: T,
        delta: T,
        hvEps: T,
        maxIter: Int,
        options: MLEOptimizationOpts<T>
    ) -> (p: [T]?, nEval: Int) {
        let k = u.count
        var nEval = 0

        // If gradient zero, return zero step
        let gNorm = norm2(g)
        if !gNorm.isFinite || gNorm == 0 {
            return (Array(repeating: 0, count: k), nEval)
        }

        // Initialize
        var p = [T](repeating: 0, count: k)
        var r = g.map { -$0 } // r0 = -g
        var d = r // d0 = r0
        var rTr = dot(r, r)

        // If radius is tiny, bail out
        if delta <= T(1e-16) {
            return (Array(repeating: 0, count: k), nEval)
        }

        for _ in 0..<maxIter {
            // Compute q = (H + λ I) d via Hv
            let (_, evals, maybeHv) = hvProduct(problem, transforms, u, d, hvEps, options)
            nEval += evals
            guard let Hv = maybeHv else {
                // If Hv failed, stop CG and return current p (could be zero)
                break
            }
            var q = Hv
            for i in 0..<k { q[i] += lambda * d[i] }

            let dTq = dot(d, q)
            if !(dTq > 0 && dTq.isFinite) {
                // Non-positive curvature or numerical issue — step to boundary along d
                let aBound = stepToBoundary(p: p, d: d, delta: delta)
                for i in 0..<k { p[i] += aBound * d[i] }
                return (p, nEval)
            }

            let alpha = rTr / dTq

            // Check if we cross the trust-region boundary
            let pNext = zip(p, d).map { $0 + alpha * $1 }
            if norm2(pNext) >= delta {
                // Compute boundary step
                let aBound = stepToBoundary(p: p, d: d, delta: delta)
                for i in 0..<k { p[i] += aBound * d[i] }
                return (p, nEval)
            }

            // Safe step
            p = pNext
            // Update residual r = r - alpha q
            for i in 0..<k { r[i] -= alpha * q[i] }
            let rTrNew = dot(r, r)

            // Convergence on residual
            if rTrNew <= (T(1e-12) * rTr) || T.sqrt(rTrNew) < options.tolGrad {
                break
            }

            let beta = rTrNew / rTr
            // d = r + beta d
            for i in 0..<k { d[i] = r[i] + beta * d[i] }
            rTr = rTrNew
        }

        return (p, nEval)
    }

    /// Compute the step length a ≥ 0 such that ||p + a d|| = delta (quadratic solve).
    @inline(__always) private static func stepToBoundary(p: [T], d: [T], delta: T) -> T {
        // Solve ||p + a d||^2 = delta^2 => (d·d) a^2 + 2 (p·d) a + (p·p - delta^2) = 0
        let dd = dot(d, d)
        let pd = dot(p, d)
        let pp = dot(p, p)
        let c = pp - delta * delta
        if dd == 0 { return T.zero }
        let disc = pd * pd - dd * c
        if disc <= 0 { return T.zero }
        let a = (-pd + T.sqrt(disc)) / dd
        return max(a, T.zero)
    }

    // Rough count of function evals for Hessian (central differences)
    @inline(__always) private static func hessEvalCost(_ k: Int) -> Int {
        // Includes f0, diagonals (2 per dim), and off-diagonals (4 per pair)
        return 1 + k * 2 + (k * (k - 1) / 2) * 4
    }

    // MARK: - Multi-start candidate generation

    /// Generate theta-candidates respecting constraints using chosen design around base theta.
    /// Enhanced: prune invalid/high-NLL starts and resample up to a cap; Sobol with scrambling.
    private static func generateCandidates(
        _ specs: [ParamSpec<T>],
        baseTheta: [T],
        count m: Int,
        design: MultiStartDesign,
        logScale: T,
        gridPerDim: Int,
        rngSeed: UInt64?,
        problem: MLEProblem<T>,
        transforms: [Transform<T>],
        options: MLEOptimizationOpts<T>
    ) -> [[T]] {
        if m <= 0 { return [] }
        let k = specs.count

        // Build per-parameter ranges or transforms to map U(0,1) -> feasible theta
        let maps: [(T, Int) -> T] = specs.enumerated().map { (i, spec) in
            let base = baseTheta[i]
            switch spec.constraint {
            case .real:
                let r: T = T(3)
                return { u01, _ in base + r * (u01 * 2 - 1) }
            case .positive:
                let s = logScale
                let b = max(base, T(1e-12))
                return { u01, _ in b * T.exp(s * (u01 * 2 - 1)) }
            case .unitInterval:
                return { u01, _ in
                    let eps: T = 1e-6
                    let v = eps + (1 - 2*eps) * u01
                    return v
                }
            case .lowerBound(let a):
                let s = logScale
                let span = max(base - a, T(1e-6))
                return { u01, _ in
                    let v = a + span * T.exp(s * (u01 * 2 - 1))
                    return max(v, a + T(1e-9))
                }
            case .upperBound(let b):
                let s = logScale
                let span = max(b - base, T(1e-6))
                return { u01, _ in
                    let v = b - span * T.exp(s * (u01 * 2 - 1))
                    return min(v, b - T(1e-9))
                }
            case .interval(let a, let b):
                return { u01, _ in
                    let eps: T = 1e-9
                    let lo = a + eps, hi = b - eps
                    return lo + (hi - lo) * u01
                }
            }
        }

        // Low-discrepancy samplers in [0,1]^k
        func sampleU01LHS(_ n: Int, _ d: Int) -> [[T]] {
            var M = Array(repeating: Array(repeating: T.zero, count: d), count: n)
            var rng = SeededRNG(seed: rngSeed)
            for j in 0..<d {
                var slots = (0..<n).map { ($0, (T($0) + rng.uniform01()) / T(n)) }
                slots.shuffle(using: &rng)
                for i in 0..<n {
                    M[i][j] = slots[i].1
                }
            }
            return M
        }

        // Owen-style digit scrambling for van-der-Corput with base 2
        func scrambledVDC(_ index: Int, scramble: UInt64) -> T {
            var i = UInt64(index)
            var bit: UInt64 = 1
            var invBase: T = T(0.5)
            var value: T = 0
            var scr = scramble
            while i > 0 || bit <= (1 << 20) {
                let di = i & 1
                let sbit = scr & 1
                let digit = di ^ sbit
                if digit == 1 {
                    value += invBase
                }
                invBase *= 0.5
                i >>= 1
                scr >>= 1
                bit <<= 1
                if i == 0 && invBase < T(1e-12) { break }
            }
            return value
        }

        // Sobol sequence using direction numbers with mild digital scrambling
        func sampleU01SobolScrambled(_ n: Int, _ d: Int) -> [[T]] {
            let baseSeed = rngSeed ?? MLESolver.sobolDefaultSeed
            var out = Array(repeating: Array(repeating: T.zero, count: d), count: n)
            let bitCount = 60
            for j in 0..<d {
                let scramble = MLESolver.sobolScrambleSeed(baseSeed: baseSeed, dimension: j)
                if let directions = MLESolver.sobolDirections(forDimension: j, bits: bitCount) {
                    for i in 0..<n {
                        out[i][j] = MLESolver.sobolComponent(index: i + 1,
                                                             directions: directions,
                                                             scramble: scramble)
                    }
                } else {
                    for i in 0..<n {
                        out[i][j] = scrambledVDC(i + 1, scramble: scramble)
                    }
                }
            }
            return out
        }

        func sampleU01Random(_ n: Int, _ d: Int) -> [[T]] {
            var rng = SeededRNG(seed: rngSeed)
            return (0..<n).map { _ in (0..<d).map { _ in rng.uniform01() } }
        }

        func sampleU01Grid(_ n: Int, _ d: Int) -> [[T]] {
            let mPerDim = max(2, gridPerDim)
            let levels: [T] = (0..<mPerDim).map { T($0 + 1) / T(mPerDim + 1) }
            if d == 1 {
                return levels.prefix(n).map { [$0] }
            }
            var grid: [[T]] = [[]]
            for _ in 0..<d {
                var next: [[T]] = []
                for p in grid {
                    for lv in levels {
                        var q = p
                        q.append(lv)
                        next.append(q)
                    }
                }
                grid = next
            }
            return Array(grid.prefix(n))
        }

        let U: [[T]]
        switch design {
        case .lhs:   U = sampleU01LHS(m, k)
        case .sobol: U = sampleU01SobolScrambled(m, k) // improved coverage vs Halton for larger k
        case .random:U = sampleU01Random(m, k)
        case .grid:  U = sampleU01Grid(m, k)
        }

        // Map U -> theta and prune bad starts by evaluating NLL in u-space, resampling if needed.
        var rngResample = SeededRNG(seed: (rngSeed ?? 0xDEAD_BEEF_F00D_BABE) ^ 0xBADC0FFE)
        var theta: [[T]] = []
        theta.reserveCapacity(U.count)

        let uPool = U // initial pool; for resamples, draw fresh randoms deterministically
        let cutoff = options.startNLLCutoff

        func drawRandomU() -> [T] {
            // uniform random in [0,1]^k for resampling
            return (0..<k).map { _ in rngResample.uniform01() as T }
        }

        for i in 0..<U.count {
            var tries = 0
            var accepted: [T]? = nil
            while tries <= max(0, options.maxStartResampleAttempts) {
                let u01 = (tries == 0) ? uPool[i] : drawRandomU()
                var th = [T](repeating: .zero, count: k)
                for j in 0..<k {
                    th[j] = maps[j](u01[j], j)
                }
                // Evaluate NLL at corresponding u to decide acceptance
                let uCandidate = zip(th, transforms).map { (theta, tf) in tf.toU(theta) }
                let f = nll(problem, transforms, uCandidate, options)
                let bad = !f.isFinite || (cutoff != nil && f > cutoff!)
                if !bad {
                    accepted = th
                    break
                }
                tries += 1
            }
            if let th = accepted {
                theta.append(th)
            } else {
                // If we failed to find a good one, still append the original mapped candidate to keep count
                var th = [T](repeating: .zero, count: k)
                for j in 0..<k {
                    th[j] = maps[j](U[i][j], j)
                }
                theta.append(th)
            }
        }

        return theta
    }

    /// Create jittered u-starts around u0 for potentially unbounded parameters to escape saddle points.
    private static func jitterAroundU(
        _ u0: [T],
        _ specs: [ParamSpec<T>],
        _ transforms: [Transform<T>],
        count: Int,
        radius: T,
        rngSeed: UInt64?
    ) -> [[T]] {
        // Currently unused (count=0). Hook kept for future.
        if count <= 0 { return [] }
        var rng = SeededRNG(seed: rngSeed)
        return (0..<count).map { _ in
            zip(u0, specs).map { (u, spec) in
                switch spec.constraint {
                case .real:
                    return u + (rng.uniform01() * 2 - 1) * radius
                default:
                    return u // leave constrained params unchanged here
                }
            }
        }
    }

    private static var sobolDefaultSeed: UInt64 { 0xA5A5_1234_9E37_F00D }
    private static var sobolNormalizer: Double { 1.0 / (Double(UInt64.max) + 1.0) }

    private static func sobolScrambleSeed(baseSeed: UInt64, dimension: Int) -> UInt64 {
        let golden: UInt64 = 0x9E37_79B1_85EB_CA87
        var scramble = baseSeed &+ UInt64(dimension + 1) &* golden
        if scramble == 0 { scramble = 0x1 }
        return scramble
    }

    private static func sobolComponent(index: Int, directions: [UInt64], scramble: UInt64) -> T {
        var idx = index
        var acc = scramble
        var bit = 1
        while idx > 0 && bit < directions.count {
            if (idx & 1) == 1 {
                acc ^= directions[bit]
            }
            idx >>= 1
            bit += 1
        }
        let value = Double(acc) * sobolNormalizer
        return T(value)
    }

    private static func sobolDirections(forDimension dimension: Int, bits: Int) -> [UInt64]? {
        guard dimension < sobolParameterTable.count else { return nil }
        let params = sobolParameterTable[dimension]
        let L = max(bits, params.s)
        var v = Array(repeating: UInt64(0), count: L + 1)
        for j in 1...params.s {
            v[j] = UInt64(params.m[j - 1]) << (64 - j)
        }
        if params.s < L {
            for j in (params.s + 1)...L {
                var value = v[j - params.s] >> params.s
                for k in 1...params.s {
                    if ((params.a >> (params.s - k)) & 1) == 1 {
                        value ^= v[j - k]
                    }
                }
                v[j] = value
            }
        }
        return v
    }

    /// Count unique u-hat up to relative Euclidean tolerance to assess multi-start diversity.
    private static func countUniqueUSolutions(_ us: [[T]], tolRel: T) -> Int {
        var uniq: [[T]] = []
        outer: for u in us {
            for v in uniq {
                let diff2: T = zip(u, v).map { (a, b) in
                    let d = a - b
                    return d * d
                }.reduce(0, +)
                let normU = T.sqrt(u.map { $0 * $0 }.reduce(0, +))
                let normV = T.sqrt(v.map { $0 * $0 }.reduce(0, +))
                let denom = T.maximum(T.one, T.maximum(normU, normV))
                let rel = T.sqrt(diff2) / denom
                if rel <= tolRel {
                    continue outer
                }
            }
            uniq.append(u)
        }
        return uniq.count
    }

    /// Count unique θ solutions by normalizing differences with parameter-range surrogates.
    private static func countUniqueThetaSolutions(_ thetas: [[T]], specs: [ParamSpec<T>], tol: T) -> Int {
        guard !thetas.isEmpty, !specs.isEmpty else { return 0 }
        let scales = specs.map { thetaScale(for: $0.constraint) }
        var uniq: [[T]] = []
        outer: for theta in thetas {
            for existing in uniq {
                var dist2: T = 0
                for i in 0..<min(theta.count, existing.count) {
                    let denom = scales[i] ?? T.maximum(T.one, T.maximum(abs(theta[i]), abs(existing[i])))
                    let delta = (theta[i] - existing[i]) / denom
                    dist2 += delta * delta
                }
                if T.sqrt(dist2) <= tol {
                    continue outer
                }
            }
            uniq.append(theta)
        }
        return uniq.count
    }

    /// θ-scale heuristic: exact range for interval/unit constraints, nil otherwise.
    private static func thetaScale(for constraint: ParamConstraint<T>) -> T? {
        switch constraint {
        case .interval(let lo, let hi):
            let span = hi - lo
            return span.isFinite && span > 0 ? span : nil
        case .unitInterval:
            return T.one
        default:
            return nil
        }
    }

    /// Gershgorin disc bound: κ ≥ λ_max/λ_min with radii derived from row sums.
    private static func estimateConditionNumberGershgorin(_ matrix: [[T]]) -> T? {
        guard !matrix.isEmpty else { return nil }
        let n = matrix.count
        guard matrix.allSatisfy({ $0.count == n }) else { return nil }
        var minBound = T.infinity
        var maxBound: T = 0
        for i in 0..<n {
            let diag = matrix[i][i]
            if !diag.isFinite { return nil }
            var radius: T = 0
            for j in 0..<n where j != i {
                radius += abs(matrix[i][j])
            }
            let lower = diag - radius
            let upper = diag + radius
            minBound = T.minimum(minBound, lower)
            maxBound = T.maximum(maxBound, upper)
        }
        if !(minBound > 0 && minBound.isFinite && maxBound.isFinite) {
            return nil
        }
        let cond = maxBound / minBound
        return cond.isFinite ? cond : nil
    }

    // MARK: - Utility: relative change in NLL (equivalently |Delta logLik|)

    @inline(__always) private static func relChangeVal(_ prev: T, _ curr: T) -> T {
        let num = abs(curr - prev)
        let den = T.maximum(T.one, abs(curr))
        return num / den
    }

    // MARK: - Linear algebra helpers (small dense)

    @inline(__always) private static func dot(_ a: [T], _ b: [T]) -> T {
        #if canImport(Accelerate)
        if T.self == Double.self {
            let res = cblas_ddot(Int32(a.count), a as! [Double], 1, b as! [Double], 1)
            return (res as Double) as! T
        }
        if T.self == Float.self {
            let res = cblas_sdot(Int32(a.count), a as! [Float], 1, b as! [Float], 1)
            return (res as Float) as! T
        }
        #endif
        var s: T = 0
        for i in 0..<a.count { s += a[i] * b[i] }
        return s
    }

    @inline(__always) private static func norm2(_ a: [T]) -> T {
        #if canImport(Accelerate)
        if T.self == Double.self {
            let res = cblas_dnrm2(Int32(a.count), a as! [Double], 1)
            return (res as Double) as! T
        }
        if T.self == Float.self {
            let res = cblas_snrm2(Int32(a.count), a as! [Float], 1)
            return (res as Float) as! T
        }
        #endif
        return T.sqrt(dot(a, a))
    }

    @inline(__always) private static func matvec(_ A: [[T]], _ x: [T]) -> [T] {
        let n = A.count
        #if canImport(Accelerate)
        if T.self == Double.self {
            let Ad = A as! [[Double]]
            var flat = toRowMajorDouble(Ad)
            let xd = x as! [Double]
            var y = [Double](repeating: 0, count: n)
            cblas_dgemv(CBLAS_ORDER(CblasRowMajor.rawValue), CBLAS_TRANSPOSE(CblasNoTrans.rawValue),
                        Int32(n), Int32(n),
                        1.0, &flat, Int32(n),
                        xd, 1,
                        0.0, &y, 1)
            return y as! [T]
        }
        if T.self == Float.self {
            let Af = A as! [[Float]]
            var flat = toRowMajorFloat(Af)
            let xf = x as! [Float]
            var y = [Float](repeating: 0, count: n)
            cblas_sgemv(CBLAS_ORDER(CblasRowMajor.rawValue), CBLAS_TRANSPOSE(CblasNoTrans.rawValue),
                        Int32(n), Int32(n),
                        1.0, &flat, Int32(n),
                        xf, 1,
                        0.0, &y, 1)
            return y as! [T]
        }
        #endif
        var y = [T](repeating: 0, count: n)
        for i in 0..<n {
            var s: T = 0
            for j in 0..<n { s += A[i][j] * x[j] }
            y[i] = s
        }
        return y
    }

    @inline(__always) private static func outer(_ a: [T], _ b: [T]) -> [[T]] {
        let n = a.count
        #if canImport(Accelerate)
        if T.self == Double.self {
            let ad = a as! [Double]
            let bd = b as! [Double]
            var y = [Double](repeating: 0, count: n * n)
            // y := a*b^T + y  (start with zeros => outer product)
            cblas_dger(CBLAS_ORDER(CblasRowMajor.rawValue),
                       Int32(n), Int32(n),
                       1.0,
                       ad, 1,
                       bd, 1,
                       &y, Int32(n))
            return fromRowMajorDouble(y, n) as! [[T]]
        }
        if T.self == Float.self {
            let af = a as! [Float]
            let bf = b as! [Float]
            var y = [Float](repeating: 0, count: n * n)
            cblas_sger(CBLAS_ORDER(CblasRowMajor.rawValue),
                       Int32(n), Int32(n),
                       1.0,
                       af, 1,
                       bf, 1,
                       &y, Int32(n))
            return fromRowMajorFloat(y, n) as! [[T]]
        }
        #endif
        var M = Array(repeating: Array(repeating: T.zero, count: n), count: n)
        for i in 0..<n {
            for j in 0..<n {
                M[i][j] = a[i] * b[j]
            }
        }
        return M
    }

    @inline(__always) private static func add(_ A: [[T]], _ B: [[T]]) -> [[T]] {
        let n = A.count
        var C = A
        for i in 0..<n {
            for j in 0..<n {
                C[i][j] += B[i][j]
            }
        }
        return C
    }

    @inline(__always) private static func scale(_ A: [[T]], _ s: T) -> [[T]] {
        let n = A.count
        var C = A
        for i in 0..<n {
            for j in 0..<n {
                C[i][j] *= s
            }
        }
        return C
    }

    @inline(__always) private static func matmul(_ A: [[T]], _ B: [[T]]) -> [[T]] {
        let n = A.count
        #if canImport(Accelerate)
        if T.self == Double.self {
            let Ad = A as! [[Double]]
            let Bd = B as! [[Double]]
            var aFlat = toRowMajorDouble(Ad)
            var bFlat = toRowMajorDouble(Bd)
            var cFlat = [Double](repeating: 0, count: n * n)
            cblas_dgemm(CBLAS_ORDER(CblasRowMajor.rawValue),
                        CBLAS_TRANSPOSE(CblasNoTrans.rawValue), CBLAS_TRANSPOSE(CblasNoTrans.rawValue),
                        Int32(n), Int32(n), Int32(n),
                        1.0,
                        &aFlat, Int32(n),
                        &bFlat, Int32(n),
                        0.0,
                        &cFlat, Int32(n))
            return fromRowMajorDouble(cFlat, n) as! [[T]]
        }
        if T.self == Float.self {
            let Af = A as! [[Float]]
            let Bf = B as! [[Float]]
            var aFlat = toRowMajorFloat(Af)
            var bFlat = toRowMajorFloat(Bf)
            var cFlat = [Float](repeating: 0, count: n * n)
            cblas_sgemm(CBLAS_ORDER(CblasRowMajor.rawValue),
                        CBLAS_TRANSPOSE(CblasNoTrans.rawValue), CBLAS_TRANSPOSE(CblasNoTrans.rawValue),
                        Int32(n), Int32(n), Int32(n),
                        1.0,
                        &aFlat, Int32(n),
                        &bFlat, Int32(n),
                        0.0,
                        &cFlat, Int32(n))
            return fromRowMajorFloat(cFlat, n) as! [[T]]
        }
        #endif
        var C = Array(repeating: Array(repeating: T.zero, count: n), count: n)
        for i in 0..<n {
            for k in 0..<n {
                let aik = A[i][k]
                if aik == 0 { continue }
                for j in 0..<n {
                    C[i][j] += aik * B[k][j]
                }
            }
        }
        return C
    }

    @inline(__always) private static func eye(_ n: Int) -> [[T]] {
        var I = Array(repeating: Array(repeating: T.zero, count: n), count: n)
        for i in 0..<n { I[i][i] = T.one }
        return I
    }

    /// Two-loop recursion for L-BFGS direction with gamma clamped.
    private static func lbfgsDirectionWithGammaClamped(
        _ grad: [T],
        _ sList: [[T]],
        _ yList: [[T]],
        _ rhoList: [T],
        floor: T,
        ceil: T
    ) -> [T] {
        var q = grad
        let m = sList.count
        var alpha = [T](repeating: 0, count: m)
        // First loop
        for i in stride(from: m - 1, through: 0, by: -1) {
            let s = sList[i]
            let y = yList[i]
            let rho = rhoList[i]
            let a = rho * dot(s, q)
            alpha[i] = a
            for j in 0..<q.count {
                q[j] -= a * y[j]
            }
        }
        // Scaling of initial Hessian (H0 = gamma I), gamma = (s^T y) / (y^T y), clamped
        var r = q
        if let lastS = sList.last, let lastY = yList.last {
            let sy = dot(lastS, lastY)
            let yy = dot(lastY, lastY)
            if sy > 0 && yy > 0 && sy.isFinite && yy.isFinite {
                var gamma = sy / yy
                if !gamma.isFinite || gamma <= 0 {
                    gamma = T.one
                }
                gamma = max(floor, min(ceil, gamma))
                for j in 0..<r.count { r[j] *= gamma }
            }
        }
        // Second loop
        for i in 0..<m {
            let s = sList[i]
            let y = yList[i]
            let rho = rhoList[i]
            let beta = rho * dot(y, r)
            for j in 0..<r.count {
                r[j] += s[j] * (alpha[i] - beta)
            }
        }
        // Direction is -H_k * grad => negative of r
        return r.map { -$0 }
    }

    /// Reconstruct a dense L-BFGS inverse Hessian approximation H by applying the two-loop
    /// recursion to each canonical basis vector e_i to compute H e_i as columns.
    private static func reconstructDenseLBFGSInverseH(
        sList: [[T]],
        yList: [[T]],
        rhoList: [T],
        dim: Int,
        floor: T,
        ceil: T
    ) -> [[T]]? {
        guard dim > 0 else { return nil }
        if sList.isEmpty || yList.isEmpty || rhoList.isEmpty { return nil }
        let m = sList.count
        guard yList.count == m, rhoList.count == m else { return nil }

        func twoLoopApply(_ v: [T]) -> [T] {
            var q = v
            var alpha = [T](repeating: 0, count: m)
            // First loop
            for i in stride(from: m - 1, through: 0, by: -1) {
                let s = sList[i]
                let y = yList[i]
                let rho = rhoList[i]
                let a = rho * dot(s, q)
                alpha[i] = a
                for j in 0..<q.count {
                    q[j] -= a * y[j]
                }
            }
            // Initial H0 scaling (gamma I)
            var r = q
            if let lastS = sList.last, let lastY = yList.last {
                let sy = dot(lastS, lastY)
                let yy = dot(lastY, lastY)
                if sy > 0 && yy > 0 && sy.isFinite && yy.isFinite {
                    var gamma = sy / yy
                    if !gamma.isFinite || gamma <= 0 {
                        gamma = T.one
                    }
                    gamma = max(floor, min(ceil, gamma))
                    for j in 0..<r.count { r[j] *= gamma }
                }
            }
            // Second loop
            for i in 0..<m {
                let s = sList[i]
                let y = yList[i]
                let rho = rhoList[i]
                let beta = rho * dot(y, r)
                for j in 0..<r.count {
                    r[j] += s[j] * (alpha[i] - beta)
                }
            }
            return r
        }

        // Build H by applying to e_i
        var H = Array(repeating: Array(repeating: T.zero, count: dim), count: dim)
        for i in 0..<dim {
            var e = [T](repeating: 0, count: dim)
            e[i] = T.one
            let col = twoLoopApply(e)
            for j in 0..<dim {
                H[j][i] = col[j]
            }
        }
        return symmetrize(H)
    }

    /// Simplex “diameter” as the maximum Euclidean distance among vertices.
    private static func simplexDiameter(_ S: [[T]]) -> T {
        var d: T = 0
        for i in 0..<S.count {
            for j in i+1..<S.count {
                let v: T = zip(S[i], S[j]).map { (a, b) in
                    let diff = a - b
                    return diff * diff
                }.reduce(0, +)
                d = T.maximum(d, T.sqrt(v))
            }
        }
        return d
    }

    /// Numerical Hessian (central differences) of the NLL in u-coordinates.
    /// Off-diagonal entries via mixed differences; diagonal via standard second differences.
    private static func numericHessian(
        _ problem: MLEProblem<T>,
        _ transforms: [Transform<T>],
        _ u: [T],
        step h: T,
        options: MLEOptimizationOpts<T>
    ) -> [[T]] {
        let k = u.count
        var H = Array(repeating: Array(repeating: T(0), count: k), count: k)
        let f0 = nll(problem, transforms, u, options)

        for i in 0..<k {
            var uiP = u; uiP[i] += h
            var uiM = u; uiM[i] -= h
            let fP = nll(problem, transforms, uiP, options)
            let fM = nll(problem, transforms, uiM, options)
            H[i][i] = (fP - T(2)*f0 + fM) / (h*h)

            for j in i+1..<k {
                var uPP = u; uPP[i] += h; uPP[j] += h
                var uPM = u; uPM[i] += h; uPM[j] -= h
                var uMP = u; uMP[i] -= h; uMP[j] += h
                var uMM = u; uMM[i] -= h; uMM[j] -= h
                let fPP = nll(problem, transforms, uPP, options)
                let fPM = nll(problem, transforms, uPM, options)
                let fMP = nll(problem, transforms, uMP, options)
                let fMM = nll(problem, transforms, uMM, options)
                let val = (fPP - fPM - fMP + fMM) / (T(4) * h * h)
                H[i][j] = val
                H[j][i] = val
            }
        }
        return H
    }

    /// Symmetrize a square matrix: 0.5 * (A + A^T)
    @inline(__always) private static func symmetrize(_ A: [[T]]) -> [[T]] {
        let n = A.count
        var S = A
        for i in 0..<n {
            for j in i+1..<n {
                let v = (A[i][j] + A[j][i]) * T.half
                S[i][j] = v
                S[j][i] = v
            }
        }
        return S
    }

    /// Matrix inversion with platform-optimized LAPACK on Apple platforms and a portable fallback elsewhere.
    private static func invertSymmetric(_ A: [[T]]) -> [[T]]? {
        let n = A.count
        guard n > 0, A.allSatisfy({ $0.count == n }) else { return nil }

        #if canImport(Accelerate)
        // Accelerate path for Double
        if T.self == Double.self {
            var aCol = toColumnMajorDouble(A as! [[Double]])
            var N = __CLPK_integer(n)
            var lda = N
            var ipiv = [__CLPK_integer](repeating: 0, count: n)
            var info: __CLPK_integer = 0
            var N1 = N
            // LU factorization
            dgetrf_(&N, &N1, &aCol, &lda, &ipiv, &info)
            if info != 0 { return nil }

            // Workspace query for dgetri
            var lwork: __CLPK_integer = -1
            var wkopt: Double = 0
            dgetri_(&N, &aCol, &lda, &ipiv, &wkopt, &lwork, &info)
            if info != 0 { return nil }
            lwork = __CLPK_integer(wkopt)
            var work = [Double](repeating: 0, count: Int(max(lwork, 1)))

            // Inversion from LU
            dgetri_(&N, &aCol, &lda, &ipiv, &work, &lwork, &info)
            if info != 0 { return nil }

            let inv = fromColumnMajorDouble(aCol, n)
            return inv as? [[T]]
        }

        // Accelerate path for Float
        if T.self == Float.self {
            var aCol = toColumnMajorFloat(A as! [[Float]])
            var N = __CLPK_integer(n)
            var lda = N
            var ipiv = [__CLPK_integer](repeating: 0, count: n)
            var info: __CLPK_integer = 0
            var N1 = N

            // LU factorization
            sgetrf_(&N, &N1, &aCol, &lda, &ipiv, &info)
            if info != 0 { return nil }

            // Workspace query for sgetri
            var lwork: __CLPK_integer = -1
            var wkopt: Float = 0
            sgetri_(&N, &aCol, &lda, &ipiv, &wkopt, &lwork, &info)
            if info != 0 { return nil }
            lwork = __CLPK_integer(wkopt)
            var work = [Float](repeating: 0, count: Int(max(lwork, 1)))

            // Inversion from LU
            sgetri_(&N, &aCol, &lda, &ipiv, &work, &lwork, &info)
            if info != 0 { return nil }

            let inv = fromColumnMajorFloat(aCol, n)
            return inv as? [[T]]
        }
        #endif

        // Portable fallback (works on Linux and for non Float/Double T)
        return invertSymmetricNaive(A)
    }

    /// Attempt Cholesky-based inversion with diagonal regularization fallback.
    /// Returns nil if Cholesky fails up to regMax or on unsupported types/platforms.
    private static func invertSymmetricCholesky(
        _ A: [[T]],
        regInit: T,
        regFactor: T,
        regMax: T
    ) -> [[T]]? {
        let n = A.count
        guard n > 0, A.allSatisfy({ $0.count == n }) else { return nil }
        let S = symmetrize(A)

        // Small diagonal jitter to reduce borderline indefiniteness
        let jitter: T = T(1e-10)

        #if canImport(Accelerate)
        if T.self == Double.self {
            var lambda = max(T.zero, regInit)
            while lambda <= regMax {
                var a = toColumnMajorDouble(S as! [[Double]])
                // Apply λI and jitter
                for i in 0..<n {
                    a[i + i*n] += Double(lambda + jitter)
                }
                var uplo = Int8(UInt8(ascii: "L"))
                var N = __CLPK_integer(n)
                var lda = N
                var info: __CLPK_integer = 0
                dpotrf_(&uplo, &N, &a, &lda, &info)
                if info == 0 {
                    // Invert from Cholesky factor
                    dpotri_(&uplo, &N, &a, &lda, &info)
                    if info == 0 {
                        var M = fromColumnMajorDouble(a, n)
                        for i in 0..<n {
                            for j in i+1..<n {
                                M[i][j] = M[j][i]
                            }
                        }
                        return M as? [[T]]
                    }
                }
                lambda *= regFactor
                if lambda == 0 { lambda = T.one * regFactor }
            }
            return nil
        }
        if T.self == Float.self {
            var lambda = max(T.zero, regInit)
            while lambda <= regMax {
                var a = toColumnMajorFloat(S as! [[Float]])
                // Apply λI and jitter
                for i in 0..<n {
                    a[i + i*n] += Float(lambda + jitter)
                }
                var uplo = Int8(UInt8(ascii: "L"))
                var N = __CLPK_integer(n)
                var lda = N
                var info: __CLPK_integer = 0
                spotrf_(&uplo, &N, &a, &lda, &info)
                if info == 0 {
                    spotri_(&uplo, &N, &a, &lda, &info)
                    if info == 0 {
                        var M = fromColumnMajorFloat(a, n)
                        for i in 0..<n {
                            for j in i+1..<n {
                                M[i][j] = M[j][i]
                            }
                        }
                        return M as? [[T]]
                    }
                }
                lambda *= regFactor
                if lambda == 0 { lambda = T.one * regFactor }
            }
            return nil
        }
        #endif

        // Non-Accelerate or non Float/Double types: cannot do Cholesky here
        return nil
    }

    /// Previous Gauss–Jordan elimination with naive pivoting kept as a portable fallback.
    private static func invertSymmetricNaive(_ A: [[T]]) -> [[T]]? {
        let n = A.count
        var M = A
        var I = (0..<n).map { i in
            (0..<n).map { j in i == j ? T(1) : T(0) }
        }

        // naive pivot choice (largest |M[p][p]| from current p)
        for p in 0..<n {
            var pivot = p
            var maxAbs = abs(M[p][p])
            for r in (p+1)..<n {
                let v = abs(M[r][p])
                if v > maxAbs {
                    maxAbs = v; pivot = r
                }
            }
            if maxAbs == 0 || !maxAbs.isFinite { return nil }

            if pivot != p {
                M.swapAt(p, pivot)
                I.swapAt(p, pivot)
            }

            let diag = M[p][p]
            let invDiag = T(1) / diag
            for j in 0..<n { M[p][j] *= invDiag; I[p][j] *= invDiag }
            for r in 0..<n where r != p {
                let factor = M[r][p]
                if factor == 0 { continue }
                for j in 0..<n {
                    M[r][j] -= factor * M[p][j]
                    I[r][j] -= factor * I[p][j]
                }
            }
        }
        return I
    }

    // PD check via Cholesky attempt (Accelerate if available, else naive test)
    private static func isSymPosDef(_ A: [[T]]) -> Bool {
        let n = A.count
        guard n > 0, A.allSatisfy({ $0.count == n }) else { return false }
        #if canImport(Accelerate)
        if T.self == Double.self {
            var a = toColumnMajorDouble(A as! [[Double]])
            var uplo = Int8(UInt8(ascii: "L"))
            var N = __CLPK_integer(n)
            var lda = N
            var info: __CLPK_integer = 0
            dpotrf_(&uplo, &N, &a, &lda, &info)
            return info == 0
        }
        if T.self == Float.self {
            var a = toColumnMajorFloat(A as! [[Float]])
            var uplo = Int8(UInt8(ascii: "L"))
            var N = __CLPK_integer(n)
            var lda = N
            var info: __CLPK_integer = 0
            spotrf_(&uplo, &N, &a, &lda, &info)
            return info == 0
        }
        #endif
        // Gershgorin-based quick check (very rough)
        for i in 0..<n {
            var sum: T = 0
            for j in 0..<n where j != i { sum += abs(A[i][j]) }
            if A[i][i] <= sum { return false }
        }
        return true
    }

    #if canImport(Accelerate)
    // MARK: - LAPACK helpers (column-major conversions)

    @inline(__always) private static func toColumnMajorDouble(_ A: [[Double]]) -> [Double] {
        let n = A.count
        var out = [Double](repeating: 0, count: n * n)
        for j in 0..<n {
            for i in 0..<n {
                out[j * n + i] = A[i][j]
            }
        }
        return out
    }

    @inline(__always) private static func fromColumnMajorDouble(_ a: [Double], _ n: Int) -> [[Double]] {
        var M = Array(repeating: Array(repeating: 0.0, count: n), count: n)
        for j in 0..<n {
            for i in 0..<n {
                M[i][j] = a[j * n + i]
            }
        }
        return M
    }

    @inline(__always) private static func toColumnMajorFloat(_ A: [[Float]]) -> [Float] {
        let n = A.count
        var out = [Float](repeating: 0, count: n * n)
        for j in 0..<n {
            for i in 0..<n {
                out[j * n + i] = A[i][j]
            }
        }
        return out
    }

    @inline(__always) private static func fromColumnMajorFloat(_ a: [Float], _ n: Int) -> [[Float]] {
        var M = Array(repeating: Array(repeating: Float(0), count: n), count: n)
        for j in 0..<n {
            for i in 0..<n {
                M[i][j] = a[j * n + i]
            }
        }
        return M
    }

    // MARK: - CBLAS helpers (row-major conversions used for gemv/gemm/ger)
    @inline(__always) private static func toRowMajorDouble(_ A: [[Double]]) -> [Double] {
        let n = A.count
        var out = [Double](repeating: 0, count: n * n)
        for i in 0..<n {
            for j in 0..<n {
                out[i * n + j] = A[i][j]
            }
        }
        return out
    }

    @inline(__always) private static func fromRowMajorDouble(_ a: [Double], _ n: Int) -> [[Double]] {
        var M = Array(repeating: Array(repeating: 0.0, count: n), count: n)
        for i in 0..<n {
            for j in 0..<n {
                M[i][j] = a[i * n + j]
            }
        }
        return M
    }

    @inline(__always) private static func toRowMajorFloat(_ A: [[Float]]) -> [Float] {
        let n = A.count
        var out = [Float](repeating: 0, count: n * n)
        for i in 0..<n {
            for j in 0..<n {
                out[i * n + j] = A[i][j]
            }
        }
        return out
    }

    @inline(__always) private static func fromRowMajorFloat(_ a: [Float], _ n: Int) -> [[Float]] {
        var M = Array(repeating: Array(repeating: Float(0), count: n), count: n)
        for i in 0..<n {
            for j in 0..<n {
                M[i][j] = a[i * n + j]
            }
        }
        return M
    }
    #endif
}

// MARK: - Tiny RNG & low-discrepancy helpers (internal)

fileprivate struct SeededRNG: RandomNumberGenerator {
    private var state: UInt64
    init(seed: UInt64? = nil) {
        if let s = seed {
            // Avoid zero seed lock; mix with golden ratio constant
            state = s ^ 0x9E3779B97F4A7C15
            if state == 0 { state = 0x9E3779B97F4A7C15 }
        } else {
            state = 0x9E3779B97F4A7C15
        }
    }
    @inline(__always) mutating func next() -> UInt64 {
        state &+= 0x9E3779B97F4A7C15
        var z = state
        z = (z ^ (z >> 30)) &* 0xBF58476D1CE4E5B9
        z = (z ^ (z >> 27)) &* 0x94D049BB133111EB
        return z ^ (z >> 31)
    }
    @inline(__always) mutating func uniform01<T: RealLike>() -> T {
        let x = next() >> 11 // 53 bits
        let d = Double(x) / Double(1 << 53)
        return T(d)
    }
}

/// Direction numbers for low-dimensional Sobol sequences (Joe–Kuo style).
private struct SobolParameter {
    let s: Int
    let a: Int
    let m: [UInt32]
}

private let sobolParameterTable: [SobolParameter] = [
    SobolParameter(s: 1, a: 0, m: [1]),
    SobolParameter(s: 2, a: 1, m: [1, 3]),
    SobolParameter(s: 3, a: 1, m: [1, 3, 1]),
    SobolParameter(s: 3, a: 2, m: [1, 1, 1]),
    SobolParameter(s: 4, a: 1, m: [1, 3, 5, 13]),
    SobolParameter(s: 4, a: 4, m: [1, 1, 5, 5]),
    SobolParameter(s: 5, a: 2, m: [1, 3, 3, 9, 7]),
    SobolParameter(s: 5, a: 4, m: [1, 1, 1, 15, 11]),
    SobolParameter(s: 5, a: 7, m: [1, 3, 5, 5, 31]),
    SobolParameter(s: 5, a: 11, m: [1, 1, 5, 5, 5]),
    SobolParameter(s: 5, a: 13, m: [1, 1, 5, 15, 41]),
    SobolParameter(s: 5, a: 14, m: [1, 3, 1, 17, 49]),
    SobolParameter(s: 5, a: 16, m: [1, 1, 1, 5, 5]),
    SobolParameter(s: 5, a: 19, m: [1, 3, 5, 15, 63]),
    SobolParameter(s: 5, a: 22, m: [1, 1, 3, 13, 27]),
    SobolParameter(s: 5, a: 25, m: [1, 3, 1, 15, 17])
]

#if swift(>=5.5)
private final class BlockingResultBox<Result>: @unchecked Sendable {
    var value: Result?
}
#endif

fileprivate func vanDerCorput<T: RealLike>(index: Int, base: Int) -> T {
    var n = index
    let b = base
    var f: T = 1
    var r: T = 0
    while n > 0 {
        f = f / T(b)
        r += f * T(n % b)
        n /= b
    }
    return r
}

fileprivate func firstPrimes(count: Int) -> [Int] {
    var primes: [Int] = []
    var x = 2
    while primes.count < count {
        var ok = true
        var d = 2
        while d * d <= x {
            if x % d == 0 { ok = false; break }
            d += 1
        }
        if ok { primes.append(x) }
        x += 1
    }
    return primes
}
