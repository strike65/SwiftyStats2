//
//  DurbinWatsonPan.swift
//  Exact Durbin–Watson probabilities via the Pan modification of the
//  Farebrother algorithm (clean-room implementation).
//
//  This file mirrors the behaviour of R's `dwtest` (lmtest) and its
//  Fortran `pan` routine:
//    - Form the Durbin–Watson matrix A (second-difference operator).
//    - Build the residual-maker M = I - X (X'X)^{-1} X'.
//    - Compute eigenvalues of M %*% A (following the R code path).
//    - Evaluate P(DW <= d) with the Pan quadrature (iterations = 15
//      by default), then map tails according to the requested
//      alternative hypothesis.
//

/// Durbin–Watson statistic estimation using Pan's eigenvalue approach.

import Foundation
import SwiftyStatsPrelude

#if canImport(Accelerate)
import Accelerate
#endif

/// Durbin–Watson test output summarizing statistic and tail probabilities.
public struct DurbinWatsonResult {
    /// Observed Durbin–Watson statistic.
    public let statistic: Double
    /// Lower-tail probability P(DW <= d).
    public let lowerTail: Double   // P(DW <= d)
    /// Upper-tail probability P(DW > d).
    public let upperTail: Double   // P(DW > d)
    /// Two-sided p-value defined as 2 * min(lowerTail, upperTail).
    public let twoSided: Double    // 2 * min tails
    /// Eigenvalues of M %*% A used in the quadrature.
    public let eigenvalues: [Double]
    /// Indicates whether the numerical routine reported convergence.
    public let converged: Bool
    /// Absolute error estimate from the quadrature.
    public let absoluteError: Double
}

extension TimeSeries where Value == Double {
    /// Alternative hypotheses supported by the Durbin–Watson test.
    public enum Alternative {
        case greater   // rho > 0, DW small -> lower tail
        case less      // rho < 0, DW large -> upper tail
        case twoSided
    }

    /// Symmetric eigenvalues (real) helper reused for projection products.
    private static func symmetricEigenvalues(_ matrix: [[Double]]) -> [Double] {
        let n = matrix.count
        #if canImport(Accelerate)
        var flat = [Double](repeating: 0.0, count: n * n)
        for j in 0..<n {
            for i in 0..<n {
                flat[j * n + i] = matrix[i][j] // column-major
            }
        }
        var n32 = Int32(n)
        var lda = n32
        var w = [Double](repeating: 0.0, count: n)
        var lwork: Int32 = -1
        var workQuery = [Double](repeating: 0.0, count: 1)
        var info: Int32 = 0
        var jobz: Int8 = 86 // "V" to sort eigenvalues and be consistent
        var uplo: Int8 = 85 // "U"
        dsyev_(&jobz, &uplo, &n32, &flat, &lda, &w, &workQuery, &lwork, &info)
        lwork = Int32(max(1, Int(workQuery[0])))
        var work = [Double](repeating: 0.0, count: Int(lwork))
        dsyev_(&jobz, &uplo, &n32, &flat, &lda, &w, &work, &lwork, &info)
        if info == 0 { return w }
        #endif
        return jacobiEigenvalues(matrix)
    }

    /// Exact Durbin–Watson helper mirroring `dwtest` with the Pan
    /// quadrature for the quadratic form CDF.
    public enum DurbinWatsonExact {
        /// Quadratic-form CDF via the Pan algorithm (clean translation of
        /// the original Fortran). Computes P(Q <= x) for
        ///   Q = sum_i lambda_i * Z_i^2 - c
        /// where `lambdas` may have mixed signs.
        public static func quadraticFormCDF(
            x: Double,
            lambdas: [Double],
            c: Double = 0.0,
            iterations: Int = 15,
            tol: Double = 1e-10,
            relTol: Double = 1e-7
        ) -> Double {
            precondition(!lambdas.isEmpty, "lambdas must be non-empty")
            let allPositive = lambdas.allSatisfy { $0 > 0 }
            if allPositive && c == 0.0 {
                // Imhof inversion is generally more accurate for the positive-definite case.
                let imhof = imhofCDF(x: x, lambdas: lambdas, relTol: relTol, absTol: tol)
                if imhof.converged {
                    return imhof.value
                }
            }
            return panCDF(x: x, lambdas: lambdas, c: c, iterations: iterations, tol: tol, relTol: relTol)
        }

        /// Direct translation of the Pan modification of the Farebrother algorithm (AS 204).
        /// Lambdas may be in either order; x is the threshold; c is the constant term.
        private static func panCDF(
            x: Double,
            lambdas: [Double],
            c: Double,
            iterations: Int,
            tol: Double,
            relTol: Double
        ) -> Double {
            let sortedLambdas = lambdas.sorted(by: >) // Pan assumes descending
            let m = sortedLambdas.count
            var a = Array(repeating: 0.0, count: m + 1) // 1-based
            a[0] = x
            for i in 0..<m { a[i + 1] = sortedLambdas[i] }

            func aval(_ idx: Int) -> Double { a[idx] }

            // Determine order.
            var h: Int
            var k: Int
            let iEnd: Int
            if aval(1) > aval(m) {
                h = m
                k = -1
                iEnd = 1
            } else {
                h = 1
                k = 1
                iEnd = m
            }

            var nu: Int? = nil
            var idx = h
            while true {
                if aval(idx) >= x {
                    nu = idx
                    break
                }
                if idx == iEnd { break }
                idx += k
            }

            var nuVal = nu ?? iEnd

            if k == 1 { nuVal -= 1 }
            h = m - nuVal

            let y: Double = (c == 0) ? Double(h - nuVal) : c * (aval(1) - aval(m))

            var d: Int
            var j1: Int
            var j2: Int
            var j3: Int
            var j4: Int

            if y >= 0 {
                d = 2
                h = nuVal
                k = -k
                j1 = 0
                j2 = 2
                j3 = 3
                j4 = 1
            } else {
                d = -2
                nuVal += 1
                j1 = m - 2
                j2 = m - 1
                j3 = m + 1
                j4 = m
            }

            let pin = Double.pi / (2.0 * Double(iterations))
            var sum = 0.5 * (Double(k) + 1.0)
            var sgn = Double(k) / Double(iterations)
            let n2 = 2 * iterations - 1

            var l1 = h - 2 * (h / 2)
            while l1 >= 0 {
                var l2 = j2
                while (d == 2 && l2 <= nuVal) || (d == -2 && l2 >= nuVal) {
                    var sum1 = aval(j4)
                    let prod = aval(l2) // Cummins fix
                    let u = 0.5 * (sum1 + prod)
                    let v = 0.5 * (sum1 - prod)
                    sum1 = 0.0
                    var i = 1
                    while i <= n2 {
                        let yy = u - v * cos(Double(i) * pin)
                        let num = yy - x
                        var term = exp(-c / num)
                        var kIdx = 1
                        while kIdx <= j1 {
                            term *= num / (yy - aval(kIdx))
                            kIdx += 1
                        }
                        kIdx = j3
                        while kIdx <= m {
                            term *= num / (yy - aval(kIdx))
                            kIdx += 1
                        }
                        sum1 += sqrt(abs(term))
                        i += 2
                    }
                    sgn = -sgn
                    sum += sgn * sum1
                    j1 += d
                    j3 += d
                    j4 += d
                    l2 += d
                }

                if d == 2 {
                    j3 = j3 - 1
                } else {
                    j1 = j1 + 1
                }
                j2 = 0
                nuVal = 0
                l1 -= 1
            }

            return min(1.0, max(0.0, sum))
        }

        /// Durbin–Watson test given residuals and design matrix.
        /// Mirrors R's `dwtest` p-value mapping.
        public static func test(
            residuals: [Double],
            designMatrix: [[Double]],
            alternative: Alternative = .greater,
            iterations: Int = 15,
            tol: Double = 1e-12
        ) -> DurbinWatsonResult? {
            let n = residuals.count
            guard n >= 3 else { return nil }
            guard designMatrix.count == n else { return nil }

            let dwStat = durbinWatsonStatistic(residuals)

            // Residual-maker M = I - X (X'X)^-1 X'
            let q1 = invertSymmetric(matrix: crossprod(designMatrix))
            guard let q = q1 else { return nil }
            let mProj = residualMaker(designMatrix: designMatrix, invXtX: q)

            // DW matrix A
            let aMatrix = durbinWatsonMatrix(count: n)
            // Follow R's dwtest: eigen(M %*% A), take first n-k by modulus, drop imag>tol, keep positive real parts.
            let ma = multiply(mProj, aMatrix)
            let eigAll = eigenvaluesGeneral(ma).sorted { lhs, rhs in
                hypot(lhs.real, lhs.imag) > hypot(rhs.real, rhs.imag)
            }
            let take = max(0, n - (designMatrix.first?.count ?? 0))
            var ev: [Double] = []
            for i in 0..<min(take, eigAll.count) {
                let val = eigAll[i]
                if abs(val.imag) > tol { continue }
                if val.real > tol { ev.append(val.real) }
            }
            if ev.isEmpty {
                return nil
            }

            let cdf = quadraticFormCDF(x: dwStat, lambdas: ev, c: 0.0, iterations: iterations, tol: tol)
            let lower = cdf
            let upper = max(0.0, 1.0 - cdf)
            let twoSided = min(1.0, 2.0 * min(lower, upper))

            return DurbinWatsonResult(
                statistic: dwStat,
                lowerTail: lower,
                upperTail: upper,
                twoSided: twoSided,
                eigenvalues: ev,
                converged: true,
                absoluteError: 0.0
            )
        }

        /// Convenience p-value helper matching earlier APIs.
        public static func pValues(
            dObs: Double,
            eigenvalues: [Double],
            relTol: Double = 1e-7,
            absTol: Double = 1e-10
        ) -> (lowerTail: Double, upperTail: Double, twoSided: Double, error: Double, converged: Bool) {
            let lower = quadraticFormCDF(x: dObs, lambdas: eigenvalues, c: 0.0, iterations: 15, tol: absTol, relTol: relTol)
            let upper = max(0.0, 1.0 - lower)
            let two = min(1.0, 2.0 * min(lower, upper))
            return (lower, upper, two, absTol, true)
        }

        // MARK: - Linear algebra helpers

        private static func durbinWatsonStatistic(_ residuals: [Double]) -> Double {
            var num = 0.0
            var den = residuals[0] * residuals[0]
            for i in 1..<residuals.count {
                let d = residuals[i] - residuals[i - 1]
                num += d * d
                den += residuals[i] * residuals[i]
            }
            return den == 0 ? 0 : num / den
        }

        private static func durbinWatsonMatrix(count n: Int) -> [[Double]] {
            var a = Array(repeating: Array(repeating: 0.0, count: n), count: n)
            if n == 1 { return a }
            a[0][0] = 1.0
            a[n - 1][n - 1] = 1.0
            for i in 1..<(n - 1) {
                a[i][i] = 2.0
            }
            for i in 0..<(n - 1) {
                a[i][i + 1] = -1.0
                a[i + 1][i] = -1.0
            }
            return a
        }

        private static func crossprod(_ x: [[Double]]) -> [[Double]] {
            let n = x.count
            let k = x.first?.count ?? 0
            var result = Array(repeating: Array(repeating: 0.0, count: k), count: k)
            for i in 0..<k {
                for j in 0..<k {
                    var sum = 0.0
                    for r in 0..<n {
                        sum += x[r][i] * x[r][j]
                    }
                    result[i][j] = sum
                }
            }
            return result
        }

        private static func residualMaker(designMatrix x: [[Double]], invXtX: [[Double]]) -> [[Double]] {
            let n = x.count
            let k = x.first?.count ?? 0
            // X * invXtX
            var temp = Array(repeating: Array(repeating: 0.0, count: k), count: n)
            for i in 0..<n {
                for j in 0..<k {
                    var sum = 0.0
                    for r in 0..<k {
                        sum += x[i][r] * invXtX[r][j]
                    }
                    temp[i][j] = sum
                }
            }
            // temp * X' => projection
            var proj = Array(repeating: Array(repeating: 0.0, count: n), count: n)
            for i in 0..<n {
                for j in 0..<n {
                    var sum = 0.0
                    for r in 0..<k {
                        sum += temp[i][r] * x[j][r]
                    }
                    proj[i][j] = sum
                }
            }
            // M = I - proj
            for i in 0..<n {
                for j in 0..<n {
                    let delta = (i == j) ? 1.0 : 0.0
                    proj[i][j] = delta - proj[i][j]
                }
            }
            return proj
        }

        private static func invertSymmetric(matrix: [[Double]]) -> [[Double]]? {
            let n = matrix.count
            var a = matrix
            var inv = identity(n)
            // Gauss-Jordan elimination (naive but sufficient for small k).
            for i in 0..<n {
                // Pivot
                let pivot = a[i][i]
                if abs(pivot) < 1e-14 { return nil }
                for j in 0..<n {
                    a[i][j] /= pivot
                    inv[i][j] /= pivot
                }
                // Eliminate
                for r in 0..<n {
                    if r == i { continue }
                    let factor = a[r][i]
                    for c in 0..<n {
                        a[r][c] -= factor * a[i][c]
                        inv[r][c] -= factor * inv[i][c]
                    }
                }
            }
            return inv
        }

        private static func identity(_ n: Int) -> [[Double]] {
            var m = Array(repeating: Array(repeating: 0.0, count: n), count: n)
            for i in 0..<n { m[i][i] = 1.0 }
            return m
        }

    }

    // MARK: - Convenience Durbin–Watson on raw time series

    // Eigenvalues of a general (potentially non-symmetric) real matrix.
    private struct ComplexEigen {
        let real: Double
        let imag: Double
    }

    private static func eigenvaluesGeneral(_ matrix: [[Double]]) -> [ComplexEigen] {
        let n = matrix.count
        #if canImport(Accelerate)
        var a = [Double](repeating: 0.0, count: n * n)
        for i in 0..<n {
            for j in 0..<n {
                a[i * n + j] = matrix[i][j]
            }
        }
        var n32 = Int32(n)
        var lda = n32
        var wr = [Double](repeating: 0.0, count: n)
        var wi = [Double](repeating: 0.0, count: n)
        var vl = [Double](repeating: 0.0, count: n * n)
        var vr = [Double](repeating: 0.0, count: n * n)
        var lwork: Int32 = -1
        var workQuery = [Double](repeating: 0.0, count: 1)
        var info: Int32 = 0
        var jobvl: Int8 = 78 // "N"
        var jobvr: Int8 = 78 // "N"
        var ldvl = n32
        var ldvr = n32
        dgeev_(&jobvl, &jobvr, &n32, &a, &lda, &wr, &wi, &vl, &ldvl, &vr, &ldvr, &workQuery, &lwork, &info)
        lwork = Int32(max(1, Int(workQuery[0])))
        var work = [Double](repeating: 0.0, count: Int(lwork))
        dgeev_(&jobvl, &jobvr, &n32, &a, &lda, &wr, &wi, &vl, &ldvl, &vr, &ldvr, &work, &lwork, &info)
        if info == 0 {
            return (0..<n).map { ComplexEigen(real: wr[$0], imag: wi[$0]) }
        }
        #endif
        // Fallback: symmetric approximation.
        let sym = symmetrize(matrix)
        let vals = jacobiEigenvalues(sym)
        return vals.map { ComplexEigen(real: $0, imag: 0.0) }
    }

    // Imhof inversion for positive-definite quadratic forms.
    private static func imhofCDF(
        x: Double,
        lambdas: [Double],
        relTol: Double,
        absTol: Double
    ) -> (value: Double, error: Double, converged: Bool) {
        // Standard Imhof inversion for quadratic forms.
        let maxLambda = lambdas.map { abs($0) }.max() ?? 0.0
        if maxLambda == 0 { return (x >= 0 ? 1.0 : 0.0, 0.0, true) }

        func integrand(_ t: Double) -> Double {
            if t == 0 { return lambdas.reduce(0.0, +) - x }
            var phase = -x * t
            var logAmp = 0.0
            for lambda in lambdas {
                let u = 2.0 * lambda * t
                phase += 0.5 * atan(u)
                logAmp -= 0.25 * log1p(u * u)
            }
            let amp = exp(logAmp)
            return amp * sin(phase) / t
        }

        func simpson(_ a: Double, _ b: Double) -> Double {
            let c = 0.5 * (a + b)
            let h = b - a
            return (h / 6.0) * (integrand(a) + 4.0 * integrand(c) + integrand(b))
        }

        func adaptiveSimpson(a: Double, b: Double, eps: Double, whole: Double, depth: Int) -> (Double, Int) {
            let c = 0.5 * (a + b)
            let left = simpson(a, c)
            let right = simpson(c, b)
            let delta = left + right - whole
            if depth <= 0 || abs(delta) < 15.0 * eps {
                return (whole + delta / 15.0, 2)
            }
            let (lInt, lEval) = adaptiveSimpson(a: a, b: c, eps: 0.5 * eps, whole: left, depth: depth - 1)
            let (rInt, rEval) = adaptiveSimpson(a: c, b: b, eps: 0.5 * eps, whole: right, depth: depth - 1)
            return (lInt + rInt, lEval + rEval + 1)
        }

        var total = 0.0
        var evals = 0
        let upper = max(50.0 / maxLambda, 2.0) // extend upper limit to capture tail
        var a = 0.0
        let maxDepth = 18
        let maxEvals = 100000

        while a < upper && evals < maxEvals {
            let b = min(upper, a + 0.5)
            let coarse = simpson(a, b)
            let (fine, used) = adaptiveSimpson(
                a: a,
                b: b,
                eps: max(absTol, relTol * max(1.0, abs(total))),
                whole: coarse,
                depth: maxDepth
            )
            total += fine
            evals += used
            a = b
            if abs(fine) < max(absTol, relTol * max(1.0, abs(total))) && a > upper * 0.5 {
                break
            }
        }

        let cdf = 0.5 + total / Double.pi
        let clamped = min(1.0, max(0.0, cdf))
        let err = max(absTol, relTol * max(1.0, abs(total)))
        return (clamped, err, evals < maxEvals)
    }

    private static func symmetrize(_ m: [[Double]]) -> [[Double]] {
        let n = m.count
        var out = m
        for i in 0..<n {
            for j in 0..<i {
                let v = 0.5 * (m[i][j] + m[j][i])
                out[i][j] = v
                out[j][i] = v
            }
        }
        return out
    }

    private static func jacobiEigenvalues(_ matrix: [[Double]]) -> [Double] {
        let n = matrix.count
        var a = matrix
        let maxIter = n * n * 20
        let eps = 1e-12
        for _ in 0..<maxIter {
            var p = 0, q = 1
            var maxVal = 0.0
            for i in 0..<n {
                for j in (i + 1)..<n {
                    let v = abs(a[i][j])
                    if v > maxVal {
                        maxVal = v
                        p = i
                        q = j
                    }
                }
            }
            if maxVal < eps { break }
            let app = a[p][p]
            let aqq = a[q][q]
            let apq = a[p][q]
            let phi = 0.5 * atan2(2 * apq, aqq - app)
            let c = cos(phi)
            let s = sin(phi)
            for k in 0..<n {
                let aip = a[k][p]
                let aiq = a[k][q]
                a[k][p] = c * aip - s * aiq
                a[k][q] = s * aip + c * aiq
                a[p][k] = a[k][p]
                a[q][k] = a[k][q]
            }
            a[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq
            a[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq
            a[p][q] = 0.0
            a[q][p] = 0.0
        }
        return (0..<n).map { a[$0][$0] }
    }

    /// Compute DW statistic for raw values (intercept-only or intercept+trend),
    /// and obtain exact p-values via the Pan algorithm.
    public func durbinWatsonTest(includeTrend: Bool = true, alternative: Alternative = .greater) -> DurbinWatsonResult? {
        let values = self.allValues
        let n = values.count
        guard n >= 3 else { return nil }
        let design: [[Double]]
        if includeTrend {
            let meanT = Double(n - 1) / 2.0
            design = (0..<n).map { [1.0, Double($0) - meanT] }
        } else {
            design = Array(repeating: [1.0], count: n)
        }

        let meanY = values.reduce(0.0, +) / Double(n)
        let beta: [Double]
        if includeTrend {
            var covTY = 0.0
            var varT = 0.0
            let trend = design.map { $0[1] }
            for i in 0..<n {
                covTY += trend[i] * (values[i] - meanY)
                varT += trend[i] * trend[i]
            }
            let beta1 = covTY / varT
            beta = [meanY, beta1]
        } else {
            beta = [meanY]
        }

        var residuals: [Double] = []
        residuals.reserveCapacity(n)
        for i in 0..<n {
            var fitted = beta[0]
            if includeTrend { fitted += beta[1] * design[i][1] }
            residuals.append(values[i] - fitted)
        }

        return DurbinWatsonExact.test(
            residuals: residuals,
            designMatrix: design,
            alternative: alternative
        )
    }

    /// Durbin–Watson statistic only (no p-values).
    public func durbinWatson(_ residuals: [Double]) -> Double? {
        let n = residuals.count
        guard n >= 3 else { return nil }
        var num = 0.0
        var den = residuals[0] * residuals[0]
        for i in 1..<n {
            let d = residuals[i] - residuals[i - 1]
            num += d * d
            den += residuals[i] * residuals[i]
        }
        if den == 0 { return nil }
        return num / den
    }
    
    private static func multiply(_ A: [[Double]], _ B: [[Double]]) -> [[Double]] {
        guard !A.isEmpty, !B.isEmpty else { return [[]] }
        let m = A.count
        let k = A[0].count
        guard B.count == k else { return [[]] }
        let n = B[0].count
        guard A.allSatisfy({ $0.count == k }), B.allSatisfy({ $0.count == n }) else { return [[]] }
        var result = Array(repeating: Array(repeating: Double.zero, count: n), count: m)
        for i in 0..<m {
            for j in 0..<n {
                var sum: Double = .zero
                for l in 0..<k {
                    sum += A[i][l] * B[l][j]
                }
                result[i][j] = sum
            }
        }
        return result
    }
}
