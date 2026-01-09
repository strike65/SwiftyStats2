//
//  Created by VT on 04.11.25.
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
//
//  GLMFit.swift
//  Minimal linear model fitting akin to R's lm.fit (OLS with pivoted QR).
//
//  Provides a least-squares solver for y = X * beta + eps.
//  Uses LAPACK pivoted QR (dgeqp3 + dormqr) when Accelerate is available,
//  falling back to normal equations with Gaussian elimination otherwise.

/// Lightweight generalized linear model fitting utilities and result containers.

import Foundation

#if canImport(Accelerate)
import Accelerate
#endif

/// Output from an OLS fit, mirroring R's `lm.fit` result fields.
public struct GLMFitResult: CustomStringConvertible {
    /// Estimated regression coefficients (length = column count of `x`).
    public let coefficients: [Double]
    /// Fitted values X * beta in the same order as the response.
    public let fitted: [Double]
    /// Residuals y - fitted.
    public let residuals: [Double]
    /// Effective rank after pivoted QR.
    public let rank: Int
    /// Residual degrees of freedom (m - rank).
    public let dfResidual: Int
    /// Pivot vector describing column permutation applied during factorization.
    public let pivot: [Int]

    /// Debug-friendly string representation summarizing key arrays.
    public var description: String {
        func preview<T>(_ array: [T], maxCount: Int = 6) -> String {
            func stringify(_ element: T) -> String {
                if let v = element as? Double { return niceNumber(v) }
                if let v = element as? Float { return niceNumber(v) }
                #if arch(i386) || arch(x86_64)
                if let v = element as? Float80 { return niceNumber(v) }
                #endif
                return String(describing: element)
            }
            if array.isEmpty { return "[]" }
            if array.count <= maxCount {
                return "[\(array.map { stringify($0) }.joined(separator: ", "))]"
            } else {
                let head = array.prefix(maxCount).map { stringify($0) }.joined(separator: ", ")
                return "[\(head), …] (count=\(array.count))"
            }
        }

        let coeffsStr = preview(coefficients)
        let fittedStr = preview(fitted)
        let residStr = preview(residuals)
        let pivotStr = preview(pivot)

        return """
        GLMFitResult:
          coefficients: \(coeffsStr)
          fitted: \(fittedStr)
          residuals: \(residStr)
          rank: \(rank)
          dfResidual: \(dfResidual)
          pivot: \(pivotStr)
        """
    }
}

/// Errors surfaced by the GLM fitting helpers.
public enum GLMFitError: Error, LocalizedError {
    case dimensionMismatch
    case singularMatrix

    /// User-presentable error descriptions.
    public var errorDescription: String? {
        switch self {
        case .dimensionMismatch:
            return "Dimension mismatch: design matrix rows must match response length and column counts must be consistent."
        case .singularMatrix:
            return "Singular system: design matrix is rank deficient or ill-conditioned; cannot compute unique OLS coefficients."
        }
    }
}

/// Namespace for least-squares fitting helpers.
public enum GLMFit {
    /// Solve min ||y - X beta||_2 via OLS.
    ///
    /// - Parameters:
    ///   - y: response vector of length m.
    ///   - x: design matrix m x n (row-major).
    /// - Returns: coefficients, fitted values, residuals, rank, and pivot.
    public static func fit(y: [Double], x: [[Double]]) throws -> GLMFitResult {
        let m = y.count
        guard m > 0 else { throw GLMFitError.dimensionMismatch }
        guard let n = x.first?.count, x.allSatisfy({ $0.count == n }), x.count == m else {
            throw GLMFitError.dimensionMismatch
        }

        #if canImport(Accelerate)
        return try fitAccelerate(y: y, x: x)
        #else
        return try fitNormalEquations(y: y, x: x)
        #endif
    }

    #if canImport(Accelerate)
    private static func fitAccelerate(y: [Double], x: [[Double]]) throws -> GLMFitResult {
        var m = Int32(y.count)
        var n = Int32(x.first?.count ?? 0)
        var a = toColumnMajor(x)
        var jpvt = [Int32](repeating: 0, count: Int(n))
        var tau = [Double](repeating: 0.0, count: Int(min(m, n)))

        // Workspace query for dgeqp3
        var lwork: Int32 = -1
        var workQuery: Double = 0.0
        var info: Int32 = 0
        var lda = m
        dgeqp3_(&m, &n, &a, &lda, &jpvt, &tau, &workQuery, &lwork, &info)
        lwork = Int32(workQuery)
        var work = [Double](repeating: 0.0, count: Int(max(1, lwork)))
        dgeqp3_(&m, &n, &a, &lda, &jpvt, &tau, &work, &lwork, &info)
        if info != 0 {
            throw GLMFitError.singularMatrix
        }

        // Apply Q^T to y
        var b = y
        var side: Int8 = 76 // 'L'
        var trans: Int8 = 84 // 'T'
        lwork = -1
        workQuery = 0.0
        var one: Int32 = 1
        var mm1 = m
        var nn1 = n
        var ldq = m
        dormqr_(&side, &trans, &mm1, &one, &nn1, &a, &lda, &tau, &b, &ldq, &workQuery, &lwork, &info)
        lwork = Int32(workQuery)
        work = [Double](repeating: 0.0, count: Int(max(1, lwork)))
        var mm2 = m
        var nn2 = n
        dormqr_(&side, &trans, &mm2, &one, &nn2, &a, &lda, &tau, &b, &ldq, &work, &lwork, &info)
        if info != 0 {
            throw GLMFitError.singularMatrix
        }

        // Solve R * beta = b[0..<n] with rank tolerance.
        var beta = [Double](repeating: 0.0, count: Int(n))
        var diagAbs: [Double] = []
        diagAbs.reserveCapacity(Int(min(m, n)))
        for i in 0..<Int(min(m, n)) {
            diagAbs.append(abs(a[i + i * Int(m)]))
        }
        let rmax = diagAbs.max() ?? 0.0
        let tolR = rmax * 1e-10
        var rank = 0
        for i in stride(from: Int(n) - 1, through: 0, by: -1) {
            let rii = a[i + i * Int(m)]
            if abs(rii) <= tolR {
                beta[i] = 0.0
                continue
            }
            rank += 1
            var sum = b[i]
            for j in (i + 1)..<Int(n) {
                sum -= a[i + j * Int(m)] * beta[j]
            }
            beta[i] = sum / rii
        }

        // Unpivot coefficients to original column order
        var coeffs = [Double](repeating: 0.0, count: Int(n))
        for (i, p) in jpvt.enumerated() {
            let idx = Int(p - 1)
            if idx >= 0 && idx < coeffs.count {
                coeffs[idx] = beta[i]
            }
        }

        let fitted = multiply(x: x, beta: coeffs)
        let residuals = zip(y, fitted).map { $0 - $1 }
        return GLMFitResult(
            coefficients: coeffs,
            fitted: fitted,
            residuals: residuals,
            rank: rank,
            dfResidual: Int(m) - rank,
            pivot: jpvt.map { Int($0) }
        )
    }
    #endif

    private static func fitNormalEquations(y: [Double], x: [[Double]]) throws -> GLMFitResult {
        let m = y.count
        let n = x.first?.count ?? 0
        // Compute XtX and XtY
        var xtx = Array(repeating: Array(repeating: 0.0, count: n), count: n)
        var xty = Array(repeating: 0.0, count: n)
        for r in 0..<m {
            for c in 0..<n {
                xty[c] += x[r][c] * y[r]
                for c2 in 0..<n {
                    xtx[c][c2] += x[r][c] * x[r][c2]
                }
            }
        }
        guard let inv = invert(matrix: xtx) else { throw GLMFitError.singularMatrix }
        var coeffs = Array(repeating: 0.0, count: n)
        for i in 0..<n {
            var sum = 0.0
            for j in 0..<n {
                sum += inv[i][j] * xty[j]
            }
            coeffs[i] = sum
        }
        let fitted = multiply(x: x, beta: coeffs)
        let residuals = zip(y, fitted).map { $0 - $1 }
        return GLMFitResult(
            coefficients: coeffs,
            fitted: fitted,
            residuals: residuals,
            rank: n,
            dfResidual: m - n,
            pivot: Array(0..<n)
        )
    }

    // MARK: - Helpers

    private static func toColumnMajor(_ x: [[Double]]) -> [Double] {
        let m = x.count
        let n = x.first?.count ?? 0
        var out = [Double](repeating: 0.0, count: m * n)
        for j in 0..<n {
            for i in 0..<m {
                out[j * m + i] = x[i][j]
            }
        }
        return out
    }

    private static func multiply(x: [[Double]], beta: [Double]) -> [Double] {
        let m = x.count
        let n = beta.count
        var out = [Double](repeating: 0.0, count: m)
        for i in 0..<m {
            var sum = 0.0
            for j in 0..<n {
                sum += x[i][j] * beta[j]
            }
            out[i] = sum
        }
        return out
    }

    private static func invert(matrix: [[Double]]) -> [[Double]]? {
        let n = matrix.count
        var a = matrix
        var inv = identity(n)
        for i in 0..<n {
            var pivot = a[i][i]
            if abs(pivot) < 1e-14 { return nil }
            pivot = 1.0 / pivot
            for j in 0..<n {
                a[i][j] *= pivot
                inv[i][j] *= pivot
            }
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

    #if canImport(Accelerate)
    private static func countRank(a: [Double], rows m: Int, cols n: Int) -> Int {
        let minmn = min(m, n)
        var rank = 0
        let tol = 1e-10
        for i in 0..<minmn {
            if abs(a[i + i * m]) > tol { rank += 1 }
        }
        return rank
    }
    #else
    private static func countRank(a: [Double], rows m: Int, cols n: Int) -> Int {
        return min(m, n)
    }
    #endif
}
