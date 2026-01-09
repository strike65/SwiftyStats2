//
//  Created by VT on 28.10.25.
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

/// Public concentration and inequality measures computed from SSExamine samples.

import SwiftyStatsPrelude
import SwiftyBoost

extension SSExamine {
    /// Gini's mean difference 2 · mean · G · n / (n - 1)
    ///
    /// - Returns: Gini's mean difference or `nil` if G or mean not available
    public var giniMeanDifference: FP? {
        return withLock { () -> FP? in
            return self.giniMeanDifferenceLocked()
        }
    }

    /// Gini coefficient measuring inequality of the sample distribution.
    ///
    /// Definition:
    /// - G = (1 / (n² · mean)) Σ_{i=1..n} (2i − n − 1) x_{(i)}, with x_{(i)} sorted ascending.
    ///
    /// Returns:
    /// - `nil` when the dataset is empty, contains negative values, or cannot be projected to numeric form.
    /// - Otherwise a value in [0, 1] where 0 indicates perfect equality and values near 1 indicate maximal concentration.
    public var giniCoefficient: FP? {
        get {
            return withLock { () -> FP? in
                return giniCoefficientLocked()
            }
        }
    }
    

    /// Normalized Gini coefficient `G_n = G · n / (n − 1)` to bound results in [0, 1].
    ///
    /// - Returns: Rescales the raw Gini so a perfectly unequal sample (all mass on one item)
    ///   maps to 1. Returns `nil` when the base Gini cannot be computed.
    public var giniNorm: FP? {
        get {
            return withLock { () -> FP? in
                return giniNormLocked()
            }
        }
    }

    /// Herfindahl–Hirschman index (HHI) computed on the sample.
    ///
    /// - Returns: Σ (xᵢ / Σ xⱼ)² for positive numeric values, or `nil` if values are invalid
    ///   or totals cannot be established.
    public var hhi : FP? {
        return withLock { () -> FP? in
            return hhiLocked()
        }
    }

    /// Concentration ratio CR\_k: sum of the k largest values.
    ///
    /// - Parameter k: Number of top observations to include (1 ≤ k ≤ n).
    /// - Returns: Σ_{i=1..k} x\_{(i)} where x\_{(i)} are order statistics sorted descending,
    ///   or `nil` when inputs are invalid or non-numeric.
    ///
    /// Throws:
    /// - `SSError(.invalidArgument, ...)` if `k` is outside [1, n].
    public func cr(k: Int) throws -> FP? {
        return try withLock { () -> FP? in
            return try crLocked(k: k)
        }
    }
    

    /// Alias for `cr(k:)`.
    ///
    /// Provided for readability in API consumers; see `cr(k:)` for mathematical definition and return semantics.
    ///
    /// Throws:
    /// - `SSError(.invalidArgument, ...)` if `k` is outside [1, n].
    public func concentrationRatio(k: Int) throws -> FP? {
        return try withLock { () -> FP? in
            return try crLocked(k: k)
        }
    }
    

    /// Computes the Lorenz curve points (cumulative unit versus cumulative value shares).
    /// - Returns: A sequence of `(u, v)` pairs where `u` is the population share and
    ///   `v` is the cumulative share of the observed values.
    /// - Note: Requires a non-empty, numeric, non-negative dataset.
    public var lorenzCurve: [(u: FP, v: FP)]? {
        get {
            return withLock{ () -> [(u: FP, v: FP)]? in
                return lorenzCurveLocked
            }
        }
    }
    
}
