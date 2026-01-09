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

/// Shape statistics, outlier detection, and related helpers for SSExamine.

import SwiftyStatsPrelude
import SwiftyBoost

extension SSExamine {
    /// Fisher–Pearson skewness coefficient based on the sample standard deviation.
    ///
    /// - Returns: Skewness, or `nil` if the variance is undefined (e.g., fewer than two items),
    ///   the sample standard deviation is zero, or the dataset is non-numeric.
    ///
    /// Sign convention:
    /// - Negative values indicate left-skewness; positive values indicate right-skewness.
    public var skewness : FP? {
        return withLock { () -> FP? in
            guard let sd = self.standardDeviationLocked(type: .sample), !sd.isZero else {
                return self.logNil("skewness requires a non-zero sample standard deviation")
            }
            guard let m3 = empiricalMomentLocked(order: 3, kind: .central) else {
                return self.logNil("skewness requires the third central moment")
            }
            return m3 / FP.pow(sd, 3)
        }
    }
    

    /// Kurtosis (fourth standardized moment) based on the sample standard deviation.
    ///
    /// - Returns: Kurtosis, or `nil` if the variance is undefined, the standard deviation is zero,
    ///   or the dataset is non-numeric.
    ///
    /// Note:
    /// - For a normal distribution, kurtosis is 3. See `kurtosisExcess` for the zero-centered form.
    public var kurtosis: FP? {
        return withLock { () -> FP? in
            guard let sd = self.standardDeviationLocked(type: .sample), !sd.isZero else {
                return self.logNil("kurtosis requires a non-zero sample standard deviation")
            }
            guard let m4 = empiricalMomentLocked(order: 4, kind: .central) else {
                return self.logNil("kurtosis requires the fourth central moment")
            }
            return m4 / FP.pow(sd, 4)
        }
    }
    

    /// Excess kurtosis (kurtosis − 3), yielding 0 for the normal distribution baseline.
    ///
    /// - Returns: Excess kurtosis, or `nil` if variance is undefined or the dataset is non-numeric.
    public var kurtosisExcess: FP? {
        return withLock { () -> FP? in
            guard let sd = self.standardDeviationLocked(type: .sample), !sd.isZero else {
                return self.logNil("kurtosisExcess requires a non-zero sample standard deviation")
            }
            guard let m4 = empiricalMomentLocked(order: 4, kind: .central) else {
                return self.logNil("kurtosisExcess requires the fourth central moment")
            }
            return m4 / FP.pow(sd, 4) - 3
        }
    }
    

    /// Indicates whether the distribution is platycurtic (flatter than normal).
    ///
    /// Returns:
    /// - `nil` if `kurtosisExcess` is unavailable; otherwise `true` if excess kurtosis < 0.
    public var isPlatyurtic: Bool? {
        if let k = self.kurtosisExcess {
            return k < FP.zero
        }
        else {
            return self.logNil("isPlatyurtic requires kurtosisExcess")
        }
    }
    

    /// Indicates whether the distribution is mesokurtic (normal-like kurtosis).
    ///
    /// Returns:
    /// - `nil` if `kurtosisExcess` is unavailable; otherwise `true` if excess kurtosis == 0.
    public var isMesokurtic: Bool? {
        if let k = self.kurtosisExcess {
            return k.isZero
        }
        else {
            return self.logNil("isMesokurtic requires kurtosisExcess")
        }
    }
    

    /// Indicates whether the distribution is leptokurtic (more peaked than normal).
    ///
    /// Returns:
    /// - `nil` if `kurtosisExcess` is unavailable; otherwise `true` if excess kurtosis > 0.
    public var isLeptokurtic: Bool? {
        if let k = self.kurtosisExcess {
            return k > FP.zero
        }
        else {
            return self.logNil("isLeptokurtic requires kurtosisExcess")
        }
    }
    

    /// Indicates whether the distribution is left-skewed (negative skewness).
    ///
    /// Returns:
    /// - `nil` if `skewness` is unavailable; otherwise `true` if skewness < 0.
    public var isLeftSkewed: Bool? {
        if let k = self.skewness {
            return k < FP.zero
        }
        else {
            return self.logNil("isLeftSkewed requires skewness")
        }
    }
    

    /// Indicates whether the distribution is right-skewed (positive skewness).
    ///
    /// Returns:
    /// - `nil` if `skewness` is unavailable; otherwise `true` if skewness > 0.
    public var isRightSkewed: Bool? {
        if let k = self.skewness {
            return k > FP.zero
        }
        else {
            return self.logNil("isRightSkewed requires skewness")
        }
    }
    

    /// Indicates whether the distribution is symmetric (zero skewness).
    ///
    /// Returns:
    /// - `nil` if `skewness` is unavailable; otherwise `true` if skewness == 0.
    public var isSymmetric: Bool? {
        if let k = self.skewness {
            return k.isZero
        }
        else {
            return self.logNil("isSymmetric requires skewness")
        }
    }
    

    /// Moors' kurtosis based on octiles (quantile-based kurtosis measure).
    ///
    /// Definition:
    /// - Uses octiles Q0, Q1, …, Q8 (with Q0 = min, Q8 = max) and is defined as
    ///   (Q7 − Q5 + Q3 − Q1) / (Q6 − Q2), when all quantiles are defined.
    ///
    /// Returns:
    /// - `nil` for empty or non-numeric datasets or when required quantiles cannot be computed.
    ///
    /// Thread-safety:
    /// - Acquires the instance lock and delegates to `moorsKurtosisLocked()`.
    ///
    /// Complexity:
    /// - O(n log n) due to quantile computations.
    public var moorsKurtosis: FP? {
        return withLock { () -> FP? in
            self.moorsKurtosisLocked()
        }
    }
    
}
