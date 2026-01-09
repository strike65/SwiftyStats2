//
//  Created by VT on 29.11.25.
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

/// Variance test wrappers for chi-square and F-based procedures.

import SwiftyStatsPrelude

extension Inferential.HypothesisTesting.VarianceTests {
    
    /// Bartlett test for equality of variances across `k` groups.
    ///
    /// The test statistic is
    /// ````text
    ///   T = ((N - k) * ln(s_p2) - sum_i ( (n_i - 1) * ln(s_i2) ))
    ///       / (1 + (1 / (3 * (k - 1))) * (sum_i 1/(n_i - 1) - 1/(N - k)))
    /// ````
    /// where
    /// ````text
    ///   s_p2 = sum_i ( (n_i - 1) * s_i2 ) / (N - k)
    /// ````
    /// and `s_i2` is the sample variance of group `i` with `n_i` observations.
    /// Under the null (equal variances and approximate normality), `T`
    /// follows `chiSquare(df: k - 1)`.
    ///
    /// - Parameters:
    ///   - samples: Columns of a `DataFrame` representing groups.
    ///   - alpha: Significance level in `(0, 1)`.
    /// - Returns: `VarianceEqualityTestResult` populated with cutoffs and p-value.
    /// - Throws: `SSError` when variances cannot be computed.
    public static func bartlettTest(samples: DataFrame<T, T>, alpha: T) throws -> VarianceEqualityTestResult<T>? {
        precondition(!samples.isEmpty, "bartlettTest: no data")
        precondition(samples.nCols >= 2, "Bartlett test requires at least two groups")
        let totalSampleSize = T(samples.totalSampleSize)
        let groupCount = T(samples.nCols)
        let pooledDegreesOfFreedom = totalSampleSize - groupCount
        var pooledVarianceNumerator: T = 0
        var weightedLogVariance: T = 0
        var inverseDegreesSum: T = 0
        for sample in samples.data {
            guard let variance = sample.sampleVariance else {
                SSLog.statisticsError("At least for one sample variance is not available.")
                throw SSError(type: .missingData, file: #fileID, line: #line, function: #function, reason: "At least for one sample variance is not available.")
            }
            let degrees = T(sample.sampleSize) - T.one
            pooledVarianceNumerator += degrees * variance
            weightedLogVariance += degrees * T.log(variance)
            inverseDegreesSum += T.one / degrees
        }
        let pooledVariance = pooledVarianceNumerator / pooledDegreesOfFreedom
        let numerator = (pooledDegreesOfFreedom * T.log(pooledVariance)) - weightedLogVariance
        let correction = T.one + (T.one / (T(3) * (groupCount - T.one))) * (inverseDegreesSum - T.one / pooledDegreesOfFreedom)
        let testStatisticValue = numerator / correction
        do {
            let dist: Distribution.ChiSquared<T> = try .init(degreesOfFreedom: groupCount - T.one)
            let cdfChiSquare = try dist.cdf(testStatisticValue)
            let cutoff90Percent = try dist.quantile(T(0.9))
            let cutoff95Percent = try dist.quantile(T(0.95))
            let cutoff99Percent = try dist.quantile(T(0.99))
            let cutoffAlpha = try dist.quantile(T.one - alpha)
            let df = groupCount - T.one
            var result = VarianceEqualityTestResult(df: nil,
                                                    cv90Pct: cutoff90Percent,
                                                    cv95Pct: cutoff95Percent,
                                                    cv99Pct: cutoff99Percent,
                                                    cvAlpha: cutoffAlpha,
                                                    pValue: T.one - cdfChiSquare,
                                                    testStatistic: testStatisticValue,
                                                    equality: !(cdfChiSquare > (T.one - alpha)),
                                                    testType: .bartlett)
            result.df = VarianceEqualityTestResult.DegreesOfFreedom(df1: df, df2: .nan)
            return result
        }
        catch {
            SSLog.statisticsError("Test not performed due to missing data.")
            throw SSError(type: .missingData, file: #fileID, line: #line, function: #function, reason: "Test not performed due to missing data.")
        }
        
    }
    
    /// Bartlett test convenience overload for an array of `SSExamine` groups.
    ///
    /// Each element is appended as one group; see ``bartlettTest(samples:alpha:)``
    /// for the test definition.
    ///
    /// - Parameters:
    ///   - samples: Collection of grouped observations.
    ///   - alpha: Significance level in `(0, 1)`.
    /// - Returns: `VarianceEqualityTestResult` populated with cutoffs and p-value.
    /// - Throws: `SSError` when variances cannot be computed.
    public static func bartlettTest(samples: [SSExamine<T, T>], alpha: T) throws -> VarianceEqualityTestResult<T>? {
        precondition(samples.count > 0, "No data")
        let dataFrame: DataFrame<T, T> = try DataFrame<T,T>(data: samples)
        return try bartlettTest(samples: dataFrame, alpha: alpha)
    }
    
    /// Levene test for equality of variances across `k` groups.
    ///
    /// Each observation `y_ij` is transformed to its absolute deviation from a
    /// group center `c_i` chosen via ``LeveneTestType``:
    /// ````text
    ///   z_ij = |y_ij - c_i|
    /// ````
    /// The test statistic is
    /// ````text
    ///   W = ((N - k) / (k - 1))
    ///       * (sum_i n_i (z_bar_i - z_bar)^2) / (sum_i sum_j (z_ij - z_bar_i)^2)
    /// ````
    /// where `n_i` is the sample size of group `i`, `z_bar_i` the mean absolute
    /// deviation within group `i`, and `z_bar` the grand mean of all `z_ij`.
    /// Under the null hypothesis of equal variances, `W` follows
    /// `F(k - 1, N - k)`.
    ///
    /// - Parameters:
    ///   - samples: Columns of a `DataFrame`, each representing one group.
    ///   - alpha: Significance level in `(0, 1)`.
    ///   - testType: Centering choice (`mean`, `median`, or 10% `trimmedMean`).
    /// - Returns: `VarianceEqualityTestResult` with critical values and p-value.
    /// - Throws: `SSError` when required statistics cannot be computed.
    public static func leveneTest(samples: DataFrame<T, T>, alpha: T, testType: LeveneTestType) throws -> VarianceEqualityTestResult<T>? {
        precondition(!samples.isEmpty, "No data")
        precondition(samples.nCols >= 2, "Levene test requires at least two groups")
        var zMeans: [T] = []
        var groupCounts: [T] = []
        var withinGroupSquares: [T] = []
        var totalCount: T = T.zero
        var totalAbsoluteDeviation: T = T.zero
        for sample in samples.data {
            let center: T
            switch testType {
            case .mean:
                if let m = sample.arithmeticMean {
                    center = m
                }
                else {
                    SSLog.statisticsError("At least for one sample-mean is not available.")
                    throw SSError(type: .missingData, file: #fileID, line: #line, function: #function, reason: "At least for one sample-mean is not available.")
                }
            case .median:
                if let m = sample.median {
                    center = m
                }
                else {
                    SSLog.statisticsError("At least for one sample-median is not available.")
                    throw SSError(type: .missingData, file: #fileID, line: #line, function: #function, reason: "At least for one sample-median is not available.")
                }
            case .trimmedMean:
                if let m = sample.trimmedMean(alpha: 0.1) {
                    center = m
                }
                else {
                    SSLog.statisticsError("At least for one sample-mean is not available.")
                    throw SSError(type: .missingData, file: #fileID, line: #line, function: #function, reason: "At least for one sample-mean is not available.")
                }
            }
            guard let values = sample.itemsAsNumericArray, !values.isEmpty else {
                SSLog.statisticsError("At least one sample is empty.")
                throw SSError(type: .missingData, file: #fileID, line: #line, function: #function, reason: "At least one sample is empty.")
            }
            var count: T = T.zero
            var meanAbsoluteDeviation: T = T.zero
            var m2: T = T.zero
            for value in values {
                let deviation = (value - center).magnitude
                count += T.one
                let delta = deviation - meanAbsoluteDeviation
                meanAbsoluteDeviation += delta / count
                let delta2 = deviation - meanAbsoluteDeviation
                m2 += delta * delta2
            }
            zMeans.append(meanAbsoluteDeviation)
            groupCounts.append(count)
            withinGroupSquares.append(m2)
            totalCount += count
            totalAbsoluteDeviation += meanAbsoluteDeviation * count
        }
        let k = T(samples.nCols)
        let N = totalCount
        let zMean = totalAbsoluteDeviation / N
        var betweenGroup: T = T.zero
        var withinGroup: T = T.zero
        var idx = 0
        while idx < zMeans.count {
            let count = groupCounts[idx]
            let diff = zMeans[idx] - zMean
            betweenGroup += count * diff * diff
            withinGroup += withinGroupSquares[idx]
            idx += 1
        }
        let numerator = (N - k) * betweenGroup
        let denominator = (k - T.one) * withinGroup
        let testStatisticValue = numerator / denominator
        let dist = try SwiftyBoost.Distribution.FisherF<T>(degreesOfFreedom1: k - T.one, degreesOfFreedom2: N - k)
        let cdfFRatio = try dist.cdf(testStatisticValue)
        let cutoffAlpha = try dist.quantile(1 - alpha)
        let cutoffAlpha90Percent = try dist.quantile(0.9)
        let cutoffAlpha95Percent = try dist.quantile(0.95)
        let cutoffAlpha99Percent = try dist.quantile(0.99)
        var result: VarianceEqualityTestResult<T> = VarianceEqualityTestResult<T>()
        result.df = VarianceEqualityTestResult.DegreesOfFreedom(df1: k - T.one, df2: N - k)
        result.cv90Pct = cutoffAlpha90Percent
        result.cv95Pct = cutoffAlpha95Percent
        result.cv99Pct = cutoffAlpha99Percent
        result.cvAlpha = cutoffAlpha
        result.pValue = 1 - cdfFRatio
        result.testStatistic = testStatisticValue
        result.equality = !(testStatisticValue > cutoffAlpha)
        switch testType {
        case .mean:
            result.testType = .levene(.mean)
        case .median:
            result.testType = .levene(.median)
        case .trimmedMean:
            result.testType = .levene(.trimmedMean)
        }
        return result
    }
    
    /// Convenience overload for the Levene variance-equality test.
    ///
    /// Wraps the provided groups into a `DataFrame` and forwards to
    /// ``leveneTest(samples:alpha:testType:)``.
    ///
    /// - Parameters:
    ///   - samples: Array of grouped observations.
    ///   - alpha: Significance level in `(0, 1)`.
    ///   - testType: Centering choice (`mean`, `median`, or 10% `trimmedMean`).
    /// - Returns: `VarianceEqualityTestResult` with critical values and p-value.
    /// - Throws: `SSError` when required statistics cannot be computed.
    public static func leveneTest(samples: [SSExamine<T, T>], alpha: T, testType: LeveneTestType) throws -> VarianceEqualityTestResult<T>? {
        precondition(!samples.isEmpty, "Levene test - no data")
        let dataFrame = try DataFrame(data: samples)
        return try leveneTest(samples: dataFrame, alpha: alpha, testType: testType)
    }
    
    /// One-sample chi-squared test for a population variance.
    ///
    /// With sample variance `s2` from `n` observations and a nominal variance
    /// `s0`, the test statistic
    /// ````text
    ///   X2 = (n - 1) * s2 / s0
    /// ````
    /// follows `chiSquare(df: n - 1)` under the null hypothesis
    /// `sigma^2 = s0`. Upper, lower, and two-tailed critical values are
    /// provided.
    ///
    /// - Parameters:
    ///   - sample: Sample container (`SSExamine`) providing `sampleVariance`.
    ///   - s0: Nominal variance to test against (must be positive).
    ///   - alpha: Significance level in `(0, 1)`.
    ///   - meanIsKnown: Flags whether the population mean is known.
    /// - Returns: `ChiSquaredVarianceTestResult` with ratios, cutoffs, and p-values.
    /// - Throws: `SSError` when variance is unavailable or inputs are invalid.
    public static func chiSquaredTest(sample: SSExamine<T, T>, nominalVariance s0: T, alpha: T, meanIsKnown: Bool = false) throws -> ChiSquaredVarianceTestResult<T>? {
        precondition(!sample.isEmpty, "No data")
        guard s0 > 0 else {
            SSLog.statisticsError("Chi-squared test - nominal variance must be positive")
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Chi-squared test - nominal variance must be positive")
        }
        let df: T = meanIsKnown ? T(sample.sampleSize) : T(sample.sampleSize) - T.one
        guard let _ = sample.sampleVariance, let _ = sample.populationVariance else {
            SSLog.statisticsError("Sample Variance is not available")
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Sample Variance is not available")
        }
        let v: T = meanIsKnown ? sample.populationVariance! : sample.sampleVariance!
        let ratio = v / s0
        let testStatistic: T = ratio * df
        let dist: Distribution.ChiSquared<T> = try .init(degreesOfFreedom: df)
        let cdfChiSquared = try dist.cdf(testStatistic)
        let lowerTail = cdfChiSquared
        let upperTail = T.one - cdfChiSquared
        let twoTailed: T = cdfChiSquared > T.half ? (T.one - cdfChiSquared) * T.two : cdfChiSquared * T.two
        let oneTailed: T = lowerTail < upperTail ? lowerTail : upperTail
        let cv = ChiSquaredVarianceTestResult<T>.ChiSquaredVarianceTestCritcalValues(oneTailedLowerCV: try dist.quantile(alpha), oneTailedUpperCV: try dist.quantile(1 - alpha), twoTailedLowerCV: try dist.quantile(alpha / T.two), twoTailedUpperCV: try dist.quantile(T.one - alpha / T.two))
        return ChiSquaredVarianceTestResult(df: df,
                                            ratio: ratio,
                                            testStatisticValue: testStatistic,
                                            p1Value: oneTailed,
                                            p2Value: twoTailed,
                                            sampleSize: T(sample.sampleSize),
                                            sigmaUEQs0: twoTailed >= alpha,
                                            sigmaLTEs0: lowerTail < alpha,
                                            sigmaGTEs0: upperTail < alpha,
                                            sd: T.sqrt(v),
                                            criticalValues: cv,
                                            meanAssumedKnown: meanIsKnown)
    }
    
    /// Convenience overload accepting a raw numeric array.
    ///
    /// The array is wrapped into an `SSExamine` container and forwarded to
    /// ``chiSquaredTest(sample:nominalVariance:alpha:meanIsKnown:)``.
    ///
    /// - Parameters:
    ///   - sample: Raw observations.
    ///   - s0: Nominal variance to test against (must be positive).
    ///   - alpha: Significance level in `(0, 1)`.
    ///   - meanIsKnown: Indicates, that the population mean is known
    /// - Returns: `ChiSquaredVarianceTestResult` with ratios, cutoffs, and p-values.
    /// - Throws: `SSError` when variance is unavailable or inputs are invalid.
    public static func chiSquaredTest(sample: Array<T>, nominalVariance s0: T, alpha: T, meanIsKnown: Bool = false) throws -> ChiSquaredVarianceTestResult<T>? {
        precondition(!sample.isEmpty, "No data")
        let e = SSExamine<T, T>.init(usingArray: sample, name: nil, characterSet: nil)
        return try chiSquaredTest(sample: e, nominalVariance: s0, alpha: alpha, meanIsKnown: meanIsKnown)
    }
    
    /// F-ratio test for comparing two variances.
    ///
    /// Tests whether the ratio of sample variances differs from 1 using the
    /// `F(df1 = n1 - 1, df2 = n2 - 1)` reference distribution. Returns two-sided
    /// and one-sided p-values plus confidence intervals for the ratio.
    ///
    /// The returned `p2Value` is the two-sided tail probability `2 * min(cdf, 1 - cdf)`,
    /// while `p1Value` is the smaller of the two one-sided tails. Directional flags
    /// (`FRatioLTE1`, `FRatioGTE1`) reflect those one-sided tails, and `FRatioEQ1`
    /// is `true` when the two-sided test does not reject at `alpha`.
    ///
    /// - Parameters:
    ///   - batch1: First sample container.
    ///   - batch2: Second sample container.
    ///   - alpha: Significance level in `(0, 1)`.
    /// - Returns: ``FRatioTestResult`` with variance estimates, p-values, and confidence intervals; `nil` if variances are undefined or zero.
    /// - Throws: `SSError` when inputs are invalid or statistics cannot be computed.
    public static func fRatioTest(batch1: SSExamine<T, T>, batch2: SSExamine<T, T>, alpha: T) throws -> FRatioTestResult<T>? {
        precondition(!batch1.isEmpty, "No data in first batch")
        precondition(!batch2.isEmpty, "No data in second batch")
        guard batch1.sampleSize >= 2 && batch2.sampleSize >= 2 else {
            SSLog.statisticsError("Batch samples must contain at least 2 observations each.")
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Batch samples must contain at least 2 observations each.")
        }
        if let s1 = batch1.sampleVariance, let s2 = batch2.sampleVariance, let m1 = batch1.arithmeticMean, let m2 = batch2.arithmeticMean {
            if s2.isZero {
                return nil
            }
            let testStat: T = s1 / s2
            let df1: T = T(batch1.sampleSize - 1)
            let df2: T = T(batch2.sampleSize - 1)
            let dist: Distribution.FisherF<T> = try .init(degreesOfFreedom1: df1, degreesOfFreedom2: df2)
            let cdf: T = try dist.cdf(testStat)
            let pGreater: T = T.one - cdf
            let pLess: T = cdf
            let pVariancesEqual: T = cdf > T.half ? T.two * (T.one - cdf) : T.two * cdf
            let pOneSided: T = pLess < pGreater ? pLess : pGreater
            let var1EQvar2: Bool = pVariancesEqual >= alpha
            let var1GTvar2: Bool = pGreater < alpha
            let var1LTvar2: Bool = pLess < alpha
            let ciEQ: ConfidenceInterval<T> = .init(
                lowerBound: try (testStat / dist.quantile(T.one - alpha / T.two)),
                upperBound: try (testStat / dist.quantile(alpha / T.two)) )
            let ciLT: ConfidenceInterval<T> = .init(
                lowerBound: .zero,
                upperBound: try (testStat / dist.quantile(alpha)) )
            let ciGT: ConfidenceInterval<T> = .init(
                lowerBound: try (testStat / dist.quantile(T.one - alpha)),
                upperBound: .infinity)
            let res: FRatioTestResult<T> = .init(sampleSize1: T(batch1.sampleSize),
                                                 sampleSize2: T(batch2.sampleSize),
                                                 dfDenominator: df2,
                                                 dfNumerator: df1,
                                                 variance1: s1,
                                                 variance2: s2,
                                                 FRatio: testStat,
                                                 p1Value: pOneSided,
                                                 p2Value: pVariancesEqual,
                                                 FRatioEQ1: var1EQvar2,
                                                 FRatioLTE1: var1LTvar2,
                                                 FRatioGTE1: var1GTvar2,
                                                 ciRatioEQ1: ciEQ,
                                                 ciRatioLTE1: ciLT,
                                                 ciRatioGTE1: ciGT,
                                                 sd1: T.sqrt(s1),
                                                 sd2: T.sqrt(s2),
                                                 mean1: m1,
                                                 mean2: m2,
                                                 postHocPower: (try powerFRatioTest(sampleSize1: batch1.sampleSize, sampleSize2: batch2.sampleSize, alpha: alpha, theta: testStat)))
            return res
        }
        else {
            return nil
        }
    }
    
    /// Convenience overload for the F-ratio variance comparison using raw arrays.
    ///
    /// Wraps the inputs into ``SSExamine`` containers and forwards to
    /// ``fRatioTest(batch1:batch2:alpha:)``.
    ///
    /// - Parameters:
    ///   - batch1: First sample.
    ///   - batch2: Second sample.
    ///   - alpha: Significance level in `(0, 1)`.
    /// - Returns: ``FRatioTestResult`` when the test can be evaluated; otherwise `nil`.
    /// - Throws: `SSError` when inputs are invalid or statistics cannot be computed.
    public static func fRatioTest(batch1: [T], batch2: [T], alpha: T) throws -> FRatioTestResult<T>? {
        let batch1Ex: SSExamine<T, T> = .init(usingArray: batch1, name: "batch1", characterSet: nil)
        let batch2Ex: SSExamine<T, T> = .init(usingArray: batch2, name: "batch2", characterSet: nil)
        return try fRatioTest(batch1: batch1Ex, batch2: batch2Ex, alpha: alpha)
    }
    
    private static func powerFRatioTest(sampleSize1 n1: Int, sampleSize2 n2: Int, alpha: T, theta: T) throws -> T? {
        guard n1 >= 2, n2 >= 2 else { return nil }
        let df1 = n1 - 1
        let df2 = n2 - 1
        let theta_hat = theta
        let dist: Distribution.FisherF<T> = try .init(degreesOfFreedom1: T(df1), degreesOfFreedom2:T(df2))
        let Flow: T = try dist.quantile(alpha / .two)
        let FHigh: T = try dist.quantile(1 - alpha / .two)
        let cdfLow = try dist.cdf(Flow / theta_hat)
        let cdfHigh = try dist.cdf(FHigh / theta_hat)
        let power_hat:T = cdfLow + (T.one - cdfHigh)
        return power_hat
    }
    
    private static func powerFRatioTest(sampleSize n: Int, alpha: T, theta: T) throws -> T? {
        guard n >= 2 else { return nil }
        let df1 = n - 1
        let df2 = n - 1
        let theta_hat = theta
        let dist: Distribution.FisherF<T> = try .init(degreesOfFreedom1: T(df1), degreesOfFreedom2:T(df2))
        let Flow: T = try dist.quantile(alpha / .two)
        let FHigh: T = try dist.quantile(1 - alpha / .two)
        let cdfLow = try dist.cdf(Flow / theta_hat)
        let cdfHigh = try dist.cdf(FHigh / theta_hat)
        let power_hat:T = cdfLow + (T.one - cdfHigh)
        return power_hat
    }
    
    /// Finds the smallest balanced sample size `n` such that the two-sided F test
    /// attains at least `targetPower` for effect size `theta = s1/s2`.
    ///
    /// - Parameters:
    ///   - alpha: Two-sided significance level for the reference F test.
    ///   - targetPower: Desired power threshold in `(0, 1)`.
    ///   - theta: Anticipated variance ratio `s1/s2`.
    ///   - nMin: Lower bound (inclusive) for the balanced group size search.
    ///   - nMax: Upper bound (inclusive) for the balanced group size search.
    /// - Returns: `(n, power)` for the first `n` meeting the target, or `nil` if no `n` in the range suffices or if power cannot be computed.
    public static func nFRatioTestBalanced(alpha: T, targetPower: T, theta: T, nMin: Int = 5, nMax: Int = 5000) throws -> (n: Int, power: T)? {
        for n in nMin...nMax {
            if let p = try powerFRatioTest(sampleSize: n , alpha: alpha, theta: theta) {
                if p >= targetPower {
                    return (n, p)
                }
            }
            else {
                return nil
            }
        }
        return nil
    }
    
    /// Holds the results of the F test for equal variances
    public struct FRatioTestResult<FP: RealLike>: CustomStringConvertible, Codable {
        /// size of batch 1
        public var sampleSize1: FP?
        /// size of batch 2
        public var sampleSize2: FP?
        /// denominator degrees of freedom
        public var dfDenominator: FP?
        /// numerator degrees of freedom
        public var dfNumerator: FP?
        /// variance of sample 1
        public var variance1: FP?
        /// variance of sample 2
        public var variance2: FP?
        /// F ratio
        public var FRatio: FP?
        /// one sided p-value
        public var p1Value: FP?
        /// two sided p-value
        public var p2Value: FP?
        /// indicates if variances are equal
        public var FRatioEQ1: Bool?
        /// indicates if var1 <= var2
        public var FRatioLTE1: Bool?
        /// indicates if var1 >= var2
        public var FRatioGTE1: Bool?
        /// confidence interval for var1 == var2
        public var ciRatioEQ1: ConfidenceInterval<FP>?
        /// confidence interval for var1 < var2
        public var ciRatioLTE1: ConfidenceInterval<FP>?
        /// confidence interval for var1 > var2
        public var ciRatioGTE1: ConfidenceInterval<FP>?
        /// Sample standard deviation of the first series.
        public var sd1: FP?
        /// Sample standard deviation of the second series.
        public var sd2: FP?
        /// Sample mean of the first series.
        public var mean1: FP?
        /// Sample mean of the second series.
        public var mean2: FP?
        /// Post-hoc (observed) power estimate for the two-sided F test based on the observed variance ratio.
        public var postHocPower: FP?
        /// Returns a textual summary of the F-ratio test.
        public var description: String {
            get {
                var descr = String()
                if let dfden = self.dfDenominator, let dfnum = self.dfNumerator, let var1 = self.variance1, let var2 = self.variance2, let f = self.FRatio, let n1 = self.sampleSize1, let n2 = self.sampleSize2, let cieq = self.ciRatioEQ1, let clt = self.ciRatioLTE1, let cig = self.ciRatioGTE1, let p1 = self.p1Value, let p2 = self.p2Value, let m1 = self.mean1, let m2 = self.mean2, let sd1 = self.sd1, let sd2 = self.sd2, let P = postHocPower {
                    descr.append("F-RATIO-TEST\n")
                    descr.append("**************************\n")
                    descr.append("variance 1: \(niceNumber(var1))\n")
                    descr.append("variance 2: \(niceNumber(var2))\n")
                    descr.append("sd1: \(niceNumber(sd1))\n")
                    descr.append("sd2: \(niceNumber(sd2))\n")
                    descr.append("mean1: \(niceNumber(m1))\n")
                    descr.append("mean2: \(niceNumber(m2))\n")
                    descr.append("n1: \(niceNumber(n1))\n")
                    descr.append("n2: \(niceNumber(n2))\n")
                    descr.append("F ratio (= effect size): \(niceNumber(f))\n")
                    descr.append("denominator df: \(niceNumber(dfden))\n")
                    descr.append("numerator df: \(niceNumber(dfnum))\n")
                    descr.append("one sided p-value: \(niceNumber(p1))\n")
                    descr.append("two sided p-value: \(niceNumber(p2))\n")
                    descr.append("post hoc power: \(niceNumber(P))\n")
                    descr.append("(1 - alpha)-CI for var1 = var2 \(cieq)\n")
                    descr.append("(1 - alpha)-CI for var1 <= var1 \(clt)\n")
                    descr.append("(1 - alpha)-CI for var1 >= var2 \(cig)\n")
                    descr.append("**************************\n")
                }
                return descr
            }
        }
    }
    
    /// Holds the results of the chi-square variance test with optional critical values.
    public struct ChiSquaredVarianceTestResult<FP: RealLike & Codable>: CustomStringConvertible, Codable {
        /// Critical values for one- and two-tailed chi-square variance tests.
        public struct ChiSquaredVarianceTestCritcalValues : Codable {
            /// Lower critical value for a one-tailed test.
            public var oneTailedLowerCV: FP
            /// Upper critical value for a one-tailed test.
            public var oneTailedUpperCV: FP
            /// Lower critical value for a two-tailed test.
            public var twoTailedLowerCV: FP
            /// Upper critical value for a two-tailed test.
            public var twoTailedUpperCV: FP
        }
        /// degrees of freedom
        public var df: FP?
        /// variance ratio
        public var ratio: FP?
        /// chi^2
        public var testStatisticValue: FP?
        /// one sided p-value
        public var p1Value: FP?
        /// two sided p-value
        public var p2Value: FP?
        /// sample size
        public var sampleSize: FP?
        /// True if sample variance == s0
        public var sigmaUEQs0: Bool?
        /// True if sample variance <= s0
        public var sigmaLTEs0: Bool?
        /// True if sample variance >= s0
        public var sigmaGTEs0: Bool?
        /// sample standard deviation
        public var sd: FP?
        /// Critical values used for the reported decision.
        public var criticalValues: ChiSquaredVarianceTestCritcalValues?
        /// Indicates whether the true mean was assumed known in the test.
        public var meanAssumedKnown: Bool?
        /// Returns a textual description of the chi-square variance test result.
        public var description: String {
            get {
                var descr = String()
                if let fr = self.ratio, let f = self.testStatisticValue, let p1 = self.p1Value, let p2 = self.p2Value, let ssd = self.sd, let n = self.sampleSize, let sdf = self.df, let cv = self.criticalValues, let sl = sigmaLTEs0, let sg = sigmaGTEs0, let su = sigmaUEQs0, let mk = self.meanAssumedKnown {
                    descr.append("CHI-SQUARE-TEST FOR EQUALITY OF VARIANCES\n")
                    descr.append("*****************************************\n")
                    descr.append("sample size: \(niceNumber(n))\n")
                    descr.append("Assuming true mean is known: \(mk ? "Yes" : "No")\n")
                    descr.append("degrees of freedom: \(niceNumber(sdf))\n")
                    descr.append("sd: \(niceNumber(ssd))\n")
                    descr.append("F ratio: \(niceNumber(fr))\n")
                    descr.append("f: \(niceNumber(f))\n")
                    descr.append("s < s0: \(sl)\n")
                    descr.append("s > s0: \(sg)\n")
                    descr.append("s <> s0: \(su)\n")
                    descr.append("lower critical value (one-tailed): \(niceNumber(cv.oneTailedLowerCV))\n")
                    descr.append("upper critical value (one-tailed): \(niceNumber(cv.oneTailedUpperCV))\n")
                    descr.append("lower critical value (two-tailed): \(niceNumber(cv.twoTailedLowerCV))\n")
                    descr.append("upper critical value (two-tailed): \(niceNumber(cv.twoTailedUpperCV))\n")
                    descr.append("one sided p-value: \(niceNumber(p1))\n")
                    descr.append("two sided p-value: \(niceNumber(p2))\n")
                }
                return descr
            }
        }
    }
    
    /// Result of tests for equality of variances (Bartlett, Levene, Brown-Forsythe)
    public struct VarianceEqualityTestResult<FP: RealLike>: CustomStringConvertible, Codable, Equatable {
        /// Codable container for degrees of freedom (replaces non-Codable tuple).
        public struct DegreesOfFreedom: Codable, Equatable {
            /// Numerator degrees of freedom.
            public var df1: FP
            /// Denominator degrees of freedom.
            public var df2: FP
            /// Creates a degrees-of-freedom pair.
            public init(df1: FP, df2: FP) {
                self.df1 = df1
                self.df2 = df2
            }
        }
        
        /// degrees of freedom
        public var df: DegreesOfFreedom? = nil
        /// critical value for alpha = 0.1
        public var cv90Pct: FP? = nil
        /// critical value for alpha = 0.05
        public var cv95Pct: FP? = nil
        /// critical value for alpha = 0.01
        public var cv99Pct: FP? = nil
        /// critical value for alpha as set by the call
        public var cvAlpha: FP? = nil
        /// p-value (= (1 - alpha)-quantile of test statistic)
        public var pValue: FP? = nil
        /// test statistic value
        public var testStatistic: FP? = nil
        /// Set to true if pValue >= alpha
        public var equality: Bool? = nil
        /// Test type
        public var testType: VarianceTestType? = nil
        /// Returns a description
        public var description: String {
            get {
                var descr = String()
                descr.append("TEST FOR EQUALITY OF VARIANCES\n")
                descr.append("******************************\n")
                if let tt = self.testType, let adf = self.df, let cv90 = self.cv90Pct, let cv95 = self.cv95Pct, let cv99 = self.cv99Pct, let eq = self.equality, let p = self.pValue, let ts = self.testStatistic {
                    descr.append("Test: \(tt)\n")
                    descr.append("variances are equal: \(eq)\n")
                    descr.append("p: \(niceNumber(p))\n")
                    descr.append("critical value for alpha = 0.01: \(niceNumber(cv99))\n")
                    descr.append("critical value for alpha = 0.05: \(niceNumber(cv95))\n")
                    descr.append("critical value for alpha = 0.10: \(niceNumber(cv90))\n")
                    descr.append("test statistic*: \(niceNumber(ts))\n")
                    descr.append("degrees of freddom (\(niceNumber(adf.df1)), \(niceNumber(adf.df2)))\n")
                }
                return descr
            }
        }
    }
    
    
    /// Defines the type of variance test
    public enum VarianceTestType: CustomStringConvertible, Codable, Equatable {
        /// Bartlett test
        case bartlett
        /// Variant of the Levene-Test
        case levene(LeveneTestType)
        
        /// Human-readable label describing the chosen variance test.
        public var description: String  {
            switch self {
            case .levene(.median):
                return "Levene/Brown-Forsythe (median)"
            case .levene(.mean):
                return "Levene (mean)"
            case .levene(.trimmedMean):
                return "Levene (trimmedMean)"
            case .bartlett:
                return "Bartlett-Test"
            }
        }
        
        private enum CodingKeys: CodingKey {
            case bartlett
            case levene
        }
        
        /// Encodes the variance test type for Codable conformance.
        public func encode(to encoder: Encoder) throws {
            var container = encoder.container(keyedBy: CodingKeys.self)
            switch self {
            case .bartlett:
                try container.encode("Bartlett", forKey: .bartlett)
            case .levene(let value):
                try container.encode(value, forKey: .levene)
            }
        }
        
        /// Decodes the variance test type for Codable conformance.
        public init(from decoder: Decoder) throws {
            let container = try decoder.container(keyedBy: CodingKeys.self)
            if let _ = try container.decodeIfPresent(String.self, forKey: .bartlett) {
                self = .bartlett
            }
            else {
                let value = try container.decode(LeveneTestType.self, forKey: .levene)
                self = .levene(value)
            }
        }
    }
    
    /// Defines the type if the Levene test.
    public enum LeveneTestType: String, Codable {
        /// Levene test using the median (Brown-Forsythe)
        case median
        /// Levene test using the mean
        case mean
        /// Levene test using the trimmed mean
        case trimmedMean
    }
}
