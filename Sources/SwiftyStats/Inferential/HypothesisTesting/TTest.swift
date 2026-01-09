//
//  Created by VT on 14.12.25.
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

import SwiftyStatsPrelude

extension Inferential.HypothesisTesting.ParametricTests {
    /// Performs a two-sample Student t test (equal-variance and Welch variants).
    ///
    /// This routine computes the classic pooled-variance t statistic assuming equal variances
    /// and the Welch t statistic allowing unequal variances. Both one- and two-tailed p-values
    /// are reported for each variant, and critical values are computed at the requested
    /// significance level `alpha` using the corresponding Student-t distributions.
    ///
    /// - Mathematical definition:
    ///   - Let `x̄₁, s₁², n₁` and `x̄₂, s₂², n₂` be the sample means, variances, and sizes.
    ///   - Equal-variance (pooled) t statistic:
    ///     `t_EQ = (x̄₁ - x̄₂) / (s_p * sqrt(1/n₁ + 1/n₂))`, where
    ///     `s_p² = ((n₁-1)s₁² + (n₂-1)s₂²) / (n₁ + n₂ - 2)` and `df_EQ = n₁ + n₂ - 2`.
    ///   - Unequal-variance (Welch) t statistic:
    ///     `t_W = (x̄₁ - x̄₂) / sqrt(s₁²/n₁ + s₂²/n₂)` with
    ///     `df_W = (s₁²/n₁ + s₂²/n₂)² / ( (s₁²/n₁)²/(n₁-1) + (s₂²/n₂)²/(n₂-1) )`.
    ///
    /// - Assumptions:
    ///   - Independent samples drawn from populations with finite variance and mean.
    ///   - For the pooled test: population variances are equal. For Welch: variances may differ.
    ///   - Approximate normality of the sampling distribution of the mean (by CLT for moderate n).
    ///
    /// - Hypotheses (two-sided): `H₀: μ₁ = μ₂` vs `H₁: μ₁ ≠ μ₂`.
    ///   One-sided decisions use the same t statistics with tail selection matching the sign of `x̄₁ - x̄₂`.
    ///
    /// - Decision rule:
    ///   - Two-sided: reject `H₀` if `|t| > t_{1-α/2, df}` (equivalently `p < α`).
    ///   - One-sided: reject `H₀` if `t > t_{1-α, df}` (upper) or `t < t_{α, df}` (lower).
    ///
    /// - Variance equality check:
    ///   A Levene median test is computed to aid in choosing between pooled and Welch variants.
    ///
    /// - Parameters:
    ///   - sample1: First sample; must have `n₁ ≥ 2` with finite mean and variance.
    ///   - sample2: Second sample; must have `n₂ ≥ 2` with finite mean and variance.
    ///   - alpha: Significance level in `(0, 1)` used for critical values and two-tailed bounds.
    /// - Returns: A `TwoSampleTTestResult` populated with test statistics, degrees of freedom,
    ///   p-values (one- and two-tailed), critical values, effect sizes `r = |t| / sqrt(t² + df)`,
    ///   and Levene test outcome; returns `nil` if inputs are invalid.
    ///
    /// - Complexity: `O(n₁ + n₂)` to compute sample summaries and statistics.
    ///
    /// - SeeAlso: ``oneSampleTTest(sample:mean:alpha:)``, ``oneWayANOVA(data:alpha:)``
    ///
    /// - Example:
    ///   ```swift
    ///   // Compare two groups at alpha = 0.05
    ///   if let res = try Parametric<Double>.twoSampleTTest(sample1: examA, sample2: examB, alpha: 0.05) {
    ///       // Two-sided decision using Welch by default if variances differ
    ///       let p = res.p2Welch ?? res.p2UEQVAR
    ///       let cv = res.CVUEQVAR
    ///       let t  = res.tWelch ?? res.tUEQVAR
    ///       // Reject H0 if |t| > cv  <=>  p < alpha
    ///       if let p, let cv, let t { print(abs(t) > cv ? "Reject H0" : "Fail to reject H0", "(p=", p, ")") }
    ///   }
    ///   ```
    ///
    /// - Throws: `SSError` if Student-t distribution CDF/quantile evaluation fails.
    public static func twoSampleTTest(sample1: SSExamine<T, T>, sample2: SSExamine<T, T>, alpha: T) throws -> TwoSampleTTestResult<T>? {
        guard sample1.sampleSize > 1, sample2.sampleSize > 1 else {
            return nil
        }
        var cdfTValueEqualVariances: T = 0
        var cdfTValueUnequalVariances: T = 0
        var dfEqualVariances: T = 0
        var dfUnequalVariances: T = 0
        var differenceInMeans: T = 0
        //        var mean1: T = 0
        //        var mean2: T = 0
        var criticalValueEqualVariances: T = 0
        var criticalValueUnequalVariances: T = 0
        var pooledStdDev: T = 0
        var pooledVariance: T = 0
        var stdDev1: T = 0
        var stdDev2: T = 0
        var tValueEqualVariances: T = 0
        var tValueUnequalVariances: T = 0
        var variancesAreEqualMedian: Bool = false
        var twoTailedEV: T = 0
        var twoTailedUEV: T = 0
        var oneTailedEV: T = 0
        var oneTailedUEV: T = 0
        //        var var1: T = 0
        //        var var2: T = 0
        var t1: T, t2: T, t3: T
        let n1 = T(sample1.sampleSize)
        let n2 = T(sample2.sampleSize)
        var dist: Distribution.StudentT<T>
        if let var1 = sample1.sampleVariance,
           let var2 = sample2.sampleVariance,
           let mean1 = sample1.arithmeticMean,
           let mean2 = sample2.arithmeticMean
        {
            t1 = var1 / n1
            t2 = var2 / n2
            t3 = t1 + t2
            let k: T = t1 / t3
            t1 = T.pow(k, 2) / (n1 - .one)
            let ktk = k * k
            let kpk = k + k
            t2 = .one - kpk + ktk
            t3 = n2 - .one
            dfUnequalVariances = ceil(.one / (t1 + (t2 / t3)))
            dfEqualVariances = n1 + n2 - .two
            stdDev1 = T.sqrt(var1)
            stdDev2 = T.sqrt(var2)
            pooledVariance = (((n1 - .one) * var1) + ((n2 - .one) * var2)) / dfEqualVariances
            pooledStdDev = T.sqrt(pooledVariance)
            differenceInMeans = mean1 - mean2
            tValueEqualVariances = differenceInMeans / (pooledStdDev * T.sqrt(T.one / n1 + T.one / n2))
            tValueUnequalVariances = differenceInMeans / T.sqrt(var1 / n1 + var2 / n2)
            dist = try .init(degreesOfFreedom: dfEqualVariances)
            cdfTValueEqualVariances = try dist.cdf(tValueEqualVariances)
            criticalValueEqualVariances = try dist.quantile(.one - alpha / .two)
            dist = try .init(degreesOfFreedom: dfUnequalVariances)
            cdfTValueUnequalVariances = try dist.cdf(tValueUnequalVariances)
            criticalValueUnequalVariances = try dist.quantile(.one - alpha / .two)
            if let leveneResult: Variance<T>.VarianceEqualityTestResult = try Variance<T>.leveneTest(samples: [sample1, sample2], alpha: alpha, testType: .median) {
                let cdfLeveneMedian = leveneResult.pValue ?? .zero
                variancesAreEqualMedian = leveneResult.equality ?? false
                if cdfTValueEqualVariances > .half {
                    twoTailedEV = (.one - cdfTValueEqualVariances) * .two
                    oneTailedEV = .one - cdfTValueEqualVariances
                }
                else {
                    twoTailedEV = cdfTValueEqualVariances * .two
                    oneTailedEV = cdfTValueEqualVariances
                }
                if cdfTValueUnequalVariances > .half {
                    twoTailedUEV = (.one - cdfTValueUnequalVariances) * .two
                    oneTailedUEV = .one - cdfTValueUnequalVariances
                }
                else {
                    twoTailedUEV = cdfTValueUnequalVariances * .two
                    oneTailedUEV = cdfTValueUnequalVariances
                }
                let effectSize_EV: T = T.sqrt((tValueEqualVariances * tValueEqualVariances) / ((tValueEqualVariances * tValueEqualVariances) + dfEqualVariances))
                let effectSize_UEV: T = T.sqrt((tValueUnequalVariances * tValueUnequalVariances) / ((tValueUnequalVariances * tValueUnequalVariances) + dfUnequalVariances))
                // Welch
                let var1OverN1 = var1 / n1
                let var2OverN2 = var2 / n2
                let sumVar = var1OverN1 + var2OverN2
                t1 = (n1 * n1 * (n1 - .one))
                t2 = T.pow(var1, 2)
                t3 = (n2 * n2 * (n2 - .one))
                let denomWelchDF: T = T.pow(var1, 2) / t1 + T.pow(var2, 2) / t3
                let welchT: T = (mean1 - mean2) / T.sqrt(var1OverN1 + var2OverN2)
                let welchDF: T = (sumVar * sumVar) / denomWelchDF
                dist = try .init(degreesOfFreedom: welchDF)
                let cdfWelch = try dist.cdf(welchT)
                var twoSidedWelch: T, oneSidedWelch: T
                if cdfWelch > .half {
                    twoSidedWelch = (.one - cdfWelch) * .two
                    oneSidedWelch = .one - cdfWelch
                }
                else {
                    twoSidedWelch = cdfWelch * .two
                    oneSidedWelch = cdfWelch
                }
                var result: TwoSampleTTestResult<T> = .init()
                result.p2EQVAR = twoTailedEV
                result.p2UEQVAR = twoTailedUEV
                result.p1EQVAR = oneTailedEV
                result.p1UEQVAR = oneTailedUEV
                result.mean1 = mean1
                result.mean2 = mean2
                result.sampleSize1 = n1
                result.sampleSize2 = n2
                result.stdDev1 = stdDev1
                result.stdDev2 = stdDev2
                result.pooledStdDev = pooledStdDev
                result.pooledVariance = pooledVariance
                result.differenceInMeans = differenceInMeans
                result.tEQVAR = tValueEqualVariances
                result.tUEQVAR = tValueUnequalVariances
                result.LeveneP = cdfLeveneMedian
                result.dfEQVAR = dfEqualVariances
                result.dfUEQVAR = dfUnequalVariances
                let lowerTwoTailed:T = alpha / .two
                let upperTwoTailed: T = .one - lowerTwoTailed
                result.mean1GTEmean2 = variancesAreEqualMedian ? ((cdfTValueEqualVariances > upperTwoTailed) ? true : false) : ((cdfTValueUnequalVariances > upperTwoTailed) ?  true : false)
                result.mean1LTEmean2 = (variancesAreEqualMedian) ? ((cdfTValueEqualVariances < lowerTwoTailed) ? true : false) : ((cdfTValueUnequalVariances < lowerTwoTailed) ? true : false)
                result.mean1EQmean2 = (variancesAreEqualMedian) ? ((cdfTValueEqualVariances >= lowerTwoTailed && cdfTValueEqualVariances <= upperTwoTailed) ? true : false) : ((cdfTValueUnequalVariances >= lowerTwoTailed && cdfTValueUnequalVariances <= upperTwoTailed) ? true : false)
                result.mean1UEQmean2 = (variancesAreEqualMedian) ? ((cdfTValueEqualVariances < lowerTwoTailed || cdfTValueEqualVariances > upperTwoTailed) ? true : false ) : ((cdfTValueUnequalVariances < lowerTwoTailed || cdfTValueUnequalVariances > upperTwoTailed) ? true : false)
                result.CVEQVAR = criticalValueEqualVariances
                result.CVUEQVAR = criticalValueUnequalVariances
                result.rUEQVAR = effectSize_UEV
                result.rEQVAR = effectSize_EV
                result.tWelch = welchT
                result.dfWelch = welchDF
                result.p2Welch = twoSidedWelch
                result.p1Welch = oneSidedWelch
                result.variancesAreEqual = variancesAreEqualMedian
                return result
            }
            else {
                return nil
            }
        }
        else {
            return nil
        }
    }
    
    /// Convenience overload creating `SSExamine` from arrays and delegating to `twoSampleTTest`.
    /// - Parameters:
    ///   - sample1: First sample values.
    ///   - sample2: Second sample values.
    ///   - alpha: Significance level.
    /// - Returns: See `TwoSampleTTestResult`.
    ///
    /// - Example:
    ///   ```swift
    ///   let a = [2.1, 2.3, 1.9, 2.5]
    ///   let b = [1.7, 1.6, 1.9]
    ///   let res = try Parametric<Double>.twoSampleTTest(sample1: a, sample2: b, alpha: 0.05)
    ///   ```
    public static func twoSampleTTest(sample1: [T], sample2: [T], alpha: T) throws -> TwoSampleTTestResult<T>? {
        let ex1: SSExamine<T, T> = .init(usingArray: sample1, name: nil, characterSet: nil)
        let ex2: SSExamine<T, T> = .init(usingArray: sample2, name: nil, characterSet: nil)
        return try twoSampleTTest(sample1: ex1, sample2: ex2, alpha: alpha)
    }
    
    /// Performs a one-sample Student t test against a hypothesised mean `μ₀`.
    ///
    /// - Mathematical definition:
    ///   - Let `x̄` be the sample mean, `s` the sample standard deviation, `n` the sample size.
    ///   - Test statistic: `t = (x̄ - μ₀) / (s / √n)` with `df = n - 1`.
    ///
    /// - Assumptions:
    ///   - Independent observations with finite mean and variance; approximate normality of the
    ///     sampling distribution of the mean (by CLT for moderate n).
    ///
    /// - Hypotheses:
    ///   - Two-sided: `H₀: μ = μ₀` vs `H₁: μ ≠ μ₀`.
    ///   - One-sided: `H₁: μ > μ₀` (upper) or `H₁: μ < μ₀` (lower); tail chosen by sign of `x̄ - μ₀`.
    ///
    /// - Decision rule:
    ///   - Two-sided: reject `H₀` if `|t| > t_{1-α/2, df}` (equivalently `p < α`).
    ///   - One-sided: reject `H₀` if `t > t_{1-α, df}` or `t < t_{α, df}`.
    ///
    /// - Parameters:
    ///   - sample: Sample to test; `n ≥ 2` with finite mean and variance.
    ///   - mean: Hypothesised population mean `μ₀`.
    ///   - alpha: Significance level in `(0, 1)`.
    ///
    /// - Complexity: `O(n)` where `n` is the sample size.
    ///
    /// - SeeAlso: ``twoSampleTTest(sample1:sample2:alpha:)-([T],_,_)``, ``oneWayANOVA(data:alpha:)``
    ///
    /// - Example:
    ///   ```swift
    ///   // Test whether the mean differs from 0 at alpha = 0.05
    ///   if let res = try Parametric<Double>.oneSampleTTest(sample: exam, mean: 0.0, alpha: 0.05) {
    ///       let p2 = res.p2Value!
    ///       let t  = res.tStat!
    ///       let cv = res.cv95Pct!  // critical value for alpha = 0.05 (two-sided)
    ///       print(abs(t) > cv ? "Reject H0" : "Fail to reject H0", "(p=", p2, ")")
    ///   }
    ///   ```
    ///
    /// - Throws: `SSError` if Student-t distribution CDF/quantile evaluation fails.
    public static func oneSampleTTest(sample: SSExamine<T, T>, mean: T, alpha: T) throws -> OneSampleTTestResult<T>? {
        guard sample.sampleSize > 1 else {
            return nil
        }
        guard let sd = sample.sampleStandardDeviation, let sampleMean = sample.arithmeticMean else {
            return nil
        }
        var tValue: T = 0
        var pValue: T = 0
        let N: T = T(sample.sampleSize)
        let diffMean: T = sampleMean - mean
        let twoTailed: T, oneTailed: T
        let dist: Distribution.StudentT<T> = try .init(degreesOfFreedom: N - 1)
        tValue = diffMean / (sd / T.sqrt(N))
        pValue = try dist.cdf(tValue)
        if pValue > .half {
            twoTailed = (.one - pValue) * .two
            oneTailed = .one - pValue
        }
        else {
            twoTailed = pValue * .two
            oneTailed = pValue
        }
        var result: OneSampleTTestResult<T> = .init()
        result.p1Value = oneTailed
        result.p2Value = twoTailed
        result.tStat = tValue
        result.cv90Pct = try dist.quantile(.one - T(0.05))
        result.cv95Pct = try dist.quantile(.one - T(0.025))
        result.cv99Pct = try dist.quantile(.one - T(0.005))
        result.cvAlpha = try dist.quantile(.one - (alpha / T.two))
        result.mean0 = mean
        result.mean = sampleMean
        result.sampleSize = N
        result.difference = diffMean
        result.stdDev = sd
        result.stdErr = sd / T.sqrt(N)
        result.df = N - 1
        result.meanEQtestValue = ((pValue < (alpha / T.two)) || (pValue > (T.one - (alpha / T.two)))) ? false : true
        result.meanLTEtestValue = (pValue < alpha) ? true : false
        result.meanGTEtestValue = (pValue > (1 - alpha)) ? true : false
        return result
    }
    
    /// Convenience overload creating `SSExamine` from an array and delegating to `oneSampleTTest`.
    /// - Parameters:
    ///   - sample: Sample values.
    ///   - mean: Hypothesised mean `μ₀`.
    ///   - alpha: Significance level.
    /// - Returns: See `OneSampleTTestResult`.
    ///
    /// - Example:
    ///   ```swift
    ///   let x = [0.1, -0.2, 0.3, 0.0]
    ///   let res = try Parametric<Double>.oneSampleTTest(sample: x, mean: 0.0, alpha: 0.05)
    ///   ```
    public static func oneSampleTTest(sample: [T], mean: T, alpha: T) throws -> OneSampleTTestResult<T>? {
        let ex1: SSExamine<T,T> = .init(usingArray: sample, name: nil, characterSet: nil)
        return try oneSampleTTest(sample: ex1, mean: mean, alpha: alpha)
    }

    /// Performs one-way ANOVA (fixed effects) to compare `k ≥ 2` group means.
    ///
    /// Use this test to determine whether there is evidence that at least one group mean
    /// differs from the others. This implementation also reports Bartlett's and Levene's
    /// tests for homoscedasticity to help assess the equal-variance assumption, and it
    /// optionally computes Tukey–Kramer post hoc comparisons when the omnibus F test is
    /// significant.
    ///
    /// - Model:
    ///   Observations follow `Y_{ij} = μ + τ_i + ε_{ij}`, where `τ_i` are fixed group effects and
    ///   `ε_{ij}` are i.i.d. with zero mean and constant variance `σ²`.
    ///
    /// - Hypotheses:
    ///   - Null: `H₀: μ₁ = μ₂ = … = μ_k` (all group means equal)
    ///   - Alternative: `H₁:` at least one mean differs
    ///
    /// - Test statistic:
    ///   - Treatment (between) sum of squares: `SS_T = Σ_i n_i (x̄_i − x̄_·)²`
    ///   - Error (within) sum of squares: `SS_E = Σ_i (n_i − 1) s_i²`
    ///   - Mean squares: `MS_T = SS_T / (k − 1)`, `MS_E = SS_E / (N − k)`
    ///   - F statistic: `F = MS_T / MS_E` with `df₁ = k − 1`, `df₂ = N − k`
    ///
    /// - Decision rule:
    ///   Reject `H₀` if `F > F_{1−α; df₁, df₂}` (equivalently if `p < α`).
    ///
    /// - Assumptions:
    ///   - Independent observations within and across groups
    ///   - Approximately normal distributions within groups (or sufficiently large samples)
    ///   - Homoscedasticity (equal variances across groups); Bartlett and Levene tests are
    ///     reported to aid assessment, but do not change the F test itself
    ///
    /// - Parameters:
    ///   - data: Array of group samples (`k ≥ 2`). Each `SSExamine` must provide finite mean and variance.
    ///   - alpha: Significance level in `(0, 1)` used to compute the F critical value.
    ///
    /// - Returns: A ``OneWayANOVAResult`` containing the F statistic, two-tailed p-value, F critical
    ///   value at `alpha`, sums and mean squares, degrees of freedom, and p-values of Bartlett and
    ///   Levene tests. If inputs are invalid (e.g., fewer than two groups, empty group, or missing
    ///   summaries), returns `nil`.
    ///
    /// - Post Hoc:
    ///   When possible, Tukey–Kramer pairwise comparisons are computed and returned in
    ///   ``OneWayANOVAResult/postHocTests``. These are appropriate when group sizes differ.
    ///
    /// - Complexity: `O(N)` where `N` is the total number of observations across all groups.
    ///
    /// - SeeAlso:
    ///   - ``twoSampleTTest(sample1:sample2:alpha:)-([T],_,_)``
    ///   - ``oneSampleTTest(sample:mean:alpha:)-([T],_,_)``
    ///
    /// - Example:
    ///   ```swift
    ///   // Compare k group means at alpha = 0.05
    ///   let groups: [SSExamine<Double, Double>] = [g1, g2, g3]
    ///   if let anova = try Parametric<Double>.oneWayANOVA(data: groups, alpha: 0.05) {
    ///       let F  = anova.FStatistic!
    ///       let cv = anova.cv!
    ///       print(F > cv ? "Reject H0 (means differ)" : "Fail to reject H0")
    ///       // Optional: inspect homoscedasticity tests and post hoc results
    ///       if let pB = anova.pBartlett, let pL = anova.pLevene { print("Bartlett:", pB, "Levene:", pL) }
    ///       if let postHoc = anova.postHocTests { for comp in postHoc { print(comp) } }
    ///   }
    ///   ```
    ///
    /// - Throws: `SSError` if F distribution CDF/quantile evaluation fails.
    public static func oneWayANOVA(data: Array<SSExamine<T, T>>, postHocTest pht: PostHoc<T>.PostHocTestType = .tukeyKramer, alpha: T) throws -> OneWayANOVAResult<T>? {
        guard data.count > 1 else { return nil }
        var overallMean: T = .zero
        var N: Int = 0
        var pBartlett: T = .zero
        var pLevene: T = .zero
        // explained sum of squares
        var SST: T = .zero
        // residual sum of squares
        var SSE: T = .zero
        var F: T = .zero
        var pValue: T = .zero
        var cdf: T = .zero
        var cutoffAlpha: T = .zero
        var meansEqual: Bool = false
        let groups: T = T(data.count)
        var sum: T = .zero
        if let bartlettTest: Variance<T>.VarianceEqualityTestResult = try Variance<T>.bartlettTest(samples: data , alpha: alpha),
           let leveneTest: Variance<T>.VarianceEqualityTestResult = try Variance<T>.leveneTest(samples: data, alpha: alpha, testType: .median) {
            pBartlett = bartlettTest.pValue ?? .nan
            pLevene = leveneTest.pValue ?? .nan
        }
        for ex in data {
            guard ex.sampleSize > 0 else { return nil }
            if let total = ex.total {
                sum += total
            }
            N += ex.sampleSize
        }
        overallMean = sum / T(N)
        var nn: T
        for ex in data {
            if let m = ex.arithmeticMean, let variance = ex.sampleVariance {
                nn = T(ex.sampleSize)
                SST += T.pow(m - overallMean, 2) * nn
                SSE += T(ex.sampleSize - 1) * variance
            }
            else {
                return nil
            }
        }
        let MSE: T = SSE / (T(N) - groups)
        let MST: T = SST / (groups - .one)
        F = MST / MSE
        let dist: Distribution.FisherF<T> = try .init(degreesOfFreedom1: groups - .one, degreesOfFreedom2: T(N) - groups)
        cdf = try dist.cdf(F)
        cutoffAlpha = try dist.quantile(.one - alpha)
        pValue = .one - cdf
        meansEqual = F <= cutoffAlpha
        let SSTotal: T = SST + SSE
        var result: OneWayANOVAResult<T> = .init()
        var postHoc: [PostHoc<T>.PostHocTestSummary]?
        switch pht {
        case .tukeyKramer:
            postHoc = try PostHoc<T>.tukeyKramerTest(data: data, alpha: alpha)
        case .scheffe: postHoc = try PostHoc<T>.scheffeTest(data: data, alpha: alpha)
        }
        result.postHocTestType = pht
        result.p2Value = pValue
        result.FStatistic = F
        result.alpha = alpha
        result.meansEQUAL = meansEqual
        result.cv = cutoffAlpha
        result.pBartlett = pBartlett
        result.pLevene = pLevene
        result.SSTotal = SSTotal
        result.SSError = SSE
        result.SSTreatment = SST
        result.MSError = MSE
        result.MSTreatment = MST
        result.dfError = T(N) - groups
        result.dfTreatment = groups - .one
        result.dfTotal = T(N) - .one
        result.postHocTests = postHoc
        return result
    }

    public static func oneWayANOVA(data: DataFrame<T,T>, postHocTest pht: PostHoc<T>.PostHocTestType = .tukeyKramer, alpha: T) throws -> OneWayANOVAResult<T>? {
        return try oneWayANOVA(data: data.data, postHocTest: pht, alpha: alpha)
    }

    /// Results for the two-sample t test (pooled and Welch variants).
    ///
    /// Includes t statistics, degrees of freedom, one- and two-tailed p-values, critical values,
    /// pooled variance/SD, and an effect size `r = |t| / sqrt(t² + df)` for both variants.
    public struct TwoSampleTTestResult<FP: RealLike>: CustomStringConvertible, Codable  {
        /// One-tailed p-value assuming equal variances (pooled t).
        public var p1EQVAR: FP?
        /// One-tailed p-value assuming unequal variances (unpooled t).
        public var p1UEQVAR: FP?
        /// Two-tailed p-value assuming equal variances (pooled t).
        public var p2EQVAR: FP?
        /// Two-tailed p-value assuming unequal variances (unpooled t).
        public var p2UEQVAR: FP?
        /// Sample mean of group 1 (x̄₁).
        public var mean1: FP?
        /// Sample mean of group 2 (x̄₂).
        public var mean2: FP?
        /// Sample size of group 1 (n₁).
        public var sampleSize1: FP?
        /// Sample size of group 2 (n₂).
        public var sampleSize2: FP?
        /// Sample standard deviation of group 1 (s₁).
        public var stdDev1: FP?
        /// Sample standard deviation of group 2 (s₂).
        public var stdDev2: FP?
        /// Pooled standard deviation s_p = sqrt( [ (n₁−1)s₁² + (n₂−1)s₂² ] / (n₁+n₂−2) ).
        public var pooledStdDev: FP?
        /// Pooled variance s_p² = [ (n₁−1)s₁² + (n₂−1)s₂² ] / (n₁+n₂−2).
        public var pooledVariance: FP?
        /// Difference in sample means (x̄₁ − x̄₂).
        public var differenceInMeans: FP?
        /// t statistic under equal-variance assumption: t = (x̄₁ − x̄₂) / ( s_p * sqrt(1/n₁ + 1/n₂) ).
        public var tEQVAR: FP?
        /// t statistic under unequal-variance assumption: t = (x̄₁ − x̄₂) / sqrt(s₁²/n₁ + s₂²/n₂).
        public var tUEQVAR: FP?
        /// p-value from Levene's median test for equal variances. Lower values suggest heteroscedasticity.
        public var LeveneP: FP?
        /// Degrees of freedom for equal-variance test (n₁ + n₂ − 2).
        public var dfEQVAR: FP?
        /// Degrees of freedom used by the unpooled variant (Satterthwaite approximation).
        public var dfUEQVAR: FP?
        /// Result of Levene test at the chosen α (true means fail to reject equal variances).
        public var variancesAreEqual: Bool? = false
        /// Decision flag: mean₁ ≥ mean₂ under the selected variant at α.
        public var mean1GTEmean2: Bool? = false
        /// Decision flag: mean₁ ≤ mean₂ under the selected variant at α.
        public var mean1LTEmean2: Bool? = false
        /// Decision flag: mean₁ = mean₂ (fail to reject) under the selected variant at α.
        public var mean1EQmean2: Bool? = false
        /// Decision flag: mean₁ ≠ mean₂ (reject) under the selected variant at α.
        public var mean1UEQmean2: Bool? = false
        /// Two-sided critical value t_{1−α/2, df} for equal-variance test.
        public var CVEQVAR: FP?
        /// Two-sided critical value t_{1−α/2, df} for unequal-variance test.
        public var CVUEQVAR: FP?
        /// Effect size r = |t| / sqrt(t² + df) for the equal-variance test.
        public var rEQVAR: FP?
        /// Effect size r = |t| / sqrt(t² + df) for the unequal-variance test.
        public var rUEQVAR: FP?
        /// Welch's t statistic using Satterthwaite degrees of freedom.
        public var tWelch: FP?
        /// Welch–Satterthwaite degrees of freedom.
        public var dfWelch: FP?
        /// Two-tailed p-value for Welch's test.
        public var p2Welch: FP?
        /// One-tailed p-value for Welch's test.
        public var p1Welch: FP?
        
        /// Human-readable multiline summary of results. Values are present only when all required fields are non-nil.
        public var description: String {
            get {
                var descr = String()
                if let m1GTEm2 = self.mean1GTEmean2, let m1LTEm2 = self.mean1LTEmean2, let m1EQm2 = self.mean1EQmean2, let m1UEqm2 = self.mean1UEQmean2, let p2EV = self.p2EQVAR, let p2UE = self.p2UEQVAR, let p2w = self.p2Welch, let p1ev = self.p1EQVAR, let p1ue = self.p1UEQVAR, let p1w = self.p1Welch, let m1 = self.mean1, let m2 = self.mean2, let n1 = self.sampleSize1, let n2 = self.sampleSize2, let sd1 = self.stdDev1, let sd2 = self.stdDev2, let psd = self.pooledStdDev, let diff = self.differenceInMeans, let teq = self.tEQVAR, let tueq = self.tUEQVAR, let cvue = self.CVUEQVAR, let cveq = self.CVEQVAR, let req = self.rEQVAR, let rueq = self.rUEQVAR, let tw = self.tWelch, let lev = self.LeveneP {
                    descr.append("TWO-SAMPLE T-TEST\n")
                    descr.append("*****************\n")
                    descr.append("mean group 1: \(niceNumber(m1))\n")
                    descr.append("mean group 2: \(niceNumber(m2))\n")
                    descr.append("sample size group 1: \(niceNumber(n1))\n")
                    descr.append("sample size group 2: \(niceNumber(n2))\n")
                    descr.append("p-value of Levene test: \(niceNumber(lev))\n")
                    descr.append("sd group 1: \(niceNumber(sd1))\n")
                    descr.append("sd group 2: \(niceNumber(sd2))\n")
                    descr.append("pooled sd: \(niceNumber(psd))\n")
                    descr.append("difference in means: \(niceNumber(diff))\n")
                    descr.append("t value for equal variances: \(niceNumber(teq))\n")
                    descr.append("t value for unequal variances: \(niceNumber(tueq))\n")
                    descr.append("t value Welch: \(niceNumber(tw))\n")
                    descr.append("critical value for equal variances: \(niceNumber(cveq))\n")
                    descr.append("critical value for unequal variances: \(niceNumber(cvue))\n")
                    descr.append("one sided p-value for equal variances: \(niceNumber(p1ev))\n")
                    descr.append("one sided p-value for unequal variances: \(niceNumber(p1ue))\n")
                    descr.append("one sided p-value (Welch): \(niceNumber(p1w))\n")
                    descr.append("two sided p-value for equal variances: \(niceNumber(p2EV))\n")
                    descr.append("two sided p-value for unequal variances: \(niceNumber(p2UE))\n")
                    descr.append("two sided p-value (Welch): \(niceNumber(p2w))\n")
                    descr.append("effect size for equal variances*: \(niceNumber(req))\n")
                    descr.append("effect size for unequal variances*: \(niceNumber(rueq))\n")
                    descr.append("mean1 >= mean2: \(m1GTEm2)\n")
                    descr.append("mean1 <= mean2: \(m1LTEm2)\n")
                    descr.append("mean1 == mean2: \(m1EQm2)\n")
                    descr.append("mean1 <> mean2: \(m1UEqm2)\n")
                }
                return descr
            }
        }
    }
    
    /// Results for the one-sample t test against a hypothesised mean.
    ///
    /// Provides the t statistic, degrees of freedom, p-values (one- and two-tailed), and critical
    /// values at common levels (0.10, 0.05, 0.01) and a user-supplied `alpha`.
    public struct OneSampleTTestResult<FP: RealLike>: CustomStringConvertible,Codable  {
        /// One-tailed p-value (tail chosen by sign of x̄ − μ₀).
        public var p1Value: FP?
        /// Two-tailed p-value.
        public var p2Value: FP?
        /// t statistic: t = (x̄ − μ₀) / ( s / √n ).
        public var tStat: FP?
        /// Critical value t_{0.95, df} (α = 0.10, two-sided).
        public var cv90Pct: FP?
        /// Critical value t_{0.975, df} (α = 0.05, two-sided).
        public var cv95Pct: FP?
        /// Critical value t_{0.995, df} (α = 0.01, two-sided).
        public var cv99Pct: FP?
        /// Critical value t_{1−α/2, df} for the supplied α (two-sided).
        public var cvAlpha: FP?
        /// Sample mean (x̄).
        public var mean: FP?
        /// Hypothesised mean (μ₀).
        public var mean0: FP?
        /// Difference x̄ − μ₀.
        public var difference: FP?
        /// Sample size (n).
        public var sampleSize: FP?
        /// Sample standard deviation (s).
        public var stdDev: FP?
        /// Standard error s/√n.
        public var stdErr: FP?
        /// Degrees of freedom (n − 1).
        public var df: FP?
        /// Decision flag: mean ≥ μ₀ at α.
        public var meanGTEtestValue: Bool? = false
        /// Decision flag: mean ≤ μ₀ at α.
        public var meanLTEtestValue: Bool? = false
        /// Decision flag: mean = μ₀ (fail to reject) at α.
        public var meanEQtestValue: Bool? = false
        
        /// Human-readable multiline summary of results.
        public var description: String {
            get {
                var descr = String()
                if let p1 = self.p1Value, let p2 = self.p2Value, let t = self.tStat, let ms = self.mean, let m0 = self.mean0, let n = self.sampleSize, let cv90 = self.cv90Pct, let cv95 = self.cv95Pct, let cv99 = self.cv99Pct, let diff = self.difference, let sd = self.stdDev, let se = self.stdErr, let sdf = self.df {
                    descr.append("ONE-SAMPLE T-TEST\n")
                    descr.append("*****************\n")
                    descr.append("sample size: \(niceNumber(n))\n")
                    descr.append("degrees of freedom: \(niceNumber(sdf))\n")
                    descr.append("mean 0: \(niceNumber(m0))\n")
                    descr.append("sample mean: \(niceNumber(ms))\n")
                    descr.append("difference: \(niceNumber(diff))\n")
                    descr.append("sample sd: \(niceNumber(sd))\n")
                    descr.append("standard error: \(niceNumber(se))\n")
                    descr.append("t value: \(niceNumber(t))\n")
                    descr.append("critical value for alpha = 0.01: \(niceNumber(cv99))\n")
                    descr.append("critical value for alpha = 0.05: \(niceNumber(cv95))\n")
                    descr.append("critical value for alpha = 0.10: \(niceNumber(cv90))\n")
                    descr.append("one sided p-value: \(niceNumber(p1))\n")
                    descr.append("two sided p-value: \(niceNumber(p2))\n")
                }
                return descr
            }
        }
    }
    
    /// Results for one-way ANOVA.
    ///
    /// Contains F statistic, p-value, critical value, homoscedasticity test p-values (Bartlett, Levene),
    /// sums of squares, mean squares, and degrees of freedom.
    public struct OneWayANOVAResult<FP: RealLike>: CustomStringConvertible {
        /// Two-tailed p-value for the omnibus F test (H₀: all group means equal).
        public var p2Value: FP?
        /// Observed F statistic = MS_Treatment / MS_Error.
        public var FStatistic: FP?
        /// p-value of Bartlett's test for homoscedasticity.
        public var pBartlett: FP?
        /// Significance level used for the F critical value.
        public var alpha: FP?
        /// Decision flag: fail to reject equal means at α (true) or reject (false).
        public var meansEQUAL: Bool?
        /// F critical value F_{1−α; df₁, df₂}.
        public var cv: FP?
        /// p-value of Levene's test for homoscedasticity.
        public var pLevene: FP?
        /// Total sum of squares (SS_Tot).
        public var SSTotal: FP?
        /// Error (within-group) sum of squares (SS_E).
        public var SSError: FP?
        /// Treatment (between-group) sum of squares (SS_T).
        public var SSTreatment: FP?
        /// Mean square error MS_E = SS_E / df₂.
        public var MSError: FP?
        /// Mean square treatment MS_T = SS_T / df₁.
        public var MSTreatment: FP?
        /// Error degrees of freedom df₂ = N − k.
        public var dfError: FP?
        /// Treatment degrees of freedom df₁ = k − 1.
        public var dfTreatment: FP?
        /// Total degrees of freedom N − 1.
        public var dfTotal: FP?
        /// Post Hoc Tests Result
        public var postHocTests: [PostHoc<FP>.PostHocTestSummary]?
        /// Post hoc test type
        public var postHocTestType: PostHoc<FP>.PostHocTestType?
        /// Human-readable multiline summary of results.
        public var description: String {
            get {
                var descr = String()
                if let f = self.FStatistic, let pB = self.pBartlett, let cv = self.cv, let p2 = self.p2Value, let pL = self.pLevene, let SST = SSTotal, let SSE = SSError, let SSF = SSTreatment, let MSE = MSError, let MST = MSTreatment, let dft = dfTotal, let dfF = dfTreatment, let dfe = dfError, let ph = self.postHocTests, let pht = self.postHocTestType {
                    descr.append("**********************************\n")
                    descr.append("*      ONE WAY ANOVA RESULTS     *\n")
                    descr.append("**********************************\n")
                    descr.append("f: \(niceNumber(f))\n")
                    descr.append("p-value Bartlett test: \(niceNumber(pB))\n")
                    descr.append("p-value Levene test: \(niceNumber(pL))\n")
                    descr.append("critical value: \(niceNumber(cv))\n")
                    descr.append("two sided p-value: \(niceNumber(p2))\n")
                    descr.append("SS_total: \(niceNumber(SST))\n")
                    descr.append("SS_error: \(niceNumber(SSE))\n")
                    descr.append("SS_Treatment: \(niceNumber(SSF))\n")
                    descr.append("df_total: \(niceNumber(dft))\n")
                    descr.append("df_error: \(niceNumber(dfe))\n")
                    descr.append("df_Treatment: \(niceNumber(dfF))\n")
                    descr.append("MS_error: \(niceNumber(MSE))\n")
                    descr.append("MS_treatment: \(niceNumber(MST))\n")
                    descr.append("**********************************\n")
                    switch pht {
                    case .tukeyKramer: descr.append(">>> Tukey-Kramer post hoc test\n")
                    case .scheffe: descr.append(">>> Scheffe post hoc test\n")
                    }
                    for phr in ph {
                        descr.append("\(phr.row), diff: \(niceNumber(phr.meanDiff)), Q: \(niceNumber(phr.testStat)), p: \(niceNumber(phr.pValue))\n")
                    }
                }
                return descr
            }
        }
    }
    
}
