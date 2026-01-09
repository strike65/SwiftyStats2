//
//  Created by VT on 13.12.25.
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

/// Public API for Anderson-Darling goodness-of-fit testing,
/// including normality testing with Stephens correction and p-value approximation.
extension Inferential.HypothesisTesting.GoodnessOfFitTests {
    
    /// Performs the Anderson–Darling normality test (with Stephens correction) on raw data.
    ///
    /// - Parameters:
    ///   - data: Sample observations.
    ///   - alpha: Significance level in (0,1).
    ///
    /// - Discussion:
    ///   Null hypothesis H₀: data ~ Normal(μ,σ²) with μ,σ estimated from data.
    ///   Alternative H₁: data not normal.
    ///
    ///   The Anderson–Darling statistic A² is defined as:
    ///   ```
    ///   A² = -n - (1/n) * sum_{i=1}^n ( (2i-1)[ln F(X_(i)) + ln(1 - F(X_(n+1-i)))] )
    ///   ```
    ///   where F is the hypothesized CDF.
    ///
    ///   Stephens (1974) adjustment for estimated parameters:
    ///   ```
    ///   A²* = A² * (1 + 0.75/n + 2.25/n²)
    ///   ```
    ///
    ///   P-value is approximated by a piecewise function from Stephens (1974).
    ///
    ///   Returns `nil` if preconditions fail: invalid alpha or insufficient data.
    ///
    /// - Returns: `ADTestResult` on success, or `nil` if preconditions not met.
    /// - Throws: If Normal distribution creation or CDF evaluation fails.
    public static func andersonDarlingNormalityTest(data: [T], alpha: T) throws -> ADTestResult<T>? {
        let ex = SSExamine<T,T>.init(usingArray: data, name: nil, characterSet: nil)
        return try andersonDarlingNormalityTest(data: ex, alpha: alpha)
    }
    
    /// Performs the Anderson–Darling normality test using an `SSExamine` summary.
    ///
    /// - Parameters:
    ///   - data: Must provide mean, sample standard deviation, and items as numeric array.
    ///   - alpha: Significance level in (0,1).
    ///
    /// - Discussion:
    ///   Uses estimated mean and sample standard deviation from `data`.
    ///   Converts items to numeric array, sorts them, maps through Normal CDF,
    ///   computes Anderson–Darling A² statistic, applies Stephens correction,
    ///   then computes p-value using piecewise approximation.
    ///
    ///   Decision rule: reject H₀ when p-value < alpha; `isNormal` = (p-value >= alpha).
    ///
    ///   Returns `nil` if sample size < 2, missing statistics, or invalid alpha.
    ///
    /// - Returns: `ADTestResult` if successful, or `nil` if preconditions not met.
    /// - Throws: If Normal distribution creation or CDF evaluation fails.
    public static func andersonDarlingNormalityTest(data: SSExamine<T, T>, alpha: T) throws -> ADTestResult<T>? {
        guard data.sampleSize >= 2 && alpha.liesInOpenRange(from: 0, to: 1) else {
            return nil
        }
        guard let estMean = data.arithmeticMean, let estStdDev = data.sampleStandardDeviation, let set: [T] = data.itemsAsNumericArray else {
            return nil
        }
        let sorted: [T] = set.sorted()
        let dist: Distribution.Normal<T> = try .init(mean: estMean, sd: estStdDev)
        let u: [T] = try sorted.map { x in
            try dist.cdf(x)
        }
        let rawAD: T = statisticFromUniform(u)
        let n = T(data.sampleSize)
        // Stephens (1974) adjustment for estimated mean and variance under normality.
        let adjustedAD: T = rawAD * (.one + T(0.75) / n + T(2.25) / (n * n))
        let pValue: T = normalityPValue(forAdjustedAD: adjustedAD)
        var res: ADTestResult<T> = ADTestResult<T>.init()
        res.pValue = pValue
        res.AD = adjustedAD
        res.sampleSize = data.sampleSize
        res.stdDev = estStdDev
        res.mean = estMean
        res.variance = estStdDev * estStdDev
        res.isNormal = pValue >= alpha
        return res
    }
    
    /// Computes the Anderson-Darling statistic A^2 for a given sample and
    /// a hypothesized continuous distribution with CDF F.
    ///
    /// - Parameters:
    ///   - sample: Observations X_1, ..., X_n
    ///   - cdf:    Cumulative distribution function F(x) under H0
    /// - Returns: Anderson-Darling statistic A^2
    private static func statistic(
        sample: [T],
        cdf: (T) -> T
    ) -> T {
        precondition(!sample.isEmpty, "Sample must not be empty")
        
        // Sort sample, then map through CDF.
        // Since F is non-decreasing, the order of F(X_(i)) matches that of X_(i).
        let sortedX = sample.sorted()
        let u = sortedX.map { x in
            cdf(x)
        }
        
        // Reuse the uniform-based implementation.
        return statisticFromUniform(u)
    }
    
    /// Computes the Anderson-Darling statistic A^2 for a given sample of
    /// values assumed to follow Uniform(0, 1) under the null hypothesis.
    ///
    /// The input array is not required to be sorted; it will be sorted internally.
    /// Small epsilon-clamping is used to avoid log(0).
    private static func statisticFromUniform(_ uniformSample: [T]) -> T {
        let n = uniformSample.count
        precondition(n > 0, "Sample must contain at least one observation")
        
        // Sort the sample (ascending)
        let sorted = uniformSample.sorted()
        let nf = T(n)
        
        // Avoid log(0) and log(1-0) by clamping
        let eps = T.ulpOfOne
        
        var sumTerms: T = .zero
        
        // i runs from 1 to n (inclusive) in the formula
        for i in 1...n {
            let uLower = T.maximum(eps, min(.one - eps, sorted[i - 1]))
            let uUpper = T.maximum(eps, min(.one - eps, sorted[n - i]))
            
            let coeff: T = T(2 * i - 1)
            let logLower: T = T.log(uLower)
            let logUpperComp: T = T.log(.one - uUpper)
            let term: T = coeff * (logLower + logUpperComp)
            sumTerms += term
        }
        
        let a2 = -nf - sumTerms / nf
        return a2
    }
    
    /// CDF of the Anderson-Darling statistic A^2 for sample size n.
    /// Returns approximately P(A_n <= a2) according to
    /// Marsaglia & Marsaglia (2004), using the quick adinf(z) approximation
    /// plus the finite-n error correction errfix(n, x).
    private static func adcdf(a2: T, sampleSize n: Int) -> T {
        precondition(n > 0, "sampleSize must be positive")
        if a2 <= .zero {
            return .zero
        }
        
        let x = asymptoticAD(a2)               // adinf(z)
        let correction = finiteSampleError(x: x, n: n)
        var p = x + correction                 // ADinf(z) + errfix(n, ADinf(z))
        
        // numerical safety: clamp to [0, 1]
        if p < .zero { p = .zero }
        if p > .one { p = .one }
        return p
    }
    
    /// P-value approximation for the normality AD test with estimated mean/variance.
    /// Uses Stephens (1974) piecewise exponential fits adopted by nortest/ad.test.
    private static func normalityPValue(forAdjustedAD a2: T) -> T {
        let p: T
        if a2 < T(0.2) {
            let exponent = T(-13.436) + T(101.14) * a2 - T(223.73) * a2 * a2
            p = .one - T.exp(exponent)
        } else if a2 < T(0.34) {
            let exponent = T(-8.318) + T(42.796) * a2 - T(59.938) * a2 * a2
            p = .one - T.exp(exponent)
        } else if a2 < T(0.6) {
            p = T.exp(T(0.9177) - T(4.279) * a2 - T(1.38) * a2 * a2)
        } else if a2 < T(13) {
            p = T.exp(T(1.2937) - T(5.709) * a2 + T(0.0186) * a2 * a2)
        } else {
            p = .zero
        }
        if p < .zero { return .zero }
        if p > .one { return .one }
        return p
    }
    
    /// Asymptotic CDF of the Anderson-Darling statistic, adinf(z).
    /// Valid for z > 0. For z <= 0, returns 0.
    private static func asymptoticAD(_ z: T) -> T {
        if z <= .zero {
            return .zero
        }
        if z < .two {
            // 0 < z < 2: z^(-1/2) * exp(-1.2337141 / z) * polynomial(z)
            let invZ: T = .one / z
            let rootZ: T = T.sqrt(z)
            
            // Horner polynomial coefficients split to help the type-checker
            let c5: T = T(-0.00168691)
            let c4: T = T(0.0116720) + z * c5
            let c3: T = T(-0.0347962) + z * c4
            let c2: T = T(-0.0649821) + z * c3
            let c1: T = T(0.247105) + z * c2
            let poly: T = T(2.00012) + z * c1
            
            let expo: T = T.exp(T(-1.2337141) * invZ)
            return expo * poly / rootZ
        } else {
            // z >= 2: exp( -exp( p(z) ) ) with p(z) as a 5th degree polynomial
            let d5: T = T(-0.0003146)
            let d4: T = T(0.008056) + z * d5
            let d3: T = T(-0.082433) + z * d4
            let d2: T = T(0.43424) + z * d3
            let d1: T = T(-2.30695) + z * d2
            let p: T = T(1.0776) + z * d1
            let inner: T = T.exp(p)
            return T.exp(-inner)
        }
    }
    
    /// Finite-sample correction errfix(n, x) as described by Marsaglia & Marsaglia
    /// and summarized e.g. in Retzer (2025).
    private static func finiteSampleError(x: T, n: Int) -> T {
        let nn: T = T(n)
        
        // Outside (0,1) there is nothing meaningful to correct.
        if x <= .zero || x >= .one {
            return .zero
        }
        
        let c = T(0.01265) + T(0.1757) / nn
        
        if x < c {
            // Region 1: x < c(n)
            let u = x / c
            // g1(u) = sqrt(u) * (1 - u) * (49 u - 102)
            let g1: T = T.sqrt(u) * (.one - u) * (T(49.0) * u - T(102.0))
            
            // scale = 0.0037/n^3 + 0.00078/n^2 + 0.00006/n
            let n2 = nn * nn
            let n3 = n2 * nn
            let scale = T(0.0037) / n3 + T(0.00078) / n2 + T(0.00006) / nn
            
            return g1 * scale
            
        } else if x < T(0.8) {
            // Region 2: c(n) <= x < 0.8
            let u = (x - c) / (T(0.8) - c)
            
            // g2(u) as 5th degree polynomial in Horner form with intermediates
            // g2(u) = -0.00022633
            //         + 6.54034 u
            //         - 14.6538 u^2
            //         + 14.458 u^3
            //         - 8.259 u^4
            //         + 1.91864 u^5
            let g2_u5: T = T(1.91864)
            let g2_u4: T = T(-8.259) + u * g2_u5
            let g2_u3: T = T(14.458) + u * g2_u4
            let g2_u2: T = T(-14.6538) + u * g2_u3
            let g2_u1: T = T(6.54034) + u * g2_u2
            let g2: T = T(-0.00022633) + u * g2_u1
            
            // scale = 0.04213/n + 0.01365/n^2
            let n2 = nn * nn
            let scale = T(0.04213) / nn + T(0.01365) / n2
            
            return g2 * scale
            
        } else {
            // Region 3: x >= 0.8
            let v = x
            
            // g3(v) = -130.2137 + 745.2337 v - 1705.091 v^2
            //         + 1950.646 v^3 - 1116.36 v^4 + 255.7844 v^5
            let g3_v5: T = T(255.7844)
            let g3_v4: T = T(-1116.36) + v * g3_v5
            let g3_v3: T = T(1950.646) + v * g3_v4
            let g3_v2: T = T(-1705.091) + v * g3_v3
            let g3_v1: T = T(745.2337) + v * g3_v2
            let g3: T = T(-130.2137) + v * g3_v1
            
            return g3 / nn
        }
    }
    
    /// Results of the Anderson–Darling test for normality.
    ///
    /// Contains the computed test statistic, p-value, sample statistics,
    /// and a flag indicating whether the null hypothesis (normality) is accepted at the given significance level.
    ///
    /// The `description` property prints a formatted human-readable report.
    public struct ADTestResult<FP: RealLike>: CustomStringConvertible, Codable  {
        /// p-value in [0,1], defined as P(A²* >= observed | H₀).
        /// Larger values indicate greater consistency with normality.
        public var pValue: FP?
        /// The adjusted Anderson–Darling statistic A²* (Stephens-corrected).
        public var AD: FP?
        /// Sample size n.
        public var sampleSize: Int?
        /// Sample standard deviation (used for Normal parameter estimate).
        public var stdDev: FP?
        /// Variance, i.e. square of stdDev.
        public var variance: FP?
        /// Sample mean.
        public var mean: FP?
        /// Decision flag for the test at the provided alpha (true when pValue >= alpha).
        public var isNormal: Bool?
        
        /// Returns a printed report of the Anderson–Darling test result,
        /// including test statistic, p-value, sample statistics, and normality decision.
        public var description: String {
            get {
                var descr = String()
                descr.append("ANDERSON-DARLING TEST FOR NORMALITY\n")
                descr.append("***********************************\n")
                if let ad = self.AD, let sd = self.stdDev, let v = self.variance, let m = self.mean, let isn = self.isNormal, let p = self.pValue, let n = self.sampleSize {
                    descr.append("normality: \(isn)\n")
                    descr.append("p: \(niceNumber(p))\n")
                    descr.append("AD: \(niceNumber(ad))\n")
                    descr.append("mean: \(niceNumber(m))\n")
                    descr.append("sd: \(niceNumber(sd))\n")
                    descr.append("var: \(niceNumber(v))\n")
                    descr.append("n: \(n)\n")
                }
                return descr
            }
        }
    }
}
