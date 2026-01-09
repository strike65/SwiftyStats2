//
//  Created by VT on 03.12.25.
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

/// Exact binomial hypothesis testing helpers and overloads.

import SwiftyStatsPrelude

extension Inferential.HypothesisTesting.ParametricTests {
    /// Exact binomial test p-values for a given number of successes.
    ///
    /// Computes lower-tail, upper-tail, and two-sided exact p-values under
    /// a Binomial(n, p0) null model using the exact distribution:
    /// - lower-tail: P(X ≤ x)
    /// - upper-tail: P(X ≥ x) = 1 − P(X ≤ x − 1)
    /// - two-sided: Sum of probabilities of all outcomes with probability ≤ P(X = x)
    ///   (the “equal-or-more-extreme by probability” definition).
    ///
    /// Implementation details:
    /// - Builds a ``SwiftyBoost/Distribution/Binomial`` with `n` and `p0`.
    /// - Uses its `cdf` for tails and `pdf` to accumulate the two-sided probability mass
    ///   over all k with pdf(k) ≤ pdf(x).
    ///
    /// - Parameters:
    ///   - x: The observed number of successes (0...n).
    ///   - n: The number of independent Bernoulli trials (n > 0).
    ///   - p0: The null success probability, must lie in (0, 1).
    /// - Returns: A tuple of p-values `(less, greater, twoSided)` when inputs are valid; otherwise `nil`.
    /// - Throws: Rethrows errors from distribution evaluation (e.g., invalid parameters).
    /// - SeeAlso: ``binomialTest(successes:numberOfTrials:p0:alternative:)``
    public static func binomialTest(successes x: Int, numberOfTrials n: Int, p0: T) throws -> (less: T, greater: T, twoSided: T)? {
        guard p0.liesInOpenRange(from: T.zero, to: T.one) else {
            return nil
        }
        guard n > 0 else {
            return nil
        }
        guard x >= 0 && x <= n else {
            return nil
        }
        var p_ValueLess: T, p_ValueGreater: T, p_ValueTwoSided: T
        var dist: Distribution.Binomial<T>
        dist = try Distribution.Binomial(numberOfTrials: T(n), probabibilityOfSuccess: p0)
        p_ValueLess = try dist.cdf(T(x))
        p_ValueGreater = try dist.cdf(T(x - 1))
        p_ValueGreater = T.one - p_ValueGreater
        p_ValueTwoSided = .zero
        let p_obs = try dist.pdf(T(x))
        var p_k: T = .zero
        for k in 0...n {
            p_k = try dist.pdf(T(k))
            if p_k <= p_obs {
                p_ValueTwoSided += p_k
            }
        }
        return (less: p_ValueLess, greater: p_ValueGreater, twoSided: p_ValueTwoSided)
    }
    
    /// Exact binomial test p-value for a chosen alternative.
    ///
    /// Convenience wrapper around
    /// ``binomialTest(successes:numberOfTrials:p0:alternative:)`` that returns only the
    /// p-value corresponding to the requested alternative hypothesis.
    ///
    /// - Parameters:
    ///   - x: The observed number of successes (0...n).
    ///   - n: The number of independent Bernoulli trials (n > 0).
    ///   - p0: The null success probability, must lie in (0, 1).
    ///   - alternative: The alternative hypothesis to test against.
    ///     - `.less`: H1: p < p0, returns P(X ≤ x).
    ///     - `.greater`: H1: p > p0, returns P(X ≥ x).
    ///     - `.twoSided`: H1: p ≠ p0, returns the two-sided exact p-value.
    /// - Returns: The requested p-value when inputs are valid; otherwise `nil`.
    /// - Throws: Rethrows errors from distribution evaluation (e.g., invalid parameters).
    public static func binomialTest(successes x: Int, numberOfTrials n: Int, p0: T, alternative: Alternative) throws -> T? {
        switch alternative {
        case .less:
            if let p = try binomialTest(successes: x, numberOfTrials: n, p0: p0) {
                return p.less
            }
            else {
                return nil
            }
        case .greater:
            if let p = try binomialTest(successes: x, numberOfTrials: n, p0: p0) {
                return p.greater
            }
            else {
                return nil
            }
        case .twoSided:
            if let p = try binomialTest(successes: x, numberOfTrials: n, p0: p0) {
                return p.twoSided
            }
            else {
                return nil
            }
        }
    }
    
    /// Exact binomial test over a categorical dataset, selecting a success code.
    ///
    /// This overload accepts a raw array and constructs an ``SSExamine`` to compute
    /// counts and proportions. It then computes:
    /// - exact p-value for the requested ``Alternative``
    /// - Jeffreys and Clopper–Pearson confidence intervals for the success probability
    ///
    /// Success is defined by `successCode`, and the success probability is estimated as
    /// `p̂ = nSuccess / n`. Confidence intervals are computed at the container’s current
    /// alpha level (as configured elsewhere in the module).
    ///
    /// - Parameters:
    ///   - data: A sample of categorical observations.
    ///   - p0: The null success probability to test against; must lie in (0, 1).
    ///   - characterSet: Optional character filter used only when `S == String`. Ignored otherwise.
    ///   - sc: The value in `data` to treat as “success”.
    ///   - alternative: The alternative hypothesis.
    /// - Returns: A ``BinomialTestResult`` aggregating counts, p-value, and confidence intervals; or `nil` if inputs are invalid.
    /// - Throws: Rethrows any errors from underlying distribution or container operations.
    /// - SeeAlso: ``binomialTest(data:p0:characterSet:successCode:alternative:)``
    public static func binomialTest<S: Codable & Comparable & Hashable & Sendable>(data: Array<S>, p0: T, characterSet: CharacterSet?, successCode sc: S, alternative: Alternative) throws -> BinomialTestResult<S, T>? {
        guard !p0.isNaN, p0.isFinite, p0.liesInOpenRange(from: .zero, to: .one) else {
            return nil
        }
        guard data.count > 2 else { return nil }
        let examine: SSExamine<S, T> = .init(usingArray: data, name: nil, characterSet: characterSet)
        return try binomialTest(data: examine, p0: p0, characterSet: characterSet, successCode: sc, alternative: alternative)
    }
    
    /// Exact binomial test using an existing ``SSExamine`` container.
    ///
    /// Interprets `successCode` as success, computes `nSuccess`, `nFailure`, `n`,
    /// and the observed success fraction `p̂ = nSuccess / n`. It then:
    /// - Computes a one-sided or two-sided exact p-value using `p̂` as the test probability argument (see note).
    /// - Computes Jeffreys and Clopper–Pearson confidence intervals using
    ///
    /// Note:
    /// - The p-value call uses `pSuccess` (the observed fraction) as the probability argument.
    ///   If you intend to test against a specific null `p0`, you likely want to pass `p0` instead.
    ///   Review your intended hypothesis to ensure this meets your needs.
    ///
    /// - Parameters:
    ///   - data: The examine container holding the sample.
    ///   - p0: The null success probability for reporting and CI context; must lie in (0, 1).
    ///   - characterSet: Optional character filter used only when `S == String`. Ignored otherwise.
    ///   - sc: The value to treat as “success”.
    ///   - alternative: The alternative hypothesis.
    /// - Returns: A ``BinomialTestResult`` aggregating counts, p-value, estimated probabilities, and confidence intervals; or `nil` if inputs are invalid.
    /// - Throws: Rethrows any errors from underlying distribution or container operations.
    public static func binomialTest<S: Codable & Comparable & Hashable & Sendable>(data: SSExamine<S, T>, p0: T, characterSet: CharacterSet?, successCode sc: S, alternative: Alternative) throws -> BinomialTestResult<S, T>? {
        guard !p0.isNaN, p0.isFinite, p0.liesInOpenRange(from: .zero, to: .one) else {
            return nil
        }
        guard !data.isEmpty else { return nil }
        let succcess = data.frequency(sc)
        let failure = data.count - succcess
        let n = data.sampleSize
        let pSuccess: T = T(succcess) / T(n)
        var cintJeffreys: ConfidenceInterval<T>
        var cintClopperPearson: ConfidenceInterval<T>
        var uB: T, lB: T
        switch alternative {
        case .less:
            uB = Distribution.Binomial<T>.findUpperBoundOnP(nTrials: n, nSuccesses: succcess, proposedSuccessFraction: pSuccess, useJeffreys: true)
            cintJeffreys = .init(lowerBound: .zero, upperBound: uB)
            uB = Distribution.Binomial<T>.findUpperBoundOnP(nTrials: n, nSuccesses: succcess, proposedSuccessFraction: pSuccess, useJeffreys: false)
            cintClopperPearson = .init(lowerBound: .zero, upperBound: uB)
        case .greater:
            lB = Distribution.Binomial<T>.findLowerBoundOnP(nTrials: n, nSuccesses: succcess, proposedSuccessFraction: pSuccess, useJeffreys: true)
            cintJeffreys = .init(lowerBound: lB, upperBound: .one)
            lB = Distribution.Binomial<T>.findLowerBoundOnP(nTrials: n, nSuccesses: succcess, proposedSuccessFraction: pSuccess, useJeffreys: false)
            cintClopperPearson = .init(lowerBound: lB, upperBound: .one)
        case .twoSided:
            lB = Distribution.Binomial<T>.findLowerBoundOnP(nTrials: n, nSuccesses: succcess, proposedSuccessFraction: pSuccess, useJeffreys: true)
            uB = Distribution.Binomial<T>.findUpperBoundOnP(nTrials: n, nSuccesses: succcess, proposedSuccessFraction: pSuccess, useJeffreys: true)
            cintJeffreys = .init(lowerBound: lB, upperBound: uB)
            lB = Distribution.Binomial<T>.findLowerBoundOnP(nTrials: n, nSuccesses: succcess, proposedSuccessFraction: pSuccess, useJeffreys: false)
            uB = Distribution.Binomial<T>.findUpperBoundOnP(nTrials: n, nSuccesses: succcess, proposedSuccessFraction: pSuccess, useJeffreys: false)
            cintClopperPearson = .init(lowerBound: .zero, upperBound: uB)
        }
        if let pValue: T = try Parametric<T>.binomialTest(successes: Int(succcess), numberOfTrials: Int(n), p0: pSuccess, alternative: alternative) {
            let res: BinomialTestResult<S, T> = .init(nTrials: n, nSuccess: succcess, nFailure: failure, pValueExact: pValue, probSuccess: pSuccess, probFailure: T(failure) / T(n), probTest: p0, successCode: sc, confIntJeffreys: cintJeffreys, confIntClopperPearson: cintClopperPearson)
            return res
        }
        else {
            return nil
        }
        
    }
    
    /// Binomial test results
    public struct BinomialTestResult<S: Comparable & Hashable & Codable, FP: RealLike>: CustomStringConvertible, Codable {
        /// number of trials
        public var nTrials: Int?
        /// number of successes
        public var nSuccess: Int?
        /// number of failures
        public var nFailure: Int?
        /// one-sided p-value (exact)
        public var pValueExact: FP?
        /// probability for success
        public var probSuccess: FP?
        /// probability for failure
        public var probFailure: FP?
        /// test probability
        public var probTest: FP?
        /// success id
        public var successCode: S?
        /// 1 - alpha confidence interval (Jeffreys)
        public var confIntJeffreys: ConfidenceInterval<FP>?
        /// 1 - alpha confidence interval (Clopper/Pearson)
        public var confIntClopperPearson: ConfidenceInterval<FP>?
        /// Returns a description
        public var description: String {
            get {
                var descr = String()
                if let nt = self.nTrials, let ns = self.nSuccess, let nf = self.nFailure, let pe = self.pValueExact, let ps = self.probSuccess, let pf = self.probFailure, let pt = self.probTest, let sc = self.successCode, let jeff = self.confIntJeffreys, let clopper = self.confIntClopperPearson {
                    descr.append("BINOMIAL TEST\n")
                    descr.append("*************\n")
                    descr.append("p-value exact: \(niceNumber(pe))\n")
                    descr.append("count of trials: \(nt)\n")
                    descr.append("count of successes: \(ns)\n")
                    descr.append("count of failures: \(nf)\n")
                    descr.append("prob. for success: \(niceNumber(ps))\n")
                    descr.append("prob. for failure: \(niceNumber(pf))\n")
                    descr.append("test prob.: \(niceNumber(pt))\n")
                    descr.append("succes coded as: \(sc)\n")
                    descr.append("(1-alpha) CI (Jeffreys): \(jeff)\n")
                    descr.append("(1-alpha) CI (Clopper-Pearson): \(clopper)\n")
                }
                return descr
            }
        }
    }
    
}
