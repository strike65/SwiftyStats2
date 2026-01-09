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

/// Public entropy and information-theoretic measures for SSExamine datasets.

import SwiftyStatsPrelude
import SwiftyBoost

// MARK: - Generic entropy class

/// Information-theoretic measures over `SSExamine` samples with interchangeable estimators.
///
/// Each API mirrors a formal definition—such as Shannon `H = -sum p log p`, Renyi `H_alpha`,
/// or KL divergence—while letting callers choose discrete plug-in, histogram, KDE, or
/// Miller-Madow corrections. The shared container keeps the probability mass function
/// and empirical distribution aligned so entropy, cross-entropy, and mutual information
/// are computed on consistent support and scaling.
extension SSExamine {
    /// Computes Shannon entropy for the receiver using the requested estimation method.
    ///
    /// - Parameters:
    ///   - type: The Shannon entropy method to use (e.g. histogram, frequency, KDE, Miller–Madow, discrete).
    ///   - base: The logarithm base for the entropy. Defaults to 2 (bits).
    ///   - bins: The binning strategy to use for histogram-based methods. Defaults to `.auto`.
    /// - Returns: The estimated Shannon entropy in the given base, or `nan` if an error occurs.
    ///
    /// Interpretation hints:
    /// - Higher values indicate more unpredictability (greater diversity in the symbol distribution).
    /// - Units depend on `base` (e.g. base 2 = bits, base e = nats).
    /// - For discrete methods, the maximum entropy is log_base(|support|) and is reached by a uniform distribution.
    /// - Histogram/KDE-based values depend on binning/bandwidth and on the units/scale of the data; compare only when methods and parameters are consistent.
    public func shannonEntropy(type: ShannonMethod<FP>, base: FP = 2, bins: BinningStrategy = .auto) -> FP {
        return withLock { () -> FP in
            do {
                return try EntropyMeasures(data: self).shannonEntropy(method: type, bins: bins, base: base)
            }
            catch _ {
                return FP.nan
            }
        }
    }

    /// Computes Rényi entropy for the receiver.
    ///
    /// - Parameters:
    ///   - alpha: The Rényi order α. Must be >= 0 and != 1.
    ///   - base: The logarithm base for the entropy. Defaults to 2 (bits).
    /// - Returns: The Rényi entropy in the given base, or `nan` if an error occurs.
    ///
    /// Interpretation hints:
    /// - α < 1 emphasizes the contribution of rare events (diversity); α > 1 emphasizes common events (concentration).
    /// - α = 0 gives Hartley entropy; α = 2 gives collision entropy; α → ∞ gives min-entropy.
    /// - Compare values only at the same α and base; different α capture different aspects of uncertainty.
    public func renyiEntropy(alpha: FP, base: FP = 2) -> FP {
        return withLock { () -> FP in
            do {
                return try EntropyMeasures(data: self).renyiEntropy(alpha: alpha, base: base)
            }
            catch _ {
                return FP.nan
            }
        }
    }

    /// Computes Hartley entropy (Rényi α = 0), using the number of unique symbols.
    ///
    /// - Parameter base: The logarithm base. Defaults to 2 (bits).
    /// - Returns: Hartley entropy in the given base, or `nan` if an error occurs.
    ///
    /// Interpretation hints:
    /// - Equals log_base(number of distinct values). Ignores the actual frequencies.
    /// - Useful as an upper bound on discrete Shannon entropy for a given support size.
    public func hartleyEntropy(base: FP = 2) -> FP {
        return withLock {
            do {
                return try EntropyMeasures(data: self).hartleyEntropy(base: base)
            } catch {
                return FP.nan
            }
        }
    }

    /// Computes collision entropy (Rényi α = 2) from discrete frequencies.
    ///
    /// - Parameter base: The logarithm base. Defaults to 2 (bits).
    /// - Returns: Collision entropy in the given base, or `nan` if an error occurs.
    ///
    /// Interpretation hints:
    /// - Related to the probability of drawing two identical symbols (index of coincidence).
    /// - Lower values indicate more concentration (dominant symbols); higher values indicate more uniformity.
    public func collisionEntropy(base: FP = 2) -> FP {
        return withLock {
            do {
                return try EntropyMeasures(data: self).collisionEntropy(base: base)
            } catch {
                return FP.nan
            }
        }
    }

    /// Computes min-entropy (Rényi α → ∞) from discrete frequencies.
    ///
    /// - Parameter base: The logarithm base. Defaults to 2 (bits).
    /// - Returns: Min-entropy in the given base, or `nan` if an error occurs.
    ///
    /// Interpretation hints:
    /// - Focuses on the most probable outcome: H_min = -log_base(max p_i).
    /// - Captures worst-case unpredictability; useful in security/cryptography contexts.
    public func minEntropy(base: FP = 2) -> FP {
        return withLock {
            do {
                return try EntropyMeasures(data: self).minEntropy(base: base)
            } catch {
                return FP.nan
            }
        }
    }

    /// Computes Tsallis entropy of order q for discrete frequencies.
    ///
    /// - Parameter q: Tsallis parameter q (> 0).
    /// - Returns: Tsallis entropy in nats (base e), or `nan` if an error occurs.
    ///
    /// Interpretation hints:
    /// - Non-additive generalization of Shannon entropy.
    /// - q < 1 gives more weight to rare events; q > 1 gives more weight to common events.
    /// - Not directly comparable to Shannon in magnitude; use consistently at the same q for comparisons.
    public func tsallisEntropy(q: FP) -> FP {
        return withLock {
            do {
                return try EntropyMeasures(data: self).tsallisEntropy(q: q)
            } catch {
                return FP.nan
            }
        }
    }

    /// Computes approximate entropy (ApEn) for a numeric time series.
    ///
    /// - Parameters:
    ///   - m: Embedding dimension (pattern length).
    ///   - r: Tolerance parameter; if `normalized == true`, interpreted as a fraction of the series' SD.
    ///   - normalized: If true, scale `r` by the sample standard deviation.
    /// - Returns: ApEn value or `nan` if not applicable or an error occurs.
    ///
    /// Interpretation hints:
    /// - Larger ApEn indicates greater irregularity/complexity in the time series; smaller values indicate more regularity.
    /// - Sensitive to record length and to the choice of (m, r); compare only when parameters and preprocessing are consistent.
    /// - Requires truly metric numeric data; categorical fallback is not meaningful here.
    public func approximateEntropy(m: Int, r: FP, normalized: Bool = true) -> FP {
        return withLock {
            do {
                return try EntropyMeasures(data: self).approximateEntropy(m: m, r: r, normalized: normalized)
            } catch {
                return FP.nan
            }
        }
    }

    /// Computes sample entropy (SampEn) for a numeric time series.
    ///
    /// - Parameters:
    ///   - m: Embedding dimension (pattern length).
    ///   - r: Tolerance parameter; if `normalized == true`, interpreted as a fraction of the series' SD.
    ///   - normalized: If true, scale `r` by the sample standard deviation.
    /// - Returns: SampEn value or `nan` if not applicable or an error occurs.
    ///
    /// Interpretation hints:
    /// - Like ApEn, higher SampEn indicates greater irregularity; unlike ApEn, it is designed to reduce bias and is less sensitive to data length.
    /// - Values can be infinite when there are no matches at scale m+1; this indicates very high irregularity or too-strict parameters.
    /// - Compare only when (m, r) and preprocessing are consistent.
    public func sampleEntropy(m: Int, r: FP, normalized: Bool = true) -> FP {
        return withLock {
            do {
                return try EntropyMeasures(data: self).sampleEntropy(m: m, r: r, normalized: normalized)
            } catch {
                return FP.nan
            }
        }
    }

    /// Computes permutation entropy using ordinal patterns of the sequence.
    ///
    /// - Parameters:
    ///   - order: The permutation order (pattern length), typically 3...7.
    ///   - delay: The delay (time step) between successive elements in the pattern.
    ///   - base: Logarithm base for the entropy. Defaults to 2 (bits).
    ///   - normalized: If true, returns entropy normalized by log(order!).
    /// - Returns: Permutation entropy or `nan` on failure.
    ///
    /// Interpretation hints:
    /// - Measures the complexity of the ordering of values; robust to monotonic transformations and amplitude scaling.
    /// - When normalized, values are in [0, 1]: near 0 indicates highly regular order patterns; near 1 indicates random-like orderings.
    /// - Results depend on `order` and `delay`; choose them based on the time scale of interest and compare consistently.
    public func permutationEntropy(order: Int, delay: Int = 1, base: FP = 2, normalized: Bool = true) -> FP {
        return withLock {
            do {
                return try EntropyMeasures(data: self).permutationEntropy(order: order, delay: delay, base: base, normalized: normalized)
            } catch {
                return FP.nan
            }
        }
    }

    /// Estimates differential entropy via the Kozachenko–Leonenko k-NN method (1D).
    ///
    /// - Parameters:
    ///   - k: Neighbor order (k >= 1).
    ///   - base: Logarithm base for the output. Defaults to 2 (bits).
    /// - Returns: Differential entropy estimate, or `nan` on failure.
    ///
    /// Interpretation hints:
    /// - Differential entropy is not bounded below and can be negative.
    /// - Not invariant to linear scaling/units (e.g., changing units shifts the value); compare only when data are in the same units and similarly scaled.
    /// - Choice of k trades bias and variance; moderate k (e.g., 3–10) is common for 1D.
    public func kozachenkoLeonenkoEntropy(k: Int = 3, base: FP = 2) -> FP {
        return withLock {
            do {
                return try EntropyMeasures(data: self).kozachenkoLeonenkoEntropy(k: k, base: base)
            } catch {
                return FP.nan
            }
        }
    }

    /// Estimates differential entropy via the Vasicek spacing estimator (1D).
    ///
    /// - Parameters:
    ///   - m: Window size; if `nil`, uses ⌊√n⌋.
    ///   - base: Logarithm base for the output. Defaults to 2 (bits).
    /// - Returns: Differential entropy estimate, or `nan` on failure.
    ///
    /// Interpretation hints:
    /// - Like other differential entropy estimators, values depend on the measurement scale and can be negative.
    /// - The window m controls the bias–variance trade-off; too small increases variance, too large increases bias.
    /// - Compare only between datasets with consistent units and similar sample sizes/parameters.
    public func vasicekEntropy(m: Int? = nil, base: FP = 2) -> FP {
        return withLock {
            do {
                return try EntropyMeasures(data: self).vasicekEntropy(m: m, base: base)
            } catch {
                return FP.nan
            }
        }
    }

    /// Computes conditional entropy H(X | condition(X) == true) for a boolean condition over elements.
    ///
    /// - Parameters:
    ///   - condition: Predicate selecting the subset.
    ///   - base: Logarithm base. Defaults to 2 (bits).
    /// - Returns: Conditional entropy of the subset, or `nan` on failure.
    ///
    /// Interpretation hints:
    /// - Measures the remaining uncertainty within the subset defined by `condition`.
    /// - Typically less than or equal to the overall entropy; can reveal heterogeneity across subpopulations.
    public func conditionalEntropy(condition: (SSElement) -> Bool, base: FP = 2) -> FP {
        return withLock {
            do {
                return try EntropyMeasures(data: self).conditionalEntropy(condition: condition, base: base)
            } catch {
                return FP.nan
            }
        }
    }

    /// Estimates mutual information I(X; Y) between two discrete transforms of the elements.
    ///
    /// - Parameters:
    ///   - transform1: Mapping from element to a discrete label for X.
    ///   - transform2: Mapping from element to a discrete label for Y.
    ///   - base: Logarithm base for the result. Defaults to 2 (bits).
    /// - Returns: Estimated mutual information, or `nan` on failure.
    ///
    /// Interpretation hints:
    /// - Non-negative and symmetric; zero indicates statistical independence (under the given discretizations).
    /// - Larger values indicate stronger dependence/shared information.
    /// - Depends on how you discretize/transform; compare only with consistent transforms and bases.
    public func mutualInformation(transform1: (SSElement) -> Int, transform2: (SSElement) -> Int, base: FP = 2) -> FP {
        return withLock {
            do {
                return try EntropyMeasures(data: self).mutualInformation(transform1: transform1, transform2: transform2, base: base)
            } catch {
                return FP.nan
            }
        }
    }

    /// Computes the Kullback–Leibler divergence D_KL(self || other) for discrete sequences.
    ///
    /// - Parameters:
    ///   - other: The comparison dataset Q.
    ///   - base: Logarithm base for the result. Defaults to 2 (bits).
    ///   - locale: Optional locale used when converting elements to numeric for fallback paths.
    /// - Returns: D_KL(self || other), `+∞` if Q assigns zero probability to a symbol with positive probability in P, or `nan` on error.
    ///
    /// Interpretation hints:
    /// - Asymmetric measure of dissimilarity: D_KL(P || Q) ≥ 0, equals 0 iff P and Q are identical distributions.
    /// - Sensitive to support mismatch: if Q assigns zero probability where P does not, the divergence is infinite.
    /// - Useful to quantify how well Q approximates P (e.g., model vs data); compare only with consistent bases.
    public func kullbackLeiblerDivergence(to other: SSExamine<SSElement, FP>, base: FP = 2, locale: Locale? = nil) -> FP {
        // Retrieve raw sequences outside EntropyMeasures to compare two datasets.
        let p = self.itemsAsArray(sorted: .original)
        let q = other.itemsAsArray(sorted: .original)
        do {
            return try EntropyMeasures.kullbackLeiblerDivergence(p: p, q: q, base: base, locale: locale)
        } catch {
            return FP.nan
        }
    }

    /// Computes cross entropy H(P, Q) where P is `self` and Q is `other`, for discrete sequences.
    ///
    /// - Parameters:
    ///   - other: The comparison dataset Q.
    ///   - base: Logarithm base for the result. Defaults to 2 (bits).
    /// - Returns: Cross entropy H(P, Q), `+∞` if Q assigns zero probability to a symbol with positive probability in P, or `nan` on error.
    ///
    /// Interpretation hints:
    /// - H(P, Q) = H(P) + D_KL(P || Q); it measures the expected code length when coding samples from P with a code optimized for Q.
    /// - Lower is better when Q is meant to model P; equals H(P) when Q = P.
    /// - Infinite when Q assigns zero probability to events that occur in P.
    public func crossEntropy(with other: SSExamine<SSElement, FP>, base: FP = 2) -> FP {
        let p = self.itemsAsArray(sorted: .original)
        let q = other.itemsAsArray(sorted: .original)
        do {
            return try EntropyMeasures.crossEntropy(p: p, q: q, base: base)
        } catch {
            return FP.nan
        }
    }

    /// Computes a profile of commonly used entropy and complexity measures.
    ///
    /// The profile may include:
    /// - Discrete Shannon entropy variants
    /// - Rényi family members
    /// - Tsallis entropy
    /// - Time-series entropies where applicable
    /// - Continuous entropy estimators when the data are numeric
    ///
    /// - Returns: A dictionary mapping metric names to values. Returns an empty dictionary on error.
    ///
    /// Interpretation hints:
    /// - Not all metrics are normalized; compare each metric only with the same method, parameters, and base.
    /// - Time-series metrics require sufficient length and truly numeric data; some entries may be absent depending on n and data type.
    /// - Use this as a quick “fingerprint”; for decisions, inspect individual metrics and their assumptions.
    public func computeEntropyProfile() -> [String: FP] {
        return withLock {
            do {
                return try EntropyMeasures(data: self).computeEntropyProfile()
            } catch {
                return [:]
            }
        }
    }

    /// Computes a set of normalized entropy and simple complexity measures.
    ///
    /// Included metrics:
    /// - normalized_entropy: Shannon entropy (nats) divided by log(|support|)
    /// - statistical_complexity: H * (1 - H)
    /// - multiscale_entropy: Average SampEn across scales (numeric-only)
    ///
    /// - Returns: A dictionary of measures, or an empty dictionary on failure.
    ///
    /// Interpretation hints:
    /// - normalized_entropy is in [0, 1] for discrete data; 0 = fully concentrated, 1 = uniform over observed support.
    /// - statistical_complexity peaks at 0.25 when defined as H*(1-H) with H in [0,1]; it is low for both very ordered and very random extremes.
    /// - multiscale_entropy summarizes irregularity across coarse-grained scales; sensitive to scale choices and data length.
    public func complexityMeasures() -> [String: FP] {
        return withLock {
            do {
                return try EntropyMeasures(data: self).complexityMeasures()
            } catch {
                return [:]
            }
        }
    }
}

