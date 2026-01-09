//
//  Created by VT on 22.11.25.
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

/// Internal entropy and information measure implementations for SSExamine.

import SwiftyStatsPrelude
import SwiftyBoost

// MARK: Entropy

/// Internal entropy estimators backing the public APIs, sharing a common container snapshot.
///
/// The routines translate discrete counts into empirical probabilities `p_i = f_i / n` and
/// continuous data into density surrogates (histograms, KDEs, kNN spacing) so that plug-in
/// estimators approximate integrals like `∫ -p log p` without leaking the lock. Bias and
/// variance controls (bin widths, bandwidths, k) are centralized here to keep the public
/// surface concise.
extension SSExamine {
    
    /// A helper class that implements a suite of entropy and information-theory
    /// measures over the data contained in an `SSExamine` instance.
    ///
    /// This type converts elements to numeric values (FP) when a continuous
    /// measure is required, and otherwise operates on the original discrete
    /// values of type `SSElement`.
    internal class EntropyMeasures {
        /// Original items as stored in the dataset, in original order.
        private let items: [SSElement]
        /// The original SSExamine dataset.
        private let orgData: SSExamine<SSElement, FP>
        /// Numeric representation of items for continuous measures.
        private let numericData: [FP]
        /// Number of observations.
        private let n: Int
        /// Sorted items (by SSElement conformance to Comparable).
        private let sortedData: [SSElement]
        /// Sorted numeric data.
        private let sortedNumeric: [FP]
        /// Optional locale influencing numeric conversion.
        private let locale: Locale?
        /// Heuristic: are the data truly convertible to a metric scale?
        private let isTrulyNumeric: Bool
        
        /// Creates an entropy measure helper for the given dataset.
        ///
        /// - Parameters:
        ///   - data: The dataset to analyze.
        ///   - locale: Optional locale used in numeric conversion of elements.
        /// - Throws: `EntropyError.emptyDataset` if data is empty or
        ///   `EntropyError.conversionFailed` if numeric conversion fails.
        init(data: SSExamine<SSElement, FP>, locale: Locale? = nil) throws {
            guard !data.isEmpty else {
                throw EntropyError<SSElement>.emptyDataset
            }
            self.orgData = data
            self.items = data.itemsAsArray(sorted: .original)
            self.n = data.count
            self.locale = locale
            self.sortedData = data.itemsAsArray(sorted: .none).sorted()
            
            // Convert to numeric values; if strings or similar cannot be interpreted
            // as numbers, the implementation falls back to stable category codes
            // (0,1,2,…). In that case, continuous measures are “not meaningful”.
            var converted = [FP]()
            converted.reserveCapacity(items.count)
            var needsCategoricalFallback = false
            var trulyNumeric = true
            
            conversionLoop: for element in items {
                do {
                    let numeric: FP = try RealConverter.to(element, locale: locale)
                    if !numeric.isFinite { throw RealConversionError.nonFinite("\(numeric)") }
                    converted.append(numeric)
                }
                catch let error as RealConversionError {
                    switch error {
                    case .invalidString, .unsupportedType:
                        needsCategoricalFallback = true
                        trulyNumeric = false
                        break conversionLoop
                    default:
                        throw EntropyError.conversionFailed(element: element)
                    }
                }
                catch _ {
                    throw EntropyError.conversionFailed(element: element)
                }
            }
            
            if needsCategoricalFallback {
                converted.removeAll(keepingCapacity: true)
                var codebook = [SSElement: FP]()
                
                for element in items {
                    if let code = codebook[element] {
                        converted.append(code)
                    } else {
                        let code = FP(codebook.count)
                        codebook[element] = code
                        converted.append(code)
                    }
                }
            }
            
            self.numericData = converted
            self.sortedNumeric = converted.sorted()
            self.isTrulyNumeric = trulyNumeric
        }
        
        /// Computes Shannon entropy using various methods.
        ///
        /// - Parameters:
        ///   - method: The estimation method (histogram, frequency, KDE with bandwidth, Miller–Madow, discrete).
        ///   - bins: Binning strategy used by histogram-based methods.
        ///   - base: The logarithm base for the entropy. Defaults to 2 (bits).
        /// - Returns: The estimated Shannon entropy in the given base.
        /// - Throws: `EntropyError` for invalid parameters or numerical issues.
        func shannonEntropy(
            method: ShannonMethod<FP> = .histogram,
            bins: BinningStrategy = .auto,
            base: FP = 2
        ) throws -> FP {
            
            switch method {
            case .histogram:
                return try shannonEntropyHistogram(bins: bins, base: base)
                
            case .frequency:
                return shannonEntropyFrequency(base: base)
                
            case .kde(let bandwidth):
                guard isTrulyNumeric else {
                    // Not numerically meaningful on categorical fallback
                    throw EntropyError<SSElement>.invalidParameter("KDE requires truly metric data")
                }
                return try shannonEntropyKDE(bandwidth: bandwidth, base: base)
                
            case .millerMadow:
                return try shannonEntropyMillerMadow(bins: bins, base: base)
                
            case .discrete:
                return shannonEntropyDiscrete(base: base)
            }
        }
        
        /// Shannon entropy for discrete values (operates directly on `SSElement`).
        ///
        /// - Parameter base: Logarithm base. Defaults to 2.
        /// - Returns: Discrete Shannon entropy in the given base.
        private func shannonEntropyDiscrete(base: FP) -> FP {
            let frequencies = Dictionary(items.map { ($0, 1) }, uniquingKeysWith: +)
            let total = FP(n)
            
            return frequencies.values.reduce(FP.zero) { (acc: FP, count: Int) in
                let p = FP(count) / total
                guard p > 0 else { return acc }
                return acc - p * FP.log(p) / FP.log(base)
            }
        }
        
        /// Shannon entropy computed from frequency counts of numeric values.
        ///
        /// Values within a small relative tolerance are grouped.
        private func shannonEntropyFrequency(base: FP) -> FP {
            // Relative tolerance against scale issues
            let relTol = FP(1e-12)
            let absTol = FP.ulpOfOne * 100
            var frequencies = [FP: Int]()
            
            for value in numericData {
                let tol = max(absTol, relTol * value.magnitude)
                let key = frequencies.keys.first { ($0 - value).magnitude <= tol } ?? value
                frequencies[key, default: 0] += 1
            }
            
            let total = FP(n)
            return frequencies.values.reduce(FP.zero) { acc, count in
                let p = FP(count) / total
                guard p > 0 else { return acc }
                return acc - p * FP.log(p) / FP.log(base)
            }
        }
        
        /// Shannon entropy using a histogram estimator.
        ///
        /// - Parameters:
        ///   - bins: Binning strategy used to build the histogram.
        ///   - base: Logarithm base. Defaults to 2.
        /// - Returns: Shannon entropy in the given base.
        /// - Throws: `EntropyError.numericalError` if min/max cannot be determined.
        private func shannonEntropyHistogram(bins: BinningStrategy, base: FP) throws -> FP {
            let binEdges = try computeBinEdges(strategy: bins)
            let binCounts = histogram(data: numericData, binEdges: binEdges)
            
            let total = FP(n)
            return binCounts.reduce(FP.zero) { (acc: FP, count: Int) in
                guard count > 0 else { return acc }
                let p = FP(count) / total
                return acc - p * FP.log(p) / FP.log(base)
            }
        }
        
        /// Shannon entropy via kernel density estimation (Gaussian kernel).
        ///
        /// - Parameters:
        ///   - bandwidth: Optional bandwidth. If `nil`, use Scott/Silverman 1D rule (~1.06·σ·n^(−1/5)).
        ///   - base: Logarithm base. Defaults to 2.
        /// - Returns: Differential Shannon entropy estimate (continuous).
        /// - Throws: `EntropyError.numericalError` if min/max cannot be determined.
        private func shannonEntropyKDE(bandwidth: FP?, base: FP) throws -> FP {
            let h = max(bandwidth ?? kdeBandwidth(), FP.leastNonzeroMagnitude)
            
            guard let minVal = numericData.min(),
                  let maxVal = numericData.max() else {
                throw EntropyError<SSElement>.numericalError("No min/max values")
            }
            
            // Couple integration range and step size to h
            let evalMin = minVal - 3 * h
            let evalMax = maxVal + 3 * h
            let span = max(evalMax - evalMin, h)
            let dx = max(h / 2, span / FP(512))  // at least ~h/2, at least 512 points across the span
            let nPoints = max(8, Int((span / dx).rounded(.up)))
            
            var entropy = FP.zero
            let logBase = FP.log(base)
            let eps = FP.leastNonzeroMagnitude
            
            for i in 0..<nPoints {
                let x = evalMin + FP(i) * dx
                let density = max(kernelDensity(at: x, bandwidth: h), eps)
                entropy -= density * FP.log(density) * dx / logBase
            }
            
            return entropy
        }
        
        /// Miller–Madow corrected Shannon entropy (bias-corrected histogram).
        ///
        /// Note: This is a bias correction on the binned distribution, not on true symbols.
        private func shannonEntropyMillerMadow(bins: BinningStrategy, base: FP) throws -> FP {
            let baseEntropy = try shannonEntropyHistogram(bins: bins, base: base)
            let nonZeroBins = try countNonZeroBins(bins: bins)
            // Miller–Madow correction (binned)
            let correction = FP(nonZeroBins - 1) / (2 * FP(n) * FP.log(base))
            return baseEntropy + correction
        }
        
        /// Computes Rényi entropy of order `alpha` using discrete frequencies of `SSElement`.
        func renyiEntropy(alpha: FP, base: FP = 2) throws -> FP {
            guard alpha >= 0 && alpha != 1 else {
                throw EntropyError<SSElement>.invalidParameter("Alpha must be >= 0 and != 1")
            }
            
            if alpha == 0 {
                return hartleyEntropy(base: base)
            }
            
            if alpha.isInfinite {
                return minEntropy(base: base)
            }
            
            let frequencies = Dictionary(items.map { ($0, 1) }, uniquingKeysWith: +)
            let total = FP(n)
            
            let sum = frequencies.values.reduce(FP.zero) { (acc: FP, count: Int) in
                let p = FP(count) / total
                return acc + FP.pow(p, alpha)
            }
            
            return FP.log(sum) / ((1 - alpha) * FP.log(base))
        }
        
        /// Hartley entropy (Rényi α = 0).
        func hartleyEntropy(base: FP = 2) -> FP {
            let uniqueValues = Set(items).count
            return FP.log(FP(uniqueValues)) / FP.log(base)
        }
        
        /// Collision entropy (Rényi α = 2).
        func collisionEntropy(base: FP = 2) throws -> FP {
            return try renyiEntropy(alpha: 2, base: base)
        }
        
        /// Min-entropy (Rényi α → ∞).
        func minEntropy(base: FP = 2) -> FP {
            let frequencies = Dictionary(items.map { ($0, 1) }, uniquingKeysWith: +)
            guard let maxFreq = frequencies.values.max() else { return FP.zero }
            let maxProb = FP(maxFreq) / FP(n)
            return -FP.log(maxProb) / FP.log(base)
        }
        
        /// Tsallis entropy of order `q` for discrete `SSElement` frequencies.
        ///
        /// Note: Output in nats (base e).
        func tsallisEntropy(q: FP) throws -> FP {
            guard q > 0 else {
                throw EntropyError<SSElement>.invalidParameter("q must be > 0")
            }
            
            let delta = q - 1
            
            if abs(delta) < FP.ulpOfOne * 100 {
                // Limit for q→1 is Shannon entropy (nats)
                return shannonEntropyDiscrete(base: FP(M_E))
            }
            
            let frequencies = Dictionary(items.map { ($0, 1) }, uniquingKeysWith: +)
            let total = FP(n)

            if abs(delta) < FP(1e-6) {
                var numerator = FP.zero
                for count in frequencies.values {
                    let p = FP(count) / total
                    guard p > 0 else { continue }
                    let logp = FP.log(p)
                    let x = delta * logp
                    // expm1-like stabilization (split up to ease the type checker)
                    let ax = abs(x)
                    let expm1Delta: FP
                    if ax < FP(1e-6) {
                        let x2 = x * x
                        let x3 = x2 * x
                        expm1Delta = x + x2 * FP(0.5) + x3 * FP(1.0 / 6.0)
                    } else {
                        expm1Delta = FP.exp(x) - FP.one
                    }
                    numerator += p * expm1Delta
                }
                return -numerator / delta
            }
            
            let sum = frequencies.values.reduce(FP.zero) { (acc: FP, count: Int) in
                let p = FP(count) / total
                return acc + FP.pow(p, q)
            }
            return (1 - sum) / (q - 1)
        }
        
        /// Approximate entropy (ApEn) for numeric time series.
        func approximateEntropy(m: Int, r: FP, normalized: Bool = true) throws -> FP {
            guard isTrulyNumeric else {
                throw EntropyError<SSElement>.invalidParameter("ApEn requires truly metric data")
            }
            guard m > 0 && m < n else {
                throw EntropyError<SSElement>.invalidParameter("m must be between 1 and n-1")
            }
            
            let tolerance = normalized ? r * standardDeviation() : r
            let phi_m = try computePhiApEn(m: m, r: tolerance)
            let phi_m1 = try computePhiApEn(m: m + 1, r: tolerance)
            return phi_m - phi_m1
        }
        
        /// Helper for approximate entropy: average log match rate of length-`m` templates.
        ///
        /// - Parameters:
        ///   - m: Template length.
        ///   - r: Chebyshev tolerance for considering two windows similar.
        /// - Returns: φ_m = (1/N) Σ_i log(C_i / N) where C_i counts r-close matches to template i.
        /// - Throws: Propagates errors from the caller when preconditions fail.
        private func computePhiApEn(m: Int, r: FP) throws -> FP {
            let N = n - m + 1
            var phi = FP.zero
            
            for i in 0..<N {
                let template = Array(numericData[i..<(i+m)])
                var matches = 0
                for j in 0..<N {
                    let window = Array(numericData[j..<(j+m)])
                    if chebyshevDistance(template, window) <= r {
                        matches += 1
                    }
                }
                phi += FP.log(FP(matches) / FP(N))
            }
            return phi / FP(N)
        }
        
        /// Sample entropy (SampEn) for numeric time series.
        func sampleEntropy(m: Int, r: FP, normalized: Bool = true) throws -> FP {
            guard isTrulyNumeric else {
                throw EntropyError<SSElement>.invalidParameter("SampEn requires truly metric data")
            }
            guard m > 0 && m < n-1 else {
                throw EntropyError<SSElement>.invalidParameter("m must be between 1 and n-2")
            }
            
            let tolerance = normalized ? r * standardDeviation() : r
            let B = countTemplateMatches(m: m, r: tolerance)
            let A = countTemplateMatches(m: m + 1, r: tolerance)
            guard A > 0 && B > 0 else { return FP.infinity }
            return -FP.log(FP(A) / FP(B))
        }
        
        /// Counts r-close matches for SampEn between all distinct template pairs of length `m`.
        ///
        /// - Parameters:
        ///   - m: Template length.
        ///   - r: Chebyshev tolerance for similarity.
        /// - Returns: Number of (i, j) pairs with Chebyshev distance ≤ r.
        private func countTemplateMatches(m: Int, r: FP) -> Int {
            let N = n - m + 1
            var totalMatches = 0
            for i in 0..<(N-1) {
                let template = Array(numericData[i..<(i+m)])
                for j in (i+1)..<N {
                    let window = Array(numericData[j..<(j+m)])
                    if chebyshevDistance(template, window) <= r {
                        totalMatches += 1
                    }
                }
            }
            return totalMatches
        }
        
        /// Permutation entropy based on ordinal patterns.
        func permutationEntropy(
            order: Int,
            delay: Int = 1,
            base: FP = 2,
            normalized: Bool = true
        ) throws -> FP {
            guard order > 1 && order <= 7 else {
                throw EntropyError<SSElement>.invalidParameter("Order should be between 2 and 7")
            }
            guard delay > 0 else {
                throw EntropyError<SSElement>.invalidParameter("Delay must be positive")
            }
            
            let maxIndex = n - (order - 1) * delay
            guard maxIndex > 0 else {
                throw EntropyError<SSElement>.insufficientData
            }
            
            var patternCounts = [String: Int]()
            for i in 0..<maxIndex {
                var subseries = [SSElement]()
                for j in 0..<order {
                    subseries.append(items[i + j * delay])
                }
                let pattern = computeOrdinalPattern(subseries)
                patternCounts[pattern, default: 0] += 1
            }
            
            let total = FP(maxIndex)
            var entropy = FP.zero
            let logBase = FP.log(base)
            for count in patternCounts.values {
                let p = FP(count) / total
                if p > 0 {
                    entropy -= p * FP.log(p) / logBase
                }
            }
            
            if normalized {
                // MaxEnt = log(order!)
                let maxEnt: FP
                do {
                    let fact = try SwiftyBoost.SpecialFunctions.factorial(UInt32(order))
                    maxEnt = FP.log(FP(fact)) / logBase
                } catch {
                    // Fallback: numerical approximation via Stirling
                    let o = FP(order)
                    let stirling = o * FP.log(o) - o + FP.log(2 * FP.pi * o) / 2
                    maxEnt = stirling / logBase
                }
                return entropy / maxEnt
            }
            return entropy
        }
        
        /// Computes the ordinal pattern string for a slice of values.
        ///
        /// Stable tie-break rule: first by value, then by index.
        private func computeOrdinalPattern(_ values: [SSElement]) -> String {
            let indices = Array(0..<values.count)
            let sortedIndices = indices.sorted {
                let lhs = values[$0]
                let rhs = values[$1]
                if lhs == rhs { return $0 < $1 }
                return lhs < rhs
            }
            return sortedIndices.map { String($0) }.joined()
        }
        
        /// Kozachenko–Leonenko differential entropy estimator (1D).
        ///
        /// Norm: 1D Euclidean (|·|). For 1D there is an additive term ln(2).
        func kozachenkoLeonenkoEntropy(k: Int = 3, base: FP = 2) throws -> FP {
            guard isTrulyNumeric else {
                throw EntropyError<SSElement>.invalidParameter("KL entropy requires truly metric data")
            }
            guard k > 0 && k < n else {
                throw EntropyError<SSElement>.invalidParameter("k must be between 1 and n-1")
            }
            
            var sum = FP.zero
            for i in 0..<n {
                var distances = [FP]()
                distances.reserveCapacity(n-1)
                for j in 0..<n where i != j {
                    distances.append(abs(numericData[i] - numericData[j]))
                }
                distances.sort()
                let kthDistance: FP = max(distances[k-1], FP.leastNonzeroMagnitude)
                sum += FP.log(kthDistance)
            }
            
            // KL estimator (1D) with digamma and ln(2)
            let psi_k: FP
            let psi_n: FP
            do {
                psi_k = try SwiftyBoost.SpecialFunctions.digamma(FP(k))
                psi_n = try SwiftyBoost.SpecialFunctions.digamma(FP(n))
            } catch {
                throw EntropyError<SSElement>.numericalError("digamma failed")
            }
            let h_nats = psi_n - psi_k + FP.log(2) + sum / FP(n)
            return h_nats / FP.log(base)
        }
        
        /// Vasicek differential entropy estimator (1D) based on order statistics spacing.
        func vasicekEntropy(m: Int? = nil, base: FP = 2) throws -> FP {
            guard isTrulyNumeric else {
                throw EntropyError<SSElement>.invalidParameter("Vasicek entropy requires truly metric data")
            }
            let window = m ?? Int(FP.sqrt(FP(n)))
            guard window > 0 && 2 * window < n else {
                throw EntropyError<SSElement>.insufficientData
            }
            
            var sum = FP.zero
            for i in window..<(n-window) {
                let upperValue = sortedNumeric[i + window]
                let lowerValue = sortedNumeric[i - window]
                let spacing = upperValue - lowerValue
                if spacing > 0 {
                    sum += FP.log(spacing * FP(n) / (2 * FP(window)))
                }
            }
            return sum / FP(n - 2 * window) / FP.log(base)
        }
        
        /// Computes the Kullback–Leibler divergence D_KL(P || Q) for discrete `SSElement` sequences.
        static func kullbackLeiblerDivergence(
            p: [SSElement],
            q: [SSElement],
            base: FP = 2,
            locale: Locale? = nil
        ) throws -> FP {
            let pFreq = Dictionary(p.map { ($0, 1) }, uniquingKeysWith: +)
            let qFreq = Dictionary(q.map { ($0, 1) }, uniquingKeysWith: +)
            let allKeys = Set(pFreq.keys).union(qFreq.keys)
            let totalP = FP(p.count)
            let totalQ = FP(q.count)
            var kl = FP.zero
            let logBase = FP.log(base)
            for key in allKeys {
                let pProb = FP(pFreq[key] ?? 0) / totalP
                let qProb = FP(qFreq[key] ?? 0) / totalQ
                if pProb > 0 {
                    if qProb > 0 {
                        kl += pProb * FP.log(pProb / qProb) / logBase
                    } else {
                        return FP.infinity
                    }
                }
            }
            return kl
        }
        
       
        /// Computes cross entropy H(P, Q) for discrete `SSElement` sequences.
        static func crossEntropy(
            p: [SSElement],
            q: [SSElement],
            base: FP = 2
        ) throws -> FP {
            let pFreq = Dictionary(p.map { ($0, 1) }, uniquingKeysWith: +)
            let qFreq = Dictionary(q.map { ($0, 1) }, uniquingKeysWith: +)
            let totalP = FP(p.count)
            let totalQ = FP(q.count)
            var ce = FP.zero
            let logBase = FP.log(base)
            for (key, pCount) in pFreq {
                let pProb = FP(pCount) / totalP
                let qProb = FP(qFreq[key] ?? 0) / totalQ
                if pProb > 0 {
                    if qProb > 0 {
                        ce -= pProb * FP.log(qProb) / logBase
                    } else {
                        return FP.infinity
                    }
                }
            }
            return ce
        }
        
        /// Computes conditional entropy H(X | condition(X) == true) for a boolean condition on elements.
        func conditionalEntropy(
            condition: (SSElement) -> Bool,
            base: FP = 2
        ) throws -> FP {
            var ss: SSExamine<SSElement, FP>
            do {
                ss = try SSExamine<SSElement, FP>.init(using: items.filter(condition), levelOfMeasurement: .ratio, name: nil, characterSet: nil)
            }
            catch _ {
                ss = SSExamine<SSElement, FP>()
            }
            guard !ss.isEmpty else {
                throw EntropyError<SSElement>.emptySubset
            }
            let subsetMeasures = try EntropyMeasures(data: ss, locale: locale)
            return subsetMeasures.shannonEntropyDiscrete(base: base)
        }
        
        /// Estimates mutual information I(X;Y) between two discrete transforms of the elements.
        func mutualInformation(
            transform1: (SSElement) -> Int,
            transform2: (SSElement) -> Int,
            base: FP = 2
        ) throws -> FP {
            let labels1 = items.map(transform1)
            let labels2 = items.map(transform2)
            // Joint distribution with tuple key (without String overhead)
            var jointCounts = [PairKey: Int]()
            var marginal1 = [Int: Int]()
            var marginal2 = [Int: Int]()
            for i in 0..<n {
                let key = PairKey(a: labels1[i], b: labels2[i])
                jointCounts[key, default: 0] += 1
                marginal1[labels1[i], default: 0] += 1
                marginal2[labels2[i], default: 0] += 1
            }
            var mi = FP.zero
            let total = FP(n)
            let logBase = FP.log(base)
            for (key, count) in jointCounts {
                let pJoint = FP(count) / total
                let p1 = FP(marginal1[key.a] ?? 0) / total
                let p2 = FP(marginal2[key.b] ?? 0) / total
                if pJoint > 0 && p1 > 0 && p2 > 0 {
                    mi += pJoint * FP.log(pJoint / (p1 * p2)) / logBase
                }
            }
            return mi
        }
        
        /// Computes a profile of commonly used entropy and complexity measures.
        func computeEntropyProfile() throws -> [String: FP] {
            var profile = [String: FP]()
            
            // Shannon entropy variants
            profile["shannon_discrete"] = shannonEntropyDiscrete(base: 2)
            profile["shannon_frequency"] = shannonEntropyFrequency(base: 2)
            
            // Rényi entropy family
            profile["hartley"] = hartleyEntropy()
            profile["renyi_0.5"] = try renyiEntropy(alpha: 0.5)
            profile["collision"] = try collisionEntropy()
            profile["min_entropy"] = minEntropy()
            
            // Tsallis entropy
            profile["tsallis_0.5"] = try tsallisEntropy(q: 0.5)
            profile["tsallis_2"] = try tsallisEntropy(q: 2)
            
            // Time-series measures (if sufficient data)
            if isTrulyNumeric && n > 100 {
                profile["approx_entropy"] = try approximateEntropy(m: 2, r: 0.2)
                profile["sample_entropy"] = try sampleEntropy(m: 2, r: 0.2)
                profile["permutation_entropy"] = try permutationEntropy(order: 3)
            } else if n > 100 {
                // Only permutation entropy is independent of metric scale
                profile["permutation_entropy"] = try permutationEntropy(order: 3)
            }
            
            // Continuous estimators
            if isTrulyNumeric && n > 50 {
                profile["kl_entropy"] = try kozachenkoLeonenkoEntropy(k: 3)
                profile["vasicek_entropy"] = try vasicekEntropy()
            }
            
            return profile
        }
        
        /// Computes a set of normalized entropy and simple complexity measures.
        ///
        /// - normalized_entropy: Shannon entropy (nats) / log(|support|)
        /// - statistical_complexity: H * (1 - H)
        /// - multiscale_entropy: Average of SampEn across multiple scales (metric-only)
        func complexityMeasures() throws -> [String: FP] {
            var measures = [String: FP]()
            
            // Normalized Shannon entropy (discrete, base e)
            let uniqueCount = FP(Set(items).count)
            let maxEntropy = uniqueCount > 1 ? FP.log(uniqueCount) : FP.zero
            let shannonNats = shannonEntropyDiscrete(base: FP(M_E))
            measures["normalized_entropy"] = maxEntropy > 0 ? shannonNats / maxEntropy : 0
            
            // Statistical complexity
            let H = measures["normalized_entropy"]!
            let D = 1 - H
            measures["statistical_complexity"] = H * D
            
            // Multiscale entropy (only for truly metric data)
            if isTrulyNumeric && n > 20 {
                let candidateScales = [1, 2, 4, 8]
                let scales = candidateScales.filter { $0 < n / 10 && $0 >= 1 }
                var mse = FP.zero
                var used = 0
                for scale in scales {
                    var coarseGrained = [FP]()
                    coarseGrained.reserveCapacity(n / scale + 1)
                    for i in stride(from: 0, to: n - scale + 1, by: scale) {
                        let segment = numericData[i..<min(i + scale, n)]
                        let mean = segment.reduce(FP.zero, +) / FP(segment.count)
                        coarseGrained.append(mean)
                    }
                    if coarseGrained.count > 10 {
                        let se = sampleEntropyForArray(series: coarseGrained, m: 2, r: 0.15, normalized: true)
                        if se.isFinite {
                            mse += se
                            used += 1
                        }
                    }
                }
                if used > 0 {
                    measures["multiscale_entropy"] = mse / FP(used)
                }
            }
            return measures
        }
        
        /// Returns the sample standard deviation or 0 on failure.
        private func standardDeviation() -> FP {
            if let sd: FP = self.orgData.sampleStandardDeviation, sd.isFinite {
                return sd
            }
            guard n > 1 else { return FP.zero }
            let mean = numericData.reduce(FP.zero, +) / FP(n)
            let variance = numericData.reduce(FP.zero) { acc, x in
                let diff = x - mean
                return acc + diff * diff
            } / FP(n - 1)
            let fallback = FP.sqrt(max(variance, FP.zero))
            return fallback.isFinite ? fallback : FP.zero
        }
        
        /// Chebyshev (L∞) distance between two equal-length vectors.
        ///
        /// - Parameters:
        ///   - a: First vector.
        ///   - b: Second vector.
        /// - Returns: max_i |a_i − b_i|, or 0 when vectors are empty.
        private func chebyshevDistance(_ a: [FP], _ b: [FP]) -> FP {
            return zip(a, b).map { abs($0 - $1) }.max() ?? FP.zero
        }
        
        /// Gaussian KDE at x.
        private func kernelDensity(at x: FP, bandwidth h: FP) -> FP {
            let factor = 1 / (FP(n) * h * FP.sqrt(2 * FP.pi))
            var acc = FP.zero
            for xi in numericData {
                let z = (x - xi) / h
                acc += factor * FP.exp(-0.5 * z * z)
            }
            return acc
        }
        
        /// Histogram Scott width (3.49·σ·n^(−1/3)).
        private func histogramScottBinWidth() -> FP {
            let sigma = standardDeviation()
            return FP(3.49) * sigma * FP.pow(FP(n), -FP.one/FP(3))
        }
        
        /// KDE bandwidth (Silverman/Scott 1D: ~1.06·σ·n^(−1/5)).
        private func kdeBandwidth() -> FP {
            let sigma = standardDeviation()
            return FP(1.06) * sigma * FP.pow(FP(n), -FP.one/FP(5))
        }
        
        /// Computes the number of bins for a given binning strategy.
        private func computeBinCount(strategy: BinningStrategy) throws -> Int {
            switch strategy {
            case .auto:
                return try computeBinCount(strategy: .sturges)
            case .count(let bins):
                return max(1, bins)
            case .sturges:
                return Int(ceil(FP.log2(FP(n)) + 1 ))
            case .sqrt:
                return Int(ceil(FP.sqrt(FP(n))))
            case .rice:
                return Int(ceil(2 * FP.pow(FP(n), 1/3)))
            case .scott, .freedmanDiaconis:
                let h = try computeBinWidth(strategy: strategy)
                guard let min = numericData.min(),
                      let max = numericData.max() else {
                    throw EntropyError<SSElement>.numericalError("No min/max values")
                }
                let range = max - min
                let countF = (range / h).isFinite ? (range / h) : FP.one
                return Int(FP.maximum(1, countF.rounded(.up)))
            }
        }
        
        /// Computes the bin width for a given strategy.
        private func computeBinWidth(strategy: BinningStrategy) throws -> FP {
            switch strategy {
            case .scott:
                // Histogram Scott; fallback if σ=0
                var w = histogramScottBinWidth()
                if !w.isFinite || w <= 0 {
                    if let minVal = numericData.min(), let maxVal = numericData.max() {
                        let range = maxVal - minVal
                        w = max(range * FP(1e-6), FP.leastNonzeroMagnitude)
                    } else {
                        w = FP.one
                    }
                }
                return w
            case .freedmanDiaconis:
                let q1 = percentile(25)
                let q3 = percentile(75)
                var iqr = q3 - q1
                if !iqr.isFinite || iqr <= 0 {
                    // Fallback: estimate directly from the sorted series
                    if sortedNumeric.count >= 4 {
                        let idx1 = sortedNumeric.count / 4
                        let idx3 = 3 * sortedNumeric.count / 4
                        iqr = sortedNumeric[idx3] - sortedNumeric[idx1]
                    }
                }
                var w = 2 * iqr * FP.pow(FP(n), -FP.one/FP(3))
                if !w.isFinite || w <= 0 {
                    if let minVal = numericData.min(), let maxVal = numericData.max() {
                        let range = maxVal - minVal
                        w = max(range * FP(1e-6), FP.leastNonzeroMagnitude)
                    } else {
                        w = FP.one
                    }
                }
                return w
            default:
                let bins = try computeBinCount(strategy: strategy)
                guard let min = numericData.min(),
                      let max = numericData.max() else {
                    throw EntropyError<SSElement>.numericalError("No min/max values")
                }
                let width = (max - min) / FP(bins)
                return FP.maximum(width, FP.leastNonzeroMagnitude)
            }
        }
        
        /// Computes the bin edges for a given binning strategy.
        private func computeBinEdges(strategy: BinningStrategy) throws -> [FP] {
            let bins = try computeBinCount(strategy: strategy)
            guard let minVal = numericData.min(),
                  let maxVal = numericData.max() else {
                throw EntropyError<SSElement>.numericalError("No min/max values")
            }
            let width = max((maxVal - minVal) / FP(bins), FP.leastNonzeroMagnitude)
            // Small buffer so that max is surely included
            let start = minVal
            var edges = [FP]()
            edges.reserveCapacity(bins + 1)
            for i in 0...bins {
                edges.append(start + FP(i) * width)
            }
            // Ensure the last edge ≥ maxVal
            if edges.last! < maxVal {
                edges[edges.count - 1] = maxVal
            }
            return edges
        }
        
        /// Builds a histogram with clamped under-/overflow and inclusive last edge.
        private func histogram(data: [FP], binEdges: [FP]) -> [Int] {
            var counts = [Int](repeating: 0, count: binEdges.count - 1)
            let lastIndex = counts.count - 1
            let lastEdge = binEdges.last!
            let firstEdge = binEdges.first!
            for value in data {
                if value <= firstEdge {
                    counts[0] += 1
                    continue
                }
                if value >= lastEdge {
                    counts[lastIndex] += 1
                    continue
                }
                // find i with e_i <= v < e_{i+1}
                if let i = binEdges.indices.dropLast().last(where: { binEdges[$0] <= value }) {
                    let idx = min(max(0, i), lastIndex)
                    counts[idx] += 1
                }
            }
            return counts
        }
        
        /// Counts the number of non-empty bins for a given strategy.
        private func countNonZeroBins(bins: BinningStrategy) throws -> Int {
            let edges = try computeBinEdges(strategy: bins)
            let counts = histogram(data: numericData, binEdges: edges)
            return counts.filter { $0 > 0 }.count
        }
        
        /// Computes a percentile using the quantile function provided by `orgData`.
        private func percentile(_ p: FP) -> FP {
            if let p = try? self.orgData.quantile(q: p / FP(100), quantileType: .type2) {
                return p
            }
            else {
                return 0
            }
        }
        
        /// Helper: SampEn on a given FP series.
        private func sampleEntropyForArray(series: [FP], m: Int, r: FP, normalized: Bool) -> FP {
            guard series.count > m + 1 else { return FP.infinity }
            let sd: FP = {
                if series.count <= 1 { return FP.zero }
                let mean = series.reduce(FP.zero, +) / FP(series.count)
                let varN1 = series.reduce(FP.zero) { $0 + ( $1 - mean ) * ( $1 - mean ) } / FP(series.count - 1)
                return FP.sqrt(max(varN1, FP.zero))
            }()
            let tol = normalized ? r * sd : r
            let N = series.count - m + 1
            var B = 0
            var A = 0
            // B: pairs in m, A: pairs in m+1
            for i in 0..<(N-1) {
                let templM = Array(series[i..<(i+m)])
                let templM1 = Array(series[i..<(i+m+1)])
                for j in (i+1)..<N {
                    let winM = Array(series[j..<(j+m)])
                    if chebyshevDistance(templM, winM) <= tol { B += 1 }
                    if j < series.count - m {
                        let winM1 = Array(series[j..<(j+m+1)])
                        if chebyshevDistance(templM1, winM1) <= tol { A += 1 }
                    }
                }
            }
            guard A > 0 && B > 0 else { return FP.infinity }
            return -FP.log(FP(A) / FP(B))
        }
        
        /// Internal hashable tuple key for MI.
        private struct PairKey: Hashable {
            let a: Int
            let b: Int
        }
    }
}
