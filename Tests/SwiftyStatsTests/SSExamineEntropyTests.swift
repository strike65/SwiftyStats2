/// Tests for entropy-related APIs in SSExamine.

import Foundation
import Testing
@testable import SwiftyStats

@Suite("Entropy and Information-Theory Measures on SSExamine")
struct SSExamineEntropyTests {

    // Helper: compute discrete Shannon entropy in base 2 for a simple frequency map
    private func shannon2(_ counts: [Int]) -> Double {
        let n = Double(counts.reduce(0, +))
        return counts.reduce(0.0) { acc, c in
            guard c > 0 else { return acc }
            let p = Double(c) / n
            return acc - p * log2(p)
        }
    }

    @Test("Shannon entropy: discrete vs expected on simple strings")
    func shannonDiscreteOnStrings() throws {
        let ex = SSExamine<String, Double>()
        #expect(try ex.append(text: "AAB", characterSet: .letters))
        // Items are ["A","A","B"]
        let h = ex.shannonEntropy(type: .discrete, base: 2)
        let expected = shannon2([2, 1]) // P(A)=2/3, P(B)=1/3
        #expect(abs(Double(h) - expected) < 1e-12)
    }

    @Test("Shannon entropy: frequency on integers matches discrete on same distribution")
    func shannonFrequencyMatchesDiscrete() {
        let ex = SSExamine<Int, Double>()
        #expect(ex.append(contentOf: [1, 1, 2, 3, 3, 3])) // counts: 1->2, 2->1, 3->3
        let hDiscrete = ex.shannonEntropy(type: .discrete, base: 2)
        let hFrequency = ex.shannonEntropy(type: .frequency, base: 2)
        #expect(abs(Double(hDiscrete - hFrequency)) < 1e-12)
    }

    @Test("Shannon entropy: histogram and Miller-Madow correction")
    func shannonHistogramAndMillerMadow() {
        let ex = SSExamine<Double, Double>()
        #expect(ex.append(contentOf: [0.0, 0.1, 0.2, 1.0, 1.1, 1.2])) // two clusters
        let hHist = ex.shannonEntropy(type: .histogram, base: 2, bins: .count(2))
        let hMM = ex.shannonEntropy(type: .millerMadow, base: 2, bins: .count(2))
        // Miller-Madow should be >= histogram estimate
        #expect(Double(hMM) >= Double(hHist) - 1e-12)
    }

    @Test("Shannon entropy: KDE returns finite value")
    func shannonKDEFinite() {
        let ex = SSExamine<Double, Double>()
        // Slightly noisy Gaussian-like sample
        let values: [Double] = [-1.2, -0.7, -0.1, 0.0, 0.3, 0.8, 1.1, 1.3, -0.5, 0.2, 0.6, -0.9, 1.0]
        #expect(ex.append(contentOf: values))
        let h = ex.shannonEntropy(type: .kde(bandwidth: nil), base: 2)
        #expect(h.isFinite)
    }

    @Test("Renyi entropy on simple distribution matches closed-form")
    func renyiEntropyMatchesClosedForm() {
        let ex = SSExamine<String, Double>()
        // Distribution: A:2, B:1 => p = [2/3, 1/3]
        #expect(ex.append(contentOf: ["A", "A", "B"]))
        let alpha: Double = 0.5
        let r = ex.renyiEntropy(alpha: alpha, base: 2)
        let p = [2.0/3.0, 1.0/3.0]
        let sum = pow(p[0], alpha) + pow(p[1], alpha)
        let expected = log(sum) / ((1 - alpha) * log(2.0))
        #expect(abs(Double(r) - expected) < 1e-12)
    }

    @Test("Hartley, collision (alpha=2), and min-entropy via EntropyMeasures")
    func hartleyCollisionMinEntropy() throws {
        let ex = SSExamine<Int, Double>()
        #expect(ex.append(contentOf: [1, 1, 2, 3, 3, 3])) // unique = 3
        let measures = try SSExamine<Int, Double>.EntropyMeasures(data: ex)

        let hartley = measures.hartleyEntropy(base: 2)
        #expect(abs(Double(hartley) - log2(3.0)) < 1e-12)

        let collision = try measures.collisionEntropy(base: 2)
        #expect(collision.isFinite)

        let minH = measures.minEntropy(base: 2)
        // max probability is 3/6 = 0.5 => min-entropy = -log2(0.5) = 1
        #expect(abs(Double(minH) - 1.0) < 1e-12)
    }

    @Test("Tsallis entropy recovers Shannon limit as q->1 for simple distribution")
    func tsallisEntropyLimit() throws {
        let ex = SSExamine<String, Double>()
        #expect(ex.append(contentOf: ["A", "A", "B"]))
        let measures = try SSExamine<String, Double>.EntropyMeasures(data: ex)
        // Use a delta that triggers the stabilized branch in the implementation
        let tsallisNear1 = try measures.tsallisEntropy(q: 1.0 + 1e-12)
        let shannonNats = ex.shannonEntropy(type: .discrete, base: Double(M_E))
        // Allow a small tolerance due to series expansion approximation
        #expect(abs(Double(tsallisNear1 - shannonNats)) < 1e-5)
    }

    @Test("Approximate entropy (ApEn) and Sample entropy (SampEn) basic sanity")
    func apEnAndSampEn() throws {
        let ex = SSExamine<Double, Double>()
        // Periodic low-complexity series
        var series: [Double] = []
        for _ in 0..<50 { series.append(contentsOf: [0.0, 1.0]) }
        #expect(ex.append(contentOf: series))

        let measures = try SSExamine<Double, Double>.EntropyMeasures(data: ex)
        let apen = try measures.approximateEntropy(m: 2, r: 0.2, normalized: true)
        let sampen = try measures.sampleEntropy(m: 2, r: 0.2, normalized: true)

        #expect(apen.isFinite && sampen.isFinite)
        // For a regular pattern, entropies should be relatively small
        #expect(Double(apen) >= 0)
        #expect(Double(sampen) >= 0)
    }

    @Test("Permutation entropy: strictly increasing sequence yields zero (normalized)")
    func permutationEntropyIncreasingIsZero() throws {
        let ex = SSExamine<Int, Double>()
        #expect(ex.append(contentOf: Array(0..<20)))
        let measures = try SSExamine<Int, Double>.EntropyMeasures(data: ex)
        let pe = try measures.permutationEntropy(order: 3, delay: 1, base: 2, normalized: true)
        #expect(abs(Double(pe)) < 1e-12)
    }

    @Test("KL divergence and cross entropy on simple discrete strings")
    func klAndCrossEntropy() throws {
        // P: A:2, B:1 ; Q: A:1, B:2
        let p = ["A", "A", "B"]
        let q = ["A", "B", "B"]
        let kl = try SSExamine<String, Double>.EntropyMeasures.kullbackLeiblerDivergence(p: p, q: q, base: 2)
        let ce = try SSExamine<String, Double>.EntropyMeasures.crossEntropy(p: p, q: q, base: 2)
        #expect(kl.isFinite)
        #expect(ce.isFinite)
        // Cross-entropy should be >= Shannon entropy of P
        let exP = SSExamine<String, Double>()
        #expect(exP.append(contentOf: p))
        let hP = exP.shannonEntropy(type: .discrete, base: 2)
        #expect(Double(ce) + 1e-12 >= Double(hP))
    }

    @Test("Conditional entropy on subset")
    func conditionalEntropySubset() throws {
        let ex = SSExamine<Int, Double>()
        #expect(ex.append(contentOf: [0, 1, 2, 3, 4, 5, 6, 7])) // uniform 0..7
        let measures = try SSExamine<Int, Double>.EntropyMeasures(data: ex)
        // Condition: even numbers only -> 0,2,4,6 => uniform of 4 items => H = log2(4) = 2
        let h = try measures.conditionalEntropy(condition: { $0 % 2 == 0 }, base: 2)
        #expect(abs(Double(h) - 2.0) < 1e-12)
    }

    @Test("Mutual information: identity transforms yield MI == H(X)")
    func mutualInformationIdentity() throws {
        let ex = SSExamine<Int, Double>()
        #expect(ex.append(contentOf: [1, 1, 2, 3, 3, 3]))
        let measures = try SSExamine<Int, Double>.EntropyMeasures(data: ex)
        let mi = try measures.mutualInformation(transform1: { $0 }, transform2: { $0 }, base: 2)
        let h = ex.shannonEntropy(type: .discrete, base: 2)
        #expect(abs(Double(mi - h)) < 1e-12)
    }

    @Test("Differential entropy estimators (KL and Vasicek) return finite values")
    func differentialEstimatorsFinite() throws {
        let ex = SSExamine<Double, Double>()
        // Smooth increasing data to avoid zero spacings
        let values = stride(from: -2.0, through: 2.0, by: 0.1).map { $0 }
        #expect(ex.append(contentOf: values))

        let measures = try SSExamine<Double, Double>.EntropyMeasures(data: ex)
        let hKL = try measures.kozachenkoLeonenkoEntropy(k: 3, base: 2)
        let hV = try measures.vasicekEntropy(m: nil, base: 2)
        #expect(hKL.isFinite)
        #expect(hV.isFinite)
    }

    @Test("Entropy profile contains expected keys and finite values")
    func entropyProfileKeysAndValues() throws {
        let ex = SSExamine<String, Double>()
        #expect(ex.append(contentOf: ["A", "A", "B", "B", "C", "C", "C"]))
        let measures = try SSExamine<String, Double>.EntropyMeasures(data: ex)
        let profile = try measures.computeEntropyProfile()

        let expectedKeys: Set<String> = [
            "shannon_discrete",
            "shannon_frequency",
            "hartley",
            "renyi_0.5",
            "collision",
            "min_entropy",
            "tsallis_0.5",
            "tsallis_2"
        ]
        for k in expectedKeys {
            #expect(profile.keys.contains(k))
            if let v = profile[k] {
                #expect(v.isFinite)
            }
        }
    }

    @Test("Complexity measures on numeric strings (SSElement == String) are finite")
    func complexityMeasuresOnNumericStrings() throws {
        // Use SSElement == String to exercise numeric conversion and multiscale path
        let ex = SSExamine<String, Double>()
        // Build 40 numeric strings convertible to Double; locale-robust parsing is handled by RealConverter
        let numericStrings = (0..<40).map { i in
            String(format: "%.3f", Double(i % 10) / 10.0)
        }
        #expect(ex.append(contentOf: numericStrings))

        let measures = try SSExamine<String, Double>.EntropyMeasures(data: ex)
        let cm = try measures.complexityMeasures()
        let expected = ["normalized_entropy", "statistical_complexity"]
        for key in expected {
            #expect(cm.keys.contains(key))
            if let v = cm[key] {
                #expect(v.isFinite)
                #expect(Double(v) >= 0.0)
            }
        }
        // multiscale_entropy may be present depending on coarse-grained length; if present, assert finite
        if let mse = cm["multiscale_entropy"] {
            #expect(mse.isFinite)
            #expect(Double(mse) >= 0.0)
        }
    }
}
