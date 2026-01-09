import Testing
import SwiftyBoost
@testable import SwiftyStats

/// Sign test regression and public API coverage tests.
struct SignTestTests {
    private func makeExamine(_ values: [Double]) -> SSExamine<Double, Double> {
        let ex = SSExamine<Double, Double>()
        _ = ex.append(contentOf: values)
        return ex
    }
    
    @Test("Sign test exact p-value matches two-sided binomial tail")
    func signTestExactMatchesBinomialTail() throws {
        let set1 = makeExamine([1, 2, 3, 4])
        let set2 = makeExamine([2, 1, 5, 6])
        
        let resultOptional = try NonParametric<Double>.signTest(set1: set1, set2: set2)
        let result = try #require(resultOptional)
        
        #expect(result.nPosDiff == 3)
        #expect(result.nNegDiff == 1)
        let expectedExact: Double = 0.625
        let pExact = try #require(result.pValueExact)
        #expect(abs(pExact - expectedExact) < 1e-12)
        #expect(result.total == 4)
    }
    
    @Test("Sign test returns finite results when all observations tie")
    func signTestHandlesAllTies() throws {
        let set1 = makeExamine([5, 5, 5])
        let set2 = makeExamine([5, 5, 5])
        
        let resultOptional = try NonParametric<Double>.signTest(set1: set1, set2: set2)
        let result = try #require(resultOptional)
        
        #expect(result.nPosDiff == 0)
        #expect(result.nNegDiff == 0)
        #expect(result.nTies == 3)
        #expect(result.total == 0)
        #expect(result.pValueExact == 1)
        #expect(result.pValueApprox == nil)
        #expect(result.zStat == nil)
    }
    
    @Test("Sign test normal approximation uses two-sided tail")
    func signTestNormalApproximationMatchesClosedForm() throws {
        let nPos = 700
        let nNeg = 301
        let sampleSize = nPos + nNeg
        
        var set1: [Double] = Array(repeating: 0, count: nPos)
        set1.append(contentsOf: Array(repeating: 1, count: nNeg))
        var set2: [Double] = Array(repeating: 1, count: nPos)
        set2.append(contentsOf: Array(repeating: 0, count: nNeg))
        
        let ex1 = makeExamine(set1)
        let ex2 = makeExamine(set2)
        
        let resultOptional = try NonParametric<Double>.signTest(set1: ex1, set2: ex2)
        let result = try #require(resultOptional)
        
        #expect(result.nPosDiff == nPos)
        #expect(result.nNegDiff == nNeg)
        let expectedExactOptional = try Parametric<Double>.binomialTest(successes: nPos, numberOfTrials: sampleSize, p0: 0.5, alternative: .twoSided)
        let expectedExact = try #require(expectedExactOptional)
        let pExact = try #require(result.pValueExact)
        #expect(abs(pExact - expectedExact) < 1e-12)
        
        let diff = Double(nPos) - 0.5 * Double(sampleSize)
        let z = (diff - 0.5) / (Double(sampleSize) * 0.25).squareRoot()
        let normal = try SwiftyBoost.Distribution.Normal<Double>()
        let expectedApprox = min(1.0, 2.0 * (1.0 - (try normal.cdf(abs(z)))))
        
        let pApprox = try #require(result.pValueApprox)
        #expect(abs(pApprox - expectedApprox) < 1e-12)
        let zStat = try #require(result.zStat)
        #expect(zStat > 0)
    }
}
