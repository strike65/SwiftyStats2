/// Tests for variance calculations and related helper functions.

import Foundation
import Testing
@testable import SwiftyStats
import SwiftyBoost

struct VarianceTests {
    private func approxEqual(_ lhs: Double, _ rhs: Double, tolerance: Double = 1e-12) -> Bool {
        let scale = max(1.0, max(abs(lhs), abs(rhs)))
        return abs(lhs - rhs) <= tolerance * scale
    }

    private func makeExamine(_ values: [Double], name: String) -> SSExamine<Double, Double> {
        let ex = SSExamine<Double, Double>()
        _ = ex.append(contentOf: values)
        ex.name = name
        return ex
    }
    
    private func leveneStatisticMeanCentered(_ groups: [[Double]]) -> Double {
        var zMeans: [Double] = []
        var counts: [Double] = []
        var withinGroup: Double = 0
        var totalCount: Double = 0
        var totalAbsDeviation: Double = 0
        for values in groups {
            guard !values.isEmpty else { return .nan }
            let center = values.reduce(0, +) / Double(values.count)
            var meanAbsDeviation = 0.0
            var m2 = 0.0
            var count = 0.0
            for value in values {
                let deviation = abs(value - center)
                count += 1
                let delta = deviation - meanAbsDeviation
                meanAbsDeviation += delta / count
                let delta2 = deviation - meanAbsDeviation
                m2 += delta * delta2
            }
            zMeans.append(meanAbsDeviation)
            counts.append(count)
            withinGroup += m2
            totalCount += count
            totalAbsDeviation += meanAbsDeviation * count
        }
        let zMean = totalAbsDeviation / totalCount
        var betweenGroup = 0.0
        var idx = 0
        while idx < zMeans.count {
            let diff = zMeans[idx] - zMean
            betweenGroup += counts[idx] * diff * diff
            idx += 1
        }
        let numerator = (totalCount - Double(zMeans.count)) * betweenGroup
        let denominator = (Double(zMeans.count) - 1.0) * withinGroup
        return numerator / denominator
    }

    @Test("Levene test matches manual F-statistic for mean-centred groups")
    func leveneMatchesManualMeanCenter() throws {
        let groups = [
            makeExamine([6, 7, 8, 9, 5], name: "g1"),
            makeExamine([3, 4, 2, 5, 1], name: "g2"),
            makeExamine([10, 12, 9, 11, 14], name: "g3")
        ]
        let dataFrame = try DataFrame(data: groups)

        let result = try #require(
            try Inferential<Double>.HypothesisTesting.VarianceTests.leveneTest(
                samples: dataFrame,
                alpha: 0.05,
                testType: .mean
            )
        )

        let centers = try groups.map {
            try #require($0.arithmeticMean)
        }
        let absoluteDeviations: [[Double]] = try zip(groups, centers).map { group, center in
            let values = try #require(group.itemsAsNumericArray)
            return values.map { abs($0 - center) }
        }
        let ni = absoluteDeviations.map { Double($0.count) }
        let zi = zip(absoluteDeviations, ni).map { values, count in
            values.reduce(0.0, +) / count
        }
        let overallMean = absoluteDeviations.flatMap { $0 }.reduce(0.0, +) / ni.reduce(0.0, +)
        let numerator = zip(ni, zi).reduce(0.0) { partial, pair in
            let (count, mean) = pair
            return partial + count * pow(mean - overallMean, 2)
        }
        let denominator = zip(absoluteDeviations, zi).reduce(0.0) { partial, pair in
            let (values, mean) = pair
            return partial + values.reduce(0.0) { $0 + pow($1 - mean, 2) }
        }
        let expectedF = ((ni.reduce(0.0, +) - Double(ni.count)) / (Double(ni.count) - 1.0)) * (numerator / denominator)

        #expect(denominator > 0)
        #expect(result.testType == .levene(.mean))
        #expect(result.testStatistic!.isFinite)
        #expect(abs(result.testStatistic! - expectedF) < 1e-10)
    }
    
    @Test("Levene test matches manual F-statistic for large mean-centred groups")
    func leveneMatchesManualLargeDataset() throws {
        let block1: [Double] = [-2, -1, 1, 2]
        let block2: [Double] = [-3, 0, 3]
        let repeats1 = 2_500
        let repeats2 = 3_334
        
        var largeGroup1: [Double] = []
        largeGroup1.reserveCapacity(block1.count * repeats1)
        for _ in 0..<repeats1 {
            largeGroup1.append(contentsOf: block1)
        }
        
        var largeGroup2: [Double] = []
        largeGroup2.reserveCapacity(block2.count * repeats2)
        for _ in 0..<repeats2 {
            largeGroup2.append(contentsOf: block2)
        }
        
        let examineGroups = [
            makeExamine(largeGroup1, name: "large1"),
            makeExamine(largeGroup2, name: "large2")
        ]
        let dataFrame = try DataFrame(data: examineGroups)
        
        let result = try #require(
            try Inferential<Double>.HypothesisTesting.VarianceTests.leveneTest(
                samples: dataFrame,
                alpha: 0.05,
                testType: .mean
            )
        )
        
        let expectedStatistic = leveneStatisticMeanCentered([largeGroup1, largeGroup2])
        #expect(result.testType == .levene(.mean))
        #expect(result.testStatistic!.isFinite)
        #expect(abs(result.testStatistic! - expectedStatistic) < 1e-10)
    }

    @Test("F-ratio test reports directional tails and confidence intervals")
    func fRatioTestDirectionalOutputs() throws {
        let batch1 = SSExamine<Double, Double>(usingArray: [1, 1, 2, 4, 7], name: "b1", characterSet: nil)
        let batch2 = SSExamine<Double, Double>(usingArray: [1, 2, 2, 2, 3], name: "b2", characterSet: nil)
        let alpha = 0.05

        let result = try #require(
            try Inferential<Double>.HypothesisTesting.VarianceTests.fRatioTest(
                batch1: batch1,
                batch2: batch2,
                alpha: alpha
            )
        )
        let resultArray = try #require(
            try Inferential<Double>.HypothesisTesting.VarianceTests.fRatioTest(
                batch1: [1, 1, 2, 4, 7],
                batch2: [1, 2, 2, 2, 3],
                alpha: alpha
            )
        )

        let s1 = try #require(batch1.sampleVariance)
        let s2 = try #require(batch2.sampleVariance)
        let fRatio = s1 / s2
        let df1 = Double(batch1.sampleSize - 1)
        let df2 = Double(batch2.sampleSize - 1)
        let dist = try Distribution.FisherF<Double>(degreesOfFreedom1: df1, degreesOfFreedom2: df2)
        let cdf = try dist.cdf(fRatio)
        let pGreater = 1.0 - cdf
        let pLess = cdf
        let pTwo = cdf > 0.5 ? 2.0 * (1.0 - cdf) : 2.0 * cdf
        let pOne = min(pGreater, pLess)

        #expect(approxEqual(result.FRatio!, fRatio))
        #expect(approxEqual(result.p2Value!, pTwo))
        #expect(approxEqual(result.p1Value!, pOne))
        #expect(result.FRatioGTE1 == (pGreater < alpha))
        #expect(result.FRatioLTE1 == (pLess < alpha))
        #expect(result.FRatioEQ1 == (pTwo >= alpha))
        #expect(result.dfDenominator != nil)
        #expect(result.dfNumerator != nil)
        #expect(result.variance1 != nil)
        #expect(result.variance2 != nil)
        #expect(result.sd1 != nil)
        #expect(result.sd2 != nil)
        #expect(result.mean1 != nil)
        #expect(result.mean2 != nil)
        #expect(!result.description.isEmpty)

        let eqCI = try #require(result.ciRatioEQ1)
        let expectedEqLower = try fRatio / dist.quantile(1.0 - alpha / 2.0)
        let expectedEqUpper = try fRatio / dist.quantile(alpha / 2.0)
        #expect(approxEqual(eqCI.lowerBound, expectedEqLower))
        #expect(approxEqual(eqCI.upperBound, expectedEqUpper))

        let ltCI = try #require(result.ciRatioLTE1)
        let expectedLtUpper = try fRatio / dist.quantile(alpha)
        #expect(ltCI.lowerBound == .zero)
        #expect(approxEqual(ltCI.upperBound, expectedLtUpper))

        let gtCI = try #require(result.ciRatioGTE1)
        let expectedGtLower = try fRatio / dist.quantile(1.0 - alpha)
        #expect(approxEqual(gtCI.lowerBound, expectedGtLower))
        #expect(gtCI.upperBound.isInfinite)
        let postHoc = try #require(result.postHocPower)
        #expect(postHoc >= 0 && postHoc <= 1)

        #expect(resultArray.sampleSize1 != nil)
        #expect(resultArray.sampleSize2 != nil)
        #expect(resultArray.FRatio != nil)
    }

    @Test("Balanced F-ratio sample size search stops at the first sufficient n")
    func balancedSampleSizeSearchStopsEarly() throws {
        let alpha = 0.05
        let theta = 1.0

        let earlyHit = try #require(
            try Inferential<Double>.HypothesisTesting.VarianceTests.nFRatioTestBalanced(
                alpha: alpha,
                targetPower: 0.049,
                theta: theta,
                nMin: 5,
                nMax: 10
            )
        )
        #expect(earlyHit.n == 5)
        #expect(earlyHit.power >= 0.049)

        let noSolution = try Inferential<Double>.HypothesisTesting.VarianceTests.nFRatioTestBalanced(
            alpha: alpha,
            targetPower: 0.2,
            theta: theta,
            nMin: 5,
            nMax: 6
        )
        #expect(noSolution == nil)
    }

    @Test("Bartlett test returns a populated variance-equality result for both overloads")
    func bartlettTestResultAccessors() throws {
        let groups = [
            makeExamine([1, 1, 2, 2, 3], name: "g1"),
            makeExamine([10, 11, 12, 13, 14], name: "g2"),
            makeExamine([0, 0, 0, 1, 10], name: "g3")
        ]
        let dataFrame = try DataFrame(data: groups)

        let res1 = try #require(try Inferential<Double>.HypothesisTesting.VarianceTests.bartlettTest(samples: dataFrame, alpha: 0.05))
        let res2 = try #require(try Inferential<Double>.HypothesisTesting.VarianceTests.bartlettTest(samples: groups, alpha: 0.05))

        for res in [res1, res2] {
            let df = try #require(res.df)
            #expect(df.df1 == 2.0)
            #expect(df.df2.isNaN)
            #expect(res.cv90Pct != nil)
            #expect(res.cv95Pct != nil)
            #expect(res.cv99Pct != nil)
            #expect(res.cvAlpha != nil)
            #expect(res.pValue != nil)
            #expect(res.testStatistic != nil)
            #expect(res.equality != nil)
            #expect(res.testType == .bartlett)
            #expect(!res.description.isEmpty)
        }
    }

    @Test("Chi-squared variance test returns a populated result for both overloads")
    func chiSquaredVarianceTestResultAccessors() throws {
        let data: [Double] = [0.8, 1.2, 1.0, 1.1, 0.9, 1.3]
        let ex = SSExamine<Double, Double>(usingArray: data, name: "x", characterSet: nil)

        let res1 = try #require(try Inferential<Double>.HypothesisTesting.VarianceTests.chiSquaredTest(sample: ex, nominalVariance: 1.0, alpha: 0.05, meanIsKnown: false))
        let res2 = try #require(try Inferential<Double>.HypothesisTesting.VarianceTests.chiSquaredTest(sample: data, nominalVariance: 1.0, alpha: 0.05, meanIsKnown: true))

        for res in [res1, res2] {
            #expect(res.df != nil)
            #expect(res.ratio != nil)
            #expect(res.testStatisticValue != nil)
            #expect(res.p1Value != nil)
            #expect(res.p2Value != nil)
            #expect(res.sampleSize != nil)
            #expect(res.sigmaUEQs0 != nil)
            #expect(res.sigmaLTEs0 != nil)
            #expect(res.sigmaGTEs0 != nil)
            #expect(res.sd != nil)
            let cvs = try #require(res.criticalValues)
            _ = cvs.oneTailedLowerCV
            _ = cvs.oneTailedUpperCV
            _ = cvs.twoTailedLowerCV
            _ = cvs.twoTailedUpperCV
            #expect(res.meanAssumedKnown != nil)
            #expect(!res.description.isEmpty)
        }
    }

    @Test("Levene trimmed-mean variant runs and tags the result")
    func leveneTrimmedMeanResultAccessors() throws {
        let groups = [
            makeExamine([6, 7, 8, 9, 5], name: "g1"),
            makeExamine([3, 4, 2, 5, 1], name: "g2"),
            makeExamine([10, 12, 9, 11, 14], name: "g3")
        ]
        let result = try #require(
            try Inferential<Double>.HypothesisTesting.VarianceTests.leveneTest(
                samples: groups,
                alpha: 0.05,
                testType: .trimmedMean
            )
        )
        #expect(result.testType == .levene(.trimmedMean))
        #expect(result.pValue != nil)
        #expect(result.testStatistic != nil)
        #expect(!result.description.isEmpty)
    }

    @Test("Variance test type enum exposes a stable description")
    func varianceTestTypeDescriptions() {
        #expect(Inferential<Double>.HypothesisTesting.VarianceTests.VarianceTestType.bartlett.description.contains("Bartlett"))
        #expect(Inferential<Double>.HypothesisTesting.VarianceTests.VarianceTestType.levene(.median).description.contains("Levene"))
    }
}
