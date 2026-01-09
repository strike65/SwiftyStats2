/// Tests for location, dispersion, and robust statistics in SSExamine.

import Foundation
import Testing
@testable import SwiftyStats
import SwiftyBoost

@Suite("SSExamine location and dispersion metrics")
struct SSExamineLocationDispersionTests {

    private let baseValues: [Double] = [1.0, 2.0, 2.0, 4.0, 5.0, 7.0]

    private func makeExamine(_ values: [Double]) -> SSExamine<Double, Double> {
        let ex = SSExamine<Double, Double>()
        #expect(ex.append(contentOf: values))
        return ex
    }

    private func approxEqual(_ lhs: Double, _ rhs: Double, tolerance: Double = 1e-9) -> Bool {
        let scale = max(1.0, max(abs(lhs), abs(rhs)))
        return abs(lhs - rhs) <= tolerance * scale
    }

    @Test("Mode, commonest, scarcest, and extrema agree")
    func modeCommonestScarcest() throws {
        let ex = makeExamine(baseValues)

        let mode = try #require(ex.mode)
        #expect(mode == [2.0])

        let common = try #require(ex.commonest)
        #expect(common == mode)

        let scarce = try #require(ex.scarcest)
        #expect(Set(scarce) == Set([1.0, 4.0, 5.0, 7.0]))

        #expect(ex.largest == 7.0)
        #expect(ex.smallest == 1.0)
        let extremes = try #require(ex.smallestAndLargest)
        #expect(extremes.smallest == 1.0)
        #expect(extremes.largest == 7.0)
    }

    @Test("Range-derived metrics are consistent")
    func rangeMetrics() throws {
        let ex = makeExamine(baseValues)

        let range = try #require(ex.range)
        #expect(approxEqual(range, 6.0))

        let mid = try #require(ex.midRange)
        #expect(approxEqual(mid, 4.0))
    }

    @Test("Quantile, quartile, and median values match Type-7 defaults")
    func quantilesAndQuartiles() throws {
        let ex = makeExamine(baseValues)

        #expect(throws: Error.self) { try ex.quantile(q: -0.1) }
        #expect(throws: Error.self) { try ex.quantile(q: 1.1) }

        let q0Opt = try ex.quantile(q: 0.0)
        let q0 = try #require(q0Opt)
        #expect(q0 == 1.0)
        let q25Opt = try ex.quantile(q: 0.25)
        let q25 = try #require(q25Opt)
        #expect(approxEqual(q25, 2.0))
        let q50Opt = try ex.quantile(q: 0.5)
        let q50 = try #require(q50Opt)
        #expect(approxEqual(q50, 3.0))
        let q75Opt = try ex.quantile(q: 0.75)
        let q75 = try #require(q75Opt)
        #expect(approxEqual(q75, 4.75))
        let q100Opt = try ex.quantile(q: 1.0)
        let q100 = try #require(q100Opt)
        #expect(q100 == 7.0)

        let quartOpt = try ex.quartile()
        let quart = try #require(quartOpt)
        #expect(approxEqual(quart.lower, 2.0))
        #expect(approxEqual(quart.median, 3.0))
        #expect(approxEqual(quart.upper, 4.75))

        let median = try #require(ex.median)
        #expect(approxEqual(median, 3.0))

        // Explicit Type-2 request should still reproduce historical expectations.
        let q75Type2Opt = try ex.quantile(q: 0.75, quantileType: .type2)
        let q75Type2 = try #require(q75Type2Opt)
        #expect(approxEqual(q75Type2, 5.0))
    }

    @Test("Hyndman-Fan quantiles track all nine reference definitions")
    func hyndmanFanQuantileFamilies() throws {
        let ex = makeExamine([1.0, 2.0, 3.0, 4.0])
        let q25 = 0.25
        let expectations: [(HyndmanFanQuantileType, Double)] = [
            (.type1, 1.0),
            (.type2, 1.5),
            (.type3, 1.0),
            (.type4, 1.0),
            (.type5, 1.5),
            (.type6, 1.25),
            (.type7, 1.75),
            (.type8, 1.4166666666666667),
            (.type9, 1.4375)
        ]
        for (type, expected) in expectations {
            let valueOpt = try ex.quantile(q: q25, quantileType: type)
            let value = try #require(valueOpt)
            #expect(approxEqual(Double(value), expected))
        }

        let q75 = 0.75
        let type7UpperOpt = try ex.quantile(q: q75, quantileType: .type7)
        let type7Upper = try #require(type7UpperOpt)
        #expect(approxEqual(Double(type7Upper), 3.25))

        let staticResultOpt = try ex.quantile(q: q25, quantileType: .type7)
        let staticResult = try #require(staticResultOpt)
        #expect(approxEqual(Double(staticResult), 1.75))
    }

    @Test("Interquartile-related measures produce expected values")
    func interquartileMetrics() throws {
        let ex = makeExamine(baseValues)

        let iqr = try #require(ex.interquartileRange)
        #expect(approxEqual(iqr, 2.75))

        let customIqrOpt = try ex.interquantileRange(lower: 0.25, upper: 0.75)
        let customIqr = try #require(customIqrOpt)
        #expect(approxEqual(customIqr, 2.75))

        // Order of arguments should not matter because result uses magnitude.
        let reversedOpt = try ex.interquantileRange(lower: 0.75, upper: 0.25)
        let reversed = try #require(reversedOpt)
        #expect(approxEqual(reversed, 2.75))

        let qDev = try #require(ex.quartileDeviation)
        #expect(approxEqual(qDev, 1.375))

        let rqDist = try #require(ex.relativeQuartileDistance)
        #expect(approxEqual(rqDist, 11.0 / 12.0))
    }

    @Test("Family of means (arithmetic, geometric, harmonic, Lehmer, power) stay coherent")
    func meanFamily() throws {
        let ex = makeExamine(baseValues)

        let arithmetic = try #require(ex.arithmeticMean)
        #expect(approxEqual(arithmetic, 3.5))

        let logSum = try #require(ex.logSum)
        let expectedLogSum = baseValues.reduce(0.0) { $0 + log($1) }
        #expect(approxEqual(logSum, expectedLogSum))

        let geom = try #require(ex.geometricMean)
        let expectedGeom = exp(expectedLogSum / Double(baseValues.count))
        #expect(approxEqual(geom, expectedGeom))
        #expect(approxEqual(try #require(ex.powerMean(0)), expectedGeom))

        let harmonic = try #require(ex.harmonicMean)
        let expectedHarm = Double(baseValues.count) / baseValues.reduce(0.0) { $0 + 1.0 / $1 }
        #expect(approxEqual(harmonic, expectedHarm))

        let contra = try #require(ex.contraharmonicMean)
        let expectedContra = baseValues.reduce(0.0) { $0 + $1 * $1 } / baseValues.reduce(0.0, +)
        #expect(approxEqual(contra, expectedContra))
        #expect(approxEqual(try #require(ex.lehmerMean(2)), expectedContra))

        #expect(approxEqual(try #require(ex.lehmerMean(1)), arithmetic))

        let lehmer3 = try #require(ex.lehmerMean(3))
        let expectedLehmer3 = baseValues.reduce(0.0) { $0 + pow($1, 3.0) } / baseValues.reduce(0.0) { $0 + pow($1, 2.0) }
        #expect(approxEqual(lehmer3, expectedLehmer3))

        let rms = try #require(ex.powerMean(2))
        let expectedRms = sqrt(baseValues.reduce(0.0) { $0 + $1 * $1 } / Double(baseValues.count))
        #expect(approxEqual(rms, expectedRms))

        let power3 = try #require(ex.powerMean(3))
        let expectedPower3 = pow(baseValues.reduce(0.0) { $0 + pow($1, 3.0) } / Double(baseValues.count), 1.0 / 3.0)
        #expect(approxEqual(power3, expectedPower3))
    }

    @Test("Trimmed, winsorized, and gastwirth means align with definitions")
    func robustMeans() throws {
        let ex = makeExamine(baseValues)

        let trimmed = try #require(ex.trimmedMean(alpha: 0.2))
        #expect(approxEqual(trimmed, 3.25))

        let winsorized = try #require(ex.winsorizedMean(alpha: 0.2))
        #expect(approxEqual(winsorized, 20.0 / 6.0))

        let winsorised = try #require(ex.winsorisedMean(alpha: 0.2))
        #expect(approxEqual(winsorised, winsorized))

        let gastwirth = try #require(ex.gastwirthEstimate())
        let expectedGastwirth = (1.0 / 3.0) * 2.0 + (2.0 / 5.0) * 3.0 + (1.0 / 3.0) * 13.0 / 3.0
        #expect(approxEqual(gastwirth, expectedGastwirth))
    }

    @Test("Median absolute deviation variants produce expected scale estimates")
    func medianAbsoluteDeviationVariants() throws {
        let ex = makeExamine(baseValues)
        let center = try #require(ex.median)

        let mad = try #require(ex.medianAbsoluteDeviation(center: center))
        #expect(approxEqual(Double(mad), 1.5))

        let sMad = try #require(ex.standardizedMedianAbsoluteDeviation(center: center))
        let expected = 1.5 * 1.4826022185056018
        #expect(approxEqual(Double(sMad), expected))

        let empty = SSExamine<Double, Double>()
        #expect(empty.medianAbsoluteDeviation(center: 0) == nil)
        #expect(empty.standardizedMedianAbsoluteDeviation(center: 0) == nil)
    }

    @Test("Box-whisker summary matches Tukey fences and notch definitions")
    func boxWhiskerSummary() throws {
        let ex = makeExamine(baseValues)
        let bwValue = try ex.boxWhisker()
        let bw = try #require(bwValue)

        #expect(approxEqual(Double(try #require(bw.median)), 3.0))
        #expect(approxEqual(Double(try #require(bw.q25)), 2.0))
        #expect(approxEqual(Double(try #require(bw.q75)), 4.75))
        #expect(approxEqual(Double(try #require(bw.iqr)), 2.75))
        #expect(approxEqual(Double(try #require(bw.lWhiskerExtreme)), 1.0))
        #expect(approxEqual(Double(try #require(bw.uWhiskerExtreme)), 7.0))

        let expectedNotch = 1.57 * 2.75 / sqrt(Double(baseValues.count))
        #expect(approxEqual(Double(try #require(bw.uNotch)), 3.0 + expectedNotch, tolerance: 1e-12))
        #expect(approxEqual(Double(try #require(bw.lNotch)), 3.0 - expectedNotch, tolerance: 1e-12))

        #expect(bw.fences?.isEmpty == true)
        #expect(bw.outliers?.isEmpty == true)
    }

    @Test("Qn robust scale matches Rousseeuw-Croux definition with finite-sample correction")
    func qnScaleEstimate() throws {
        let ex = makeExamine(baseValues)
        let qn = try #require(ex.Qn)
        let expectedQn = 2.717115016
        #expect(approxEqual(Double(qn), expectedQn, tolerance: 1e-12))

        let constant = SSExamine<Double, Double>()
        #expect(constant.append(contentOf: Array(repeating: 5.0, count: 4)))
        #expect(constant.Qn == 0.0)
    }

    @Test("Sn robust scale matches Rousseeuw-Croux definition with finite-sample correction")
    func snScaleEstimate() throws {
        let ex = makeExamine(baseValues)
        let sn = try #require(ex.Sn)
        let expectedSn = 2.3685036
        #expect(approxEqual(Double(sn), expectedSn, tolerance: 1e-12))

        let constant = SSExamine<Double, Double>()
        #expect(constant.append(contentOf: [3.0, 3.0, 3.0]))
        #expect(constant.Sn == 0.0)
    }

    @Test("Log sum handles zero, negative, and empty datasets")
    func logSumEdgeCases() throws {
        let ex = makeExamine(baseValues)
        let positiveLogSum = try #require(ex.logSum)
        let expected = baseValues.reduce(0.0) { $0 + log($1) }
        #expect(approxEqual(positiveLogSum, expected))

        let zeroEx = makeExamine([2.0, 0.0, 3.0])
        let zeroLog = try #require(zeroEx.logSum)
        #expect(zeroLog == -Double.infinity)

        let negativeEx = makeExamine([1.0, -2.0])
        #expect(negativeEx.logSum == nil)

        let emptyEx = SSExamine<Double, Double>()
        if let emptyLog = emptyEx.logSum {
            #expect(emptyLog == 0.0)
        } else {
            #expect(Bool(false), "Expected logSum to be 0 for an empty examine")
        }
    }

    @Test("Unnormalized moments and variances match manual calculations")
    func momentsAndVariance() throws {
        let ex = makeExamine(baseValues)

        let arithmetic = try #require(ex.arithmeticMean)
        let moment2 = try #require(ex.unnormalizedMoment(order: 2, about: arithmetic))
        #expect(approxEqual(moment2, 25.5))

        let rawMoment0 = try #require(ex.unnormalizedMoment(order: 0, about: 0.0))
        #expect(approxEqual(rawMoment0, Double(baseValues.count)))

        let sampleVar = try #require(ex.sampleVariance)
        #expect(approxEqual(sampleVar, 5.1))

        let populationVar = try #require(ex.populationVariance)
        #expect(approxEqual(populationVar, 4.25))
    }

    @Test("Variance aliases, semivariances, and dispersion derivatives stay in sync")
    func varianceAliasesAndSemivariances() throws {
        let ex = makeExamine(baseValues)
        let popVar = try #require(ex.populationVariance)
        let sampleVar = try #require(ex.sampleVariance)

        #expect(ex.biasedVariance == popVar)
        #expect(ex.unbiasedVariance == sampleVar)

        let popSd = try #require(ex.populationStandardDeviation)
        let sampleSd = try #require(ex.sampleStandardDeviation)
        #expect(approxEqual(Double(popSd * popSd), Double(popVar)))
        #expect(approxEqual(Double(sampleSd * sampleSd), Double(sampleVar)))

        let mean = baseValues.reduce(0.0, +) / Double(baseValues.count)
        let lowerTerms = baseValues.filter { $0 < mean }.map { pow($0 - mean, 2.0) }
        let expectedSemiLower = lowerTerms.isEmpty ? 0.0 : lowerTerms.reduce(0.0, +) / Double(lowerTerms.count)
        let upperTerms = baseValues.filter { $0 > mean }.map { pow($0 - mean, 2.0) }
        let expectedSemiUpper = upperTerms.isEmpty ? 0.0 : upperTerms.reduce(0.0, +) / Double(upperTerms.count)
        let downsideTerms = baseValues.map { min(0.0, $0 - mean) }.map { pow($0, 2.0) }
        let expectedDownside = downsideTerms.reduce(0.0, +) / Double(baseValues.count)

        let semiLower = try #require(ex.semiVarianceLower)
        let semiUpper = try #require(ex.semiVarianceUpper)
        let semiDown = try #require(ex.semiVarianceDownside)

        #expect(approxEqual(Double(semiLower), expectedSemiLower))
        #expect(approxEqual(Double(semiUpper), expectedSemiUpper))
        #expect(approxEqual(Double(semiDown), expectedDownside))

        let sem = try #require(ex.standardError)
        let expectedSem = Double(sampleSd) / sqrt(Double(baseValues.count))
        #expect(approxEqual(Double(sem), expectedSem))

        let cv = try #require(ex.coefficientOfVariation)
        let cvAlias = try #require(ex.cv)
        let expectedCv = Double(sampleSd) / mean
        #expect(approxEqual(Double(cv), expectedCv))
        #expect(cv == cvAlias)
    }

    @Test("Dispersion metrics gracefully handle empty, singleton, and non-numeric data")
    func dispersionNilCases() throws {
        let empty = SSExamine<Double, Double>()
        #expect(empty.range == nil)
        #expect(empty.interquartileRange == nil)
        let emptyIqr = try? empty.interquantileRange(lower: 0.25, upper: 0.75)
        #expect(emptyIqr == nil)
        #expect(empty.midRange == nil)
        #expect(empty.quartileDeviation == nil)
        #expect(empty.relativeQuartileDistance == nil)
        #expect(empty.sampleVariance == nil)
        #expect(empty.populationVariance == nil)
        #expect(empty.populationStandardDeviation == nil)
        #expect(empty.sampleStandardDeviation == nil)
        #expect(empty.semiVarianceUpper == nil)
        #expect(empty.semiVarianceLower == nil)
        #expect(empty.semiVarianceDownside == nil)
        #expect(empty.standardError == nil)
        #expect(empty.coefficientOfVariation == nil)
        #expect(empty.cv == nil)
        #expect(empty.centralMoment(order: 2) == nil)
        #expect(empty.originMoment(order: 3) == nil)
        #expect(empty.standardizedMoment(order: 3) == nil)
        #expect(empty.rawMoment(order: 2, about: 0.0) == nil)
        #expect(empty.skewness == nil)
        #expect(empty.kurtosis == nil)
        #expect(empty.kurtosisExcess == nil)

        let singleton = SSExamine<Double, Double>()
        singleton.append(4.2)
        #expect(singleton.sampleSize == 1)
        #expect(singleton.sampleVariance == nil)
        #expect(singleton.populationVariance == nil)
        #expect(singleton.sampleStandardDeviation == nil)
        #expect(singleton.populationStandardDeviation == nil)
        #expect(singleton.standardError == nil)
        #expect(singleton.semiVarianceUpper == nil)
        #expect(singleton.semiVarianceLower == nil)
        #expect(singleton.semiVarianceDownside == nil)
        #expect(singleton.coefficientOfVariation == nil)
        #expect(singleton.cv == nil)
        #expect(singleton.centralMoment(order: 3) == 0)
        #expect(singleton.standardizedMoment(order: 3) == nil)
        #expect(singleton.skewness == nil)
        #expect(singleton.kurtosis == nil)
        #expect(singleton.kurtosisExcess == nil)

        let strings = SSExamine<String, Double>()
        #expect(strings.append(contentOf: ["a", "b", "c"]))
        #expect(strings.isNumeric == false)
        #expect(strings.range == nil)
        #expect(strings.interquartileRange == nil)
        let stringIqr = try? strings.interquantileRange(lower: 0.25, upper: 0.75)
        #expect(stringIqr == nil)
        #expect(strings.midRange == nil)
        #expect(strings.quartileDeviation == nil)
        #expect(strings.relativeQuartileDistance == nil)
        #expect(strings.sampleVariance == nil)
        #expect(strings.populationVariance == nil)
        #expect(strings.sampleStandardDeviation == nil)
        #expect(strings.populationStandardDeviation == nil)
        #expect(strings.standardError == nil)
        #expect(strings.semiVarianceUpper == nil)
        #expect(strings.semiVarianceLower == nil)
        #expect(strings.semiVarianceDownside == nil)
        #expect(strings.coefficientOfVariation == nil)
        #expect(strings.cv == nil)
        #expect(strings.centralMoment(order: 2) == nil)
        #expect(strings.originMoment(order: 3) == nil)
        #expect(strings.standardizedMoment(order: 3) == nil)
        #expect(strings.rawMoment(order: 2, about: 0.0) == nil)
        #expect(strings.skewness == nil)
        #expect(strings.kurtosis == nil)
        #expect(strings.kurtosisExcess == nil)
    }

    @Test("Mean absolute deviations and Gini mean difference match hand calculations")
    func meanDeviationsAndGini() throws {
        let values: [Double] = [1.0, 2.0, 3.0, 4.0, 7.0]
        let ex = makeExamine(values)

        let center: Double = 3.0
        let upperMad = try #require(ex.upperMeanAbsoluteDeviation(from: center))
        let lowerMad = try #require(ex.lowerMeanAbsoluteDeviation(from: center))

        // Upper deviations: {4, 7} -> |4-3|=1, |7-3|=4 -> mean = 2.5
        #expect(approxEqual(Double(upperMad), 2.5))
        // Lower deviations: {1,2} -> |1-3|=2, |2-3|=1 -> mean = 1.5
        #expect(approxEqual(Double(lowerMad), 1.5))

        // Gini mean difference = 2/(n^2) * sum_{i<j} |x_i - x_j| * n/(n-1) (SwiftyStats implementation)
        var sumPairs = 0.0
        for i in 0..<(values.count - 1) {
            for j in (i + 1)..<values.count {
                sumPairs += abs(values[i] - values[j])
            }
        }
        let n = Double(ex.sampleSize)
        let expectedGiniMeanDifference = (2 * sumPairs) / (n * n) * n / (n - 1)
        let giniMD = try #require(ex.giniMeanDifference)
        #expect(approxEqual(Double(giniMD), expectedGiniMeanDifference))
    }

    @Test("Relative quartile distance returns nil when the median is zero")
    func relativeQuartileDistanceZeroMedian() throws {
        let symmetric = SSExamine<Double, Double>()
        #expect(symmetric.append(contentOf: [-2.0, -1.0, 0.0, 1.0, 2.0]))
        #expect(symmetric.relativeQuartileDistance == nil)
    }

    @Test("Interquantile range rejects quantiles outside the unit interval")
    func interquantileRangeInvalidBounds() throws {
        let ex = makeExamine(baseValues)
        #expect(throws: Error.self) { try ex.interquantileRange(lower: -0.1, upper: 0.9) }
        #expect(throws: Error.self) { try ex.interquantileRange(lower: 0.25, upper: 1.1) }
    }

    @Test("Moment helpers and shape statistics line up with manual calculations")
    func momentHelpersAndShapeStatistics() throws {
        let values = baseValues
        let ex = makeExamine(values)
        let n = Double(values.count)
        let mean = values.reduce(0.0, +) / n
        let central2 = values.reduce(0.0) { $0 + pow($1 - mean, 2.0) } / n
        let central3 = values.reduce(0.0) { $0 + pow($1 - mean, 3.0) } / n
        let central4 = values.reduce(0.0) { $0 + pow($1 - mean, 4.0) } / n
        let origin3 = values.reduce(0.0) { $0 + pow($1, 3.0) } / n
        let rawAboutMean = values.reduce(0.0) { $0 + pow($1 - mean, 2.0) }
        let sampleVar = try #require(ex.sampleVariance)
        let sampleSd = sqrt(Double(sampleVar))

        let centralMoment2 = try #require(ex.centralMoment(order: 2))
        #expect(approxEqual(Double(centralMoment2), central2))

        let centralMoment3 = try #require(ex.centralMoment(order: 3))
        #expect(approxEqual(Double(centralMoment3), central3))

        let centralMoment4 = try #require(ex.centralMoment(order: 4))
        #expect(approxEqual(Double(centralMoment4), central4))

        let originMoment3 = try #require(ex.originMoment(order: 3))
        #expect(approxEqual(Double(originMoment3), origin3))

        let standardizedMoment3 = try #require(ex.standardizedMoment(order: 3))
        let expectedStandardized3 = central3 / pow(sampleSd, 3.0)
        #expect(approxEqual(Double(standardizedMoment3), expectedStandardized3))

        let standardizedMoment4 = try #require(ex.standardizedMoment(order: 4))
        let expectedStandardized4 = central4 / pow(sampleSd, 4.0)
        #expect(approxEqual(Double(standardizedMoment4), expectedStandardized4))

        let rawMomentMean = try #require(ex.rawMoment(order: 2, about: mean))
        #expect(approxEqual(Double(rawMomentMean), rawAboutMean))

        let rawMomentZero = try #require(ex.rawMoment(order: 0, about: 0.0))
        #expect(approxEqual(Double(rawMomentZero), n))

        let skewness = try #require(ex.skewness)
        #expect(approxEqual(Double(skewness), expectedStandardized3))

        let kurtosis = try #require(ex.kurtosis)
        let expectedKurtosis = central4 / pow(sampleSd, 4.0)
        #expect(approxEqual(Double(kurtosis), expectedKurtosis))

        let kurtosisExcess = try #require(ex.kurtosisExcess)
        #expect(approxEqual(Double(kurtosisExcess), expectedKurtosis - 3.0))
    }

    @Test("Confidence interval helpers reproduce analytic intervals and errors")
    func confidenceIntervalHelpers() throws {
        let ex = makeExamine(baseValues)
        let n = Double(baseValues.count)
        let mean = baseValues.reduce(0.0, +) / n

        let alpha: Double = 0.1
        let zQuantile = try SwiftyBoost.Distribution.Normal().quantile(1.0 - alpha / 2.0)
        let knownSD: Double = 2.0

        let ciKnownOptional = try ex.normalCI(alpha: alpha, sd: knownSD)
        let ciKnown = try #require(ciKnownOptional)
        let halfKnown = Double(ciKnown.upperBound - ciKnown.lowerBound) / 2.0
        let expectedHalfKnown = zQuantile * knownSD / n.squareRoot()
        #expect(approxEqual(halfKnown, expectedHalfKnown))
        #expect(approxEqual(Double(ciKnown.lowerBound), mean - expectedHalfKnown))
        #expect(approxEqual(Double(ciKnown.upperBound), mean + expectedHalfKnown))

        let populationSd = try #require(ex.populationStandardDeviation)
        let ciEstimatedOptional = try ex.normalCI(alpha: alpha, sd: nil)
        let ciEstimated = try #require(ciEstimatedOptional)
        let expectedHalfEstimated = zQuantile * Double(populationSd) / n.squareRoot()
        let halfEstimated = Double(ciEstimated.upperBound - ciEstimated.lowerBound) / 2.0
        #expect(approxEqual(halfEstimated, expectedHalfEstimated))

        let tAlpha: Double = 0.1
        let tQuantile = try SwiftyBoost.Distribution.StudentT(degreesOfFreedom: n - 1.0).quantile(1.0 - tAlpha / 2.0)
        let studentOptional = try ex.studentTCI(alpha: tAlpha)
        let student = try #require(studentOptional)
        let sampleSd = try #require(ex.sampleStandardDeviation)
        let expectedHalfStudent = tQuantile * Double(sampleSd) / n.squareRoot()
        let halfStudent = Double(student.upperBound - student.lowerBound) / 2.0
        #expect(approxEqual(halfStudent, expectedHalfStudent, tolerance: 1e-9))

        let meanCI = try #require(ex.meanCI)
        let student95Optional = try ex.studentTCI(alpha: 0.05)
        let student95 = try #require(student95Optional)
        #expect(approxEqual(Double(meanCI.lowerBound), Double(student95.lowerBound), tolerance: 1e-12))
        #expect(approxEqual(Double(meanCI.upperBound), Double(student95.upperBound), tolerance: 1e-12))

        #expect(throws: SSError.self) {
            _ = try ex.normalCI(alpha: 1.5, sd: knownSD)
        }

        #expect(throws: SSError.self) {
            _ = try ex.studentTCI(alpha: -0.1)
        }
    }
}
