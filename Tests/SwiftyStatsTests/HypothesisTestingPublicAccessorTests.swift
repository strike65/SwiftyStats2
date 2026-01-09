import Testing
@testable import SwiftyStats

/// Exercises every public HypothesisTesting accessor at least once.
struct HypothesisTestingPublicAccessorTests {
    private func makeExamine(_ values: [Double], name: String? = nil) -> SSExamine<Double, Double> {
        SSExamine<Double, Double>(usingArray: values, name: name, characterSet: nil)
    }

    private func makeStringExamine(_ values: [String]) -> SSExamine<String, Double> {
        SSExamine<String, Double>(usingArray: values, name: nil, characterSet: nil)
    }

    private func makeUniforms(n: Int, seed: UInt64) -> [Double] {
        precondition(n >= 0)
        var state = seed
        var out: [Double] = []
        out.reserveCapacity(n)
        for _ in 0..<n {
            state = state &* 6364136223846793005 &+ 1442695040888963407
            let top53 = state >> 11
            let u = Double(top53) / Double(1 << 53)
            out.append(u)
        }
        return out
    }

    @Test("Binomial test overloads return consistent tails")
    func binomialTestOverloadsAreConsistent() throws {
        let x = 7
        let n = 20
        let p0 = 0.3

        let tails = try #require(try Parametric<Double>.binomialTest(successes: x, numberOfTrials: n, p0: p0))
        #expect(tails.less >= 0 && tails.less <= 1)
        #expect(tails.greater >= 0 && tails.greater <= 1)
        #expect(tails.twoSided >= 0 && tails.twoSided <= 1)

        let pLess = try #require(try Parametric<Double>.binomialTest(successes: x, numberOfTrials: n, p0: p0, alternative: .less))
        let pGreater = try #require(try Parametric<Double>.binomialTest(successes: x, numberOfTrials: n, p0: p0, alternative: .greater))
        let pTwoSided = try #require(try Parametric<Double>.binomialTest(successes: x, numberOfTrials: n, p0: p0, alternative: .twoSided))
        #expect(abs(pLess - tails.less) < 1e-12)
        #expect(abs(pGreater - tails.greater) < 1e-12)
        #expect(abs(pTwoSided - tails.twoSided) < 1e-12)
    }

    @Test("Binomial test categorical overload populates result accessors")
    func binomialTestCategoricalResultAccessors() throws {
        let data = ["s", "f", "s", "s", "f", "s", "f", "s", "f", "s"]
        let ex = makeStringExamine(data)

        let arrayRes = try #require(
            try Parametric<Double>.binomialTest(
                data: data,
                p0: 0.5,
                characterSet: nil,
                successCode: "s",
                alternative: .twoSided
            )
        )
        let exRes = try #require(
            try Parametric<Double>.binomialTest(
                data: ex,
                p0: 0.5,
                characterSet: nil,
                successCode: "s",
                alternative: .twoSided
            )
        )

        for res in [arrayRes, exRes] {
            #expect(res.nTrials != nil)
            #expect(res.nSuccess != nil)
            #expect(res.nFailure != nil)
            #expect(res.pValueExact != nil)
            #expect(res.probSuccess != nil)
            #expect(res.probFailure != nil)
            #expect(res.probTest != nil)
            #expect(res.successCode != nil)
            #expect(res.confIntJeffreys != nil)
            #expect(res.confIntClopperPearson != nil)
            #expect(!res.description.isEmpty)
        }
    }

    @Test("Anderson-Darling normality test populates all result fields")
    func andersonDarlingNormalityTestResultAccessors() throws {
        let data: [Double] = [-1.2, -0.8, -0.5, -0.2, 0.0, 0.1, 0.4, 0.7, 1.1, 1.4]
        let ex = makeExamine(data)

        let arrayRes = try #require(try GoodnessOfFit<Double>.andersonDarlingNormalityTest(data: data, alpha: 0.05))
        let exRes = try #require(try GoodnessOfFit<Double>.andersonDarlingNormalityTest(data: ex, alpha: 0.05))

        for res in [arrayRes, exRes] {
            #expect(res.pValue != nil)
            #expect(res.AD != nil)
            #expect(res.sampleSize == data.count)
            #expect(res.stdDev != nil)
            #expect(res.variance != nil)
            #expect(res.mean != nil)
            #expect(res.isNormal != nil)
            #expect(!res.description.isEmpty)
        }
    }

    @Test("Grubbs test returns a fully-populated result")
    func grubbsTestResultAccessors() throws {
        let values: [Double] = [0, 0, 0, 10]
        let ex = makeExamine(values)
        let mean = try #require(ex.arithmeticMean)
        let sd = try #require(ex.sampleStandardDeviation)

        let res = try #require(
            try Outliers<Double>.grubbsTest(
                values: values,
                mean: mean,
                standardDeviation: sd,
                alpha: 0.05
            )
        )
        #expect(res.criticalValue != nil)
        #expect(res.largest != nil)
        #expect(res.smallest != nil)
        #expect(res.sampleSize == values.count)
        #expect(res.maxDiff != nil)
        #expect(res.mean != nil)
        #expect(res.G != nil)
        #expect(res.stdDev != nil)
        #expect(res.hasOutliers != nil)
        #expect(!res.description.isEmpty)
    }

    @Test("Rosner ESD result exposes expected accessors")
    func rosnerESDResultAccessors() throws {
        let values: [Double] = [1, 2, 2, 3, 100, 3, 2, 1]
        let res = try #require(
            try Outliers<Double>.rosnerESD(
                data: values,
                alpha: 0.05,
                maxOutliers: 2,
                testType: .bothTails
            )
        )
        #expect(res.stdDeviations != nil)
        #expect(res.itemsRemoved != nil)
        #expect(res.testStatistics != nil)
        #expect(res.lambdas != nil)
        #expect(res.countOfOutliers != nil)
        #expect(res.outliers != nil)
        #expect(res.alpha != nil)
        #expect(res.maxOutliers != nil)
        #expect(res.testType != nil)
        #expect(res.means != nil)
        #expect(!res.description.isEmpty)
    }

    @Test("Two-sample Kolmogorov-Smirnov test runs for both overloads")
    func ksTwoSampleTestOverloads() throws {
        let a = makeUniforms(n: 30, seed: 0x1111_1111)
        let b = makeUniforms(n: 30, seed: 0x2222_2222).map { $0 + 0.5 }
        let ex1 = makeExamine(a)
        let ex2 = makeExamine(b)

        let res1 = try #require(try NonParametric<Double>.ksTest(set1: ex1, set2: ex2))
        let res2 = try #require(try NonParametric<Double>.ksTest(set1: a, set2: b))

        for res in [res1, res2] {
            #expect(res.dMaxPos != nil)
            #expect(res.dMaxNeg != nil)
            #expect(res.dMaxAbs != nil)
            #expect(res.zStatistic != nil)
            #expect(res.p2Value != nil)
            #expect(res.sampleSize1 != nil)
            #expect(res.sampleSize2 != nil)
            #expect(!res.description.isEmpty)
        }
    }

    @Test("One-sample KS statistic and bootstrap helpers return coherent values")
    func ksOneSampleStatisticAndBootstrapHelpers() throws {
        let data = makeUniforms(n: 25, seed: 0xABCDEF)

        let stat = try GoodnessOfFit<Double>.ksStatistic(data: data) { x in
            x <= 0 ? 0 : (x >= 1 ? 1 : x)
        }
        #expect(stat.d >= 0)
        #expect(stat.d == max(stat.dPlus, stat.dMinus))

        let fit: ([Double]) -> [Double] = { xs in
            [xs.min() ?? 0.0, xs.max() ?? 1.0]
        }
        let cdfFrom: ([Double]) -> (Double) -> Double = { theta in
            let a = theta[0]
            let b = theta[1]
            let denom = max(b - a, Double.leastNonzeroMagnitude)
            return { x in
                if x <= a { return 0 }
                if x >= b { return 1 }
                return (x - a) / denom
            }
        }
        let sampler: (Int, [Double], UInt64?) -> [Double] = { n, theta, seed in
            let u = makeUniforms(n: n, seed: seed ?? 0x1234)
            let a = theta[0]
            let b = theta[1]
            return u.map { a + (b - a) * $0 }
        }
        var seedState: UInt64 = 0xB00B_1E5
        let samplerSeed: () -> UInt64? = {
            seedState &+= 1
            return seedState
        }

        let boot = try GoodnessOfFit<Double>.bootstrapOneSampleKS(
            data: data,
            targetDistribution: .uniform,
            fit: fit,
            cdfFrom: cdfFrom,
            sampler: sampler,
            samplerSeed: samplerSeed,
            reestimateEachReplicate: false,
            B: 50
        )
        #expect(boot.n == data.count)
        #expect(boot.b == 50)
        #expect(boot.d.isFinite)
        #expect(boot.dPlus.isFinite)
        #expect(boot.dMinus.isFinite)
        #expect(boot.pBootstrap > 0 && boot.pBootstrap <= 1)
        #expect(!boot.thetaHat.isEmpty)
        #expect(!boot.distribution.isEmpty)
    }

    @Test("One-sample KS entry point returns a populated result for Bernoulli data")
    func ksOneSampleEntryPointBernoulliReturnsFallbackOnEmptyInput() {
        let res = GoodnessOfFit<Double>.ksTestOneSample(data: [], testTarget: .bernoulli)
        #expect(res.n == 0)
        #expect(res.b == 0)
        #expect(res.d.isNaN)
        #expect(res.dPlus.isNaN)
        #expect(res.dMinus.isNaN)
        #expect(res.pBootstrap.isNaN)
        #expect(res.thetaHat.isEmpty)
        #expect(!res.distribution.isEmpty)
    }

    @Test("One-sample KS entry point is reproducible with a fixed seed")
    func ksOneSampleEntryPointRespectsSeed() {
        let data = makeUniforms(n: 40, seed: 0x5151_5151)

        let res1 = GoodnessOfFit<Double>.ksTestOneSample(
            data: data,
            testTarget: .uniform,
            boostrapCount: 40,
            seed: 0xCAFE_BABE
        )
        let res2 = GoodnessOfFit<Double>.ksTestOneSample(
            data: data,
            testTarget: .uniform,
            boostrapCount: 40,
            seed: 0xCAFE_BABE
        )

        #expect(abs(res1.pBootstrap - res2.pBootstrap) < 1e-12)
        #expect(abs(res1.d - res2.d) < 1e-12)
        #expect(res1.b == res2.b)
        #expect(res1.thetaHat == res2.thetaHat)
    }

    @Test("KS bootstrap helpers surface public errors")
    func ksErrorsAreThrown() throws {
        #expect(throws: GoodnessOfFit<Double>.KSError.self) {
            _ = try GoodnessOfFit<Double>.ksStatistic(data: []) { _ in 0 }
        }

        let data = makeUniforms(n: 10, seed: 0xD00D)
        let fit: ([Double]) -> [Double] = { _ in [0.0, 1.0] }
        let cdfFrom: ([Double]) -> (Double) -> Double = { _ in { x in x <= 0 ? 0 : (x >= 1 ? 1 : x) } }
        let sampler: (Int, [Double], UInt64?) -> [Double] = { n, _, seed in makeUniforms(n: n, seed: seed ?? 0xBEEF) }

        #expect(throws: GoodnessOfFit<Double>.KSError.self) {
            _ = try GoodnessOfFit<Double>.bootstrapOneSampleKS(
                data: data,
                targetDistribution: .uniform,
                fit: fit,
                cdfFrom: cdfFrom,
                sampler: sampler,
                samplerSeed: nil,
                reestimateEachReplicate: false,
                B: 0
            )
        }
    }

    @Test("Kruskal-Wallis test returns a fully-populated result")
    func kruskalWallisResultAccessors() throws {
        let groups: [[Double]] = [
            [1, 2, 3, 4],
            [10, 11, 12, 13],
            [20, 21, 22, 23]
        ]
        let res = try #require(try NonParametric<Double>.kruskalWallisTest(data: groups, alpha: 0.05))
        #expect(res.H_value != nil)
        #expect(res.H_value_corrected != nil)
        #expect(res.pValue != nil)
        #expect(res.nGroups == groups.count)
        #expect(res.df != nil)
        #expect(res.nObservations != nil)
        #expect(res.meanRanks != nil)
        #expect(res.sumRanks != nil)
        #expect(res.cv != nil)
        #expect(res.nTies != nil)
        #expect(res.alpha != nil)
        #expect(!res.description.isEmpty)
    }

    @Test("Friedman test returns a fully-populated result")
    func friedmanResultAccessors() throws {
        let data: [[Double]] = [
            [1, 2, 3],
            [2, 3, 1],
            [3, 1, 2],
            [1, 3, 2]
        ]
        let res = try #require(try NonParametric<Double>.friedmanTest(data: data, alpha: 0.05, tieCorrect: true))
        #expect(res.numberOfBlocks == data.count)
        #expect(res.numberOfTreatments == data[0].count)
        #expect(res.rankSums != nil)
        #expect(res.meanRanks != nil)
        #expect(res.chiSquare != nil)
        #expect(res.chiSquareTieCorrected != nil)
        #expect(res.degreesOfFreedom != nil)
        #expect(res.criticalValueChiSquare != nil)
        #expect(res.pValueChiSquare != nil)
        #expect(res.tieCorrectionFactor != nil)
        #expect(res.tieSum != nil)
        #expect(res.kendallsW != nil)
        #expect(res.imanDavenportF != nil)
        #expect(res.imanDavenportDf1 != nil)
        #expect(res.imanDavenportDf2 != nil)
        #expect(res.criticalValueImanDavenport != nil)
        #expect(res.pValueImanDavenport != nil)
        #expect(!res.description.isEmpty)
    }

    @Test("Mann-Whitney U test returns a populated result for both overloads")
    func mannWhitneyResultAccessors() throws {
        let g1: [Double] = [1, 2, 3, 4, 5]
        let g2: [Double] = [10, 11, 12, 13, 14]
        let ex1 = makeExamine(g1, name: "group 1")
        let ex2 = makeExamine(g2, name: "group 2")

        let res1 = try #require(try NonParametric<Double>.mannWhitneyUTest(group1: ex1, group2: ex2, continuityCorrection: true))
        let res2 = try #require(try NonParametric<Double>.mannWhitneyUTest(group1: g1, group2: g2, continuityCorrection: false))

        for res in [res1, res2] {
            #expect(res.sumRanks1 != nil)
            #expect(res.sumRanks2 != nil)
            #expect(res.meanRank1 != nil)
            #expect(res.meanRank2 != nil)
            #expect(res.sumTiedRanks != nil)
            #expect(res.nTies != nil)
            #expect(res.U != nil)
            #expect(res.W != nil)
            #expect(res.zStat != nil)
            #expect(res.pValueAsymTwoSided != nil)
            #expect(res.pValueExactTwoSided != nil)
            #expect(res.pValueAsym != nil)
            #expect(res.pValueExact != nil)
            #expect(res.effectSize != nil)
            #expect(res.continuityCorrected != nil)
            #expect(!res.description.isEmpty)
        }
    }

    @Test("Mann-Whitney power and sample size helpers return non-nil tuples")
    func mannWhitneyPowerHelpers() throws {
        let planned = try #require(try NonParametric<Double>.nMannWhitney(effectSize: 0.65, alpha: 0.05, alternative: .twoSided, power: 0.8))
        #expect(planned.n1 > 0)
        #expect(planned.n2 > 0)
        #expect(planned.power == 0.8)

        let achieved = try #require(try NonParametric<Double>.powerMannWhitney(effectSize: 0.65, alpha: 0.05, alternative: .twoSided, N: 40))
        #expect(achieved.n1 > 0)
        #expect(achieved.n2 > 0)
        #expect(achieved.power >= 0 && achieved.power <= 1)
    }

    @Test("Exact Mann-Whitney p-value helper returns probabilities in [0, 1]")
    func mannWhitneyExactPValueHelper() throws {
        let pLess = try NonParametric<Double>.MMannWhithneyDistribution.pValue(
            U: 0,
            n1: 3,
            n2: 4,
            alternative: .less
        )
        let pGreater = try NonParametric<Double>.MMannWhithneyDistribution.pValue(
            U: 0,
            n1: 3,
            n2: 4,
            alternative: .greater
        )
        let pTwo = try NonParametric<Double>.MMannWhithneyDistribution.pValue(
            U: 0,
            n1: 3,
            n2: 4,
            alternative: .twoSided
        )
        #expect(pLess >= 0 && pLess <= 1)
        #expect(pGreater >= 0 && pGreater <= 1)
        #expect(pTwo >= 0 && pTwo <= 1)
    }

    @Test("Wilcoxon signed-rank test returns a populated result for both overloads")
    func wilcoxonSignedRankResultAccessors() throws {
        let a: [Double] = [1, 2, 3, 4, 5, 6]
        let b: [Double] = [1, 1, 4, 2, 5, 9]
        let ex1 = makeExamine(a, name: "tier1")
        let ex2 = makeExamine(b, name: "tier2")

        let res1 = try #require(
            try NonParametric<Double>.wilcoxonSignedRankTest(
                tier1: ex1,
                tier2: ex2,
                zeroHandling: .pratt,
                continuityCorrection: true,
                exact: true,
                exactMaxN: 30
            )
        )
        let res2 = try #require(
            try NonParametric<Double>.wilcoxonSignedRankTest(
                tier1: a,
                tier2: b,
                zeroHandling: .discardZeroDifferences,
                continuityCorrection: false,
                exact: true,
                exactMaxN: 30
            )
        )

        for res in [res1, res2] {
            #expect(res.p2Value != nil)
            #expect(res.sampleSize != nil)
            #expect(res.nPosRanks != nil)
            #expect(res.nNegRanks != nil)
            #expect(res.nTies != nil)
            #expect(res.nZeroDiff != nil)
            #expect(res.sumNegRanks != nil)
            #expect(res.sumPosRanks != nil)
            #expect(res.meanNegRank != nil)
            #expect(res.meanPosRank != nil)
            #expect(res.zStat != nil)
            #expect(res.cohenD != nil)
            #expect(res.Tplus != nil)
            #expect(res.zeroHandling != nil)
            #expect(res.p2ValueExact != nil)
            #expect(!res.description.isEmpty)
        }
    }

    @Test("Runs test returns a populated result for both overloads")
    func runsTestResultAccessors() throws {
        let data: [Double] = [1, 2, 3, 4, 5, 6, 7, 8]
        let ex = makeExamine(data)

        let median: Randomness<Double>.RunsTestCuttingPoint<Double> = .median
        let userDefined: Randomness<Double>.RunsTestCuttingPoint<Double> = .userDefined(cuttingPoint: 4.5)
        let res1 = try #require(try Randomness<Double>.runsTest(data: ex, cuttingPoint: median, alpha: 0.05, withContinuityCorrection: false))
        let res2 = try #require(try Randomness<Double>.runsTest(data: data, cuttingPoint: userDefined, alpha: 0.05, withContinuityCorrection: true))

        for res in [res1, res2] {
            #expect(res.continuityCorrection != nil)
            #expect(res.alpha != nil)
            #expect(res.cp != nil)
            #expect(res.nAboveCP != nil)
            #expect(res.nBelowcp != nil)
            #expect(res.nRuns != nil)
            #expect(res.meanRuns != nil)
            #expect(res.varianceRuns != nil)
            #expect(res.zStatApprox != nil)
            #expect(res.criticalValue != nil)
            #expect(res.pValueApprox != nil)
            #expect(res.pTwoSidedExact != nil)
            #expect(res.pTwoSidedExactSymmetric != nil)
            #expect(res.pUpperExact != nil)
            #expect(res.pLowerExact != nil)
            #expect(res.signs != nil)
            #expect(res.randomness != nil)
            #expect(!res.description.isEmpty)
        }
    }

    @Test("Runs test supports all cutting point modes")
    func runsTestCuttingPointModes() throws {
        let dataWithMode: [Double] = [1, 1, 2, 1, 2,2,2,2,3,3,2, 3, 4, 5, 6, 7]
        let ex = makeExamine(dataWithMode)
        let alpha = 0.05

        let median: Randomness<Double>.RunsTestCuttingPoint<Double> = .median
        let mean: Randomness<Double>.RunsTestCuttingPoint<Double> = .mean
        let mode: Randomness<Double>.RunsTestCuttingPoint<Double> = .mode
        let user: Randomness<Double>.RunsTestCuttingPoint<Double> = .userDefined(cuttingPoint: 3.5)

        let r1 = try Randomness<Double>.runsTest(data: ex, cuttingPoint: median, alpha: alpha, withContinuityCorrection: false)
        let r2 = try Randomness<Double>.runsTest(data: ex, cuttingPoint: mean, alpha: alpha, withContinuityCorrection: false)
        let r3 = try Randomness<Double>.runsTest(data: ex, cuttingPoint: mode, alpha: alpha, withContinuityCorrection: false)
        let r4 = try Randomness<Double>.runsTest(data: ex, cuttingPoint: user, alpha: alpha, withContinuityCorrection: false)

        #expect(r1 != nil)
        #expect(r2 != nil)
        #expect(r3 != nil)
        #expect(r4 != nil)
    }

    @Test("Wald-Wolfowitz two-sample runs test returns a populated result for both overloads")
    func waldWolfowitzTwoSampleResultAccessors() throws {
        let a: [Double] = [1, 2, 3, 4, 5]
        let b: [Double] = [10, 11, 12, 13, 14]
        let ex1 = makeExamine(a)
        let ex2 = makeExamine(b)

        let res1 = try #require(try NonParametric<Double>.waldWolfowitzTwoSampleTest(set1: ex1, set2: ex2, alpha: 0.05, withContinuityCorrection: true))
        let res2 = try #require(try NonParametric<Double>.waldWolfowitzTwoSampleTest(set1: a, set2: b, alpha: 0.05, withContinuityCorrection: false))

        for res in [res1, res2] {
            #expect(res.nRuns != nil)
            #expect(res.zStat != nil)
            #expect(res.pValueExact != nil)
            #expect(res.pValueAsymp != nil)
            #expect(res.continuityCorrection != nil)
            #expect(res.alpha != nil)
            #expect(res.criticalValue != nil)
            #expect(res.rejectNull != nil)
            #expect(res.mean != nil)
            #expect(res.variance != nil)
            #expect(res.nTiesIntergroup != nil)
            #expect(res.nTiedCases != nil)
            #expect(res.sampleSize1 != nil)
            #expect(res.sampleSize2 != nil)
            #expect(!res.description.isEmpty)
        }
    }

    @Test("Sign test overloads return populated results")
    func signTestOverloadsReturnResults() throws {
        let set1 = makeExamine([1, 2, 3, 4], name: "s1")
        let set2 = makeExamine([2, 1, 5, 6], name: "s2")
        let oneSample = makeExamine([-1, 0, 2, -3, 4], name: "s")

        let pairedEx = try #require(try NonParametric<Double>.signTest(set1: set1, set2: set2))
        let pairedArray = try #require(try NonParametric<Double>.signTest(set1: [1, 2, 3, 4], set2: [2, 1, 5, 6]))
        let oneEx = try #require(try NonParametric<Double>.signTest(data: oneSample))
        let oneArray = try #require(try NonParametric<Double>.signTest(data: [-1, 0, 2, -3, 4]))

        for res in [pairedEx, pairedArray, oneEx, oneArray] {
            #expect(res.pValueExact != nil)
            #expect(res.pValueApprox != nil)
            #expect(res.nPosDiff != nil)
            #expect(res.nNegDiff != nil)
            #expect(res.nTies != nil)
            #expect(res.total != nil)
            #expect(res.zStat != nil)
            #expect(!res.description.isEmpty)
        }
    }

    @Test("Friedman test error cases are surfaced via the public error enum")
    func friedmanTestErrorCases() throws {
        #expect(throws: NonParametric<Double>.FriedmanTestError.self) {
            _ = try NonParametric<Double>.friedmanTest(data: [[Double]](), alpha: 0.05, tieCorrect: true)
        }

        #expect(throws: NonParametric<Double>.FriedmanTestError.self) {
            _ = try NonParametric<Double>.friedmanTest(data: [[Double]](repeating: [1.0], count: 1), alpha: 0.05, tieCorrect: true)
        }

        #expect(throws: NonParametric<Double>.FriedmanTestError.self) {
            _ = try NonParametric<Double>.friedmanTest(data: [[Double]]([[1.0, 2.0], [1.0]]), alpha: 0.05, tieCorrect: true)
        }

        #expect(throws: NonParametric<Double>.FriedmanTestError.self) {
            _ = try NonParametric<Double>.friedmanTest(data: [[1.0, 1.0], [2.0, 2.0]], alpha: 0.05, tieCorrect: true)
        }
    }

    @Test("Wilcoxon zero-handling enum exposes a description")
    func wilcoxonZeroHandlingDescriptions() {
        #expect(NonParametric<Double>.WilcoxonZeroHandling.discardZeroDifferences.description.contains("discardZeroDifferences"))
        #expect(NonParametric<Double>.WilcoxonZeroHandling.pratt.description == "Pratt")
    }

    @Test("T test entry points and one-way ANOVA populate all result accessors")
    func tTestAndAnovaResultAccessors() throws {
        let a: [Double] = [2, 2, 3, 3, 4]
        let b: [Double] = [10, 11, 12, 13, 14]
        let ex1 = makeExamine(a, name: "a")
        let ex2 = makeExamine(b, name: "b")

        let twoEx = try #require(try Parametric<Double>.twoSampleTTest(sample1: ex1, sample2: ex2, alpha: 0.05))
        let twoArray = try #require(try Parametric<Double>.twoSampleTTest(sample1: a, sample2: b, alpha: 0.05))

        for res in [twoEx, twoArray] {
            #expect(res.p1EQVAR != nil)
            #expect(res.p1UEQVAR != nil)
            #expect(res.p2EQVAR != nil)
            #expect(res.p2UEQVAR != nil)
            #expect(res.mean1 != nil)
            #expect(res.mean2 != nil)
            #expect(res.sampleSize1 != nil)
            #expect(res.sampleSize2 != nil)
            #expect(res.stdDev1 != nil)
            #expect(res.stdDev2 != nil)
            #expect(res.pooledStdDev != nil)
            #expect(res.pooledVariance != nil)
            #expect(res.differenceInMeans != nil)
            #expect(res.tEQVAR != nil)
            #expect(res.tUEQVAR != nil)
            #expect(res.LeveneP != nil)
            #expect(res.dfEQVAR != nil)
            #expect(res.dfUEQVAR != nil)
            #expect(res.mean1GTEmean2 != nil)
            #expect(res.mean1LTEmean2 != nil)
            #expect(res.mean1EQmean2 != nil)
            #expect(res.mean1UEQmean2 != nil)
            #expect(res.CVEQVAR != nil)
            #expect(res.CVUEQVAR != nil)
            #expect(res.rEQVAR != nil)
            #expect(res.rUEQVAR != nil)
            #expect(res.tWelch != nil)
            #expect(res.dfWelch != nil)
            #expect(res.p2Welch != nil)
            #expect(res.p1Welch != nil)
            #expect(res.variancesAreEqual != nil)
            #expect(!res.description.isEmpty)
        }

        let oneEx = try #require(try Parametric<Double>.oneSampleTTest(sample: ex1, mean: 0.0, alpha: 0.05))
        let oneArray = try #require(try Parametric<Double>.oneSampleTTest(sample: a, mean: 0.0, alpha: 0.05))

        for res in [oneEx, oneArray] {
            #expect(res.p1Value != nil)
            #expect(res.p2Value != nil)
            #expect(res.tStat != nil)
            #expect(res.cv90Pct != nil)
            #expect(res.cv95Pct != nil)
            #expect(res.cv99Pct != nil)
            #expect(res.cvAlpha != nil)
            #expect(res.mean != nil)
            #expect(res.mean0 != nil)
            #expect(res.difference != nil)
            #expect(res.sampleSize != nil)
            #expect(res.stdDev != nil)
            #expect(res.stdErr != nil)
            #expect(res.df != nil)
            #expect(res.meanGTEtestValue != nil)
            #expect(res.meanLTEtestValue != nil)
            #expect(res.meanEQtestValue != nil)
            #expect(!res.description.isEmpty)
        }

        let groups: [SSExamine<Double, Double>] = [
            makeExamine([1, 2, 3, 4], name: "g1"),
            makeExamine([2, 3, 4, 5], name: "g2"),
            makeExamine([10, 11, 12, 13], name: "g3")
        ]
        let anova = try #require(try Parametric<Double>.oneWayANOVA(data: groups, alpha: 0.05))
        #expect(anova.p2Value != nil)
        #expect(anova.FStatistic != nil)
        #expect(anova.pBartlett != nil)
        #expect(anova.alpha != nil)
        #expect(anova.meansEQUAL != nil)
        #expect(anova.cv != nil)
        #expect(anova.pLevene != nil)
        #expect(anova.SSTotal != nil)
        #expect(anova.SSError != nil)
        #expect(anova.SSTreatment != nil)
        #expect(anova.MSError != nil)
        #expect(anova.MSTreatment != nil)
        #expect(anova.dfError != nil)
        #expect(anova.dfTreatment != nil)
        #expect(anova.dfTotal != nil)
        #expect(!anova.description.isEmpty)
    }
}
