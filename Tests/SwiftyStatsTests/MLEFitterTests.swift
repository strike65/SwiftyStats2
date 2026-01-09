/// Tests that validate MLE fitters against analytic solutions and simulations.

import Testing
import SwiftyStats
import SwiftyBoost
import Foundation

@Suite("MLEFitter public surfaces")
struct MLEFitterTests {

    // MARK: - Analytic fitters

    @Test("Exponential analytic estimator matches reciprocal mean (seeded 0xFCEFF)")
    func exponentialFitMatchesAnalytic() async throws {
        let rate = 2.5
        let seed: UInt64 = 0x0000_0000_000F_CEFF
        let data = RNGSampler<Double>.randoms(n: 512, dist: .exponential(rate: rate), seed: seed)

        let result = MLEFitter<Double>.fitExponential(data: data)
        let lambdaHat = try #require(result.thetaHat.first)
        let manual = Double(data.count) / data.reduce(0, +)

        // Exact formula match in floating arithmetic; absolute tolerance OK.
        #expect(abs(lambdaHat - manual) <= 1e-12)
        #expect(result.converged)
    }

    @Test("Rayleigh analytic fit recovers sigma^2 ~= sum(x^2)/(2n) (seeded 0xABCDEF)")
    func rayleighFitRecoversSigmaSquared() async throws {
        let sigma = 1.2
        let seed: UInt64 = 0x00AB_CDEF
        let u = RNGSampler<Double>.randoms(n: 5000, dist: .uniform01, seed: seed)
        let data = u.map { sigma * sqrt(-2.0 * log(max($0, .leastNonzeroMagnitude))) }

        let result = MLEFitter<Double>.fitRayleigh(data: data)
        #expect(result.thetaHat.count == 1)
        let vHat = result.thetaHat[0] // sigma
        let trueV = sigma
        let relErr = abs(vHat - trueV) / trueV

        #expect(result.converged)
        #expect(relErr < 0.05)
    }

    @Test("Laplace analytic fit recovers location and scale (seeded 0xBEEFBEEF)")
    func laplaceFitRecoversParameters() async throws {
        let mu = -0.4
        let b = 0.8
        let seed: UInt64 = 0xBEEF_BEEF
        let u = RNGSampler<Double>.randoms(n: 6000, dist: .uniform01, seed: seed)
        let centered = u.map { $0 - 0.5 }
        let data = centered.map { t -> Double in
            mu - b * (t >= 0 ? 1 : -1) * log(1 - 2 * abs(t))
        }

        let result = MLEFitter<Double>.fitLaplace(data: data)
        #expect(result.thetaHat.count == 2)
        let muHat = result.thetaHat[0]
        let bHat = result.thetaHat[1]

        let muErr = abs(muHat - mu)
        let bRelErr = abs(bHat - b) / b

        #expect(result.converged)
        #expect(muErr < 0.05)
        #expect(bRelErr < 0.08)
    }

    @Test("Uniform analytic fit returns sample min/max (seeded 0x12345678)")
    func uniformFitReturnsOrderStats() async throws {
        let a = -2.0
        let b = 5.0
        let seed: UInt64 = 0x1234_5678
        let u = RNGSampler<Double>.randoms(n: 4000, dist: .uniform01, seed: seed)
        let data = u.map { a + (b - a) * $0 }
        let minX = try #require(data.min())
        let maxX = try #require(data.max())

        let result = MLEFitter<Double>.fitUniform(data: data)
        #expect(result.thetaHat.count == 2)
        let aHat = result.thetaHat[0]
        let bHat = result.thetaHat[1]

        // Exact equality with order stats for the sample
        #expect(aHat == minX)
        #expect(bHat == maxX)
        #expect(result.converged)
        // And bounds relative to true params
        #expect(aHat >= a)
        #expect(bHat <= b)
    }

    @Test("Wald (Inverse Gaussian) analytic fit recovers mean and shape (seeded 0xCAFEFEED)")
    func waldFitRecoversParameters() async throws {
        let mu = 1.0
        let lambda = 2.5
        let seedZ: UInt64 = 0xCAFE_FEED
        let seedU: UInt64 = 0xFACE_FEED

        let z = RNGSampler<Double>.randoms(n: 6000, dist: .gaussian(mean: 0, standardDeviation: 1), seed: seedZ)
        let u = RNGSampler<Double>.randoms(n: 6000, dist: .uniform01, seed: seedU)
        let data = zip(z, u).map { (zVal, uVal) -> Double in
            let y = zVal * zVal
            let muOverLam = mu / lambda
            let term = mu + (mu * mu * y) * 0.5 * (1.0 / lambda) - (muOverLam * 0.5) * sqrt(4 * mu * lambda * y + mu * mu * y * y)
            return (uVal <= mu / (mu + term)) ? term : (mu * mu / term)
        }

        let result = MLEFitter<Double>.fitWald(data: data)
        #expect(result.thetaHat.count == 2)
        let muHat = result.thetaHat[0]
        let lambdaHat = result.thetaHat[1]

        let muErr = abs(muHat - mu)
        let lambdaRelErr = abs(lambdaHat - lambda) / lambda

        #expect(result.converged)
        #expect(muErr < 0.03)
        #expect(lambdaRelErr < 0.15)
    }

    @Test("Arcsine analytic fit returns sample min/max (seeded 0xA11CE)")
    func arcsineFitReturnsOrderStats() async throws {
        let a = -1.5
        let b = 2.0
        let seed: UInt64 = 0x0A11_CE
        let data = RNGSampler<Double>.randoms(n: 5000, dist: .arcsine(a: a, b: b), seed: seed)
        let minX = try #require(data.min())
        let maxX = try #require(data.max())

        let result = MLEFitter<Double>.fitArcsine(data: data)
        #expect(result.thetaHat.count == 2)
        let aHat = result.thetaHat[0]
        let bHat = result.thetaHat[1]

        #expect(aHat == minX)
        #expect(bHat == maxX)
        #expect(result.converged)
        #expect(aHat >= a)
        #expect(bHat <= b)
    }

    @Test("Bernoulli analytic fit recovers success probability (seeded 0xB30B)")
    func bernoulliFitRecoversProbability() async throws {
        let p = 0.37
        let seed: UInt64 = 0x0B30_BEEF
        let data = RNGSampler<Double>.randoms(n: 4000, dist: .bernoulli(p: p), seed: seed)

        let result = MLEFitter<Double>.fitBernoulli(data: data)
        #expect(result.thetaHat.count == 1)
        let pHat = result.thetaHat[0]
        let relErr = abs(pHat - p) / p

        #expect(result.converged)
        #expect(relErr < 0.02)
    }

    @Test("Binomial numerical fit recovers success probability (seeded 0xB1B1)")
    func binomialFitRecoversProbability() async throws {
        let nTrials = 12
        let p = 0.58
        let seed: UInt64 = 0x0B1B_1E
        let data = RNGSampler<Double>.randoms(n: 5000, dist: .binomial(nTrials: nTrials, p: p), seed: seed)

        let result = MLEFitter<Double>.fitBinomial(data, n: nTrials)
        #expect(result.thetaHat.count == 1)
        let pHat = result.thetaHat[0]
        let relErr = abs(pHat - p) / p

        #expect(result.converged)
        #expect(relErr < 0.03)
    }

    @Test("Geometric numerical fit recovers success probability (seeded 0xC0DEC0DE)")
    func geometricFitRecoversProbability() async throws {
        let p = 0.42
        let seed: UInt64 = 0xC0DE_C0DE
        let data = RNGSampler<Double>.randoms(n: 5000, dist: .geometric(p: p), seed: seed)

        let result = MLEFitter<Double>.fitGeometric(data, optimizer: .lbfgs)
        #expect(result.thetaHat.count == 1)
        let pHat = result.thetaHat[0]
        let relErr = abs(pHat - p) / p

        #expect(result.converged)
        #expect(relErr < 0.03)
    }

    @Test("Poisson analytic fit recovers rate (seeded 0xAA55AA55)")
    func poissonFitRecoversLambda() async throws {
        let lambda = 3.6
        let seed: UInt64 = 0xAA55_AA55
        let data = RNGSampler<Double>.randoms(n: 4000, dist: .poisson(lambda: lambda), seed: seed)

        let result = MLEFitter<Double>.fitPoisson(data)
        #expect(result.thetaHat.count == 1)
        let lambdaHat = result.thetaHat[0]
        let relErr = abs(lambdaHat - lambda) / lambda

        #expect(result.converged)
        #expect(relErr < 0.02)
    }

    @Test("Negative Binomial numerical fit recovers r and p (seeded 0xF00D)")
    func negativeBinomialFitRecoversParameters() async throws {
        let rTrue = 5.0
        let pTrue = 0.4
        let seed: UInt64 = 0x0F00_D00D
        let data = RNGSampler<Double>.randoms(n: 5000, dist: .negativeBinomial(r: rTrue, p: pTrue), seed: seed)

        let result = MLEFitter<Double>.fitNegativeBinomial(data, optimizer: .lbfgs)
        #expect(result.thetaHat.count == 2)
        let rHat = result.thetaHat[0]
        let pHat = result.thetaHat[1]

        let rRel = abs(rHat - rTrue) / rTrue
        let pRel = abs(pHat - pTrue) / pTrue

        #expect(result.converged)
        #expect(rRel < 0.08)
        #expect(pRel < 0.05)
    }

    @Test("Pareto numerical fit recovers shape (seeded 0xBAD5EED)")
    func paretoFitRecoversParameters() async throws {
        let scale = 1.3
        let shape = 2.4
        let seed: UInt64 = 0x0BAD_5EED
        let data = RNGSampler<Double>.randoms(n: 6000, dist: .pareto(scale: scale, shape: shape), seed: seed)

        let result = MLEFitter<Double>.fitPareto(data, optimizer: .lbfgs)
        #expect(result.thetaHat.count == 1)
        let shapeHat = result.thetaHat[0]
        let rel = abs(shapeHat - shape) / shape

        #expect(result.converged)
        #expect(rel < 0.05)
    }

    @Test("Chi-Squared numerical fit recovers degrees of freedom (seeded 0xC1C1)")
    func chiSquaredFitRecoversDegreesOfFreedom() async throws {
        let df = 9.5
        let seed: UInt64 = 0x0C1C_1E
        let data = RNGSampler<Double>.randoms(n: 4000, dist: .chiSquare(df: df), seed: seed)

        let opts = centralOptions(optimizer: .lbfgs, seed: seed)
        let result = MLEFitter<Double>.fitChiSquared(data, optimizer: .lbfgs, options: opts)
        #expect(result.thetaHat.count == 1)
        let dfHat = result.thetaHat[0]
        let relErr = abs(dfHat - df) / df

        #expect(result.converged)
        #expect(relErr < 0.08)
    }

    @Test("Inverse Chi-Squared numerical fit recovers parameters (seeded 0x1C1C)")
    func inverseChiSquaredFitRecoversParameters() async throws {
        let nu = 15.0
        let tau2 = 1.2
        let seed: UInt64 = 0x1C1C_1C1C
        let chi = RNGSampler<Double>.randoms(n: 4000, dist: .chiSquare(df: nu), seed: seed)
        let data = chi.map { (nu * tau2) / $0 }

        let opts = centralOptions(optimizer: .lbfgs, seed: seed)
        let result = MLEFitter<Double>.fitInverseChiSquared(data, optimizer: .lbfgs, options: opts)
        #expect(result.thetaHat.count == 2)
        let nuHat = result.thetaHat[0]
        let tauHat = result.thetaHat[1]

        let nuRel = abs(nuHat - nu) / nu
        let tauRel = abs(tauHat - tau2) / tau2

        #expect(result.converged)
        #expect(nuRel < 0.08)
        #expect(tauRel < 0.08)
    }

    @Test("Inverse Normal numerical fit recovers mean and shape (seeded 0x1N1N)")
    func inverseNormalFitRecoversParameters() async throws {
        let mu = 0.8
        let lambda = 2.4
        let seed: UInt64 = 0x1A1A_0001
        let data = RNGSampler<Double>.randoms(n: 5000, dist: .inverseNormal(mu: mu, sigma: lambda), seed: seed)

        let opts = centralOptions(optimizer: .lbfgs, seed: seed)
        let result = MLEFitter<Double>.fitInverseNormal(data, optimizer: .lbfgs, options: opts)
        #expect(result.thetaHat.count == 2)
        let muHat = result.thetaHat[0]
        let lambdaHat = result.thetaHat[1]

        #expect(result.converged)
        #expect(abs(muHat - mu) < 0.03)
        #expect(abs(lambdaHat - lambda) / lambda < 0.12)
    }

    @Test("Log-Normal numerical fit recovers mu and sigma (seeded 0x10G10G)")
    func logNormalFitRecoversParameters() async throws {
        let mu = -0.2
        let sigma = 0.45
        let seed: UInt64 = 0x10A1_0A10
        let data = RNGSampler<Double>.randoms(n: 5000, dist: .lognormal(lmean: mu, lsigma: sigma), seed: seed)

        let result = MLEFitter<Double>.fitLogNormal(data)
        #expect(result.thetaHat.count == 2)
        let muHat = result.thetaHat[0]
        let sigmaHat = result.thetaHat[1]

        #expect(result.converged)
        #expect(abs(muHat - mu) < 0.03)
        #expect(abs(sigmaHat - sigma) / sigma < 0.05)
    }

    @Test("Landau numerical fit recovers location and scale (seeded 0x1A2B)")
    func landauFitRecoversParameters() async throws {
        let mu = 0.1
        let c = 0.6
        let seed: UInt64 = 0x1A2B_3C4D
        let data = RNGSampler<Double>.randoms(n: 6000, dist: .landau(location: mu, scale: c), seed: seed)

        var opts = centralOptions(optimizer: .lbfgs, seed: seed)
        opts.multiStartCount = 6
        opts.randomRestartCount = 2

        let result = MLEFitter<Double>.fitLandau(data, optimizer: .lbfgs, options: opts)
        #expect(result.thetaHat.count == 2)
        let muHat = result.thetaHat[0]
        let cHat = result.thetaHat[1]

        #expect(result.converged)
        #expect(abs(muHat - mu) < 0.15)
        #expect(abs(cHat - c) / c < 0.25)
    }

    @Test("Map Airy numerical fit recovers location and scale (seeded 0xAA11)")
    func mapAiryFitRecoversParameters() async throws {
        let mu = 0.25
        let sigma = 0.7
        let seed: UInt64 = 0xAA11_5511
        let data = RNGSampler<Double>.randoms(n: 2500, dist: .mapAiry(location: mu, scale: sigma), seed: seed)

        var opts = centralOptions(optimizer: .lbfgs, seed: seed)
        opts.multiStartCount = 4
        opts.randomRestartCount = 1
        opts.warmStartTheta = [mu, sigma]

        let result = MLEFitter<Double>.fitMapAiry(data, optimizer: .lbfgs, options: opts)
        #expect(result.thetaHat.count == 2)
        let muHat = result.thetaHat[0]
        let sigmaHat = result.thetaHat[1]

        #expect(result.converged)
        #expect(abs(muHat - mu) < 0.5)
        #expect(abs(sigmaHat - sigma) / sigma < 0.2)
    }

    @Test("SAS 0.5 numerical fit recovers location and scale (seeded 0x55AA)")
    func sasPoint5FitRecoversParameters() async throws {
        let mu = 1.2
        let sigma = 0.6
        let dist = try SwiftyBoost.Distribution.SASPoint5<Double>(location: mu, scale: sigma)
        let n = 64
        var data = [Double]()
        data.reserveCapacity(n)
        for k in 0..<n {
            let u = (Double(k) + 0.5) / Double(n + 1)
            data.append(try dist.quantile(u))
        }

        var opts = centralOptions(optimizer: .lbfgs)
        opts.warmStartTheta = [mu, sigma]

        let result = MLEFitter<Double>.fitSASPoint5(data: data, optimizer: .lbfgs, options: opts)
        #expect(result.thetaHat.count == 2)
        let muHat = result.thetaHat[0]
        let sigmaHat = result.thetaHat[1]

        #expect(result.converged)
        #expect(abs(muHat - mu) < 0.05)
        #expect(abs(sigmaHat - sigma) / sigma < 0.05)
    }

    @Test("Triangular numerical fit recovers lower, mode, upper (seeded 0x1357)")
    func triangularFitRecoversParameters() async throws {
        let lower = -1.0
        let upper = 2.5
        let mode = 0.6
        let dist = try SwiftyBoost.Distribution.Triangular<Double>(lower: lower, mode: mode, upper: upper)
        let n = 64
        var data = [Double]()
        data.reserveCapacity(n)
        for k in 0..<n {
            let u = (Double(k) + 0.5) / Double(n + 1)
            data.append(try dist.quantile(u))
        }

        var opts = centralOptions(optimizer: .lbfgs)
        opts.warmStartTheta = [mode, lower, upper]
        let result = MLEFitter<Double>.fitTriangular(data, optimizer: .lbfgs, options: opts)
        #expect(result.thetaHat.count == 3)
        let modeHat = result.thetaHat[0]
        let lowerHat = result.thetaHat[1]
        let upperHat = result.thetaHat[2]

        #expect(result.converged)
        #expect(abs(modeHat - mode) < 0.08)
        #expect(abs(lowerHat - lower) < 0.08)
        #expect(abs(upperHat - upper) < 0.08)
    }


    // MARK: - Numerical fitters (central families)

    // Options helper for central continuous families that benefit from stronger optimization
    private func centralOptions(optimizer: OptimizerKind = .lbfgs, seed: UInt64? = nil) -> MLEOptimizationOpts<Double> {
        var opt = MLEOptimizationOpts<Double>()
        opt.optimizer = optimizer
        opt.computeCovariance = true

        // More robust multi-starts for tougher surfaces
        opt.multiStartCount = 10
        opt.multiStartDesign = .lhs
        opt.randomRestartCount = 3
        opt.rngSeed = seed
        opt.enableParallelStarts = false
        opt.parallelTaskPriority = TaskPriority.low

        // Relative steps in θ for better scaling
        opt.initialStepStrategy = .relativeTheta(0.35)

        // Tuned finite differences and Hessian steps
        opt.gradStep = 5e-6
        opt.hessianStep = 1e-3

        // Slightly relaxed tolerances
        opt.relTolLogLik = 1e-7
        opt.tolSimplexNM = 1e-6
        opt.relTolFSpreadNM = 1e-8

        opt.diagnosticsEnabled = true
        return opt
    }

    @Test("Beta numerical fit recovers alpha and beta (seeded 0xBADA55)")
    func betaFitRecoversParameters() async throws {
        let a = 2.2
        let b = 0.6
        let seed: UInt64 = 0x00BA_DA55
        let data = RNGSampler<Double>.randoms(n: 5000, dist: .beta(a: a, b: b), seed: seed)

        let result = MLEFitter<Double>.fitBeta(data)
        #expect(result.thetaHat.count == 2)

        let aHat = result.thetaHat[0]
        let bHat = result.thetaHat[1]
        let relA = abs(aHat - a) / a
        let relB = abs(bHat - b) / b

        #expect(result.converged)
        #expect(relA < 0.08)
        #expect(relB < 0.08)
    }

    @Test("Gamma numerical fit recovers shape and scale (seeded 0xFACEB00C)")
    func gammaFitRecoversParameters() async throws {
        let shape = 2.2
        let scale = 1.35
        let seed: UInt64 = 0xFACE_B00C
        let data = RNGSampler<Double>.randoms(n: 4000, dist: .gamma(shape: shape, scale: scale), seed: seed)

        let result = MLEFitter<Double>.fitGamma(data)
        #expect(result.thetaHat.count == 2)

        let shapeHat = result.thetaHat[0]
        let scaleHat = result.thetaHat[1]
        let relShapeErr = abs(shapeHat - shape) / shape
        let relScaleErr = abs(scaleHat - scale) / scale

        #expect(result.converged)
        #expect(relShapeErr < 0.08)
        #expect(relScaleErr < 0.08)
    }

    @Test("Weibull numerical fit recovers shape and scale (seeded 0xDEADC0DE)")
    func weibullFitRecoversParameters() async throws {
        let k = 1.4
        let lambda = 0.9
        let seed: UInt64 = 0xDEAD_C0DE
        // Increase n to reduce variance; Weibull MLE can be noisy for shape/scale
        let data = RNGSampler<Double>.randoms(n: 10_000, dist: .weibull(k: k, lambda: lambda), seed: seed)

        let opts = centralOptions(seed: seed)
        let result = MLEFitter<Double>.fitWeibull(data, optimizer: .lbfgs, options: opts)
        #expect(result.thetaHat.count == 2)
        let kHat = result.thetaHat[0]
        let lamHat = result.thetaHat[1]
        let relK = abs(kHat - k) / k
        let relLam = abs(lamHat - lambda) / lambda

        #expect(result.converged)
        #expect(relK < 0.20)
        #expect(relLam < 0.20)
    }

    @Test("Cauchy numerical fit recovers location and scale (seeded 0xC0C0A)")
    func cauchyFitRecoversParameters() async throws {
        let mu = 0.3
        let gamma = 1.1
        let seed: UInt64 = 0x0C0C_0A
        let data = RNGSampler<Double>.randoms(n: 6000, dist: .cauchy(mu: mu, gamma: gamma), seed: seed)

        let result = MLEFitter<Double>.fitCauchy(data)
        #expect(result.thetaHat.count == 2)
        let muHat = result.thetaHat[0]
        let gamHat = result.thetaHat[1]

        let muErr = abs(muHat - mu)
        let gamRel = abs(gamHat - gamma) / gamma

        #expect(result.converged)
        #expect(muErr < 0.05)
        #expect(gamRel < 0.05)
    }

    @Test("Logistic numerical fit recovers location and scale (seeded 0xC0FFEE)")
    func logisticFitRecoversParameters() async throws {
        let location = -0.75
        let scale = 0.65
        let seed: UInt64 = 0x00C0_FFEE
        let data = RNGSampler<Double>.randoms(n: 4500, dist: .logistic(mu: location, s: scale), seed: seed)

        let result = MLEFitter<Double>.fitLogistic(data)
        #expect(result.thetaHat.count == 2)

        let locationHat = result.thetaHat[0]
        let scaleHat = result.thetaHat[1]
        let locErr = abs(locationHat - location)
        let scaleErr = abs(scaleHat - scale) / scale

        #expect(result.converged)
        #expect(locErr < 0.05)
        #expect(scaleErr < 0.05)
    }

    @Test("Extreme Value (Gumbel) numerical fit recovers location and scale (seeded 0xE17E)")
    func extremeValueFitRecoversParameters() async throws {
        let mu = 3.0
        let beta = 5.0
        let seed: UInt64 = 0x0E17_E
        let data = RNGSampler<Double>.randoms(n: 5000, dist: .gumbel(mu: mu, beta: beta), seed: seed)

        let result = MLEFitter<Double>.fitExtremeValue(data)
        #expect(result.thetaHat.count == 2)
        let muHat = result.thetaHat[0]
        let betaHat = result.thetaHat[1]

        let muErr = abs(muHat - mu)
        let betaRel = abs(betaHat - beta) / beta

        #expect(result.converged)
        #expect(muErr < 0.10)
        #expect(betaRel < 0.10)
    }

    @Test("Student's t (standard) numerical fit recovers nu (seeded 0x0707_E7E7)")
    func studentsTStandardFitRecoversNu() async throws {
        let nu = 7.5
        let seed: UInt64 = 0x0707_E7E7
        let data = RNGSampler<Double>.randoms(n: 6000, dist: .studentT(mu: 0, sigma: 1, nu: nu), seed: seed)

        let result = MLEFitter<Double>.fitStudentsT(data)
        #expect(result.thetaHat.count == 1)
        let nuHat = result.thetaHat[0]
        let rel = abs(nuHat - nu) / nu

        #expect(result.converged)
        #expect(rel < 0.20)
    }

    @Test("Fisher F numerical fit recovers d1 and d2 (seeded 0xF17F17)")
    func fisherFFitRecoversParameters() async throws {
        let d1 = 10.0
        let d2 = 14.0
        let seed: UInt64 = 0x0F17_F17
        let data = RNGSampler<Double>.randoms(n: 6000, dist: .f(d1: d1, d2: d2), seed: seed)

        var opts = centralOptions(seed: seed)
        opts.multiStartCount = 6
        opts.randomRestartCount = 2
        opts.multiStartDesign = .sobol
        opts.initialStepStrategy = .relativeTheta(Double(0.25))
        opts.optimizer = .lbfgs

        let result = MLEFitter<Double>.fitFisherF(data, optimizer: .lbfgs, options: opts)
        #expect(result.thetaHat.count == 2)
        let d1Hat = result.thetaHat[0]
        let d2Hat = result.thetaHat[1]
        let rel1 = abs(d1Hat - d1) / d1
        let rel2 = abs(d2Hat - d2) / d2

        #expect(result.converged)
        #expect(rel1 < 0.08)
        #expect(rel2 < 0.08)
    }

    @Test("Inverse-Gamma numerical fit recovers alpha and beta (seeded 0x01AF_E0)")
    func inverseGammaFitRecoversParameters() async throws {
        let alpha = 3.0
        let beta = 2.0
        let seed: UInt64 = 0x01AF_E0
        // Increase n for better stability
        let y = RNGSampler<Double>.randoms(n: 10_000, dist: .gamma(shape: alpha, scale: 1), seed: seed)
        let data = y.map { beta / $0 }

        let opts = centralOptions(seed: seed)
        let result = MLEFitter<Double>.fitInverseGamma(data, optimizer: .lbfgs, options: opts)
        #expect(result.thetaHat.count == 2)
        let aHat = result.thetaHat[0]
        let bHat = result.thetaHat[1]
        let relA = abs(aHat - alpha) / alpha
        let relB = abs(bHat - beta) / beta

        #expect(result.converged)
        #expect(relA < 0.20)
        #expect(relB < 0.20)
    }

    @Test("Skew-Normal numerical fit recovers location, scale, shape (seeded 0x05AE_E5)")
    func skewNormalFitRecoversParameters() async throws {
        let xi = -0.2
        let omega = 1.1
        let alpha = 3.0
        let seed: UInt64 = 0x05AE_E5
        // Reduce n and lighten options to avoid timeouts while keeping accuracy reasonable.
        let data = RNGSampler<Double>.randoms(n: 1001, dist: .skewNormal(xi: xi, omega: omega, alpha: alpha), seed: seed)

        var opts = centralOptions(seed: seed)
        // Lighten the optimization load for this test
        opts.multiStartCount = 4
        opts.randomRestartCount = 1
        opts.diagnosticsEnabled = false
        // Slightly larger grad step stabilizes finite differences on this surface
        opts.gradStep = 1e-6

        let result = MLEFitter<Double>.fitSkewNormal(data, optimizer: .lbfgs, options: opts)
        #expect(result.thetaHat.count == 3)
        let xiHat = result.thetaHat[0]
        let omegaHat = result.thetaHat[1]
        let alphaHat = result.thetaHat[2]

        let xiErr = abs(xiHat - xi)
        let omegaRel = abs(omegaHat - omega) / omega
        let alphaAbs = abs(alphaHat - alpha)

        #expect(result.converged)
        #expect(xiErr < 0.1)
        #expect(omegaRel < 0.1)
        #expect(alphaAbs < 0.5)
    }

    // MARK: - Options helper for noncentral families

    private func ncOptions(optimizer: OptimizerKind = .nelderMead, seed: UInt64? = nil) -> MLEOptimizationOpts<Double> {
        var opt = MLEOptimizationOpts<Double>()
        opt.optimizer = optimizer
        opt.computeCovariance = true

        // Multi-start to avoid poor local optima
        opt.multiStartCount = 8
        opt.multiStartDesign = .lhs
        opt.randomRestartCount = 2
        opt.rngSeed = seed

        // Relative simplex steps in θ to stabilize scaling
        opt.initialStepStrategy = .relativeTheta(0.25)

        // Gentler finite differences and Hessian steps for rough likelihoods
        opt.gradStep = 1e-5
        opt.hessianStep = 5e-4

        // Slightly relaxed tolerances for noisy/flat regions
        opt.relTolLogLik = 1e-7
        opt.tolSimplexNM = 1e-6
        opt.relTolFSpreadNM = 1e-8
        opt.stallIterationsNM = 300
        opt.stallRelTolNM = 1e-9

        // Keep diagnostics on for visibility in failures
        opt.diagnosticsEnabled = true

        return opt
    }

    // MARK: - Numerical fitters (noncentral and specialized)

    @Test("Noncentral Chi-Squared numerical fit recovers k and lambda (seeded 0x0AC2_C2)")
    func ncChiSquaredFitRecoversParameters() async throws {
        let k = 6.0
        let lambda = 5.0
        let seed: UInt64 = 0x0AC2_C2
        let data = RNGSampler<Double>.randoms(n: 1000, dist: .noncentralChiSquare(k: k, lambda: lambda), seed: seed)

        // Estimate λ with known k; use tuned options
        let opts = ncOptions(optimizer: .lbfgs, seed: seed)
        let result = MLEFitter<Double>.fitNCChiSquared(data, degreesOfFreedom: k, optimizer: .lbfgs, options: opts)

        #expect(result.thetaHat.count == 1)
        let lamHat = result.thetaHat[0]
        let relLam = abs(lamHat - lambda) / lambda

        #expect(result.converged)
        #expect(relLam < 0.25)
    }

    @Test("Noncentral Student's t numerical fit recovers nu and delta (seeded 0x00AC_E7)")
    func ncStudentsTFitRecoversParameters() async throws {
        let nu = 41.0
        let delta = 2.5
        let seed: UInt64 = 0x00AC_E7
        let data = RNGSampler<Double>.randoms(n: 10_000, dist: .noncentralT(mu: 0, sigma: 1, nu: nu, delta: delta), seed: seed)

        // Warm start δ from sample mean (for ν not too small, E[T] ≈ δ)
        let mean = data.reduce(0.0, +) / Double(data.count)

        var opts = ncOptions(optimizer: .lbfgs, seed: seed)
        opts.warmStartTheta = [mean]
        // Slightly larger grad step and larger relative theta step help leave flat regions
        opts.gradStep = 2e-5
        opts.initialStepStrategy = .relativeTheta(0.5)

        let result = MLEFitter<Double>.fitNCStudentsT(data, degreesOfFreedom: nu, optimizer: .lbfgs, options: opts)

        #expect(result.thetaHat.count == 1)
        let deltaHat = result.thetaHat[0]
        let absDelta = abs(deltaHat - delta)

        #expect(result.converged)
        #expect(absDelta < 0.30)
    }

    @Test("Noncentral Fisher F numerical fit recovers d1, d2, lambda (seeded 0x0ACF_001/002)")
    func ncFisherFFitRecoversParameters() async throws {
        let d1 = 8.0, d2 = 16.0, lambda = 3.5
        let seedX: UInt64 = 0x0ACF_001
        let seedY: UInt64 = 0x0ACF_002
        let x1 = RNGSampler<Double>.randoms(n: 9000, dist: .noncentralChiSquare(k: d1, lambda: lambda), seed: seedX)
        let x2 = RNGSampler<Double>.randoms(n: 9000, dist: .chiSquare(df: d2), seed: seedY)
        let data = zip(x1, x2).map { (a, b) in (a / d1) / (b / d2) }

        var opts = ncOptions(optimizer: .lbfgs, seed: seedX ^ seedY)

        // Warm start using mean-based λ heuristic: λ = (d1(d2−2)/d2)·E[F] − d1
        let ex = try SSExamine<Double, Double>(using: data, levelOfMeasurement: .ratio, name: nil, characterSet: nil)
        let m = max(ex.arithmeticMean ?? 1.0, 1.000001)
        let lamWarm = max(((d1 * (d2 - 2.0)) / d2) * m - d1, 1e-6)
        opts.warmStartTheta = [d1, d2, lamWarm]

        let result = MLEFitter<Double>.fitNCFisherF(data, df1: d1, df2: d2, optimizer: .lbfgs, options: opts)
        #expect(result.converged)

        let lamHat = result.thetaHat.last ?? .nan
        let relLam = abs(lamHat - lambda) / lambda
        #expect(relLam < 0.30)
    }

    @Test("Noncentral Beta numerical fit recovers alpha, beta; lambda is checked weakly (seeded 0x0ACB_001/002)")
    func ncBetaFitRecoversParameters() async throws {
        let a = 3.0, b = 2.5, lambda = 1.2
        // X ~ χ²(2a, 2λ), Y ~ χ²(2b) => X/(X+Y) ~ NCBeta(α=a, β=b, λ=lambda)
        let seedX: UInt64 = 0x0ACB_001
        let seedY: UInt64 = 0x0ACB_002

        let X = RNGSampler<Double>.randoms(n: 10_000, dist: .noncentralChiSquare(k: 2 * a, lambda: 2 * lambda), seed: seedX)
        let Y = RNGSampler<Double>.randoms(n: 10_000, dist: .chiSquare(df: 2 * b), seed: seedY)
        let data = zip(X, Y).map { $0 / ($0 + $1) }

        var opts = ncOptions(optimizer: .lbfgs, seed: seedX ^ seedY)
        // Slightly larger grad step helps with very flat regions in λ
        opts.gradStep = 2e-5

        let result = MLEFitter<Double>.fitNCBeta(data, optimizer: .lbfgs, options: opts)
        #expect(result.thetaHat.count == 3)
        let aHat = result.thetaHat[0]
        let bHat = result.thetaHat[1]
        let lamHat = result.thetaHat[2]

        let relA = abs(aHat - a) / a
        let relB = abs(bHat - b) / b

        #expect(result.converged)
        #expect(relA < 0.25)
        #expect(relB < 0.25)
        // λ is weakly identified; ensure it's finite and nonnegative, and allow a very loose bound.
        #expect(lamHat.isFinite && lamHat >= 0)
        let relLam = abs(lamHat - lambda) / max(lambda, 1e-9)
        #expect(relLam < 1.50)
    }

    @Test("Holtsmark numerical fit recovers location and scale (seeded 0x0A17_0005)")
    func holtsmarkFitRecoversParameters() async throws {
        let mu = -1.5
        let c = 1.4
        let seed: UInt64 = 0x0A17_0005
        // Heavier tails: increase n for stability
        let data = RNGSampler<Double>.randoms(n: 10_000, dist: .holtsmark(mu: mu, c: c), seed: seed)

        let opts = centralOptions(seed: seed)
        let result = MLEFitter<Double>.fitHoltsmark(data, optimizer: .lbfgs, options: opts)
        #expect(result.thetaHat.count == 2)
        let muHat = result.thetaHat[0]
        let cHat = result.thetaHat[1]

        let muErr = abs(muHat - mu)
        let cRel = abs(cHat - c) / c

        #expect(result.converged)
        #expect(muErr < 0.1)
        #expect(cRel < 0.1)
    }

    // MARK: - Hyperexponential fitter tests

    // Simple local helper: sample an exponential(rate) via inverse CDF using a provided U(0,1) variate.
    private func inverseExp(rate: Double, u: Double) -> Double {
        let uu = max(u, .leastNonzeroMagnitude)
        return -log(uu) / rate
    }

    // Generate hyperexponential samples by mixture: choose phase by probabilities, then draw exponential with that rate.
    private func sampleHyperexponential(n: Int, probs: [Double], rates: [Double], seed: UInt64) -> [Double] {
        precondition(n > 0)
        precondition(probs.count == rates.count && !probs.isEmpty)
        let k = probs.count
        var cum = [Double](repeating: 0, count: k)
        var s: Double = 0
        for i in 0..<k { s += probs[i]; cum[i] = s }
        for i in 0..<k { cum[i] /= s }

        // Use our RNGSampler uniform stream for determinism.
        let u1 = RNGSampler<Double>.randoms(n: n, dist: .uniform01, seed: seed ^ 0x9E37)
        let u2 = RNGSampler<Double>.randoms(n: n, dist: .uniform01, seed: seed ^ 0xA5A5_5A5A)

        var out = [Double](repeating: 0, count: n)
        for i in 0..<n {
            let u = u1[i]
            var idx = 0
            while idx < k-1 && u > cum[idx] { idx += 1 }
            out[i] = inverseExp(rate: rates[idx], u: u2[i])
        }
        return out
    }

    // Compute best permutation match error between estimated (p, λ) and true (p, λ) to handle label switching.
    private func bestPermutationErrors(phat: [Double], lhat: [Double], p: [Double], l: [Double]) -> (wL1: Double, rateRelMax: Double) {
        let k = p.count
        var used = [Bool](repeating: false, count: k)
        var wL1Best = Double.infinity
        var rRelBest = Double.infinity

        func backtrack(_ perm: inout [Int], _ depth: Int) {
            if depth == k {
                var wErr = 0.0
                var rErrMax = 0.0
                for i in 0..<k {
                    wErr += abs(phat[i] - p[perm[i]])
                    rErrMax = max(rErrMax, abs(lhat[i] - l[perm[i]]) / max(l[perm[i]], 1e-12))
                }
                wL1Best = min(wL1Best, wErr)
                rRelBest = min(rRelBest, rErrMax)
                return
            }
            for j in 0..<k where !used[j] {
                used[j] = true
                perm.append(j)
                backtrack(&perm, depth + 1)
                perm.removeLast()
                used[j] = false
            }
        }
        var perm: [Int] = []
        backtrack(&perm, 0)
        return (wL1Best, rRelBest)
    }

    @Test("Hyperexponential (2-phase) numerical fit recovers probabilities and rates")
    func hyperexponential2PhaseFitRecoversParameters() async throws {
        let probs = [0.1, 0.9]
        let rates = [1.5, 5.0]
        let seed: UInt64 = 0x00AB_C012

        let n = 20_000
        let data = sampleHyperexponential(n: n, probs: probs, rates: rates, seed: seed)

        var opts = centralOptions(optimizer: .lbfgs, seed: seed)
        // Tougher surface: increase multi-starts slightly
        opts.multiStartCount = 32
        opts.randomRestartCount = 5
        opts.gradStep = 1e-6
        opts.hessianStep = 5e-4
        opts.initialStepStrategy = .relativeTheta(0.35)

        let result = MLEFitter<Double>.fitHyperexponential(data, phaseCount: 2, optimizer: .lbfgs, options: opts)

        #expect(result.converged)
        #expect(result.thetaHat.count == 4)

        let phat = Array(result.thetaHat[0..<2])
        let lhat = Array(result.thetaHat[2..<4])

        // Probabilities should sum to ~1 and be in (0,1); rates positive
        let sumP = phat.reduce(0, +)
        #expect(abs(sumP - 1.0) < 1e-6)
        #expect(phat.allSatisfy { $0 > 0 && $0 < 1 })
        #expect(lhat.allSatisfy { $0 > 0 })

        // Handle label switching by taking best permutation match
        let (wL1, rRelMax) = bestPermutationErrors(phat: phat, lhat: lhat, p: probs, l: rates)

        // Mixture models are harder; use reasonable tolerances
        #expect(wL1 < 0.2)
        #expect(rRelMax < 0.35)
    }

}
