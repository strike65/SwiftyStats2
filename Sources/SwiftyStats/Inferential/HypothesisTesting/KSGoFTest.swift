/// Kolmogorov–Smirnov test utilities, distributions, and bootstrap helpers.

import SwiftyStatsPrelude

extension Inferential.HypothesisTesting.GoodnessOfFitTests {
    /// Performs a one-sample KS test by fitting the target distribution and bootstrapping D.
    ///
    /// - Parameters:
    ///   - data: Sample observations.
    ///   - testTarget: Distribution to fit/compare against.
    ///   - bootstrapCount: Max. count of samples to draw
    ///   - seed: Optional base seed for reproducible bootstrap sampling.
    /// - Returns: `KSTestResult` containing θ̂, D statistics, and bootstrap p-value.
    public static func ksTestOneSample(data: [T], testTarget: GOFTestTarget, boostrapCount bsCount: Int = 1000, seed: UInt64? = nil) ->  KSTestResult {
        let samplerSeed = makeBootstrapSamplerSeed(seed: seed)
        switch testTarget {
        case .arcsine:
            func cdf(from theta: [T]) -> (T) -> T {
                let dist = try? SwiftyBoost.Distribution.Arcsine<T>(minX: theta[0], maxX: theta[1])
                return { x in
                    guard let d = dist else { return .nan }
                    do {
                        return try d.cdf(x)
                    } catch {
                        return .nan
                    }
                }
            }
            
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .arcsine(a: theta[0], b: theta[1]), seed: seed)
            }
            
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .arcsine,
                    fit: { x in MLEFitter<T>.fitArcsine(data: x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return KSTestResult(
                    distribution: GOFTestTarget.arcsine.description,
                    n: 0,
                    d: .nan,
                    dPlus: .nan,
                    dMinus: .nan,
                    pBootstrap: .nan,
                    b: 0,
                    thetaHat: []
                )
            }
        case .bernoulli:
            func cdf(from theta: [T]) -> (T) -> T {
                let p = theta[0]
                let dist = try? SwiftyBoost.Distribution.Bernoulli<T>(p: p)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .bernoulli(p: theta[0]), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .bernoulli,
                    fit: { x in MLEFitter<T>.fitBernoulli(data: x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.bernoulli.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .beta:
            func cdf(from theta: [T]) -> (T) -> T {
                let dist = try? SwiftyBoost.Distribution.Beta<T>(alpha: theta[0], beta: theta[1])
                return { x in
                    guard let d = dist else { return .nan }
                    do {
                        return try d.cdf(x)
                    } catch {
                        return .nan
                    }
                }
            }
            
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .beta(a: theta[0], b: theta[1]), seed: seed)
            }
            
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .beta,
                    fit: { x in MLEFitter<T>.fitBeta(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return KSTestResult(
                    distribution: GOFTestTarget.beta.description,
                    n: 0,
                    d: .nan,
                    dPlus: .nan,
                    dMinus: .nan,
                    pBootstrap: .nan,
                    b: 0,
                    thetaHat: []
                )
            }
        case .binomial:
            func estimateTrials(_ x: [T]) -> Int {
                let maxX = max(0, x.map { Int($0.rounded(.toNearestOrAwayFromZero)) }.max() ?? 0)
                let ex = SSExamine<T, T>(usingArray: x, name: nil, characterSet: nil)
                let mean = ex.arithmeticMean ?? T.zero
                let variance = ex.populationVariance ?? T.zero
                
                var momGuess: Int? = nil
                if mean > T.zero && variance >= T.zero {
                    let pMom = T.one - variance / mean
                    if pMom > T.zero && pMom < T.one {
                        let nMom = mean / pMom
                        if nMom.isFinite && nMom > T.zero {
                            momGuess = max(maxX, Int(nMom.rounded(.toNearestOrAwayFromZero)))
                        }
                    }
                }
                
                func logLik(n: Int) -> T? {
                    if n < maxX { return nil }
                    let nT = T(n)
                    let p = mean / nT
                    if !(p > T.zero && p < T.one && nT.isFinite) { return nil }
                    guard let dist = try? SwiftyBoost.Distribution.Binomial<T>(numberOfTrials: nT, probabibilityOfSuccess: p) else {
                        return nil
                    }
                    var ll: T = 0
                    for value in x {
                        do {
                            ll += try dist.logPdf(value)
                        } catch {
                            return nil
                        }
                    }
                    return ll.isFinite ? ll : nil
                }
                
                let minSpan = 10
                var lower = maxX
                var upper = max(lower + minSpan, momGuess ?? (lower + minSpan))
                var bestN = maxX
                var bestLL = -T.infinity
                var expandCount = 0
                let maxExpand = 4
                
                while true {
                    if lower <= upper {
                        for n in lower...upper {
                            if let ll = logLik(n: n), ll > bestLL {
                                bestLL = ll
                                bestN = n
                            }
                        }
                    }
                    if bestN == upper && expandCount < maxExpand {
                        expandCount += 1
                        lower = upper + 1
                        let grow = max(minSpan, upper / 2)
                        upper = upper + grow
                        continue
                    }
                    break
                }
                return max(bestN, maxX)
            }
            func cdf(from theta: [T]) -> (T) -> T {
                let nTrials = Int(theta[0].rounded())
                let p = theta[1]
                let dist = try? SwiftyBoost.Distribution.Binomial<T>(numberOfTrials: T(nTrials), probabibilityOfSuccess: p)
                return { x in
                    guard let d = dist else { return .nan }
                    do {
                        return try d.cdf(x)
                    } catch {
                        return .nan
                    }
                }
            }
            
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let nTrials = Int(theta[0].rounded())
                let p = theta[1]
                return RNGSampler.randoms(n: n, dist: .binomial(nTrials: nTrials, p: p), seed: seed)
            }
            
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .binomial,
                    fit: { x in
                        let inferredN = estimateTrials(x)
                        let pHat = MLEFitter<T>.fitBinomial(x, n: inferredN).thetaHat[0]
                        return [T(inferredN), pHat]
                    },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return KSTestResult(
                    distribution: GOFTestTarget.binomial.description,
                    n: 0,
                    d: .nan,
                    dPlus: .nan,
                    dMinus: .nan,
                    pBootstrap: .nan,
                    b: 0,
                    thetaHat: []
                )
            }
        case .cauchy:
            func cdf(from theta: [T]) -> (T) -> T {
                let mu = theta[0]
                let s = max(theta[1], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.Cauchy<T>(location: mu, scale: s)
                return { x in
                    guard let d = dist else { return .nan }
                    do {
                        return try d.cdf(x)
                    } catch {
                        return .nan
                    }
                }
            }
            
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let mu = theta[0], s = max(theta[1], T.leastNonzeroMagnitude)
                return RNGSampler.randoms(n: n, dist: .cauchy(mu: mu, gamma: s), seed: seed)
            }
            
            do {
                let ksCauch = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .cauchy,
                    fit: { x in MLEFitter<T>.fitCauchy(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return ksCauch
            } catch {
                return KSTestResult(
                    distribution: GOFTestTarget.cauchy.description,
                    n: 0,
                    d: .nan,
                    dPlus: .nan,
                    dMinus: .nan,
                    pBootstrap: .nan,
                    b: 0,
                    thetaHat: []
                )
            }
        case .chisquared:
            func cdf(from theta: [T]) -> (T) -> T {
                let dist = try? SwiftyBoost.Distribution.ChiSquared<T>(degreesOfFreedom: theta[0])
                return { x in
                    guard let d = dist else { return .nan }
                    do {
                        return try d.cdf(x)
                    } catch {
                        return .nan
                    }
                }
            }
            
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .chiSquare(df: theta[0]), seed: seed)
            }
            
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .chisquared,
                    fit: { x in MLEFitter<T>.fitChiSquared(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return KSTestResult(
                    distribution: GOFTestTarget.chisquared.description,
                    n: 0,
                    d: .nan,
                    dPlus: .nan,
                    dMinus: .nan,
                    pBootstrap: .nan,
                    b: 0,
                    thetaHat: []
                )
            }
        case .exponential:
            func cdf(from theta: [T]) -> (T) -> T {
                let dist = try? SwiftyBoost.Distribution.Exponential<T>(lambda: theta[0])
                return { x in
                    guard let d = dist else { return .nan }
                    do {
                        return try d.cdf(x)
                    } catch {
                        return .nan
                    }
                }
            }
            
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .exponential(rate: theta[0]), seed: seed)
            }
            
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .exponential,
                    fit: { x in MLEFitter<T>.fitExponential(data: x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return KSTestResult(
                    distribution: GOFTestTarget.exponential.description,
                    n: 0,
                    d: .nan,
                    dPlus: .nan,
                    dMinus: .nan,
                    pBootstrap: .nan,
                    b: 0,
                    thetaHat: []
                )
            }
        case .gumbel:
            func cdf(from theta: [T]) -> (T) -> T {
                let dist = try? SwiftyBoost.Distribution.ExtremeValueGumbel<T>(location: theta[0], scale: theta[1])
                return { x in
                    guard let d = dist else { return .nan }
                    do {
                        return try d.cdf(x)
                    } catch {
                        return .nan
                    }
                }
            }
            
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .gumbel(mu: theta[0], beta: theta[1]), seed: seed)
            }
            
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .gumbel,
                    fit: { x in MLEFitter<T>.fitExtremeValue(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return KSTestResult(
                    distribution: GOFTestTarget.gumbel.description,
                    n: 0,
                    d: .nan,
                    dPlus: .nan,
                    dMinus: .nan,
                    pBootstrap: .nan,
                    b: 0,
                    thetaHat: []
                )
            }
        case .fisherF:
            func cdf(from theta: [T]) -> (T) -> T {
                let d1 = max(theta[0], T.leastNonzeroMagnitude)
                let d2 = max(theta[1], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.FisherF<T>(degreesOfFreedom1: d1, degreesOfFreedom2: d2)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .f(d1: theta[0], d2: theta[1]), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .fisherF,
                    fit: { x in MLEFitter<T>.fitFisherF(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.fisherF.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .gamma:
            func cdf(from theta: [T]) -> (T) -> T {
                let k = max(theta[0], T.leastNonzeroMagnitude)
                let th = max(theta[1], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.Gamma<T>(shape: k, scale: th)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .gamma(shape: theta[0], scale: theta[1]), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .gamma,
                    fit: { x in MLEFitter<T>.fitGamma(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.gamma.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .geometric:
            func cdf(from theta: [T]) -> (T) -> T {
                let p = theta[0]
                let dist = try? SwiftyBoost.Distribution.Geometric<T>(probabibilityOfSuccess: p)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .geometric(p: theta[0]), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .geometric,
                    fit: { x in MLEFitter<T>.fitGeometric(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.geometric.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .holtsmark:
            func cdf(from theta: [T]) -> (T) -> T {
                let mu = theta[0]
                let c = max(theta[1], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.Holtsmark<T>(loc: mu, scale: c)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .holtsmark(mu: theta[0], c: theta[1]), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .holtsmark,
                    fit: { x in MLEFitter<T>.fitHoltsmark(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.holtsmark.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .inverseChiSquared:
            func cdf(from theta: [T]) -> (T) -> T {
                let df = max(theta[0], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.InverseChiSquared<T>(degreesOfFreedom: df, scale: nil)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let df = theta[0]
                let y = RNGSampler.randoms(n: n, dist: .chiSquare(df: df), seed: seed)
                return y.map { v in v > 0 ? 1 / v : .leastNonzeroMagnitude }
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .inverseChiSquared,
                    fit: { x in MLEFitter<T>.fitInverseChiSquared(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.inverseChiSquared.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .inverseGamma:
            func cdf(from theta: [T]) -> (T) -> T {
                let a = max(theta[0], T.leastNonzeroMagnitude)
                let b = max(theta[1], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.InverseGamma<T>(shape: a, scale: b)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let a = theta[0], b = theta[1]
                let y = RNGSampler.randoms(n: n, dist: .gamma(shape: a, scale: 1), seed: seed)
                return y.map { yy in b / max(yy, T.leastNonzeroMagnitude) }
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .inverseGamma,
                    fit: { x in MLEFitter<T>.fitInverseGamma(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.inverseGamma.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .inverseNormal:
            func cdf(from theta: [T]) -> (T) -> T {
                let mu = theta[0]
                let lam = max(theta[1], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.InverseNormal<T>(mean: mu, shape: lam)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let mu = theta[0]
                let lambda = theta[1]
                return RNGSampler.randoms(n: n, dist: .inverseNormal(mu: mu, sigma: lambda), seed: seed)
            }
            
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .inverseNormal,
                    fit: { x in MLEFitter<T>.fitInverseNormal(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.inverseNormal.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .landau:
            func cdf(from theta: [T]) -> (T) -> T {
                let loc = theta[0]
                let sc = max(theta[1], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.Landau<T>(location: loc, scale: sc)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .landau(location: theta[0], scale: theta[1]), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .landau,
                    fit: { x in MLEFitter<T>.fitLandau(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.landau.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .laplace:
            func cdf(from theta: [T]) -> (T) -> T {
                let mu = theta[0]
                let b = max(theta[1], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.Laplace<T>(location: mu, scale: b)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let mu = theta[0], b = max(theta[1], T.leastNonzeroMagnitude)
                let u = RNGSampler<T>.randoms(n: n, dist: .uniform01, seed: seed)
                let centered = u.map { $0 - 0.5 }
                return centered.map { t in mu - b * (t >= 0 ? 1 : -1) * T.log(1 - 2 * t.magnitude) }
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .laplace,
                    fit: { x in MLEFitter<T>.fitLaplace(data: x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.laplace.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .logistic:
            func cdf(from theta: [T]) -> (T) -> T {
                let mu = theta[0], s = max(theta[1], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.Logistic<T>(location: mu, scale: s)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let mu = theta[0], s = max(theta[1], T.leastNonzeroMagnitude)
                return RNGSampler.randoms(n: n, dist: .logistic(mu: mu, s: s), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .logistic,
                    fit: { x in MLEFitter<T>.fitLogistic(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.logistic.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .lognormal:
            func cdf(from theta: [T]) -> (T) -> T {
                let mu = theta[0]
                let b = max(theta[1], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.LogNormal<T>(location: mu, scale: b)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let mu = theta[0], b = max(theta[1], T.leastNonzeroMagnitude)
                return RNGSampler<T>.randoms(n: n, dist: .lognormal(lmean: mu, lsigma: b), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .lognormal,
                    fit: { x in MLEFitter<T>.fitLogNormal(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.lognormal.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .mapAiry:
            func cdf(from theta: [T]) -> (T) -> T {
                let loc = theta[0]
                let scale = theta[1]
                let dist = try? SwiftyBoost.Distribution.MapAiry<T>(location: loc, scale: scale)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .mapAiry(location: theta[0], scale: theta[1]), seed: seed)
            }
            
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .mapAiry,
                    fit: { x in MLEFitter<T>.fitMapAiry(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.mapAiry.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .negativeBinomial:
            func cdf(from theta: [T]) -> (T) -> T {
                let r = theta[0]
                let p = theta[1]
                let dist = try? SwiftyBoost.Distribution.NegativeBinomial<T>(successes: r, probabilityOfSuccess: p)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .negativeBinomial(r: theta[0], p: theta[1]), seed: seed)
            }
            
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .negativeBinomial,
                    fit: { x in MLEFitter<T>.fitNegativeBinomial(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.negativeBinomial.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .noncentralBeta:
            func cdf(from theta: [T]) -> (T) -> T {
                let a = max(theta[0], T.leastNonzeroMagnitude)
                let b = max(theta[1], T.leastNonzeroMagnitude)
                let lam = max(theta[2], T.zero)
                let dist = try? SwiftyBoost.Distribution.NonCentralBeta<T>(alpha: a, beta: b, lambda: lam)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let a = theta[0], b = theta[1], lam = theta[2]
                let X = RNGSampler<T>.randoms(n: n, dist: .noncentralChiSquare(k: 2 * a, lambda: 2 * lam), seed: seed.map { $0 ^ 0x0ACB_001 })
                let Y = RNGSampler<T>.randoms(n: n, dist: .chiSquare(df: 2 * b), seed: seed.map { $0 ^ 0x0ACB_002 })
                return zip(X, Y).map { $0 / ($0 + $1) }
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .noncentralBeta,
                    fit: { x in MLEFitter<T>.fitNCBeta(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.noncentralBeta.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .noncentralChiSquared:
            func cdf(from theta: [T]) -> (T) -> T {
                let k = max(theta[0], T.leastNonzeroMagnitude)
                let lam = max(theta[1], T.zero)
                let dist = try? SwiftyBoost.Distribution.NonCentralChiSquared<T>(degreesOfFreedom: k, nonCentrality: lam)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .noncentralChiSquare(k: theta[0], lambda: theta[1]), seed: seed)
            }
            do {
                func thetaGuess(_ x: [T]) -> (k: T, lambda: T) {
                    let ex = SSExamine<T, T>(usingArray: x, name: nil, characterSet: nil)
                    guard let mean = ex.arithmeticMean,
                          let variance = ex.populationVariance,
                          mean.isFinite, variance.isFinite, mean > .zero else {
                        return (T(2.0), T.zero)
                    }
                    let lambdaGuess = max(T.zero, variance / T.two - mean)
                    let kGuess = max(T(0.5), mean - lambdaGuess)
                    return (kGuess, lambdaGuess)
                }
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .noncentralChiSquared,
                    fit: { x in
                        let guess = thetaGuess(x)
                        let result = MLEFitter<T>.fitNCChiSquared(x, degreesOfFreedom: guess.k)
                        let lambdaHat = result.thetaHat.first ?? guess.lambda
                        return [guess.k, lambdaHat]
                    },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.noncentralChiSquared.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .noncentralStudentT:
            func cdf(from theta: [T]) -> (T) -> T {
                let nu = max(theta[2], T.leastNonzeroMagnitude)
                let delta = theta[3]
                let dist = try? SwiftyBoost.Distribution.NonCentralStudentT<T>(degreesOfFreedom: nu, nonCentrality: delta)
                return { x in
                    guard let d = dist else { return .nan }
                    let mu = theta[0], sigma = max(theta[1], T.leastNonzeroMagnitude)
                    let z = (x - mu) / sigma
                    do { return try d.cdf(z) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let mu = theta[0], sigma = max(theta[1], T.leastNonzeroMagnitude), nu = max(theta[2], T.leastNonzeroMagnitude), delta = theta[3]
                return RNGSampler.randoms(n: n, dist: .noncentralT(mu: mu, sigma: sigma, nu: nu, delta: delta), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .noncentralStudentT,
                    fit: { x in
                        let ex = try? SSExamine<T, T>(using: x.map { T($0) }, levelOfMeasurement: .ratio, name: nil, characterSet: nil)
                        let nuGuess = T(max(5.0, (ex?.sampleStandardDeviation ?? T.one)))
                        let r = MLEFitter<T>.fitNCStudentsT(x, degreesOfFreedom: nuGuess)
                        return [0, 1, nuGuess, r.thetaHat[0]]
                    },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.noncentralStudentT.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .noncentralFisherF:
            func cdf(from theta: [T]) -> (T) -> T {
                let d1 = max(theta[0], T.leastNonzeroMagnitude)
                let d2 = max(theta[1], T.leastNonzeroMagnitude)
                let lam = max(theta[2], T.zero)
                let dist = try? SwiftyBoost.Distribution.NonCentralF<T>(degreesOfFreedom1: d1, degreesOfFreedom2: d2, nonCentrality: lam)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let d1 = theta[0], d2 = theta[1], lam = theta[2]
                let x1 = RNGSampler<T>.randoms(n: n, dist: .noncentralChiSquare(k: d1, lambda: lam), seed: seed.map { $0 ^ 0x0ACF_001 })
                let x2 = RNGSampler<T>.randoms(n: n, dist: .chiSquare(df: d2), seed: seed.map { $0 ^ 0x0ACF_002 })
                return zip(x1, x2).map { (a, b) in (a / d1) / (b / d2) }
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .noncentralFisherF,
                    fit: { x in
                        let d1 = T(8), d2 = T(16)
                        let r = MLEFitter<T>.fitNCFisherF(x, df1: d1, df2: d2)
                        return [d1, d2, r.thetaHat.last ?? .nan]
                    },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.noncentralFisherF.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .normal:
            func cdf(from theta: [T]) -> (T) -> T {
                let mu = theta[0]
                let sigma = max(theta[1], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.Normal<T>(mean: mu, sd: sigma)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .gaussian(mean: theta[0], standardDeviation: theta[1]), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .normal,
                    fit: { x in
                        let exLocal: SSExamine<T, T> = SSExamine<T, T>.init(usingArray: x, name: nil, characterSet: nil)
                        let muLocal = exLocal.arithmeticMean ?? T.zero
                        let sigmaLocal = max(exLocal.sampleStandardDeviation ?? T.leastNonzeroMagnitude, T.leastNonzeroMagnitude)
                        return [muLocal, sigmaLocal]
                    },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.normal.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .normalSTD:
            func cdf(from theta: [T]) -> (T) -> T {
                let dist = try? SwiftyBoost.Distribution.Normal<T>(mean: 0, sd: 1)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .gaussian(mean: 0, standardDeviation: 1), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .normalSTD,
                    fit: { _ in [] },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: false,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.normalSTD.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .pareto:
            func cdf(from theta: [T]) -> (T) -> T {
                let scale = theta[0]
                let shape = theta[1]
                let dist = try? SwiftyBoost.Distribution.Pareto<T>(scale: scale, shape: shape)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .pareto(scale: theta[0], shape: theta[1]), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .pareto,
                    fit: { x in
                        // Pareto-I fitter returns [alpha] with scale fixed to sample minimum.
                        // KS expects [scale, shape] = [xm, alpha]. Inject xm here to avoid mismatch.
                        let xm = x.filter { $0.isFinite && $0 > 0 }.min() ?? (x.min() ?? T.leastNonzeroMagnitude)
                        let alpha = MLEFitter<T>.fitPareto(x).thetaHat[0]
                        return [xm, alpha]
                    },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.pareto.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .poisson:
            func cdf(from theta: [T]) -> (T) -> T {
                let lambda = theta[0]
                let dist = try? SwiftyBoost.Distribution.Poisson<T>(lambda: lambda)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .poisson(lambda: theta[0]), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .poisson,
                    fit: { x in MLEFitter<T>.fitPoisson(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.poisson.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .rayleigh:
            func cdf(from theta: [T]) -> (T) -> T {
                let sigma = theta[0]
                let dist = try? SwiftyBoost.Distribution.Rayleigh<T>(scale: sigma)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .rayleigh(scale: theta[0]), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .rayleigh,
                    fit: { x in MLEFitter<T>.fitRayleigh(data: x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.rayleigh.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .sasPoint5:
            func cdf(from theta: [T]) -> (T) -> T {
                let loc = theta[0]
                let sc = max(theta[1], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.SASPoint5<T>(location: loc, scale: sc)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let mu = theta[0], sc = max(theta[1], T.leastNonzeroMagnitude)
                return RNGSampler<T>.randoms(n: n, dist: .sasPoint5(location: mu, sigma: sc), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .sasPoint5,
                    fit: { x in MLEFitter<T>.fitSASPoint5(data: x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.sasPoint5.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .skewNormal:
            func cdf(from theta: [T]) -> (T) -> T {
                let xi = theta[0]
                let omega = max(theta[1], T.leastNonzeroMagnitude)
                let alpha = theta[2]
                let dist = try? SwiftyBoost.Distribution.SkewNormal<T>(location: xi, scale: omega, shape: alpha)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .skewNormal(xi: theta[0], omega: theta[1], alpha: theta[2]), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .skewNormal,
                    fit: { x in MLEFitter<T>.fitSkewNormal(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.skewNormal.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .StudentT:
            func standardize(_ x: [T]) -> [T] {
                let ex: SSExamine<T, T> = SSExamine<T, T>.init(usingArray: x, name: nil, characterSet: nil)
                let mu = ex.arithmeticMean ?? T.zero
                let sd = max(ex.sampleStandardDeviation ?? T.leastNonzeroMagnitude, T.leastNonzeroMagnitude)
                return x.map { ($0 - mu) / sd }
            }
            func cdf(from theta: [T]) -> (T) -> T {
                let nu = max(theta[0], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.StudentT<T>(degreesOfFreedom: nu)
                return { x in
                    // x is expected to be a z-score already when used inside ksStatistic.
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let nu = max(theta[0], T.leastNonzeroMagnitude)
                // Draw from standard t (mu=0, sigma=1)
                let raw = RNGSampler.randoms(n: n, dist: .studentT(mu: 0, sigma: 1, nu: nu), seed: seed)
                return standardize(raw)
            }
            do {
                let zData = standardize(data)
                let res = try bootstrapOneSampleKS(
                    data: zData,
                    targetDistribution: .StudentT,
                    fit: { x in
                        return MLEFitter<T>.fitStudentsT(x).thetaHat // [nu]
                    },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.StudentT.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .triangular:
            func cdf(from theta: [T]) -> (T) -> T {
                let a = theta[1], b = theta[2], c = theta[0]
                let dist = try? SwiftyBoost.Distribution.Triangular<T>(lower: a, mode: c, upper: b)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let a = theta[1], b = theta[2], c = theta[0]
                return RNGSampler<T>.randoms(n: n, dist: .triangular(lower: a, upper: b, mode: c), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .triangular,
                    fit: { x in MLEFitter<T>.fitTriangular(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.triangular.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .uniform:
            func cdf(from theta: [T]) -> (T) -> T {
                let a = theta[0], b = theta[1]
                let dist = try? SwiftyBoost.Distribution.Uniform<T>(lower: a, upper: b)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                let a = theta[0], b = theta[1]
                let u = RNGSampler<T>.randoms(n: n, dist: .uniform01, seed: seed)
                return u.map { a + (b - a) * $0 }
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .uniform,
                    fit: { x in MLEFitter<T>.fitUniform(data: x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.uniform.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        case .weibull:
            func cdf(from theta: [T]) -> (T) -> T {
                let k = max(theta[0], T.leastNonzeroMagnitude)
                let lam = max(theta[1], T.leastNonzeroMagnitude)
                let dist = try? SwiftyBoost.Distribution.Weibull<T>(shape: k, scale: lam)
                return { x in
                    guard let d = dist else { return .nan }
                    do { return try d.cdf(x) } catch { return .nan }
                }
            }
            func sample(_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T] {
                RNGSampler.randoms(n: n, dist: .weibull(k: theta[0], lambda: theta[1]), seed: seed)
            }
            do {
                let res = try bootstrapOneSampleKS(
                    data: data,
                    targetDistribution: .weibull,
                    fit: { x in MLEFitter<T>.fitWeibull(x).thetaHat },
                    cdfFrom: cdf(from:),
                    sampler: sample,
                    samplerSeed: samplerSeed,
                    reestimateEachReplicate: true,
                    B: bsCount
                )
                return res
            } catch {
                return .init(distribution: GOFTestTarget.weibull.description, n: 0, d: .nan, dPlus: .nan, dMinus: .nan, pBootstrap: .nan, b: 0, thetaHat: [])
            }
        }
    }
    
    /// Computes the empirical KS statistics for a dataset against a supplied CDF.
    ///
    /// - Parameters:
    ///   - data: Observations to evaluate.
    ///   - cdf: The hypothesised CDF evaluated at arbitrary points.
    /// - Returns: The tuple `(D, D⁺, D⁻)`.
    /// - Throws: `KSError.emptySample` when the dataset is empty.
    /// - Throws: `KSError.invalidCDF` when the supplied CDF yields NaN or infinity.
    internal static func ksStatistic(
        data: [T],
        cdf: (T) -> T
    ) throws -> (d: T, dPlus: T, dMinus: T) {
        let n = data.count
        guard n > 0 else { throw KSError.emptySample }
        let x = data.sorted()
        
        var dPlus: T = 0
        var dMinus: T = 0
        
        var index = 0
        while index < n {
            let value = x[index]
            var count = 1
            while index + count < n && x[index + count] == value {
                count += 1
            }
            let cumBefore = index
            let edfLeft = T(cumBefore) / T(n)
            let edfRight = T(cumBefore + count) / T(n)
            let cdfRaw = T(cdf(value))
            guard cdfRaw.isFinite, !cdfRaw.isNaN else { throw KSError.invalidCDF }
            let cdfAtValue = min(T.one, max(T.zero, cdfRaw))
            
            let dPlusCandidate = edfRight - cdfAtValue
            if dPlusCandidate > dPlus { dPlus = dPlusCandidate }
            
            let dMinusCandidate = cdfAtValue - edfLeft
            if dMinusCandidate > dMinus { dMinus = dMinusCandidate }
            
            index += count
        }
        return (max(dPlus, dMinus), dPlus, dMinus)
    }

    private static func makeBootstrapSamplerSeed(seed: UInt64?) -> (() -> UInt64?)? {
        guard let seed else { return nil }
        let normalised = seed == 0 ? 0x9E3779B185EBCA87 : seed
        var state = normalised
        var first = true
        return {
            if first {
                first = false
                return normalised
            }
            state &*= 2862933555777941757
            state &+= 3037000493
            return state
        }
    }
    
    /// Bootstraps the one-sample KS statistic for a fitted distribution.
    ///
    /// - Parameters mirror `bootstrapCore` and allow callers to plug in custom fitters/CDFs.
    /// - Returns: A `KSTestResult` describing the observed statistic and bootstrap p-value.
    public static func bootstrapOneSampleKS(
        data: [T],
        targetDistribution: GOFTestTarget,
        fit: ([T]) -> [T],
        cdfFrom thetaToCDF: @escaping ([T]) -> (T) -> T,
        sampler: @escaping (_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T],
        samplerSeed: (() -> UInt64?)? = nil,
        reestimateEachReplicate: Bool = true,
        B: Int
    ) throws -> KSTestResult {
        try bootstrapCore(
            data: data,
            targetDistribution: targetDistribution,
            fit: fit,
            cdfFrom: thetaToCDF,
            sampler: sampler,
            samplerSeed: samplerSeed,
            reestimateEachReplicate: reestimateEachReplicate,
            B: B
        )
    }
    
    private static func bootstrapCore(
        data: [T],
        targetDistribution: GOFTestTarget,
        fit: ([T]) -> [T],
        cdfFrom thetaToCDF: @escaping ([T]) -> (T) -> T,
        sampler: @escaping (_ n: Int, _ theta: [T], _ seed: UInt64?) -> [T],
        samplerSeed: (() -> UInt64?)?,
        reestimateEachReplicate: Bool,
        B: Int
    ) throws -> KSTestResult {
        let n = data.count
        guard n > 0 else { throw KSError.emptySample }
        guard B > 0 else { throw KSError.invalidReplicates }
        
        let thetaHat = fit(data)
        let cdfHat = thetaToCDF(thetaHat)
        
        let (dObs, dPlusObs, dMinusObs) = try ksStatistic(data: data, cdf: cdfHat)
        
        var ge = 0
        var successes = 0
        for _ in 0..<B {
            let seed = samplerSeed?()
            let y = sampler(n, thetaHat, seed)
            let thetaStar = reestimateEachReplicate ? fit(y) : thetaHat
            let cdfStar = thetaToCDF(thetaStar)
            let stats = try? ksStatistic(data: y, cdf: cdfStar)
            guard let dStar = stats?.d, dStar.isFinite else { continue }
            successes += 1
            if dStar >= dObs { ge &+= 1 }
        }
        
        let p: T
        if dObs.isFinite, successes > 0 {
            p = T(ge + 1) / T(successes + 1)
        } else {
            p = .nan
        }
        
        return .init(
            distribution: targetDistribution.description,
            n: n,
            d: dObs,
            dPlus: dPlusObs,
            dMinus: dMinusObs,
            pBootstrap: p,
            b: successes,
            thetaHat: thetaHat
        )
    }
    /// Summary of a one-sample Kolmogorov–Smirnov test run.
    public struct KSTestResult {
        /// Display name of the fitted distribution.
        public let distribution: String
        /// Sample size of the observed dataset.
        public let n: Int
        /// Two-sided KS statistic `D = max(D⁺, D⁻)`.
        public let d: T
        /// Upper deviation `D⁺ = sup_x (F_n(x) − F(x; θ̂))`.
        public let dPlus: T
        /// Lower deviation `D⁻ = sup_x (F(x; θ̂) − F_n(x⁻))`.
        public let dMinus: T
        /// Bootstrap p-value based on `b` replicates.
        public let pBootstrap: T
        /// Number of successful bootstrap replicates.
        public let b: Int
        /// MLE parameters θ̂ under the null distribution.
        public let thetaHat: [T]
    }
    
    /// Errors thrown by the bootstrap-based KS harness.
    public enum KSError: Error {
        case emptySample
        case invalidReplicates
        case invalidCDF
    }
    
}
