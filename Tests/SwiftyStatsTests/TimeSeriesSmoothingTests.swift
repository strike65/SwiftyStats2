/// Tests for time series smoothing utilities and forecast helpers.

import Testing
@testable import SwiftyStats

struct TimeSeriesSmoothingTests {
    private let values = (1...25).map { Double($0) }
    private let times = (0..<25).map { Double($0) }

    @Test
    func movingAverageWindowThreeMatchesExpected() throws {
        let ts = TimeSeries(times: times, values: values)
        let expected = (2...24).map { Double($0) }

        let result = try #require(ts.movingAverage(window: 3))
        #expect(result.count == expected.count)
        zip(result, expected).forEach { actual, expected in
            #expect(abs(actual - expected) < 1e-12)
        }
    }

    @Test
    func singleExponentialSmoothingAlphaHalf() throws {
        let ts = TimeSeries(times: times, values: values)
        var expected = Array(repeating: 0.0, count: values.count)
        expected[0] = values[0]
        for idx in 1..<values.count {
            expected[idx] = 0.5 * values[idx] + 0.5 * expected[idx - 1]
        }

        let result = try #require(ts.singleExponentialSmoothing(alpha: 0.5))
        #expect(result.count == expected.count)
        zip(result, expected).forEach { actual, expected in
            #expect(abs(actual - expected) < 1e-12)
        }
    }

    @Test
    func autocorrelationMatchesReferenceComputation() throws {
        let ts = TimeSeries(times: times, values: values)
        let expected = referenceACF(values: values, maxLag: 4)

        let result = try #require(ts.autocorrelationFunction(maxLag: 4))
        #expect(result.count == expected.count)
        zip(result, expected).forEach { actual, expected in
            #expect(abs(actual - expected) < 1e-12)
        }
    }

    // MARK: - Helpers

    private func referenceACF(values: [Double], maxLag: Int) -> [Double] {
        let n = values.count
        let mean = values.reduce(0, +) / Double(n)

        func autocovariance(lag: Int) -> Double {
            let limit = n - lag
            let sum = (0..<limit).reduce(0.0) { partial, idx in
                let x0 = values[idx] - mean
                let xk = values[idx + lag] - mean
                return partial + x0 * xk
            }
            return sum / Double(n)
        }

        let gamma0 = autocovariance(lag: 0)
        var result: [Double] = [1.0]
        for k in 1...maxLag {
            result.append(autocovariance(lag: k) / gamma0)
        }
        return result
    }
}
