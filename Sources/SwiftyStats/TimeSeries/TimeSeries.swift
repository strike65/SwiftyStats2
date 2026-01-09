//
//  Created by VT on 23.11.25.
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

/// Generic time series container with smoothing and autocorrelation helpers.

import SwiftyStatsPrelude

/// A single observation in a time series.
///
/// - Parameters:
///   - Time: Comparable, codable timestamp type.
///   - Value: Codable payload type.
public struct TimeSeriesPoint<Time: Comparable & Codable, Value: Codable>: Codable {
    /// The observation time.
    public let time: Time
    /// The observed value.
    public let value: Value
    
    /// Creates a point with a time and value.
    ///
    /// - Parameters:
    ///   - time: The timestamp.
    ///   - value: The observed value.
    @inlinable
    public init(time: Time, value: Value) {
        self.time = time
        self.value = value
    }
}

/// An ordered collection of time-stamped values.
///
/// Points are maintained in non-decreasing order by `time` when constructed
/// from unsorted inputs or via `insert(_:)`. Use `appendSorted(_:)` to
/// append in already-sorted scenarios with a precondition check.
///
/// - Parameters:
///   - Time: Comparable, codable timestamp type.
///   - Value: Codable payload type.
public struct TimeSeries<Time: Comparable & Codable, Value: Codable>: Codable {
    /// The underlying storage of points, expected to be sorted by time.
    public var points: [TimeSeriesPoint<Time, Value>]
    
    /// Creates an empty time series.
    @inlinable
    public init() {
        self.points = []
    }
    
    /// Creates a time series from points assumed to be already sorted by time.
    ///
    /// - Important: No sorting is performed; caller guarantees order.
    /// - Parameter sortedPoints: Points in non-decreasing time order.
    @inlinable
    public init(sortedPoints: [TimeSeriesPoint<Time, Value>]) {
        self.points = sortedPoints
    }
    
    /// Creates a time series from unsorted points.
    ///
    /// The points are sorted by ascending `time`.
    /// - Parameter unsortedPoints: Unordered points.
    @inlinable
    public init(unsortedPoints: [TimeSeriesPoint<Time, Value>]) {
        self.points = unsortedPoints.sorted(by: { $0.time < $1.time })
    }
    
    /// Creates a time series from parallel arrays of times and values.
    ///
    /// - Precondition: `times.count == values.count`.
    /// - Parameters:
    ///   - times: Timestamps.
    ///   - values: Values corresponding to `times`.
    @inlinable
    public init(times: [Time], values:[Value]) {
        precondition(times.count == values.count, "times and values must have the same length")
        
        var temp: [TimeSeriesPoint<Time, Value>] = []
        for i in 0..<times.count {
            temp.append(TimeSeriesPoint(time: times[i], value: values[i]))
        }
        self.points = temp.sorted(by: { $0.time < $1.time })
    }
    
    /// The number of points.
    @inlinable
    public var count: Int {
        return self.points.count
    }
    
    /// The earliest time, or `nil` if empty.
    @inlinable
    public var startTime: Time? {
        return self.points.first?.time
    }
    
    /// The latest time, or `nil` if empty.
    @inlinable
    public var endTime: Time? {
        return self.points.last?.time
    }
    
    /// All times in order.
    @inlinable
    public var allTimes: [Time] {
        return self.points.map{ $0.time }
    }
    /// All values in time order.
    @inlinable
    public var allValues: [Value] {
        return self.points.map{ $0.value }
    }
    
    /// Accesses the point at `index`.
    @inlinable
    public subscript(index: Int) -> TimeSeriesPoint<Time, Value> {
        get {
            return self.points[index]
        }
        set {
            self.points[index] = newValue
        }
    }
    
    /// Returns the index of `queryTime` or the next greater time.
    ///
    /// Uses binary search on the sorted points.
    ///
    /// - Parameter queryTime: The time to locate.
    /// - Returns: The index of the first point with `time >= queryTime`, or `nil` if none.
    @inlinable
    public func indexOfTimeOrNext(_ queryTime: Time) -> Int? {
        guard !self.points.isEmpty else {
            return nil
        }
        var low = 0
        var high = self.points.count
        while low < high {
            let mid = (low + high) / 2
            if self.points[mid].time < queryTime {
                low = mid + 1
            }
            else {
                high = mid
            }
        }
        return low < self.points.count ? low : nil
    }
    
    /// Appends a point that is not earlier than the current last point.
    ///
    /// - Precondition: `point.time >= last.time` if the series is non-empty.
    /// - Parameter point: The point to append.
    @inlinable
    public mutating func appendSorted(_ point: TimeSeriesPoint<Time, Value>) {
        if let lastTime = self.points.last?.time {
            precondition(point.time >= lastTime, "appendSorted: points must be in non-decreasing order by time")
        }
        self.points.append(point)
    }
    
    /// Inserts a point while maintaining sort order by time.
    ///
    /// Uses binary search to find the insertion index.
    /// - Parameter point: The point to insert.
    @inlinable
    public mutating func insert(_ point: TimeSeriesPoint<Time, Value>) {
        if self.points.isEmpty {
            self.points.append(point)
            return
        }
        
        var low = 0
        var high = self.points.count
        while low < high {
            let mid = (low + high) / 2
            if self.points[mid].time < point.time {
                low = mid + 1
            }
            else {
                high = mid
            }
        }
        self.points.insert(point, at: low)
    }
}

extension TimeSeries where Value: RealLike {
    // MARK: - Moving averages

    /// Simple moving average (trailing window) of size `window`.
    /// Returns an array of length `count - window + 1`, or `nil` if window is invalid.
    @inlinable
    public func movingAverage(window: Int) -> [Value]? {
        let n = points.count
        guard window > 0 && window <= n else { return nil }
        if window == 1 { return valuesAsArray() }

        var result: [Value] = []
        result.reserveCapacity(n - window + 1)

        var acc: Value = 0
        for i in 0..<window {
            acc += points[i].value
        }
        result.append(acc / Value(window))

        if n == window { return result }

        for i in window..<n {
            acc += points[i].value
            acc -= points[i - window].value
            result.append(acc / Value(window))
        }
        return result
    }

    /// Centered moving average of window size `window` (must be odd).
    /// Returns an array aligned to the center; leading/trailing points without a full window are omitted.
    @inlinable
    public func centeredMovingAverage(window: Int) -> [Value]? {
        let n = points.count
        guard window > 0 && window % 2 == 1 && window <= n else { return nil }

        var result: [Value] = []
        result.reserveCapacity(n - window + 1)

        var acc: Value = 0
        for i in 0..<window {
            acc += points[i].value
        }
        result.append(acc / Value(window))

        if n == window { return result }

        for i in window..<n {
            acc += points[i].value
            acc -= points[i - window].value
            result.append(acc / Value(window))
        }
        return result
    }

    // MARK: - Exponential smoothing

    /// Single exponential smoothing (SES). Returns smoothed series of length `count`,
    /// initialized at the first observation.
    @inlinable
    public func singleExponentialSmoothing(alpha: Value) -> [Value]? {
        let n = points.count
        guard n > 0 else { return nil }
        guard alpha >= 0 && alpha <= 1 else { return nil }

        var result = [Value](repeating: 0, count: n)
        result[0] = points[0].value
        if n == 1 { return result }

        for i in 1..<n {
            let x = points[i].value
            result[i] = alpha * x + (1 - alpha) * result[i - 1]
        }
        return result
    }

    /// Double exponential smoothing (Holt's linear trend).
    /// Returns smoothed level/trend forecasts of length `count`.
    @inlinable
    public func doubleExponentialSmoothing(alpha: Value, beta: Value) -> [Value]? {
        let n = points.count
        guard n >= 2 else { return nil }
        guard alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1 else { return nil }

        var level = points[0].value
        var trend = points[1].value - points[0].value
        var result = [Value](repeating: 0, count: n)
        result[0] = level

        for i in 1..<n {
            let x = points[i].value
            let prevLevel = level
            level = alpha * x + (1 - alpha) * (level + trend)
            trend = beta * (level - prevLevel) + (1 - beta) * trend
            result[i] = level + trend
        }
        return result
    }

    /// Triple exponential smoothing (Holt-Winters additive seasonality).
    /// - Parameters:
    ///   - alpha: level smoothing (0...1)
    ///   - beta: trend smoothing (0...1)
    ///   - gamma: seasonal smoothing (0...1)
    ///   - seasonLength: length of season (>1)
    /// - Returns: Smoothed forecasts of length `count`.
    @inlinable
    public func tripleExponentialSmoothing(
        alpha: Value,
        beta: Value,
        gamma: Value,
        seasonLength: Int
    ) -> [Value]? {
        let n = points.count
        guard n >= seasonLength * 2 else { return nil }
        guard alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 else { return nil }
        guard seasonLength > 1 else { return nil }

        var level: Value = 0
        var trend: Value = 0
        var season = [Value](repeating: 0, count: seasonLength)

        // Initialize level and trend from first season
        for i in 0..<seasonLength {
            level += points[i].value
        }
        level /= Value(seasonLength)
        trend = (points[seasonLength].value - points[0].value) / Value(seasonLength)

        // Initialize seasonals
        for i in 0..<seasonLength {
            season[i] = points[i].value - level
        }

        var result = [Value](repeating: 0, count: n)
        for t in 0..<n {
            let s = season[t % seasonLength]
            result[t] = level + trend + s

            // Update components only if we have a full season ahead
            if t >= seasonLength {
                let x = points[t].value
                let lastLevel = level
                level = alpha * (x - season[t % seasonLength]) + (1 - alpha) * (level + trend)
                trend = beta * (level - lastLevel) + (1 - beta) * trend
                season[t % seasonLength] = gamma * (x - level) + (1 - gamma) * season[t % seasonLength]
            }
        }
        return result
    }

    /// The arithmetic mean of all values, or `nil` if empty.
    ///
    /// Uses a numerically stable summation.
    public func meanValue() -> Value? {
        guard !self.points.isEmpty else {
            return nil
        }
        var sum: Value = Value.zero
        var values: [Value] = []
        for p in points {
            values.append(p.value)
        }
        sum = Helpers.sum(&values)
        return sum / Value(self.points.count)
    }
    
    @inlinable
    func variance(unbiased: Bool) -> Value? {
        let examine: SSExamine<Value, Value> = SSExamine<Value, Value>.init(usingArray: self.valuesAsArray(), name: nil, characterSet: nil)
        if unbiased {
            return examine.unbiasedVariance ?? nil
        }
        else {
            return examine.biasedVariance ?? nil
        }
    }
    
    func covariance(with other: TimeSeries<Time, Value>, unbiased: Bool) -> Value? {
        let n = other.points.count
        guard n >= 2 else { return nil }
        guard n == other.points.count else { return nil }
        if let meanX = self.meanValue(), let meanY = other.meanValue() {
            var sum: Value = 0
            var terms : [Value] = []
            for i in 0..<n {
                let dx = self.points[i].value - meanX
                let dy = other.points[i].value - meanY
                terms.append(dx * dy)
            }
            sum = Helpers.sum(&terms)
            let denom = unbiased ? Value(n - 1) : Value(n)
            return sum / denom
        }
        else {
            return nil
        }
    }
    
    /// Returns all values as an array in time order.
    @inlinable
    public func valuesAsArray() -> [Value] {
        return allValues
    }
    
    // MARK: - Autocovariance and autocorrelation

    /// Biased sample autocovariance at integer lag k (k >= 0).
    /// Returns nil if series too short or lag invalid.
    @inlinable
    public func autocovariance(lag k: Int) -> Value? {
        let n = points.count
        guard n >= 2 else { return nil }
        guard k >= 0 && k < n else { return nil }
        guard let m = meanValue() else { return nil }

        if k == 0 {
            // gamma(0) = variance * n / (n-1) if you use unbiased variance,
            // but we compute directly with division by n:
            var sum: Value = 0
            for p in points {
                let d = p.value - m
                sum += d * d
            }
            return sum / Value(n)
        }

        let limit = n - k
        var sum: Value = 0
        for i in 0..<limit {
            let x0 = points[i].value - m
            let xk = points[i + k].value - m
            sum += x0 * xk
        }
        return sum / Value(n)
    }

    /// Sample autocorrelation at integer lag k (k >= 0).
    /// Uses biased autocovariance gamma_hat(k) / gamma_hat(0).
    @inlinable
    public func autocorrelation(lag k: Int) -> Value? {
        guard let gamma0 = autocovariance(lag: 0) else { return nil }
        if gamma0 == 0 {
            return nil // constant series -> undefined correlation
        }
        guard let gammak = autocovariance(lag: k) else { return nil }
        return gammak / gamma0
    }

    /// Autocorrelation function from lag 0 to maxLag (inclusive).
    /// Returns array rho[0...maxLag].
    @inlinable
    public func autocorrelationFunction(maxLag: Int) -> [Value]? {
        let n = points.count
        guard n >= 2 else { return nil }
        guard maxLag >= 0 && maxLag < n else { return nil }

        var result = [Value]()
        result.reserveCapacity(maxLag + 1)

        guard let gamma0 = autocovariance(lag: 0), gamma0 != 0 else {
            return nil
        }

        result.append(1) // lag 0

        for k in 1...maxLag {
            guard let gammak = autocovariance(lag: k) else {
                return nil
            }
            result.append(gammak / gamma0)
        }

        return result
    }
}
