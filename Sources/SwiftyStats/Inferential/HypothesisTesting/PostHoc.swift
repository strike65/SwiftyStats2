//
//  Created by VT on 21.12.25.
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

import SwiftyStatsPrelude

extension Inferential.HypothesisTesting.PostHocTests {
    
    public static func scheffeTest(data: DataFrame<T,T>, alpha: T) throws -> Array<PostHocTestSummary>? {
        return try scheffeTest(data: data.data, alpha: alpha)
    }
    
    public static func scheffeTest(data: [SSExamine<T,T>], alpha: T) throws -> Array<PostHocTestSummary>? {
        guard let tukeyResult = try tukeyKramerTest(data: data, alpha: alpha) else { return nil }
        var scheffeResults: [PostHocTestSummary] = []
        var n_total: T = .zero
        let k: T = T(data.count)
        var df_error: T
        var tempRes: PostHocTestSummary
        for ex in data {
            n_total += T(ex.sampleSize)
        }
        df_error = n_total - k
        let dist = try Distribution.FisherF(degreesOfFreedom1: k - .one, degreesOfFreedom2: df_error)
        for tk in tukeyResult {
            tempRes.meanDiff = tk.meanDiff
            tempRes.row = tk.row
            tempRes.testType = .scheffe
            tempRes.testStat = tk.testStat / T.sqrt(.two)
            tempRes.pValue = T.one - (try dist.cdf(T.pow(tempRes.testStat, 2)))
            scheffeResults.append(tempRes)
        }
        return scheffeResults
    }
    
    public static func tukeyKramerTest(data: DataFrame<T, T>, alpha: T) throws -> Array<PostHocTestSummary>? {
        return try tukeyKramerTest(data: data.data, alpha: alpha)
    }
    
    public static func tukeyKramerTest(data: Array<SSExamine<T, T>>, alpha: T) throws -> Array<PostHocTestSummary>? {
        guard data.count > 2 else { return nil }
        var means: [T] = []
        var i: Int = 1
        var names: Set<String> = Set<String>()
        for e in data {
            if let m = e.arithmeticMean {
                means.append(m)
            }
            else {
                return nil
            }
            if let n = e.name {
                if !names.insert(n).inserted {
                    throw SSError.init(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: "Duplicate group name: \(n)")
                }
            }
            else {
                e.name = "Group \(i)"
            }
            i += 1
        }
        var Q: [[T]] = Array<Array<T>>()
        var differences: [[T]] = Array<Array<T>>()
        var N_total: T = .zero
        for i in stride(from: 0, through: data.count - 1, by: 1) {
            differences.append(Array<T>())
            for j in stride(from: 0, through: data.count - 1, by: 1) {
                differences[i].append(.nan)
                if j >= i + 1 {
                    let diff: T = means[i] - means[j]
                    differences[i][j] = diff
                }
            }
        }
        var sse: T = .zero
        for i in stride(from: 0, through: data.count - 1, by: 1) {
            N_total += T(data[i].sampleSize)
            if let variance = data[i].sampleVariance {
                sse += T(data[i].sampleSize - 1) * variance
            }
            else {
                return nil
            }
        }
        let dfError: T = N_total - T(data.count)
        guard dfError > .zero else { return nil }
        let s: T = T.sqrt(sse / dfError)
        var n_i: T
        var n_j: T
        for i in stride(from: 0, through: data.count - 1, by: 1) {
            Q.append(Array<T>())
            n_i = T(data[i].sampleSize)
            for j in stride(from: 0, through: data.count - 1, by: 1) {
                Q[i].append(T.nan)
                if j >= i + 1 {
                    n_j = T(data[j].sampleSize)
                    let se: T = s * T.sqrt((.one / n_i + .one / n_j) / .two)
                    Q[i][j] = differences[i][j].magnitude / se
                }
            }
        }
        var pValues: [[T]] = Array<Array<T>>()
        var temp: T = .nan
        for i in stride(from: 0, through: data.count - 1, by: 1) {
            pValues.append(Array<T>())
            for j in stride(from: 0, through: data.count - 1, by: 1) {
                pValues[i].append(.nan)
                if j >= i + 1 {
                    temp = .nan
                    temp = try Inferential.HypothesisTesting.ParametricTests.ptukey(q: Q[i][j], nranges: 1, numberOfMeans: T(data.count), df: N_total - T(data.count), tail: .upper, returnLogP: false)
                }
                pValues[i][j] = temp
            }
        }
        var confIntv: Array<Array<ConfidenceInterval<T>>> = Array<Array<ConfidenceInterval<T>>>()
        var lb: T, ub: T, hwidth: T
        let criticalQ: T = try Inferential.HypothesisTesting.ParametricTests.qtukey(p: .one - alpha, nranges: .one, numberOfMeans: T(data.count), df: N_total - T(data.count), tail: .lower, log_p: false)
        for i in stride(from: 0, through: data.count - 1, by: 1) {
            confIntv.append(Array<ConfidenceInterval<T>>())
            n_i = T(data[i].sampleSize)
            for j in stride(from: 0, through: data.count - 1, by: 1) {
                confIntv[i].append(.init(lowerBound: .nan, upperBound: .nan))
                n_j = T(data[j].sampleSize)
                if j >= i + 1 {
                    hwidth = criticalQ * s * T.sqrt((.one / n_i + .one / n_j) / .two)
                    lb = differences[i][j] - hwidth
                    ub = differences[i][j] + hwidth
                    confIntv[i][j].lowerBound = lb
                    confIntv[i][j].upperBound = ub
                    confIntv[i][j].width = .two * hwidth
                }
            }
        }
        var summary: Array<PostHocTestSummary> = Array<PostHocTestSummary>()
        for i in stride(from: 0, through: data.count - 2, by: 1) {
            for j in stride(from: i + 1, through: data.count - 1, by: 1) {
                if let ni = data[i].name, let nj = data[j].name {
                    summary.append((row: ni + " - " + nj + ":", meanDiff: differences[i][j], testStat: Q[i][j], pValue: pValues[i][j], testType: .tukeyKramer))
                }
                else {
                    return nil
                }
            }
        }
        return summary
    }
    
    public static func nicePostHocTestResults(results: [PostHocTestSummary]) -> String {
        var descr: String = ""
        descr.append("****************************\n")
        descr.append("POST HOC TEST RESULT SUMMARY\n")
        descr.append("****************************\n")
        for phr in results {
            descr.append("\(phr.row)\tdiff: \(phr.meanDiff)\tQ: \(phr.testStat)\tp: \(phr.pValue)\n")
        }
        descr.append("****************************\n")
        return descr
    }
    
    /// A tuple containing the results of one out of multiple comparisons.
    public typealias PostHocTestSummary = (row: String, meanDiff: T, testStat: T, pValue: T, testType: PostHocTestType)
    
    public enum PostHocTestType: Int, Codable, CustomStringConvertible {
        case tukeyKramer, scheffe
        public var description: String {
            switch self {
            case .tukeyKramer:
                return "Tukey/Kramer"
            case .scheffe:
                return "Scheffé"
            }
        }

    }


}
