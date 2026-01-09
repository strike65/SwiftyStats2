//
//  Created by VT on 14.12.25.
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

extension Inferential.HypothesisTesting.NonParametricTests {
    //
    //  FriedmanTest.swift
    //
    //  Generic Friedman test (one-way repeated measures ANOVA by ranks).
    //
    //  Requires:
    //  - GroupedData<Value, Group> (stable sort)
    //  - Rank<Value, Group, FP> (average ranks + tie metadata)
    //
    //  Notes:
    //  - Input must be a complete block design: n blocks, k treatments, exactly 1 value per cell.
    //  - Ranks are computed within each block (row).
    //  - Optional tie correction uses C = 1 - sum(t^3 - t) / (n * (k^3 - k)) and Q_adj = Q / C.
    //
    
    // MARK: - Friedman test
    
    /// Performs the Friedman test (one-way repeated-measures ANOVA by ranks).
    ///
    /// The input is treated as a complete block design: each row is a block (e.g. subject),
    /// and each column is a treatment/condition. Values are ranked within each block and
    /// compared via rank sums across treatments.
    ///
    /// The chi-square statistic is
    /// `Q = 12/(n*k*(k+1)) * sum_j R_j^2 - 3*n*(k+1)`,
    /// where `R_j` is the rank sum of treatment `j`, with `df = k - 1`.
    ///
    /// When `tieCorrect` is `true`, the statistic is adjusted by the correction factor
    /// `C = 1 - sum(t^3 - t) / (n * (k^3 - k))` and `Q_adj = Q / C`.
    ///
    /// - Parameters:
    ///   - data: Row-major matrix of observations (blocks x treatments).
    ///   - alpha: Significance level in `(0, 1)` used to compute critical values.
    ///   - tieCorrect: Whether to apply the standard tie correction.
    /// - Returns: A ``FriedmanTestResult`` when the design is valid.
    /// - Throws: ``FriedmanTestError`` if the design is invalid or the tie correction is undefined.
    public static func friedmanTest<Value: Comparable & Hashable & Sendable & Codable>(
        data: [[Value]],
        alpha: T = 0.05,
        tieCorrect: Bool = true
    ) throws -> FriedmanTestResult<T>?
    {
        guard !data.isEmpty else { throw FriedmanTestError.emptyData }
        let n = data.count
        let k = data[0].count
        guard k >= 2 else { throw FriedmanTestError.invalidNumberOfTreatments }
        
        for (i, row) in data.enumerated() {
            if row.count != k {
                throw FriedmanTestError.unequalRowLengths(expected: k, got: row.count, rowIndex: i)
            }
        }
        
        var rankSums = Array(repeating: T.zero, count: k)
        var tieSum = T.zero
        
        // Prebuild treatment indices [0, 1, ..., k-1] once.
        var treatmentIndex: [Int] = []
        treatmentIndex.reserveCapacity(k)
        for j in 0..<k { treatmentIndex.append(j) }
        
        // Rank within each block
        for i in 0..<n {
            let row = data[i]
            
            // Stable sort the row and keep original treatment indices as "groups"
            let gd = GroupedData<Value, Int>(data: row, groups: treatmentIndex)
            let sorted = gd.sorted()
            let sortedValues = sorted.sortedData
            let sortedGroups = sorted.sortedGroups // original treatment index per sorted position
            
            // Compute average ranks on sorted values
            let rinfo = Rank<Value, Int, T>(sortedData: sortedValues)
            let sortedRanks = rinfo.ranks
            
            // Map ranks back to original treatment order, then accumulate
            for p in 0..<k {
                let originalJ = sortedGroups[p]
                rankSums[originalJ] += sortedRanks[p]
            }
            
            // Sum tie blocks across all blocks (each block tie contributes t^3 - t)
            tieSum += rinfo.totalTieCorrection
        }
        
        // Mean ranks
        let nFP = T(n)
        var meanRanks = Array(repeating: T.zero, count: k)
        for j in 0..<k {
            meanRanks[j] = rankSums[j] / nFP
        }
        
        // Friedman chi-square statistic (uncorrected)
        // Q = 12/(n*k*(k+1)) * sum(Rj^2) - 3*n*(k+1)
        let kFP = T(k)
        let kp1 = kFP + T.one
        
        var sumR2 = T.zero
        for j in 0..<k {
            sumR2 += rankSums[j] * rankSums[j]
        }
        
        let q = (T(12) / (nFP * kFP * kp1)) * sumR2 - T(3) * nFP * kp1
        
        // Tie correction factor: C = 1 - sum(t^3 - t) / (n * (k^3 - k))
        // Then Q_adj = Q / C
        let k3mk = (kFP * kFP * kFP) - kFP
        var c = T.one
        var qAdj = q
        
        if tieCorrect {
            if k3mk == T.zero { throw FriedmanTestError.invalidNumberOfTreatments }
            c = T.one - (tieSum / (nFP * k3mk))
            // If all values are tied within every block, C can hit 0 (test not informative).
            if c <= T.zero {
                throw FriedmanTestError.invalidTieCorrectionFactor
            }
            qAdj = q / c
        }
        
        let df = k - 1
        
        // Kendall's W: W = Q_adj / (n * (k - 1))
        let w = qAdj / (nFP * T(df))
        
        // Iman-Davenport approximation:
        // F = ((n - 1) * Q_adj) / (n*(k - 1) - Q_adj)
        // df1 = k - 1, df2 = (k - 1)*(n - 1)
        let denomF = (nFP * T(df)) - qAdj
        var fID: T? = nil
        var df1: Int? = nil
        var df2: Int? = nil
        if denomF > T.zero && n >= 2 {
            fID = ((nFP - T.one) * qAdj) / denomF
            df1 = df
            df2 = df * (n - 1)
        }
        var result: FriedmanTestResult<T> = FriedmanTestResult<T>()
        result.numberOfBlocks = n
        result.numberOfTreatments = k
        result.rankSums = rankSums
        result.meanRanks = meanRanks
        result.chiSquare = q
        result.chiSquareTieCorrected = qAdj
        result.degreesOfFreedom = df
        result.criticalValueChiSquare = try criticalValueChiSquare(df: df, alpha: alpha)
        result.tieCorrectionFactor = c
        result.tieSum = tieSum
        result.kendallsW = w
        result.imanDavenportF = fID
        result.imanDavenportDf1 = df1
        result.imanDavenportDf2 = df2
        result.pValueChiSquare = try pValueChiSquare(df: df, chi: qAdj)
        if let ff = fID, let d1 = df1, let d2 = df2 {
            result.pValueImanDavenport = try pValueImanDavenport(f: ff, df1: d1, df2: d2, alpha: alpha)
            result.criticalValueImanDavenport = try criticalValueImanDavenport(df1: d1, df2: d2, alpha: alpha)
        }
        else {
            result.pValueImanDavenport = nil
            result.criticalValueImanDavenport = nil
        }
        return result
    }
    // Convenience: compute p-values via injected survival functions (sf = 1 - cdf).
    fileprivate static func pValueChiSquare(df: Int, chi x: T) throws -> T {
        let chi: Distribution.ChiSquared<T> = try .init(degreesOfFreedom: T(df))
        let cdf = try chi.cdf(x)
        return T.one - cdf
    }
    
    fileprivate static func criticalValueChiSquare(df: Int, alpha: T) throws -> T {
        // Reject H0 if chiSquareTieCorrected >= critical value
        precondition(alpha > T.zero && alpha < T.one)
        let chi = try Distribution.ChiSquared<T>(degreesOfFreedom: T(df))
        return try chi.quantile(T.one - alpha)
    }
    
    fileprivate static func pValueImanDavenport(f: T, df1: Int, df2: Int, alpha: T) throws -> T? {
        let fd: Distribution.FisherF<T> = try .init(degreesOfFreedom1: T(df1), degreesOfFreedom2: T(df2))
        let cdf = try fd.cdf(f)
        return .one - cdf
    }
    
    fileprivate static func criticalValueImanDavenport(df1: Int, df2: Int, alpha: T) throws -> T? {
        precondition(alpha > T.zero && alpha < T.one)
        let f = try Distribution.FisherF<T>(degreesOfFreedom1: T(df1), degreesOfFreedom2: T(df2))
        return try f.quantile(T.one - alpha)
    }
    
    /// Holds the results of the Friedman test.
    ///
    /// Values are optional so the result can be incrementally populated by the implementation.
    public struct FriedmanTestResult<FP: RealLike>: CustomStringConvertible, Codable {
        
        /// Number of blocks (rows) in the complete block design.
        public var numberOfBlocks: Int?
        /// Number of treatments/conditions (columns) in the design.
        public var numberOfTreatments: Int?
        
        // Rank aggregates
        /// Rank sums `R_j` per treatment.
        public var rankSums: [FP]?
        /// Mean ranks per treatment (`R_j / n`).
        public var meanRanks: [FP]?
        
        // Test statistics
        /// Uncorrected Friedman chi-square statistic `Q`.
        public var chiSquare: FP?
        /// Tie-corrected statistic `Q_adj` (equals ``chiSquare`` when no correction is applied).
        public var chiSquareTieCorrected: FP?
        /// Degrees of freedom (`k - 1`).
        public var degreesOfFreedom: Int?
        /// Chi-square critical value at `1 - alpha` and ``degreesOfFreedom``.
        public var criticalValueChiSquare: FP?
        /// Right-tail p-value for the chi-square approximation, evaluated at ``chiSquareTieCorrected``.
        public var pValueChiSquare: FP?
        
        // Tie correction diagnostics
        /// Tie correction factor `C` used to compute ``chiSquareTieCorrected``.
        public var tieCorrectionFactor: FP?
        /// Sum of tie contributions `sum(t^3 - t)` across all blocks.
        public var tieSum: FP?
        
        // Effect size
        /// Kendall's coefficient of concordance `W = Q_adj / (n * (k - 1))`.
        public var kendallsW: FP?
        
        // Iman-Davenport approximation (optional, may be nil if denominator <= 0)
        /// Iman-Davenport `F` approximation (optional).
        public var imanDavenportF: FP?
        /// Numerator degrees of freedom for ``imanDavenportF`` (equals `k - 1` when available).
        public var imanDavenportDf1: Int?
        /// Denominator degrees of freedom for ``imanDavenportF`` (equals `(k - 1) * (n - 1)` when available).
        public var imanDavenportDf2: Int?
        /// Fisher-F critical value at `1 - alpha` for ``imanDavenportDf1`` and ``imanDavenportDf2``.
        public var criticalValueImanDavenport: FP?
        /// Right-tail p-value for the Iman-Davenport approximation.
        public var pValueImanDavenport: FP?
        
        /// Multi-line, human-readable summary of the result.
        public var description: String {
            var descr = String()
            descr.append("FRIEDMAN TEST FOR RELATED SAMPLES\n")
            descr.append("*********************************\n")
            let blocksStr = numberOfBlocks.map(String.init) ?? "nil"
            let treatmentsStr = numberOfTreatments.map(String.init) ?? "nil"
            let dfStr = degreesOfFreedom.map(String.init) ?? "nil"
            let rankSumsStr = rankSums.map({ $0.map { niceNumber($0) }.joined(separator: ", ") }) ?? "nil"
            let meanRanksStr = meanRanks.map({ $0.map { niceNumber($0) }.joined(separator: ", ") }) ?? "nil"
            let tieSumStr = tieSum.map({ niceNumber($0) }) ?? "nil"
            let tieCorrStr = tieCorrectionFactor.map({ niceNumber($0) }) ?? "nil"
            let kendallStr = kendallsW.map({ niceNumber($0) }) ?? "nil"
            let chiStr = chiSquare.map({ niceNumber($0) }) ?? "nil"
            let chiTieStr = chiSquareTieCorrected.map({ niceNumber($0) }) ?? "nil"
            let cvChiStr = criticalValueChiSquare.map({ niceNumber($0) }) ?? "nil"
            let pChiStr = pValueChiSquare.map({ niceNumber($0) }) ?? "nil"
            
            descr.append("blocks: \(blocksStr)\n")
            descr.append("treatments: \(treatmentsStr)\n")
            descr.append("degrees of freedom: \(dfStr)\n")
            descr.append("sum of ranks: \(rankSumsStr)\n")
            descr.append("mean ranks: \(meanRanksStr)\n")
            descr.append("tie sum: \(tieSumStr)\n")
            descr.append("tie correction factor: \(tieCorrStr)\n")
            descr.append("Kendall's W: \(kendallStr)\n")
            descr.append("chi-square (uncorrected): \(chiStr)\n")
            descr.append("chi-square (tie-corrected): \(chiTieStr)\n")
            descr.append("critical value (chi-square): \(cvChiStr)\n")
            descr.append("p-value (chi-square): \(pChiStr)\n")
            
            if let f = imanDavenportF,
               let df1 = imanDavenportDf1,
               let df2 = imanDavenportDf2 {
                descr.append("Iman-Davenport F: \(niceNumber(f))\n")
                descr.append("Iman-Davenport df1: \(df1)\n")
                descr.append("Iman-Davenport df2: \(df2)\n")
                if let cvF = criticalValueImanDavenport {
                    descr.append("critical value (Iman-Davenport F): \(niceNumber(cvF))\n")
                }
                if let pF = pValueImanDavenport {
                    descr.append("p-value (Iman-Davenport): \(niceNumber(pF))\n")
                }
            }
            return descr
        }
    }
    
    /// Errors that may be thrown when evaluating the Friedman test.
    public enum FriedmanTestError: Error, Sendable {
        /// The input matrix is empty.
        case emptyData
        /// The design has fewer than two treatments/conditions.
        case invalidNumberOfTreatments
        /// At least one row has a different length than the first row.
        case unequalRowLengths(expected: Int, got: Int, rowIndex: Int)
        /// The tie correction factor is non-positive (e.g. all values tied within each block).
        case invalidTieCorrectionFactor
    }
}
