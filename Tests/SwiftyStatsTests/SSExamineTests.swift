/// Core SSExamine behavior and thread-safety test coverage.

import Foundation
import Testing
@testable import SwiftyStats

@Test("SSExamine maintains frequencies and ordering")
func ssExamineFrequenciesAndOrdering() {
    let examine = SSExamine<Int, Double>()
    #expect(examine.isEmpty)

    examine.append(1)
    let appended = examine.append(repeating: 2, item: 2)
    #expect(appended)

    #expect(examine.sampleSize == 3)
    #expect(examine.length == 2)
    #expect(examine.frequency(1) == 1)
    #expect(examine.frequency(2) == 2)
    #expect(examine.frequency(3) == 0)

    let relFrequency = examine.relativeFrequency(2)
    #expect(abs(relFrequency - (2.0 / 3.0)) < 1e-9)

    #expect(examine.itemsAsArray(sorted: .original) == [1, 2, 2])
    #expect(examine.itemsAsArray(sorted: .ascending) == [1, 2, 2])
    #expect(examine.itemsAsArray(sorted: .descending) == [2, 2, 1])
    #expect(examine.uniqueItems(sorted: .ascending) == [1, 2])

    examine.name = "Sample"
    let row = examine.itemsAsString(delimiter: ",", asRow: true)
    #expect(row == "1,2,2")

    let column = examine.itemsAsString(delimiter: nil, asRow: false)
    #expect(column == "Sample\n1\n2\n2")
}

@Test("SSExamine removal semantics respect raw order")
func ssExamineRemovalSemantics() {
    let examine = SSExamine<Int, Double>()
    let appended = examine.append(contentOf: [1, 2, 2, 3, 2])
    #expect(appended)
    #expect(examine.itemsAsArray(sorted: .original) == [1, 2, 2, 3, 2])

    examine.remove(2, allOccurences: false)
    #expect(examine.itemsAsArray(sorted: .original) == [1, 2, 3, 2])
    #expect(examine.sampleSize == 4)
    #expect(examine.frequency(2) == 2)

    examine.remove(2, allOccurences: true)
    #expect(examine.itemsAsArray(sorted: .original) == [1, 3])
    #expect(examine.sampleSize == 2)
    #expect(examine.frequency(2) == 0)

    examine.removeAll()
    #expect(examine.isEmpty)
    #expect(examine.sampleSize == 0)
}

@Test("SSExamine supports string ingestion and codable round trip")
func ssExamineStringCodableRoundTrip() throws {
    let examine = SSExamine<String, Double>()
    examine.name = "letters"
    examine.tag = "id-42"
    examine.alpha = 0.1
    examine.levelOfMeasurement = .ordinal

    let appended = try examine.append(text: "Hello, 42!", characterSet: .letters)
    #expect(appended)
    #expect(examine.itemsAsArray(sorted: .original) == ["H", "e", "l", "l", "o"])
    #expect(examine.sampleSize == 5)
    #expect(examine.length == 4)

    let csv = examine.itemsAsString(delimiter: ",", asRow: true, encloseElementsBy: "\"")
    #expect(csv == "\"H\",\"e\",\"l\",\"l\",\"o\"")

    let copy = examine.copy() as? SSExamine<String, Double>
    #expect(copy != nil)
    #expect(copy?.isEqual(examine) == true)

    let encoder = JSONEncoder()
    let data = try encoder.encode(examine)
    let decoder = JSONDecoder()
    let decoded = try decoder.decode(SSExamine<String, Double>.self, from: data)

    #expect(decoded.isEqual(examine))
    #expect(decoded.name == examine.name)
    #expect(decoded.tag == examine.tag)
    #expect(examine.alpha == 0.1)
    #expect(decoded.alpha == 0.05)
    #expect(decoded.levelOfMeasurement == .nominal)
    #expect(decoded.itemsAsArray(sorted: .original) == examine.itemsAsArray(sorted: .original))
}

@Test("SSExamine exposes mutable metadata properties")
func ssExamineMetadataProperties() {
    let examine = SSExamine<Int, Double>()
    #expect(examine.rootObject == nil)
    #expect(examine.tag == nil)
    #expect(examine.name == nil)
    #expect(examine.levelOfMeasurement == .nominal)
    #expect(examine.descriptionString == "SSExamine instance - standard")
    #expect(abs(Double(examine.alpha) - 0.05) < 1e-12)
    #expect(examine.hasChanges == false)
    #expect(examine.isNumeric)
    #expect(examine.isEmpty)
    #expect(examine.sampleSize == 0)
    #expect(examine.length == 0)

    examine.rootObject = ["meta": 7]
    examine.tag = "tag-1"
    examine.name = "dataset"
    examine.levelOfMeasurement = .ordinal
    examine.descriptionString = "Custom description"
    examine.alpha = 0.01
    examine.hasChanges = true
    examine.isNumeric = false

    let meta = examine.rootObject as? [String: Int]
    #expect(meta?["meta"] == 7)
    #expect(examine.tag == "tag-1")
    #expect(examine.name == "dataset")
    #expect(examine.levelOfMeasurement == .ordinal)
    #expect(examine.descriptionString == "Custom description")
    #expect(abs(Double(examine.alpha) - 0.01) < 1e-12)
    #expect(examine.hasChanges)
    #expect(examine.isNumeric == false)
}

@Test("SSExamine throwing initializer variants populate data")
func ssExamineInitializerVariants() throws {
    let numbers = try SSExamine<Int, Double>(using: [3, 1, 2, 1], levelOfMeasurement: .ratio, name: "numbers", characterSet: nil)
    #expect(numbers.name == "numbers")
    #expect(numbers.tag == nil)
    // Initialiser resets level to nominal through initializeLocked.
    #expect(numbers.levelOfMeasurement == .nominal)
    #expect(numbers.isNumeric)
    #expect(numbers.sampleSize == 4)
    #expect(numbers.length == 3)
    #expect(numbers.itemsAsArray(sorted: .original) == [3, 1, 2, 1])
    #expect(numbers.itemsAsArray(sorted: .none).sorted() == numbers.itemsAsArray(sorted: .ascending))
    #expect(numbers.contains(2))
    let relativeThree = numbers.relativeFrequency(3)
    #expect(abs(Double(relativeThree) - 0.25) < 1e-12)
    #expect(numbers[0] == 3)
    #expect(numbers[4] == nil)
    #expect(numbers[-1] == nil)
    #expect(numbers.uniqueItems(sorted: .descending) == [3, 2, 1])
    let uniqueNone = numbers.uniqueItems(sorted: .none)
    #expect(Set(uniqueNone) == Set([1, 2, 3]))

    let letters = try SSExamine<String, Double>(using: "a1BcA", levelOfMeasurement: .interval, name: nil, characterSet: .letters)
    #expect(letters.name != nil)
    #expect(letters.tag == letters.name)
    #expect(letters.levelOfMeasurement == .nominal)
    #expect(letters.isNumeric == false)
    #expect(letters.itemsAsArray(sorted: .original) == ["a", "B", "c", "A"])
    #expect(letters.sampleSize == 4)
    #expect(letters.length == 4)
    let relB = letters.relativeFrequency("B")
    #expect(abs(Double(relB) - 0.25) < 1e-12)

    #expect(throws: SSError.self) {
        _ = try SSExamine<Int, Double>(using: "abc", levelOfMeasurement: .interval, name: nil, characterSet: nil)
    }
}

@Test("SSExamine append variants handle edge cases and errors")
func ssExamineAppendVariants() throws {
    let numeric = SSExamine<Int, Double>()
    #expect(numeric.append(repeating: 0, item: 9) == false)
    #expect(numeric.append(contentOf: []) == false)
    #expect(numeric.isEmpty)
    #expect(numeric.hasChanges == false)

    #expect(numeric.append(contentOf: [1, 2]))
    #expect(numeric.sampleSize == 2)
    #expect(numeric.hasChanges)
    #expect(numeric.append(repeating: 3, item: 2))
    #expect(numeric.sampleSize == 5)
    #expect(numeric.frequency(2) == 4)
    let relTwo = numeric.relativeFrequency(2)
    #expect(abs(Double(relTwo) - (4.0 / 5.0)) < 1e-12)

    #expect(throws: SSError.self) {
        _ = try numeric.append(text: "abc", characterSet: nil)
    }

    let strings = SSExamine<String, Double>()
    #expect(try strings.append(text: "", characterSet: nil) == false)
    #expect(try strings.append(text: "123", characterSet: .letters) == false)
    let filtered = try strings.append(text: "ab1", characterSet: .letters)
    #expect(filtered)
    #expect(strings.itemsAsArray(sorted: .original) == ["a", "b"])
    #expect(strings.isNumeric)
}

@Test("SSExamine renders arrays and strings across orders")
func ssExamineArrayAndStringRendering() {
    let empty = SSExamine<Int, Double>()
    #expect(empty.itemsAsString(delimiter: ",") == nil)
    #expect(empty.itemsAsArray(sorted: .ascending).isEmpty)
    #expect(empty.uniqueItems(sorted: .ascending).isEmpty)

    let examine = SSExamine<Int, Double>()
    #expect(examine.append(contentOf: [3, 1, 2, 1]))
    examine.name = "numbers"

    #expect(examine.itemsAsArray(sorted: .ascending) == [1, 1, 2, 3])
    #expect(examine.itemsAsArray(sorted: .descending) == [3, 2, 1, 1])
    let noneOrder = examine.itemsAsArray(sorted: .none).sorted()
    #expect(noneOrder == [1, 1, 2, 3])

    #expect(examine.uniqueItems(sorted: .ascending) == [1, 2, 3])
    #expect(examine.uniqueItems(sorted: .descending) == [3, 2, 1])
    let uniqueNone = examine.uniqueItems(sorted: .none)
    #expect(Set(uniqueNone) == Set([1, 2, 3]))

    let row = examine.itemsAsString(delimiter: "|", asRow: true)
    #expect(row == "3|1|2|1")
    let column = examine.itemsAsString(delimiter: nil, asRow: false, encloseElementsBy: "'")
    #expect(column == "numbers\n'3'\n'1'\n'2'\n'1'")
}

@Test("SSExamine empirical CDF matches cumulative frequencies")
func ssExamineEmpiricalCDF() throws {
    let examine = SSExamine<Int, Double>()
    #expect(examine.empiricalCDF(item: 42) == nil)

    #expect(examine.append(contentOf: [1, 1, 2, 3, 5]))

    let belowMinimum = examine.empiricalCDF(item: 0)
    #expect(belowMinimum == 0)

    let cdfAtOne = try #require(examine.empiricalCDF(item: 1))
    #expect(abs(Double(cdfAtOne) - 0.4) < 1e-12)

    let cdfAtTwo = try #require(examine.empiricalCDF(item: 2))
    #expect(abs(Double(cdfAtTwo) - 0.6) < 1e-12)

    let betweenThreeAndFive = try #require(examine.empiricalCDF(item: 4))
    #expect(abs(Double(betweenThreeAndFive) - 0.8) < 1e-12)

    let cdfAtFive = try #require(examine.empiricalCDF(item: 5))
    #expect(abs(Double(cdfAtFive) - 1.0) < 1e-12)

    let aboveMaximum = try #require(examine.empiricalCDF(item: 6))
    #expect(abs(Double(aboveMaximum) - 1.0) < 1e-12)
}

@Test("SSExamine produces frequency and cumulative tables")
func ssExamineFrequencyTables() {
    let examine = SSExamine<Int, Double>()
    #expect(examine.append(contentOf: [3, 1, 2, 1]))

    let expectedCounts: [Int: Int] = [1: 2, 2: 1, 3: 1]
    let tableByValue = examine.frequencyTable(sorted: .valueAscending)
    #expect(tableByValue.map(\.item) == [1, 2, 3])
    #expect(tableByValue.map(\.frequency) == [2, 1, 1])
    let relativeValues = tableByValue.map { Double($0.relativeFrequency) }
    let expectedRelatives: [Double] = [0.5, 0.25, 0.25]
    for (value, expected) in zip(relativeValues, expectedRelatives) {
        #expect(abs(value - expected) < 1e-12)
    }

    let tableByValueDesc = examine.frequencyTable(sorted: .valueDescending)
    #expect(tableByValueDesc.map(\.item) == [3, 2, 1])

    let tableByFrequencyAsc = examine.frequencyTable(sorted: .frequencyAscending)
    #expect(tableByFrequencyAsc.first?.frequency == 1)
    #expect(tableByFrequencyAsc.last?.frequency == 2)

    let tableByFrequencyDesc = examine.frequencyTable(sorted: .frequencyDescending)
    #expect(tableByFrequencyDesc.first?.frequency == 2)

    let tableNone = examine.frequencyTable(sorted: .none)
    let countsDict = Dictionary(uniqueKeysWithValues: tableNone.map { ($0.item, $0.frequency) })
    #expect(countsDict == expectedCounts)

    let cumulativeItem = examine.cumulativeFrequencyTable(format: .eachItem)
    #expect(cumulativeItem.count == 1)
    if let entry = cumulativeItem.first {
        #expect(entry.item == 1)
        #expect(entry.cumulativeCount == 2)
        #expect(abs(Double(entry.cumulativeRelativeCount) - 0.5) < 1e-12)
    }

    let cumulativeUnique = examine.cumulativeFrequencyTable(format: .eachUniqueItem)
    #expect(cumulativeUnique.count == 1)
    if let entry = cumulativeUnique.first {
        #expect(entry.item == 1)
        #expect(entry.cumulativeCount == 2)
        #expect(abs(Double(entry.cumulativeRelativeCount) - 0.5) < 1e-12)
    }
}

@Test("SSExamine cumulative caches and numeric projection track mutations")
func ssExamineCumulativeCachesAndNumericProjection() {
    let examine = SSExamine<Int, Double>()
    #expect(examine.cumulativeFrequencies.isEmpty)
    #expect(examine.cumulativeRelativeFrequencies.isEmpty)
    #expect(examine.itemsAsNumericArray?.isEmpty == true)

    #expect(examine.append(contentOf: [1, 2, 2, 3]))

    let cumulative = examine.cumulativeFrequencies
    #expect(cumulative == [1: 1, 2: 3, 3: 4])

    let cumulativeRelative = examine.cumulativeRelativeFrequencies
    #expect(cumulativeRelative.keys.count == 3)
    #expect(abs(Double(cumulativeRelative[1] ?? -1) - 0.25) < 1e-12)
    #expect(abs(Double(cumulativeRelative[2] ?? -1) - 0.75) < 1e-12)
    #expect(abs(Double(cumulativeRelative[3] ?? -1) - 1.0) < 1e-12)

    if let numericValues = examine.itemsAsNumericArray {
        #expect(numericValues.sorted() == [1.0, 2.0, 2.0, 3.0])
    } else {
        #expect(Bool(false), "Expected numeric projection for integer dataset")
    }

    examine.remove(2, allOccurences: true)

    let updatedCumulative = examine.cumulativeFrequencies
    #expect(updatedCumulative == [1: 1, 3: 2])

    let updatedRelative = examine.cumulativeRelativeFrequencies
    #expect(abs(Double(updatedRelative[1] ?? -1) - 0.5) < 1e-12)
    #expect(abs(Double(updatedRelative[3] ?? -1) - 1.0) < 1e-12)

    if let compactNumeric = examine.itemsAsNumericArray {
        #expect(Set(compactNumeric) == Set([1.0, 3.0]))
    } else {
        #expect(Bool(false), "Expected numeric projection after removal")
    }
}

@Test("SSExamine hash reflects raw sequence identity")
func ssExamineHashFollowsRawSequence() throws {
    let examine = SSExamine<Int, Double>()
    #expect(examine.hash == 0)

    #expect(examine.append(contentOf: [1, 2, 2]))
    let copy = try #require(examine.copy() as? SSExamine<Int, Double>)
    #expect(copy.hash == examine.hash)

    let differentOrder = SSExamine<Int, Double>()
    #expect(differentOrder.append(contentOf: [2, 1, 2]))
    #expect(differentOrder.hash != examine.hash)

    examine.append(3)
    copy.append(3)
    #expect(examine.hash == copy.hash)
}

@Test("SSExamine removal, contains, and subscript behaviours")
func ssExamineRemovalAndSubscript() {
    let examine = SSExamine<Int, Double>()
    #expect(examine.append(contentOf: [5, 6, 7]))

    #expect(examine[0] == 5)
    #expect(examine[2] == 7)
    #expect(examine[3] == nil)
    #expect(examine[-1] == nil)

    #expect(examine.contains(6))
    examine.hasChanges = false
    examine.remove(99, allOccurences: true)
    #expect(examine.itemsAsArray(sorted: .original) == [5, 6, 7])
    #expect(examine.hasChanges)

    examine.hasChanges = false
    examine.remove(6, allOccurences: false)
    #expect(examine.itemsAsArray(sorted: .original) == [5, 7])
    #expect(examine.contains(6) == false)
    #expect(examine.sampleSize == 2)
    #expect(examine.length == 2)

    examine.removeAll()
    #expect(examine.isEmpty)
    #expect(examine.sampleSize == 0)
    #expect(examine.length == 0)
    #expect(examine.hasChanges)
    #expect(examine.relativeFrequency(5) == 0)
    #expect(examine.itemsAsArray(sorted: .ascending).isEmpty)
}

@Suite("SSExamine Stats: totals and powers over array sizes 2...51")
struct SSExamineStatsTests {

    // Helper for magnitude-aware comparisons
    private func approxEqual(_ a: Double, _ b: Double, ulpScaledEps: Double = 1e-12, minAbsEps: Double = 1e-12) -> Bool {
        let scale = max(1.0, max(abs(a), abs(b)))
        return abs(a - b) <= max(minAbsEps, ulpScaledEps * scale)
    }

    // Generates strictly positive non-zero values to avoid inverse/odd-power issues.
    private func positiveValues(count n: Int) -> [Double] {
        // 0.5, 1.5, 2.5, ... ensures no zeros and all positive
        return (0..<n).map { Double($0) + 0.5 }
    }

    // Generates alternating sign values with no zeros to stress compensation paths.
    private func alternatingValues(count n: Int) -> [Double] {
        // ±(k + 1.0) to avoid zeros, starting negative for variety
        return (0..<n).map { i in
            let base = Double(i + 1)
            return (i % 2 == 0) ? -base : base
        }
    }

    @Test("Totals and means for positive sequences (2...51)")
    func totalsAndMeansPositive() {
        for n in 2...51 {
            let values = positiveValues(count: n)
            let expectedTotal = values.reduce(0, +)
            let expectedSquares = values.reduce(0) { $0 + $1 * $1 }
            let expectedMean = expectedTotal / Double(n)

            let ex = SSExamine<Double, Double>()
            #expect(ex.append(contentOf: values))

            guard
                let totalOpt = ex.total,
                let squareTotalOpt = ex.squareTotal,
                let p2Opt = ex.poweredTotal(power: 2),
                let meanOpt = ex.arithmeticMean,
                let p0Opt = ex.poweredTotal(power: 0),
                let p1Opt = ex.poweredTotal(power: 1)
            else {
                #expect(Bool(false), "Expected non-nil totals for numeric dataset of size \(n)")
                return
            }

            let total = Double(totalOpt)
            let squareTotal = Double(squareTotalOpt)
            let p2 = Double(p2Opt)
            let mean = Double(meanOpt)

            #expect(approxEqual(total, expectedTotal))
            #expect(approxEqual(squareTotal, expectedSquares))
            #expect(approxEqual(p2, expectedSquares))
            #expect(approxEqual(mean, expectedMean))

            // Fast paths
            let p0 = Double(p0Opt)
            let p1 = Double(p1Opt)
            #expect(p0 == Double(n))
            #expect(approxEqual(p1, expectedTotal))
        }
    }

    @Test("Totals and powers for alternating-sign sequences (2...51)")
    func totalsAndPowersAlternating() {
        for n in 2...51 {
            let values = alternatingValues(count: n)
            let expectedTotal = values.reduce(0, +)
            let expectedSquares = values.reduce(0) { $0 + $1 * $1 }

            let ex = SSExamine<Double, Double>()
            #expect(ex.append(contentOf: values))

            // total and p1
            guard
                let totalOpt = ex.total,
                let p1Opt = ex.poweredTotal(power: 1),
                let squareTotalOpt = ex.squareTotal,
                let p2Opt = ex.poweredTotal(power: 2),
                let p0Opt = ex.poweredTotal(power: 0)
            else {
                #expect(Bool(false), "Expected non-nil totals for alternating dataset of size \(n)")
                return
            }

            let total = Double(totalOpt)
            let p1 = Double(p1Opt)
            #expect(approxEqual(total, expectedTotal))
            #expect(approxEqual(p1, expectedTotal))

            // squares
            let squareTotal = Double(squareTotalOpt)
            let p2 = Double(p2Opt)
            #expect(approxEqual(squareTotal, expectedSquares))
            #expect(approxEqual(p2, expectedSquares))

            // p0
            let p0 = Double(p0Opt)
            #expect(p0 == Double(n))
        }
    }

    @Test("Powered totals with non-integer power on positive data (2...51)")
    func poweredTotalsNonIntegerPowerPositive() {
        let power: Double = 1.5
        for n in 2...51 {
            let values = positiveValues(count: n)
            let expected = values.reduce(0.0) { $0 + pow($1, power) }

            let ex = SSExamine<Double, Double>()
            #expect(ex.append(contentOf: values))

            guard let pOpt = ex.poweredTotal(power: power) else {
                #expect(Bool(false), "Expected non-nil powered total (\(power)) for dataset size \(n)")
                return
            }
            let p = Double(pOpt)
            #expect(approxEqual(p, expected, ulpScaledEps: 1e-11, minAbsEps: 1e-12))
        }
    }

    @Test("Inverse totals for positive non-zero data (2...51)")
    func inverseTotalsPositive() {
        for n in 2...51 {
            let values = positiveValues(count: n) // strictly > 0
            let expected = values.reduce(0.0) { $0 + 1.0 / $1 }

            let ex = SSExamine<Double, Double>()
            #expect(ex.append(contentOf: values))

            guard let invOpt = ex.inverseTotal else {
                #expect(Bool(false), "Expected non-nil inverse total for dataset size \(n)")
                return
            }
            let inv = Double(invOpt)
            #expect(approxEqual(inv, expected, ulpScaledEps: 1e-11, minAbsEps: 1e-12))
        }
    }

    @Test("Consistency: squareTotal equals poweredTotal(2) (2...51)")
    func squareTotalConsistency() {
        for n in 2...51 {
            // Use a mix to exercise both code paths and numeric stability
            let values = (n % 2 == 0) ? positiveValues(count: n) : alternatingValues(count: n)

            let ex = SSExamine<Double, Double>()
            #expect(ex.append(contentOf: values))

            guard
                let squareTotalOpt = ex.squareTotal,
                let p2Opt = ex.poweredTotal(power: 2)
            else {
                #expect(Bool(false), "Expected non-nil square totals for dataset size \(n)")
                return
            }

            let s = Double(squareTotalOpt)
            let p2 = Double(p2Opt)
            #expect(approxEqual(s, p2))
        }
    }
}
