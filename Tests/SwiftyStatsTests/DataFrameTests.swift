/// Tests covering SSDataFrame creation, mutation, and derived statistics.

import Foundation
import Testing
@testable import SwiftyStats

struct DataFrameTests {
    private func makeExamine(_ values: [Double], name: String? = nil, tag: String? = nil) -> SSExamine<Double, Double> {
        let ex = SSExamine<Double, Double>()
        _ = ex.append(contentOf: values)
        ex.name = name
        ex.tag = tag
        return ex
    }

    @Test("Initialization assigns unique column names and preserves counts")
    func initAndNaming() throws {
        let colA = makeExamine([1.0, 2.0, 3.0], name: "colA")
        let colADuplicate = makeExamine([4.0, 5.0], name: "colA")

        let df = try DataFrame(data: [colA, colADuplicate])

        #expect(df.nCols == 2)
        #expect(df.isEmpty == false)
        #expect(df.columnNames[0] == "colA")
        #expect(df.columnNames[1].hasPrefix("colA"))
        #expect(df.columnNames[1] != "colA") // uniqueness enforced
        #expect(df.sampleSize(column: 0) == 3)
        #expect(df.sampleSize(column: 1) == 2)
        #expect(df.totalSampleSize == 5)
        #expect(df[columnName: "does-not-exist"] == nil)
    }

    @Test("Append, remove, and subscript behave consistently")
    func appendAndRemove() throws {
        let df = DataFrame<Double, Double>()

        let first = makeExamine([1.0, 1.0], name: "first", tag: "t1")
        try df.append(first, name: nil)
        #expect(df.nCols == 1)
        #expect(df.columnNames[0] == "first")
        #expect(df.tags[0] == "t1")

        let second = makeExamine([2.0], name: "first", tag: "t2")
        try df.append(second, name: "custom")
        #expect(df.nCols == 2)
        #expect(df.columnNames == ["first", "custom"])
        #expect(df.tags == ["t1", "t2"])
        #expect(df[1].sampleSize == 1)

        let removed = df.remove(name: "first")
        #expect(removed?.sampleSize == 2)
        #expect(df.nCols == 1)
        #expect(df.columnNames == ["custom"])

        df.removeAll()
        #expect(df.isEmpty)
        #expect(df.nCols == 0)
    }

    @Test("Copying and Codable round-trips preserve structure")
    func copyAndCodable() throws {
        let c1 = makeExamine([1.0, 2.0], name: "x", tag: "tx")
        let c2 = makeExamine([3.0], name: "y", tag: "ty")
        let original = try DataFrame(data: [c1, c2])

        let copied = try #require(original.copy() as? DataFrame<Double, Double>)
        #expect(copied.columnNames == original.columnNames)
        #expect(copied.tags == original.tags)
        #expect(copied.nCols == original.nCols)

        let encoded = try JSONEncoder().encode(original)
        let decoded = try JSONDecoder().decode(DataFrame<Double, Double>.self, from: encoded)
        #expect(decoded.columnNames == original.columnNames)
        #expect(decoded.tags == original.tags)
        #expect(decoded.nCols == original.nCols)
        #expect(decoded.sampleSize(column: 0) == 2)
        #expect(decoded.sampleSize(column: 1) == 1)
    }

    @Test("saveTo and dataframe(from:) round-trip JSON and overwrite semantics")
    func saveToRoundTrip() throws {
        let fm = FileManager.default
        let tmpRoot = fm.temporaryDirectory.appendingPathComponent("SwiftyStats-DataFrame-\(UUID().uuidString)")
        try fm.createDirectory(at: tmpRoot, withIntermediateDirectories: true)
        defer { try? fm.removeItem(at: tmpRoot) }

        let path = tmpRoot.appendingPathComponent("frame.json").path
        let df = try DataFrame(data: [
            makeExamine([1.0, 2.0], name: "c1", tag: "t1"),
            makeExamine([3.0], name: "c2", tag: "t2")
        ])

        #expect(try df.saveTo(fileName: path, overwrite: true))
        let loadedOpt = try DataFrame<Double, Double>.dataframe(fromFile: path)
        let loaded = try #require(loadedOpt)
        #expect(loaded.nCols == 2)
        #expect(loaded.columnNames == df.columnNames)
        #expect(loaded.tags == df.tags)
        #expect(loaded.sampleSize(column: 0) == 2)
        #expect(loaded.sampleSize(column: 1) == 1)

        #expect(throws: SSError.self) {
            _ = try df.saveTo(fileName: path, overwrite: false)
        }

        #expect(try df.saveTo(fileName: path, overwrite: true))
    }

    @Test("saveTo supports ZIP compression and dataframe(from:) inflates it")
    func saveToCompressedRoundTrip() throws {
        let fm = FileManager.default
        let tmpRoot = fm.temporaryDirectory.appendingPathComponent("SwiftyStats-DataFrame-ZIP-\(UUID().uuidString)")
        try fm.createDirectory(at: tmpRoot, withIntermediateDirectories: true)
        defer { try? fm.removeItem(at: tmpRoot) }

        let path = tmpRoot.appendingPathComponent("frame.zip").path
        let df = try DataFrame(data: [
            makeExamine([10.0, 20.0, 30.0], name: "zip1", tag: "tz1"),
            makeExamine([5.0], name: "zip2", tag: "tz2")
        ])

        #expect(try df.saveTo(fileName: path, overwrite: true, compressAsZip: true, zipEntryFileName: "data.json"))
        let loadedOpt = try DataFrame<Double, Double>.dataframe(fromFile: path)
        let loaded = try #require(loadedOpt)
        #expect(loaded.nCols == 2)
        #expect(loaded.columnNames == df.columnNames)
        #expect(loaded.tags == df.tags)
        #expect(loaded.sampleSize(column: 0) == 3)
        #expect(loaded.sampleSize(column: 1) == 1)
    }
}
