/// Thread-safety and concurrent mutation tests for SSExamine.

import Foundation
import Testing
@testable import SwiftyStats

@Test("SSExamine concurrent append maintains integrity")
func ssExamineConcurrentAccess() async throws {
    let examine = SSExamine<Int, Double>()
    let taskCount = 32
    let valuesPerTask = 100
    try await withThrowingTaskGroup(of: Void.self) { group in
        for taskID in 0..<taskCount {
            group.addTask {
                let base = taskID * valuesPerTask
                for offset in 0..<valuesPerTask {
                    examine.append(base + offset)
                }
            }
        }
        try await group.waitForAll()
    }

    let expectedSize = taskCount * valuesPerTask
    #expect(examine.sampleSize == expectedSize)
    for value in 0..<expectedSize {
        #expect(examine.frequency(value) == 1)
    }
}

@Test("SSExamine removal and encoding remain consistent under load")
func ssExamineRemovalAndEncodingStress() async throws {
    let examine = SSExamine<Int, Double>()
    let totalValues = 1_000
    let chunkSize = 100

    try await withThrowingTaskGroup(of: Void.self) { group in
        for offset in stride(from: 0, to: totalValues, by: chunkSize) {
            group.addTask {
                let upperBound = min(offset + chunkSize, totalValues)
                for value in offset..<upperBound {
                    examine.append(value)
                }
            }
        }
        try await group.waitForAll()
    }

    #expect(examine.sampleSize == totalValues)

    try await withThrowingTaskGroup(of: Void.self) { group in
        group.addTask {
            for value in 0..<(totalValues / 2) {
                examine.remove(value, allOccurences: true)
            }
        }

        group.addTask {
            let encoder = JSONEncoder()
            let decoder = JSONDecoder()
            for _ in 0..<75 {
                let data = try encoder.encode(examine)
                let decoded = try decoder.decode(SSExamine<Int, Double>.self, from: data)
                let decodedSize = decoded.sampleSize
                #expect(decodedSize <= totalValues)
                #expect(decodedSize >= totalValues / 2)
            }
        }

        try await group.waitForAll()
    }

    #expect(examine.sampleSize == totalValues / 2)
    for value in 0..<(totalValues / 2) {
        #expect(examine.contains(value) == false)
    }
    for value in (totalValues / 2)..<totalValues {
        #expect(examine.contains(value) == true)
    }
}
