/// Tests ensuring stable product calculations and sign handling in SSExamine.

import Testing
import SwiftyStats
import Foundation

@Suite("Product calculations across methods")
struct ProductTests {

    // Helper to build an SSExamine from an array of FP values.
    private func makeExamine(_ values: [Double]) throws -> SSExamine<Double, Double> {
        return try SSExamine(using: values, levelOfMeasurement: .ratio, name: nil, characterSet: nil)
    }

    @Test("Logarithmic method handles negatives via abs-logs + sign")
    func logarithmicHandlesNegatives() async throws {
        // Product: (-2) * 3 * (-4) = 24
        let values = [-2.0, 3.0, -4.0]
        let ex = try makeExamine(values)

        let res = try #require(ex.product(method: .logarithmic))
        #expect(res.value.isFinite)
        #expect(res.logValue != nil)
        #expect(abs(res.value - 24.0) <= 1e-12)
        // log(24)
        #expect(abs((res.logValue ?? .nan) - log(24.0)) <= 1e-12)
        #expect(res.overflowRisk == false || res.overflowRisk == nil)
        #expect(res.underflowRisk == false || res.underflowRisk == nil)
    }

    @Test("Zero short-circuit preserves signed zero")
    func signedZeroWhenZeroPresent() async throws {
        // Overall sign negative (one negative), contains zero → -0.0
        let values = [-2.0, 0.0, 5.0]
        let ex = try makeExamine(values)

        // Check multiple methods for consistency
        for method in [ProductMethod.naive, .logarithmic, .scaled, .exponentMantissa, .compensated] {
            let res = try #require(ex.product(method: method))
            #expect(res.value == 0.0)
            // Signed zero check: 1/(-0.0) is -∞
            let inv = 1.0 / res.value
            #expect(inv.isInfinite && inv.sign == .minus)
        }
    }

    @Test("Logarithmic overflow produces signed infinity with overflowRisk")
    func logarithmicOverflow() async throws {
        // sqrt(max)^3 overflows
        let s = Double.greatestFiniteMagnitude.squareRoot()
        let values = [s, s, s] // positive sign
        let ex = try makeExamine(values)

        let res = try #require(ex.product(method: .logarithmic))
        #expect(res.value.isInfinite && res.value.sign == .plus)
        #expect(res.overflowRisk == true)
        #expect(res.underflowRisk == false || res.underflowRisk == nil)
    }

    @Test("Logarithmic underflow to signed zero with underflowRisk")
    func logarithmicUnderflowSignedZero() async throws {
        // Smallest positive subnormal squared is 0. Ensure odd negatives → -0.0
        let a = Double.leastNonzeroMagnitude
        let values = [-a, a] // product is -a^2 -> definitely underflows to -0.0
        let ex = try makeExamine(values)

        let res = try #require(ex.product(method: .logarithmic))
        #expect(res.value == 0.0)
        let inv = 1.0 / res.value
        #expect(inv.isInfinite && inv.sign == .minus) // -0.0
        #expect(res.underflowRisk == true)
    }

    @Test("Agreement across methods on moderate range")
    func agreementModerateRange() async throws {
        // Product: 1.5 * 2.0 * 0.5 * (-3.0) * (-4.0) = 18.0
        let values = [1.5, 2.0, 0.5, -3.0, -4.0]
        let ex = try makeExamine(values)

        var results: [Double] = []
        for method in [ProductMethod.naive, .logarithmic, .scaled, .exponentMantissa, .compensated] {
            let r = try #require(ex.product(method: method))
            results.append(r.value)
        }

        // All methods should be very close to 18
        for v in results {
            #expect(abs(v - 18.0) <= 1e-10)
        }
    }

    @Test("ExponentMantissa returns exponent and mantissa")
    func exponentMantissaFieldsPresent() async throws {
        let values = [1.25, 2.5, 3.5]
        let ex = try makeExamine(values)

        let res = try #require(ex.product(method: .exponentMantissa))
        // Value should be finite
        #expect(res.value.isFinite)
        // Exponent/mantissa should be populated for scaled/exponentMantissa
        #expect(res.exponent != nil)
        #expect(res.mantissa != nil)
        // Reconstruct from mantissa * 2^exponent (sign included in mantissa via product method)
        if let m = res.mantissa, let e = res.exponent {
            let recon = m * Double.pow(2.0, Double(e))
            // Allow some rounding error from normalization steps
            #expect(abs(recon - abs(res.value)) <= 1e-12 || abs(recon - res.value) <= 1e-12)
        }
    }
}
