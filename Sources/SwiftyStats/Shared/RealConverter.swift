//
//  Created by VT on 22.12.25.
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
import Foundation

/// Errors thrown during conversion in `RealConverter`.
public enum RealConversionError: Error, CustomStringConvertible {
    /// The input type is not supported by the converter.
    case unsupportedType(Any.Type)
    /// The provided string cannot be parsed as a finite number.
    case invalidString(String)
    /// The numeric value is out of the representable range for the target type.
    case outOfRange(String)
    /// The value represents a non-finite number (NaN or ±Infinity).
    case nonFinite(String)
    /// A generic conversion failure with context.
    case conversionFailed(String)

    /// Human-readable description of the error.
    public var description: String {
        switch self {
        case .unsupportedType(let t): return "Unsupported input type: " + String(describing: t)
        case .invalidString(let s):   return "Invalid numeric representation: \"" + s + "\""
        case .outOfRange(let s):      return "Value out of representable range: " + s
        case .nonFinite(let s):       return "Non-finite value (NaN/Inf): " + s
        case .conversionFailed(let s):return "Conversion failed: " + s
        }
    }
}

/// A locale-aware, robust converter for turning heterogenous values into concrete
/// floating-point types (`Float`, `Double`, and `Float80` on Intel).
public final class RealConverter {

    /// Convert an arbitrary value to a concrete floating-point type `T`.
    ///
    /// Supported inputs include numeric primitives, `NSNumber`, `Decimal`, `String`
    /// (with locale-aware parsing), and `CGFloat` (where available).
    ///
    /// - Parameters:
    ///   - value: The source value to convert.
    ///   - locale: Optional locale hint used for parsing strings; falls back to
    ///             reasonable candidates (autoupdating, current, de_DE, en_US_POSIX).
    /// - Returns: The converted value as type `T`.
    /// - Throws: `RealConversionError` if conversion fails or yields a non-finite value.
    public static func to<T: RealLike>(_ value: Any, locale: Locale? = nil) throws -> T {
        if let t = value as? T { return t }

        switch value {
        case let s as String:
            return try parseString(s, as: T.self, localeHint: locale)

        case let dec as Decimal:
            return try fromDecimal(dec, as: T.self)

        case let num as NSNumber:
            // Treat NSNumber(boolean) as 0/1
            if CFGetTypeID(num) == CFBooleanGetTypeID() {
                return try fromInteger(num.boolValue ? 1 : 0, as: T.self)
            }
            let dec = num.decimalValue
            return try fromDecimal(dec, as: T.self)

        case let d as Double: return try fromFP(d, as: T.self)
        case let f as Float:  return try fromFP(f, as: T.self)
        #if arch(i386) || arch(x86_64)
        case let f80 as Float80: return try fromFP(f80, as: T.self)
        #endif

        #if canImport(CoreGraphics)
        case let cg as CGFloat:
            return try fromFP(Double(cg), as: T.self)
        #endif

        case let i as Int:     return try fromInteger(i, as: T.self)
        case let i as Int8:    return try fromInteger(i, as: T.self)
        case let i as Int16:   return try fromInteger(i, as: T.self)
        case let i as Int32:   return try fromInteger(i, as: T.self)
        case let i as Int64:   return try fromInteger(i, as: T.self)
        case let u as UInt:    return try fromInteger(u, as: T.self)
        case let u as UInt8:   return try fromInteger(u, as: T.self)
        case let u as UInt16:  return try fromInteger(u, as: T.self)
        case let u as UInt32:  return try fromInteger(u, as: T.self)
        case let u as UInt64:  return try fromInteger(u, as: T.self)

        default:
            throw RealConversionError.unsupportedType(type(of: value))
        }
    }

    /// Convenience: Convert to `Float`.
    @inline(__always) public static func toFloat(_ v: Any, locale: Locale? = nil) throws -> Float  { try to(v, locale: locale) as Float }
    /// Convenience: Convert to `Double`.
    @inline(__always) public static func toDouble(_ v: Any, locale: Locale? = nil) throws -> Double { try to(v, locale: locale) as Double }
    #if arch(i386) || arch(x86_64)
    /// Convenience: Convert to `Float80` (Intel architectures only).
    @inline(__always) public static func toFloat80(_ v: Any, locale: Locale? = nil) throws -> Float80 { try to(v, locale: locale) as Float80 }
    #endif
}

// MARK: - Private

private extension RealConverter {

    /// Build a list of locale candidates to try when parsing strings.
    static func localeCandidates(from hint: Locale?) -> [Locale] {
        var list: [Locale] = []
        if let h = hint { list.append(h) }
        list.append(Locale.autoupdatingCurrent)
        let curr = Locale.current
        if curr.identifier != Locale.autoupdatingCurrent.identifier { list.append(curr) }
        list.append(contentsOf: [Locale(identifier: "de_DE"), Locale(identifier: "en_US_POSIX")])
        var seen = Set<String>()
        return list.filter { seen.insert($0.identifier).inserted }
    }

    /// Parse a `String` into `T` using multiple locale strategies.
    /// On Intel, attempts to preserve Float80 precision via canonicalization and `strtold`.
    static func parseString<T: RealLike>(_ s: String, as _: T.Type, localeHint: Locale?) throws -> T {
        let trimmed = s.trimmingCharacters(in: .whitespacesAndNewlines)

        switch trimmed.lowercased() {
        case "nan": throw RealConversionError.nonFinite("NaN")
        case "inf", "+inf": throw RealConversionError.nonFinite("+Inf")
        case "-inf": throw RealConversionError.nonFinite("-Inf")
        default: break
        }

        // Intel-only: try canonicalization + strtold for Float80 first.
        #if arch(x86_64) || arch(i386)
        if T.self == Float80.self {
            for loc in localeCandidates(from: localeHint) {
                if let canon = canonicalizeNumericString(trimmed, locale: loc),
                   let f80 = parseFloat80CanonicalString(canon),
                   f80.isFinite {
                    return f80 as! T
                }
            }
            if let f80 = parseFloat80CanonicalString(trimmed), f80.isFinite {
                return f80 as! T
            }
        }
        #endif

        // 1) NumberFormatter
        for loc in localeCandidates(from: localeHint) {
            if let d = parseWithFormatter(trimmed, locale: loc) {
                return try fromFP(d, as: T.self)
            }
        }

        // 2) Decimal -> T
        if let dec = Decimal(string: trimmed) {
            return try fromDecimal(dec, as: T.self)
        }

        // 3) Plain Double
        if let d = Double(trimmed) {
            return try fromFP(d, as: T.self)
        }

        // 4) Heuristic
        let runtime = localeCandidates(from: localeHint).first ?? .autoupdatingCurrent
        if let d = coerceLocaleAwareDecimal(trimmed, locale: runtime) {
            return try fromFP(d, as: T.self)
        }

        throw RealConversionError.invalidString(s)
    }

    /// Parse using `NumberFormatter` after normalizing common grouping characters.
    static func parseWithFormatter(_ s: String, locale: Locale) -> Double? {
        let nf = NumberFormatter()
        nf.locale = locale
        nf.numberStyle = .decimal
        let grpSep = nf.groupingSeparator ?? ","
        var cleaned = s
            .replacingOccurrences(of: "\u{00A0}", with: " ")
            .replacingOccurrences(of: "\u{202F}", with: " ")
            .replacingOccurrences(of: "\u{2009}", with: " ")
        let quoteLikes = ["'", "\\u{2019}", "\\u{201B}", "\\u{02BC}"]
        for q in quoteLikes {
            cleaned = cleaned.replacingOccurrences(of: q, with: grpSep)
        }
        return nf.number(from: cleaned)?.doubleValue
    }

    /// Heuristically coerce a localized decimal string to a `Double` using likely separators.
    static func coerceLocaleAwareDecimal(_ s: String, locale: Locale) -> Double? {
        let nf = NumberFormatter()
        nf.locale = locale
        nf.numberStyle = .decimal
        let decSep = nf.decimalSeparator ?? "."
        let grpSep = nf.groupingSeparator ?? ","

        var t = s.replacingOccurrences(of: "\u{00A0}", with: "")
                  .replacingOccurrences(of: "\u{202F}", with: "")
                  .replacingOccurrences(of: " ", with: "")

        if t.contains(decSep) && t.contains(grpSep) {
            t = t.replacingOccurrences(of: grpSep, with: "")
                 .replacingOccurrences(of: decSep, with: ".")
            return Double(t)
        }

        if t.contains(decSep) {
            t = t.replacingOccurrences(of: decSep, with: ".")
            return Double(t)
        }

        let altPairs: [(grp: String, dec: String)] = [
            (",", "."), (".", ","),
            ("\\u{2019}", ","), ("\\u{2019}", "."),
            ("'", ","), ("'", "."),
            ("\u{00A0}", ","), ("\u{00A0}", "."),
            ("\u{202F}", ","), ("\u{202F}", "."),
            ("\u{2009}", ","), ("\u{2009}", ".")
        ]
        for (g, d) in altPairs where t.contains(g) || t.contains(d) {
            let u = t.replacingOccurrences(of: g, with: "")
                     .replacingOccurrences(of: d, with: ".")
            if let val = Double(u) { return val }
        }
        return nil
    }

    /// Convert `Decimal` to `T`. On Intel, preserves `Float80` via `strtold` on a canonical string.
    static func fromDecimal<T: RealLike>(_ dec: Decimal, as _: T.Type) throws -> T {
        #if arch(x86_64) || arch(i386)
        if T.self == Float80.self {
            let s = withUnsafePointer(to: dec) { ptr in
                NSDecimalString(ptr, Locale(identifier: "en_US_POSIX"))
            }
            if let f80 = parseFloat80CanonicalString(s), f80.isFinite {
                return f80 as! T
            }
            throw RealConversionError.nonFinite(s)
        }
        #endif
        let d = (dec as NSDecimalNumber).doubleValue
        return try fromFP(d, as: T.self)
    }

    /// Convert an integer value to `T`, validating finiteness.
    static func fromInteger<T: RealLike, I: BinaryInteger>(_ i: I, as _: T.Type) throws -> T {
        let t = T(i)
        if t.isFinite { return t }
        throw RealConversionError.outOfRange("\(i)")
    }

    /// Convert a floating-point value to `T`. On Intel, preserves `Float80` precision.
    static func fromFP<T: RealLike, F: BinaryFloatingPoint>(_ x: F, as _: T.Type) throws -> T {
        #if arch(x86_64) || arch(i386)
        if T.self == Float80.self {
            let f80: Float80
            if let v = x as? Float80 {
                f80 = v
            } else if let v = x as? Double {
                f80 = Float80(v)
            } else if let v = x as? Float {
                f80 = Float80(v)
            } else {
                f80 = Float80(Double(x))
            }
            if f80.isFinite {
                return f80 as! T
            }
            throw RealConversionError.nonFinite("\(x)")
        }
        #endif

        let d = Double(x)
        let t = T(d)
        if t.isFinite { return t }
        throw RealConversionError.nonFinite("\(x)")
    }

    // MARK: Float80 helpers (Intel)

    /// Produces a canonical C-locale numeric string ('.' decimal, no grouping) from a localized representation.
    #if arch(x86_64) || arch(i386)
    static func canonicalizeNumericString(_ s: String, locale: Locale) -> String? {
        let nf = NumberFormatter()
        nf.locale = locale
        nf.numberStyle = .decimal
        let decSep = nf.decimalSeparator ?? "."
        let grpSep = nf.groupingSeparator ?? ","

        var t = s.trimmingCharacters(in: .whitespacesAndNewlines)
        t = t
            .replacingOccurrences(of: "\u{00A0}", with: " ")
            .replacingOccurrences(of: "\u{202F}", with: " ")
            .replacingOccurrences(of: "\u{2009}", with: " ")

        let quoteLikes = ["'", "\\u{2019}", "\\u{201B}", "\\u{02BC}"]
        for q in quoteLikes {
            t = t.replacingOccurrences(of: q, with: grpSep)
        }

        t = t.replacingOccurrences(of: " ", with: "")

        if t.contains(decSep) && t.contains(grpSep) {
            t = t.replacingOccurrences(of: grpSep, with: "")
                 .replacingOccurrences(of: decSep, with: ".")
        } else if t.contains(decSep) {
            t = t.replacingOccurrences(of: decSep, with: ".")
        } else if t.contains(grpSep) {
            t = t.replacingOccurrences(of: grpSep, with: "")
        }

        // Handle swapped notation if present
        if t.contains(",") || t.contains("'") || t.contains("\\u{2019}") {
            let candidates: [(grp: String, dec: String)] = [
                (",", "."), (".", ","),
                ("\\u{2019}", ","), ("\\u{2019}", "."),
                ("'", ","), ("'", ".")
            ]
            for (g, d) in candidates where t.contains(g) || t.contains(d) {
                let u = t.replacingOccurrences(of: g, with: "")
                         .replacingOccurrences(of: d, with: ".")
                if isCanonicalNumericCandidate(u) {
                    t = u
                    break
                }
            }
        }

        guard isCanonicalNumericCandidate(t) else { return nil }
        return t
    }

    /// Quick validation that a string looks like a C-locale numeric literal.
    static func isCanonicalNumericCandidate(_ s: String) -> Bool {
        if s.isEmpty { return false }
        let allowed = CharacterSet(charactersIn: "+-0123456789.eE")
        if s.unicodeScalars.contains(where: { !allowed.contains($0) }) {
            return false
        }
        let dotCount = s.filter { $0 == "." }.count
        if dotCount > 1 { return false }
        let eCount = s.filter { $0 == "e" || $0 == "E" }.count
        if eCount > 1 { return false }
        if let eIndex = s.firstIndex(where: { $0 == "e" || $0 == "E" }) {
            let before = s[..<eIndex]
            let after = s[s.index(after: eIndex)...]
            if before.isEmpty || !before.contains(where: { $0.isNumber }) {
                return false
            }
            let afterStr = String(after)
            let afterBody = afterStr.drop(while: { $0 == "+" || $0 == "-" })
            if afterBody.isEmpty || !afterBody.contains(where: { $0.isNumber }) {
                return false
            }
        }
        return true
    }

    /// Parse canonical numeric string to Float80 via strtold.
    static func parseFloat80CanonicalString(_ s: String) -> Float80? {
        var result: Float80 = 0
        s.withCString { cstr in
            var endPtr: UnsafeMutablePointer<CChar>? = nil
            errno = 0
            let val = strtold(cstr, &endPtr)
            if let end = endPtr, end == UnsafeMutablePointer(mutating: cstr) {
                result = Float80.nan
                return
            }
            if let end = endPtr {
                let remainder = String(cString: end)
                if !remainder.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
                    result = Float80.nan
                    return
                }
            }
            result = val
        }
        guard result.isFinite else { return nil }
        return result
    }
    #endif
}

