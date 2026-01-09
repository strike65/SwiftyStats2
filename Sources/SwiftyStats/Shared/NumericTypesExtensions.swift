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
#if canImport(CoreGraphics)
// import CoreGraphics
#endif

#if canImport(RealModule)
import RealModule
/// A “Real-like” floating-point protocol used by the converter and helpers.
///
/// On platforms where RealModule is available, this refines both `Real` and
/// `BinaryFloatingPoint`. Otherwise it falls back to `BinaryFloatingPoint` only.
///
/// Conforming types:
/// - Float
/// - Double
/// - Float80 (where available)
public protocol RealLike: Real & BinaryFloatingPoint & Sendable & Comparable & Codable & RuntimeNumeric where Self.RawSignificand: FixedWidthInteger {}
extension Float: RealLike {}
extension Double: RealLike {}
#if arch(i386) || arch(x86_64)
/// Add Codable protocol conformance to Float80
extension Float80: @retroactive Decodable {}
extension Float80: @retroactive Encodable {}
extension Float80: RealLike {
    public init(from decoder: any Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        let stringRep: String = try container.decode(String.self, forKey: .string)
        if let f = Float80.init(stringRep) {
            self = f
        }
        else {
            self = Float80.nan
        }
    }
    
    private enum CodingKeys: String, CodingKey {
        case string = "stringValue"
    }
    public func encode(to encoder: any Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        let stringRep: String = "\(self)"
        try container.encode(stringRep, forKey: .string)
    }
}
#endif
#endif
// #else
#if !canImport(RealModule)
/// A “Real-like” floating-point protocol used by the converter and helpers.
///
/// When RealModule is not available, this is simply an alias for
/// `BinaryFloatingPoint` with the concrete conformances used in this project.
///
/// Conforming types:
/// - Float
/// - Double
/// - Float80 (on Intel architectures)
public protocol RealLike: BinaryFloatingPoint & Sendable & Comparable & Codable & RuntimeNumeric where Self.RawSignificand: FixedWidthInteger {}
extension Float: RealLike {}
extension Double: RealLike {}
#if arch(i386) || arch(x86_64)
extension Float80: RealLike {
    public init(from decoder: any Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        let stringRep: String = try container.decode(String.self, forKey: .string)
        if let f = Float80.init(stringRep) {
            self = f
        }
        else {
            self = Float80.nan
        }
    }
    
    private enum CodingKeys: String, CodingKey {
        case string = "stringValue"
    }
    public func encode(to encoder: any Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        let stringRep: String = "\(self)"
        try container.encode(stringRep, forKey: .string)
    }
}
#endif
#endif

/// Marker protocol identifying numeric types recognised by runtime checks.
public protocol RuntimeNumeric  {}

extension    Int: RuntimeNumeric {}
extension   Int8: RuntimeNumeric {}
extension  Int16: RuntimeNumeric {}
extension  Int32: RuntimeNumeric {}
extension  Int64: RuntimeNumeric {}
extension   UInt: RuntimeNumeric {}
extension  UInt8: RuntimeNumeric {}
extension UInt16: RuntimeNumeric {}
extension UInt32: RuntimeNumeric {}
extension UInt64: RuntimeNumeric {}
extension  Float: RuntimeNumeric {}
extension Double: RuntimeNumeric {}
#if canImport(CoreGraphics)
extension CGFloat: RuntimeNumeric {}
#endif
#if arch(x86_64) || arch(i386)
extension Float80: RuntimeNumeric {}
#endif
extension NSNumber: RuntimeNumeric {}
extension Decimal: RuntimeNumeric {}

/// Range convenience methods for all numeric types recognized at runtime.
///
/// These helpers are available for any type that is both `Comparable` and `RuntimeNumeric`
/// (i.e., all standard integer/unsigned types, Float/Double/Float80, CGFloat, Decimal, NSNumber).
internal extension Comparable where Self: RuntimeNumeric {
    
    /// Returns true if `self` lies in the open interval `(from, to)`.
    func liesInOpenRange(from: Self, to: Self) -> Bool {
        return self > from && self < to
    }

    /// Returns true if `self` lies in the left-open interval `(from, to]`.
    func liesInLeftOpenRange(from: Self, to: Self) -> Bool {
        return self > from && self <= to
    }

    /// Returns true if `self` lies in the right-open interval `[from, to)`.
    func liesInRightOpenRange(from: Self, to: Self) -> Bool {
        return self >= from && self < to
    }

    /// Returns true if `self` lies in the closed interval `[from, to]`.
    func liesInClosedRange(from: Self, to: Self) -> Bool {
        return self >= from && self <= to
    }
}

internal extension BinaryFloatingPoint {
    @inlinable static var half: Self { 0.5 }
    @inlinable static var quarter: Self { 0.25 }
    @inlinable static var one: Self { 1 }
    @inlinable static var two: Self { 2 }
    @inlinable static var pi: Self { 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651e+00 }
    @inlinable static var twoPi: Self { 6.28318530717958647692528676655900576839433879875021164194988918461563281257241799725606965068423413596429617303e+00 }
    @inlinable static var halfPi: Self { 1.57079632679489661923132169163975144209858469968755291048747229615390820314310449931401741267105853399107404326e+00 }
    @inlinable static var quarterPi: Self { 0.785398163397448309615660845819875721049292349843776455243736148076954101571552249657008706335529266995537021628320576661773 }
    @inlinable static var thirdPi: Self { 1.04719755119659774615421446109316762806572313312503527365831486410260546876206966620934494178070568932738269550e+00 }
    @inlinable static var twoThirdsPi: Self { 2.09439510239319549230842892218633525613144626625007054731662972820521093752413933241868988356141137865476539101e+00 }
    @inlinable static var threeQuartersPi: Self { 2.35619449019234492884698253745962716314787704953132936573120844423086230471465674897102611900658780098661106488e+00 }
    @inlinable static var sixthPi: Self { 5.23598775598298873077107230546583814032861566562517636829157432051302734381034833104672470890352844663691347752e-01 }
    @inlinable static var piSqr: Self { 9.86960440108935861883449099987615113531369940724079062641334937622004482241920524300177340371855223182402591377e+00 }
    @inlinable static var rootTwo: Self { 1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753432764157273501384623e+00 }
    @inlinable static var rootThree: Self { 1.73205080756887729352744634150587236694280525381038062805580697945193301690880003708114618675724857567562614142e+00 }
    @inlinable static var rootPi: Self { 1.77245385090551602729816748334114518279754945612238712821380778985291128459103218137495065673854466541622682362e+00 }
    @inlinable static var oneDivPi: Self { 0.31830988618379067153776752674502872406891929148091289749533468811779359526845307018022760553250617191214568545351 }
    @inlinable static var oneDivTwoPi: Self { 1.59154943091895335768883763372514362034459645740456448747667344058896797634226535090113802766253085956072842727e-01 }
    @inlinable static var neDivRootPi: Self { 5.64189583547756286948079451560772585844050629328998856844085721710642468441493414486743660202107363443028347906e-01 }
    @inlinable static var twoDivPi: Self { 6.36619772367581343075535053490057448137838582961825794990669376235587190536906140360455211065012343824291370907e-01 }
    @inlinable static var twoDivRootPi: Self { 1.12837916709551257389615890312154517168810125865799771368817144342128493688298682897348732040421472688605669581272 }
    @inlinable static var lnTwo: Self { 6.93147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687542001481021e-01 }
    @inlinable static var lnTen: Self { 2.30258509299404568401799145468436420760110148862877297603332790096757260967735248023599720508959829834196778404e+00 }
    @inlinable static var lnLnTwo: Self { -3.66512920581664327012439158232669469454263447837105263053677713670561615319352738549455822856698908358302523045e-01 }
    @inlinable static var e: Self { 2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193e+00 }
    @inlinable static var oneDivE: Self { 3.67879441171442321595523770161460867445811131031767834507836801697461495744899803357147274345919643746627325277e-01 }
    @inlinable static var euler: Self { 5.77215664901532860606512090082402431042159335939923598805767234884867726777664670936947063291746749514631447250e-01 }
    @inlinable static var catalan: Self { 9.15965594177219015054603514932384110774149374281672134266498119621763019776254769479356512926115106248574422619e-01 }
    @inlinable static var zetaThree: Self { 1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933525814619915780e+00 }
    @inlinable static var phi: Self { 1.61803398874989484820458683436563811772030917980576286213544862270526046281890244970720720418939113748475408808e+00 }
    @inlinable static var holtsmarkEntropy: Self { 2.06944850513462440031558003845421663807675255168165702483694991535648437010785015575051452878363155526878657629 }
    @inlinable static var sqrt2PiInv: Self { 0.39894228040143267793994605993438186847585863116493465766592582967065792589930183850125233390730693643030255886 }
}
