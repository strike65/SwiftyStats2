/// Prelude that reexports core dependencies used across SwiftyStats modules.
///
/// This target centralizes common imports (Foundation, Numerics, SwiftyBoost, and C libs)
/// so downstream modules can rely on a single prelude import. Nothing in this file alters
/// behavior; it only exposes dependencies for convenience.

@_exported import SwiftyBoost
@_exported import Numerics
//@_exported import RealModule
@_exported import Foundation
#if canImport(CoreGraphics)
import CoreGraphics
#endif
#if canImport(Darwin)
import Darwin
#elseif canImport(Glibc)
import Glibc
#endif

