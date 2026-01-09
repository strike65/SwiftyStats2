//
//  Created by VT on 28.10.25.
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

/// Lightweight logging helpers for SwiftyStats and SSExamine operations.

import SwiftyStatsPrelude
#if canImport(OSLog)
import OSLog
#endif

/// Encapsulates logging endpoints for different SwiftyStats subsystems.
///
/// This type centralizes the definition of loggers used throughout the library.
/// When `OSLog` is available (Apple platforms, appropriate SDKs), it exposes
/// strongly-typed `Logger` instances. Otherwise, it falls back to a simple
/// stdout-based logger to keep logging available in all build environments.
struct SSLogger {
    #if canImport(OSLog)
    /// Statistics- and parameter-related messages (errors, warnings).
    ///
    /// Category: "functions_parameters"
    /// Subsystem: Defaults to the app's bundle identifier, or
    /// "de.volker-thieme.SwiftyStats" if unavailable.
    @available(macOS 11.0, iOS 14.0, tvOS 14.0, watchOS 7.0, *)
    static let stat = Logger(
        subsystem: Bundle.main.bundleIdentifier ?? "de.volker-thieme.SwiftyStats",
        category: "functions_parameters"
    )

    /// Internal developer-facing diagnostics for severe or unexpected conditions.
    ///
    /// Category: "severe_bugs"
    /// Subsystem: Defaults to the app's bundle identifier, or
    /// "de.volker-thieme.SwiftyStats" if unavailable.
    @available(macOS 11.0, iOS 14.0, tvOS 14.0, watchOS 7.0, *)
    static let dev = Logger(
        subsystem: Bundle.main.bundleIdentifier ?? "de.volker-thieme.SwiftyStats",
        category: "severe_bugs"
    )

    /// Filesystem-related diagnostics (I/O failures, permissions, etc.).
    ///
    /// Category: "filesystem"
    /// Subsystem: Defaults to the app's bundle identifier, or
    /// "de.volker-thieme.SwiftyStats" if unavailable.
    @available(macOS 11.0, iOS 14.0, tvOS 14.0, watchOS 7.0, *)
    static let fs = Logger(
        subsystem: Bundle.main.bundleIdentifier ?? "de.volker-thieme.SwiftyStats",
        category: "filesystem"
    )
    #else
    /// Lightweight stdout-based logger used when `OSLog` is unavailable.
    ///
    /// This is intended for non-Apple platforms or environments where `OSLog`
    /// cannot be imported. Messages are formatted with a level and category.
    struct SimpleLogger {
        /// Logical grouping for the messages emitted by this logger.
        let category: String

        /// Emits an error-level message.
        ///
        /// - Parameter message: The message to print.
        func error(_ message: String) {
            print("[ERROR][\(category)] \(message)")
        }

        /// Emits an informational message.
        ///
        /// - Parameter message: The message to print.
        func info(_ message: String) {
            print("[INFO][\(category)] \(message)")
        }

        /// Emits a debug-level message.
        ///
        /// - Parameter message: The message to print.
        func debug(_ message: String) {
            print("[DEBUG][\(category)] \(message)")
        }
    }

    /// Statistics- and parameter-related messages (errors, warnings).
    static let stat = SimpleLogger(category: "functions_parameters")
    /// Internal developer-facing diagnostics for severe or unexpected conditions.
    static let dev = SimpleLogger(category: "severe_bugs")
    /// Filesystem-related diagnostics (I/O failures, permissions, etc.).
    static let fs = SimpleLogger(category: "filesystem")
    #endif
}

/// Convenience façade for writing log messages without exposing logging backends.
///
/// Use these static functions to record errors or diagnostics from call sites
/// without importing or depending on `OSLog` directly. On supported Apple
/// platforms at suitable OS versions, messages are routed to `OSLog`. Otherwise,
/// they are printed to stdout as a fallback.
public enum SSLog {
    /// Logs a statistics-related error, prioritising `OSLog` when available.
    ///
    /// Intended for issues around parameter validation, domain errors, or
    /// statistical computation failures that are relevant to end users.
    ///
    /// - Parameter message: The message to log.
    public static func statisticsError(_ message: String) {
#if canImport(OSLog)
#if os(macOS) || os(iOS) || os(tvOS) || os(watchOS) || os(visionOS)
        if #available(macOS 11.0, iOS 14.0, tvOS 14.0, watchOS 7.0, visionOS 1.0, *) {
            SSLogger.stat.error("\(message, privacy: .public)")
        } else {
            print("[ERROR][functions_parameters] \(message)")
        }
#else
        // Non-Apple platforms that still import OSLog (rare). Use a simple print fallback.
        print("[ERROR][functions_parameters] \(message)")
#endif
#else
        print("[ERROR][functions_parameters] \(message)")
#endif
    }

    /// Logs an internal developer-facing error to aid debugging.
    ///
    /// Use this for unexpected states, assertions, or severe conditions that
    /// indicate a bug rather than a user-facing error.
    ///
    /// - Parameter message: The message to log.
    public static func developerError(_ message: String) {
#if canImport(OSLog)
#if os(macOS) || os(iOS) || os(tvOS) || os(watchOS) || os(visionOS)
        if #available(macOS 11.0, iOS 14.0, tvOS 14.0, watchOS 7.0, visionOS 1.0, *) {
            SSLogger.dev.error("\(message, privacy: .public)")
        } else {
            print("[ERROR][severe_bugs] \(message)")
        }
#else
        // Non-Apple platforms that still import OSLog (rare). Use a simple print fallback.
        print("[ERROR][severe_bugs] \(message)")
#endif
#else
        print("[ERROR][severe_bugs] \(message)")
#endif
    }

    /// Logs filesystem access errors for storage diagnostics.
    ///
    /// Use this for file read/write failures, missing directories, permission
    /// problems, or other I/O-related issues.
    ///
    /// - Parameter message: The message to log.
    public static func filesystemError(_ message: String) {
#if canImport(OSLog)
#if os(macOS) || os(iOS) || os(tvOS) || os(watchOS) || os(visionOS)
        if #available(macOS 11.0, iOS 14.0, tvOS 14.0, watchOS 7.0, visionOS 1.0, *) {
            SSLogger.fs.error("\(message, privacy: .public)")
        } else {
            print("[ERROR][filesystem] \(message)")
        }
#else
        // Non-Apple platforms that still import OSLog (rare). Use a simple print fallback.
        print("[ERROR][filesystem] \(message)")
#endif
#else
        print("[ERROR][filesystem] \(message)")
#endif
    }
}
