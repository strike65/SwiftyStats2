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
/// Error modeling and entropy error definitions shared across SwiftyStats.

import SwiftyStatsPrelude

/// A rich, contextual error used by SwiftyStats.
///
/// `SSError` is an `NSError` subclass that captures additional source metadata
/// (file, line, and function) alongside a categorised `ErrorType`.
/// It conforms to `LocalizedError` for user-facing descriptions and is
/// marked `@unchecked Sendable` as it is a reference type with mutable
/// properties. Only pass instances across concurrency domains if you
/// control their mutation.
///
/// - Note: The `domain` is fixed to `de.volker-thieme.SwiftyStats`.
/// - SeeAlso: ``SSError/ErrorType``
public class SSError: NSError, LocalizedError, @unchecked Sendable {
   // MARK: - Classification

   /// Categorisation of the underlying fault.
   ///
   /// This value determines both the error code of the underlying `NSError` and
   /// the text returned by `description` and `localizedDescription`.
   public var type: ErrorType = .none

   // MARK: - Source Context

   /// Source line number where the error originated.
   public var line: Int = 0

   /// Function in which the error was detected.
   ///
   /// This is typically populated via `#function` at the call site.
   public var function: String = ""

   /// Fully-qualified file name associated with the failure.
   ///
   /// This is typically populated via `#fileID` or `#filePath` at the call site.
   public var file: String = ""

   /// Optional, human-readable explanation describing the failure.
   public private(set) var reason: String?

   /// Error severity and domain codes emitted by SwiftyStats.
   ///
   /// The raw value of each case maps to the `code` of the produced `NSError`.
   public enum ErrorType: Int {
       /// No error.
       case none
       /// Invalid arguments were provided to a function or initializer.
       case invalidArgument
       /// A collection index was out of range.
       case indexOutOfRange
       /// A row index was out of range.
       case rowIndexOutOfRange
       /// A column index was out of range.
       case columnIndexOutOfRange
       /// A referenced row name was not found.
       case unknownRowName
       /// A referenced column name was not found.
       case unknownColumnName
       /// The function is not defined for the provided domain or parameters.
       case functionNotDefinedInDomainProvided
       /// Required data was missing.
       case missingData
       /// The data was in an unexpected or unsupported format.
       case wrongDataFormat
       /// Dimensions or sizes did not match the required constraints.
       case sizeMismatch
       /// The maximum number of iterations was reached before convergence.
       case maxNumberOfIterationReached
       /// The operation is only available for numeric elements.
       case availableOnlyForNumbers
       /// An underlying POSIX error occurred.
       case posixError
       /// The file could not be written.
       case fileNotWriteable
       /// The specified directory does not exist.
       case directoryDoesNotExist
       /// The specified file was not found.
       case fileNotFound
       /// The specified file already exists.
       case fileExists
       /// A required object could not be created.
       case errorCreatingObject
       /// An unexpected internal error occurred.
       case internalError
       /// A mathematical singularity was encountered.
       case singularity
       /// The exponent exceeded the maximum supported value (e.g., for hypergeometric pFq).
       case maxExponentExceeded
   }
   
   // MARK: - NSError Overrides

   /// Developer-focused summary including file, line, and function metadata.
   ///
   /// This is intended for diagnostics and logging, and may reveal file paths.
   override public var description: String {
       switch self.type {
       case .none:
           return "No error"
       case .invalidArgument:
           return "Invalid argument in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .indexOutOfRange:
           return "Index out of range in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .rowIndexOutOfRange:
           return "Row index out of range in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .columnIndexOutOfRange:
           return "Column index out of range in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .unknownRowName:
           return "Unknown row name in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .unknownColumnName:
           return "Unknown column name in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .functionNotDefinedInDomainProvided:
           return "Called function is not defined in that domain in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .missingData:
           return "Missing data in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .wrongDataFormat:
           return "Wrong data format in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .sizeMismatch:
           return "Size mismatch in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .maxNumberOfIterationReached:
           return "Maximum number of iterations reached in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .availableOnlyForNumbers:
           return "Only available for numbers in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .posixError:
           return "POSIX error in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .fileNotWriteable:
           return "File not writable error in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .fileExists:
           return "File not writable because the file already exists in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .directoryDoesNotExist:
           return "Directory does not exist in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .fileNotFound:
           return "File does not exist in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .errorCreatingObject:
           return "Unable to create examine object: " + self.file + " Line: \(self.line) in function: " + self.function
       case .internalError:
           return "Fatal internal error: " + self.file + " Line: \(self.line) in function: " + self.function + ". Contact the developer."
       case .singularity:
           return "Argument singularity: " + self.file + " Line: \(self.line) in function: " + self.function + ". Function returns INF."
       case .maxExponentExceeded:
           return "Value of exponent required for summation (pFq):" + self.file + " Line: \(self.line) in function: " + self.function + ". Hints: (1) try using lnpfq = 1 or (2) use Float80."
           
       }
   }

   /// Localized, user-facing description of the error.
   ///
   /// This mirrors the developer `description` but is suitable for presentation
   /// to users. It still includes minimal context to aid support.
   override public var localizedDescription: String {
       if let reason, !reason.isEmpty {
           return reason
       }
       switch self.type {
       case .none:
           return "No error"
       case .invalidArgument:
           return "Invalid argument in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .indexOutOfRange:
           return "Index out of range in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .rowIndexOutOfRange:
           return "Row index out of range in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .columnIndexOutOfRange:
           return "Column index out of range in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .unknownRowName:
           return "Unknown row name in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .unknownColumnName:
           return "Unknown column name in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .functionNotDefinedInDomainProvided:
           return "Called function is not defined in that domain in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .missingData:
           return "Missing data in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .wrongDataFormat:
           return "Wrong data format in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .sizeMismatch:
           return "Size mismatch in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .maxNumberOfIterationReached:
           return "Maximum number of iterations reached in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .availableOnlyForNumbers:
           return "Only available for numbers in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .posixError:
           return "POSIX error in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .fileNotWriteable:
           return "File not writable error in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .fileExists:
           return "File not writable because the file already exists in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .directoryDoesNotExist:
           return "Directory does not exist in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .fileNotFound:
           return "File does not exist in: " + self.file + " Line: \(self.line) in function: " + self.function
       case .errorCreatingObject:
           return "Unable to create examine object: " + self.file + " Line: \(self.line) in function: " + self.function
       case .internalError:
           return "Fatal internal error: " + self.file + " Line: \(self.line) in function: " + self.function + ". Contact the developer."
       case .singularity:
           return "Argument singularity: " + self.file + " Line: \(self.line) in function: " + self.function + ". Function returns INF."
       case .maxExponentExceeded:
           return "Value of exponent required for summation (pFq):" + self.file + " Line: \(self.line) in function: " + self.function + ". Hints: (1) try using lnpfq = 1 or (2) use Float80."

       }
   }
   
   /// Additional, user-facing reason string surfaced through `NSError`.
   override public var localizedFailureReason: String? {
       return reason
   }
   
   // MARK: - Initialization

   /// Creates a new `SSError` with concrete source metadata.
   ///
   /// - Parameters:
   ///   - type: The categorised error type. Also used as the `NSError` code.
   ///   - file: The file path at which the error originated. Prefer `#fileID`/`#filePath`.
   ///   - line: The line number at which the error originated. Prefer `#line`.
   ///   - function: The function name in which the error was detected. Prefer `#function`.
   ///   - reason: The reason for failing.
   public init(type: ErrorType, file: String, line: Int, function: String, reason: String? = nil) {
       var info: [String: Any] = [:]
       if let reason, !reason.isEmpty {
           info[NSLocalizedFailureReasonErrorKey] = reason
       }
       super.init(domain: "de.volker-thieme.SwiftyStats", code: type.rawValue, userInfo: info)
       self.type = type
       self.line = line
       self.function = function
       self.file = file
       self.reason = reason
   }
   
   /// Decoder-based initializer to support `Codable` interoperability via `NSCoding`.
   ///
   /// - Parameter aDecoder: The decoder used to reconstruct the object.
   required public init?(coder aDecoder: NSCoder) {
       super.init(coder: aDecoder)
       self.reason = aDecoder.decodeObject(forKey: "SSError.reason") as? String
   }

   /// Persists the custom `SSError` fields for `NSCoding` serialization.
   ///
   /// - Parameter aCoder: Encoder receiving the stored error metadata.
   override public func encode(with aCoder: NSCoder) {
       super.encode(with: aCoder)
       if let reason {
           aCoder.encode(reason, forKey: "SSError.reason")
       }
   }
   
}

/// Errors thrown by entropy computations.
enum EntropyError<SSElement: Hashable & Comparable & Codable & Sendable>: Error, LocalizedError {
    case emptyDataset
    case conversionFailed(element: SSElement)
    case invalidParameter(String)
    case insufficientData
    case numericalError(String)
    case emptySubset
    
    var errorDescription: String? {
        switch self {
        case .emptyDataset:
            return "Dataset is empty"
        case .conversionFailed(let element):
            return "Conversion failed for element: \(element)"
        case .invalidParameter(let msg):
            return "Invalid parameter: \(msg)"
        case .insufficientData:
            return "Insufficient amount of data"
        case .numericalError(let msg):
            return "Numerical error: \(msg)"
        case .emptySubset:
            return "Empty subset"
        }
    }
}
