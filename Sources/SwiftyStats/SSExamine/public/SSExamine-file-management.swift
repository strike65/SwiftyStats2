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

/// Persistence, import/export, and file management helpers for SSExamine content.

import SwiftyStatsPrelude
import ZIPFoundation
import SwiftyBoost

/// File import/export helpers that preserve the statistical meaning of an `SSExamine` sample.
///
/// Routines here treat on-disk data as realizations of the same random variable, whether
/// stored plainly, in CSV-like delimited text, or zipped archives. Parsing maintains the
/// original insertion order (raw sequence) while reconstructing frequency counts, so any
/// computed moments or entropy metrics remain invariant to the persistence round-trip.
extension SSExamine {
    
    /// Import a delimited text file and parse each separated token with a client-supplied parser.
    ///
    /// Safety and robustness improvements:
    /// - Uses URL-based APIs and expands `~`.
    /// - Validates separator non-emptiness to avoid pathological splits.
    /// - Optionally strips enclosure strings; if empty, ignores to avoid O(n²) replace loops on empty pattern.
    /// - Provides actionable `SSError` on I/O failures and parser failures.
    /// - Avoids silent nil returns except when the file clearly does not contain the separator (original behavior preserved).
    /// - Detects ZIP archives and transparently reads the first non-directory entry as text using `stringEncoding`.
    ///
    /// - Parameters:
    ///   - path: File system path; `~` will be expanded.
    ///   - separator: Non-empty separator string to split the file content.
    ///   - elementsEnclosedBy: Optional enclosure to strip from each token.
    ///   - stringEncoding: Encoding of the file on disk.
    ///   - parser: A closure that converts a token into an `SSElement`; return `nil` to signal a parse failure.
    /// - Returns: An `SSExamine` instance or `nil` if the file does not contain the separator or parsing fails.
    public class func examine(
        fromFile path: String,
        separator: String,
        elementsEnclosedBy: String? = nil,
        stringEncoding: String.Encoding = .utf8,
        _ parser: (String?) -> SSElement?
    ) throws -> SSExamine<SSElement, FP>? {
        // Validate separator
        guard !separator.isEmpty else {
            let msg = "Separator must not be empty."
            SSLog.filesystemError(msg)
            throw SSError(type: .wrongDataFormat, file: #fileID, line: #line, function: #function, reason: msg)
        }
        
        // Expand ~ and form URL
        let expanded = NSString(string: path).expandingTildeInPath
        // Basic NUL guard
        guard !expanded.contains("\0") else {
            let msg = "Path contains NUL character."
            SSLog.filesystemError(msg)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: msg)
        }
        let url = URL(filePath: expanded)
        
        // Pre-check readability for clearer error, but still rely on IO errors below
        let fm = FileManager.default
        var isDir: ObjCBool = false
        guard fm.fileExists(atPath: url.path, isDirectory: &isDir), !isDir.boolValue else {
            SSLog.filesystemError("File not found or is a directory: \(url.path)")
            throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: "File not found or is a directory: \(url.path)")
        }
        guard fm.isReadableFile(atPath: url.path) else {
            SSLog.filesystemError("File is not readable: \(url.path)")
            throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: "File is not readable: \(url.path)")
        }
        
        // Read entire file or unzip and read entry text.
        let content: String
        do {
            if try isZip(url: url) {
                content = try readFirstZipEntryString(url: url, encoding: stringEncoding)
            } else {
                content = try String(contentsOf: url, encoding: stringEncoding)
            }
        } catch let e as SSError {
            // Already logged and wrapped
            throw e
        } catch {
            SSLog.filesystemError("Failed reading file: \(url.path). Error: \(error.localizedDescription)")
            throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: "Failed reading file: \(error.localizedDescription)")
        }
        
        // If separator not found, preserve original behavior and return nil
        guard content.contains(separator) else {
            SSLog.statisticsError("File \(url.lastPathComponent) does not contain the separator \"\(separator)\"")
            return nil
        }
        guard !content.isEmpty else {
            SSLog.statisticsError("File \(url.lastPathComponent) is empty; cannot build SSExamine")
            return nil
        }
        
        // Split and optionally strip enclosures
        let rawTokens = content.components(separatedBy: separator)
        let enclosure = (elementsEnclosedBy?.isEmpty == false) ? elementsEnclosedBy! : nil
        let tokens: [String] = {
            guard let e = enclosure else { return rawTokens }
            // Avoid pathological performance by skipping replace when enclosure not present
            if rawTokens.allSatisfy({ !$0.contains(e) }) { return rawTokens }
            return rawTokens.map { $0.replacingOccurrences(of: e, with: "") }
        }()
        
        // Parse tokens
        var values: [SSElement] = []
        values.reserveCapacity(tokens.count)
        for t in tokens {
            // Skip empty tokens (compatibility with previous code)
            if t.isEmpty { continue }
            guard let v = parser(t) else {
                // Stop parsing on first failure and return nil (preserve original behavior),
                // but log a helpful message for diagnostics.
                SSLog.statisticsError("Parsing failed for token: \"\(t)\" from file: \(url.lastPathComponent)")
                return nil
            }
            values.append(v)
        }
        
        // Construct result
        let name = url.lastPathComponent
        let result = SSExamine<SSElement, FP>(usingArray: values, name: name, characterSet: nil)
        // Heuristically infer and set level of measurement
        result.levelOfMeasurement = inferLevelOfMeasurement(for: result)
        return result
    }
    
    /// Load an `SSExamine` instance from a JSON file.
    ///
    /// Safety improvements:
    /// - Uses Data-based decoding directly from file URL (no String round-trip).
    /// - Propagates decoding errors verbosely.
    /// - Validates existence, non-directory, and readability.
    public class func examine(
        fromJSON path: String,
        stringEncoding: String.Encoding = .utf8
    ) throws -> SSExamine<SSElement, FP>? {
        let expanded = NSString(string: path).expandingTildeInPath
        guard !expanded.contains("\0") else {
            let msg = "Path contains NUL character."
            SSLog.filesystemError(msg)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: msg)
        }
        let url = URL(filePath: expanded)
        
        let fm = FileManager.default
        var isDir: ObjCBool = false
        guard fm.fileExists(atPath: url.path, isDirectory: &isDir), !isDir.boolValue else {
            SSLog.filesystemError("File \(path) not found or is a directory.")
            throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: "File \(path) not found or is a directory.")
        }
        guard fm.isReadableFile(atPath: url.path) else {
            SSLog.filesystemError("File \(path) is not readable.")
            throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: "File \(path) is not readable.")
        }
        
        // Read raw bytes and decode
        let data: Data
        do {
            data = try Data(contentsOf: url, options: [.mappedIfSafe])
        } catch {
            SSLog.filesystemError("Failed reading JSON file: \(url.path). Error: \(error.localizedDescription)")
            throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: "Failed reading JSON file: \(error.localizedDescription)")
        }
        
        do {
            let decoded = try JSONDecoder().decode(SSExamine<SSElement, FP>.self, from: data)
            // Heuristically infer and set level of measurement
            decoded.levelOfMeasurement = inferLevelOfMeasurement(for: decoded)
            return decoded
        } catch {
            SSLog.filesystemError("JSON decode failed for \(url.lastPathComponent): \(error.localizedDescription)")
            throw error
        }
    }
    
    /// Save the dataset to a delimited text file or a ZIP containing the text.
    ///
    /// Safety improvements:
    /// - Uses a write-to-temp-then-atomic-move strategy to avoid partial writes and TOCTOU issues.
    /// - Validates parent directory (creates if missing); rejects when parent exists but is not a directory.
    /// - Honors `overwrite` by using replace semantics when appropriate.
    /// - Uses URL-based APIs and avoids race-prone pre-check-only logic.
    /// - Logs and throws `SSError` with underlying reasons.
    ///
    /// - Parameters:
    ///   - path: Destination path; `~` is expanded.
    ///   - atomically: If true, write string atomically to the temp file (still followed by atomic move).
    ///   - overwrite: If false and file exists, throw `.fileExists`.
    ///   - separator: Delimiter between elements.
    ///   - encloseElementsBy: Optional enclosure string for each element.
    ///   - asRow: If true, write all items on one line; otherwise write each on its own line.
    ///   - stringEncoding: Encoding for the output file.
    ///   - compressAsZip: When true, write a ZIP archive containing one entry with the text data.
    ///   - zipEntryFileName: Optional file name for the inner entry (defaults to basename + ".txt").
    /// - Returns: True on success.
    public func saveTo(
        fileName path: String,
        atomically: Bool = true,
        overwrite: Bool,
        separator: String = ",",
        encloseElementsBy: String? = nil,
        asRow: Bool = true,
        stringEncoding: String.Encoding = .utf8,
        compressAsZip: Bool = false,
        zipEntryFileName: String? = nil
    ) throws -> Bool {
        // Validate separator (allow empty only when single-column line-per-item)
        if separator.isEmpty && asRow {
            let msg = "Separator must not be empty when writing as a single row."
            SSLog.filesystemError(msg)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: msg)
        }
        
        let fm = FileManager.default
        let expanded = NSString(string: path).expandingTildeInPath
        guard !expanded.contains("\0") else {
            let msg = "Path contains NUL character."
            SSLog.filesystemError(msg)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: msg)
        }
        var destURL = URL(filePath: expanded)
        
        // Append ".zip" to the destination filename if compression is requested and it’s missing.
        if compressAsZip {
            let ext = destURL.pathExtension.lowercased()
            if ext != "zip" {
                // Keep existing extension as part of the base name to avoid losing information.
                let baseName = destURL.lastPathComponent
                let parent = destURL.deletingLastPathComponent()
                let zippedName = baseName + ".zip"
                destURL = parent.appending(path: zippedName, directoryHint: .notDirectory)
            }
        }
        
        let parentURL = destURL.deletingLastPathComponent()
        
        // Ensure parent directory exists and is a directory
        var isDir: ObjCBool = false
        if !fm.fileExists(atPath: parentURL.path, isDirectory: &isDir) {
            do {
                try fm.createDirectory(at: parentURL, withIntermediateDirectories: true, attributes: nil)
            } catch {
                SSLog.filesystemError("Unable to create directory for export: \(parentURL.path). Error: \(error.localizedDescription)")
                throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: "Unable to create directory for export: \(error.localizedDescription)")
            }
        } else if !isDir.boolValue {
            SSLog.filesystemError("Parent path is not a directory: \(parentURL.path)")
            throw SSError(type: .directoryDoesNotExist, file: #fileID, line: #line, function: #function, reason: "Parent path is not a directory: \(parentURL.path)")
        }
        
        // Materialize content
        guard let content = self.itemsAsString(delimiter: separator, asRow: asRow, encloseElementsBy: encloseElementsBy) else {
            SSLog.filesystemError("Failed to serialize items to string.")
            throw SSError(type: .wrongDataFormat, file: #fileID, line: #line, function: #function, reason: "Failed to serialize items to string.")
        }
        
        // Prepare a secure temporary file in the same directory for atomic replace
        let tmpName = UUID().uuidString + (compressAsZip ? ".zip.tmp" : ".tmp")
        let tmpURL = parentURL.appending(path: tmpName, directoryHint: .notDirectory)
        
        do {
            if compressAsZip {
                // Create ZIP at tmpURL with a single entry containing the text data
                try writeZip(at: tmpURL,
                             entryName: makeZipEntryName(for: destURL, override: zipEntryFileName),
                             text: content,
                             encoding: stringEncoding)
            } else {
                // Write plain text to temp file
                try content.write(to: tmpURL, atomically: atomically, encoding: stringEncoding)
            }
        } catch {
            SSLog.filesystemError("Failed writing temporary file: \(tmpURL.lastPathComponent). Error: \(error.localizedDescription)")
            throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: "Failed writing temporary file: \(error.localizedDescription)")
        }
        
        // Move temp into place atomically, handling overwrite semantics
        do {
            if fm.fileExists(atPath: destURL.path) {
                if !overwrite {
                    // Clean up temp file
                    try? fm.removeItem(at: tmpURL)
                    SSLog.filesystemError("File already exists and overwrite is false: \(destURL.path)")
                    throw SSError(type: .fileExists, file: #fileID, line: #line, function: #function, reason: "File already exists: \(destURL.lastPathComponent)")
                } else {
                    _ = try fm.replaceItemAt(destURL, withItemAt: tmpURL, backupItemName: nil, options: [.usingNewMetadataOnly])
                }
            } else {
                try fm.moveItem(at: tmpURL, to: destURL)
            }
            return true
        } catch {
            // Attempt cleanup of temp file if it still exists
            try? fm.removeItem(at: tmpURL)
            SSLog.filesystemError("Failed moving file into place: \(error.localizedDescription)")
            throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: "Failed moving file into place: \(error.localizedDescription)")
        }
    }
    
    /// Exports this dataset as JSON to the given file path using an atomic, replace-in-place strategy.
    ///
    /// Behavior:
    /// - Expands `~` in `path` and rejects NUL characters.
    /// - Ensures the parent directory exists (creates it if missing) and is a directory.
    /// - Encodes `self` with `JSONEncoder()` to UTF-8 `Data`.
    /// - Writes to a temporary file in the destination directory, then atomically moves or replaces it at `path`.
    /// - Honors `overwrite`: if `false` and the destination exists, throws `.fileExists`; if `true`, uses `FileManager.replaceItemAt` to swap atomically.
    /// - Logs failures via `SSLog.filesystemError` and throws `SSError` with actionable reasons.
    ///
    /// Parameters:
    /// - path: Destination path for the JSON file. `~` will be expanded.
    /// - atomically: When `true`, writes the temporary file with `.atomic` in addition to the final atomic move/replace.
    /// - overwrite: When `true`, replaces an existing file atomically; when `false`, throws if the file exists.
    ///
    /// Returns:
    /// - `true` on success.
    ///
    /// Throws:
    /// - `SSError(.invalidArgument)` if the path contains a NUL character.
    /// - `SSError(.directoryDoesNotExist)` if the parent directory cannot be created or is not a directory.
    /// - `SSError(.wrongDataFormat)` if JSON encoding fails.
    /// - `SSError(.fileNotWriteable)` if writing the temporary file or the final move/replace fails.
    /// - `SSError(.fileExists)` if `overwrite == false` and the destination already exists.
    public func exportJSON(fileName path: String, atomically: Bool = true, overwrite: Bool) throws -> Bool {
        let fm = FileManager.default
        let expanded = NSString(string: path).expandingTildeInPath
        // Guard against embedded NUL which Foundation paths do not support
        guard !expanded.contains("\0") else {
            let msg = "Path contains NUL character."
            SSLog.filesystemError(msg)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: msg)
        }
        let destURL = URL(filePath: expanded, directoryHint: .notDirectory)
        let parentURL = destURL.deletingLastPathComponent()
        
        // Ensure parent exists and is a directory
        var isDir = ObjCBool(false)
        if !fm.fileExists(atPath: parentURL.path, isDirectory: &isDir) {
            do {
                try fm.createDirectory(at: parentURL, withIntermediateDirectories: true, attributes: nil)
            } catch {
                SSLog.filesystemError("Unable to create directory for export: \(parentURL.path). Error: \(error.localizedDescription)")
                throw SSError(type: .directoryDoesNotExist, file: #fileID, line: #line, function: #function, reason: "Unable to create directory for export: \(error.localizedDescription)")
            }
        } else if !isDir.boolValue {
            SSLog.filesystemError("Parent path is not a directory: \(parentURL.path)")
            throw SSError(type: .directoryDoesNotExist, file: #fileID, line: #line, function: #function, reason: "Parent path is not a directory: \(parentURL.path)")
        }
        
        // Encode JSON to Data using specified stringEncoding where possible
        // JSON data is UTF-8 by definition; we keep Data and write it, avoiding lossy String round-trips.
        let encoder = JSONEncoder()
        let jsonData: Data
        do {
            jsonData = try encoder.encode(self)
        } catch {
            SSLog.filesystemError("JSON encoding failed: \(error.localizedDescription)")
            throw SSError(type: .wrongDataFormat, file: #fileID, line: #line, function: #function, reason: "JSON encoding failed: \(error.localizedDescription)")
        }
        
        // Prepare a temporary file in the same directory for atomic replace
        let tmpURL = parentURL.appending(path: UUID().uuidString + ".json.tmp", directoryHint: .notDirectory)
        do {
            // Write bytes. atomically parameter is honored by writing to a temp path first anyway.
            // If 'atomically' is requested, we can still just write to tmpURL; the final replace is atomic.
            if atomically {
                try jsonData.write(to: tmpURL, options: [.atomic])
            } else {
                try jsonData.write(to: tmpURL, options: [])
            }
        } catch {
            SSLog.filesystemError("Unable to write temporary JSON file: \(error.localizedDescription)")
            // Best-effort cleanup
            try? fm.removeItem(at: tmpURL)
            throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: "Unable to write temporary JSON file: \(error.localizedDescription)")
        }
        
        // Move temp into place atomically; honor overwrite semantics
        do {
            if fm.fileExists(atPath: destURL.path) {
                if !overwrite {
                    try? fm.removeItem(at: tmpURL)
                    SSLog.filesystemError("File already exists and overwrite is false: \(destURL.lastPathComponent)")
                    throw SSError(type: .fileExists, file: #fileID, line: #line, function: #function, reason: "File already exists")
                } else {
                    _ = try fm.replaceItemAt(destURL, withItemAt: tmpURL, backupItemName: nil, options: [.usingNewMetadataOnly])
                }
            } else {
                try fm.moveItem(at: tmpURL, to: destURL)
            }
            return true
        } catch {
            // Cleanup temp on failure
            try? fm.removeItem(at: tmpURL)
            SSLog.filesystemError("Unable to finalize JSON write: \(error.localizedDescription)")
            throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: "Unable to finalize JSON write: \(error.localizedDescription)")
        }
    }

    
    // MARK: - ZIP helpers
    
    /// Builds a reasonable default inner entry name for the zip file.
    /// Defaults to basename(without extension) + ".txt" unless `override` is provided.
    private func makeZipEntryName(for destURL: URL, override: String?) -> String {
        if let o = override, !o.isEmpty { return o }
        let base = destURL.deletingPathExtension().lastPathComponent
        return base.isEmpty ? "data.txt" : base + ".txt"
    }
    
    /// Writes a ZIP archive at `zipURL` with a single entry named `entryName` containing the given text.
    private func writeZip(at zipURL: URL, entryName: String, text: String, encoding: String.Encoding) throws {
        // Convert text to data
        guard let data = text.data(using: encoding) else {
            let msg = "Failed to encode text using \(encoding)."
            SSLog.filesystemError(msg)
            throw SSError(type: .wrongDataFormat, file: #fileID, line: #line, function: #function, reason: msg)
        }
        
        // Use ZIPFoundation’s throwing initializer to create the archive
        let archive: Archive
        do {
            archive = try Archive(url: zipURL, accessMode: .create)
        } catch {
            let msg = "Unable to create ZIP archive: \(error.localizedDescription)"
            SSLog.filesystemError(msg)
            throw SSError(type: .errorCreatingObject, file: #fileID, line: #line, function: #function, reason: msg)
        }
        
        do {
            try archive.addEntry(with: entryName,
                                 type: .file,
                                 uncompressedSize: Int64(data.count),
                                 compressionMethod: .deflate,
                                 bufferSize: 32 * 1024,
                                 progress: nil) { (position: Int64, size: Int) -> Data in
                // Ensure bounds: position within [0, data.count]
                if position < 0 || position > Int64(data.count) {
                    return Data()
                }
                let start = Int(position)
                let end = min(start + size, data.count)
                if start >= end { return Data() }
                return data.subdata(in: start..<end)
            }
        } catch {
            SSLog.filesystemError("Failed to add entry to ZIP: \(error.localizedDescription)")
            throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: "Failed to add entry to ZIP: \(error.localizedDescription)")
        }
    }
    
    
    
    /// Detects whether the file at `url` appears to be a ZIP archive.
    /// Uses both extension-based heuristics and magic bytes "PK\x03\x04".
    private class func isZip(url: URL) throws -> Bool {
        // Quick extension check
        if url.pathExtension.lowercased() == "zip" {
            return true
        }
        // Magic bytes check
        do {
            let fh = try FileHandle(forReadingFrom: url)
            defer { try? fh.close() }
            let sig = try fh.read(upToCount: 4) ?? Data()
            if sig.count == 4 {
                // 50 4B 03 04
                let magic: [UInt8] = [0x50, 0x4B, 0x03, 0x04]
                return Array(sig) == magic
            }
        } catch {
            // If we cannot read header, fall back to extension only
            SSLog.filesystemError("Unable to read file header for ZIP detection: \(error.localizedDescription)")
        }
        return false
    }
    
    /// Reads the first non-directory entry from a ZIP file as a String using `encoding`.
    /// Prefers an entry with a text-like extension when available; otherwise uses the first file entry.
    private class func readFirstZipEntryString(url: URL, encoding: String.Encoding) throws -> String {
        let archive: Archive
        do {
            archive = try Archive(url: url, accessMode: .read)
        } catch {
            let msg = "Unable to open ZIP archive: \(error.localizedDescription)"
            SSLog.filesystemError(msg)
            throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: msg)
        }
        
        // Collect candidate file entries
        let fileEntries = archive.compactMap { entry -> Entry? in
            switch entry.type {
            case .file: return entry
            default: return nil
            }
        }
        guard !fileEntries.isEmpty else {
            let msg = "ZIP archive contains no file entries."
            SSLog.filesystemError(msg)
            throw SSError(type: .wrongDataFormat, file: #fileID, line: #line, function: #function, reason: msg)
        }
        
        // Prefer text-like extensions
        let textExts: Set<String> = ["txt", "csv", "tsv", "json", "dat", "log"]
        let preferred = fileEntries.first { e in
            let ext = URL(filePath: e.path).pathExtension.lowercased()
            return textExts.contains(ext)
        } ?? fileEntries.first!
        
        // Extract into memory
        var data = Data()
        data.reserveCapacity(Int(preferred.uncompressedSize))
        do {
            _ = try archive.extract(preferred, bufferSize: 32 * 1024, progress: nil) { chunk in
                data.append(chunk)
            }
        } catch {
            let msg = "Failed to extract entry '\(preferred.path)': \(error.localizedDescription)"
            SSLog.filesystemError(msg)
            throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: msg)
        }
        
        // Decode to String
        if let s = String(data: data, encoding: encoding) {
            return s
        } else if let sUTF8 = String(data: data, encoding: .utf8) {
            // Fallback to UTF-8 if requested encoding fails
            return sUTF8
        } else {
            let msg = "Failed to decode ZIP entry text using \(encoding)."
            SSLog.filesystemError(msg)
            throw SSError(type: .wrongDataFormat, file: #fileID, line: #line, function: #function, reason: msg)
        }
    }
    
    // MARK: - Level-of-measurement heuristic
    
    /// Infers a conservative level-of-measurement classification from the dataset.
    ///
    /// Rules:
    /// - Empty dataset: .nominal
    /// - Numeric convertible (itemsAsNumericArray != nil):
    ///   • any negative values -> .interval
    ///   • else (all >= 0) -> .ratio
    /// - Non-numeric: .nominal
    private class func inferLevelOfMeasurement(for ex: SSExamine<SSElement, FP>) -> LevelOfMeasurement {
        guard !ex.isEmpty else { return .nominal }
        if let nums = ex.itemsAsNumericArray, !nums.isEmpty {
            let hasNegative = nums.contains(where: { $0 < 0 })
            return hasNegative ? .interval : .ratio
        } else {
            return .nominal
        }
    }
    
    
    /// Loads an array of `RealLike` values from a text file using a configurable separator.
    ///
    /// - Parameters:
    ///   - path: File system path to the text file.
    ///   - separator: Character used to separate tokens (e.g., ',', ';', ' ', '\t' ...).
    ///   - encoding: Text encoding (defaults to UTF-8).
    ///   - allowNaNAndInfinity: If true, tokens "nan", "inf", "+inf", "-inf" (case-insensitive) are allowed and mapped to corresponding Double values.
    /// - Returns: Array of parsed doubles (empty tokens are skipped).
    /// - Throws: An error if the file cannot be read or a token cannot be parsed as a RealLike.
    public class func loadReals(from path: String,
                     separator: Character,
                     encoding: String.Encoding = .utf8,
                     allowNaNAndInfinity: Bool = true) -> SSExamine<FP, FP>? {
        let text: String
        do {
            text = try String(contentsOfFile: path, encoding: encoding)
        }
        catch {
            SSLog.filesystemError("File at \(path) could not be read.")
            return nil
        }
        if text.isEmpty { return SSExamine<FP, FP>() }
        
        // Split strictly by the provided separator; then trim and filter empties
        let tokens = text.split(separator: separator, omittingEmptySubsequences: true)
        
        var result: [FP] = []
        result.reserveCapacity(tokens.count)
        
        for (_, tok) in tokens.enumerated() {
            let trimmed = tok.trimmingCharacters(in: .whitespacesAndNewlines)
            if trimmed.isEmpty { continue }
            
            if allowNaNAndInfinity {
                let lower = trimmed.lowercased()
                if lower == "nan" {
                    result.append(.nan)
                    continue
                } else if lower == "inf" || lower == "+inf" || lower == "infinity" || lower == "+infinity" {
                    result.append(.infinity)
                    continue
                } else if lower == "-inf" || lower == "-infinity" {
                    result.append(-.infinity)
                    continue
                }
            }
            do {
                let value: FP = try RealConverter.to(trimmed)
                result.append(value)
            }
            catch {
                return nil
            }
        }
        let ex: SSExamine<FP, FP> = .init(usingArray: result, name: nil, characterSet: nil)
        return ex
    }
}

