//
//  Created by VT on 23.11.25.
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

/// DataFrame-style view over SSExamine rows for tabular operations.

import SwiftyStatsPrelude
import ZIPFoundation

/// A column-oriented container for statistical examinations.
///
/// `DataFrame` holds an ordered collection of `SSExamine` columns along with
/// associated column names and tags. It provides safe accessors, unique name
/// generation on insertion, and Codable/Copying conformances.
///
/// - Generics:
///   - SSElement: The element type contained by each `SSExamine` column. Must be `Hashable`, `Comparable`, `Codable`, and `Sendable`.
///   - FP: The floating-point-like numeric type used by `SSExamine`. Must conform to `RealLike`.
public class DataFrame<SSElement, FP>: NSObject, NSCopying, Codable, NSMutableCopying, SSDataFrameContainer
where SSElement: Hashable & Comparable & Codable & Sendable, FP: RealLike {

    /// Removes the column with the given name.
    ///
    /// - Parameter name: The column name to remove.
    /// - Returns: The removed `SSExamine` if found; otherwise `nil`.
    /// - Precondition: The data frame must contain at least one column.
    func remove(name: String!) -> SSExamine<SSElement, FP>? {
        precondition(self.nCols > 0)
        if let idx = self.columnNames.firstIndex(of: name) {
            self._columnNames.remove(at: idx)
            self._tags.remove(at: idx)
            self._nCols -= 1
            return self._data.remove(at: idx)
        }
        else {
            return nil
        }
    }
    
    /// Convenience alias for the `SSExamine` column type.
    typealias Examine = SSExamine<SSElement, FP>
    
    
//    public typealias Examine = SSExamine<SSElement, FP>

    /// Backing storage for column names.
    private var _columnNames: Array<String>
    /// Backing storage for column tags.
    private var _tags: Array<String>
    /// Backing storage for `SSExamine` columns.
    private var _data: Array<SSExamine<SSElement, FP>>
    /// Number of columns.
    private var _nCols: Int
    
    /// The names of all columns in order.
    public var columnNames: Array<String> {
        get {
            return self._columnNames
        }
        set {
            self._columnNames = newValue
        }
    }
    
    /// The tags associated with columns in order.
    public var tags: Array<String> {
        get {
            return self._tags
        }
    }
    
    /// The number of columns in the data frame.
    public var nCols: Int {
        get {
            return self._nCols
        }
    }
    
    /// The array of `SSExamine` columns.
    public var data: Array<SSExamine<SSElement, FP>> {
        return self._data
    }
    
    /// A Boolean value indicating whether the data frame has no columns.
    public var isEmpty: Bool {
        return self._nCols == 0
    }
    
    /// Checks whether a column index is within valid bounds.
    ///
    /// - Parameter idx: The zero-based column index to validate.
    /// - Returns: `true` if valid, otherwise `false`.
    private func isValidColumnIndex(_ idx: Int) -> Bool {
        return idx >= 0 && idx < self._nCols
    }
    
    /// Returns the sample size for a given column index.
    ///
    /// - Parameter column: The zero-based column index.
    /// - Returns: The sample size for the column, or `0` if the index is invalid.
    public func sampleSize(column: Int) -> Int {
        if self.isValidColumnIndex(column) {
            let e: SSExamine<SSElement, FP> = self.data[column]
            return e.sampleSize
        }
        else {
            return 0
        }
    }
    
    /// The sum of sample sizes across all columns.
    public var totalSampleSize: Int {
        get {
            var ss: Int = 0
            for e in self._data {
                ss += e.sampleSize
            }
            return ss
        }
    }

    /// Accesses the column at the specified index.
    ///
    /// - Parameter idx: The zero-based column index.
    /// - Returns: The `SSExamine` at the given index.
    /// - Precondition: `idx` must be within bounds.
    public subscript (_ idx: Int) -> SSExamine<SSElement, FP> {
        precondition(self.isValidColumnIndex(idx), "Index out of bounds")
        return self._data[idx]
    }
    
    /// Accesses the column with the specified column name.
    ///
    /// - Parameter col: The column name.
    /// - Returns: The matching `SSExamine` if found; otherwise `nil` and logs an error.
    public subscript (columnName col: String) -> SSExamine<SSElement, FP>? {
        if let idx: Int = self.columnNames.firstIndex(of: col) {
            return self._data[idx]
        }
        else {
            SSLog.statisticsError("No object found for column name: \(col)")
            return nil
        }
    }
    
    
    
    /// Creates an empty data frame.
    override public init() {
        _columnNames = []
        _tags = []
        _data = []
        _nCols = 0
        super.init( )
    }
    
    /// Generates a unique column name from a base and proposed index.
    ///
    /// Behavior:
    /// - If `base` is unused, returns `base`.
    /// - If used, tries `base + "%02d"` (with `i`), then `base + "_%02d"` with increasing suffix until unique.
    ///
    /// - Parameters:
    ///   - base: The preferred base name.
    ///   - proposedIndex: The 1-based index to embed into the name if needed.
    /// - Returns: A unique column name not present in `columnNames`.
    private func makeUniqueColumnName(base: String, proposedIndex i: Int) -> String {
        if self._columnNames.firstIndex(of: base) == nil {
            return base
        }
        var temp = base + String(format: "%02d", arguments: [i as CVarArg])
        var k = 1
        while self._columnNames.firstIndex(of: temp) != nil {
            k += 1
            temp = base + "_" + String(format: "%02d", arguments: [k as CVarArg])
        }
        return temp
    }
    
    /// Creates a data frame from an array of `SSExamine` columns.
    ///
    /// Column names are derived from `SSExamine.name` when available; otherwise
    /// they are auto-generated as `Sample_###`. All names are made unique.
    /// Tags are taken from `SSExamine.tag` or default to an empty string.
    ///
    /// - Parameter data: The columns to include in the data frame.
    /// - Throws: Rethrows any error originating from underlying operations.
    init(data: Array<SSExamine<SSElement, FP>>) throws {
        self._data = []
        self._data.reserveCapacity(data.count)
        self._columnNames = []
        self._columnNames.reserveCapacity(data.count)   // reserve column name capacity
        self._tags = []
        self._tags.reserveCapacity(data.count)          // reserve tag capacity
        self._nCols = data.count
        super.init()
        var i: Int = 1
        for e in data {
            self._data.append(e)
            if let name = e.name {
                let unique = makeUniqueColumnName(base: name, proposedIndex: i)
                self._columnNames.append(unique)
            }
            else {
                // Auto name "Sample_###" and ensure uniqueness as well (rare but safe)
                let auto = "Sample_" + String(format:"%03d", arguments: [i as CVarArg])
                let unique = makeUniqueColumnName(base: auto, proposedIndex: i)
                self._columnNames.append(unique)
            }
            if let t = e.tag {
                self._tags.append(t)
            }
            else {
                self._tags.append("")
            }
            i += 1
        }
    }
    
    /// Encodes this instance into the given encoder.
    ///
    /// - Parameter encoder: The encoder to write data to.
    /// - Throws: Any error encountered during encoding.
    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        try container.encodeIfPresent(self._data, forKey: .data)
        try container.encodeIfPresent(self._tags, forKey: .tags)
        try container.encodeIfPresent(self._columnNames, forKey: .columnNames)
        try container.encodeIfPresent(self._nCols, forKey: .nCols)
    }
    
    /// Creates a new instance by decoding from the given decoder.
    ///
    /// - Parameter decoder: The decoder to read data from.
    /// - Throws: Any error encountered during decoding.
    public required init(from decoder: Decoder) throws {
        // Initialize all stored properties with safe defaults first
        self._data = []
        self._tags = []
        self._columnNames = []
               self._nCols = 0

        let container = try decoder.container(keyedBy: CodingKeys.self)

        // Decode with defaults to ensure properties are always initialized
        if let d = try container.decodeIfPresent(Array<SSExamine<SSElement, FP>>.self, forKey: .data) {
            self._data = d
        }
        if let t = try container.decodeIfPresent(Array<String>.self, forKey: .tags) {
            self._tags = t
        }
        if let n = try container.decodeIfPresent(Array<String>.self, forKey: .columnNames) {
            self._columnNames = n
        }
        if let c = try container.decodeIfPresent(Int.self, forKey: .nCols) {
            self._nCols = c
        } else {
            // Fallback: make nCols consistent with decoded data if present
            self._nCols = self._data.count
        }

        super.init()
    }
    
    /// Coding keys for `Codable` conformance.
    private enum CodingKeys: String, CodingKey {
        case data = "ExamineArray"
        case tags = "tags"
        case columnNames = "columnNames"
        case nRowsMax = "nRowsMax"
        case nCols = "nCols"
    }
    
    /// Creates an immutable copy of the data frame.
    ///
    /// - Parameter zone: The zone identifies an area of memory from which to allocate for the new instance.
    /// - Returns: A new `DataFrame` with duplicated columns, names, and tags.
    public func copy(with zone: NSZone? = nil) -> Any {
        do {
            let res = try DataFrame<SSElement, FP>.init(data: self.data)
            res._tags = self.tags
            res._columnNames = self.columnNames
            return res;
        }
        catch {
            SSLog.developerError("copy(with:) failed")
            return DataFrame()
        }
    }
    
    /// Creates a mutable copy of the data frame.
    ///
    /// - Parameter zone: The zone identifies an area of memory from which to allocate for the new instance.
    /// - Returns: A new `DataFrame` identical to `copy(with:)`.
    public func mutableCopy(with zone: NSZone? = nil) -> Any {
        return self.copy(with: zone)
    }
    
    /// Creates a mutable copy of the data frame.
    ///
    /// - Returns: A new `DataFrame` identical to `copy(with:)`.
    public override func mutableCopy() -> Any {
        return self.copy(with: nil)
    }

    /// Appends a column to the data frame.
    ///
    /// The column name is chosen in the following order and made unique:
    /// 1. The provided `name` parameter if non-nil.
    /// 2. The `examine.name` if non-nil.
    /// 3. An auto-generated name `Sample_###`.
    ///
    /// Tags are taken from `examine.tag` or default to an empty string.
    ///
    /// - Parameters:
    ///   - examine: The `SSExamine` column to append.
    ///   - name: An optional preferred column name.
    /// - Throws: Rethrows any error originating from underlying operations.
    public func append(_ examine: SSExamine<SSElement, FP>, name: String?) throws {
        self._data.append(examine)
        if let t = examine.tag {
            self._tags.append(t)
        }
        else {
            self._tags.append("")
        }
        let i = self._nCols + 1
        if let base = name ?? examine.name {
            let unique = makeUniqueColumnName(base: base, proposedIndex: i)
            self._columnNames.append(unique)
        }
        else {
            let auto = "Sample_" + String(format:"%03d", arguments: [i as CVarArg])
            let unique = makeUniqueColumnName(base: auto, proposedIndex: i)
            self._columnNames.append(unique)
        }
        self._nCols += 1
    }
    
    /// Removes all columns, names, and tags from the data frame.
    public func removeAll() {
        self._data.removeAll()
        self._tags.removeAll()
        self._columnNames.removeAll()
        self._nCols = 0
    }

    /// Saves the data frame to a JSON file, optionally wrapped in a ZIP archive.
    ///
    /// - Parameters:
    ///   - path: Destination path; `~` is expanded. If `compressAsZip == true` and the path
    ///     lacks a `.zip` extension, one is appended.
    ///   - atomically: When `true`, writes the temporary file atomically before the final move/replace.
    ///   - overwrite: When `false`, throws `.fileExists` if the destination already exists.
    ///   - compressAsZip: When `true`, stores the JSON payload inside a single-entry ZIP.
    ///   - zipEntryFileName: Optional inner entry name; defaults to basename-without-extension + `.json`.
    /// - Returns: `true` on success.
    /// - Throws: `SSError` for invalid paths, directory issues, encoding failures, or I/O errors.
    public func saveTo(
        fileName path: String,
        atomically: Bool = true,
        overwrite: Bool,
        compressAsZip: Bool = false,
        zipEntryFileName: String? = nil
    ) throws -> Bool {
        let fm = FileManager.default
        let expanded = NSString(string: path).expandingTildeInPath
        guard !expanded.contains("\0"), !expanded.isEmpty else {
            let msg = "Path is empty or contains NUL."
            SSLog.filesystemError(msg)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: msg)
        }

        var destURL = URL(filePath: expanded, directoryHint: .notDirectory)
        if compressAsZip && destURL.pathExtension.lowercased() != "zip" {
            destURL.appendPathExtension("zip")
        }
        let parentURL = destURL.deletingLastPathComponent()

        // Ensure parent exists (create if needed)
        var isDir = ObjCBool(false)
        if !fm.fileExists(atPath: parentURL.path, isDirectory: &isDir) {
            do {
                try fm.createDirectory(at: parentURL, withIntermediateDirectories: true)
            } catch {
                SSLog.filesystemError("Unable to create directory: \(parentURL.path)")
                throw SSError(type: .directoryDoesNotExist, file: #fileID, line: #line, function: #function, reason: error.localizedDescription)
            }
        } else if !isDir.boolValue {
            let msg = "Parent path is not a directory: \(parentURL.path)"
            SSLog.filesystemError(msg)
            throw SSError(type: .directoryDoesNotExist, file: #fileID, line: #line, function: #function, reason: msg)
        }

        let encoder = JSONEncoder()
        encoder.outputFormatting = []
        #if DEBUG
        encoder.outputFormatting.insert(.prettyPrinted)
        #endif
        let jsonData: Data
        do {
            jsonData = try encoder.encode(self)
        } catch {
            let msg = "JSON encoding failed: \(error.localizedDescription)"
            SSLog.filesystemError(msg)
            throw SSError(type: .wrongDataFormat, file: #fileID, line: #line, function: #function, reason: msg)
        }

        // Prepare temporary file alongside destination
        let tmpName = UUID().uuidString + (compressAsZip ? ".zip.tmp" : ".tmp")
        let tmpURL = parentURL.appending(path: tmpName, directoryHint: .notDirectory)

        do {
            if compressAsZip {
                let entry = makeZipEntryName(for: destURL, override: zipEntryFileName)
                try writeZip(at: tmpURL, entryName: entry, data: jsonData)
            } else {
                let opts: Data.WritingOptions = atomically ? [.atomic] : []
                try jsonData.write(to: tmpURL, options: opts)
            }
        } catch {
            SSLog.filesystemError("Failed writing temporary file: \(error.localizedDescription)")
            try? fm.removeItem(at: tmpURL)
            throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: error.localizedDescription)
        }

        do {
            if fm.fileExists(atPath: destURL.path) {
                guard overwrite else {
                    try? fm.removeItem(at: tmpURL)
                    let msg = "File exists and overwrite is false: \(destURL.lastPathComponent)"
                    SSLog.filesystemError(msg)
                    throw SSError(type: .fileExists, file: #fileID, line: #line, function: #function, reason: msg)
                }
                _ = try fm.replaceItemAt(destURL, withItemAt: tmpURL, backupItemName: nil, options: [.usingNewMetadataOnly])
            } else {
                try fm.moveItem(at: tmpURL, to: destURL)
            }
            return true
        } catch {
            try? fm.removeItem(at: tmpURL)
            SSLog.filesystemError("Unable to finalize write: \(error.localizedDescription)")
            throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: "Finalize failed: \(error.localizedDescription)")
        }
    }

    /// Loads a data frame from a JSON or single-entry ZIP file produced by `saveTo`.
    ///
    /// - Parameters:
    ///   - path: Source path; `~` is expanded. ZIP archives are detected via extension or magic bytes.
    ///   - stringEncoding: Encoding used when reading zipped text; defaults to UTF-8.
    /// - Returns: A decoded `DataFrame` on success.
    /// - Throws: `SSError` for I/O errors, invalid paths, or decoding failures.
    public class func dataframe(
        fromFile path: String,
        stringEncoding: String.Encoding = .utf8
    ) throws -> DataFrame<SSElement, FP>? {
        let fm = FileManager.default
        let expanded = NSString(string: path).expandingTildeInPath
        guard !expanded.contains("\0"), !expanded.isEmpty else {
            let msg = "Path is empty or contains NUL."
            SSLog.filesystemError(msg)
            throw SSError(type: .invalidArgument, file: #fileID, line: #line, function: #function, reason: msg)
        }
        let url = URL(filePath: expanded, directoryHint: .notDirectory)

        var isDir = ObjCBool(false)
        guard fm.fileExists(atPath: url.path, isDirectory: &isDir), !isDir.boolValue else {
            let msg = "File not found or is a directory: \(url.path)"
            SSLog.filesystemError(msg)
            throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: msg)
        }
        guard fm.isReadableFile(atPath: url.path) else {
            let msg = "File not readable: \(url.path)"
            SSLog.filesystemError(msg)
            throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: msg)
        }

        let data: Data
        if try Self.isZip(url: url) {
            let text = try Self.readFirstZipEntryString(url: url, encoding: stringEncoding)
            guard let d = text.data(using: stringEncoding) ?? text.data(using: .utf8) else {
                let msg = "Unable to transcode ZIP entry text to data."
                SSLog.filesystemError(msg)
                throw SSError(type: .wrongDataFormat, file: #fileID, line: #line, function: #function, reason: msg)
            }
            data = d
        } else {
            do {
                data = try Data(contentsOf: url, options: [.mappedIfSafe])
            } catch {
                let msg = "Unable to read file: \(error.localizedDescription)"
                SSLog.filesystemError(msg)
                throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: msg)
            }
        }

        do {
            let decoder = JSONDecoder()
            let decoded = try decoder.decode(DataFrame<SSElement, FP>.self, from: data)
            return decoded
        } catch {
            let msg = "JSON decode failed: \(error.localizedDescription)"
            SSLog.filesystemError(msg)
            throw SSError(type: .wrongDataFormat, file: #fileID, line: #line, function: #function, reason: msg)
        }
    }
    
    /// Legacy JSON archiver (non-zipped) preserved for compatibility.
    ///
    /// Prefer ``saveTo(fileName:atomically:overwrite:compressAsZip:zipEntryFileName:)`` for
    /// atomic writes, optional ZIP compression, and clearer error signalling.
    ///
    /// - Parameters:
    ///   - path: Destination file path; `~` is expanded.
    ///   - overwrite: When `false`, throws if the file exists.
    /// - Returns: `true` on success.
    /// - Throws: `SSError` on invalid paths, directory issues, or I/O failures.
    public func archiveTo(filePath path: String!, overwrite: Bool!) throws -> Bool {
        // Validate and normalize inputs
        guard let rawPath = path?.trimmingCharacters(in: .whitespacesAndNewlines), !rawPath.isEmpty else {
            #if os(macOS) || os(iOS)
            if #available(macOS 10.12, iOS 13, *) {
                SSLog.filesystemError("Empty file path")
            }
            #endif
            throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: "Empty file path")
        }
        let overwrite = overwrite ?? false
        
        let fm = FileManager.default
        let expandedPath = NSString(string: rawPath).expandingTildeInPath
        let fileURL = URL(filePath: expandedPath)
        let dirURL = fileURL.deletingLastPathComponent()
        
        // Ensure directory exists and is a directory
        var isDir = ObjCBool(false)
        guard fm.fileExists(atPath: dirURL.path, isDirectory: &isDir), isDir.boolValue else {
            #if os(macOS) || os(iOS)
            if #available(macOS 10.12, iOS 13, *) {
                SSLog.filesystemError("No writable directory found: \(dirURL.path)")
            }
            #endif
            throw SSError(type: .directoryDoesNotExist, file: #fileID, line: #line, function: #function, reason: "Directory does not exist: \(dirURL.path)")
        }
        
        // If the target exists
        if fm.fileExists(atPath: fileURL.path) {
            if overwrite {
                // Ensure deletable (covers permissions on parent directory)
                guard fm.isDeletableFile(atPath: fileURL.path) else {
                    #if os(macOS) || os(iOS)
                    if #available(macOS 10.12, iOS 13, *) {
                        SSLog.filesystemError("File is not deletable: \(fileURL.path)")
                    }
                    #endif
                    throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: "File not deletable")
                }
                do {
                    try fm.removeItem(at: fileURL)
                } catch {
                    #if os(macOS) || os(iOS)
                    if #available(macOS 10.12, iOS 13, *) {
                        SSLog.filesystemError("Unable to remove existing file: \(error.localizedDescription)")
                    }
                    #endif
                    throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: "Remove failed: \(error.localizedDescription)")
                }
            } else {
                #if os(macOS) || os(iOS)
                if #available(macOS 10.12, iOS 13, *) {
                    SSLog.filesystemError("File exists: \(fileURL.path)")
                }
                #endif
                throw SSError(type: .fileExists, file: #fileID, line: #line, function: #function, reason: "File exists at path")
            }
        } else {
            // Ensure parent directory is writable
            guard fm.isWritableFile(atPath: dirURL.path) else {
                #if os(macOS) || os(iOS)
                if #available(macOS 10.12, iOS 13, *) {
                    SSLog.filesystemError("Directory not writable: \(dirURL.path)")
                }
                #endif
                throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: "Directory not writable")
            }
        }
        
        // Encode JSON safely
        let encoder = JSONEncoder()
        encoder.outputFormatting = []
        #if DEBUG
        encoder.outputFormatting.insert(.prettyPrinted)
        #endif
        
        let data = try encoder.encode(self)
        
        // Write atomically
        do {
            try data.write(to: fileURL, options: [.atomic])
            return true
        } catch {
            #if os(macOS) || os(iOS)
            if #available(macOS 10.12, iOS 13, *) {
                SSLog.filesystemError("Unable to write data: \(error.localizedDescription)")
            }
            #endif
            throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: "Write failed: \(error.localizedDescription)")
        }
    }

    // MARK: ZIP helpers (mirrors SSExamine)

    /// Builds a reasonable default inner entry name for the zip file.
    /// Defaults to basename(without extension) + ".json" unless `override` is provided.
    private func makeZipEntryName(for destURL: URL, override: String?) -> String {
        if let o = override, !o.isEmpty { return o }
        let base = destURL.deletingPathExtension().lastPathComponent
        return base.isEmpty ? "dataframe.json" : base + ".json"
    }

    /// Writes a ZIP archive at `zipURL` with a single entry named `entryName` containing the given data.
    private func writeZip(at zipURL: URL, entryName: String, data: Data) throws {
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
                if position < 0 || position > Int64(data.count) { return Data() }
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
        if url.pathExtension.lowercased() == "zip" { return true }
        do {
            let fh = try FileHandle(forReadingFrom: url)
            defer { try? fh.close() }
            let sig = try fh.read(upToCount: 4) ?? Data()
            if sig.count == 4 {
                let magic: [UInt8] = [0x50, 0x4B, 0x03, 0x04]
                return Array(sig) == magic
            }
        } catch {
            SSLog.filesystemError("Unable to read file header for ZIP detection: \(error.localizedDescription)")
        }
        return false
    }

    /// Reads the first non-directory entry from a ZIP file as a String using `encoding`.
    private class func readFirstZipEntryString(url: URL, encoding: String.Encoding) throws -> String {
        let archive: Archive
        do {
            archive = try Archive(url: url, accessMode: .read)
        } catch {
            let msg = "Unable to open ZIP archive: \(error.localizedDescription)"
            SSLog.filesystemError(msg)
            throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: msg)
        }

        let fileEntries = archive.compactMap { entry -> Entry? in
            entry.type == .file ? entry : nil
        }
        guard let preferred = fileEntries.first else {
            let msg = "ZIP archive contains no file entries."
            SSLog.filesystemError(msg)
            throw SSError(type: .wrongDataFormat, file: #fileID, line: #line, function: #function, reason: msg)
        }

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

        if let s = String(data: data, encoding: encoding) {
            return s
        } else if let sUTF8 = String(data: data, encoding: .utf8) {
            return sUTF8
        } else {
            let msg = "Failed to decode ZIP entry text."
            SSLog.filesystemError(msg)
            throw SSError(type: .wrongDataFormat, file: #fileID, line: #line, function: #function, reason: msg)
        }
    }

    /// Legacy JSON unarchiver preserved for compatibility.
    ///
    /// Prefer ``dataframe(fromFile:stringEncoding:)`` for ZIP-awareness and clearer diagnostics.
    ///
    /// - Parameter path: Source path; `~` is expanded.
    /// - Returns: Decoded data frame or `nil` if parsing fails.
    /// - Throws: `SSError` on invalid paths or I/O/decoding errors.
    public class func unarchiveFrom(filePath path: String!) throws -> DataFrame<SSElement, FP>? {
        // Validate and normalize input
        guard let rawPath = path?.trimmingCharacters(in: .whitespacesAndNewlines), !rawPath.isEmpty else {
            #if os(macOS) || os(iOS)
            if #available(macOS 10.12, iOS 13, *) {
                SSLog.filesystemError("Empty file path")
            }
            #endif
            throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: "Empty file path")
        }
        
        let fm = FileManager.default
        let expandedPath = NSString(string: rawPath).expandingTildeInPath
        let fileURL = URL(filePath: expandedPath)
        
        // Ensure file exists and is readable
        var isDir = ObjCBool(false)
        guard fm.fileExists(atPath: fileURL.path, isDirectory: &isDir), !isDir.boolValue else {
            #if os(macOS) || os(iOS)
            if #available(macOS 10.12, iOS 13, *) {
                SSLog.filesystemError("Path does not point to a regular file: \(fileURL.path)")
            }
            #endif
            throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: "Not a regular file")
        }
        guard fm.isReadableFile(atPath: fileURL.path) else {
            #if os(macOS) || os(iOS)
            if #available(macOS 10.12, iOS 13, *) {
                SSLog.filesystemError("File not readable: \(fileURL.path)")
            }
            #endif
            throw SSError(type: .fileNotWriteable, file: #fileID, line: #line, function: #function, reason: "File not readable")
        }
        
        // Read and decode
        do {
            let data = try Data(contentsOf: fileURL, options: [.mappedIfSafe])
            let decoder = JSONDecoder()
            let result = try decoder.decode(DataFrame<SSElement, FP>.self, from: data)
            return result
        } catch let err as DecodingError {
            #if os(macOS) || os(iOS)
            if #available(macOS 10.12, iOS 13, *) {
                SSLog.filesystemError("JSON decoding failed: \(err)")
            }
            #endif
            throw SSError(type: .wrongDataFormat, file: #fileID, line: #line, function: #function, reason: "Decoding error: \(err)")
        } catch {
            #if os(macOS) || os(iOS)
            if #available(macOS 10.12, iOS 13, *) {
                SSLog.filesystemError("Unable to read file: \(error.localizedDescription)")
            }
            #endif
            throw SSError(type: .fileNotFound, file: #fileID, line: #line, function: #function, reason: "Read failed: \(error.localizedDescription)")
        }
    }
}

extension DataFrame where FP: RealLike, SSElement: RealLike {
    /// Mean of the per-row means across the entire data frame, or `nil` if any row has no mean.
    public func totalMean() -> FP? {
        var terms: [FP] = []
        for e in self.data {
            if let m = e.arithmeticMean {
                terms.append(m)
            }
            else {
                return nil
            }
        }
        let totalMean = Helpers.sum(&terms) / FP(terms.count)
        return totalMean
    }
    /// Median of the per-row medians across the entire data frame, or `nil` if any row has no median.
    public func totalMedian() -> FP? {
        var terms: [FP] = []
        for e in self.data {
            if let m = e.median {
                terms.append(m)
            }
            else {
                return nil
            }
        }
        
        let e: SSExamine<FP, FP> = SSExamine(usingArray: terms, name: "totalMedian", characterSet: nil)
        return e.median
    }

}
