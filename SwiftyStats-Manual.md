# SwiftyStats Manual

SwiftyStats is a concurrency-safe statistics toolkit for Swift packages, apps, and scripts. This manual covers installation, module layout, and common usage patterns so you can move from ingestion to inference without hunting through source files.

## Table of Contents

1. [Purpose and Scope](#purpose-and-scope)
2. [Requirements](#requirements)
3. [Installing and Building](#installing-and-building)
4. [Quick Start](#quick-start)
5. [Module Overview](#module-overview)
6. [CLI Demo](#cli-demo)
7. [Documentation](#documentation)
8. [Testing and Quality Assurance](#testing-and-quality-assurance)

---

## Purpose and Scope

- Highlight key entry points across SSExamine, Inferential, MLE, and TimeSeries.
- Provide runnable snippets for Playgrounds, tests, or scripts.
- Summarise preconditions and error surfaces so you can compose APIs safely.

---

## Requirements

- Swift 6.2 or later
- macOS 13 / iOS 16 (or newer) SDKs

---

## Installing and Building

### Add the package to `Package.swift`

```swift
// swift-tools-version: 6.2
import PackageDescription

let package = Package(
    name: "MyApp",
    platforms: [.macOS(.v13), .iOS(.v16)],
    products: [
        .library(name: "MyApp", targets: ["MyApp"])
    ],
    dependencies: [
        .package(url: "https://github.com/volker/SwiftyStatsGH.git", from: "2.0.0")
    ],
    targets: [
        .target(name: "MyApp", dependencies: [
            .product(name: "SwiftyStats", package: "SwiftyStatsGH")
        ])
    ]
)
```

Then run:

```bash
swift package resolve
swift build
```

---

## Quick Start

### SSExamine basics

```swift
import SwiftyStats

let examine = SSExamine<Int, Double>()
examine.append(42)
examine.append(repeating: 3, item: 7)

let items = examine.itemsAsArray(sorted: .original)
let freq = examine.frequency(7)
let rel = examine.relativeFrequency(7)
```

### Hypothesis testing

```swift
import SwiftyStats

let a = [1.0, 2.0, 3.0, 4.0, 5.0]
let b = [10.0, 11.0, 12.0, 13.0, 14.0]

let mw = try NonParametric<Double>.mannWhitneyUTest(
    group1: a,
    group2: b,
    continuityCorrection: true
)
let pValue = mw?.pValueAsymTwoSided
```

### MLE fitters

```swift
import SwiftyStats

let data = [1.2, 0.9, 1.1, 1.3, 0.95]
let result = MLEFitter<Double>.fitGamma(data: data)
let thetaHat = result.thetaHat
```

### Time series helpers

```swift
import SwiftyStats

let values = (1...10).map { Double($0) }
let times = (0..<10).map { Double($0) }
let ts = TimeSeries(times: times, values: values)

let smoothed = ts.singleExponentialSmoothing(alpha: 0.5)
let acf = ts.autocorrelationFunction(maxLag: 4)
```

---

## Module Overview

- `SSExamine`: exploratory statistics, frequency tables, entropy, and `SSDataFrame` utilities.
- `Inferential`: hypothesis tests (parametric, non-parametric, goodness-of-fit, outliers, variance tests).
- `MLE`: maximum-likelihood fitters plus an internal solver (`MLEProblem`, `MLESolver`).
- `TimeSeries`: smoothing, autocorrelation, and Durbin-Watson helpers.
- `LinearModel`: minimal OLS helper via `GLMFit.fit`.

Note: `SwiftyStatsPrelude` is an internal target that reexports shared dependencies for SwiftyStats sources. It is not exposed as a product dependency for external packages.

---

## CLI Demo

Run the demo executable to smoke-test behaviour:

```bash
swift run swiftystats-demo
```

---

## Documentation

- Task-oriented references live in `README.md`, `SwiftyStats-Manual.md`, and `Using-MLE.md`.
- DocC output is generated from symbol documentation in `Sources/SwiftyStats` using `make doc`.

---

## Testing and Quality Assurance

```bash
swift test
swift test --enable-code-coverage
```

Use `make test` if you prefer the Makefile shortcuts.
