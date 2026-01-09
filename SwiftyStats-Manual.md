# SwiftyStats Manual

SwiftyStats is a concurrency-safe statistics toolkit for Swift packages, apps, and scripts. This manual explains how to install the library, understand the module layout, operate the `SSExamine` container, call the inferential and MLE helpers, drive the CLI demo, and extend the package with your own models. The intent is to provide a single, task-oriented reference so you can move from ingestion to inference without hunting through source files.

## Table of Contents

1. [Purpose and Scope](#purpose-and-scope)
2. [Installing and Building](#installing-and-building)
3. [How to Use SwiftyStats (Cheat Sheet)](#how-to-use-swiftystats-cheat-sheet)
4. [Package Tour](#package-tour)
5. [Working with SSExamine](#working-with-ssexamine)
6. [Derived Views and Sample Inspection](#derived-views-and-sample-inspection)
7. [Summary Statistics](#summary-statistics)
8. [Entropy and Complexity](#entropy-and-complexity)
9. [Outliers and Hypothesis Testing](#outliers-and-hypothesis-testing)
10. [Maximum-Likelihood Fitters](#maximum-likelihood-fitters)
11. [Custom MLE Problems and Solvers](#custom-mle-problems-and-solvers)
12. [Random Sampling Utilities](#random-sampling-utilities)
13. [Score Functions for Custom Gradients](#score-functions-for-custom-gradients)
14. [CLI Demo and Automation](#cli-demo-and-automation)
15. [Testing and Quality Assurance](#testing-and-quality-assurance)
16. [Troubleshooting and Recipes](#troubleshooting-and-recipes)

---

## Purpose and Scope

- Highlight every public entry point exposed by the Swift Package.
- Provide runnable snippets you can paste into Swift Playgrounds, unit tests, or scripts.
- Describe preconditions, concurrency guarantees, and error surfaces so you can compose APIs safely.
- Outline best practices for persistence, logging, and troubleshooting.

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
        .package(url: "https://github.com/volker/SwiftyStatsGH.git", branch: "main")
    ],
    targets: [
        .target(name: "MyApp", dependencies: [
            .product(name: "SwiftyStats", package: "SwiftyStatsGH")
        ])
    ]
)

