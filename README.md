# SwiftyStats

SwiftyStats is a Swift Package that delivers a concurrency-safe toolkit for inspecting univariate samples - making heavy use of generics. The core `SSExamine` type tracks raw insertion order, absolute and relative frequencies, and produces cumulative views, making it useful for exploratory statistics in apps, scripts, or server tooling. A lightweight CLI (`swiftystats-demo`) ships with the package for quick experimentation.

## Release Status

SwiftPM versions are published via git tags. Replace `2.0.1` with the tag you want to depend on.

## Features

- Thread-safe sample container (`SSExamine`) with raw-order reconstruction, frequency tables, and Codable support.
- `SSDataFrame` for lightweight tabular storage with JSON persistence and optional ZIP compression.
- Inferential tests (parametric, non-parametric, variance, and outlier checks) under `Inferential`.
- MLE fitters and a shared solver for common distributions, with analytic and numerical options.
- Time series smoothing, autocorrelation, and exact Durbin-Watson helpers.
- Integration with [SwiftyBoost](https://github.com/strike65/SwiftyBoost.git) for advanced numeric routines.
- Ready-to-run demo executable plus a Swift Testing suite covering SSExamine, MLE, and time series helpers.

## Modules

- `SSExamine`: exploratory statistics, entropy, frequency tables, and `SSDataFrame` utilities.
- `Inferential`: hypothesis tests (parametric, non-parametric, goodness-of-fit, outliers, variance tests).
- `MLE`: maximum-likelihood fitters (`MLEFitter`) plus an internal solver (`MLEProblem`, `MLESolver`).
- `TimeSeries`: smoothing, autocorrelation, and Durbin-Watson helpers.
- `LinearModel`: minimal OLS helper via `GLMFit.fit`.

## Requirements

- Swift 6.2 or later
- macOS 13 / iOS 16 or newer SDKs

## Getting Started

Add SwiftyStats to your `Package.swift`:

```swift
.package(url: "https://github.com/volker/SwiftyStatsGH.git", from: "2.0.1")
```

Then add `SwiftyStats` to your target dependencies.

### Quick Sample

```swift
import SwiftyStats

let examine = SSExamine<Int, Double>()
examine.append(42)
examine.append(repeating: 3, item: 7) // add seven three times

print(examine.itemsAsArray(sorted: .original))     // [42, 7, 7, 7]
print(examine.frequency(7))                        // 3
print(examine.relativeFrequency(7))                // 0.75
print(examine.frequencyTable(sorted: .valueAscending))
```

### Hypothesis Testing (Example)

```swift
import SwiftyStats

let a = [1.0, 2.0, 3.0, 4.0, 5.0]
let b = [10.0, 11.0, 12.0, 13.0, 14.0]

let mw = try NonParametric<Double>.mannWhitneyUTest(group1: a, group2: b, continuityCorrection: true)
print(mw?.pValueAsymTwoSided)
```

### Demo Executable

```bash
swift run swiftystats-demo
```

## Development

The repository includes a `Makefile` to streamline common workflows:

```bash
make debug        # swift build
make release      # swift build -c release
make test         # swift test
make doc          # generate static DocC site under ./docs for static hosting
make doc-archive  # emit SwiftyStats.doccarchive into ./docs for IDE previews
make clean-docs   # remove generated docs, archives, and symbol graphs
make clean        # swift package clean
```

`make doc` produces a static DocC site ready for GitHub Pages (`docs/`), while `make doc-archive` keeps the `.doccarchive` format for local browsing in Xcode. Use `make clean-docs` to clear generated documentation assets.

## Documentation

- Generate static docs with `make doc`; the site is written to `docs/`. The hosting base path is configured via `HOSTING_BASE_PATH` in `Makefile` and should match your GitHub Pages project name.
- Publish by committing `docs/` and enabling GitHub Pages for the `main/docs` source, or preview locally with `python3 -m http.server --directory docs 8000`.
- Task-oriented references live in `SwiftyStats-Manual.md` and `Using-MLE.md`, while DocC output is generated from symbol documentation under `Sources/SwiftyStats`.

## Testing

All tests live under `Tests/SwiftyStatsTests` and use Swift's `Testing` framework (`@Test`, `#expect`). Run the entire suite with:

```bash
swift test
# or
make test
```

## Contributing

Issues and pull requests are welcome. Please include relevant tests and follow the Swift API Design Guidelines described in the project documentation.

## Licence

SwiftyStats is available under the MIT Licence. See [LICENSE.md](LICENSE.md) for details.
