# Changelog

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [2.0.1] - 2026-01-14
### Added
- Optional `continuityCorrection` flag for `mannWhitneyUTest` to enable continuity-corrected normal approximations.
- Optional seed parameter for `ksTestOneSample` to make KS bootstrap sampling reproducible.
- Rank now surfaces per-tie-block metadata (`TieBlock`) containing start index, length, value, and exact `t^3 - t` correction.
- Minimal OLS helper `GLMFit.fit` (pivoted QR with Accelerate, normal-equation fallback) akin to R's `lm.fit`, returning coefficients, fitted values, residuals, rank, and pivot indices.
- Pan-modified Farebrother implementation for Durbin-Watson p-values in `Sources/SwiftyStats/TimeSeries/DurbinWatsonPan.swift`, exposing quadratic-form CDFs and raw `TimeSeries` entry points under the MIT license.
- Durbin-Watson statistic helpers and exact Imhof-based p-values with SwiftyBoost quadrature integration.
- Swift Testing coverage for HypothesisTesting public entry points and result accessors.
- DocC-style documentation for the Wald–Wolfowitz runs test in `Randomness.swift`, covering cutting-point selection, tie handling, and the approximate versus exact p-values returned.
- DocC coverage for the two-sample Kolmogorov-Smirnov test, clarifying inputs, statistic computation, and invalid-data guards.
- Swift Testing coverage for `TimeSeries` smoothing and autocorrelation using a 25-point deterministic dataset.
- DocC-style documentation for the Bartlett variance homogeneity test, including the statistic and chi-square reference distribution.
- DocC coverage for variance tests (Bartlett, Levene, F-ratio, chi-squared) and the sign test overloads to clarify inputs, statistics, and return values.
- R reference script `R-Scripts/DurbinWatsonReference.R` to reproduce Durbin-Watson statistics and Imhof-based p-values via CompQuadForm and lmtest, including dwtest support for arbitrary time series vectors passed on the command line.
- DocC-style coverage for the Levene variance homogeneity test, detailing the absolute-deviation transformation and F-statistic.
- DocC-style coverage for the one-sample chi-squared variance test, describing the `(n - 1) * s2 / s0` statistic and chi-square reference distribution.
- R script `R-Scripts/DurbinWatsonValidation.R` to mirror `DurbinWatsonValidation.runValidation()` output using CompQuadForm’s `davies`.
- DurbinWatsonTest convenience helpers to derive residuals and design matrices from raw time series (intercept-only or intercept + trend) and to run the DW test directly on a series.
- R script `R-Scripts/TukeyKramerReference.R` to compute Tukey-Kramer post-hoc comparisons for a fixed eleven-group dataset.

### Changed
- Updated `README.md`, `SwiftyStats-Manual.md`, and `Using-MLE.md` to align module descriptions, installation guidance, and MLE options naming with current sources.
- Corrected MLE fitter documentation to reference `MLEOptimizationOpts` and fix `OptimizerKind` typos.
- Updated `.gitignore` to drop tooling-specific ignore entries while keeping local workspace and secret exclusions.
- Added file-level DocC summaries across all Swift sources and completed public symbol documentation (TimeSeries, inferential result structs, score helpers, SSDataFrame) for consistent API coverage.
- Completed DocC coverage for remaining public hypothesis-testing entry points and result accessors (Kruskal-Wallis, Friedman, Mann-Whitney exact p-values, rank helpers, and standard deviation aliases).
- Documented `Parametric.twoSampleTTest` with DocC summary covering inputs, tail thresholds, and the embedded Levene variance check.
- Bartlett and Levene variance-equality helpers now stream pooled variance and absolute-deviation statistics to avoid per-observation buffers on large datasets, with a large-sample Levene Swift Testing case for regression coverage.
- Split the monolithic `SSExamine.swift` into semantic slices (`SSExamine-core`, `-location-measures`, `-scatter-measures`, `-shape-parameters`, `-moments`, `-concentration-measures`, `-entropy-measures`, and `-file-management`) to keep accessors grouped by category.
- Split `SSExamine-internal.swift` into focused impl files (`SSExamine-internal-moments-concentration`, `-entropy-measures`, `-items`, `-frequencies`, `-dispersion`, `-accumulators`) under `SSExamine/impl` for clearer maintenance.
- `R-Scripts/DurbinWatsonValidation.R` now mirrors the Swift validation output, uses `dwtest` for cross-checks, and falls back to `imhof` when `davies` reports an error for indefinite forms.
- Wilcoxon signed-rank test accepts configurable zero handling (`discardZeroDifferences` or Pratt-style zero retention) to align statistics with R when zero differences are present.
- Added concise DocC summaries for Wilcoxon signed-rank overloads documenting zero handling, continuity correction, and exact p-value options.

### Fixed
- KS bootstrap now honours the caller-supplied `bsCount` across all distributions.
- KS bootstrap now reports the number of successful replicates and derives p-values from that count.
- KS statistic now validates CDF outputs, clamping to [0, 1] and rejecting NaN/inf values.
- Noncentral Fisher F KS bootstrap now re-estimates parameters per replicate.
- KS binomial trial-count estimation now uses a MoM seed with discrete log-likelihood search instead of `max(data)`.
- KS one-sample statistic now handles tied values by evaluating EDF left/right limits per unique value.
- StudentT KS bootstrap no longer applies hidden standardisation in `bootstrapCore`; standardisation now happens explicitly in the StudentT path.
- KS noncentral chi-squared bootstrap now derives a data-driven df guess and returns `[k, lambda]` so the KS CDF/sampler stay consistent.
- Tukey-Kramer post-hoc testing now uses the within-group mean square error and symmetric standard error in test statistics and confidence intervals.
- Documentation now matches current API names and signatures across DocC, README, and the manuals (SSExamine sorting, entropy helpers, MLE surface, RNG sampler, and inferential types).
- Sign test now reports the correct two-sided exact binomial p-value (including large samples), applies the symmetric normal approximation with continuity correction, guards all-tie inputs, and uses the effective non-tied sample size in its results.
- Rank tie-correction terms now computed in integer space before conversion to maintain exact `t^3 - t` values.
- Corrected Mann–Whitney U distribution PMF/CDF dynamic programming to normalise by the total combinations and initialise base cases properly.
- Anderson-Darling normality test now applies the Stephens adjustment for estimated mean/variance, uses the normal-specific Stephens p-value approximation (matching R's `ad.test`) instead of the generic AD CDF, and clamps the result to [0, 1] so normality decisions reflect the correct tail behaviour.
- Durbin-Watson p-value computation now evaluates upper and lower tails separately (upper tail via `4 - d` symmetry) to align with dwtest tail probabilities.
- Corrected the `SSExamineLocationDispersionTests` Gini mean difference expectation by computing pairwise absolute differences directly so the expected value matches the implementation.
- Levene variance-equality test now accumulates within-group squared absolute deviations correctly and is covered by a mean-centred Swift Testing case.
- F-ratio variance test now reports two-sided and one-sided p-values consistently with the CDF tails, sets directional/equality flags accordingly, and derives two-sided/one-sided confidence intervals from the correct quantiles.
- Balanced F-ratio power search now returns the first sample size meeting the target within bounds (or `nil` when unreachable), chi-squared variance equality flags now reflect non-rejection, and Bartlett/Levene tests guard against being called with fewer than two groups.

## [2.0.0] - 2025-11-22
### Added
- Swift Testing coverage for box-whisker summaries plus the Rousseeuw-Croux `Qn` and `Sn` robust scale estimators.
- Swift Testing coverage for `DataFrame` creation, mutation, copy/Codable round-trips, and archive/unarchive file management.
- `DataFrame.saveTo` and `DataFrame.dataframe(from:)` for JSON persistence with optional ZIP compression, mirroring the `SSExamine` file-management API.
- Convenience `Makefile` with developer targets: `doc`, `doc-archive`, `clean-docs`, `test`, `debug`, `release`, and `clean`.
- Repository documentation: `README.md`, `LICENSE.md`, and `CHANGELOG.md`.
- Comprehensive Swift Testing coverage for `SSExamine`, including frequency tables, cumulative caches, metadata, and append variants.
- Additional `SSExamine` tests covering cumulative frequency snapshots, numeric projections, hashing, dispersion aliases, semivariances, confidence intervals, coefficient-of-variation helpers, dispersion edge cases, and moment/shape metrics.
- `HyndmanFanQuantileType` enumeration plus `SSExamine.quantileHyndmanFan` helpers implementing the nine Hyndman–Fan quantile definitions.
- Python utility `entropy_measures.py` that reproduces the public entropy APIs for a reference dataset and persists the results to `entropy_results.json`.
- RNG-driven Swift Testing coverage for the public `MLEFitter`, exercising exponential, gamma, and logistic estimators via deterministic samples from `RNGSampler`.
- `MLEFitter.fitHyperexponential`, reusing the analytic Hyperexponential score with a logits/log-rate parameterisation to estimate mixture probabilities and phase rates under the simplex constraint.
- Richer `MLEResult` diagnostics: u-space uniqueness counts remain, θ-space uniqueness now scales by feasible ranges, Hessian condition numbers are estimated via Gershgorin discs, and convergence reasons are surfaced alongside the updated CLI/docs for reproducibility.
- Multi-start orchestration upgrades: deterministic Sobol (scrambled VDC) sampling for coverage, TaskGroup-backed parallel local solves on Apple platforms, and configurable pruning/resampling for invalid starts.
- Parameter spec heuristics tightened for Beta/Gamma/Fisher families, and the CLI demo now uses shared option presets (central, heavy-tail, noncentral) with moment-based warm starts for inverse-gamma and related examples.
- `MLEFitter.fitBernoulli` plus `RNGSampler.RandomNumbers.bernoulli` so Bernoulli data can be sampled and fit consistently alongside existing discrete families.
- Log-parameter score helpers (`scoreWithLog`/`totalScoreWithLog`) for the Landau, Map Airy, and Negative Binomial distributions so callers can work in unconstrained space without hand-coded chain rules.
- Dispersion internals now emit `SSLog.statisticsError` messages whenever computations return `nil`, mirroring the existing Gini coefficient diagnostics.

### Changed
- Completed DocC coverage for remaining public surfaces (DataFrame legacy archive helpers, inferential enums, experimental Welford accumulator, and score helpers) to keep symbol docs exhaustive.
- Cleaned demo output and source strings to keep the compiled code ASCII-only ahead of the 2.0.0 release.
- Completed DocC coverage for SSExamine inequality metrics and entropy helpers, documenting Gini/HHI/CR surfaces and the Approximate/Sample entropy kernels for clearer mathematical intent.
- Added concise, math-focused DocC overviews across `SSExamine/public` and `SSExamine/impl` to explain how dispersion, entropy, summation, and container helpers map to the empirical distribution.
- Documented Inferential result structs and demo harness helpers so previously undocumented properties and sample runners now surface DocC guidance.
- Completed DocC coverage for SSExamine cumulative frequency snapshots and dispersion internals (normalized Gini, concentration ratios, HHI) to describe their mathematical definitions and return conventions.
- Added DocC-style summaries for every remaining public declaration spanning GOF tests, MLE optimisation options, RNG helpers, score functions, and SSExamine accessors so symbol documentation is now exhaustive.
- Completed DocC symbol coverage for remaining public enums and error overrides so the library's API documentation now spans every public function and property.
- Documentation generation now supports both static hosting (`make doc`) and `.doccarchive` output (`make doc-archive`), and the README describes each target.
- `SSExamine` internals are refactored into focused impl files under `Sources/SwiftyStats/SSExamine/impl` to improve readability and future maintenance.
- Split the monolithic `SSExamine.swift` into semantic slices (`SSExamine-location*.swift`, `SSExamine-dispersion*.swift`, `SSExamine-items*.swift`, `SSExamine-summation*.swift`, and `SSExamine-shared.swift`) with dedicated internal helpers.
- Test target now depends on `SwiftyBoost` directly to assert analytic confidence interval bounds.
- Expanded inline documentation throughout `SSExamine-stats-dispersion.swift` to clarify range, semi-variance, and moment helpers.
- `SSError` now captures the logged failure reason and surfaces it via `localizedDescription`/`localizedFailureReason`, giving downstream callers actionable error messages.
- Updated `SSExamineLocationDispersionTests` expectations (quartiles, IQR, Gastwirth) to align with the new Hyndman–Fan Type 7 default quantile behavior while keeping explicit Type 2 coverage.
- `Inferential.HypothesisTesting.bootstrapOneSampleKS` now accepts sampler closures with the `(n: Int, dist: RandomNumbers, seed: UInt64?) -> [T]` shape used by `RNGSampler.randoms`, eliminating the explicit `inout RNG` requirement and making it easier to pipe in `MLEFitter.fit…` results and deterministic seeds.
- Normalized the `KSTest` helpers by renaming the local CDF closures and `ksStatistic` parameters to lowerCamelCase, making the file's identifiers consistent with Swift naming conventions.
- Ensured every fitter under `Sources/SwiftyStats/MLE/Fitters` uses the same doc-comment structure (Landau, MapAiry, Inverse Normal, Pareto, Triangular) and refreshed the narrative descriptions of their parameterizations.
- Added DocC-friendly documentation for every score extension plus the modified Bessel helpers in `Score-functions.swift`, clarifying how the distribution gradients and supporting derivatives are computed for the noncentral families.
- `MLEFitter.fitBernoulli` now emits the analytic sample-mean estimator directly and the triangular param-spec builder exposes all three parameters, improving solver stability for those shapes.
- Expanded `MLEFitterTests` to cover Bernoulli, Binomial, Geometric, Poisson, Negative Binomial, Pareto, inverse Chi-Squared, inverse Normal, Log-Normal, Landau, Map Airy, SAS 0.5, and Triangular fitters so every public MLE surface is exercised.
- Quantile and concentration-ratio helpers now throw `SSError(.invalidArgument, ...)` for out-of-range inputs, with call sites updated to use `try`/`try?` so invalid parameters no longer silently return `nil`.
- Wald-Wolfowitz two-sample runs test now counts runs over the sorted combined series (including inter-group ties), applies the requested continuity correction, uses the supplied significance level to derive critical values and decisions, and reports inter-group tie blocks and counts.

### Fixed
- Swift Testing `#require` calls in `VarianceTests` now wrap throwing helpers with `try` so the test target builds under SwiftPM.
- One-sample t-test now populates the correct 99% critical value, sets the alpha-specific critical value, and documents the entry point.
- Two-sample t-test now applies two-tailed alpha thresholds when setting equality/direction flags so decisions align with the reported p-values and critical values.
- Stripped non-ASCII characters from compiled sources and test names to avoid toolchain issues.
- `SSExamine.isNumeric` now preserves its stored heuristic flag for empty datasets and text-append workflows, keeping metadata defaults and `logSum`'s empty-case return of `0` intact.
- `SSExamine` append helpers now correctly report success by returning `true` when data is appended.
- `SSExamine` dispersion helpers now guard zero-median quartile distances, compute standard error via `s / √n`, evaluate skewness/kurtosis using sample standard deviation, and stabilise raw/high-order moment calculations to avoid NaNs for odd powers.
- Restored Hyndman–Fan type 6 quantile behaviour after moving the location helpers into split files.
- `SSExamine.empiricalCDF` now consults cumulative frequency caches so queries between observed values return the correct cumulative probability.
- `SSExamine.EntropyMeasures` now falls back to categorical coding when numeric conversion fails and recomputes sample deviation from the embedded numeric data, ensuring multiscale entropy stays finite for string-based samples.
- Tsallis entropy near the Shannon limit uses a numerically stable series expansion, preventing precision loss for `q ≈ 1`.
- Numerical MLE fitters now declare their `logPDF`/`gradLogPDF` closures as `@Sendable`, quieting Swift 6 data-race diagnostics across every distribution.
- Hyperexponential MLE wrapper now forwards the robust covariance and diagnostic metadata (θ-uniqueness, condition estimates, convergence reasons) reported by the core solver.
- Reimplemented `lorenzCurveLocked` so Lorenz curves are only produced for finite, non-negative samples and always terminate exactly at `(1, 1)`, logging diagnostics instead of trapping on invalid data.

## [0.1.0] - 2025-10-31
### Added
- Initial `SSExamine` reference implementation with concurrency-safe sampling, frequency tracking, and Codable support.
- `swiftystats-demo` executable showcasing sample ingestion and concurrent appends.
- Stress tests covering concurrent mutation and encoding/decoding behaviour.
