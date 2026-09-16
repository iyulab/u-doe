# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.14.0] - 2026-09-16

### Added

- **`optimize_desirability(specs, candidates)`** (WASM and Rust) -- the search
  over candidate settings that Derringer & Suich's method ends in. The crate
  scored one response vector and stopped there, so every consumer wrote the loop
  itself, and the loop has a trap the scores do not point out.

  **An overall desirability of 0 is not an optimum.** D is a weighted geometric
  mean, so one response that misses its limits at every candidate makes every
  candidate score exactly 0; a loop keeping the maximum then returns whichever
  candidate it scored first -- a setting that satisfies no specification,
  presented as the recommended one. The result is therefore two shapes rather
  than one with a sentinel:

  ```js
  { kind: "best", index, overall, individual }
  { kind: "infeasible", unreachable: [number] }
  ```

  `unreachable` names the responses that score 0 everywhere, and is empty when
  each response is met somewhere but never all at once -- the distinction
  between moving a limit and relaxing the trade-off. Ties go to the first
  strictly best candidate; a response that is not finite is refused rather than
  scored as 0.

- **`rsm_predict(coefficients, factor_count, coded)`** (WASM) -- evaluates a
  fitted response surface at coded factor levels, batched: `coded` is a flat
  `n x factor_count` row-major array and the result is `n` responses. The crate
  has always had `RsmModel::predict`; only the binding was missing, so a JS
  consumer wanting a contour, an optimum or a confirmation-run estimate had to
  rebuild the polynomial and copy a term order that is private to the crate.

- **`fit_rsm`'s result carries `terms`** -- what each coefficient is, in the
  same order: `["Intercept", "A", "B", "A^2", "B^2", "A:B"]`. Without it the
  order was knowable only from the crate's source, so a caller displaying
  coefficients copied a layout it could not check against the value it was
  handed.

- **`analysis::rsm::model_terms`** (Rust) -- the same list, built beside the
  function that builds the model row so the labels and the order they describe
  cannot drift apart.

### Fixed

- **`doe_anova` takes every effect name `estimate_effects` returns.** It
  estimated contrasts up to order 2 whatever was asked for, so a caller who
  requested three-factor interactions got `A:B:C` back from `estimate_effects`
  and had it refused here as `unknown_effect` -- with a message naming
  `:`-joined factor names as exactly the accepted form. Every `2^3` and larger
  full factorial has such a term, so the natural pipeline (estimate effects at
  an order, then ANOVA on those effects) ended in an error at order 3.

  The order is now read off the requested names. An over-long or misspelled
  name is still `unknown_effect`; `fit_least_squares` already accepted any
  order and is unchanged.

## [0.13.0] - 2026-09-15

### Changed

- **Breaking:** every WASM function throws an `Error` carrying `code` and the
  values behind the refusal (`err.terms`, `err.runs`, `err.first`,
  `err.supported`, ...) instead of a bare string; `err.message` keeps the
  readable text. The codes and fields are listed in the README.
- **Breaking:** `DoeError` carries fields instead of text. The catch-all
  `UnsupportedDesign(String)`, `InvalidSpecification(String)` and
  `MatrixError(String)` are gone, replaced by `UnsupportedFraction`,
  `UnknownArray`, `EmptyDesign`, `ParameterOutOfRange`, `TooFewContrasts`,
  `EffectOutsideDesign`, `InvalidCodingRange`, `EmptyResponses`,
  `TooFewReplicates`, `NonPositiveResponse` and `SingularModel`.
  `InsufficientResponses` is renamed `ResponseCountMismatch` (too many
  responses are refused too), and `UnknownEffect { name }` is
  `UnknownEffect { effect }`. `DoeError::code()` returns the same `code` the
  WASM error carries.
- **Breaking:** `analysis::rsm::steepest_ascent` and `RsmModel::predict` return
  a `Result`.
- **Breaking:** `power::two_level_factorial_power` returns `Result<f64, _>`,
  `power_curve` returns `Result<Vec<_>, _>`, and `required_replicates` returns
  `Result<Option<usize>, _>` -- `None` when no replicate count up to
  `max_replicates` reaches the target.
- **Breaking:** `ResponseSpec` can only be built valid. Its fields are private
  (read them with `goal()`, `lower()`, `target_value()`, `upper()`, `s1()`,
  `s2()`, `importance()`); `ResponseSpec::new`, `maximize`, `minimize`,
  `target` and `with_importance` return a `Result`, and
  `overall_desirability` returns `Result<f64, DoeError>`.
- **Breaking:** `EffectRow::f_statistic` and `LeastSquaresTerm`'s
  `f_statistic`, `std_error` and `t_statistic` are `Option<f64>`, `None`
  where they used to be `NaN`: with no residual degrees of freedom, or (for
  the F and t statistics) a residual of exactly zero -- the same convention
  `CurvatureTest::f_statistic` and every `p_value` already followed. Over WASM
  they are `null` instead of `NaN`.
- **Breaking:** WASM results carry `null` for an absent value (`curvature`,
  `pure_error`, `lenth`, `p_value`), as documented, instead of `undefined`.
  Test with `== null`.

### Added

- `design::factorial::standard_fractions()` and the WASM function
  `standard_fractions()`: every 2^(k-p) fraction the crate can build, with its
  resolution, defining relation and generators.
- `DesignMatrix::check_shape()`.

### Fixed

- `desirability` accepted a specification whose ramp has no width or runs
  backwards and scored it anyway: `Maximize` with `target <= lower` gave a
  step at `lower` or a constant 1, `Target` with its limits reversed gave 0
  everywhere -- scores that look valid and rank every candidate the same. The
  limits a goal uses must now be finite with `lower < target` (Maximize),
  `target < upper` (Minimize) or `lower < target < upper` (Target), the
  exponents it uses finite and positive, and `importance` finite and not
  negative; otherwise `InvalidDesirabilityLimits` or `ValueOutOfDomain`
  (WASM: `invalid_desirability_limits` / `value_out_of_domain` with the spec's
  `index`).
- `overall_desirability` returned 0.0 for no specifications, a response count
  that did not match, or weights summing to zero -- indistinguishable from a
  design with a response out of range. These are now `EmptyResponses`,
  `ResponseCountMismatch` and `NoWeightedResponse`.
- `power::two_level_factorial_power` returned 0.0 for inputs that have no power
  to report (`k = 0`, `p >= k`, no replicates, a non-positive effect size or
  sigma) and did not check `alpha` at all (`alpha = 2` gave power 1). It now
  refuses them with `ParameterOutOfRange` or `ValueOutOfDomain`; the WASM
  function throws instead of returning a number. Valid inputs give the same
  values as before.
- `2^(k-p)` in the power calculation was an integer shift, which overflows once
  `k - p` reaches the width of `usize` -- 32 bits in WebAssembly, where the
  shift was silently masked and the power computed for the wrong run count.
- `power::required_replicates` returned `max_replicates` both when the target
  power was first reached there and when it was not reached at all.
- `RsmModel::predict` predicted from part of the model when the point or the
  coefficient list had the wrong length. It now refuses both.
- A run whose length differs from the number of factor names is refused
  (`DesignShapeMismatch`). It used to index past the row -- a panic, which a
  WASM caller receives as `RuntimeError: unreachable` -- or, for a short row of
  zeros in `doe_anova`, pass for a centre point.
- `fit_rsm` refuses an empty design and a design with fewer runs than the
  quadratic model has coefficients, instead of panicking or reporting a
  matrix error; `steepest_ascent` refuses a coefficient list of the wrong
  length instead of panicking.
- `fit_least_squares` reported `OverSpecifiedModel` with the mean counted as a
  term ("model has 8 terms but the design has only 6 degrees of freedom" for
  seven terms on seven runs). `terms` now counts the terms, as `doe_anova`
  always did.
- `fit_least_squares` named a term constant over the runs as
  `effects 'I' and 'C'` rather than "effect 'C' is aliased with the overall
  mean": the mean now goes second, as everywhere else in the crate.

## [0.12.1] - 2026-09-15

### Fixed

- Power and sample-size calculations use `u-numflow` 0.6's normal quantile and
  CDF: the quantile was accurate only to 4.5e-4 (Abramowitz & Stegun 26.2.23)
  and is now exact to double precision (Wichura AS 241), so computed power
  changes in its trailing digits.

## [0.12.0] - 2026-09-13

### Added

- `analysis::least_squares::fit_least_squares` -- ordinary least squares on
  two-level coded data with Type III (partial) sums of squares, for the runs
  the contrast formula cannot analyse: a run left out, corners replicated
  unequally, terms that are correlated. Each term reports its coefficient,
  effect (twice the coefficient), standard error, t and F statistics and
  p-value; the fit reports residual, R², fitted values and residuals. On
  balanced data it returns the same values as `doe_anova` and
  `estimate_effects`. Requested columns that coincide over the runs kept, or a
  term constant over them, are refused as aliased with each other or with the
  mean. Also exposed as the WASM function `fit_least_squares`.

## [0.11.0] - 2026-09-12

### Added

- `analysis::lenth::lenth` -- Lenth's (1989) pseudo standard error and margin
  of error for unreplicated and saturated designs, where there are no residual
  degrees of freedom to test against. Terms that share a contrast column are
  counted once: in a 2^(5-2), `A`, `B:D` and `C:E` are one contrast, and
  counting each would repeat the same magnitude and inflate the degrees of
  freedom. Only the individual margin is reported; the simultaneous margin is
  quoted with different quantiles by different sources.
- `DoeAnovaResult::fitted` and `DoeAnovaResult::residuals`, per run, for
  residuals-versus-fitted and normal probability plots. With orthogonal
  contrasts each coefficient is half its effect, so no matrix solve is needed.
- WebAssembly: `doe_anova` returns `fitted` and `residuals`; `estimate_effects`
  returns `lenth: { pse, margin_of_error, df, distinct_contrasts }` (or `null`
  with fewer than 3 distinct contrasts) -- the line a half-normal plot is read
  against.
- `doe_anova` accepts centre points -- runs with every factor at 0. The effects
  come from the factorial runs alone, and `DoeAnovaResult::curvature` reports
  the curvature test the centre points exist for: a one-degree-of-freedom sum of
  squares, `n_F · n_C · (ȳ_F − ȳ_C)² / (n_F + n_C)`, tested against the pure
  error. 0.10.0 refused such designs, and `analysis::rsm::fit_rsm`, which its
  error message pointed to, cannot fit them either: with only three levels the
  pure quadratic columns coincide.
- `DoeAnovaResult::pure_error`, from every group of identical design rows --
  centre points and replicated factorial runs alike.
- WebAssembly: `doe_anova` returns `curvature` and `pure_error` (or `null`).

### Changed

- **Breaking:** `doe_anova` refuses effects whose contrasts are correlated
  without being identical, with `DoeError::PartiallyAliasedEffects`. It built
  the model sum of squares by adding per-term sums of squares, which is exact
  only for orthogonal contrasts; a Plackett-Burman two-factor interaction is
  correlated with other main effects, so the table was not an ANOVA of the data
  and the residual could be clamped to zero without notice.
- **Breaking:** `estimate_effects` applies the same check, to every term up to
  `max_order`: a contrast that does not sum to zero is refused as partially
  aliased with the mean (`I`), and two that are correlated without coinciding
  are refused as partially aliased with each other. The contrast formula it
  uses is exact only for balanced, orthogonal columns, and a design with a run
  left out or unequal replication breaks both -- the effects, sums of squares
  and percent contributions it returned looked plausible and were not. For a
  replicated 2^3 with one run removed, the sums of squares added up to 3.6
  times the total. Terms whose columns coincide exactly, as aliases in a regular
  fraction do, are still returned and `lenth` counts them once; a 12-run
  Plackett-Burman is refused at `max_order` 2 and accepted at 1. `lenth` checks
  the effects it is given the same way.
- **Breaking:** `DoeAnovaResult` has four new public fields (`fitted`,
  `residuals`, `curvature`, `pure_error`). With centre points, `r_squared`
  counts the curvature term as part of the model.

### Fixed

- Several error messages contained runs of spaces in mid-sentence, left by line
  continuations that had lost their backslash. `NotTwoLevelCoded` now also
  tells a caller with centre points to use `doe_anova` rather than `fit_rsm`.

### Removed

- The dependency on `u-analytics`. It was declared, and kept current through
  several version bumps, but nothing in the crate used it; it only added a
  crate to every build and tied this crate's releases to that one's.

## [0.10.0] - 2026-09-10

### Changed

- **Breaking:** `analysis::effects::estimate_effects` and `analysis::anova::doe_anova`
  now reject a design that is not two-level coded, with the new
  `DoeError::NotTwoLevelCoded`. Both build the contrast for a term as the product
  of its factor columns, which estimates an effect only when every column is
  `-1` or `+1`. A centre point zeroes that product for every term it touches and
  an axial point scales it, so the `2/n` divisor stops matching the runs that
  contributed — the functions returned numbers that looked like effects and were
  not. Designs carrying centre or axial points (central composite, Box-Behnken,
  definitive screening, three-level Taguchi arrays such as L9/L18/L27) belong in
  `analysis::rsm::fit_rsm`; two-level arrays (L4/L8/L12/L16), full and fractional
  factorials and Plackett-Burman designs are unaffected.
- **Breaking:** `doe_anova` rejects two requested effects that share a contrast
  column, with the new `DoeError::AliasedEffects`. In a fractional design this is
  the ordinary case rather than an exotic one — in a 2^(5-2), `A`, `B:D` and
  `C:E` are one column — and admitting both counted the same sum of squares
  twice, which previously surfaced only as a residual silently clamped to zero.
- **Breaking:** `doe_anova` rejects a model with more terms than the design has
  degrees of freedom, with the new `DoeError::OverSpecifiedModel`. An exactly
  saturated model stays legal: its residual degrees of freedom are zero, which is
  the normal outcome for a screening design.
- `doe_anova` computes residual degrees of freedom by subtraction rather than
  saturating subtraction. With the three guards above in place the value can no
  longer go negative, so saturating the result would only hide a bug.

### Added

- `design::DesignMatrix::two_level_violation` reports the first entry that is
  neither `-1` nor `+1`, as `(run, factor, value)`, so a caller can test a design
  before handing it to effect estimation.

## [0.9.0] - 2026-09-07

### Changed

- **Breaking:** `design::mixture::simplex_lattice` and `simplex_centroid` now
  return `Result<DesignMatrix, DoeError>` instead of `DesignMatrix`, matching
  every other design generator in the crate. Their doc comments already stated
  the contract (`q >= 2`, `m >= 1`); nothing enforced it.
  - `simplex_lattice(0, 2)` used to abort — reaching WebAssembly callers as a
    raw `unreachable` trap rather than a domain error.
  - `simplex_lattice(q, 0)` used to divide by zero and return a design matrix of
    `NaN`, which serialises to `null` and reached the caller looking
    structurally valid. This was the more dangerous of the two.
  - `simplex_centroid(q)` computes `2^q - 1` runs by shifting, which overflowed
    for large `q`.
- Component count is now bounded to `2..=12` and lattice degree to `1..=10`,
  reported through the existing `InvalidFactorCount` and `UnsupportedDesign`
  errors.

### Added

- Documentation for both mixture designs in the README, including the
  WebAssembly signatures and returned JSON shape.
- Boundary regression tests, plus an invariant that accepted designs never
  contain non-finite values.

## [0.8.1] - 2026-09-07

### Changed

- **`rand` is now 0.10** and **`getrandom` 0.4** on WebAssembly targets. This
  crate does not name `rand` types in its public signatures, so the change is
  internal and the API is unaffected. The
  `RUSTFLAGS --cfg getrandom_backend="wasm_js"` that `getrandom` 0.3 required is
  no longer needed.
- **`u-numflow` is now required at 0.4 and `u-analytics` at 0.7** (previously 0.3
  and 0.6), following those crates' own releases.
- **The minimum supported Rust version is now declared as 1.85** and is verified
  by building on that exact toolchain; 1.84 and below fail. The crate previously
  declared no `rust-version` at all.

## [0.8.0] - 2026-07-17

### Changed

- **BREAKING: `fractional_factorial(7, 3)` now emits the published standard
  fraction.** The 2^(7-3) generators change from the non-standard
  `E=ABD, F=ACD, G=BCD` to the textbook `E=ABC, F=BCD, G=ACD`
  (Montgomery 2019, Table 8.14; NIST/SEMATECH e-Handbook §5.3.3.4.7),
  matching every other supported (k, p). Both fractions are valid
  resolution IV designs, but the previous one silently diverged from the
  published alias tables consumers pair these designs with. Runs stored
  from an earlier 7-3 design do not match the new matrix.

### Fixed

- **`fractional_factorial_info(7, 3)` returned a mathematically invalid
  defining relation.** The stored string
  `I=ABDE=ACDF=AEFG=BCFG=BDEG=CDEF=ABCG` was not closed under word
  products (e.g. ABDE·ACDF = BCEF, absent), so it described *no*
  2^(7-3) fraction — neither the emitted one nor the standard one. The
  relation is now `I=ABCE=BCDF=ACDG=ADEF=BDEG=ABFG=CEFG`, and new
  table-wide invariant tests pin every entry's defining relation to the
  emitted design matrix (each word's column product is +1 on every run,
  2^p−1 distinct words, resolution = shortest word length), so a
  generators ↔ relation mismatch can no longer pass CI.

### Added

- **`FractionalInfo.generators`** — generator equations of the emitted
  fraction (e.g. `["E=ABC", "F=BCD", "G=ACD"]`), derived from the same
  numeric table that builds the design matrix so the strings cannot
  drift. Consumers can now derive the alias structure of the design they
  actually received instead of pairing it with an external published
  table.
- **WASM `fractional_factorial_info(k, p)`** — exposes the metadata to
  npm consumers: `{ k, p, resolution: "III"|"IV"|"V", defining_relation,
  generators: [str] }`.
- **`fractional_factorial(6, 3)`** — the standard 8-run resolution III
  fraction (`D=AB, E=AC, F=BC`, Montgomery Table 8.14) was missing from
  the generator table and rejected as unsupported.
- `Resolution` now implements `Display` (`"III"` / `"IV"` / `"V"`).

## [0.7.0] - 2026-07-15

### Added

- **`desirability` — response importance weights.** `ResponseSpec` gains an
  `importance` field (default `1.0`) and a `with_importance()` builder, and
  `overall_desirability` now computes the canonical Derringer-Suich
  importance-weighted geometric mean D = (∏ dᵢ^rᵢ)^(1/Σrᵢ). With all weights at
  the default this is identical to the previous unweighted geometric mean. The
  WASM `desirability` spec accepts an optional `importance` field. This separates
  the response weight from the curve-shape exponents `s1`/`s2`, which were
  previously the only per-response knob and were easy to mistake for importance.

### Changed

- **`half_normal_plot_data` / WASM `estimate_effects` — term identity.**
  `half_normal_plot_data` now returns `Vec<HalfNormalPoint>` (`term_index`,
  `abs_effect`, `quantile`) instead of bare `(f64, f64)` tuples, and the WASM
  `half_normal` output is now `[{ term_index, abs_effect, quantile }]`. The points
  are sorted by `|effect|`, a different order from `effects` (model-term order);
  the new `term_index` lets callers label each point correctly instead of pairing
  positionally, which mislabelled every point.

## [0.6.2] - 2026-07-05

### Fixed

- npm: expose the `./package.json` subpath in the `exports` map so tools
  that `require('<pkg>/package.json')` (license scanners, version
  reporters) keep working alongside the conditional exports introduced in
  the previous release (`ERR_PACKAGE_PATH_NOT_EXPORTED`).

## [0.6.1] - 2026-07-05

### Fixed

- **npm packaging — Node-compatible entry.** The npm package previously
  shipped only the wasm-bindgen *bundler*-target output, whose static
  `.wasm` import fails on Node's CJS path (`tsx`/`ts-node` in non-ESM
  packages) with an opaque `SyntaxError: Invalid or unexpected token`.
  The package now additionally ships the *nodejs*-target CJS glue under
  `node/` and routes Node consumers to it via a conditional `exports`
  map (`node` → CJS with filesystem wasm loading, `default` → bundler
  ESM). `require()`, native ESM `import`, and CJS TS runners all work
  without loader hooks. A pre-publish smoke test (CJS `require` + ESM
  `import`) now guards this path in CI. Rust API unchanged.

### Changed

- `u-numflow` dependency `^0.2` → `^0.3` (compatible; 0.3.0 publishes the
  previously-unreleased `wasm` feature and input-validation hardening —
  no API used by this crate changed).

## [0.6.0] - 2026-06-12

### Changed — BREAKING (WASM)

- WASM `desirability` response specs now **reject unknown keys** with an
  explicit `unknown field` error instead of silently ignoring them
  (`serde(deny_unknown_fields)`). Other exports take positional arrays/strings
  and are unaffected.

### Changed

- Dependency: `u-analytics` `^0.5` → `^0.6`.

## [0.5.1] - 2026-06-11

### Changed

- Dependency refresh: `u-analytics` `^0.4` → `^0.5` (upstream adds
  `BinMethod::Fixed`; no u-doe API or behavior change).

## [0.5.0] - 2026-06-10

### Changed — BREAKING

- **Interaction effect names now join factor names with `":"`** (R/JMP convention):
  `estimate_effects` and `doe_anova` emit `"A:B"` instead of the bare concatenation
  `"AB"`. Bare concatenation is ambiguous for multi-character factor names
  (`"Injection pressure" + "Hold pressure"` → unsplittable token).
- **`doe_anova` errors on unknown effect names** (`DoeError::UnknownEffect`) instead
  of silently skipping them. Silently dropping a term produced a structurally
  different model with no signal to the caller. Callers using the old `"AB"`-style
  names get an actionable error instead of an empty/smaller model.
- **`EffectEstimate` gains a `columns: Vec<usize>` field** (design-matrix column
  indices of the term's factors). Code constructing `EffectEstimate` literally must
  add the field. The WASM `estimate_effects` response includes `columns` per effect —
  consumers can format display labels (e.g. `"A × B"`) without parsing `name`.

### Changed — WASM / npm

- Analysis exports (`doe_anova`, `estimate_effects`, `signal_to_noise`, `fit_rsm`,
  `steepest_ascent`, `desirability`) drop the legacy `_json` parameter-name suffix
  in the generated `.d.ts` — they have taken native JS arrays/objects (not JSON
  strings) since 0.4.0, matching the mixture exports.
- Passing a JSON **string** to an analysis export now fails fast with a descriptive
  error naming the parameter (`"design: expected a native JS array/object, got a
  string — pass the value directly, not JSON.stringify(...)"`) instead of an opaque
  `TypeError: Reflect.get called on non-object` from deep inside serde.

## [0.4.0] - 2026-05

- Mixture designs (Simplex Lattice + Centroid) + WASM bindings.
- Analysis exports migrated from JSON-string to native-value (serde-wasm-bindgen)
  convention.
