# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
