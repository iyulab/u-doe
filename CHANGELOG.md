# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
