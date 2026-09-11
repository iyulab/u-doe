# u-doe

Design of Experiments (DOE) library for Rust — classical design generation, effect estimation, RSM, and multi-response optimization.

## Features

- **Design Generation**: Full factorial (2^k), fractional factorial (2^(k-p)), Plackett-Burman screening, Central Composite Design (CCD), Box-Behnken, Taguchi orthogonal arrays, Definitive Screening Design (DSD)
- **Analysis**: Effect estimation (main effects + 2FI), half-normal plot data, DOE ANOVA with F-statistics and p-values, Taguchi S/N ratio analysis
- **Response Surface Methodology**: Second-order OLS model fitting, steepest ascent path
- **Power & Sample Size**: Two-level factorial power, required replicates, power curves
- **Multi-response Optimization**: Derringer-Suich desirability functions
- **Coding**: Actual ↔ coded value transformation (Encoder)

## Quick Start

```toml
[dependencies]
u-doe = "0.2"
```

## Examples

### 2^4 Full Factorial Design

```rust
use u_doe::design::factorial::full_factorial;
use u_doe::analysis::effects::estimate_effects;

let design = full_factorial(4).unwrap(); // 16 runs, Yates order
let responses = vec![
    45.0, 71.0, 48.0, 65.0, 68.0, 60.0, 80.0, 65.0,
    43.0, 100.0, 45.0, 104.0, 75.0, 86.0, 70.0, 96.0,
];

let effects = estimate_effects(&design, &responses, 2).unwrap();
// Prints A = 21.625, AC = -18.125, ...
for e in &effects {
    println!("{}: {:.3}", e.name, e.estimate);
}
```

### Response Surface Methodology

```rust
use u_doe::design::ccd::{ccd, AlphaType};
use u_doe::analysis::rsm::{fit_rsm, steepest_ascent};

let design = ccd(2, AlphaType::Rotatable, 3).unwrap();
let responses = vec![/* measured values */];
// let model = fit_rsm(&design, &responses).unwrap();
// let path = steepest_ascent(&model, 5, 0.5);
```

### Multi-response Desirability

```rust
use u_doe::optimization::desirability::{ResponseSpec, overall_desirability};

let specs = vec![
    ResponseSpec::maximize(50.0, 100.0, 100.0, 1.0),
    ResponseSpec::minimize(0.0, 0.0, 30.0, 1.0),
];
let d = overall_desirability(&specs, &[80.0, 10.0]);
println!("Overall desirability: {d:.3}");
```

### Taguchi S/N Ratio

```rust
use u_doe::analysis::taguchi_sn::{signal_to_noise, sn_factor_effects, SnGoal};
use u_doe::design::factorial::full_factorial;

// Replicated responses per run (e.g., 3 replicates each)
let responses = vec![
    vec![10.0, 11.0, 10.5],  // run 1
    vec![20.0, 21.0, 19.5],  // run 2
    vec![15.0, 14.5, 15.5],  // run 3
    vec![25.0, 26.0, 24.5],  // run 4
];
let sn_values = signal_to_noise(&responses, SnGoal::LargerIsBetter).unwrap();

let design = full_factorial(2).unwrap();
let effects = sn_factor_effects(&design, &sn_values).unwrap();
for e in &effects {
    println!("Factor {}: optimal level {}, Δ = {:.2} dB",
             e.factor_index, e.optimal_level, e.delta);
}
```

### Power & Sample Size

```rust
use u_doe::power::{two_level_factorial_power, required_replicates};

// Power for a 2^3 full factorial, detecting effect of 2σ with α=0.05
let power = two_level_factorial_power(3, 0, 2, 2.0, 1.0, 0.05);
println!("Power = {power:.3}");  // ~0.95 at n=2 replicates

// Minimum replicates to achieve 80% power
let n = required_replicates(3, 0, 2.0, 1.0, 0.05, 0.80, 10);
println!("Required replicates: {n}");
```

## Design Types

| Design | Function | Use Case |
|--------|----------|----------|
| Full Factorial | `full_factorial(k)` | k = 2..7, all factor combinations |
| Fractional Factorial | `fractional_factorial(k, p)` | Screening, k=4..7, p=1..3 |
| Plackett-Burman | `plackett_burman(k)` | Screening, k ≤ 19, N = 8/12/16/20 |
| CCD | `ccd(k, alpha, n_center)` | RSM, k = 2..6 |
| Box-Behnken | `box_behnken(k, n_center)` | RSM, k = 3/4/5 |
| Taguchi | `taguchi_array(name, k)` | Robust design, L4–L27 |
| DSD | `definitive_screening(k)` | Screening + quadratic, k = 2..12, 2k+1 runs |
| Simplex Lattice | `simplex_lattice(q, m)` | Mixture, q = 2..12, degree m = 1..10 |
| Simplex Centroid | `simplex_centroid(q)` | Mixture, q = 2..12, 2^q − 1 subset centroids |

## References

- Montgomery, D.C. (2019). *Introduction to Statistical Quality Control*, 8th ed. Wiley.
- Derringer, G. & Suich, R. (1980). Simultaneous Optimization of Several Response Variables. *Journal of Quality Technology* 12(4), 214–219.
- Plackett, R.L. & Burman, J.P. (1946). The Design of Optimum Multifactorial Experiments. *Biometrika* 33(4), 305–325.
- Jones, B. & Nachtsheim, C.J. (2011). A Class of Three-Level Designs for Definitive Screening in the Presence of Second-Order Effects. *Journal of Quality Technology* 43(1), 1–15.
- Taguchi, G. (1986). *Introduction to Quality Engineering*. Asian Productivity Organization.

## WebAssembly / npm

Available as an npm package via [wasm-pack](https://rustwasm.github.io/wasm-pack/).

```bash
npm install @iyulab/u-doe
```

### Quick Start

```javascript
import init, { full_factorial, doe_anova } from '@iyulab/u-doe';

await init();
const design = full_factorial(3); // 2^3 = 8 runs
```

> **Note:** All analysis functions take **native JS arrays/objects** — pass values
> directly, not `JSON.stringify(...)` strings. A string argument is rejected with a
> descriptive error naming the offending parameter.

### Functions

#### `full_factorial(k) -> DesignMatrix`

Generate a 2^k full factorial design (k = 1..7).

**Output:**
```json
{ "data": [[...]], "factor_names": ["X1","X2"], "run_count": 4, "factor_count": 2 }
```

#### `fractional_factorial(k, p) -> DesignMatrix`

Generate a 2^(k-p) fractional factorial design (k = 4..7, p = 1..3).

#### `fractional_factorial_info(k, p) -> FractionalInfo`

Metadata for the fraction `fractional_factorial(k, p)` emits: resolution, defining relation, and generator equations. Use it to derive the alias structure of the design you actually received instead of pairing the matrix with an external published table.

**Output:**
```json
{
  "k": 7,
  "p": 3,
  "resolution": "IV",
  "defining_relation": "I=ABCE=BCDF=ACDG=ADEF=BDEG=ABFG=CEFG",
  "generators": ["E=ABC", "F=BCD", "G=ACD"]
}
```

#### `plackett_burman(k) -> DesignMatrix`

Generate a Plackett-Burman screening design (k = 1..19).

#### `ccd(k, design_type, n_center) -> DesignMatrix`

Generate a Central Composite Design. `design_type`: `"FaceCentered"` | `"Rotatable"` | `"Inscribed"`. k = 2..6.

#### `box_behnken(k, n_center) -> DesignMatrix`

Generate a Box-Behnken design (k = 3, 4, or 5).

#### `taguchi_array(name, k) -> DesignMatrix`

Get a Taguchi orthogonal array. `name`: `"L4"` | `"L8"` | `"L9"` | `"L12"` | `"L16"` | `"L18"` | `"L27"`.

#### `definitive_screening(k) -> DesignMatrix`

Generate a Definitive Screening Design (k = 2..12).


#### `simplex_lattice(q, m) -> DesignMatrix`

Simplex lattice design for a `q`-component mixture, degree `m`. Rows are the
compositions whose components are multiples of `1/m` and sum to 1.

`q` must be in `2..=12` and `m` in `1..=10`; anything else is rejected with an
error rather than returning a degenerate design.

```js
simplex_lattice(3, 2)
// { data: [[0,0,1],[0,0.5,0.5],[0,1,0],[0.5,0,0.5],[0.5,0.5,0],[1,0,0]],
//   factor_names: ["X1","X2","X3"], run_count: 6, factor_count: 3 }
```

#### `simplex_centroid(q) -> DesignMatrix`

Simplex centroid design: the centroid of every non-empty subset of the `q`
components, giving `2^q − 1` runs (q=3 → 7, q=4 → 15). `q` must be in `2..=12` —
the run count doubles with every additional component.

#### `doe_anova(design, responses, factor_names, effect_names) -> AnovaResult`

Perform DOE ANOVA on a coded design matrix.

**Input:** `design`: `[[f64]]`, `responses`: `Float64Array`, `factor_names`: `["A","B"]`, `effect_names`: `["A","B","A:B"]`

Interaction effect names join factor names with `":"` (e.g. `"A:B"`). An unknown entry in `effect_names` is an error (since 0.5.0; previously silently skipped).

Since 0.10.0 three further inputs are errors rather than misleading numbers:

- the design is not two-level coded (any value other than `-1`/`+1`) — use `fit_rsm` for designs with centre or axial points;
- two requested effects share a contrast column, so they are aliased and would contribute the same sum of squares twice;
- the model asks for more terms than the design has degrees of freedom. An exactly saturated model stays legal, with `residual_df: 0`.
- two requested effects are partially aliased — their contrasts are correlated without being identical, as a Plackett-Burman two-factor interaction is with other main effects. Their sums of squares would overlap, so the table would not be an ANOVA of the data.

**Output:**
```json
{ "effects": [{ "name": "A", "sum_of_squares": 10.0, "df": 1, "mean_square": 10.0, "f_statistic": 5.0, "p_value": 0.03 }], "residual_ss": 4.0, "residual_df": 2, "total_ss": 14.0, "r_squared": 0.71, "r_squared_adj": 0.57, "fitted": [12.1, 14.3, ...], "residuals": [0.4, -0.3, ...], "curvature": null, "pure_error": null }
```

#### `signal_to_noise(responses, goal) -> [f64]`

Compute Taguchi S/N ratios. `responses`: `[[f64]]` (replicates per run). `goal`: `"LargerIsBetter"` | `"SmallerIsBetter"` | `"NominalIsBest"`.

#### `estimate_effects(design, responses, factor_names, max_order) -> EffectsResult`

Estimate main effects and interactions for a 2-level factorial design.

**Output:**
```json
{ "effects": [
    { "name": "A",   "columns": [0],    "estimate": 21.6,  "sum_of_squares": 1870.6, "percent_contribution": 45.2 },
    { "name": "A:C", "columns": [0, 2], "estimate": -18.1, "sum_of_squares": 1314.1, "percent_contribution": 31.7 }
  ],
  "half_normal": [{ "term_index": 1, "abs_effect": 18.1, "quantile": 0.57 },
                  { "term_index": 0, "abs_effect": 21.6, "quantile": 1.15 }] }
```

Interaction names join factor names with `":"` (since 0.5.0; previously bare concatenation `"AC"`). `columns` holds the design-matrix column indices of the term's factors — use it for display formatting (e.g. `"A × C"`) instead of parsing `name`.

Since 0.10.0 a design that is not two-level coded is an error. The contrast for a term is the product of its factor columns, which estimates an effect only when every column is `-1` or `+1`: a centre point zeroes that product and an axial point scales it. Central composite, Box-Behnken, definitive screening and three-level Taguchi arrays (L9/L18/L27) therefore belong in `fit_rsm`; two-level arrays (L4/L8/L12/L16), full and fractional factorials and Plackett-Burman designs are unaffected.

Each `half_normal` point carries `term_index` — an index into `effects` — so a point can be labelled directly (e.g. `effects[point.term_index].name`). The points are sorted by `|effect|`, a **different order** from `effects` (model-term order), so pairing them positionally (`effects[i]` ↔ `half_normal[i]`) mislabels every point; always use `term_index`.

The result also carries `lenth: { pse, margin_of_error, df, distinct_contrasts }` — Lenth's (1989) pseudo standard error and margin of error, the line to draw on the half-normal plot: effects beyond `margin_of_error` are judged active. It needs no replication, which is the point of it for saturated screening designs. Terms that share a contrast column (in a 2^(5-2), `A`, `B:D` and `C:E`) are counted once. `lenth` is `null` when there are fewer than 3 distinct contrasts.

`max_order` goes up to the factor count (`4` on a 2^4 includes `A:B:C:D`); a higher order adds nothing, and `0` is an error.

`doe_anova` returns per-run `fitted` and `residuals`, in run order, for residuals-versus-fitted and normal probability plots.

`doe_anova` accepts **centre points** — runs with every factor at 0. They carry no information about the effects, which come from the factorial runs alone; what they measure is curvature. The result then carries `curvature: { sum_of_squares, df: 1, f_statistic, p_value }`, tested against `pure_error: { sum_of_squares, df }` — the spread among runs made at the same design point, centre replicates and repeated factorial runs alike. Both are `null` when the design gives nothing to compute them from. Axial and three-level designs still belong in `fit_rsm`.

#### `fit_rsm(design, responses, factor_names) -> RsmModel`

Fit a second-order Response Surface Model via OLS.

**Output:**
```json
{ "coefficients": [10.0, 2.5, -1.3, 0.5, 0.8, -0.2], "r_squared": 0.95, "factor_count": 2 }
```

#### `steepest_ascent(coefficients, factor_count, n_steps, step_size) -> AscentResult`

Compute the steepest ascent path from a fitted RSM model.

**Output:**
```json
{ "steps": [{ "coded": [0.5, 0.3], "step_number": 1 }] }
```

#### `desirability(specs, responses) -> DesirabilityResult`

Compute Derringer-Suich desirability for multiple responses.

**Input:** `specs`: `[{ "goal": "Maximize"|"Minimize"|"Target", "lower": 0, "target": 100, "upper": 100, "s1": 1, "s2": 1, "importance": 1 }]`, `responses`: `Float64Array`

`s1`/`s2` are curve-shape exponents (`s = 1` linear, `> 1` convex/stricter, `< 1` concave). `importance` (optional, default `1`) is the Derringer-Suich response weight rᵢ — it is **distinct from** the shape exponents: raise `importance` to make a response count more in the aggregate, not to reshape its curve.

**Output:**
```json
{ "individual": [0.8, 0.6], "overall": 0.69 }
```

`overall` is the importance-weighted geometric mean D = (∏ dᵢ^rᵢ)^(1/Σrᵢ); with all weights at the default `1` this is the plain geometric mean.

#### `two_level_factorial_power(k, p, n_replicates, effect_size, sigma, alpha) -> f64`

Compute statistical power of a 2^(k-p) factorial design. Returns power in [0, 1].

## npm (WebAssembly)

```bash
npm install @iyulab/u-doe
```

The package resolves per environment via a conditional `exports` map:

| Environment | Entry |
|---|---|
| Bundlers (webpack, Vite, …) | ESM + WebAssembly ESM-integration (`default` condition) |
| Node.js — `require()`, ESM `import`, CJS TS runners (`tsx`, `ts-node`) | CJS glue loading the wasm from the filesystem (`node` condition) — no loader hooks or flags |

## Related

- [`u-analytics`](https://crates.io/crates/u-analytics) — SPC, process capability, statistical analysis
- [`u-numflow`](https://crates.io/crates/u-numflow) — Math primitives (used internally)
