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
| Fractional Factorial | `fractional_factorial(k, p)` | Screening, k=3..7, p=1..3 |
| Plackett-Burman | `plackett_burman(k)` | Screening, k ≤ 19, N = 8/12/16/20 |
| CCD | `ccd(k, alpha, n_center)` | RSM, k = 2..6 |
| Box-Behnken | `box_behnken(k, n_center)` | RSM, k = 3/4/5 |
| Taguchi | `taguchi_array(name, k)` | Robust design, L4–L27 |
| DSD | `definitive_screening(k)` | Screening + quadratic, k = 2..12, 2k+1 runs |

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

### Functions

#### `full_factorial(k) -> DesignMatrix`

Generate a 2^k full factorial design (k = 1..7).

**Output:**
```json
{ "data": [[...]], "factor_names": ["X1","X2"], "run_count": 4, "factor_count": 2 }
```

#### `fractional_factorial(k, p) -> DesignMatrix`

Generate a 2^(k-p) fractional factorial design (k = 4..7, p = 1..3).

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

#### `doe_anova(design, responses, factor_names, effect_names) -> AnovaResult`

Perform DOE ANOVA on a coded design matrix.

**Input:** `design`: `[[f64]]`, `responses`: `Float64Array`, `factor_names`: `["A","B"]`, `effect_names`: `["A","B","AB"]`

**Output:**
```json
{ "effects": [{ "name": "A", "sum_of_squares": 10.0, "df": 1, "mean_square": 10.0, "f_statistic": 5.0, "p_value": 0.03 }], "residual_ss": 4.0, "residual_df": 2, "total_ss": 14.0, "r_squared": 0.71, "r_squared_adj": 0.57 }
```

#### `signal_to_noise(responses, goal) -> [f64]`

Compute Taguchi S/N ratios. `responses`: `[[f64]]` (replicates per run). `goal`: `"LargerIsBetter"` | `"SmallerIsBetter"` | `"NominalIsBest"`.

#### `estimate_effects(design, responses, factor_names, max_order) -> EffectsResult`

Estimate main effects and interactions for a 2-level factorial design.

**Output:**
```json
{ "effects": [{ "name": "A", "estimate": 21.6, "sum_of_squares": 1870.6, "percent_contribution": 45.2 }], "half_normal": [[21.6, 1.15]] }
```

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

**Input:** `specs`: `[{ "goal": "Maximize"|"Minimize"|"Target", "lower": 0, "target": 100, "upper": 100, "s1": 1, "s2": 1 }]`, `responses`: `Float64Array`

**Output:**
```json
{ "individual": [0.8, 0.6], "overall": 0.69 }
```

#### `two_level_factorial_power(k, p, n_replicates, effect_size, sigma, alpha) -> f64`

Compute statistical power of a 2^(k-p) factorial design. Returns power in [0, 1].

## Related

- [`u-analytics`](https://crates.io/crates/u-analytics) — SPC, process capability, statistical analysis
- [`u-numflow`](https://crates.io/crates/u-numflow) — Math primitives (used internally)
