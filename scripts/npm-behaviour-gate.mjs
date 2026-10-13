// Binding behaviour gate: errors and absent values as a JS caller sees them.
//
// An import-only smoke test proves the package loads, not that any binding
// answers with the documented shape. The host-target tests never build a JS
// `Error` and serialize with serde_json, so what is only observable against
// the built package is whether a refusal arrives as an `Error` carrying
// `code` and its fields, and an absent value as `null` rather than
// `undefined` or `NaN`.
//
//   PKG_NAME=<package-name> node scripts/npm-behaviour-gate.mjs
//
// Run it from a directory where the package is installed -- the publish
// workflow's smoke directory, or any directory with a local build under
// node_modules. The package is resolved from the working directory, not from
// this file. The publish workflow runs it after the README examples and before
// `npm publish`; keeping the assertions in a file rather than in the workflow
// means the same checks can run against a local build before the push that
// publishes, so an assertion the engine no longer satisfies is found there.

import { createRequire } from 'node:module'
import { join } from 'node:path'
import { pathToFileURL } from 'node:url'

const pkg = process.env.PKG_NAME
if (!pkg) {
  console.error('usage: PKG_NAME=<package-name> node npm-behaviour-gate.mjs')
  process.exit(2)
}
const requireFromCwd = createRequire(join(process.cwd(), 'noop.js'))
const m = await import(pathToFileURL(requireFromCwd.resolve(pkg)).href)
let checks = 0
const ok = (c, msg) => { checks++; if (!c) throw new Error(msg) }
const refusal = (f) => { try { f() } catch (e) { return e } return null }

// A refusal is an Error with a stable code and the values behind it.
const pb = m.plackett_burman(7)
let e = refusal(() => m.fit_least_squares(pb.data.slice(1), [1, 2, 3, 4, 5, 6, 7], pb.factor_names, pb.factor_names))
ok(e instanceof Error, 'a refusal must be an Error')
ok(e.code === 'over_specified_model' && e.terms === 7 && e.runs === 7, `code=${e.code} terms=${e.terms} runs=${e.runs}`)
ok(typeof e.message === 'string' && e.message.includes('degrees of freedom'), `message=${e.message}`)

const f3 = m.full_factorial(3)
e = refusal(() => m.fit_least_squares(f3.data.filter((r) => r[2] > 0), [1, 2, 3, 5], f3.factor_names, f3.factor_names))
ok(e.code === 'aliased_effects' && e.first === 'C' && e.second === 'I', `code=${e.code} first=${e.first} second=${e.second}`)

e = refusal(() => m.fractional_factorial(4, 2))
const table = m.standard_fractions()
ok(e.code === 'unsupported_fraction' && JSON.stringify(e.supported) === JSON.stringify(table.map((t) => [t.k, t.p])),
  'unsupported_fraction must list the table standard_fractions() returns')

const spec = (goal, lower, target, upper) => ({ goal, lower, target, upper, s1: 1, s2: 1 })
e = refusal(() => m.desirability([spec('Maximize', 0, 10, 10), spec('Target', 10, 5, 0)], [5, 5]))
ok(e.code === 'invalid_desirability_limits' && e.index === 1 && e.goal === 'Target', `code=${e.code} index=${e.index}`)

e = refusal(() => m.two_level_factorial_power(3, 0, 2, 1, 1, 2))
ok(e.code === 'parameter_out_of_range' && e.parameter === 'alpha' && e.exclusive_min === 0 && e.exclusive_max === 1 && e.got === 2,
  `code=${e.code} parameter=${e.parameter} exclusive_min=${e.exclusive_min} exclusive_max=${e.exclusive_max} got=${e.got}`)

// Power is the noncentral t test with the model's error df (0.648 here;
// the normal approximation said 0.760), and an unreplicated full model
// has no error df to test against.
const pw = m.two_level_factorial_power(3, 0, 2, 2.0, 1.5, 0.05)
ok(Math.abs(pw - 0.6480) < 1e-3, `two_level_factorial_power 2^3 x 2 = ${pw}, want 0.648`)
e = refusal(() => m.two_level_factorial_power(3, 0, 1, 2.0, 1.5, 0.05))
ok(e.code === 'no_error_degrees_of_freedom' && e.terms === 7 && e.runs === 8, `code=${e.code}`)

// Inputs that once trapped the module are refused like any other.
e = refusal(() => m.fit_rsm([], [], []))
ok(e instanceof Error && e.code === 'empty_input' && e.parameter === 'design', `fit_rsm([]): ${e.code} ${e.parameter}`)
ok(typeof m.full_factorial(2).run_count === 'number', 'the module must still answer after a refusal')

// Too few values is the shared reason, at the element that falls short.
e = refusal(() => m.signal_to_noise([[1, 2], [3]], 'NominalIsBest'))
ok(e.code === 'insufficient_data' && e.parameter === 'responses' && e.index === 1 && e.min === 2 && e.got === 1,
  `signal_to_noise short run: ${e.code} ${e.parameter}@${e.index} min=${e.min} got=${e.got}`)

// A fitted surface can be evaluated, and says what its coefficients are.
// The oracle is worked by hand: k=2, coded (0.5, -1),
// row [1, 0.5, -1, 0.25, 1, -0.5] . [10, 2, -3, -1.5, 0.5, 4] = 12.125
const rsmC = [10, 2, -3, -1.5, 0.5, 4]
const p1 = m.rsm_predict(rsmC, 2, Float64Array.from([0.5, -1]))
ok(p1.length === 1 && Math.abs(p1[0] - 12.125) < 1e-12, `rsm_predict: ${p1}`)
const pN = m.rsm_predict(rsmC, 2, Float64Array.from([0.5, -1, 0, 0, 1, 1]))
ok(pN.length === 3 && Math.abs(pN[1] - 10) < 1e-12, `rsm_predict batched: ${pN}`)
e = refusal(() => m.rsm_predict(rsmC, 2, Float64Array.from([0.5])))
ok(e.code === 'point_dimension_mismatch', `ragged coded: ${e.code}`)
e = refusal(() => m.rsm_predict([1, 2, 3], 2, Float64Array.from([0, 0])))
ok(e.code === 'coefficient_count_mismatch', `wrong coefficient count: ${e.code}`)

// A search that cannot report a 0 as an optimum.
const maxSpec = (lower, target) => [{ goal: 'Maximize', lower, target, upper: target, s1: 1, s2: 1 }]
let opt = m.optimize_desirability(maxSpec(0, 100), Float64Array.from([10, 50, 90]))
ok(opt.kind === 'best' && opt.index === 2 && Math.abs(opt.overall - 0.9) < 1e-12, `optimum: ${JSON.stringify(opt)}`)
opt = m.optimize_desirability(maxSpec(200, 300), Float64Array.from([10, 50, 90]))
ok(opt.kind === 'infeasible' && JSON.stringify(opt.unreachable) === '[0]', `infeasible: ${JSON.stringify(opt)}`)
e = refusal(() => m.optimize_desirability(maxSpec(0, 100), Float64Array.from([10, NaN])))
ok(e.code === 'value_not_finite' && e.index === 1, `a NaN is not a score: ${e.code}@${e.index}`)
// A number array is read as sent: null is refused at its index, not
// scored as 0 by the glue's typed-array copy.
e = refusal(() => m.optimize_desirability(maxSpec(0, 100), [10, null]))
ok(e.code === 'malformed_input' && e.index === 1 && e.parameter === 'candidates', `null candidate: ${e.code}@${e.index}`)

// An absent value is null, with its key present.
const f2 = m.full_factorial(2)
const sat = m.doe_anova(f2.data, [1, 2, 3, 5], f2.factor_names, ['A', 'B', 'A:B'])
ok('curvature' in sat && sat.curvature === null, 'curvature must be null')
ok(sat.effects.every((r) => r.f_statistic === null && r.p_value === null), 'a saturated F and p must be null, not NaN')
const fx = m.estimate_effects(f2.data, [1, 2, 3, 5], f2.factor_names, 1)
ok('lenth' in fx && fx.lenth === null, 'lenth must be null with two contrasts')

console.log(`binding behaviour gate: ok (${checks} checks)`)
