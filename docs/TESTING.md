# Testing & Quality Harnesses

GRASS is a single-binary Fortran research solver with no network surface,
no database, no UI, and no deployment pipeline. The harness inventory
below reflects that scope: web-service-flavoured categories (E2E, canary,
migration, SAST) are explicitly **not applicable by design** and are
listed at the bottom so future readers know they are omitted on
purpose, not by oversight.

## Quick reference

```bash
# Full suite (unit + integration + regression)
ctest --test-dir build --output-on-failure

# Single test
ctest --test-dir build -R test_eos --output-on-failure
```

CI runs the full suite on every push to `dev` / `main` in three
configurations (release, debug, sanitizer); see
[`.github/workflows/ci.yml`](../.github/workflows/ci.yml).

---

## Test harnesses

### Unit (`tests/unit/`)

Isolated correctness checks for individual modules. Each test is a
standalone Fortran program that calls a fixed API and reports
pass/fail per assertion via `tests/unit/test_utils.f90`.

| Test | Target | Notable assertions |
|---|---|---|
| `test_ad` | `ad_mod` (dual-number AD) | derivative chain rule, log/exp consistency |
| `test_brent` | `brent_mod` (1D root find) | bracketing, convergence, edge cases |
| `test_eos` | `eos_mod` (PCHIP EOS) | round-trips at 1e-10, `pe_at_h` consistency at 1e-12, dual derivative vs FD at 1e-4, table-edge stability |
| `test_spectral` | `spectral_hub_mod` | basis orthogonality, derivative recovery |
| `test_spline` | `spline_mod` | knot-by-knot interpolation, end-condition handling |
| `test_uryu_inc` | `rotation_law_mod` (`uryu_inc`) | first-integral term G vs Gauss-Legendre quadrature of the law at 1e-12 (p = 0.5, 1, 1.3, 2); equator identity G(lambda-1) = F_e Omega_e c(lambda,p); bulk Omega root vs an analytic flat-metric case at 1e-10, incl. the first-root selection and the axis |
| `test_uryu_fmax` | `rotation_law_mod` (`uryu`, analytic Fmax) | y* = Fmax/F_e solves h(y) = 0 to 1e-13 for five (lambda1, lambda2); with A, B from the value conditions, Omega(Fmax) = lambda1, Omega(F_e) = lambda2 (1e-12), dOmega/dF(Fmax) = 0 and Omega <= lambda1 Omega_c on a dense F grid, for three F_e (y* independent of F_e); `uryu_omega` vs a brute-force scan in a flat metric (root nearest the guess), and the axis |

Run individually: `./build/tests/test_<name>`.

### Integration (`tests/integration/`)

End-to-end equilibrium solves driven from program entry, with reference
output gated by `tests/check_regression.sh` at **1e-4 relative
tolerance** against `tests/reference/test_<name>.ref`.

| Test | Configuration |
|---|---|
| `test_gr_uniform` | GR + uniform rotation, MPA1, SDIV=401, MDIV=41, r_ratio = 0.7 |
| `test_gr_constj` | GR + constant-J differential rotation |
| `test_gr_uryu` | GR + Uryu rotation law |
| `test_st_uniform` | Scalar–tensor + uniform rotation |
| `test_st_uniform_r07` | ST + uniform rotation, r_ratio = 0.7 (locked-axis variant) |
| `test_restart` | round-trip of `Res/res.rst` (writer ↔ reader unit consistency) |

Each test sets its own `eos_file`, `active_theory`, `solver_type`,
`SDIV`, `MDIV`, etc. inside the test program — the global compile-time
defaults in `src/para_panel.f90` are overridden so a single binary can
exercise multiple configurations.

### Regression gate (`tests/check_regression.sh`)

Reads each numeric line from a reference `.ref` file, finds the
corresponding label in `Cont/properties.dat`, and fails if relative
error exceeds **1e-4** (4 significant digits). Catches reintroduced
solver bugs and silent EOS-table changes. The gate is invoked
automatically by every integration test via CTest.

### Runtime config overrides (`GRASS_TIMING`, `GRASS_E_CENTER`)

The integration tests pin two knobs that double as developer
sweep / profiling state:

- `e_center` default-branch literal (`src/core/starting_model_mod.f90`),
  pulled from `e_center_default` in `src/para_panel.f90`.
- `timing` flag (`src/para_panel.f90`).

Both are now read at startup from environment variables in
`apply_runtime_overrides` (called from `initialize_theory`):

| Env var | Effect | Accepted values |
|---|---|---|
| `GRASS_TIMING` | overrides `timing` | `0/1`, `true/false`, `yes/no`, `on/off`, `t/f` |
| `GRASS_E_CENTER` | overrides `e_center_default` (cgs g/cm³) | float literal, e.g. `0.8e15` |

`tests/CMakeLists.txt` sets `ENVIRONMENT "GRASS_TIMING=0;GRASS_E_CENTER=0.8e15"`
on every integration test, so `ctest` is hermetic regardless of the
developer's local compile-time defaults. No source patching, no rebuild
between dev work and tests.

Manual examples:

```bash
GRASS_TIMING=1 ./build/grass                    # profiling run
GRASS_E_CENTER=1.0e15 ./build/grass             # different central density
ctest --test-dir build -R test_gr_uryu          # baseline (env set by CTest)
```

---

## Quality harnesses

### Lint / static analysis

`gfortran -Wall -Wextra -Wno-maybe-uninitialized` is enabled by the
debug and sanitizer presets. There is no separate lint pass; the
compiler's diagnostics are the gate.

### Type checking

Fortran is statically and strongly typed at compile time; no separate
type-check step is needed.

### Sanitizers

The `sanitizer` CMake preset enables AddressSanitizer + UndefinedBehavior
sanitizer. CI runs the full unit + integration suite under sanitizer
on every push (with `test_brent` excluded due to a gfortran-runtime
ASan crash — see commit `d5e8ac3`).

### Coverage

Not currently enforced. `gfortran -fprofile-arcs -ftest-coverage` plus
`gcov` produces line coverage if needed; treat it as ad-hoc.

---

## Performance harnesses

There is no recurring perf gate in CI today. One-shot microbenchmarks
have been written to support specific perf-audit questions and live
under `.omc/research/` as research artifacts (not regression gates):

- `.omc/research/bench_dgemm_stack.f90` — wall-clock evidence for the
  proposed batched-derivative-DGEMM speedup. Used to verify the
  audit's speculation; the optimization was *not* adopted, so this
  file remains as methodology reference only and is not run by any
  automation.

If a recurring perf gate is needed in the future, a microbench should
be promoted to `tests/perf/` and wired into CTest with a tolerance-
gated baseline (same pattern as the integration regression).

### Production timing breakdown

`src/para_panel.f90 :: timing = .true.` enables in-binary CPU-time
buckets that print every 5 outer iterations and then early-exit:

```
precompute avg over 5:
  bessel              5.40e-06
  first_deriv_s       4.96e-04
  first_deriv_m       6.19e-03   ← summed across BLAS threads
  …
get_all_targets avg over 5:
  precompute        1.06e-02
  build_source_terms 3.56e-03
  …
```

Note: cpu_time is the *sum* across all BLAS-internal threads. With
OpenBLAS at default thread count, the CPU/wall ratio for these buckets
is roughly 3–6×; divide by that ratio to recover wall-clock impact.

Integration tests force `GRASS_TIMING=0` so the solver runs to
convergence regardless of the source default; for profiling, set
`GRASS_TIMING=1` in the shell before invoking the binary.

### Audit reports

- `.omc/research/theory_perf_audit.md` — ranked-recommendation perf audit.
- `.omc/research/theory_perf_audit_evidence.md` — wall-clock-vs-cpu_time
  evidence pass on top recommendations.

---

## Categories not applicable to GRASS

These appear in production-codebase checklists but have no analogue here.
Listed so the omission is explicit:

| Category | Why not applicable |
|---|---|
| End-to-end (Playwright / Cypress / Selenium) | no UI |
| API contract tests | no HTTP / RPC surface |
| Database integration | no persistent store |
| Dependency audit (npm audit, pip-audit, …) | system-library deps only (BLAS, LAPACK, FFTW); no package manifest |
| SAST / semgrep / CodeQL | no untrusted input, no network surface, no auth code |
| Load / stress tests | single-threaded by design (no parallelization, per CLAUDE.md) |
| Smoke test post-deploy | no deploy pipeline |
| Canary gate | no rollout |
| Migration dry-run | no schema |

---

## Test minimum-viable set (what we run on every CI push)

- Unit (5 tests) ✅
- Integration (6 tests) ✅
- Regression gate (1e-4) ✅
- Sanitizer build (ASan + UBSan) ✅
- Compiler-warning gate (`-Wall -Wextra`) ✅

This satisfies the "minimum viable production set" for a research-HPC
codebase: correctness + boundary contracts + reintroduction protection.
