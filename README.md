# GRASS

**General Relativistic Axisymmetric Spacetime Solver**

[![CI](https://github.com/Hao-Jui/GRASS/actions/workflows/ci.yml/badge.svg?branch=dev)](https://github.com/Hao-Jui/GRASS/actions/workflows/ci.yml)
[![Fortran](https://img.shields.io/badge/Fortran-2003-734F96)](https://wg5-fortran.org/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

GRASS computes rotating neutron-star equilibria in either General Relativity (GR)
or Scalar–Tensor (ST) gravity. The code integrates the field equations on a
compactified meridional grid, reads tabulated equations of state (EOS),
and outputs diagnostic profiles that characterise the equilibrium model.

- Latest release: `66dfd4e`
- v2: `8a7b201`
- v1: `6139a9f`

See [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md) for the development workflow and
[docs/TESTING.md](docs/TESTING.md) for the full test/regression matrix.

---

## Quickstart

```bash
git clone https://github.com/Hao-Jui/GRASS.git
cd GRASS

# Build (CMake is the primary build system)
cmake --preset release
cmake --build build -j

# Run a single equilibrium model
./build/grass

# Run the full test suite
ctest --test-dir build --output-on-failure
```

Outputs land in `Cont/` (diagnostics) and `Res/res.rst` (restart binary).

---

## Features

- Uniformly and differentially rotating, axisymmetric neutron-star equilibria.
- Tabulated EOS support (`p`, `e`, `h`, `n0` columns), monotone PCHIP
  interpolation in log-space with phase-transition handling.
- Compactified meridional grid with spectral-type interpolation.
- Newton 1D / 2D shooting solvers to match target bulk properties (M, M_b, J, χ, axis ratio).
- Adaptive 4-stage relaxation: Picard → Chebyshev → Anderson → Aitken δ².
- Restart binary (`Res/res.rst`, magic `GRASSRST01`) for warm starts.
- M–R sequence builder, scalar-burning helper for ST gravity.

---

## Requirements

| Component | Version | Notes |
|---|---|---|
| Fortran compiler | gfortran ≥ 11 (or any F2003) | tested gfortran 13 / 14 on macOS + Linux |
| CMake | ≥ 3.21 | primary build; presets need 3.21 |
| BLAS / LAPACK | OpenBLAS, Accelerate, MKL | linked dynamically |
| FFTW3 | optional | only if FFT-based diagnostics are enabled |

The legacy GNU Make build is preserved in `Makefile`; CMake is recommended.

---

## Build

### CMake (recommended)

Three presets defined in `CMakePresets.json`:

```bash
cmake --preset release      # -O3, native arch, no debug
cmake --preset debug        # -O0 -g
cmake --preset sanitizer    # -fsanitize=address,undefined

cmake --build build -j
```

Custom flags:

```bash
cmake -B build -S . \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_FLAGS="-O3 -march=native"
cmake --build build -j
```

### Make (legacy)

```bash
mkdir -p Cont Res
make
```

Binary lands at `build/bin/a.out`.

### Reproducibility

Builds are reproducible given identical compiler, BLAS implementation, and
flags. CI runs in three configurations (release, debug, sanitizer); see
[`.github/workflows/ci.yml`](.github/workflows/ci.yml).

---

## Run

```bash
./build/grass
```

The binary reads its run configuration from compile-time defaults in
`src/para_panel.f90` (theory, EOS file, rotation law, grid resolution,
shooting target). Edit and rebuild to change the active configuration.

Diagnostics:

| File | Contents |
|---|---|
| `Cont/properties.dat` | bulk quantities (M, M_b, J, χ, Ω_c, Ω_e, R, …) |
| `Cont/Omega.dat` | equatorial profile: r, Ω, h, c_s², F_j |
| `Cont/velocity.dat` | orbital velocity diagnostics (GR) |
| `Cont/hamiltonian.dat` | Hamiltonian-constraint field |
| `Cont/moment_tail.dat` | far-field multipoles |
| `Res/res.rst` | restart binary |

---

## Configuration

All compile-time knobs live in `src/para_panel.f90`. Common edits:

```fortran
integer :: active_theory  = THEORY_GR        ! or THEORY_ST
integer :: run_task       = shoot            ! shoot | OneModel | MRbuild
character(len=128) :: eos_file = "MPA1"      ! file under eos/<name>.dat
integer :: SDIV = 2 * res + 1                ! radial grid; res ~ 200
integer :: MDIV = 41                         ! angular grid
```

See `src/para_panel.f90` for the full list (rotation law, shooting target,
relaxation thresholds, output verbosity).

Runtime overrides: if `./para.nml` exists in the working directory, its
`&para` namelist is read first in `initialize_theory`, so a parameter
search needs no rebuild. Recognised: `solver_type`, `A_diff`, `lambda1`,
`lambda2`, `uryu_p`, `uryu_q`, `Mb_goal`, `J_goal`, `M_goal`, `chi_goal`,
`omc_goal`, `run_task`, `run_mode`, `shooting`, `FIX1`, `FIX2`, `SDIV`,
`MDIV`, `eos_file`, `e_center_default`, `output`, `timing`, `r_ratio`
and `regrid_e_center`. Anything absent keeps its compiled-in default.
`regrid_e_center` (g/cm^3, > 0) replaces the restart's central energy
density in `MODE_REGRID`, e.g. to step a warm-started density chain.

Rotation law `solver_type = "uryu_inc"`: Uryu et al. 2017 (PRD 96,
103011) eq. (5) with the plus sign, `j = A^2 Omega (Omega/Omega_c - 1)^p`
for `Omega >= Omega_c`, so `Omega` rises monotonically outward. Set by
`inc_p` (real, `> 0`) and `inc_lambda = Omega_e/Omega_c` (`> 1`), both in
`para.nml`; `A` follows from the equator condition. The first integral
adds `-A^2 Omega_c^2 [x^(p+1)/(p+1) + x^(p+2)/(p+2)]`, `x = Omega/Omega_c - 1`,
to the log-enthalpy.

```fortran
&para
  solver_type = "uniform"
  run_task = 3
  run_mode = 1
  r_ratio = 0.85
  regrid_e_center = 7.157e14
/
```

---

## Repository Layout

```
src/
  main.f90                     program entry point
  para_panel.f90               global parameters, grid arrays, field allocations
  core/                        EOS, grid, shoot, analysis, constraint, restart
    eos_mod.f90                PCHIP-interpolated tabulated EOS
    grid_mod.f90               compactified meridional grid
    shoot_mod.f90              Newton driver
    shoot_solver_*_mod.f90     1D / 2D / r_ratio shooters
    analysis_mod.f90           ADM mass, M_b, J, I, Love number
    constraint_mod.f90         Hamiltonian constraint
    starting_model_mod.f90     initial guess + single-model driver
    MRcurve_mod.f90            M–R sequence builder
    regrid_mod.f90             restart binary I/O
    sphere_mod.f90             spherical pre-iteration
    miscellaneous_mod.f90      converged-block printer
    scalar_burning_mod.f90     ST scalar-mass annealer
  theory/                      rotation solver + relaxation
    rotation_solver_mod.f90    outer loop on equatorial radius
    spin_integration_mod.f90   Green's-function targets, multipole projection
    spin_relaxation_mod.f90    4-stage relaxation orchestrator
    relaxation_mod.f90         Picard / Chebyshev / Anderson / Aitken kernels
    spin_workspace_mod.f90     geometry + workspace caching
    spin_updates_mod.f90       metric / scalar field updates
    spin_derivatives_mod.f90   grid derivatives (DGEMM-based when collocated)
    rotational_law_mod.f90     uniform / const-J / Uryu / increasing-Uryu rotation laws
  tool/                        numerical utilities
    toolkit_mod.f90            interp, derivatives, integral accumulation
    ad_mod.f90                 forward-mode AD (dual numbers)
    cheb_mod.f90               Chebyshev / Legendre bases
    spectral_hub_mod.f90       spectral basis dispatcher
    brent_mod.f90              1D root finding
    nag_compat_mod.f90         adaptive quadrature (NAG-D01 drop-in)
    spline_mod.f90             1D spline (utility)
    lapack_interfaces_mod.f90  explicit BLAS / LAPACK interfaces

eos/                           EOS table files (user-supplied)
tests/                         unit + integration + regression + harness scripts
docs/                          supplementary notes
.omc/research/                 perf audit reports + microbenchmarks
```

---

## EOS Format

Plain-text table:

```
num_tab
e(1)   p(1)   h(1)   n0(1)
…
e(N)   p(N)   h(N)   n0(N)
```

Columns are physical quantities (CGS energy density and pressure;
dimensionless relativistic enthalpy `h = 1 + (e + p) / (ρ_0 c²)`; baryon
number density). The loader (`eos_mod::loadEos`) converts to log-space
and detects phase-transition rows automatically. PCHIP slopes are
precomputed once per EOS load for every direction the public API uses.

---

## Solver Data Flow

1. `initialize_theory()` — set GR / ST and compile-time constants.
2. `loadEos` — read EOS, build PCHIP slopes.
3. `make_grid` + `GridTrig` — compactified `s ∈ [0, SMAX]`, `μ ∈ [0, 1]`.
4. Task dispatch:
   - `OneModel` → single equilibrium via `initialize_starting_model`.
   - `shoot` → 1D / 2D Newton on bulk targets.
   - `MRbuild` → mass-radius sequence sweep.
5. `rotation_solver` outer loop on equatorial radius `r_e`:
   - Compute Green's-function targets (`spin_integration_mod`).
   - Apply 4-stage relaxation to drive metric + scalar to targets.
   - Convergence: `|r_e_new / r_e_old − 1| < 1e-7` (configurable).
6. `solution_properties` — global diagnostics.
7. `hamiltonian` — constraint check (Ham L2 ≲ 1e-10 at convergence).

---

## 4-Stage Relaxation Strategy

Implemented in `src/theory/relaxation_mod.f90` and orchestrated by
`spin_relaxation_mod.f90`.

| Stage | Regime | Method | Threshold |
|---|---|---|---|
| 1. Picard damping | `dif ≥ 0.1` | `x ← (1 − w) · x + w · target`, `w = 0.7` | `PICARD_THRESH` |
| 2. Chebyshev | `0.01 ≤ dif < 0.1` | 3-term recurrence on EMA-smoothed spectral radius | `CHEB_THRESH` |
| 3. Anderson | `dif < 0.01` | rank-3 Anderson, per-field | `N_ANDERSON` |
| 4. Aitken δ² | exponential tail | quadratic extrapolation when `ρ_obs` is locked | `N_MONOTONE` |

Typical run: 40–90 iterations total (5–7 Picard, 2–3 Chebyshev, 3–5
Anderson, 30–75 Aitken at ρ ≈ 0.85).

| Parameter | Default | Effect |
|---|---|---|
| `PICARD_THRESH` | 1.0e-1 | Picard → Chebyshev cutover |
| `CHEB_THRESH` | 1.0e-3 | Chebyshev → Anderson cutover |
| `w_picard` | 0.7 | Picard damping; lower for tighter coupling |
| `N_MONOTONE` | 5 | Aitken activation gate |
| `MAX_AITKEN_JUMP` | 10 | per-element Aitken safety bound |

---

## Validation

- Hamiltonian constraint norm `O(1e-10)` at convergence.
- M–R sequences cross-checked against published references (MPA1, APR4, …).
- `tests/integration/test_*` exercise GR-uniform, GR-constJ, GR-Uryu, and
  ST-uniform configurations with `1e-4` reference-tolerance gates.
- Restart round-trip (`test_restart`) validates writer / reader unit consistency.

See [docs/TESTING.md](docs/TESTING.md) for the full harness inventory.

---

## Citing GRASS

If you use GRASS in published work, please cite:

```bibtex
@software{grass,
  author = {Kuan, Hao-Jui},
  title  = {GRASS: General Relativistic Axisymmetric Spacetime Solver},
  url    = {https://github.com/Hao-Jui/GRASS},
  year   = {2026}
}
```

---

## License

MIT. See [LICENSE](LICENSE).

---

## Support

Open an issue with:

- compiler version + OS + build flags,
- a `Cont/properties.dat` excerpt or console log,
- the EOS file used (or its first 5 rows).

Pull requests welcome — see [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md).
