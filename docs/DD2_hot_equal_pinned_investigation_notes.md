# DD2_hot_equal_pinned in GRASS: TOV star, Uryu-law bugs, and the unresolved rotating-Ham-L2 mystery

Investigation notes from building a TOV star and then a differentially-rotating
star in GRASS from `DD2_hot_equal_pinned` (PyCompOSE's hot-merger-remnant
extension of the DD2 EOS). Covers everything found and fixed, in
chronological order, plus what remains open.

## ⚠️ Critical methodology finding (read this first)

**A single big jump in `r_ratio` (`MODE_REGRID` straight from the `r_ratio=1.0`
TOV restart to the target, e.g. `0.995`, in one solve) can under-converge by
4+ orders of magnitude relative to a gradual spin-up (`1.0→0.9995→0.999→
…→0.995`, each step re-relaxed before the next).** Demonstrated on MPA1:
target `r_ratio=0.995`, same EOS, same shape parameters, same everything —
single jump gives `Ham L2=153.12` (confirmed **exactly** reproducible 4/4
runs from an identical frozen restart — this is not run-to-run noise, it is
a genuine, deterministic, differently-converged fixed point of the
relaxation); gradual spin-up (6 steps) gives `Ham L2=4.91e-3` at the *same*
target `r_ratio`. A ~31,000× difference from convergence path alone, nothing
physical.

**This retroactively confounds the MPA1-decimation ("EOS knot-density")
result reported earlier in this file.** Re-run with gradual spin-up:
`MPA1_coarse` (the decimated table) converges to `Ham L2=4.54e-3` —
statistically indistinguishable from full-resolution MPA1's own `4.91e-3`.
**The "113,000× degradation from table sparsity" claim below does not survive
proper convergence and should be treated as retracted** — it compared a
gradually-converged historical baseline against a single-jump result for the
decimated table, an apples-to-oranges comparison, not a genuine EOS effect.

**DD2 is the one thing that survives this correction.** Re-run with the
*same* gradual 6-step spin-up (using the bulk-fixed, seam-smooth,
standard-surface, 2441-row table): `Ham L2` moves smoothly through
`1.20→1.16→1.14→1.13→1.13`, landing at `1.1256` — **matching its own
single-jump result to 4 significant figures.** Unlike MPA1, DD2's single-jump
and gradual-spin-up paths converge to the *same* fixed point. This means (a)
DD2's poor `Ham L2` is real, not a convergence-path artifact, and (b) every
single-jump DD2 measurement earlier in this file (all the hypothesis tests)
was very likely already representative of DD2's true converged behavior,
unlike the MPA1 comparisons.

**Not yet re-verified under gradual spin-up** (flagging rather than
re-running, given time already spent establishing the above): the "MPA1
elevated-surface → `Ham L2=153`" surface-level result. It was almost
certainly measured via single-jump (matching every other test's
methodology, this session and prior) and needs the same gradual-spin-up
check before it can be trusted. Until re-verified, treat the "confirmed
surface-level mechanism" claim below as **unconfirmed**, not disproven.

**Practical rule going forward**: any new EOS/config A/B test must use a
matched convergence path for both arms (ideally gradual spin-up for both,
or at minimum single-jump for both) — never compare a historical
gradually-converged number against a fresh single-jump one.

---

## Summary (superseded in part — see the box above)

- **Phase 1 (TOV/non-rotating star): fully working**, `Ham L2≈5e-3`.
- **Two real, confirmed GRASS bugs found and fixed** in the Uryu
  differential-rotation law's implementation (independent of any EOS
  choice — both also affect MPA1 in the same regime).
- **Surface-level mechanism (MPA1 `1e-3→153`)**: status downgraded to
  **unconfirmed** — likely measured via single-jump regrid against a
  gradually-converged baseline; needs re-verification (see box above).
- **First spatial localization of the residual** (via `Cont/hamiltonian.dat`,
  enabled this round): DD2's poor `Ham L2` is grid-scale noise concentrated
  at **30–45% of the stellar radius** (`s_gp≈0.30–0.45`), not at the surface.
  Rotation amplifies exactly these shells by 175–450× relative to the TOV
  baseline; the near-surface shell's *share* of the norm actually *shrinks*
  under rotation. This localization used DD2's single-jump result, which is
  now confirmed representative of its true converged state (see box above)
  — this finding stands.
- **RETRACTED**: the "EOS table knot-sparsity" mechanism (MPA1 decimation
  `1e-3→113`, and the subsequent DD2 "bulk fix" attempt). Does not survive
  gradual-spin-up re-verification (`MPA1_coarse` gradual: `4.54e-3`, no
  degradation). See box above.
- **The one finding that has survived every test in this investigation,
  including this session's most rigorous check**: `DD2_hot_equal_pinned`
  gives `Ham L2≈1.05-1.2` for *any* differential rotation, *any* convergence
  path (single-jump or gradual spin-up), and now also *any* EOS-table
  smoothness/knot-density (13 knots or 258, linear or cubic lookup) tried
  against it. Root cause still not found. The `alpha`/`ricci` term-
  decomposition finding below (a real second-derivative-operator noise-
  amplification mechanism, specific to how `alpha` is constructed) is the
  strongest still-standing lead, and — unlike the knot-density finding — was
  measured on a run now confirmed to be properly converged, so it is not
  subject to the same retraction.

---

## Phase 1: TOV (non-rotating) star — SOLVED

**Problem**: `sphere_mod::TOV`'s radial integrator has no iteration cap, and
`DD2_hot_equal_pinned` (built from `PyCompOSE`'s `fit_appendixA.py`, a hot
merger-remnant profile) doesn't extend anywhere near GRASS's default vacuum
surface convention (`E_SURFACE=7.8 g/cm³`, `P_SURFACE=1.01e8 dyn/cm²` —
literally iron density). The table's own low-density floor sat ~13 orders
of magnitude higher in pressure at the time, so the integrator ran
effectively forever trying to reach an unreachable target.

**Fix (at the time)**: re-pointed `E_SURFACE`/`P_SURFACE` in
`src/para_panel.f90` at the table's own low-density floor instead of the
standard convention. (Superseded later — see "Absolute surface level"
below — but this got Phase 1 unblocked.)

**Also fixed**: `r_ratio` was being unconditionally clobbered to `0.7` by
`single_model()` in `src/core/starting_model_mod.f90` regardless of what was
configured, so `OneModel` always built a rotator, never a true TOV star.
Removed the hard-coded overwrite (and an earlier `solver_type`-based guess
in the `case default` branch of `initialize_starting_model`) so `r_ratio` is
honored as configured — `1.0` for Phase 1, a genuine rotator value for
Phase 2.

**Environment note**: the build links MKL's `libmkl_rt.so.2`, which lazily
dispatches to `libmkl_intel_thread.so.2` on first real MKL call — that
threading layer needs Intel's OpenMP runtime, not present in this shell.
Run with `MKL_THREADING_LAYER=GNU` to route through `libgomp` instead (no
source change needed).

**Result**: converges cleanly, `Ham L2≈5e-3` — actually *better* than
GRASS's own MPA1 reference test achieves for a rotating model.

---

## Phase 2: Uryu-law differential rotation — two real bugs found and fixed

### Bug 1: `brent_core: failed to bracket root` (NaN in `AA_h`)

**Symptom**: crashed when spinning up from the TOV restart to
`r_ratio=0.995` with the Uryu rotation law.

**Root cause**: `AA_h`'s defining exponent, `1/(uryu_p+uryu_q) = 1/4`, is
fractional. Its bracketed base — `(F_e·Fmax_h)·tmp1/tmp2` — goes negative
for the *entire band* `0 ≤ F_e < (lambda1/lambda2)^(1/uryu_q)·Fmax_h` (not
just `F_e<0`). Fortran's `real**noninteger` on a negative base silently
returns `NaN`. `find_omega_e`'s geometric bracket search starts right in
this near-spherical, `F_e≈0` regime (exactly where a gradual spin-up seeds
it from a TOV restart) and can't escape the NaN wall.

**Fix**: `src/theory/rotational_law_mod.f90` — added a shared `safe_F_e`
function that clamps `F_e` comfortably above the whole unsafe band before
calling `AA_h`/`BB_h`. Applied at **two** call sites:
- `diff_rotation_uryu` (the original crash site).
- `cache_uryu_ab` (a second, initially-missed call site using the raw,
  unclamped `F_equator_h` global — this was poisoning `cached_aa`/`cached_bb`
  with `NaN` even after the first fix, breaking the downstream peak-search
  `zbrent_rot`/`rotation_law_uryu` calls with a different crash signature,
  `Cont/rotation_law.dat`).

**Verified universal, not EOS-specific**: MPA1 hits the *identical* crash
at the same `r_ratio`/default shape parameters — this is a general Uryu-law
implementation bug, unrelated to `DD2_hot_equal_pinned`.

### Bug 2 (not a crash, a stop): `negative F_equator_h`

**Symptom**: after fixing Bug 1, spin-up to `r_ratio=0.995` with the
*default* Uryu shape (`lambda1=1.5, lambda2=0.3`) ran cleanly through the
whole solve, then hit `STOP negative F_equator_h; L120 in uryu` (a clean
stop, not a crash) — the found `Omega_e` came out less than the local
frame-dragging value.

**Root cause**: confirmed via diagnostic instrumentation — the found root
tracks smoothly for hundreds of outer `Fmax_h`-shooting iterations, then at
a specific `Fmax_h` threshold the root **jumps discontinuously** to a
spurious branch (`Omega_e` collapsing toward `0`). Two fix attempts (a
physically-motivated `Fmax_h` reseed; a bounded-step continuation method on
the shooting loop) both failed to help — the bounded-step attempt failed
*faster*, which rules out a numerical-overshoot explanation and confirms
this is a genuine branch discontinuity in the equation for extreme
differential rotation this close to spherical.

**Resolution**: not a bug fix — milder Uryu shape parameters
(`lambda1=1.1, lambda2=0.9` instead of `1.5/0.3`) avoid the discontinuity
entirely and converge cleanly. **Confirmed universal**: MPA1 hits the same
stop at the same `r_ratio` with the default shape too — this is inherent to
extreme differential rotation near `r_ratio=1` for *any* EOS, not specific
to our table.

---

## The `Ham L2` investigation: a large real effect found, but DD2's own problem remains unexplained

With the crash fixed and milder shape parameters avoiding the stop, DD2
rotates without crashing — but `Ham L2≈1.05-1.3`, 1-3 orders of magnitude
worse than MPA1 achieves under the identical solver/shape/rotation rate
(MPA1: `≈1e-3` at `r_ratio=0.995`, `≈0.09` at `r_ratio=0.7`).

### Hypotheses tested and eliminated (in order)

1. **Vacuum-surface placement within the table's own narrow range** — moved
   the cutoff from `e=2.4e3` to `e=1.9e3 g/cm³`: no change.
2. **Table resolution** — resampled 326→6000 rows (PCHIP, same curve): no
   change.
3. **Warm-vs-cold-crust construction** — reworked `fit_appendixA.py`'s
   low-density leg to prefer staying at the table's coldest tabulated `T`
   (using DD2's own physically self-consistent NSE crust) instead of
   warming via entropy-pinning: fixed a real physics issue (the warmed
   version was unphysically soft, `Gamma1` ~0.5-0.85× MPA1's), but `Ham L2`
   unchanged.
4. **A sharp local `Gamma1` discontinuity** at the `i_min` boundary
   (`T` jumping `2.32→5.67 MeV` in one table row) — fixed with a smooth
   log-linear `T` glide window; `Gamma1` spike reduced `9.37→2.15`, but
   `Ham L2` **unchanged**, and — checked directly in the 2D constraint
   field — the exact same dominant grid shell contributed the exact same
   magnitude before and after.
5. **First-law/Euler-relation thermodynamic consistency** — checked
   directly against the table's own `Q1,Q2,Q3,Q4,Q7` at and around the
   implicated density: smooth, no violation, actually *better* than a
   random baseline point elsewhere in the table.
6. **GRASS's own grid resolution** — doubled `SDIV` (`801→1601`): `Ham L2`
   got **worse** (`1.067→1.288`), ruling out simple under-resolution (a
   discretization error should improve with resolution, not degrade).
7. **Rotation-law-specific bug** — tried `const_j` (a completely different
   formulation, no shared code with Uryu): identical `Ham L2≈1.06`,
   confirming the issue has nothing to do with the Uryu implementation.

### The hypothesis that WAS confirmed — on MPA1, not on DD2 (⚠️ status downgraded, see top-of-file box)

8. **Absolute vacuum-surface pressure level** (not just "which row within
   DD2's own table"): raised **MPA1's own** surface cutoff to match DD2's
   absolute level (`e≈2000 g/cm³`, an in-table, non-extrapolated MPA1 row).
   Result: MPA1's `Ham L2` exploded from `~1e-3` to **`153`** — a
   ~150,000× degradation. This is a large, real, confirmed effect: an
   elevated vacuum-surface cutoff, independent of which EOS underlies it,
   badly destabilizes a rotating solve.

   **⚠️ Update, this session**: the *same* number, `Ham L2=153`, is exactly
   what a plain, non-elevated-surface MPA1 gives when reached via a
   single-jump `r_ratio` regrid instead of gradual spin-up (see the
   top-of-file box) — and MPA1's `~1e-3` baseline was almost certainly
   itself established via gradual spin-up. This raises real doubt about
   whether finding 8 measured the surface-level effect at all, versus simply
   reproducing the single-jump convergence artifact under a different name.
   Not retracted outright (unlike the knot-density finding below, which
   *was* directly re-checked and failed) because it hasn't been directly
   re-run with matched convergence paths yet — but it should be treated as
   **unconfirmed** until that re-check happens.

### Fixing DD2's own surface level — still didn't help

Given finding 8, extended `DD2_hot_equal_pinned` down to the standard
near-vacuum convention using a cold degenerate polytropic append in
`fit_appendixA.py` (matched continuously at the table's own floor,
`n0=1e-12 fm⁻³`):

- First attempt used `Gamma=1.584` (matched to *MPA1's* edge slope) — hit a
  catastrophic-cancellation bug (absolute `MeV/fm³` energy-density
  arithmetic subtracting two near-equal `O(mn·n0)` terms), then, once fixed
  numerically, introduced a **new** `Gamma1` discontinuity at the seam
  (`1.58→0.95`) because DD2's own edge slope (`~0.95`, confirming the
  "warm-crust softness" finding) doesn't match MPA1's.
- Fixed by matching `Gamma` to **DD2's own edge slope** instead (via the
  same 3-point one-sided log-log derivative as PyCompOSE's
  `Table.get_polytrope(nb_idx=0)`, `Gamma≈0.9495`) — confirmed smooth
  `Gamma1` continuity at the seam (`0.946→0.951`, no jump).
- Result with the surface now at `e≈8.8 g/cm³` (matching the standard
  convention almost exactly) **and** a smooth seam: **`Ham L2` still
  `≈1.05`, completely unchanged.**

### The puzzle this left (surface level)

MPA1 is *highly* sensitive to surface level across the tested range
(`8.9→2000 g/cm³`: `1e-3→153`). DD2 shows **zero** sensitivity across the
*same* range (`8.8→2400 g/cm³`: always `~1.05-1.08`). If surface level were
DD2's dominant problem, some sensitivity should show up somewhere in that
range — it doesn't. This means DD2 has its own, apparently much larger,
independent problem that the surface-level effect (real, proven on MPA1)
doesn't explain and fixing the surface can't touch.

---

## Localizing the residual: `Cont/hamiltonian.dat` field analysis

Enabled the previously dead-code-gated 2D Hamiltonian field dump
(`constraint_mod.f90`, `if (.false.)` → `if (.true.)`), added `r_phys` and
`enthalpy` columns, and later `ww`/`velocity_sq` — and read the field
directly rather than only the scalar `Ham L2` norm.

### Where the residual actually lives

Shell-by-shell (proxy for the norm's own `s_gp²/(1-s_gp)⁴`-weighted integrand,
validated against GRASS's own reported `Ham L2` — proxy ratio `205.5` vs
actual `210.9` for TOV-vs-rotating): comparing a genuine TOV run
(`Ham L2=4.977e-3`) against the seam-smoothed rotating run
(`Ham L2=1.0498`) on the *same* table:

| `s_gp` window | TOV contrib | rotating contrib | ratio |
|---|---|---|---|
| 0.00–0.20 | tiny | tiny | 1–2500× (both negligible) |
| 0.20–0.30 | 0.65 | 233 | ~385× |
| 0.30–0.40 | 91 | 20,627 | ~230× |
| 0.40–0.45 | 108 | 44,255 | ~408× |
| 0.45–0.50 (surface) | 69 | 165 | **2.4×** |

Rotation amplifies the 30–45% shells by two-to-three orders of magnitude
more than it amplifies the near-surface shell, and the surface's *fractional
share* of the total actually **shrinks** under rotation (`0.217→0.0025`).
This directly explains, retroactively, why all nine surface/low-density
hypotheses above came back negative or ineffective: none of them touch the
region that actually dominates.

On the enthalpy field itself the metric solve is clean: 4th-difference
smoothness relative to the field's own range is `8.7e-9` (machine noise). The
`log e(h)` reconstructed from the table via PCHIP at those same grid points
was `1.5e-2` before the fix below — six orders of magnitude rougher than its
own input, with the roughness spikes landing exactly on the `ham` spikes.

This also retroactively explains hypothesis 6's anomaly (doubling `SDIV`
made `Ham L2` *worse*, `1.067→1.288`): halving `DS` amplifies a fixed-
amplitude curvature kink by `~1/DS²` in the Laplacian operator — the
opposite of what a genuine discretization error would do.

### Mechanism found: `eval_given_rtx`'s default is piecewise-linear

`nb_orig` (326 points, `fit_appendixA.py`) is CompOSE's own native table
grid. `eos.eval_given_rtx(..., method="linear")` (the default;
`compose/eos.py:284`, `scipy.interpolate.RegularGridInterpolator`) is
piecewise-**trilinear** — between two adjacent `nb` knots, `Q1(nb)` (hence
pressure) is exactly linear, so `d(Q1)/d(log nb)` (`Gamma1`) jumps
discontinuously at every knot. DD2's core has only **13 such knots** across
the whole stellar density range (`log_h≈0.09–0.32`); MPA1 (5931 rows, a
smooth analytic-fit EOS) has **116–182** there and shows no such noise.

### Confirming it: two independent A/B tests (⚠️ RETRACTED — see below)

1. **PCHIP-resampling the output curve does nothing** (settles the earlier
   hypothesis-2 result, which had never been checked against the field
   directly): resampling `DD2_hot_equal_pinned.dat` 446→6000 rows (same
   curve, denser sampling, `resample_grass_table.py`) changes the rotating
   field's shell-by-shell noise by a ratio of `1.00` everywhere — `Ham L2`
   `1.0498→1.0536` (no improvement, marginally worse). Resampling a curve
   through the same knots cannot remove curvature discontinuities baked
   into those knots; more points on a piecewise-linear function are still
   piecewise-linear. (This piece is unaffected by the retraction below —
   both numbers were measured the same way, single-jump, so the comparison
   is internally valid; it's the conclusion drawn from item 2 that failed.)
2. **Decimating MPA1 reproduces the effect from scratch** — the decisive
   test, independent of anything DD2-specific: subsampled `MPA1.dat` to
   DD2-like knot sparsity (458 rows, every 13th row, 8 core knots vs DD2's
   13 — real data removed, not resampled) with the **standard surface level
   held fixed** throughout (`E_SURFACE=7.8`, `P_SURFACE=1.01e8`, unchanged
   from MPA1's normal working configuration). TOV converges cleanly
   (`Ham L2=7.4e-6`); rotating (single-jump): `Ham L2` `1e-3 → 113`, a
   ~113,000× degradation from table sparsity **alone**. The noise localizes
   to the same kind of mid-radius shells (`s_gp≈0.15–0.35` and `0.40–0.50`,
   ~97% of the norm) as DD2's own case.

   **⚠️ RETRACTED**: re-run with a matched gradual 6-step spin-up (see the
   top-of-file box) — `MPA1_coarse` converges to `Ham L2=4.54e-3`,
   indistinguishable from full-resolution MPA1's own gradual result
   (`4.91e-3`). The `1e-3→113` comparison silently mixed a gradually-
   converged baseline (`1e-3`) against a single-jump result for the
   decimated case (`113`) — the "~113,000× degradation" was almost entirely
   convergence-path artifact, not an EOS-knot-density effect. Decimating
   MPA1's table has, at most, a small, unremarkable effect on `Ham L2` once
   both cases are converged the same way.

This mechanism is **not** established — the MPA1 sparsifying-direction test
that was supposed to confirm it does not survive proper reconvergence.

### The fix attempted — and its negative result on DD2 (this part still stands)

Modified `fit_appendixA.py`: for `i≥i_min` (the bulk/core region, where
`T_full`/`Yl_full` are the smooth Appendix-A fit curves `T_model()`/
`Yl_model()`, not the discrete per-row entropy-pin/T-glide optimization that
only runs in the crust below `i_min`) — refined each inter-knot interval
20×, then re-looked-up `Q1`/`Q2`/`Q7` via `eval_given_rtx(..., method="cubic")`
(tricubic reconstruction of the native CompOSE table, using curvature from
all surrounding grid points, not just the 2 nearest) instead of `"linear"`.
Crust logic (`i<i_min`) and the polytropic low-density extension are
untouched. Result: 13→**258** core knots (exceeding MPA1's 116), confirmed
monotonic end-to-end, confirmed smooth `Gamma1` (no seam jump), confirmed
`e(h)` mapping unchanged from the old table at matched `log_h` (sanity
check that the fix didn't silently distort the physics).

- TOV: `Ham L2=4.957e-3` — unchanged, as expected (this only touches the
  interior; TOV was never the problem).
- Rotating (single-jump): `Ham L2=1.1256` — **no improvement**, marginally
  worse than the pre-fix `1.0498`.
- **Re-verified with gradual 6-step spin-up** (this check is *not* subject
  to the retraction above, since it was run specifically to test whether
  DD2's own numbers hold up — see the top-of-file box): `Ham L2` moves
  smoothly `1.20→1.16→1.14→1.13→1.13`, landing at `1.1256` — matching the
  single-jump result to 4 significant figures. DD2's poor `Ham L2` is
  confirmed real and convergence-path-independent, using the exact same
  table this section describes.
- The field itself did move, just not net-favorably: shells at `s_gp≈0.20–
  0.30` improved ~2×, but `s_gp≈0.30–0.35` got ~2.1× *worse*, netting to a
  wash (proxy total ratio `1.07`).
- Checked `ww` (frame-dragging) and `velocity_sq` directly in the
  still-noisy `s_gp≈0.28–0.47` window: both are smooth and monotonic, no
  discontinuity — rules out a rotation-law-profile artifact as the residual
  driver.

### The sharpened, still-open puzzle (revised)

The EOS-knot-density mechanism itself is **retracted** (see above) — it
does not explain why DD2 differs from MPA1. What *does* still stand,
independent of that retraction:

- DD2's poor `Ham L2` (`~1.05–1.2`) is now confirmed real on the most
  rigorous footing available in this investigation: unaffected by table
  knot density (13 vs 258), unaffected by linear-vs-cubic lookup, and
  — critically — identical whether reached via single-jump or gradual
  6-step spin-up. None of the convergence-methodology concerns that
  retracted the MPA1 comparisons apply to this specific finding.
- The spatial localization (30–45% of stellar radius, not the surface) and
  the `alpha`/`ricci` term-decomposition mechanism (next section) were both
  measured on DD2 single-jump runs, which are now confirmed representative
  of DD2's true converged state — so neither is subject to the same
  retraction.
- What's now open, cleanly: DD2 has a real, structure-specific problem;
  EOS-table properties (surface level pending re-check, knot density
  retracted) don't explain it. The `alpha`-construction mechanism below is
  the strongest remaining lead.

---

## Term decomposition: `alpha`'s construction is the mechanism

Added `rho`, `gama`, `alpha`, `rhoH`, `ricci`, `KK` to `Cont/hamiltonian.dat`
(all already available inside `hamiltonian()`) to decompose
`ham = ricci + KK - 16π·rhoH·acoup4` directly, rather than infer which term
carries the noise.

### The finding

In the noisy window (`s_gp≈0.28–0.47`), relative 4th-difference smoothness:

| field | rel\|4th-diff\| |
|---|---|
| `rho` | `3.2e-8` |
| `gama` | `1.7e-8` |
| `KK` | `1.2e-5` |
| `alpha` | **`4.5e-6`** |
| `rhoH` (source term) | `2.3e-4` |
| `ricci` | **`2.8e-2`** |
| `ham` | `4.6e-2` |

`rho`/`gama` are essentially machine-smooth. `alpha` is ~500× rougher than
`rho`/`gama` despite feeding from the same underlying physics. `ricci`
(`= -2·gurr·(op%lap2(alpha) + 0.5·op%laplacian(logPsi4))`) is a further
~6,000× rougher than `alpha` alone. `ricci` and `16π·rhoH` are both `O(20–30)`
in this window and nearly cancel to produce `ham` — but `ricci`'s own
intrinsic roughness (not a cancellation artifact) propagates through
essentially unchanged in absolute terms (`ricci`'s roughness ≈ `0.84`,
`ham`'s ≈ `1.0`, same window) and accounts for nearly all of it.

### Why `alpha` specifically

Read `spin_integration_mod.f90::update_alpha_potential`. Unlike `rho`/
`gama`/`ww` (built from a spectral Legendre-multipole Green's-function sum,
`dgemm` against `P_2n`, infinitely smooth by construction), `alpha` is built
completely differently: a per-point nonlinear, **rational** (division-
containing) combination of *already-once-differentiated* cached quantities
(`dg_s_cache`, `dr_s_cache`, `d2g_ss_cache`, etc. — first and second radial
derivatives of `gama`/`rho`/`ww`, computed via the same high-order
finite-difference stencils used elsewhere), followed by a trapezoidal
**angular integration** (not a spectral sum) to build up `alpha(:,m)` from
`alpha(:,1)=0`. Two effects compound: (1) taking derivatives of already-
tiny-roughness fields amplifies that roughness (a generic property of
high-order finite-difference stencils — `deriv_s_and_ss_1d`'s interior
stencil is 6th-order/7-point, with coefficients up to `490` over `180·DS²`,
and `DS≈1/800` at this resolution gives a `~6.4e5` scale factor sensitive to
near-Nyquist content); (2) the nonlinear/rational combination in `da_dm`'s
formula can amplify further. `constraint_mod.f90`'s `op%lap2(alpha)` then
applies a *second* such derivative on top, which is where the `6,000×`
alpha→ricci jump comes from.

### Confirming this is not DD2-specific — but the resulting *severity* might be

Checked `alpha`'s own roughness for a normal, full-resolution MPA1 rotating
run (single-jump, so this specific number should be re-checked under
gradual spin-up too, per the methodology box — but the qualitative
construction-mechanism point does not depend on convergence quality):
`rel|4th-diff|(alpha) = 1.6e-4` — actually **~36× rougher than DD2's own
alpha** (`4.5e-6`) in this same diagnostic! Yet MPA1 converges far better
once properly (gradually) converged. This means the `alpha`-construction
amplification mechanism is universal (present in both EOS), and by itself
does not explain the MPA1/DD2 gap — DD2's problem is not simply "worse
`alpha` roughness." What differs between them is still open: possibly how
the *second* derivative operator (`ricci`'s `op%lap2`) interacts with each
star's specific radial `alpha` *profile* (not just its high-frequency
content) — e.g. a near-cancellation in `inv_denom` (the rational
combination's denominator, `m1·(1+s1·gs)² + numer_m²`) at a specific radius
that happens to coincide with DD2's stellar structure. Not yet checked.

---

## Separate issue, now fixed: low-density polytropic Γ was sub-linear

Distinct from the Ham L2 investigation above (that problem is localized to
the stellar *core*, `s_gp≈0.30-0.45`, far from this table's low-density
tail) — raised directly by the user while reviewing this EOS's TOV M-R
curve: the low-density polytropic extension's `Γ` (`d ln P/d ln n₀`,
confirmed as the physically correct basis over `d ln P/d ln E`, which only
coincides with it in the rest-mass-dominated limit) was **`0.9495`**, matched
at construction time to DD2's own soft edge slope. `Γ<1` is sub-linear —
pressure falls off *slower* than density as the star thins toward the
surface — the classical Newtonian Lane-Emden `n<0` regime, which has no
finite radius. Concretely: an M-R sweep (`GRASS_E_CENTER`-driven, 45
log-spaced points) showed 8 low-mass points (`M≲0.75 M☉`) failing to
converge within an 8s timeout, and the lowest-converged point (`R=25.2 km`)
still visibly growing as `M→0`.

For comparison, MPA1's own low-density plateau was measured the same way:
`Γ=1.58425`, exactly constant across ~4.5 decades in density
(`e≈1e3–2.4×10⁷ g/cm³`) — which turns out to be precisely the Read, Lackey,
Owen & Friedman (2009) shared low-density-crust polytropic index used
across their whole piecewise-polytrope EOS framework, not a coincidence
specific to MPA1.

**Fix**: in `fit_appendixA.py`'s Part 4, replaced the edge-slope-matching
`GAMMA_CRUST` computation with a fixed constant, `Γ=1.58425` (user-selected
from three literature-motivated options). Same stitch point (DD2's own
native CompOSE floor, `n0≈1e-12 fm⁻³`), same value-continuity matching (the
`Q7`-relative-shift formula from the earlier cancellation-bug fix, unchanged
— still needed regardless of which `Γ` is used), same extension range. The
seam now has a **genuine `Γ₁` discontinuity by design**
(`1.579→0.954` at the join) — standard piecewise-polytrope crust practice
(Read et al. 2009's own framework stitches this way), not a defect.

**Result** (same 45-point `GRASS_E_CENTER` sweep, `res=100`, before vs
after):

| | before (`Γ≈0.945`) | after (`Γ=1.58425`) |
|---|---|---|
| lowest converged `M` | `0.752 M☉` | `0.568 M☉` |
| `R` at lowest converged `M` | `25.20 km` | `17.48 km` |
| `R` at matched `M≈0.75 M☉` | `25.20 km` | `15.85 km` (−37%) |
| `M_max` | `2.3743 M☉` | `2.3743 M☉` (unchanged) |
| `R` at `M_max` | `12.38 km` | `11.86 km` |
| low-mass `R(M)` trend | still climbing at cutoff | monotonic, no upturn, across the whole converged range |

The radius-divergence problem is fixed: `R(M)` is now well-behaved and
bounded across the whole converged range, with no upturn as `M` decreases.
`M_max` is unchanged (untouched core physics, as expected — the change only
touches density below the table's own native floor). `R` at `M_max` shifted
by a small, understood, and expected amount (~4%): the table's global
`MN_GRASS` rescale (needed to keep `h[0]>1`) is smaller now
(`-1.28%` vs the old `-4.14%`, since the stiffer tail needs less correction),
and since that constant divides GRASS's internal enthalpy `h` at *every*
density, not just the tail, it shifts the whole `h(e)` mapping slightly —
the physical `e`/`p`/`n₀` columns outside the tail are untouched.

7 points below `e_c≈3.4×10¹⁴ g/cm³` (`M<0.57 M☉`) still fail to converge,
confirmed even with a 60s timeout (not a slowness issue) — a separate,
likely starting-guess-related limitation for very diffuse configurations,
not the `Γ<1` divergence this fix targeted, and not investigated further
here.

`E_SURFACE`/`P_SURFACE` in `para_panel.f90` updated to match: row 52 of the
regenerated table (`e=8.95 g/cm³`, still directly tabulated, no
extrapolation), `P=1.307030e18 dyn/cm²` — about 27× lower than the old
`3.589e19`, exactly as expected from a much steeper pressure drop-off with
the stiffer `Γ`.

Interactive before/after M-R comparison:
https://claude.ai/code/artifact/e4c51da3-4105-4ed8-8642-d5566e815fe0

### Confirmed: the Γ-tail fix does not touch the rotating Ham L2 problem

Built a fresh differentially-rotating star on the `Γ=1.58425` table to check
directly (not just assume, since the tail fix genuinely could have shifted
something table-wide via the `MN_GRASS` rescale — see above): standard TOV
(`solver_type="uniform"`, `run_mode=MODE_DEFAULT`, `r_ratio=1.0`, default
`e_center=0.8e15`, production `res=400`) → `Ham L2=4.957e-3`, clean. Spun up
via a single `MODE_REGRID` jump to `solver_type="uryu"`, `r_ratio=0.995`,
`lambda1=1.1`/`lambda2=0.9` (this investigation's established mild shape,
avoiding the `negative F_equator_h` stop) → **`Converged`, `Ham L2=1.1406`**.

No crash, no stop — and the number lands squarely in the same `~1.05-1.2`
range every other DD2 rotating test in this file has given, on every table
variant tried (446-row, 2441-row bulk-fixed, now the `Γ=1.58425`-tail
version). Exactly as expected: the low-density-tail fix targeted a genuine,
separate radius-divergence problem in the TOV branch and had no reason to
touch the core-localized (`s_gp≈0.30-0.45`) rotating residual, and it
didn't. This also re-confirms DD2's convergence-path-independence (single
big `r_ratio` jump, same as always used here) carries over to the new
table — no gradual-spin-up tie-breaker needed this time.

---

## A post-merger-like configuration, and a new numerical wall found

### Literature targets

Researched concrete parameter choices for a genuinely post-merger-remnant-
like Uryu-law configuration (not just the mild, stability-driven shape used
everywhere above). Choices made, with the specific reference and the exact
number/quote supporting each:

1. **Rotation-law shape** (`lambda1`=`λ₁=Ω_max/Ω_c`, `lambda2`=`λ₂=Ω_e/Ω_c`
   in this code — the parametrization is Uryū et al. (2017)'s own
   4-parameter U8 law, `p=1,q=3` fixed, matching this code's
   `uryu_p`/`uryu_q` defaults exactly):
   - **Origin of the law and its first post-merger-motivated numbers**:
     Uryū, Tsokaros, Baiotti, Galeazzi, Taniguchi & Yoshida (2017),
     *"Modeling differential rotations of compact stars in
     equilibriums,"* Phys. Rev. D **96**, 103011
     ([arXiv:1709.02643](https://arxiv.org/abs/1709.02643)) — introduces
     the U8/U9 rotation laws this code implements.
   - **Chosen value, `λ₁=1.5, λ₂=0.5`**: Cassing & Rezzolla (2024),
     *"Realistic models of general-relativistic differentially rotating
     stars,"* MNRAS **532**(1), 945
     ([arXiv:2405.06609](https://arxiv.org/abs/2405.06609)) — use this as
     their standard/reference shape, calibrated directly on the DD2 EOS
     (their reference Model Y: `M=2.29 M☉` on DD2, `p=1,q=3`).
   - **Alternative literature value, `λ₁=2.0, λ₂=0.5`** (more extreme,
     not the one chosen): Iosif & Stergioulas (2021), *"Equilibrium
     sequences of differentially rotating stars with post-merger-like
     rotational profiles,"* MNRAS **503**(1), 850-866
     ([arXiv:2011.10612](https://arxiv.org/abs/2011.10612)) — quote
     ("we adopt the choice {2.0, 0.5} as in Uryū et al. 2017") attributes
     this specific pair back to the original 2017 paper as *its own*
     post-merger-motivated example. Their Model C5 (`rp/re≈0.5`) gives a
     concrete `χ≈0.715` example.
2. **Mass/spin target**: Ma, Jiang, Wang, Jin, Fan & Wei (2018), *"GW170817
   and the Prospect of Forming Supramassive Remnants in Neutron Star
   Mergers,"* ApJ **858**(2), 74
   ([DOI:10.3847/1538-4357/aabafe](https://iopscience.iop.org/article/10.3847/1538-4357/aabafe)) —
   total binary mass `2.74 M☉`; remnant (after ~0.04 M☉ ejecta) widely
   cited as hypermassive, `~2.7-2.8 M☉`. Dimensionless remnant spin
   (distinct from the low pre-merger component spins) commonly quoted
   `χ≈0.6` for cases that collapse; Iosif & Stergioulas's Model C5
   (above) gives the concrete `χ≈0.715` example used as the upper end of
   the target range.
3. **Axis ratio**: `rp/re≈0.5-0.8` in both Iosif & Stergioulas (2021) and
   Cassing & Rezzolla (2024) above — this code's `r_ratio` is exactly this
   quantity.

### What was actually reachable

Deliberately avoided GRASS's `shoot`/`SHOOT_2D` mass+spin targeting mode
for this — no test anywhere in the repo actually drives `shoot_v2`'s
Newton loop, and an unrecognized `FIX2` string silently falls through to
targeting Keplerian break-up instead of erroring (found via code
inspection this session). Used only the established `OneModel`+
`MODE_REGRID`+gradual-stepping mechanism instead, accepting whatever
mass/spin comes out rather than shooting to hit a target exactly.

**New finding: `r_ratio≈0.80` is a wall for the mild shape, independent of
central density.** Tried gradual `r_ratio` spin-down at `lambda1=1.1,
lambda2=0.9` (the established-safe shape) at two different central
densities:
- Elevated (`e_c=2.0e15`, near this EOS's own TOV `M_max` density): reached
  `r_ratio=0.80` cleanly (`M=2.588 M☉`, `M₀=3.013 M☉` — already
  hypermassive, `χ=0.504`), but **every** further step tried
  (`0.78,0.76,0.75,0.74,0.72,0.70`) failed identically
  (`brent_core: failed to bracket root`) — not a step-size issue, since
  even the smallest step (`0.80→0.78`) failed the same way.
- Default (`e_c=0.8e15`): reached `r_ratio=0.85` cleanly (`χ=0.394`), but
  `r_ratio=0.80` failed (`negative F_equator_h`) at *this* density too —
  same wall, one step further out, confirming it's a `r_ratio`-driven
  limit rather than a mass/compactness-driven one.

**New finding: `lambda1≈1.3` is a wall for the shape strength, independent
of exact `r_ratio`.** At the default-density, `r_ratio=0.85` configuration,
stepped `(λ₁,λ₂)` toward the literature `(1.5,0.5)` target:
`(1.1,0.9)→(1.2,0.8)→(1.25,0.75)→(1.28,0.72)` all converged, each with only
marginal further increase in `χ` (`0.394→0.406→0.409→0.4091` — visibly
saturating); every attempt to push past this
(`(1.3,0.7)`, `(1.3,0.65)`, `(1.35,0.65)`) failed with `negative
F_equator_h`, regardless of step size or which of `λ₁`/`λ₂` moved. This
matches the *type* of failure already documented earlier in this file
(a genuine branch discontinuity in the Uryu equations for strong
differential rotation, previously only characterized near `r_ratio→1`) —
this session extends that finding to a boundary along the `λ₁` direction
as well, at a fixed, already-oblate `r_ratio`.

### The delivered configuration

Best reached, fully converged:

| quantity | value | literature target |
|---|---|---|
| `r_ratio` | 0.85 | 0.5–0.8 |
| `(λ₁, λ₂)` | (1.28, 0.72) | (1.5–2.0, 0.5) |
| ADM Mass `M` | 1.861 M☉ | — |
| Baryonic Mass `M₀` | 2.039 M☉ | ~2.7–2.8 M☉ (not hypermassive here) |
| `χ` (dimensionless spin) | 0.409 | ~0.6–0.7 |
| Central Ω | 679.0 Hz | — |
| Equatorial Ω | 488.9 Hz | — |
| `Ω_c/Ω_e` (`=1/λ₂`) | 1.389 | — |
| `T/W` | 0.0579 | — |
| Areal R | 14.54 km | — |
| Ham L2 | 0.356 | (rotating-branch residual, per this file) |

`Central Ω > Equatorial Ω` with `λ₁=1.28>1` confirms a genuinely
non-monotonic, off-center-peaked `Ω(r)` profile — qualitatively the right
*shape* for a post-merger remnant, even though the *magnitude* of
differential rotation and the spin `χ` fall short of the literature
targets. Built at the default (non-hypermassive) central density rather
than the elevated one, since the `r_ratio`-wall problem made the elevated-
density, genuinely-hypermassive branch unreachable past `r_ratio=0.80`
with this rotation-law implementation.

**Not resolved**: this is a genuine, currently unfixed limitation of the
Uryu-law implementation for strong differential rotation — both along
`r_ratio` (near `0.80`, any density) and along `λ₁` (near `1.3`, at
`r_ratio=0.85`). A literature-faithful post-merger profile (`λ₁≥1.5`,
`r_ratio≤0.8`, hypermassive) is not currently reachable via gradual
stepping from a TOV seed with this implementation — reaching one would
need either fixing the underlying branch-discontinuity bug (root cause
still uncharacterized beyond "near `r_ratio→1`," now also seen at fixed
`r_ratio=0.85` along `λ₁`) or attempting the untested `SHOOT_2D` targeting
path.

### References (this section)

- Uryū, K., Tsokaros, A., Baiotti, L., Galeazzi, F., Taniguchi, K., &
  Yoshida, S. (2017). *Modeling differential rotations of compact stars in
  equilibriums.* Phys. Rev. D, 96, 103011.
  [arXiv:1709.02643](https://arxiv.org/abs/1709.02643)
- Iosif, P., & Stergioulas, N. (2021). *Equilibrium sequences of
  differentially rotating stars with post-merger-like rotational
  profiles.* MNRAS, 503(1), 850-866.
  [arXiv:2011.10612](https://arxiv.org/abs/2011.10612)
- Cassing, M., & Rezzolla, L. (2024). *Realistic models of
  general-relativistic differentially rotating stars.* MNRAS, 532(1), 945.
  [arXiv:2405.06609](https://arxiv.org/abs/2405.06609)
- Ma, P.-X., Jiang, J.-L., Wang, H., Jin, Z.-P., Fan, Y.-Z., & Wei, D.-M.
  (2018). *GW170817 and the Prospect of Forming Supramassive Remnants in
  Neutron Star Mergers.* ApJ, 858(2), 74.
  [DOI:10.3847/1538-4357/aabafe](https://iopscience.iop.org/article/10.3847/1538-4357/aabafe)

---

## Files changed

### GRASS (`/u/tlam/GRASS`)

- `src/para_panel.f90` — `r_ratio` default, `E_SURFACE`/`P_SURFACE`
  (currently `8.82 g/cm³` / matching pressure — the near-standard,
  seam-smoothed table's own row), `eos_file` (currently
  `"DD2_hot_equal_pinned"`).
- `src/core/starting_model_mod.f90` — removed the `r_ratio` hard-coded
  overwrites (case-default seed guess, `single_model`'s `0.7` clobber);
  `MODE_REGRID` branch now correctly re-imposes a *target* `r_ratio`
  captured before `regrid_read` can clobber it, instead of the old
  hard-coded `e_center=0.74e15` override.
- `src/theory/rotational_law_mod.f90` — new shared `safe_F_e` function;
  applied in `diff_rotation_uryu` and `cache_uryu_ab`.
- `src/core/constraint_mod.f90` — enabled the (previously dead-code-gated)
  2D Hamiltonian-field dump to `Cont/hamiltonian.dat`, added `r_phys`,
  `enthalpy`, `ww`, and `velocity_sq` columns to it (diagnostic only, still
  enabled; `ww`/`velocity_sq` were already imported into the module, so this
  was a one-line addition to the `write` statement).

### PyCompOSE (`/sakura/ptmp/tlam/PyCompOSE/examples/DD2_hot`)

- `fit_appendixA.py`:
  - Reverse-pass low-density construction re-ordered to prefer staying at
    the table's coldest tabulated `T` (using DD2's own physical NSE crust),
    falling back to warming only as a last resort (previously the reverse:
    entropy-pin/warm first).
  - Added a smooth `T`-glide transition window (`WIN=20` rows) approaching
    the `i_min` boundary, replacing an abrupt jump to the untouched,
    original Appendix-A-fit point.
  - Added a cold degenerate polytropic low-density extension (120 rows)
    below the table's own floor, matched continuously in `Q7` (not
    absolute energy density — avoids catastrophic cancellation) with
    `Gamma` fit to the table's own edge slope, extending the GRASS-format
    output down to `n0=1e-16 fm⁻³` (`e≈0.165 g/cm³`) so a standard
    near-vacuum surface cutoff is directly tabulated, no extrapolation
    needed.
  - `write_grass_table` call now writes the full extended (446-row, later
    2441-row) table.
  - **Bulk (`i>=i_min`) knot-density fix**: for the core/bulk region (where
    `T_full`/`Yl_full` are the smooth `T_model()`/`Yl_model()` fit curves,
    not the crust's discrete entropy-pin/T-glide output), refined each
    inter-knot interval 20× (`BULK_REFINE=20`) and re-looked-up
    `Q1`/`Q2`/`Q7` via `eos.eval_given_rtx(..., method="cubic")` instead of
    the default `"linear"`. 13→258 core knots; confirmed monotonic, seam-
    smooth, `e(h)`-preserving. Did not improve `Ham L2` on DD2 (see
    "Localizing the residual" above) — kept in the pipeline anyway since it
    is a genuine, unrelated improvement in table quality, just not the fix
    for this problem.
  - **`GAMMA_CRUST` fixed at `1.58425`** (was: matched to DD2's own soft
    edge slope, `≈0.9495`) — see "Separate issue, now fixed" above. Same
    stitch point and `Q7`-shift continuity formula; only the constant's
    value and its derivation (now a literature/MPA1-anchored fixed value,
    not a per-run 3-point derivative of DD2's own data) changed.
- New diagnostic scripts (all still present, not part of the regular
  pipeline): `resample_grass_table.py` (re-run this session against the
  current table each time it changes — it is a snapshot, not live),
  `compare_eos_derivatives.py`, `check_first_law.py`.

### GRASS `eos/` directory (untracked test tables)

- `MPA1_coarse.dat` (new) — MPA1.dat decimated to DD2-like knot sparsity
  (458 rows, every 13th row of the original 5931) for the mirror test
  above; real subsampled data, not resampled.
- `DD2_hot_equal_pinned_fine.dat` — PCHIP-resampled 6000-row version of
  whatever `DD2_hot_equal_pinned.dat` was at the time it was last generated
  (currently the *old*, pre-bulk-fix 446-row table — stale relative to the
  current 2441-row `DD2_hot_equal_pinned.dat`; re-run
  `resample_grass_table.py` before reusing it for anything).

## Current repo state

- `/u/tlam/GRASS`: git-dirty with the four source files below
  (`constraint_mod.f90` now also carries `ww`, `velocity_sq`, `rho`, `gama`,
  `alpha`, `rhoH`, `ricci`, `KK` — 13 columns total); `git diff --stat` shows
  the exact changes. `eos/DD2_hot_equal_pinned.dat` (2441 rows, untracked,
  regenerated with the `Γ=1.58425` tail — same row counts as the earlier
  bulk-fixed version, different tail content), `eos/DD2_hot_equal_pinned_fine.dat`
  (stale — resampled from the *pre-Γ-fix* table; re-run
  `resample_grass_table.py` before reusing it), `eos/MPA1_coarse.dat`
  (untracked), and `eos/DD2_hot_equal.dat` (the older, pre-investigation
  table, untracked) all present.
- `para_panel.f90` last left at: `eos_file="DD2_hot_equal_pinned"`,
  `solver_type="uniform"`, `run_mode=MODE_DEFAULT`, `r_ratio=1.0`,
  `E_SURFACE=8.95`/`P_SURFACE=1.307030e18` (current, matching the `Γ=1.58425`
  table), `res=400` (production resolution, restored after the M-R sweep's
  `res=100` speed hack) — a TOV, not rotating, configuration. Switch to
  `solver_type="uryu"`, `run_mode=MODE_REGRID`, `r_ratio=0.995` for any
  further rotating test, same as throughout this investigation.
- `Res/res.rst` holds whatever the M-R sweep's last point left
  (`e_c=9e15 g/cm³`, `M≈2.03 M☉`) — not a deliberately-prepared restart.
  Regenerate a fresh TOV restart (`GRASS_E_CENTER` unset, uses
  `e_center_default=0.8e15`, or set explicitly) before chaining into a
  rotating test. Frozen copies of both EOSs' clean TOV restarts from
  earlier this investigation are saved in the scratchpad (see below) —
  note `res_DD2_TOV_frozen.rst` predates the `Γ` fix and is now stale for
  DD2; `res_MPA1_TOV_frozen.rst` is still valid (MPA1 untouched).
- Run with `MKL_THREADING_LAYER=GNU` from `/u/tlam/GRASS` (relative
  `./eos/`, `./Cont/`, `./Res/` paths require this cwd).
- Saved field/restart snapshots for every run this session, in the
  scratchpad: `hamiltonian_TOV_baseline.dat`, `hamiltonian_ROT_baseline.dat`,
  `hamiltonian_TOV_fine.dat`, `hamiltonian_ROT_fine.dat`,
  `hamiltonian_MPA1coarse_ROT.dat`, `hamiltonian_TOV_bulkfix.dat`,
  `hamiltonian_ROT_bulkfix.dat`, `hamiltonian_ROT_decomp.dat` (the 13-column
  term-decomposition run), `res_MPA1_TOV_frozen.rst`,
  `res_DD2_TOV_frozen.rst` — useful for follow-up analysis without
  re-running GRASS, and the two frozen restarts specifically enable a
  matched-convergence-path re-test without rebuilding a TOV star from
  scratch each time.

## AthenaK evolution tables: 3D trapped-neutrino DD2 + 1D initial-data slice

Two new `.athtab` files were built for AthenaK's `pgen/grass` branch, so a
GRASS-built star can actually be *evolved* with a genuine finite-temperature,
composition-dependent EOS (`<mhd> dyn_eos = compose`) instead of just reusing
GRASS's own barotropic `(e,p,h,n0)` pair as static initial data:

- `/u/tlam/GRASS/eos/DD2_hot_trapped_3d.athtab` (96 MB) — real 3D
  `(nb,yq,t)` trapped-neutrino DD2 table (`examples/DD2` with electrons +
  `add_trapped_neutrinos()`), for the runtime evolution EOS.
- `/u/tlam/GRASS/eos/DD2_hot_slice.athtab` (91 KB) — 1D `(nb)`-only slice
  (`Q1`, `Q7`, `temp`, `Y[e]`) along the SAME `(nb,T,Yl)` trajectory that
  produced `DD2_hot_equal_pinned.dat`, used to seed the passive-scalar `Y[e]`
  when installing initial data from a GRASS restart, consistently with the
  pressure/energy profile GRASS actually solved its stellar structure
  against.

Generating scripts (in `/sakura/ptmp/tlam/PyCompOSE/examples/DD2_hot/`):
`build_dd2_hot_trapped_3d_athtab.py`, `build_dd2_hot_slice_1d_athtab.py`,
`verify_mn_scalar.py`.

**Baryon-mass (`mn`) consistency**: both tables use `mn=mp=MN_GRASS =
939.5660386746146 MeV` (the same effective mass GRASS's own
`MB=1.6749286e-24 g`, `para_panel.f90:164`, implies) — rescaled from the raw
3D table's native amu convention (`mn_old=931.494061` MeV) via the standard
`Q_new=(Q_old+offset)*mn_old/mn_new-offset` identity (keeping every physical
quantity invariant). Verified, not just asserted: `verify_mn_scalar.py`
greps `kGrassMB` directly out of AthenaK's `grass_units.hpp` and confirms
it's byte-identical to the Python side's `GRASS_MB_G` — no change was needed
there.

**dfloor / minimum-n0 consistency**: the 1D slice's `nb` range is clipped to
exactly `[1e-12, 10]` fm⁻³ — the real 3D table's own native floor/ceiling,
*not* GRASS's own sub-floor `Γ=1.58425` cold-polytrope crust extension
(that crust exists only in `DD2_hot_equal_pinned.dat`, for building the
star's low-density envelope, and has no grounding in the real 3D EOS). The
corresponding AthenaK `<mhd> dfloor`, derived from AthenaK's own
`unit_system.hpp` constants (`<mhd> units=geometric_solar`, the default):

```
GetBaryonMass()-equivalent(mn=MN_GRASS), code units = 0.002711994197437188
dfloor_code = nb_min_3D(=1e-12 fm^-3) * GetBaryonMass()-equivalent
            = 2.711994197437188e-15
```

confirmed empirically by an actual smoke-test run: AthenaK's own
`tov::TabulatedEOS` printed `rho = [2.71199e-15, 0.0271199]` on loading
`DD2_hot_slice.athtab` — an exact match. Use `dfloor >= 2.8e-15` (or the
raw number above) in any `.athinput` combining these tables with
`<mhd> units=geometric_solar`; redo the derivation with the corresponding
`Make*()` constants for any other unit choice.

**Known, accepted approximation**: `add_trapped_neutrinos()`'s Yq-axis
reinterpolation touches `Q3/Q4/Q5` (`mu_b`,`mu_q`,`mu_l`) but adds no
explicit neutrino chemical-potential term to them (only `Q1,Q2,Q6,Q7` get
one, via `add_equilibrium_neutrinos()`). AthenaK's `EOSCompOSE` reader does
consume `Q3/Q4/Q5` as real chemical potentials — accepted as-is since
nothing in the current `dyngr_grass` pipeline uses them for source terms.

**AthenaK-side wiring** (branch `pgen/grass`): `grass_reader.hpp` gained a
`Yq` field on `Point` and a `const tov::TabulatedEOS &slice_eos` parameter
on `Interpolate`, set via `slice_eos.GetYeFromRho<Host>(out->rho0)`;
`dyngr_grass.cpp` constructs `slice_eos` from `<problem> table=...` and
writes it to `IYF` when `<mhd> nscalars>0` — same pattern
`kadath_bns.cpp` uses. Compiles clean; two smoke tests
(`inputs/dyn_grmhd/dyngr_grass_{tov,rot}_compose_test.athinput`, run from
`/sakura/ptmp/tlam/athenak_run/pgen_grass_smoke/`) both pass with
`rho-max`/`alpha-min` matching the pre-existing ideal-EOS baseline almost
exactly (regression check) and a nonzero seeded `Y[e]` scalar in both.

## Next steps, in priority order

1. **Re-verify the MPA1 surface-level result (finding 8) under gradual
   spin-up**, the same way the knot-density finding was re-checked and
   retracted. Use `res_MPA1_TOV_frozen.rst`, raise `E_SURFACE`/`P_SURFACE`
   to DD2's cutoff level (`e≈2000 g/cm³`, row 130 of `MPA1.dat`, as before),
   and run the same 6-step gradual chain instead of a single jump. If it
   still lands near `~150`, the surface-level mechanism survives and is
   genuinely confirmed. If it lands near MPA1's normal `~5e-3`, it joins
   the knot-density finding as retracted, and every surface-level-based
   conclusion in this file (including DD2's own extensive surface-level
   remediation work) needs re-examination.
2. **Push on the `alpha`-construction mechanism specifically for DD2.** The
   term decomposition (previous section) is the strongest lead that has
   *not* been undermined by the convergence-methodology finding. Concretely:
   dump `da_dm` and `inv_denom` (both local variables in
   `update_alpha_potential`, not currently exposed) at the same radii,
   across both EOS, to check whether `inv_denom`'s denominator
   (`m1·(1+s1·gs)² + numer_m²`) comes anomalously close to zero for DD2's
   specific stellar structure in the `s_gp≈0.30–0.45` band — a near-
   singular divide would explain why DD2 (not just "an EOS with somewhat
   rougher alpha," which MPA1 already is) is uniquely bad.
3. If (2) doesn't localize it further, consider whether `alpha`'s
   trapezoidal (not spectral) angular integration is itself the weak link
   — e.g. temporarily switching it to use the same Legendre-basis machinery
   already used for `rho`/`gama`/`ww` (if that's physically valid for this
   equation) as a diagnostic, not necessarily a real fix.
