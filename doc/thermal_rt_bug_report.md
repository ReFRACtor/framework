# Thermal Radiative Transfer: Bugs Found and Fixed

Branch: `thermal_rt_followup`, a follow-up to `plans/thermal-rt-coverage-handoff.md`.

This branch added tests that check thermal radiances against independent
solutions and check Jacobians against finite differences (FD). Those tests
turned up a series of bugs in the ReFRACtor drivers and in the third-party
LIDORT 3.8 / VLIDORT 2.8 Fortran code. Every Fortran change is marked in the
source with a `! ReFRACtor fix:` comment.

## Commits

| Commit | Scope |
|---|---|
| `a2d3b90` | Handoff item 1: the VLIDORT vector "multiple-scattering shortfall" (turned out to be a bad reference) |
| `963f0a8` | Item 2a: downwelling thermal output; LIDORT/VLIDORT/FO thermal Jacobian fixes |
| `3301a77` | Items 2b+3: Cox-Munk thermal tests; FO emissivity and BRDF emissivity linearization fixes |
| `b22c557` | Item 2c: Python integration tests; fixes to the SpurrRt chain rule, sphericity setup, VLIDORT Pi matrix and FO |

## Important deployment note

**LIDORT is not built from `thirdparty/lidort-3.8` by default.** CMake finds a
prebuilt `liblidort.so` in the conda/pixi environment instead. The LIDORT
Fortran fixes below only take effect if one of these is done:

- build with `BUILD_LIDORT=ON`, or
- release an updated LIDORT package.

For testing on this branch, LIDORT was built from source and copied into the
environment. The original library was kept as `liblidort.so.3.8.conda-orig`.

VLIDORT *is* built from source (as an ExternalProject), so the VLIDORT fixes
apply as soon as it is rebuilt.

---

## 1. Not a bug: the VLIDORT vector thermal reference was wrong

The handoff reported that refractor's VLIDORT vector (nstokes=3) result for the
Rayleigh thermal scenario was missing the multiply-scattered field. The old
offline reference gave I=0.014753 and Q=-8.15e-4, while refractor gave
I≈0.013937 and Q≈0.

The reference was wrong, not refractor:

- **Same tester, scalar mode:** the 2.7-based offline tester run at nstokes=1
  gives I=0.0139377, which agrees with refractor.
- **Energy:** a conservative layer with τ=0.02 cannot brighten the surface
  emission by 6%.
- **Polarization:** at Fourier order m=0, single scattering of an isotropic
  upwelling hemisphere produces Q≡0.
- **Independent solver:** a new vector successive-orders-of-scattering (SOS)
  solver (`reference/rt/vlidort_thermal_test/sos_thermal_check.py`) agrees with
  refractor to about 1e-4 in I and gives |Q|≈1e-9.

Test changes:

- The vector thermal test now checks I, Q and U for all scenarios against the
  SOS solver.
- The scalar thermal tolerances were tightened from 7e-3 absolute (about 50% of
  the signal) to 1e-5.

A separate observation: stock VLIDORT 2.8.3, using its internal Lambertian
surface path, returned direct transmission only for this case. Refractor does
not use that path (it always uses the BRDF supplement), so this was not
investigated further.

## 2. Fortran bugs

### LIDORT 3.8

| # | File | Bug | Effect |
|---|---|---|---|
| L1 | `lidort_lbbf_jacobians.f90` | `UBRDF_F` was declared `(MAXSTREAMS, MAX_USER_STREAMS)`, transposed relative to the caller's `USER_BRDF_F` | Upwelling atmospheric blackbody (BB) Jacobians lost most of the surface-reflected term; about 22% error with albedo 0.3 |
| L2 | `lidort_lbbf_jacobians.f90` (whole layers and partial layers) | In the 2-layer BVP, the bottom (BOA) level used the upper-level `Gp1` coefficient instead of `Gp2` | Wrong BOA-level BB Jacobian whenever there is scattering; about 0.3% error |
| L3 | `FO_Thermal_RTCalcs_ILPS/ILCS.f90` | In the enhanced-pseudo-spherical path, the transmittance derivative used the vertical distance `zj` instead of the line-of-sight distance `dj` | Wrong downwelling (and off-nadir upwelling) atmospheric Jacobians; about 5% error |
| L4 | `brdf_lin_sup_masters.f90` | The emissivity derivative with respect to the kernel factor, `-LOCAL/FF`, is NaN when FF=0 | NaN surface Jacobians for Cox-Munk with a zero Lambertian component. Now uses exactly -1 for the Lambertian kernel and guards FF=0 |

### VLIDORT 2.8

| # | File | Bug | Effect |
|---|---|---|---|
| V1 | `VFO_Thermal_RTCalcs_ILPS/ILCS.f90` | Same as L3 | Same as L3 |
| V2 | `vlidort_lbbf_jacobians_vector.f90` | Same as L2 | Same as L2 |
| V3 | `vlidort_lps_masters.f90`, `vlidort_lcs_masters.f90` (all 8 call sites) | The BB-Jacobian (LBBF) routine was passed `PI_XQP_PRE` where it expects `PI_XQM_PRE` | Wrong upwelling BB Jacobians with anisotropic scattering; up to about 17% error |
| V4 | `vlidort_lp_bvproblem.f90` | Thermal particular-solution linearizations were copied for *all* layers instead of only the layer being varied (LIDORT does this correctly) | Spurious boundary-condition terms in the profile Jacobians of every layer above the bottom; about 17% error with a reflecting surface |
| V5 | `vlidort_ls_wfsurface.f90` | The user-stream and thermal-transmittance-only emissivity linearization was subtracted instead of added | Wrong surface Jacobians in thermal mode; about 0.06·B error |
| V6 | `vbrdf_lin_sup_masters.f90` | The emissivity linearization with respect to the kernel factor and kernel parameters was off by a factor of FF (upstream 2.8.3 assumed the local emissivity excludes FF; it includes it) | Wrong surface Jacobians. Restored the form LIDORT uses, plus the FF=0 guard from L4 |
| V7 | `vbrdf_sup_routines.f90`, `vbrdf_lin_sup_routines.f90` | The Lambertian emissivity branch set only Stokes component 1; components 2..n stayed uninitialized stack memory | Garbage Q/U in vector thermal results, with values that depended on what ran before (failures depended on test order) |
| V8 | `vlidort_lps_masters.f90` | `DO_RTSOL_VARY` and `NPARAMS_VARY` were uninitialized when atmospheric (e.g. BB) linearization was on but profile linearization was off | Out-of-range LAPACK calls or memory corruption |
| V9 | `vlidort_pars.f90.in` | `MAX_SURFACEWFS = 1` | Any Cox-Munk surface Jacobian (3 weighting functions) indexed out of bounds. Raised to 7 to match LIDORT |

## 3. ReFRACtor C++ driver bugs

| # | Location | Bug | Effect |
|---|---|---|---|
| R1 | `SpurrRt::stokes_and_jacobian_single_wn` | The atmospheric BB chain rule looped over *layers*, but the blackbody is defined on *levels*, so the surface-level term was dropped | About 25% error in the temperature-offset Jacobian; the surface-pressure Jacobian was also wrong (affects every thermal retrieval) |
| R2 | `MultiScattRtDriver::setup_sphericity` | In thermal-only mode, the first-order (FO) correction was left on with neither FO type selected | LIDORT rejects this input, so every non-plane-parallel thermal `LidortRt`/`VlidortRt` failed. Thermal-only now turns FO correction off, because the FO thermal code has no BB Jacobians |
| R3 | `MultiScattRtDriver` / `VLidortRtDriver::clear_linear_inputs` | The BB-Jacobian (LBBF) flags and VLIDORT's `do_linearization` were left on after a Jacobian call | A later radiance-only call ran a half-configured linearized solution (this is what triggers V8) |
| R4 | `FirstOrderDriver` | The BRDF supplement ran in direct-bounce-only mode, which skips the Kirchhoff emissivity calculation | FO thermal emissivity was 1 for every surface. A Lambertian special case had hidden this (handoff item 3). Fixed, and the special cases were removed |
| R5 | `FirstOrderDriver::setup_thermal_inputs` | Emissivity was read before `SpurrRt` runs the BRDF supplement | The emissivity lagged one call behind, so the first call was wrong. It is now read in `calculate_rt` |
| R6 | `FirstOrderDriver::copy_jacobians` | BB Jacobians were left empty, but `SpurrRt` still indexes them | Undefined behaviour. They are now computed using the fact that the radiance is linear in the BB inputs |
| R7 | Tests | The FO-vs-LIDORT thermal radiance tolerance was 6e-3 absolute | This hid R4. Tightened to 1e-7 |

## 4. New capability

Downwelling output (an up-looking observer at the bottom of the atmosphere):

- `MultiScattRtDriver::set_downwelling()` / `downwelling()` for LIDORT and
  VLIDORT.
- `FirstOrderDriver::set_downwelling()` / `downwelling()`, thermal-only (it
  throws if solar sources are on).
- Both are exposed in Python. The setting is serialized as part of the RT
  interface state.

## 5. New tests

- `lib/RadiativeTransfer/thermal_jacobian_check.h` checks Jacobians against FD
  for the atmospheric BB, surface BB, optical depth and surface parameters. It
  covers:
  - up and down output;
  - Lambertian and Cox-Munk surfaces;
  - absorbing and scattering atmospheres;
  - 2 and 5 layers;
  - an anisotropic Henyey-Greenstein phase function.

  It runs for LIDORT, VLIDORT (nstokes 1 and 3) and the First Order driver.
- Downwelling radiances are checked against the SOS solver and a closed-form
  result, along with LIDORT↔VLIDORT and FO↔LIDORT consistency.
- Cox-Munk emissivity is checked between FO and LIDORT.
- `python/test/thermal_rt_test.py` runs integration-level state-vector Jacobian
  FD checks through `LidortRt`, `VlidortRt` and `FirstOrderRt`.

## 6. Known remaining limitations

- **Solar+thermal mode:** FO correction stays on for the solar part, and the
  FO thermal code has no BB Jacobians. As a result, BB Jacobians in combined
  mode miss the direct-emission term.
- **FF=0 for non-Lambertian kernels:** in the BRDF supplement emissivity
  linearization, a non-Lambertian kernel with kernel factor FF=0 now gets a
  derivative of 0 instead of NaN, but that zero is not the correct derivative.
- **Out of scope** (per the handoff): 2Stream BB Jacobians and FO
  pure-Rayleigh thermal.
