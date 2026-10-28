# Batch 15 — B/C UV-running scalar-mass feasibility, before any more fits

## Scope

The latest Batch 14 live replay for B(-1) / C(-1) found

| Model | mS,max(100 TeV) | MF,min(100 TeV) | mS,max/MF,min | Former exploratory chi2 |
|---|---:|---:|---:|---:|
| B(-1) | 145.31 TeV | 84.395 TeV | 1.7218 | 18,376.7 |
| C(-1) | 180.64 TeV | 104.30 TeV | 1.7319 | 20,279.3 |

The same fitted Yukawa point and a zero-extra-scalar-quartic boundary give
relative scalar mass-squared changes of only ~10^-12, so adding those
small quartics at the boundary does not cause the hierarchy inversion.

**Verified source detail (GitHub main a05b3b7):** In the native projected
UV mass-squared beta expressions, the terms quadratic in the heavy-fermion
mass and Yukawas occur with coefficient **-8** on each of two `Tr[...]`
structures for B(-1), and **-12** for C(-1). The Yukawa-dependent `mS_i^2`
term has coefficient 2 (B) / 3 (C) multiplying `Tr[y_i.Trans[Bar[y_i]]]`.
The convention is `16*pi^2 dm_i^2/dln(mu)=beta^(1)`. Because running is
from 10^7 **down** to 10^5 GeV, negative MF²Y² beta contributions can
*increase* positive running scalar mass-squared. This is why scanning the
Yukawas' overall scale is more informative than changing small quartic inputs.

## What the new script genuinely does

**New source:** `scripts/StudyT3BCUVScalarMassFeasibility.py`.

1. Read saved Batch 13 / Batch 12 general-scalar UV-running best Yukawas from
   `general_fresh_v1_B_m1_cont` / `general_fresh_v1_C_m1_cont` using the
   corrected Batch 14 loader. Verify every saved source SHA-256.
2. Copy the historical initial state *in memory*. Apply a uniform factor
   `f_y` to all 18 real fitted Yukawas, and a positive factor `f_S` to both
   historical initial scalar masses: `mS_i^2(UV) -> f_S^2*mS_i^2(UV)`.
   No negative scalar mass input, new fermion masses, or modified beta is used.
   Each general quartic remains `0.1*abs(lambdaT3)` at the UV boundary,
   or the actual source-fit value if it differs.
3. Reintegrate native candidate UV RGBeta (`10^7 -> 10^5 GeV`) for each point,
   including running of the full Yukawas and MF. Diagonalise the candidate
   two-doublet 2x2 scalar quadratic matrix including its actual generated
   `mu2` and compare to the smallest heavy-fermion singular mass at `muF`.
4. **If and only if** the candidate has positive scalar eigenvalues and
   `mS,max(muF) < MF,min(muF)`, run the existing candidate scalar-only EFT1
   (`10^5 -> 10^3 GeV`), retaining the additional quartic and quadratic
   evolution. Otherwise mark the inversion and skip the invalid EFT1 stage.
5. Record two independent screening conditions:
   - conservative F-first mass ratio at muF: `mS,max/MF,min <= 0.5`;
   - both positive running scalar masses at muS inside `[0.5,5] TeV`.
   Both thresholds are adjustable *screening conventions*, not physical pole
   threshold prescriptions.
6. At the baseline (`f_y=f_S=1`) separately evaluate the **actual native UV**
   `mS1Sq` and `mS2Sq` beta expressions at the UV start and `muF`, substituting
   `(y1,y2) -> 0` and `MF -> 0` at the *same frozen state*. Inclusion–exclusion
   gives exact local contributions: no-(Y/MF), Yukawa-only at MF=0,
   MF-only at Y=0, and joint Y–MF dependence. They sum algebraically to the
   full beta. These are not unique diagrams and do not integrate by themselves
   to the net mass change. The baseline must reproduce the previous fitted
   MF(muF), and will independently verify the Batch 14 scalar mass at muF
   when `uv_mass_hierarchy_v1/diagnostic.json` is present.

**Not done:** no new C5 evaluation, no new chi-square minimisation, no
physical pole masses, no sequential thresholds, no new certification of the
Matchete/RGBeta coupling map. These are **potential starting points** for
future constrained Yukawa refits, not accepted neutrino fits. A uniform
Yukawa reduction generally changes the neutrino observables drastically.

If the scan has zero successes, it does not establish a no-go theorem: it
only excludes the finite positive input-factor grid under this candidate RGE.
For tests reporting numerical failures, inspect the saved reason rather than
interpreting these points as a physical exclusion.

## Files to extract (complete files, no patches)

- `scripts/StudyT3BCUVScalarMassFeasibility.py`
- `tests/python/test_t3_bc_uv_mass_feasibility.py`
- `docs/T3_BC_UV_SCALAR_MASS_FEASIBILITY_BATCH15.md`

The prior Batch 7, 9, 11, 12 and **fixed** Batch 14 scripts must already be
installed, and the existing source beta JSON exports must be available.

## PowerShell, project root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_BC_UV_Scalar_Mass_Feasibility_Batch15.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_bc_uv_mass_feasibility.py tests\python\test_t3_bc_uv_scalar_hierarchy.py

# Default: B and C, 4 Yukawa factors x 3 scalar-mass factors = 24 UV replays.
# Only positive, F-first ordered endpoints proceed into EFT1.
python -m scripts.StudyT3BCUVScalarMassFeasibility --run-id uv_mass_feasibility_v1

Get-Content -Encoding UTF8 output\studies\general_scalar_extended\diagnostics\uv_mass_feasibility_v1\summary.md
```

For a smaller B-only pilot with exact historical baseline included:

```powershell
python -m scripts.StudyT3BCUVScalarMassFeasibility `
  --model B_-1 --yukawa-factors "0.1,0.5,1" `
  --uv-scalar-factors "0.5,1" `
  --run-id uv_mass_feasibility_b_pilot_v1
```

## Outputs (new directories only)

`output/studies/general_scalar_extended/diagnostics/<new-run-id>/`

- `screen.json`: all candidates, source hashes, running-mass screens and
  local UV mass beta dependence, with no invented chi2 scores.
- `scan.csv`: one row per trial, ready for inspection/plotting.
- `summary.md`: threshold screen results and beta-sector table.

Existing historical files, previous fits and Batch 14 results are never written.
Run ID collisions are refused. The script does not edit GitHub.

## What to do with the output

If the screen finds well-separated, TeV-scale scalar candidates, pick the
candidate closest to the saved fit Yukawa direction and construct a **new**
constrained 18-real-Yukawa fit with those UV masses fixed or included in a
controlled optimisation. Require both `chi2 <= 5.86` and mass ordering on
recomputed trajectories. This requires a separate follow-up batch; do not
read the old chi2 as the new candidate's fit quality.

If no positive-UV-mass configuration reaches both screens, a more general
scan of Yukawa textures and/or heavy-fermion/scalar input masses might help,
but those are different benchmark definitions. If the running mass hierarchy
remains inverted near decoupling, a proper scalar-first or near-degenerate
matching construction is needed; more optimisation alone cannot cure the EFT.
