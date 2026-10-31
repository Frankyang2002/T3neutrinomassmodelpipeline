# B/C running-input one-scale leading triangle — v1

## Intended improvement over overnight v1

The previous overnight campaign independently passed (1) 16-model UV RGBeta
integration, (2) B/C five-mass UV→EFT1 integration, and (3) **frozen** B/C
Matchete v23 matching. Those three independent successes did *not* prove a
single coupled trajectory was used in the matching.

This package takes the **same 1e5 GeV UV endpoint** from the B/C full-scalar-mass
runner for *all* mass/Yukawa/quartic inputs and evaluates the already cross-
checked **480-term direct all-heavy ordered Matchete coefficient** in that
common-scale input state. It does not substitute EFT1 1 TeV scalar masses
into the UV Yukawas or reconstruct a missing EFT Wilson boundary.

## New complete files, relative to local repository root

- `Numerical/experimental/T3BCRunningOneScaleTriangle.py`: input and symmetry
gates, complex symmetric MF Takagi rotation, full three-heavy-column Yukawa
lift, scalar-mass eigenbasis and action checks, and candidate source C5/mν.
- `scripts/RunT3BCRunningOneScaleTriangle.py`: load and hash-check B/C
overnight numerical source; import and verify archived direct v22 Matchete
source and numerical coefficient; compare the **frozen** v22 C5 to a freshly
recomputed C5; then evaluate the same 480-term expression at the running
UV endpoint. It also recomputes the independent F/Y/X triangle kernel at
those **running** masses, quartics and U, testing direct-vs-sequential
ordered coefficient parity for each heavy column. Optional scan reruns B/C full coupled UV→EFT1 RGBeta with three
additional mixed quadratic UV settings.
- `scripts/RunT3BCRunningOneScaleTriangle.ps1`: tests and one full B/C run
with mixed-bilinear numerical scenarios.
- `tests/test_t3_bc_running_one_scale_triangle_v1.py`: source-independent
complex matrix/Takagi/C5 contracts and live archived matching parity when
installed in the complete repository.
- `docs/T3_BC_RUNNING_ONE_SCALE_TRIANGLE_V1.md`: this document.

## Installation and execution (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_running_one_scale_triangle_v1.zip" -DestinationPath . -Force
python -u -X utf8 -m unittest discover -s tests -p "test_t3_bc_running_one_scale_triangle_v1.py" -v

# Fast: evaluate only existing B/C overnight UV running endpoints
python -u -X utf8 -m scripts.RunT3BCRunningOneScaleTriangle --all

# Full numerical mixed-quadratic scan: recommended *instead of* fast command
# if you want all scenarios without needing --resume.
& .\scripts\RunT3BCRunningOneScaleTriangle.ps1
```

**Do not run both fast and full commands into the same output directory**
unless the full invocation includes `-Resume`, or select a distinct `--output-id`.
For an already-completed full scan, resume with:

```powershell
& .\scripts\RunT3BCRunningOneScaleTriangle.ps1 -Resume
```

The default B/C numerical inputs come from the successful overnight runner:

`output/studies/t3_full_stack_overnight_v1/bc_full_mass/B_m1_uv_to_eft1_coupled_five_mass_v1.json`

`output/studies/t3_full_stack_overnight_v1/bc_full_mass/C_m1_uv_to_eft1_coupled_five_mass_v1.json`

Both include a grouped (not truly sequential) F boundary at 1e5 GeV.
The source SHA-256 inventory is checked, so changing RGBeta during Workstream A
may require rerunning the B/C numerical source with new beta exports.

## Outputs

`output/studies/t3_bc_running_one_scale_triangle_v1/`:

- `summary.json`, `progress.json`
- `B_m1_running_one_scale.json`, `C_m1_running_one_scale.json`
  (including per-column **running-input direct/sequential residuals**)
- with `--scan`: `bc_scenarios/{B_m1,C_m1}/{label}/*` original native
  B/C numeric trajectories, plus `scenarios/{B_m1,C_m1}/{label}.json`

The reports give the running `MF` Takagi singular masses, transformed y1/y2
flavour coupling columns, running scalar U/MX/MY, the three matched candidate
Matchete quartics, 3×3 candidate C5 and neutrino mass matrices, and singular
neutrino masses **explicitly marked NOT PHYSICAL**.

### Coefficient interpretation

The archived source is an *ordered* two-lepton coefficient A, with three
heavy columns. For each column r, its three distinct spurion coefficients
`(a11, a12, a22)` give

\[
 C_5^{\mathrm{candidate}}=\sum_{r=1}^3
 (\bar y_{1r},\bar y_{2r})
 \begin{pmatrix}2a_{11}^{(r)}&a_{12}^{(r)}\\
                 a_{12}^{(r)}&2a_{22}^{(r)}\end{pmatrix}
 (\bar y_{1r},\bar y_{2r})^T.
\]

For the complex symmetric running MF, the chosen heavy basis is the Takagi
basis `V^T MF V = diag(MF_r)`. The inspected original T3Topology
builds `y1` with `CConj[F]` but `y2` with `F`. The candidate spinor
orientation therefore uses `y1_mass = y1_UV conjugate(V)` and
`y2_mass = y2_UV V`, and **also explicitly computes** the alternative
`y1_UV V` result and its relative C5 difference. Because the full native
Matchete spinor-phase bridge is not certified, **neither selection is a
physically certified conversion**. Source v22 already
verified this exact ordered-source C5 convention against the v18 frozen
reference, and the **new runner independently replays that frozen parity**.

For the three Matchete scalar couplings, the **existing Stage17 *candidate***
dictionary is kept unchanged:

- `lambdaT3 (Matchete) = lambdaT3 (RGBeta)`
- `lambdaT3ExtraQ020020I1 = -lambdaHHdagS2S2 (RGBeta)`
- `lambdaT3ExtraQ020200I1 = -lambdaHHdagS1barS1bar (RGBeta)`

A source and numerical identity in the declared common basis is *not* the
physical cross-package SU(2)/spinor-phase certificate. Running-input parity
therefore tests matching **within the restricted source**, not the full
renormalized EFT equivalence. The other B/C quartics
enter the running masses/couplings but are absent in this restricted leading
`O(y^2 lambda)` `LLHH` source term.

### Mixed quadratic scenarios

Optional `--scan` repeats the *original* B/C coupled UV→EFT1 numerical
integration with (a) zero independent mixed quadratic, (b) phase reversed,
(c) 2.5× nominal mixed quadratic, bounded at `fraction<=0.45`. It evaluates
**every** scenario only at its common 1e5 GeV UV endpoint, never combining
Yukawas from one scenario with scalar masses of another.

A changed independent mixed bilinear changes the propagator mass basis.
Because the archived v22 source was generated for the restricted UV operator
inventory **without that independent bilinear**, scanning it through the
running mass eigenstates is a candidate extrapolation, **not a proof** that
this is the physical general-model matching.

## Remaining physical blockers (not addressed by this batch)

1. Workstream A/B common physically normalized Matchete↔RGBeta scalar,
   fermion, Yukawa and quartic dictionary.
2. Recompute/revalidate native matching with the independently nonzero
   general mixed quadratic in the full UV Lagrangian.
3. Correct EFT1/EFT2 renormalized matching, including finite field,
   parameter, counterterm, EOM and Wilson-operator effects.
4. True sequential decoupling of heavy fermions and both running scalar
   eigenstates, with matched EFT evolution between their mass scales.
5. Full 16-model physical C5/mν including the six other independent
   Yukawa matrix types where gauge symmetry allows them.

Thus: **a source-controlled, running-input one-scale leading coefficient has
been implemented, but `physical_C5_certified` remains `false`.**
