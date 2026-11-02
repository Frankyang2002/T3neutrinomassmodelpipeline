# Batch 4 — A_-4 dense fitted common input and joint scalar matching

**Status: CONDITIONAL / EXPLORATORY. Not a certified physical one-loop neutrino-mass prediction.**

## Scope

This new isolated workstream uses the *same frozen fitted A_-4 UV coupling state*
for the first-threshold running and LLSS calculation instead of the previous
illustrative diagonal Yukawas `y1=y2=4e-4 I`. It uses the existing files and
code, not an invented differential equation or Clebsch–Gordan coefficient.

1. Read `Numerical/benchmarks/fitted_reference/T3_dS1_1_dS2_3_dF_2_alpha_m4.json`
   and enforce SHA256 `3cb8dc3e94f2bd6f60571a7a85477638e8fc02109bccfcec0e04055852d81247`.
2. Require the frozen UV state to supply **all 20** couplings of the canonical
   `full16_general_v1` UV RGBeta export. No missing values, borrowed initial
   state, RG-interpolated Yukawas or implicit zeros are permitted. Original
   frozen SM `gY,g2,g3,lambdaH,yu,yd,ye` are included, not only BSM couplings.
3. Reuse `Numerical.orchestration.GeneralRunningBenchmark` for the two numerical
   ODE stages, `muUV=1e7` to `muF=1e5` to `muS=1e3` GeV. **All BSM fermions
   disappear together at muF**; EFT1 retains **both** BSM scalar multiplets.
4. Reuse the original native Matchete F-only tree Wilson *source* and exact CG:
   `output/general_observable/fermion_first_all16_v1/A_m4/data/eft1_after_F_wilson_seed.json`.
   The existing `prepare_heavy_flavour` rotates the dense Yukawas into the
   *evolved* heavy mass basis, and `flavour_matrix` computes the two **barred**
   Yukawa contractions. An independent einsum checks the contractions.
   `project_ordered_term` then evaluates the native `T3MixCG` into real scalar
   components. `operator_rows` applies the established lepton flavour order,
   canonical orbit and symmetry convention. PR Hermitian conjugates are not
   counted as extra neutrino-mass operators.
5. Write a **new provisional** F-tree LLSS JSON with the same adapter-compatible
   source grammar, different frozen benchmark, full off-diagonal flavour
   support and explicit unknown GammaCC normalization. Existing releases
   are never overwritten.
6. Reuse `build_source_channels` on the original Matchete scalar-quartic seed,
   then the existing `build_ordinary_background` and
   `evaluate_ordinary_response` at the **same** muS for both scalars. This is
   the previously implemented fixed-order LLSS-induced leading log plus one
   simultaneous finite scalar bubble. It does NOT add second-order LLSS running.
7. Output the complex 3x3 conditional `C5`, `mnu=v² C5/2`, three singular
   values in eV per unknown `|N_GammaCC|`, scalar mass eigenvalues, flavour
   support counts, source hashes, Ward/matching-scale residuals, and explicit
   physical-certification flags.

The source-derived scalar HHSS aliases are **candidate** Matchete/RGBeta
identifications, not a physical certified dictionary. The ordinary scalar
basis is appropriate to `A_-4`; B/C mixed and self-split scalar cases are **not**
supported by this runner. This batch does not invent or evaluate native
Matchete spinor signs, hard one-loop fermion matching, full EFT1 Wilson mixing,
EWSB-resummed mass shifts or PMNS fits.

## Package contents (project-relative paths)

```
Numerical/lo_certification/benchmark/common_input_batch4.py   (new)
Numerical/lo_certification/benchmark/README_BATCH4.md         (new)
scripts/RunLOCommonInputBatch4.py                              (new)
tests/test_lo_common_input_batch4.py                           (new)
output/general_observable/lo_certification_benchmark/
    batch4_pinned_UV_flavour_witness.json                      (new, static source arithmetic)
```

Only **new files** are supplied. No existing code, benchmark, frozen release,
matching summary, GitHub file or other workstream should be replaced.

## Windows PowerShell

Download `T3_Batch4_Common_Fitted_Input.zip` into Downloads. From the
repository root and activated `.venv`:

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = "$env:USERPROFILE\Downloads\T3_Batch4_Common_Fitted_Input.zip"
if (-not (Test-Path -LiteralPath $zip)) { throw "ZIP not found: $zip" }
Expand-Archive -LiteralPath $zip -DestinationPath $repo -Force
Set-Location $repo
$env:PYTHONPATH = (Get-Location).Path

python -X utf8 -m unittest discover -s tests -p 'test_lo_common_input_batch4.py' -v
python -X utf8 scripts\RunLOCommonInputBatch4.py --check-only
python -X utf8 scripts\RunLOCommonInputBatch4.py --run-id source01
```

The last command runs the existing native-source projector and RGEs locally.
The generated directory is:

```
output/general_observable/lo_certification_benchmark/batch4_common_input_source01/
    fitted_running.json                       # one full frozen UV state + muF + muS snapshots
    A_m4_fitted_LLSS_provisional.json          # new dense-flavour source-resolved LLSS
    numeric_backgrounds/A_m4.json              # joint scalar matching metadata
    numeric_backgrounds/A_m4.npz               # source HHSS channels and masses
    conditional_scalar_matching.json           # log and one joint scalar finite bubble
    summary.json                               # conditional C5/mnu and diagnostics
    report.md                                  # concise report
```

The output directory is **immutable**. If it exists, run with e.g.
`--run-id source02` (do not delete earlier runs). Missing native or canonical
symbolic files cause an explicit error, not substituted inputs.

The direct UV fitted-reference output is included as a **separate diagnostic**
when available. Its source is the frozen UV state, whereas the new conditional
EFT matching uses the RGE-evolved `muF` Yukawas and `muS` scalar parameters.
The two matrices therefore must **not** be compared as a scale-matched
physical equality and no best-fit GammaCC phase is extracted from them.

## Validated and not yet validated

- Verified from current GitHub: pinned fitted benchmark, full 20-coupling
  original UV export, original F-only Matchete tree/CG source, correct
  `A_-4` ordinary source/real-scalar topology, and existing numerical RGE,
  CG, LLSS orbit and joint matching functions.
- Locally tested for this package: pure source-contract, complex-matrix,
  conjugated heavy-flavour contraction, benchmark support, fail-closed and
  immutable ownership tests. Full Windows-repository integration is skipped
  when those source files are not present in the preparation environment.
- **Not yet executed in this environment:** the full local canonical ODE
  evolution and 1-TeV matching, which require your entire checkout.
  Do not report numerical `C5` values until the Windows run completes.

## EFT matching convention

`UV -> integrate out ALL heavy fermions at muF -> EFT1(S1,S2) -> integrate
out BOTH scalar multiplets at muS -> SMEFT`.
