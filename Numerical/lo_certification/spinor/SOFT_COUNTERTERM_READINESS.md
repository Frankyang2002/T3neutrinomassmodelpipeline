# Workstream 1 — F-only matching: EFT1 soft / counterterm readiness

## Source and status

This batch **does not implement EFT1 soft diagrams, MS-bar counterterms or a physical C5**. It follows the existing immutable native F-only, raw LLHH, and pole/log outputs. It recomputes the upstream 176 Laurent triplets from the original source, checks all stored `trios.json` and four expression files per model **byte-for-byte**, and issues 176 source-resolved soft/CT requirements. The original 16-model Wilson release is untouched.

The existing `RunT3GeneralFThreshold.wl` and `RunGeneralT3FermionFirstEFT1OneLoop.wl` use Matchete `Match` to return a **matched EFT**. It is not independently established whether the saved native loop expression includes the hard/soft subtraction. **It must not be relabeled the unsubtracted UV amplitude or subtracted a second time by assumption.**

Each kernel record retains its source index map, exact signed native pole coefficient, logarithm `Log[μbar2/MF^2]`, and chirality. The unknown EFT1-soft amplitude and CT remain `null` / `false`. In a separately proven convention in which the input is an un-subtracted UV amplitude, the *conditional* pole balance would be

\[
R_{\rm UV}-R_{\rm EFT1,soft}+R_{\rm CT}=0.
\]

**This equation has not been applied to Matchete's already-matched output.** The output does not generate an invented counterterm `-R_native`.

## Actual numerical building block

The scalar master integral in a standard dimensionally regulated Passarino–Veltman convention is

\[
B_0(0;a,b)=\frac{1}{\bar\epsilon}+B_{0,\rm fin},\quad
B_{0,\rm fin}=1-\frac{a\log(a/\bar\mu^2)-b\log(b/\bar\mu^2)}{a-b},
\]

with `a,b` scalar mass-squared eigenvalues and the equal-mass limit `-Log(a/μbar²)`. The finite function is evaluated using a numerically stable `log1p` implementation. In this convention the log-scale derivative is `d B0_fin / d ln μbar² = +1`. The master has **no** included overall `i/(16π²)` factor, SU(2) tensor contraction, fermionic spinor phase, or model-specific vertex symmetry factor. Neither the `+1/epsilon-bar` normalization nor the finite result is by itself an EFT1 operator coefficient.

The two mixed-doublet runs (`B_-1`, `C_-1`) have *hashed* numerical mass-basis input files in the native batch. The script validates their SHA-256 against the original Matchete inputs and exports the three mass-eigenstate masters `B0(XX)`, `B0(XY)`, `B0(YY)` at `μ_F = 100000 GeV`. For the other 14 models it exports the symbolic `B0` formula and leaves scalar masses null, because the original native batch does not bind their masses numerically. The B/C mass inputs are the existing experimental matching inputs and are **not asserted to be synchronized with the 16-model neutrino benchmark**.

## Installation paths

- `Numerical/lo_certification/spinor/soft_counterterm.py`
- `Numerical/lo_certification/spinor/SOFT_COUNTERTERM_READINESS.md`
- `scripts/BuildLOSpinorSoftSubtractionLedger.py`
- `scripts/ProbeMatcheteMatchingSource.wl` (optional unexecuted native diagnostic)
- `tests/test_lo_spinor_soft_counterterm.py`
- `output/general_observable/lo_certification_spinor/github_soft_counterterm_reference_20261103.json`

## PowerShell (repository root)

```powershell
$env:PYTHONPATH = (Get-Location).Path
python -X utf8 -m unittest discover -s tests -p 'test_lo_spinor_*.py' -v
python -X utf8 scripts\BuildLOSpinorSoftSubtractionLedger.py --check-only
python -X utf8 scripts\BuildLOSpinorSoftSubtractionLedger.py --run-id soft01
```

Output root: `output/general_observable/lo_certification_spinor/soft_ct_e515d82d118284ba_soft01/`.
For each of 16 models: `pole_and_soft_ct_requirements.json` and `scalar_B0_master.json`; top-level `summary.json`. The directory must be new (immutable). No Matchete rerun is required.

Optional Matchete **source-level diagnostic**, not matching execution:

```powershell
wolframscript -file scripts\ProbeMatcheteMatchingSource.wl
```

It writes a uniquely named `matchete_match_semantics_<UUID>.json` inside `output/general_observable/lo_certification_spinor/`. Its text markers and `Match` downvalues are not, without independent interpretation, evidence of correct hard/soft subtraction.

## Remaining physics

1. Inspect Matchete's `Match`/one-loop subtraction semantics before adding another EFT1 soft loop (avoid double counting).
2. Calculate the actual EFT1 `LLSS × HHSS` scalar loop with **both scalars active** and correct CG/flavour contractions, plus any other required EFT1 diagrams.
3. Establish counterterms and `1/epsilon` cancellation in a **declared** renormalization scheme, including wavefunctions and coupling renormalization where needed.
4. Independently determine the native GammaCC to master Weyl normalization with signed connected amplitudes, including Dirac/Majorana cases.
5. Only after a verified physical fermion boundary: evolve EFT1 to `mu_S`, integrate out **both BSM scalar multiplets together**, then match to SMEFT.
