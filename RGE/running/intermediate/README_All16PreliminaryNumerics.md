# All-16 source-backed numerical comparison — conditional λT3-only batch

**Scope:** Read-only GitHub `main` inspected on 2026-11-02, HEAD `7f998f2ae4afb38e3864be07c4f628079b00fbc7`. This batch is **new files only** and does not overwrite earlier outputs or alter any original matching/RGE/mass-basis source code.

## Installed project paths (four NEW files)

1. `RGE/running/intermediate/All16PreliminaryNumerics.py` — source-backed all-model numerical adapter, provenance checks, one common joint scalar threshold, full SU(2) residuals, per-PL-source decomposition.
2. `scripts/RunAll16PreliminaryNumerics.py` — 16-model CLI, fail-closed input checks, strict output ownership/no overwrites.
3. `tests/test_all16_preliminary_numerics.py` — 12 unit tests covering basis rotations, scalar dimensions, provenance, cancellation/Weinberg interpretation, CLI and output ownership.
4. `RGE/running/intermediate/README_All16PreliminaryNumerics.md` — this document.

Existing dependencies (not replaced):

- `RGE/running/intermediate/B1PreliminaryNumerics.py` (previous successful B−1 numerical pilot)
- `Numerical/joint_scalar_matching/{core,fermion_boundary_bridge,source_quartic_v4,adapter}.py`
- `output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/` (16 formal F-tree LLSS coefficients)
- `output/full16_symbolic/` (Matchete scalar quartic seeds)
- `output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels/` (16 native quartic exports)
- `output/general_observable/all16_v1/{slug}/running.json` (illustrative EFT1 scalar running; used in 8 ordinary cases)
- `output/general_observable/mixed_f_synchronised_v1/{B_m1,C_m1}/` and `output/general_t3_mass_basis/{B_m1,C_m1}_final_mass_basis.json` (two verified *algebraic*, not physically certified, mixed-doublet mass basis handoffs)

## What the batch can and cannot compute

The pipeline always retains **both** scalars in EFT1 until one common scalar threshold:

\[
\text{UV} \longrightarrow \text{integrate ALL fermions at }\mu_F
\longrightarrow \text{EFT1}(H,S_1,S_2)
\longrightarrow \text{integrate ALL scalars at }\mu_S
\longrightarrow \text{SMEFT}.
\]

For a **single native `lambdaT3` source vertex** with all other quartic channels **excluded (NOT set to physical zero)**, construct the formal unsummed complex scalar source Wilson tensor and its one-channel response:

\[
\widehat\beta_{LLHH}^{(\lambda T3)}
=\sum_{a,b} C_{LLab}^{\rm F-tree,formal}(\mu_F)
             \lambda_{HHab}^{\lambda T3}(\mu_S),
\quad
C_{5,\log}^{\rm master}/N_\Gamma
=\frac{\ln(\mu_S/\mu_F)}{16\pi^2}\,\Pi_W^{\rm master}[\widehat\beta],
\]

\[
C_{5,\rm bubble}^{\rm master}/N_\Gamma
= 2\,\Pi_W^{\rm joint}\left[-\frac{1}{2(16\pi^2)}
  \sum_{a,b}\widehat C_{LLab}\,\widehat\lambda_{HHab}\,
  \overline B_0(0;m_a^2,m_b^2;\mu_S^2)\right].
\]

The master coefficient is **twice** the joint matching module's projected coefficient; the implementation carries that factor in both beta/log/finite parts. It diagonalizes the **whole heavy real scalar mass matrix at once**. It verifies full SU(2) Weinberg geometry and fixed-coupling cancellation of the explicit matching-scale dependence. When the Ward residual fails, no mass singular values are reported. PL native source sums are audited and PR conjugates are not counted twice. Neutrino mass singular values use `v=246.22 GeV` and `mnu = v² C5_master/2`; they are **per unit unknown |N_GammaCC| only**.

### Why numerical coverage is 10 candidate models, NOT all 16

| Model class | Models | Input treatment |
|---|---:|---|
| Original scalar real basis matches LLSS source | **8** | Source-derived quartic `lambdaT3` channel; masses and couplings at muS from **illustrative** `all16_v1` running records. These are **not independently certified as the same fermion-boundary benchmark**. |
| Mixed B−1/C−1 | **2** | Exact source-recorded complex-epsilon 8-real-component basis rotations at F and S, matched F Matchete parameter provenance, S eigenmasses and recorded candidate lambdaT3. |
| Split-real self-bilinear | **6** | **Blocked**. Matchete original complex S basis ≠ released minus/plus real LLSS basis. Native phases/CG map and F/S mass basis have not been independently established. No identity/permutation/reinterpretation of unverified spherical fields is allowed. |

For B−1 this extends, rather than overwrites, the prior numerical pilot. Normalized physical effects may be significantly different from the one-channel conditional answers. **Do not compare the conditional masses as physical predictions or as fitted oscillation data.**

The B/C masses and scalar `lambdaT3` were obtained from one later scale; the 8 ordinary records come from an independent illustrative running workstream. No model is certified to have all scalar couplings, full 7-term EFT1 mixing, hard F-loop matching, or signed `GammaCC` conventions. This is an incomplete first one-loop sector, not the full physical SMEFT threshold.

## Commands (PowerShell, local repository root)

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_All16_Preliminary_OneChannel_Batch.zip" -DestinationPath .

python -m unittest discover -s tests -p 'test_all16_preliminary_numerics.py' -v

# Validate the actual source data without writing any new output.
python scripts\RunAll16PreliminaryNumerics.py `
  --output-root 'output\general_observable\lo_eft1_mixing\all16_lambdaT3_preliminary_v1' `
  --check-only

# First complete attempt. Creates one summary and per-model detail files
# for models whose actual inputs were compatible.
python scripts\RunAll16PreliminaryNumerics.py `
  --output-root 'output\general_observable\lo_eft1_mixing\all16_lambdaT3_preliminary_v1'

Get-Content 'output\general_observable\lo_eft1_mixing\all16_lambdaT3_preliminary_v1\summary.json' -Raw
```

If you only want a C−1 result or one ordinary model, use `--model C_-1` or `--model E_-4` with a **fresh** output subdirectory. `--require-all` exits 2 unless every selected model has a Ward-valid one-channel result; absent split-real basis evidence makes `--require-all` fail as it should. `InputError` returns exit 2 regardless of this option; the known six split cases alone do not force a nonzero exit status.

Expected output layout:

```text
output/general_observable/lo_eft1_mixing/all16_lambdaT3_preliminary_v1/
    summary.json
    A_m4/lambdaT3_one_channel.json
    B_m1/lambdaT3_one_channel.json
    ... [other model files only where inputs are valid]
```

Each model result includes source and Matchete quartic seed SHA256, muF/muS, provenance of scalar input, source operator count, all native quartic labels, the **only included** `lambdaT3` channel, a per-PL-source beta breakdown, quartic log, finite joint bubble, full SU(2) residual, matching-scale cancellation residual and conditional mass singular values only when Ward-valid. Known split-real blockers appear in `summary.json`, not as fake zeros. Nothing writes outside `output/general_observable/lo_eft1_mixing/<new_run>/`.

## Next physics work if this batch succeeds

For an all-16 conditional table, complete the **six split-real source-matched basis transformations** from installed Matchete's actual CG/field definitions and use their correct muS quadratic spectra, then reassess. For genuinely physical results, certify native Matchete/RGBeta quartic normalization and signed `GammaCC` -> Weyl factor, plus complete EFT1 operator and one-loop matching. The current B−1 one-channel pilot remains an independent regression reference, not a physical normalization certificate.
