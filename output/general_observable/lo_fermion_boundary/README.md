# Fermion-first T3 LLSS numerical boundary — fail-closed handoff

**Intended repository directory:** `output/general_observable/lo_fermion_boundary/`.

## What is implemented

`build_boundary.py` reads three *existing* local, previously verified reference stages:

- `output/general_observable/ordered_wilson_tensor_v4/`
- `output/general_observable/general_wilson_polynomial_batch_v2/`
- `output/general_observable/general_wilson_normalization_audit_v1/`

It validates all 16 models and matches the numerical complex three-family LLSS tensor back to every independent **Matchete PL source operator**. The mixed `B_-1`/`C_-1` model source inventories must have exactly **10/12 PL plus 10/12 PR** terms, respectively, with native `Proj[-1]` / `Proj[1]` chirality. It never adds the PR conjugate to the PL coefficient. The source-order flavour permutations, exact SU(2) CG projections already present in the verified real-component tensors, and the 1/2/4 distinct permutation-orbit divisions are preserved. Source-resolved coefficients are summed and compared numerically against the existing polynomial bridges. Input hashes and the saved physical-normalization audits must agree.

**It does not regenerate Matchete data** or write to existing shared stages. If the local `v4` / `v2` / `v1` source artifacts are absent or differ, the script refuses to synthesize them from historical sources. In particular, it does not read the wrong-chirality `v3` mixed tensors.

## The unresolved physical issue

The reference source and audits establish a *formal commuting-spinor polynomial* with basis coordinates

`[lepton_SU2_i, flavour_p, lepton_SU2_j, flavour_q, real_scalar_a, real_scalar_b]`,

all 1-based, in GeV^-1. The independent representative expands over the distinct exchange orbit under `(i,p)<->(j,q)` and `a<->b`. **No extra identical-field orbit multiplicity should be applied to this formal polynomial.**

They do **not** independently demonstrate the absolute four-component Matchete `GammaCC` to master-RGE left-handed Weyl-spinor bilinear mapping. That mapping needs an actual signed amplitude / spinor convention calculation, ideally with an independent matched amplitude. We therefore leave `gamma_cc_to_master_weyl_factor` and `gamma_cc_to_master_weyl_sign` as JSON `null`, never `1` or `+`. A genuine normalized physical boundary **has not been certified**.

The original all-16 F-first output also explicitly sets `raw_one_loop_F_matching_included: false`. Only **tree-level F-matched LLSS operators** (which subsequently produce the leading one-loop neutrino-mass mechanism through EFT1 scalar propagation/mixing) are included by this tool. It does not establish absence or inclusion of suppressed hard one-loop matching contributions and must **not** be presented as complete one-loop UV matching. Neither EFT1 RGE nor the scalar threshold at mu_S is implemented here.

## Install and run (PowerShell)

From `C:\thesistemp\T3neutrinomassmodelpipeline\` after extracting the ZIP, run:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 output\general_observable\lo_fermion_boundary\test_boundary.py
python -X utf8 output\general_observable\lo_fermion_boundary\test_boundary.py --live
python -X utf8 output\general_observable\lo_fermion_boundary\build_boundary.py --check-only
python -X utf8 output\general_observable\lo_fermion_boundary\build_boundary.py
```

The final command writes a **new** immutable release under

`output/general_observable/lo_fermion_boundary/releases/source_<input_sha_prefix>/`.

Inside each of the 16 model folders (`A_m4`, `B_m1`, etc.) a `boundary.json` contains the native source term, source index, full complex three-generation symmetric LLSS coordinates, scalar-basis metadata and a verified summed relative tensor. `summary.json` lists the 16 model results and explicit certification blockers. The script **refuses to overwrite** a prior release; use `--release-id new_unique_identifier` for a new run. All output is contained within the one allowed directory.

Downstream EFT1-mixing/scalar-matching workers may inspect these values as *uncalibrated numerical LLSS data*. They must not accept these as physically normalized master Wilson coefficients until the GammaCC-to-Weyl sign/factor is independently certified. The handoff schema explicitly flags this.

## Prerequisites / limitation

All three reference directories must exist in your **current local checkout** and correspond to the same matching, polynomial, normalization audit and mu_F numerical run. The GitHub copy inspected when preparing this package reports the 16/16 structural audit, but it is read-only and may lag your local files. This package was syntax-checked and its offline unit tests were run in an isolated environment; the real 16-model test and numerical generation must be performed from your Windows checkout, since your `C:` drive and saved Matchete JSON are not mounted here.
