# T3 fermion-first EFT1: source-resolved full-SU(2) quartic response

**Scope:** New standalone EFT1-mixing workstream files only. No modification to matching, scalar-threshold, existing master-RGE or previously certified outputs.

## Current source state (read-only GitHub `main` at `a687cb45263a6f46cf44537861df79a7dd07a0df`, 2 Nov 2026)

- `output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/` has all 16 **provisional**, source-resolved fermion-tree LLSS boundaries. PL sources and the exact SU(2) CG/real-scalar contractions are retained, and the PR Hermitian conjugates have been excluded. The numerical sums have been checked, but the Matchete `GammaCC` → Weyl sign/factor **is unresolved**.
- `output/general_observable/lo_eft1_mixing/source_kernel_v1/` supplies a **conditional neutrino/H0R pivot** for each heavy scalar pair. This new module extends that calculation to *all* lepton and Higgs SU(2) components and checks the entire Weinberg tensor. It does not change or regenerate the original pivot stage.
- `output/general_observable/lo_fermion_boundary/matchete_gamma_cc_evidence_v2.txt` confirms Matchete was loaded and package sources inspected, but explicitly reports `EvidenceOnly_NoPhysicalSignInferred`. Do not silently set `N_GammaCC_to_Weyl=+1` or `-1`.
- No normalized, fully certified real-scalar EFT1 quartic dictionary is available to this workstream. The Stage-22 candidate Matchete/RGBeta tensor alignment is not the same as a certified **numeric** hhss vertex map at `mu_F`.

## What the implementation does

Input: the *real* Chat-1 provisional boundary JSON, including each separately labelled Matchete PL operator and all complex three-generation coefficients. The input coordinates are `(lepton_SU2_i, flavour_p, lepton_SU2_j, flavour_q, real_scalar_a, real_scalar_b)` with indices 1-based. Each independent coordinate is expanded over its **distinct** fermion and scalar exchange orbit, without extra factorials. Missing/modified certifications, inconsistent source sums and invalid indices fail closed. The source recomposition check uses **relative** tolerance scaled to the sizes of actual contributions (no absolute zero cutoff for suppressed channels).

The implemented scalar-pair master RGE is

\[
 16\pi^2 \mu\frac{d C_{(ip)(jq)AB}}{d\mu}\bigg|_{\lambda}
 =\sum_{a,b\in\text{heavy}} \lambda_{ABab} C_{(ip)(jq)ab}
 =\sum_{a\le b}(2-\delta_{ab})\lambda_{ABab} C_{(ip)(jq)ab}.
\]

Here `A,B = 1..4` are the four Higgs real coordinates, `a,b=5..n` contain **both** heavy scalar multiplets, and the full real-scalar potential convention is `V4 = lambda_abcd phi_a phi_b phi_c phi_d/4!`. This is the quartic term in the repository's `RGE/general/WilsonTensorRGE.py`. The other **six** terms of that master equation remain in the previous `GeneralEFT1C5Mixing.py`; this source-specific diagnostic does not declare them zero in an unverified EFT1 coupling inventory.

For the *symmetric* project convention

\[
 \mathcal L_5=\tfrac12\sum_{pq} C_{5,pq}(L_p\cdot H)(L_q\cdot H)+\mathrm{h.c.},\quad
 C_{5,pq}=C_{5,qp}, \qquad H=(H^+,H^0),
\]

and `H+ = (phi_1+i phi_2)/sqrt(2)`, `H0=(phi_3+i phi_4)/sqrt(2)`, let

\[
 z_{\nu,A}=\tfrac1{\sqrt2}(0,0,1,i),\qquad
 z_{e,A}=-\tfrac1{\sqrt2}(1,i,0,0),\qquad
 T_{ijAB}=\tfrac14(z_{iA}z_{jB}+z_{iB}z_{jA}).
\]

Then the canonical master tensor must have

\[
 C_{(ip)(jq)AB}=C_{5,pq}\,T_{ijAB},\quad
 T_{\nu\nu,H0_RH0_R}=1/4.
\]

We project the **complete** LLHH beta onto `T` using the complex Hermitian inner product and report the entire orthogonal residual. The neutrino/H0R pivot `4 beta_(nu_p nu_q, H0R H0R)` is independently reported. A nonzero residual means the projected number is **not** an independent physical Weinberg beta, even if the pivot is nonzero. A zero residual is a *conditional* SU(2) tensor-algebra check, **not** a Matchete normalization or quartic source certificate.

### Without a quartic input

By default, the output is a factorized list of complex `(3,3,2,2)` lepton-SU(2) coefficient blocks per Matchete PL source and per **unordered** heavy scalar pair. The real scalar pair multiplicity `1` or `2` is explicit. No beta or physical C5 is assigned.

### With a separately supplied diagnostic quartic

An optional `--quartic-root` supplies a **source-hash-bound, explicitly diagnostic** hhss real-component tensor. It may contain several **independently named quartic operators**, each with its own component list; even if two operators have the same hhss component, their values are **added** and kept separately. It reports every `(Matchete source, quartic operator, scalar pair)` beta projection and residual, every sum by source and quartic operator, and the total. Every complex `3×3` matrix is encoded as `{"real":[...],"imag":[...]}`.

A quartic input file should be placed at `output/general_observable/lo_eft1_mixing/quartic_diagnostics/A_m4/quartic.json`, etc., with content such as:

```json
{
  "format": "t3-eft1-diagnostic-hhss-quartic-v1",
  "model": "A_-4",
  "authority": "diagnostic_only",
  "source_boundary_sha256": "REPLACE_WITH_SHA256_OF_THE_EXACT_A_m4_BOUNDARY_JSON",
  "mu_F_gev": 100000.0,
  "coverage": "partial_diagnostic",
  "quartic_terms": [
    {
      "label": "test_lambda_1_NOT_PHYSICAL",
      "components": [
        {"indices": [3, 3, 5, 7], "value": 0.1}
      ]
    },
    {
      "label": "test_lambda_2_NOT_PHYSICAL",
      "components": [
        {"indices": [3, 3, 5, 7], "value": -0.02}
      ]
    }
  ]
}
```

All four indices are **1-based**, with exactly two Higgs and two heavy-scalar indices in each row. Values are real fully symmetric potential components, **not** unnormalized coupling-symbol coefficients. The sample coefficients above are arbitrary mathematical test values, **not** a source Matchete/RGBeta coupling dictionary and should not be interpreted as physical T3 input. A legacy single `quartic_components` list is also accepted as one unlabelled diagnostic channel; do **not** mix the two representations. The SHA and scale must match the boundary file.

For `0 < mu_S < mu_F` the optional fixed-order log is exported **per unresolved overall spinor factor**:

\[
 \Delta C_5^{\mathrm{EFT1,log}}/N_\Gamma
 =\frac{\ln(\mu_S/\mu_F)}{16\pi^2}
   [\beta_{C_5}^{(1)}/N_\Gamma].
\]

**Scale/scheme bookkeeping:** This is only the MSbar one-loop RGE logarithm implied by the master beta. A finite simultaneous-two-scalar MSbar threshold evaluated at the *same* `mu_S` must be added **once** by its owner. Its explicit scalar bubble logarithms must not be added here a second time. The matching-scale cancellations and source counterterms are not fully verified in this stage.

## Install/run (PowerShell, project root)

Extract `T3_EFT1_FullSU2_Next.zip` into the existing project (no `-Force` is needed because all files are new):

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_EFT1_FullSU2_Next.zip" -DestinationPath .
python -m unittest discover -s tests -p 'test_full_su2_source_response.py' -v
```

First inspect the real 16-model source inputs **without writing any files**:

```powershell
python scripts\RunFullSU2SourceResponse.py `
  --release-root 'output\general_observable\lo_fermion_boundary\releases\source_021dcc1e8fe8ba0f' `
  --output-root 'output\general_observable\lo_eft1_mixing\source_su2_response_v1' `
  --check-only
```

Then build the factorized actual-source SU(2) kernel **without assumed quartics**:

```powershell
python scripts\RunFullSU2SourceResponse.py `
  --release-root 'output\general_observable\lo_fermion_boundary\releases\source_021dcc1e8fe8ba0f' `
  --output-root 'output\general_observable\lo_eft1_mixing\source_su2_response_v1'
```

Expected: `output/general_observable/lo_eft1_mixing/source_su2_response_v1/summary.json` plus `A_m4/full_su2_response.json` through `E_p2/full_su2_response.json`. The scripts write only to the specified new workstream output root and refuse overwriting files. If the root already exists, choose a fresh `source_su2_response_v2` destination. Existing 16-model source boundaries are **read-only**.

For future *diagnostic-only* component quartics, use a **different** output root:

```powershell
python scripts\RunFullSU2SourceResponse.py `
  --release-root 'output\general_observable\lo_fermion_boundary\releases\source_021dcc1e8fe8ba0f' `
  --quartic-root 'output\general_observable\lo_eft1_mixing\quartic_diagnostics' `
  --mu-S-gev 1000 `
  --output-root 'output\general_observable\lo_eft1_mixing\source_su2_quartic_probe_v1'
```

Requires **all 16** `quartic.json` files when running without `--model`. Start with `--model A_-4` to investigate one case. A quartic tensor missing any genuine coupling must remain flagged as `partial_diagnostic` and is *not* equivalent to the physical EFT1 beta.

## Verified and unresolved

- **Locally tested:** 11 independent unit tests including exact Weinberg tensor factors, complex flavour, distinguish pivot from full-SU(2) projection, scalar pair multiplicities, independent master beta function agreement for a nontrivial complex case, distinct quartic labels on the same real component, suppression-preserving summation, strict input checks and CLI no-overwrite behavior.
- **From read-only GitHub:** source all-16 release has verified model inventory; prior `source_kernel_v1` results exist for all 16. This package's all-16 CLI must still be run on the local Windows checkout. No new 16-model numerical results are claimed before that run.
- **Unresolved:** actual Matchete `GammaCC` sign/factor; verified complete hhss quartic map; other dimension-five operator-class closure and one-loop hard fermion matching; complete physical UV+EFT1+simultaneous scalar threshold matching. `physical_C5_computed` **always** remains `false`.
