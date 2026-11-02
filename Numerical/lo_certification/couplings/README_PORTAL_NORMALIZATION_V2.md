# Workstream 2 — exact portal normalization and full-16 HHSS inventory (v2)

This is an **additive** audit layered on the existing Workstream-2 v1/v2 batch; it does not edit shared `RGE/`, `Lagrangian/`, `Numerical/joint_scalar_matching/`, source-generated registries, prior reference results, or GitHub. All outputs are under `output/general_observable/lo_certification_couplings/`.

## Files in this ZIP (complete files, not patches)

- `Numerical/lo_certification/couplings/portal_normalization.py`: derives exact real-component HHSS portal tensors from source-declared RGBeta singlet/adjoint invariants, matches the full Matchete native component tensors by exact SymPy span solving, and converts numerical muF/muS portal couplings only when the map is square/full-rank and exact.
- `Numerical/lo_certification/couplings/quartic_audit.py`: checks the full current 16-model Matchete/release/cache and RGBeta source, enumerates the generated 24 additional HHSS candidate invariant entries, builds portal-only numerical arrays at both scales, and records independent blocking dependencies.
- `scripts/RunLOQuarticNormalization.py`: all-16 CLI runner.
- `tests/test_lo_couplings_normalization.py`, `tests/test_lo_couplings_normalization_integration.py`: analytic and full-repository regressions.
- `Numerical/lo_certification/couplings/README_PORTAL_NORMALIZATION_V2.md`: this guide.
- `output/general_observable/lo_certification_couplings/PORTAL_MAP_SCHEMA_V2.json`: machine-readable output contract (not a numerical input).

## Running in local PowerShell

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WS2_Portal_Normalization_V2.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_lo_couplings_*.py' -v
python -m scripts.RunLOQuarticNormalization --repo-root . --run-name ws2_portal_normalization_v1
Get-Content 'output\general_observable\lo_certification_couplings\ws2_portal_normalization_v1\summary.json' -Raw
```

Choose a **new** `--run-name` on subsequent runs. The runner refuses to overwrite an existing result directory. An optional repeatable `--model` argument selects a subset for debugging; default is all 16. To inspect exactly verified directions in PowerShell:

```powershell
$j = Get-Content 'output\general_observable\lo_certification_couplings\ws2_portal_normalization_v1\summary.json' -Raw | ConvertFrom-Json
$j.models | Format-Table model,S1_span_status,S2_span_status,RGBeta_portal_directions_missing_from_native_seed,extra_HHSS_candidate_count -AutoSize
```

### Physics and exact mathematical method

The current `RGE/running/rgbeta/T3RGBetaModel.wl` defines the **two** ordinary heavy-scalar portal directions for a nontrivial SU(2) multiplet:

\[
\mathcal L_{4,\mathrm{portal}} =
\sum_{i=1}^2\big[\lambda_{Hi}(H^\dagger H)(S_i^\dagger S_i)
+\mathbf{1}_{d_i>1}\lambda_{Hi\mathrm{Adj}}(H^\dagger t^A H)(S_i^\dagger T_i^A S_i)\big].
\]

For a singlet \(d_i=1\), only \(\lambda_{Hi}\) exists. With canonical real fields \(z_p=(x_{2p}+ix_{2p+1})/\sqrt2\), and real Hessian \(K(G)_{ab}=\partial^2(z^\dagger Gz)/\partial x_a\partial x_b\), each portal invariant has the unique quartic tensor

\[
\lambda_{ABab}^{(G,U)} =
\frac{\partial^4[-(H^\dagger G H)(S^\dagger U S)]}{\partial h_A\partial h_B\partial s_a\partial s_b}
=-K(G)_{AB}K(U)_{ab}.
\]

The **minus sign** is the source's \(V_4=-\mathcal L_{4,\mathrm{int}}\); derivatives enforce the required \(4!\) tensor convention. The symmetric tensor comparison includes *every* Higgs and scalar real-component coordinate, including off-diagonal components and exact irrational SU(2) CG factors.

For each scalar species, the native source tensors \(N_j\) are compared with independent RGBeta portal tensors \(R_i\) by solving

\[
N_j=\sum_i R_i A_{ij},\qquad
\lambda_\mathrm{RGBeta}=A\lambda_\mathrm{Matchete}.
\]

The solver checks **exact** SymPy equality for *every* component. It does not certify a map based only on a projection coefficient or a single neutral component. Only when \(A\) is square, full-rank, and exactly closes is the inverse used to derive numerical native portal couplings at the independently recorded \(\mu_F\) and \(\mu_S\). A separate numerical real-tensor reconstruction checks the result. Otherwise the corresponding native couplings remain `null` and the rank/index mismatch is a blocker. The assumption that both packages use the *same* standard complex SU(2) highest-weight component coordinates is explicitly recorded; a failure is not concealed by inventing a field rotation.

This adds an **independent basis-spanning check** absent from the earlier provisional alias dictionary. In particular, it will detect a native seed with only one portal direction for a non-singlet whose RGBeta potential has two independent directions. It does not assume any result in advance for individual models.

### Outputs at `output/general_observable/lo_certification_couplings/<run-name>/`

`summary.json` records 16 model statuses, total missing ordinary portal directions, the prior 24 additional HHSS candidate entries, source hashes, and blocking dependencies. Each `<slug>/` contains:

- `portal_normalization.json`: exact symbolic native-to-RGBeta basis conversion matrices per scalar, inverse coupling conversion *when proven*, numerical conversions at \(\mu_F\) and \(\mu_S\), source-hashed generated additional invariant inventory, and individual physics blockers.
- `rgbeta_portals_only_muF_muS.npz`: `muF_portal_channel_tensors`, `muF_portal_sum_tensor`, `muF_portal_labels`, with matching `muS_...` counterparts. A channel is created only when the corresponding numerical coupling is present at that scale. The metadata explicitly lists every omitted required coupling; absence is never treated as a physical zero. Arrays are in the assumed original interleaved RGBeta complex highest-weight component frame; **do not pass them into complete EFT1 Wilson or joint scalar solvers**.
- `failure.json` *only* when current local source integrity or a model-specific conversion check fails. The CLI returns nonzero if any model fails.

The generated extra quartic files `RGE/general/generated/<model>_extra_quartics.wl` are checked read-only against their registry entries. Additional channels with unique phase/intertwiner structure, e.g. the B/C mixed doublet entries and special split-real cases, are **not** assigned arbitrary real CG factors or spectral mass rotations.

### What remains blocked

1. Independent Matchete-to-RGBeta **topology mixing `lambdaT3`** normalization, complex phase and genuine representation intertwiners.
2. The full **24 additional HHSS registry candidate entries**, including independent channels absent from the numeric running state. Their counts are checked, but their physical fourth-derivative tensors and scheme-matched muF/muS coefficients remain unresolved.
3. The complete renormalizable SM Weyl Yukawa and gauge-component embeddings, signed Matchete `GammaCC` to master Weyl factor, and independent one-loop hard-F/operator closure.
4. The verified scalar frame transformation in six split-real and B/C mixed cases, including a running \(\mu_F\to\mu_S\) basis transformation where needed.

**No physical `C5` is produced**. The established fixed hierarchy remains exactly:

`UV -> integrate out ALL heavy fermions together at muF -> EFT1(H,S1,S2) -> integrate out BOTH heavy scalar multiplets together at muS -> SMEFT`.

As in the prior workstream, no sequential heavy-fermion or individual-scalar decoupling is introduced.
