# T3 EFT1: direct one-loop Weinberg-source mixing

**Scope:** newly introduced files only. Read-only inspection base: `Frankyang2002/T3neutrinomassmodelpipeline`, default branch `main`, 2026-11-02. No GitHub writes. No edits to matching, threshold or pre-existing RGE modules.

## Physics conventions verified in repository

- `RGE/general/WilsonTensorRGE.py` implements the seven terms of arXiv:2512.16890v1 Eq. (4.85), using **symmetric** Weyl `C[i,j,a,b]` and real-scalar components. This new NumPy evaluator reproduces those seven terms, rather than modifying that source.
- `RGE/matching/WeinbergTensorAdapter.py` defines a polynomial `kappa*(nu*H0 - e*H+)**2`, uses `H=(Re+i Im)/sqrt(2)`, and divides four polynomial derivatives by 4 to build `C[i,j,a,b]`. In this repository the physical symmetric Weinberg convention is `L5 = (1/2) C5[p,q] (L_p.H)(L_q.H) + h.c.` and `m_nu = -v**2*C5/2` (`physics/NeutrinoMass.py` and `Lagrangian/RunMatching.wl`). **Hence `C[nu_p,nu_q,ReH0,ReH0]=C5[p,q]/4`.** The tensor is transpose-symmetric in flavour, without complex conjugation of the off-diagonals.
- `Numerical/running/GeneralT3OrderedTensor.py` and `GeneralT3WilsonNormalizationAudit.py` label the existing source-factorised, ordered Matchete outputs **not canonical**, and explicitly leave the `GammaCC` to two-component Weyl sign/factor unresolved. `output/general_observable/ordered_wilson_tensor_v4/` is therefore NOT a physical input for this code.

## The retained master RGE

With `dot(C) = 16*pi**2 * dC/dln(mu)`, `i,j,k,l` Weyl components and `a,b,c,d` real scalar components, the contributions are (summed internal indices and gauge-generator labels implicit):

```text
1  yukawa_wavefunction:
       y_ijd y*_{kld} C_klab
2  scalar_pair:
      -(1/2) sum_{swap a,b} [2 g^2 theta_ac theta_bd - lambda_abcd] C_ijcd
3  mixed_yukawa_gauge:
       sum_{swap a,b; swap i,j} [2 y_jkc y*_{klb} + y_jkb y*_{klc}
                                + 4 g^2 t_lj theta_bc] C_ilac
4  crossed_yukawa:
       sum_{simultaneous swap (a,b),(i,j)} y_jla y_ikc C_lkbc
5  conjugate_coefficient_yukawa:
       [2 y_ild y_jkd + y_ijd y_kld] C*_{klab}
6  scalar_anomalous_dimension:
       sum_{swap a,b} gamma_s(b,d) C_ijad
7  fermion_anomalous_dimension:
       sum_{swap i,j} gamma_f(j,k) C_ikab
```

Here `lambda_abcd = d^4 V/(dphi_a dphi_b dphi_c dphi_d)` for `V_4 = lambda_abcd phi_a phi_b phi_c phi_d / 4!`. `gamma_s` and `gamma_f` are computed in the *same* sign/generator convention as `RGE/general/AnomalousDimensions.py` (not guessed from a simplified model). All gauge, quartic and Yukawa components are read explicitly and hermitian generator/symmetric-tensor checks are applied. No suppression-based filtering is applied to any operator or RGE term. Individual terms may be non-Weinberg-shaped; the **total sum over all sources** must project into the complete `LLHH` Weinberg subspace before an output is labeled physically certified (individual operators also report their Ward residuals).

For a source that is strictly `LL S_i S_j` (both scalars odd under an exact Z2 and Higgs even), all the non-quartic direct `LLSS -> LLHH` contractions vanish **if** the EFT1 scalar kinetic/Yukawa/coupling tensors respect that block structure. The scalar-pair quartic term can receive contributions from *every* allowed scalar invariant and representation component, including suppressed ones. The evaluator nevertheless retains all seven terms; no such zeros are hard-coded. Other fermion-first operator classes must be independently audited before physical certification.

## Fixed overall one-loop running and threshold handoff

At leading overall one-loop order, with `C^(0)` the tree-level matched and canonical dimension-five tensor at `mu_F`,

\[
\beta_{C_5}^{(1)} = 4 \;\left[\beta_{C_{ijab}}^{(1)}\right]_{i=\nu_p,j=\nu_q,a=b=\operatorname{Re}H^0},
\qquad
\Delta C_5(\mu_S) = \frac{\beta_{C_5}^{(1)}(\mu_F)}{16\pi^2}\ln\frac{\mu_S}{\mu_F}.
\]

It is mathematically sufficient at this fixed order to freeze `C^(0)` and renormalizable couplings at `mu_F`; their one-loop changes inside the mixing term are order `(16*pi**2)^-2`. A 21-point analytic-log trajectory is provided for numerical plotting/diagnostics. **This does not perform LLSS self-running**. There is no scalar-by-scalar decoupling, and no intermediate scalar threshold: both heavy scalars belong to EFT1 until the common `mu_S`.

The output records `beta_C5_16pi2`, `delta_C5_log`, complex flavour matrices per *operator source × seven RGE terms*, their sum, and the logarithmic scale derivative. Scheme is **MSbar**. The finite **joint two-scalar threshold matching** is *not* calculated by this workstream. When adding it, the scalar-threshold provider must subtract the EFT1 infrared logarithm **once** and use the same `C5` and `C_ijab` conventions; do not add an un-subtracted hard logarithm on top of this evolution. Any direct hard contribution already matched at `mu_F` must be tracked as a **separate** one-loop boundary, not assigned to this log. A scheme change requires an explicit finite counterterm. The scalar mass/splitting dependence beyond that in the renormalizable boundary is in the finite joint-scalar matching, not in the UV mixing anomalous dimension.

## Certified-boundary JSON format

A fully normalized input is explicitly identified by:

```json
{
  "format": "t3-eft1-certified-master-weyl-v1",
  "authority": "certified",
  "model": "A_-4",
  "scheme": "MSbar",
  "mu_F_gev": 100000.0,
  "mu_S_gev": 2000.0,
  "n_fermions": 6,
  "n_real_scalars": 12,
  "lepton_components": [[1,2],[3,4],[5,6]],
  "higgs_real_components": [1,2,3,4],
  "scalar_field_blocks": [
    {"name":"H", "role":"Higgs", "indices":[1,2,3,4]},
    {"name":"S1", "role":"heavy", "indices":[5,6]},
    {"name":"S2", "role":"heavy", "indices":[7,8,9,10,11,12]}
  ],
  "master_lagrangian_convention": "C_ijab psi_i psi_j phi_a phi_b + hc",
  "weinberg_lagrangian_convention": "(1/2) C5_ij (L_i.H)(L_j.H) + hc",
  "certificates": {"... nine required explicit boolean keys ...": true},
  "yukawa_components": [],
  "gauge_sectors": [],
  "quartic_components": [
    {"indices": [3,3,5,7], "value": 0.01}
  ],
  "operators": [
    {"source_id": "F_tree_1", "components": [
      {"indices": [1,1,5,7], "value": {"real": 1e-11, "imag": 2e-12}}
    ]}
  ]
}
```

**This is a shape illustration, NOT a certified physical input.** A real `certified` object MUST have all nine true certificates and explicit, complete nonempty gauge-generator sectors (SU(2), U(1)). The adapter requires the future Chat 1 canonical Matchete-to-Weyl output in `t3-chat1-canonical-master-weyl-boundary-v1` and an independently certified EFT1 quartic/Yukawa/gauge context in `t3-eft1-renormalizable-couplings-MSbar-v1`. There is no automatic factor fitting or inferred coupling map. Every independent `C[i,j,a,b]` is specified directly in master convention; symmetric i/j, a/b partners are automatically populated, **without adding symmetry factors**. For the quartic, every permutation of a real-scalar index quadruple is populated. For the Yukawa, `y[i,j,a]=y[j,i,a]` is enforced.

Generator entry format is `{"fermion": [[...],...], "scalar": [[...],...]}`, or each matrix as `{"real": [...], "imag": [...]}`; each sector has `{"coupling": 0.65, "generators": [...]}`. Tensor entries have 1-based indices, as in the repository. The synthetic fixture generator deliberately supplies zero gauge/Yukawa sectors and fake quartics, and cannot pass the certified gate.

## Source/verification distinctions

| Aspect | Status |
|---|---|
| 7-term one-loop master `psi^2 phi^2` beta evaluator | Implemented, component-sum tested |
| Complex 3×3 flavour and source/term decomposition | Implemented, tested |
| Higgs/Weinberg normalization `C5=4 C_nunu,R0R0` | Checked against repository adapter; tested with exact doublet quartic |
| All-16 schema/numerical run | Synthetic only; no physical parameter scan |
| Canonical normalization of CURRENT ordered Matchete boundary | **Not certified** in checked source |
| Full physical RGBeta-to-Matchete quartic / Yukawa dictionary | **Not supplied**, must be certified by responsible workstream |
| General EFT1 operator-class completeness beyond `psi² phi²` | **Not established**; required external audit |
| Finite simultaneous-scalar matching | **Not implemented** (separate owner's task) |
| Complete leading-one-loop neutrino masses | **Not yet available** |

Source files inspected include the repository `RGE/general/{WilsonTensorRGE,AnomalousDimensions,ScalarBasis,GaugeGenerators}.py`, `RGE/matching/WeinbergTensorAdapter.py`, `physics/NeutrinoMass.py`, `Numerical/running/{GeneralT3OrderedTensor,GeneralT3WilsonEvolution,GeneralT3WilsonNormalizationAudit}.py`, plus the older scalar-only intermediate-EFT modules. References: arXiv:2512.16890v1 Eq. (4.85), arXiv:1507.06314v3 and arXiv:0601225v1 for the scotogenic limit.

## Exact local installation/run commands (PowerShell)

Save the ZIP as `C:\Users\<your-user>\Downloads\T3_EFT1_C5_Mixing.zip`, then run:

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = Join-Path $HOME 'Downloads\T3_EFT1_C5_Mixing.zip'
Expand-Archive -LiteralPath $zip -DestinationPath .
python -m unittest discover -s tests -p 'test_general_eft1_c5_mixing.py' -v

# Generate ALL 16 synthetic-only input files (no Matchete normalization assumed).
python scripts\BuildSyntheticGeneralEFT1Fixtures.py `
  --output-root 'output\general_observable\lo_eft1_mixing\synthetic_input'

# Compute ALL 16 synthetic-only complex C5 one-loop logarithms.
python scripts\RunGeneralEFT1C5Mixing.py `
  --boundary-root 'output\general_observable\lo_eft1_mixing\synthetic_input' `
  --output-root 'output\general_observable\lo_eft1_mixing\synthetic_results' `
  --allow-synthetic
```

The synthetic run writes `synthetic_results\summary.json` and, for each model such as `A_m4`, `synthetic_results\A_m4\c5_eft1_log.json`. Its status is `SyntheticDiagnostic_NotPhysicalC5`. Repeating a command that targets existing files **raises `FileExistsError` rather than overwriting**; use a new output directory for subsequent studies.

After Chat 1 supplies **certified** canonical master-Weyl boundary tensors, AND a complete independently certified physical EFT1 coupling/gauge context is available, adapt one model with:

```powershell
python scripts\AdaptCertifiedEFT1Boundary.py `
  --chat1 'PATH_TO_CHAT1_FINAL_A_m4.json' `
  --eft1-couplings 'PATH_TO_EFT1_COUPLINGS_A_m4.json' `
  --output 'output\general_observable\lo_eft1_mixing\certified_boundaries\A_m4\boundary.json'
```

For all 16, place one resulting `boundary.json` under each of the tags `A_m4,A_m2,A_p0,B_m3,B_m1,B_p1,C_m3,C_m1,C_p1,D_m2,D_p0,D_p2,E_m4,E_m2,E_p0,E_p2`, and run:

```powershell
python scripts\RunGeneralEFT1C5Mixing.py `
  --boundary-root 'output\general_observable\lo_eft1_mixing\certified_boundaries' `
  --output-root 'output\general_observable\lo_eft1_mixing\certified_results'
```

There is intentionally **no `--allow-synthetic`** on the certified command. A physical log may be combined with finite joint-scalar matching only after the provider fixes compatible normalization and logarithm-subtraction conditions. No automatic finite matching or neutrino mass conversion is performed here.
