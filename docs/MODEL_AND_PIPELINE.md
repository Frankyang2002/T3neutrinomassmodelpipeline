# T3 model and calculation contract

## Symmetries and field content

The ordinary T3 model extends the Standard Model by heavy fermion multiplets `F` and two scalar multiplets `S1` and `S2`, odd under the imposed Z2. Standard Model fields are even. The package's hypercharge convention is

\[
Q=T_3+Y,\qquad Y(S_1)=\alpha/2,\quad Y(S_2)=(\alpha+2)/2,\quad Y(F)=(\alpha+1)/2.
\]

Some source literature uses `Y_literature=2Y`. The `alpha` parameter and all charge checks must use **one** convention within a calculation. The ordinary SU(2) representation classes in `common/T3Model.py`/`Lagrangian/model/T3ModelCatalog.wl` are:

| Class | `(dim S1, dim S2, dim F)` |
|---|---|
| A | `(1,3,2)` |
| B | `(2,2,1)` |
| C | `(2,2,3)` |
| D | `(3,1,2)` |
| E | `(3,3,2)` |

Across allowed hypercharges, the project uses a 16-model ordinary catalogue. A shared-scalar/scotogenic model has one **physical** scalar: introducing two formal topology roles must not produce two independent physical scalar multiplets.

## Definition of the general Lagrangian

Start from all Lorentz-, gauge- and Z2-invariant renormalisable operators for the chosen field representations. Schematically,

\[
\mathcal L_{\mathrm{UV}}=\mathcal L_{\mathrm{SM}}+\mathcal L_{\mathrm{kin,mass}}^{\mathrm{BSM}}
+\mathcal L_{\mathrm{Yukawa}}^{\mathrm{allowed}}-V_{\mathrm{scalar}}^{\mathrm{allowed}}.
\]

This equation specifies the **inclusion rule**, not an explicit sum over Clebsch-Gordan tensors. In code, the representation-dependent CGs, conjugations, independent contractions, flavour indices, and Hermitian conjugates are generated rather than replaced by a schematic topology vertex.

The operator basis must include, wherever symmetry allows:

- The two ordinary T3 Yukawa topology vertices **and** independent conjugate-scalar Yukawas. The latter are present for a hypercharge-zero scalar, subject to the model's field and flavour content.
- General permitted heavy-fermion mass structures, with correct Dirac/Majorana representation constraints and flavour indices.
- Higgs, scalar self and portal quartics, with *every independent* SU(2) invariant tensor (not merely one product of norms).
- The T3 Higgs-scalar mixing quartic and additional holomorphic/non-Hermitian quartics, with independent complex coefficients as appropriate.
- Every additional symmetry-allowed scalar bilinear, including the off-diagonal B(-1)/C(-1) scalar mixing when permitted.

The existing saved 16-model **structural** inventory records 50 additional complex-quartic model-slots, eight additional complex-bilinear model-slots, and six additional Yukawa model-slots (sums over different models, **not** interactions in one model). Their count is an operator-inventory cross-check, not a certification of all physical normalisations and flavour structures. Sources for the stored inventory: `docs/general_scalar_trial_operator_manifest.json`, `RGE/general/generated/general_rgbeta_quartic_registry.json`, `scripts/AuditFull16GeneralSymbolicPhysics.py` and the actual model-builder files.

A benchmark may set some allowed parameters to zero, but their operators remain part of the UV definition. Such boundary values should be listed alongside the matching scale, and RG evolution must use the relevant general-theory beta functions. In particular, zero extra bilinears do **not** follow from Z2 symmetry alone.

## Matching and running

The required hierarchy is

\[
\mathcal L_{\rm UV}\ \xrightarrow{\;\text{integrate out all }F\;}
\mathcal L_{\rm EFT,scalars}\ \xrightarrow{\;\text{integrate out all }S\;}
\mathcal L_{\rm SMEFT}.
\]

1. Run the complete UV couplings while all heavy fields are dynamical. Integrate out all heavy fermions before any scalar. Flavour-dependent fermion masses may produce multiple fermion thresholds.
2. Retain the scalar degrees of freedom and every operator needed at the chosen EFT dimension and loop order. Use matching conditions obtained from the same general Lagrangian and compatible field/mass bases.
3. Run the scalar EFT between fermionic and scalar thresholds. Integrate out the complete scalar sector, respecting mass mixing, operator inheritance, hard contributions and RG-induced terms without double counting.
4. Below the last scalar threshold, run SMEFT Wilson coefficients and form the light neutrino mass matrix. In the repository's documented Weinberg-operator convention, `m_nu = -(v^2/2) C5`; the Wilson coefficient must first be symmetrised and expressed in the consistent charged-lepton flavour basis.

The one-loop running convention is `16 pi^2 dX/d(ln mu) = beta_X^(1)` when the exporter supplies the one-loop numerator. The `16 pi^2` factor must be applied exactly once.

A complete RG/EFT prediction requires the renormalisable and nonrenormalisable operators relevant at the specified order, not just the terms appearing explicitly in a single neutrino-mass diagram. B(-1) and C(-1) can have additional independent flavour contributions to C5. A mass-mixing term can change the **leading** one-loop answer through scalar eigenstates and is not automatically higher-loop suppressed.

## Source ownership and what is established

| Area | Source location | Verified source-level observation |
|---|---|---|
| Model definitions | `common/T3Model.py`, `Lagrangian/model/` | Ordinary classes and representation tests |
| General UV interaction candidates | `Lagrangian/interactions/`, `Lagrangian/model/LagrangianBuilder.wl` | General scalar/extra Yukawa candidates exist, but some are opt-in |
| Matching | `Lagrangian/RunMatching.wl`, general matching runners | Distinct matching routes and basis checks exist |
| EFT ordering | `common/Thresholds.py`, `common/PipelinePlan.py` | Plan structures support multiple orders; execution support is narrower |
| RGEs | `RGE/general/`, `RGE/running/` | General templates/exports coexist with restricted paths |
| Numerics | `Numerical/core/`, `Numerical/running/`, `Numerical/fitting/` | Extra couplings are conditionally included in numerical states |
| Observables | `physics/` | Low-energy mass and mixing extraction |

In the inspected source, `Lagrangian/interactions/T3Topology.wl` uses the environment variables `T3_EXPERIMENTAL_GENERAL_QUARTICS`, `T3_EXPERIMENTAL_GENERAL_QUADRATICS` and `T3_EXPERIMENTAL_GENERAL_FERMIONS`. `Lagrangian/model/LagrangianBuilder.wl` itself returns `FullSymmetryCompletenessCertified -> False` and `FermionFlavorMassGeneralityCertified -> False`. `pipeline.py` also retains `--restricted-production` alongside general-symbolic operation. These are **actual source conditions**, not evidence that a fully general end-to-end computation has been demonstrated.

Before treating an output as a result of the general Lagrangian, verify all of the following against the generated files: the model's full allowed-interaction list; fermion flavour and scalar mass bases; coefficient normalisation between Matchete, RGBeta and numerical adapters; the fermion-first threshold plan; one-loop EFT Wilson mixing and matching; and independent numerical cross-checks. If any condition fails, the result must report the missing component rather than silently substitute restricted output.

## Numerical result interpretation

An optimizer fit to oscillation data is not a parameter-free prediction. Keep model inputs, fitted Yukawas, source hash/manifests, solver accuracy, accepted/failed cases and numerical output separate. Reports with a name suggesting completion or validation do not independently certify the general UV physics. Saved restricted-potential results may be retained in an **external historical archive** for reproducibility but should not be described as outputs of this general theory.

## References

- `T3_Honours_Thesis` (user-supplied November 2026 thesis draft): ordinary model catalogue, one-loop T3 analysis and B(-1)/C(-1) additional contributions. The thesis describes calculations and should not override present source limitations.
- Restrepo, Zapata and Yaguna, arXiv:1308.3655v2: radiative T3 model classifications (note their hypercharge convention).
- Ma, arXiv:hep-ph/0601225: scotogenic example.
- Merle and Platscher, arXiv:1507.06314v3: scotogenic running and matching issues.
- Aebischer, Bresciani and Selimović, arXiv:2512.16890v1: general EFT anomalous-dimension framework, requiring specialization and convention matching before use.
