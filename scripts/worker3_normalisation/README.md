# Worker 3 — Matchete-to-Weyl normalisation, read-only audit

All paths below are relative to `C:\thesistemp\T3neutrinomassmodelpipeline`. No production source is replaced. The `audit.py` script is independent; it additionally checks local production-source markers, and the unit suite calls the existing `project_c5` and `_weinberg_template` functions when available.

## Execute

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
python -m unittest discover -s scripts/worker3_normalisation -p 'test_*.py' -v
python scripts/worker3_normalisation/audit.py --root . --output output/audits/worker3_normalisation
```

Optional **explicitly supplied** reference matrices (not an implicit data search):

```powershell
python scripts/worker3_normalisation/audit.py --root . --output output/audits/worker3_normalisation --fixtures 'path\to\your_worker3_trajectories.json'
```

Input schema:

```json
{"format":"worker3-c5-trajectories-v1","cases":[{"name":"fixture_A","C5_start":{"real":[[1,0,0],[0,2,0],[0,0,3]],"imag":[[0,0,0],[0,0,0],[0,0,0]]},"channel_A":{"real":[[1,0,0],[0,0,0],[0,0,0]],"imag":[[0,0,0],[0,0,0],[0,0,0]]},"channel_B":{"real":[[0,0,0],[0,2,0],[0,0,3]],"imag":[[0,0,0],[0,0,0],[0,0,0]]}}]}
```

Each matrix is a complex, symmetric 3×3 coefficient **in the same charged-lepton mass basis**. `C5_start` is the initial matrix; `channel_A`, `channel_B` are independent contribution matrices for sensitivity testing. The default synthetic input is 32 deterministic nondegenerate complex-symmetric matrices, and is **not** a physical benchmark replay. Stage norm ratios are tested using synthetic stage matrices, not full ODE solutions. This distinction is preserved in generated JSON.

## Convention table and physical caveats

`output/audits/worker3_normalisation/convention_table.csv` contains the machine-readable source reference ledger. Source references point to exact file paths, symbols, and GitHub blobs in `source_manifest.json`; the audit does not claim that the installed Matchete package was executed. The source file `Numerical/lo_certification/physical_normalisation/B1_MATCHETE_WEYL_DERIVATION.md` already records a restricted chosen-spinor-basis identity and its remaining native-package blocker.

| Item | Derivation/status |
|---|---|
| Real scalar components | `H0=(hR+i hI)/sqrt(2)` gives a factor 1/2 for two H0 fields. |
| Symmetric Weinberg definition | `(1/2)C5_ij(L_i.H)(L_j.H)+h.c.` gives `C_master[nu_i,nu_j,hR,hR]=C5_ij/4`, assuming the **explicit unsuppressed master polynomial**. |
| Master Wilson symmetry | `C_ijab=C_jiab=C_ijba` is enforced by the boundary parser; never multiply by both orbit sizes after summing the full tensor. |
| Scalar quartic | `lambda_abcd phi_a phi_b phi_c phi_d / 4!`; source portal CG and component factors require separate verification. |
| Tree LLSS | Fermions integrated out, both heavy scalar multiplets retained; not a direct Weinberg coefficient. Separate `XX`, `XY`, `YY` and conjugate sectors before projection. |
| Direct LLHH | Distinct source class and loop order; its normalisation cannot be inherited from LLSS. |
| Matchete `GammaCC`, `Proj[-1]` | Source PL chirality identified in existing audit; numerical native four-component C and phase not established by isolated algebra. |
| Matchete `Proj[1]` | Conjugate barred-lepton sector; must not be summed as independent holomorphic PL. |
| Dirac/Majorana | Chosen-basis B1 derivation supports canonical mass `-M chi chi/2+h.c.`; no extra Majorana factor. Native source basis may still differ. |
| SU(2), triplet CG | Changing epsilon sign, doublet basis, or triplet `G=U.T@U` requires transforming **every vertex and mass** consistently. |
| Universal source map | **Uncertified**. Channel/model/flavour dependence cannot be excluded by a successful pipeline fit. |

## Exact scaling theorem and its limits

For a real Wilson normalisation `r != 0` with fixed renormalisable couplings, an RGE/matching system **homogeneous and real-linear** in the full vector of Wilson coefficients obeys `F(rC)=rF(C)` even when terms in `C*` occur. By uniqueness of the IVP, if every boundary, matching injection, and independent source is also multiplied by `r`, then `C(t;rC0)=r C(t;C0)`. Real-linear thresholds, sequential stages and the linear Higgs projection retain this property. The theorem is conditional: it fails for additive unscaled sources, Wilson feedback into running renormalisable couplings, quadratic Wilson terms, or inconsistent rescaling of some channels.

For `m_nu=-(v^2/2)C5` and fixed v, this implies

- `m_i -> |r| m_i`, `Delta m^2_ij -> r^2 Delta m^2_ij`;
- `m_i/m_j`, nonzero mass-squared ratios, ordering, `|U_PMNS|`, and three mixing angles stay unchanged;
- a stage ratio `m_i(mu_b)/m_i(mu_a)` is unchanged if the **same nonzero r** multiplies both stage matrices (likewise homogeneous matrix norms); the difference of masses and absolute stages are not protected;
- eigenvector conventions/angles at exact mass degeneracy and ratios with vanishing denominators are undefined.

For a complex constant `z`, **matrix-level** `m -> z m` retains `|U|`, angles, ratios, and gives masses `|z|m`; Majorana phases may change by a common convention-dependent shift. But along an RGE containing `B(C*)`,

`beta(z C) - z beta(C) = (conj(z) - z) B(C*)`.

Thus general complex rescaling **does not commute with evolution**, unless `z` is real, the conjugate term vanishes on the trajectory, or other special cancellations apply. For the current `GeneralEFT1C5Mixing.py`, `conjugate_coefficient_yukawa` is explicitly present. The regression tests include nonzero complex obstruction and a zero-conjugate-term control.

For independent source channels `C5=sum_a r_a M_a`, only a common `r_a=r` guarantees exact protection. Unequal factors generally rotate Takagi eigenvectors and change mass ratios, mixing angles, and stage-running ratios; if all channel matrices are proportional, or a special accidental alignment holds, shape may remain invariant. A fitted normalisation does not prove universal `r`.

## Certification missing — concrete closure experiment

Obtain from the **same installed Matchete version and identical fixed SU(2) basis**: (1) exact `PlusHc` source operator `InputForm` for each LLSS (`XX,XY,YY`, PL and PR) and direct LLHH channel; (2) the native `GammaCC`, `Proj[-1]`, `Proj[1]` charge-conjugation and two-component epsilon definitions, including their signs; (3) a one-component amputated amplitude/functional derivative computed directly from the Matchete operator and independently from `C_ijab psi_i psi_j phi_a phi_b` with canonical external spinors; (4) full CG contraction including doublet epsilon, triplet metric, rotations and scalar real-field `1/sqrt(2)` per field; (5) matching of all identical-particle permutations and the symmetric flavour map, with dimensional factors and h.c. handled once; (6) repeat for a complex off-diagonal flavour component in each model class and each operator channel, including two non-collinear contributions; (7) directly compare computed amplitudes. Only after these pass for every class can a universal conversion be claimed. No factor is invented by this audit.

**Scope**: This suite is intentionally a convention regression suite with synthetic or user-exported matrices, not a substitute for running the local sixteen-model coupled-IVP pipeline. It does not supply certified Wilson matching inputs or physical neutrino predictions.

## Synthetic coupled-IVP stage regression

```powershell
python scripts/worker3_normalisation/ivp_regression.py --output output/audits/worker3_normalisation
```

Uses SciPy `solve_ivp` on an independent two-channel, two-stage **real-linear/conjugate-linear** test system with cross-channel coupling, a threshold injection, and separate stage running. Checks 16 initial conditions against common real scaling (`-2.4`), PMNS moduli-derived mixing angles, mass ratios, three stage norm ratios, complex-scaling obstruction and unequal channel factor sensitivity. This **does not** execute or validate the physical coupled T3 IVP: numerical source and an installed physical conversion factor are not available in this execution environment. All such output fields are explicitly labelled synthetic.
