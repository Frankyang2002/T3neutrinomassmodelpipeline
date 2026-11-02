# T3 fermion-threshold spinor worker — fail-closed implementation

Repository path: `Numerical/lo_certification/spinor/`.

## Fixed EFT sequence

`UV -> integrate ALL heavy BSM fermions together at mu_F -> EFT1 with BOTH BSM scalar multiplets and SM fermions still active -> integrate BOTH scalars together at mu_S -> SMEFT`.

This workstream only handles the **first** threshold. Never split fermion thresholds, remove scalars at `mu_F`, or treat a PR Hermitian conjugate as another PL Wilson coefficient. A scalar complex multiplet may be represented as two charged `Plus`/`Minus` real-coordinate blocks; these are not additional scalar multiplets or thresholds.

## Inputs and provenance

The numerical first-threshold source data are not regenerated; they are the immutable **provisional** reference at `output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/`. The runner requires the 16 `boundary.json` files, `summary.json`, exact source Matchete SHA256, and the original `ordered_wilson_tensor_v4`, `general_wilson_polynomial_batch_v2`, and `general_wilson_normalization_audit_v1` files with matching SHA256. It refuses to synthesize any missing input or to overwrite an existing output with changed content.

All 16 model labels are checked. The full 3-generation complex components are carried source-by-source, including zero benchmark sources; chirality is inspected in the actual native InputForm operator, and each source's CG identity is retained as supplied by the previously computed Matchete tensor. Duplicate sources, mismatched source inventory, an incomplete scalar basis, wrong orbital multiplicities, wrong scalar charge-block spans, altered provenance, non-finite coefficients, and mismatched sums all fail.

Each independent six-coordinate row has indices `[lepton SU2 i, lepton flavour p, lepton SU2 j, lepton flavour q, real scalar a, real scalar b]` (1-based). The independent entry is the **coefficient per member of its distinct permutation orbit**, not a coefficient to be multiplied again by 1/2/4. Source-resolved and total tensors are **formal commuting-placeholder Wilson data**, *not* a calibrated physical Weyl Wilson tensor. The code preserves the original native spinor operator text in the new interface.

## Conditional signed spinor algebra

In an explicitly declared Weyl basis with `epsilon=i sigma2`, `gamma5=diag(-I,+I)`, `PL=diag(I,0)` and `C=i gamma2 gamma0`, matrix multiplication proves

`(chi,0)^T C PL (eta,0) = chi^T epsilon eta`

and the same expression with `C^-1` has a **minus** sign. An exact rational exterior algebra verifies from the **declared** UV actions:

- Majorana: `-M/2 h^T eps h + sum_i y_i h^T eps l_i` gives `+1/(2M) sum_ij y_i y_j l_i^T eps l_j`, with cross terms `+1/M`.
- Dirac: `-M u^T eps v + u^T eps A + v^T eps B` gives `+1/M A^T eps B`.

These exact identities validate Grassmann ordering, crossings in this declared action, and symmetry factors. **They do not establish** that Matchete's symbolic `GammaCC` equals this chosen C rather than C inverse, nor that the native external-leg ordering and the master Wilson normalization have a matching phase.

The saved installed Matchete source evidence provides `CConj`, `GammaCC` and projector *symbolic* rules and shows that `Proj[-1]` denotes left chirality while `Proj[1]` denotes right chirality. It does not fix the numerical charge-conjugation matrix. The A_-4 native zero-momentum and B/C native CG probes in the historical directory also explicitly leave independent full-amplitude / master-phase checks open. Nothing in this package may turn this open factor into `+1`.

## Handoff for EFT1 and scalar matching

Each output `source_<summary SHA prefix>/<model_tag>/boundary.json` contains a fully indexed source-resolved complex numerical LLSS tensor, its sum, all source InputForm text, scalar-basis metadata, and strict provenance. The JSON schema is `Numerical/lo_certification/spinor/handoff.schema.json`. The stable `normalization.matchete_GammaCC_to_master_Weyl_complex_factor` and `physical_master_weyl_boundary` fields are `null` by design. The downstream `safe_for_physical_EFT1_or_scalar_matching` flag is `false`; only **structural debugging or conditional algebra** may use the data. A physical worker must fail closed until independently calibrated.

The source data are full 3-family complex **benchmark tensors**; they are not an all-parameter symbolic re-evaluation for arbitrary new Yukawa matrices. Any change of Yukawa input requires a fresh native matching run and refreshed checked source tensors. **No hard one-loop F matching** is present, so these data must not be advertised as a complete LO one-loop neutrino-mass result.

## Commands

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$env:PYTHONPATH = (Get-Location).Path
python -X utf8 -m unittest discover -s tests -p 'test_lo_spinor_*.py' -v
python -X utf8 scripts\RunLOSpinorCertification.py --check-only
python -X utf8 scripts\RunLOSpinorCertification.py
python -X utf8 scripts\RunLOSpinorCertification.py --require-physical
```

The first two runner commands require the local original source release and reference output directories. `--require-physical` is **expected to exit 2** until all three blocked physical verifications are completed; it will not deliver plausible-looking but uncertified coefficients. `--check-only` avoids writes. The normal runner writes exclusively within `output/general_observable/lo_certification_spinor/`. No RGE, mass-matching or SMEFT file is modified.

### Exact remaining blockers

1. Installed Matchete `GammaCC`, `CConj`, native `PL/PR`, and fermion-flow rules must be combined with the numerical definition of `C`/`C^-1` and matched against the master's **signed** Weyl bilinear (including operator prefactor).
2. A **connected** four-point UV-to-EFT1 native amplitude must be computed for all channels, with Dirac and Majorana exchange and both crossings, native CG contractions (including C_-1 bridge tensor), and external Grassmann ordering; saved diagnostics are not sufficient.
3. The full first-threshold **hard one-loop fermion** contribution has not been included in the saved source-resolved F-tree release. It must be computed before declaring the complete leading-one-loop mass certified (or rigorously shown unnecessary at the selected order).

The release created here is therefore a usable **provisional** machine-readable first-threshold boundary, not a physically normalized result.
