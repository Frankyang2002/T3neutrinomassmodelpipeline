# T3 B(-1)/C(-1): native Y/X finite-term audit v1

## Scope

This is an **exact textual partition of the native local Matchete operator
expressions** already saved from the independent-canonical Y-only matching
and fresh-kernel X-only matching. It does not rerun Matchete, change any input
matching source, run RGE, or substitute benchmark parameters.

The starting amplitudes have the form

\[
\Delta\mathcal L^{\mathrm{raw}}_{\mathrm{LL}\phi\phi}
= \frac{\mathcal P}{\epsilon}
+ \mathcal L\,\log\frac{\bar\mu^2}{M_{X,Y}^2}
+ \mathcal K,
\]

where the output records the **unsimplified Matchete InputForm** for:

- `pole_residue`: the expression multiplying the formal `1/ϵ`;
- `pole_free_unrenormalized`: logarithmic plus scale-independent terms;
- `at_mu_equal_heavy_mass_unrenormalized`: the scale-independent terms alone, after formally setting the **specific Matchete** `Log[μbar2/MSX^2]` or `Log[μbar2/MSY^2]` to zero.

The splitter preserves the operators' flavour, SU(2) Clebsch-Gordan,
chiral, Yukawa, mixing-spurion, hbar and normalization structures.
It separates additive terms at zero bracket depth and **only** removes a
single multiplicative `ϵ*` factor from a known denominator when computing
pole residues. It aborts for unknown poles/log arguments, mixed logs and
poles in one monomial, incorrect field inventory, or inconsistent source
reports. No algebraic cancellations or SU(2)/flavour projections are
claimed. Counts are not counts of independent invariant operators.

## Reference results, verified with archived Matchete input files

| Model | Threshold and operator | 1/ϵ terms | log terms | constant terms |
|---|---|---:|---:|---:|
| B(-1) | Y: LLHH | 50 | 50 | 32 |
| C(-1) | Y: LLHH | 50 | 50 | 32 |
| B(-1) | Y: LLXX | 240 | 240 | 186 |
| C(-1) | Y: LLXX | 198 | 198 | 186 |
| Both | Y: LLHX | 0 | 0 | 0 |
| Both | X: LLHH | 18 | 18 | 0 |

The X LLHH raw expression has no explicit scale-independent additive terms,
so the **formal** X LLHH `at_mu_equal_heavy_mass_unrenormalized` text is `0`.
This is **not proof of zero physical one-loop X matching**: the X-only
Matchete run starts from tree-level Y EFT2, and complete EFT subtraction,
field renormalization, counterterms, and the inherited Y contribution are
not yet combined. The Y constant terms are not individually guaranteed
nonzero after tensor simplification and numerical flavour contraction.

### Inputs consumed locally

- `output/general_t3_mass_eigen_matchete/eft1_scalar_Y_only_canonical_probe/B_m1_final_y_only_canonical_match_raw_llhh.inputform.txt` (likewise `C_m1`, and `raw_llhx`/`raw_llxx`).
- `output/general_t3_mass_eigen_matchete/eft1_scalar_Y_only_canonical_probe/B_m1_final_y_only_canonical_match.json` (likewise C).
- `output/general_t3_mass_eigen_matchete/eft1_scalar_xy_fresh_kernel_v3/B_m1_final_x_from_saved_y_match_x_followup_x_raw_llhh.inputform.txt` (likewise C).
- `output/general_t3_mass_eigen_matchete/eft1_scalar_xy_fresh_kernel_v3/B_m1_final_x_from_saved_y_match_x_followup.json` and matching `.log` (likewise C).

The original X v3 report may remain stuck at `XFollowupOneLoopMatchingInProgress`
due to the **post-export `AssociateTo` reporting bug**. The runner does not
promote that status; it independently demands the positive source/handoff
flags, two completed Matchete pipelines, the post-export checkpoint, and the
36-term LLHH file. For stronger provenance, first run
`scripts.RecoverT3BCXRawMatching --all` from the previous v4 package.

### Install and run (Windows PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_finite_threshold_term_audit_v1.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p "test_t3_bc_matchete_finite_terms_v1.py" -v
python -u -X utf8 -m scripts.ExtractT3BCFiniteScalarTerms --all
```

Outputs are saved to
`output/audits/t3_bc_scalar_finite_term_extraction_v1/`.
Each model subfolder contains 4 sector audit JSON files and 12 exact
InputForm expressions (`pole_residue`, `pole_free_unrenormalized`,
`at_mu_equal_heavy_mass_unrenormalized` for each Y/LLHH, Y/LLHX,
Y/LLXX, X/LLHH); top-level `summary.json` and `summary.csv` record
counts and source SHA256 digests.

The included `docs/T3_BC_FINITE_THRESHOLD_REFERENCE_V1.json` records a
reference audit computed **from the files provided in the B/C v3+v4
archives**, independent of the user's next local execution. If a hash differs
locally, compare the specific sources before using any extracted expression.

## Outstanding physical work — do not skip

1. **Renormalization/subtraction:** establish the correct MS-bar pole
   counterterms and EFT loop subtraction in precisely the Matchete scheme;
   do not equate deleting pole monomials with physical matching.
2. **Operator projection:** project the resulting finite tensor-valued
   LLHH and LLXX structures into the numerical `(3,3,8,8,6)` Wilson basis,
   with all SU(2), Majorana, flavour and hbar conventions.
3. **Sequential inheritance:** propagate the full Y-threshold correction,
   not only its tree-level EOMEFT, into the X threshold.
4. **Running parameter map:** substitute actual couplings/masses at the
   model-dependent muY and muX, not the preliminary Matchete parameter JSON.
5. **Equivalence proof and scale checks:** resolve original F-first EOMEFT
   versus independently reconstructed EFT1, and verify matching-scale
   dependence against the EFT1/EFT2 RGE and known limits.

`physical_sequential_matching_certified=false` is set deliberately.
Do not plug these formal expressions directly into certified neutrino masses.
