# T3 scalar relevance — Batch 2: leading one-loop Wick selection

## Question

Which of the symmetry-allowed scalar quartics can appear at **one scalar-quartic vertex** in the **leading unbroken-phase one-loop Weinberg-sector** amplitude with two original T3 Yukawa insertions?

This is more restrictive than counting the number of Higgs legs, but **not** a calculation of the full general scalar-potential neutrino mass. The result is conditional on diagonal complex scalar propagators and the selected chirality/field basis.

## Sources

- Original Yukawas and the topology vertex: `Lagrangian/interactions/T3Topology.wl`, inspected at `a02e04c6def5cf8522b2f612ccd6bff7d9963b91`.
- SU(2), hypercharge and Hermitian-conjugation quartic basis: `Lagrangian/interactions/GeneralScalarOperators.wl` and `GeneralScalarVerifiedQuartics.wl`.
- The input `output/scalar_relevance/operator_classification.json` is the strictly complete Batch 1 operator classification, tied to original repository artifacts by five SHA256 source hashes.
- Saved independent symbolic C5 observations come from `output/general_matching_symbolic_report.json` and `output/general_t3_physics_summary.json`, **indirectly** via Batch 1. This batch does **not** invoke Matchete.

## Derivation: necessary propagator and external-field matching

The source file defines the **lepton-side** vertices schematically by

\[
 y_1\,\bar L F^c S_1, \qquad y_2\,\bar L F S_2^\dagger,
\]

plus their Hermitian conjugates. Thus, the scalar ends attached to the two chosen lepton-side Yukawa vertices can be `S1` and/or `S2dag`. In the **complex, diagonal, gauge-unbroken scalar basis** the only scalar propagators are

\[
\langle S_1 S_1^\dagger\rangle, \qquad
\langle S_2 S_2^\dagger\rangle.
\]

A scalar quartic that closes those two lines must therefore contain one of three scalar-end multisets:

| Scalar ends of quartic | Source Yukawa pair |
|---|---|
| `(S1dag)^2` | `y1 y1` |
| `S1dag S2` | `y1 y2` |
| `S2^2` | `y2 y2` |

For the holomorphic `LLHH` convention used in the existing Matchete matching, the corresponding quartic orientation must contain `Hdag^2`. This orientation is **calibrated by the Hermitian conjugate of the original** `HH S1 S2dag` topology vertex. The script tests both the named quartic and its Hermitian conjugate.

Consequently, an additional quartic with `H Hdag` cannot provide the required Higgs pair at that one vertex, and a quartic with `Hdag^2` but the wrong scalar conjugation cannot close the two scalar lines with the assumed propagators. Matching these external/propagator legs is a **necessary**, not sufficient, criterion: fermion chirality, charge conjugation and SU(2) CG factors still matter. The saved matching result independently supports that the four Wick-eligible **additional** quartic slots enter the extracted symbolic C5.

### Verified field-flow inventory from the uploaded Batch 1 output

| Additional quartic category | Across 16 models | Single insertion at fixed one-loop C5 |
|---|---:|---|
| `Hdag^2` with compatible scalar pair | 4 | Wick-eligible and observed in saved symbolic C5 |
| `H Hdag` with new scalar pair | 10 | Excluded by the specified Higgs orientation |
| `Hdag^2` with incompatible scalar pair | 10 | Excluded by diagonal Wick contractions |
| Four new BSM scalar legs, no H | 26 | No direct single-vertex contribution |

The four eligible vertices are, separately in `B_-1` and `C_-1`:

- `lambdaT3ExtraQ020020I1`: `Hdag^2 S2^2` (`y2 y2` pairing);
- `lambdaT3ExtraQ020200I1`: `Hdag^2 S1dag^2` (`y1 y1` pairing).

`lambdaT3` is also eligible in all 16 models through its Hermitian conjugate (`y1 y2`).

The ordinary `Hdag H S_i^dag S_i` portal vertices (two baseline portal families per model) similarly fail the **direct** holomorphic one-vertex test. This does not show their physical couplings are negligible: portal-induced scalar masses after EWSB, vacuum properties, wavefunction effects and scalar RGEs are separate questions.

## Hard limitations — do not overinterpret the 20 exclusions

1. **Quadratic mixing:** Eight allowed additional quadratic slots occur across the 16 configurations. If a holomorphic/mixed mass coefficient is nonzero, propagators such as `\langle S_i S_i\rangle` or off-diagonal contractions can exist in an appropriate representation/mass basis. The diagonal field-flow test then **does not apply**. These are unsuppressed tree-level effects; their zero values are boundary/model assumptions, not loop-order conclusions.
2. **Two-threshold calculations:** The graph test does not calculate or renormalize a full fermion-first EFT, scalar-threshold matching, or the running/mixing of the `LLSS` Wilson operators.
3. **Loop order:** Additional four-BSM-scalar quartics cannot be the *only* quartic providing the two external Higgs fields. Extra insertions can matter at higher order; running and scalar self-energies are not erased.
4. **Higgs VEV effects:** In broken phase, portal interactions can change scalar mass matrices at order `v^2`. An expansion of dimension-five C5 in the unbroken gauge phase does not automatically include every higher-order `v^2/M^2` correction.
5. **RG generation:** Zero boundary conditions for omitted allowed couplings are generally not invariant under running. A coupling generated at one-loop inserted in a one-loop C5 is formally higher fixed order; for long scale hierarchies, logs require a numerical check.
6. **Evidence:** The four *observed* quartic contributions are grounded in saved symbolic matching reports, not a newly performed Matchete calculation or a certified Matchete–RGBeta physical coupling map.

## Execution (PowerShell, repository root)

Install the Batch 2 ZIP by extracting into the project root (no existing source needs to be overwritten):

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Scalar_Relevance_Batch2.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_scalar_vertex_batch2.py
python -m scripts.AuditT3OneLoopScalarVertices
```

A normal run checks Batch 1's recorded SHA256 hashes against all five existing local sources and validates the `T3Topology.wl` Yukawa convention before producing new outputs. If the sources changed, rerun Batch 1 first. `--skip-source-hash-check` exists **only** for reproducing this uploaded-report analysis in an isolated offline environment; it marks the report as not locally source-verified.

Outputs:

- `output/scalar_relevance/leading_c5_vertex_audit.json`
- `output/scalar_relevance/leading_c5_vertex_audit.csv`
- `output/scalar_relevance/leading_c5_vertex_summary.md`

The ZIP includes a **preview copy** of these outputs produced from the user's uploaded Batch 1 report with `source_provenance_verified_on_local_checkout=false`. Running the script in the actual repository regenerates the same classification and asserts the local source hashes.

## Next batch

Batch 3: quantify the eight allowed quadratic terms and their influence on the one-loop mass. Test zero-versus-nonzero quadratic mixing/splitting, identify when the diagonal-propagator selection rule fails, and record the benchmark restriction explicitly. A full general `C5` mass-eigenstate matching calculation is a separate, larger goal.
