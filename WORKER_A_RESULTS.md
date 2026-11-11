# Worker A — source-backed direct-channel audit

## Verified from the uploaded local source archive

- A(-4): 135 raw one-loop entries, including 6 direct LLHH terms (3 Laurent terms in each chirality), 48 LLSS, and 81 other operators.
- D(-2): 1200 raw one-loop entries, including 12 direct LLHH terms (6 baseline Y1/Y2, 6 extra Y1/y2Conj), 352 LLSS, and 836 other operators.
- Both Laurent v2 audits explicitly set `uv_hard_minus_eft1_soft_subtraction_performed=false`, `local_MSbar_counterterms_added=false`, and `renormalised_MSbar_Wilson_boundary_computed=false`.
- All three heavy fermions are integrated together at 1e5 GeV in the source setup, while S1 and S2 are retained until 1e3 GeV.
- D(-2) y2Conj is numerically zero in the supplied benchmark at muF. Its raw extra channel is not evidence for an independent coupling under the restricted Y1/Y2 theory.
- Relative off-diagonal heavy mass Frobenius norm: A(-4) 5.13e-10, D(-2) 5.23e-10. The audit therefore uses the existing documented diagonal-MF numerical approximation, not a certified generic chiral rotation.

## Computed numbers (GeV^-1)

| Raw/formal diagnostic | A(-4) | D(-2) |
|---|---:|---:|
| Norm of raw direct hard 1/epsilon residue | 2.44523139e-16 | 2.55026756e-16 |
| Norm of raw direct hard regular term at muF | 1.44552295e-16 | 1.92877387e-16 |
| Norm of A4 fixed-coupling formal hard/soft remainder using scalar masses at muS | 2.01015905e-15 | not evaluated |

**None of these is a physical one-loop EFT1 boundary.** The formal A4 remainder takes a formula from the existing `a4_direct_hard_soft.py`; it is a diagnostic that mixes Yukawa couplings at F with candidate scalar masses at S. It must not be added to the tree-level EFT1 Wilson tensor.

## Outstanding physics

The uploaded archive contains first-threshold UV Matchete expressions, formal Laurent files, and benchmark numerical inputs, but no independently certified scalar-soft final EFT1 matching result in a common GammaCC-to-Weyl basis. Full-basis operator projection, EFT1-soft/UV counterterms, sign/CG normalization, and (for D(-2)) the scalar self-bilinear mass basis are still unresolved. No production code was modified.

## Installation

Put `scripts/WorkerAFiniteChannelAudit.py` into the project path of the same name. Run in activated venv from project root:

```powershell
$runId = Get-Date -Format 'yyyyMMdd_HHmmss'
python scripts\WorkerAFiniteChannelAudit.py --repo . --output "output\audits\worker_a_direct_channels_$runId.json"
```

The output path is fresh and the script uses exclusive creation to avoid overwriting prior results.
