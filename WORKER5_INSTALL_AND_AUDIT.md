# Worker 5 — joint scalar threshold completeness (audit only)

This ZIP contains additive files only. No `core.py`, `worker_c_joint.py`, solver, matched coefficient, or Matchete source is replaced.

## Physics statement

The verified existing kernel computes

\[
C_5^{\rm SMEFT}(\mu_S)=C_5^{\rm EFT1}(\mu_S)+\Delta C_5^{\rm LLSS\times HHSS}_{\rm bubble}(\mu_S), \qquad \mu_S=10^3\,{\rm GeV}.
\]

The bubble is the full one-insertion LLSS–HHSS zero-momentum finite MSbar bubble, evaluated in the joint *real mass eigenbasis* of both scalar multiplets. The `core.py` implementation uses

\[
\Delta {\cal A}_{fgijAB}=-\frac{1}{2(16\pi^2)}\sum_{ab}C_{fgijab}\lambda_{ABab}\overline B_0(0;m_a^2,m_b^2;\mu_S),
\]

with the off-diagonal eigenpair multiplicity 2. Projection onto its Weinberg tensor yields the matrix that is added exactly once to the inherited coefficient. Scalar species S1 and S2 are not decoupled sequentially. Quartic channels, when supplied, sum to the same bubble.

This is *not* certified as the complete finite one-loop matching of the general EFT1 operator inventory. Additional dimension-five insertions, counterterm/wavefunction effects, or scalar topologies must be treated as **potential** sources until the renormalized, scheme-matched, bubble-subtracted finite coefficient is established. A different Matchete expression is not by itself a new physical diagram.

## Verified saved D(-2) evidence in public repository (as inspected 2026-11-12)

`output/audits/worker_c_restricted_D_m2_20261112_052721_aaa5ec32/result.json` reports exact total pole cancellation for NF=1 after setting the separately introduced `y2Conj` source channel to zero; it does not certify finite renormalization, physical mass relations, operator normalization, or a nonbubble remainder.

`output/audits/worker_c_numeric_overlap_D_m2_20261112_060707_96a2df24/result.json` reports provisional **hard finite-part** norm `9.359924597942799e-16 GeV^-1` and **total finite-part** norm `4.52810332477528e-18 GeV^-1`. The hard residue has maximum `4.021724619722921e-15 GeV^-1`; total residue is about `4.733165431326071e-30 GeV^-1`. For these two saved matrices the numerical Frobenius alignment is `1.0`, the total/hard norm ratio is `0.004837756199201119`, and the least-squares complex-line residual is `1.8553800501533808e-12`. This tests source-level flavour proportionality only. The total explicitly includes inherited and running terms. The hard/total finite-part difference is **not** a missing physical threshold coefficient.

Crucially that saved benchmark uses **3000 GeV**, fermion 4000 GeV, S1 3200 GeV, S2 2400 GeV, while the production joint threshold is **1000 GeV**. The archive reports `NotAttemptedNoCanonicalTensorInput` for bubble overlap. Thus its numbers must never be compared as a same-scale correction to production C5.

A source `y2Conj` term can be an artifact of a generalized intermediate source. Restricting it to zero is not proof that every other term has been consistently restricted, nor proof that it is an independent coupling in the actual Y1/Y2-only T3 theory.

## Flavour implications

The bubble inherits flavour from the symmetric EFT1 LLSS Wilson tensor and does not create arbitrary new flavour tensors. Other source-level contractions might have different flavour structures, but no independently normalized and subtracted remainder has been established. This audit computes normalized complex Frobenius overlaps only for available *diagnostic* matrices, and labels the comparison explicitly as **hard vs total**, never **new correction vs bubble**. It cannot determine an additional shift of neutrino masses or angles.

For a genuinely certified matrix \(\delta C_5\), one would compare \(\|\delta C_5\|_F/\|C_5^{\rm current}\|_F\) and the complex Frobenius overlap, then propagate both complex-symmetric matrices through the same SMEFT evolution and Takagi diagonalization. Until that exists, even the sign of a missing neutrino mass shift is not fixed. The conservative conditional norm bound at the matching scale is \(\|\delta m_\nu\|_2\leq v^2\|\delta C_5\|_2/2\), for the convention \(m_\nu=-v^2 C_5/2\); verify actual project normalization first. No numerical bound is asserted.

## Installation (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$HOME\Downloads\T3_Worker_5_Scalar_Matching.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
```

## Test and run (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
python -m unittest discover -s tests -p 'test_worker5_scalar_matching.py' -v
python -m Numerical.diagnostics.worker5_scalar_matching --root .
Get-ChildItem output\audits\worker5_scalar_matching -Recurse -File
```

The script uses only existing Worker C gate and endpoint code, performs zero Matchete evaluations, and creates timestamped non-overwriting audit files under `output/audits/worker5_scalar_matching/`. It will fail closed on invalid endpoint metadata. A missing local saved profile remains missing rather than fabricated. No production coefficient is updated.

## Complete installed paths

* `Numerical/diagnostics/worker5_scalar_matching.py` — exhaustive read-only inventory, hashes, existing Worker C source gates and canonical endpoint checks; diagnostic finite-part flavour alignments
* `tests/test_worker5_scalar_matching.py` — linear algebra, scale-log and inherited-plus-bubble regression checks
* `WORKER5_INSTALL_AND_AUDIT.md` — physics audit, source status, limitations, installation and execution

Existing Worker C files are deliberately reused without replacement. GitHub was read-only; local repository was not accessible in this execution environment, so actual Windows checkout tests must be performed by the commands above. Archive facts were checked against public GitHub `main`, not against any uncommitted local changes.
