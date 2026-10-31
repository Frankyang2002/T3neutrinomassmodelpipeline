# Workstream C v23 — final consolidated leading one-loop benchmark

## Purpose

One local-only batch for the archived T3 B(-1), C(-1), **two original Yukawas / no independently switched-on mixed scalar quadratic**, with three heavy fermion flavours. It consolidates v18, v19, v21, v22 and runs the repo's *existing* final SM+Weinberg RGE as an explicitly **nonphysical SM-only background proxy**.

This is intended to end further triangle topology/pole/fixed-order numerics batches for **rough estimates**. It does not establish a physical renormalized T3 Wilson coefficient. It never changes the UV Lagrangian, existing matching code, RGE modules, git repository or GitHub.

## Install at project root

From PowerShell:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive `
  -LiteralPath "$HOME\Downloads\T3_WorkstreamC_FinalConsolidated_v23.zip" `
  -DestinationPath . -Force
& .\scripts\RunT3BCFinalConsolidated.ps1
```

Full files delivered:

- `scripts/FinalizeT3BCLeadingBenchmark.experimental.py`
- `scripts/RunT3BCFinalConsolidated.ps1`
- `docs/WorkstreamC_FinalConsolidated_v23_README.md`

Local output (never included as an input to source code):

- `output/audits/t3_bc_consolidated_leading_reference_v23/summary.json`
- `output/audits/t3_bc_consolidated_leading_reference_v23/summary.md`
- `output/audits/t3_bc_consolidated_leading_reference_v23/masses.csv`

Optional source-only self-test and comparison without RG evolution:

```powershell
python .\scripts\FinalizeT3BCLeadingBenchmark.experimental.py --self-test
python .\scripts\FinalizeT3BCLeadingBenchmark.experimental.py --skip-running
```

`--require-physical` is deliberately rejected. The normal default computes *all* checks, including illustrative RGE, and writes a report. Self-test only does not need project evidence.

## Physics and numerical definitions

Audited route (v22): full-mass all-heavy Matchete ordered coefficient \(A_5\) = full-mass F/Y/X assembled ordered coefficient, tested in 12 configurations/model for each of three flavour structures. The project convention is

\[
 C_5=A_5+A_5^T,\qquad m_\nu=-\frac{v^2}{2}C_5,\qquad \hbar\equiv\frac{1}{16\pi^2}.
\]

v23 reads v22's complex symmetric 3x3 coefficient and v18's frozen-reference 3x3 coefficient, enforces equality and source SHA-256 provenance, and uses SVD/Takagi singular masses in eV. Its separate numerical loop cross-check evaluates for all 3 fermion masses and all scalar pairs XX/XY/YY:

\[
 \Phi(M_F^2,M_a^2,M_b^2)
 =-M_F^2\int_0^\infty \frac{t\,dt}
 {(t+M_F^2)(t+M_a^2)(t+M_b^2)}.
\]

This Euclidean *mass integral* is independent of the Matchete algebra, but **the overall minus sign follows the signed comparison convention introduced in v17**, not an independent Minkowski UV vertex/fermion-loop sign calculation. Do not claim v23 independently fixes the physical sign or the Majorana symmetry factor.

### Running used by this batch

The fitted reference contains \(y_1,y_2,M_F,\lambda_{T3}\) and SM inputs **defined at \(10^7\) GeV**. v18/v22 freeze the BSM values and evaluate the loop at the current Matchete \(M_X,M_Y\), while keeping the two additional direct quartics explicitly zero. v23 uses the existing `Numerical.running.SMWeinbergEvolution` with the same full flavour matrices to:

1. Run the SM-only parameters (with `K=0`) from \(10^7\) GeV to \(M_X\).
2. Insert the exact full-mass v22 source-normalized \(C_5\) at \(\mu=M_X\).
3. Run SM+Weinberg (including matrix \(K\)) from \(M_X\) to **100 GeV**.
4. Produce low-scale illustrative singular masses and squared differences.

**Crucial limitation:** Step 1 is *not* the true T3 UV/EFT1/EFT2 parameter trajectory, and the Wilson matrix insertion in step 2 is not a certified physical matching boundary. This is a transparent reference to estimate the size of SM+Weinberg evolution only. Running one-loop C5 at one loop is beyond strict \(O(\hbar)\) order, and BSM corrections can be important. Do not use the running output as a precision prediction.

## What is and is not verified

| Statement | Result in v23 |
|---|---|
| Source files unchanged since v18/19/21/22 | Fail-closed SHA-256 checks |
| Both available leading \(y_i y_j\lambda\) triangle topologies complete for archived source | v21 evidence required |
| Restricted triangle hard/soft \(\overline{\mathrm{MS}}\) cancellation | v19 evidence required |
| Direct one-scale all-heavy vs F/Y/X full-mass Matchete \(C_5\) agreement | v22 12 cases/model plus v18 3x3 cross-check |
| Independent convergent Euclidean triangle *mass integral* | 3 fermions × 3 scalar pairs × 2 models |
| Approximate SM-only running at/under \(M_X\) | Yes, project-native solver |
| Absolute UV Minkowski Feynman-rule sign and vertex/combinatorial normalization independent of Matchete | **NOT certified** |
| UV→EFT1→EFT2 fully matched BSM running and finite field/parameter/EOM matching | **NOT certified** |
| Physical \(C_5\) / measured neutrino masses | **NOT certified** |

The report explicitly records the outstanding blockers. Achieving strict physical certification would require a new independent UV amplitude derivation and full parameter/Wilson threshold transport, which are not available from these archived outputs alone. It cannot be honestly guaranteed within one batch.
