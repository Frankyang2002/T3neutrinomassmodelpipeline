# Workstream C v21 — B/C(-1) renormalizable one-loop LLHH topology closure

**Status:** Conditional topological closure for the **archived, two-Yukawa UV source**, not full renormalized physical Weinberg matching. All operations are local, read-only with respect to existing project sources, and produce one experimental audit report.

## Executables (repo-relative)

- `scripts/CloseT3BCOneLoopTopologies.experimental.py` — full replacement independent v21 source/proof gate.
- `scripts/RunT3BCTopologyClosure.ps1` — self-test plus local B/C run.
- `docs/WorkstreamC_OneLoopTopologyClosure_v21_README.md` — this documentation.

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_OneLoopTopologyClosure_v21.zip" -DestinationPath . -Force
& .\scripts\RunT3BCTopologyClosure.ps1
```

Or:

```powershell
python scripts\CloseT3BCOneLoopTopologies.experimental.py --self-test
python scripts\CloseT3BCOneLoopTopologies.experimental.py --output output\audits\t3_bc_topology_closure_v21\summary.json
```

Output: `output/audits/t3_bc_topology_closure_v21/summary.json`.

## One-loop graph budget

For connected one-particle-irreducible four-point `LLHH`, with four external legs, one loop and renormalizable vertices, `sum_v (valence(v)-2)=4`. Two original T3 Yukawa vertices each contribute 1; the remaining budget is exactly 2. Apart from arbitrary two-point mass insertions, the only partitions are **one quartic vertex** or **two extra cubic vertices**.

The archived B/C UV source contains two lepton-side T3 Yukawas, `y1 Bar[l] CConj[F] S1` and `y2 Bar[l] F Bar[S2]` (plus Hermitian conjugates). For `alpha=-1`, the relevant new scalar species and original Yukawa scalar ends are SU(2) doublets with `2Y=-1`; the Majorana F is either an SU(2) singlet (B) or triplet (C) with `2Y=0`. The only quartic vertices with matching holomorphic scalar ends and two conjugate Higgses are:

| Yukawa channel | Source quartic |
|---|---|
| `y1 y1` | `lambdaT3ExtraQ020200I1` |
| `y1 y2` | `lambdaT3` |
| `y2 y2` | `lambdaT3ExtraQ020020I1` |

v20.1 has already checked that the complete F/Y/X saved matched-EFT `LLHH` GammaCC sectors contain only these three quartics, with 204/132/36 additive terms per model; v21 checks those sources' hashes remain unchanged.

## Two extra cubic vertices: why no Weinberg contribution in this source

1. `H S S`: for any three SU(2) doublets `2⊗2⊗2`, no singlet exists. Also `2Y(H)=+1`, while a pair of B/C scalars has total twice-hypercharge in `{-2,0,+2}`, so no neutral `HSS` cubic exists.
2. `H F F`: `F` is an SU(2) singlet (B) or triplet (C); two integer-isospin F's cannot combine with a Higgs doublet to make an SU(2) singlet. The F fields also have `Y=0`, so `HFF` is hypercharge-forbidden.
3. `H` plus one SM and one BSM fermion: for F singlet or triplet the `H l F` seesaw Yukawa is **gauge invariant**, but is excluded by the defining T3 `Z2` assumption and verified **absent from the archived UV interaction inventory**. This model assumption is essential: without it the leading neutrino-mass topology need not be a loop.
4. Gauge–Higgs cubic `H† H V` contains **one H and one H†**, not the holomorphic `HH`/`H†H†` pair; gauge–scalar/gauge–fermion cubics lack an external H. SM Yukawa `H f f` cubics do not supply BSM scalar ends to connect to the two T3 Yukawa vertices, at the available four-vertex budget. Unbroken SU(2)×U(1) has no pure-Higgs cubic.

Thus the two-extra-cubic graph class is excluded by explicit gauge representation, Z2, Higgs orientation and Wick connectivity restrictions. This is a conditional proof at the stated order, **not** an automatic full renormalized EFT coefficient.

## Mixed quadratic limitation

The general B/C operator manifest also admits `mu2_extra_000101_I1` (counts `S1dag S2dag` plus conjugate). A **nonzero** value changes propagators already at tree level, so arbitrarily many mass insertions are permitted without changing the loop count. It must then be included in the mass/eigenbasis and the original triangle weights **recalculated**. The **archived mass-eigenstate Matchete UV interaction lists contain only 11 baseline plus 7 extra quartic entries; no `ExtraQuadratic` or extra Yukawa**. v21 validates this for both models and uses the earlier zero-extra-quadratic boundary assumption. The result does not apply to arbitrary nonzero `mu2_extra` without a new matching calculation. The existing scalar eigenbasis may still incorporate other mixing already present in its source mass parameters; this is not an assertion that the scalar mixing matrix is diagonal.

## Evidence and safety gates

v21 requires: the complete 16-model operator inventory and source-verified Wick audit; v20.1 restricted source-census report with its exact F/Y/X source hashes; unchanged current UV files; the v19 restricted hard-soft subtraction report; and the archived B/C mass-eigenstate source lists and SU(2) tensor checks. No Matchete run or GitHub writes occur.

`--require-physical` fails deliberately. A successful result certifies the absence of extra renormalizable one-loop **graph topologies** within the archived source, not full finite operator/parameter/EOM matching, RG transport or the physical Wilson coefficient.

### Next physics task

With the diagram classes closed for the archived benchmark, the remaining work is the scheme-consistent UV→EFT1→EFT2 finite operator matching (including the original-vs-canonical EFT1 equality at the operator level), followed by scale-appropriate RG transport and benchmark evaluation. Avoid further duplicative loop-topology audits unless the UV source changes.
