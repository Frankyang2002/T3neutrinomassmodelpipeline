# Workstream C — Native spinor/Warsaw identity (v11)

## Concrete obstruction

The B(-1)/C(-1) native Y/LLHH expressions are **not** normalized Wilson coefficients. Before promoting a pole-subtracted neutral coefficient to the Weinberg matrix, its `GammaCC` spinor chain and chiral orientation must be mapped to the exact installed Warsaw operator. In the independently extracted installed Matchete v0.6.0 source (`installed_source.json`), `SMEFT_Warsaw.m` uses `cllHH[p,r]` without an explicit 1/2 in `Lagd5`; the native raw sector has **PL/direct** `Transp[l] GammaCC Proj[-1] l` with `H H`, and **PR/conjugate** `Bar[l] GammaCC Proj[1] Transp[Bar[l]]` with `Bar[H] Bar[H]`. The earlier v9.1/v10.1 checks audited metadata but did not test the operator identity.

## What is now independently proved

`ResolveT3BCWeinbergOperator.experimental.py` uses an *exact* chiral representation, with

- `C = i gamma^2 gamma^0`, `C^T = -C`, `C^2 = -1`, `C^dagger = -C`;
- `P_L = diag(1,1,0,0)` and `P_R = diag(0,0,1,1)`;
- explicit **Grassmann anti-commutation**, flavour transposition, and Hermitian conjugation;
- source anchors for `CConj`, `GammaCC`, `Proj` and the installed Warsaw `Lagd5` definition;
- direct audits of all 132 native Y/LLHH InputForm terms per B/C model.

The physical spinor identity in this convention has **positive** relative sign:

`(l_p^T C P_L l_r)^dagger = Bar[l_r] C P_R Transp[Bar[l_p]]`.

The Majorana bilinear is symmetric under `p <-> r`. A conventional half-normalized `L5 = (1/2) C5_ij Q_ij + h.c.` therefore corresponds to `C5 = A + A^T`, provided `A` is the complete correctly subtracted coefficient of Matchete's Warsaw `Q`; `m_nu = -v^2/2 C5` in the stated mass convention. The two conjugate PR/PL sectors must NOT be added together as independent contributions.

This identifies the **spinor sign** at the representation/source level. The absolute installed-Matchete algebraic rewrite is **not certified** unless the optional live Wolfram probe returns `ExactPlusEquality`; such a result would still not establish the threshold counterterms.

## Files (project-relative)

1. `scripts/ResolveT3BCWeinbergOperator.experimental.py` — full Python implementation and audits.
2. `Lagrangian/ProbeT3BCWarsawDirectIdentity.experimental.wl` — optional native Matchete test using only built-in `SM`, no `SMEFT_Warsaw` loading.
3. `docs/WorkstreamC_WarsawDirectIdentity_README.md` — these instructions.

No production file is modified. No GitHub write.

## Commands — PowerShell from project root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_WarsawDirectIdentity_v11.zip" -DestinationPath . -Force
python scripts\ResolveT3BCWeinbergOperator.experimental.py --self-test
python scripts\ResolveT3BCWeinbergOperator.experimental.py --output output\audits\t3_bc_y_warsaw_direct_identity_v11\summary.json

$dir = 'output\audits\t3_bc_y_warsaw_direct_identity_v11'
wolframscript -file Lagrangian\ProbeT3BCWarsawDirectIdentity.experimental.wl "$dir\native_identity.json"
if ($LASTEXITCODE -ne 0) { Get-Content "$dir\native_identity.json" }
python scripts\ResolveT3BCWeinbergOperator.experimental.py --matchete-identity "$dir\native_identity.json" --output "$dir\summary_live.json"
```

The first Python execution uses the previous local `installed_source.json` and saved B/C Y raw `LLHH` files. It can succeed independently of the optional Mathematica live identity check. The final command is **deliberately fail-closed** if the native simplifier cannot prove exact plus equality; upload `native_identity.json` if that happens. No re-running the costly Y or X matching is necessary.

## Outstanding, explicitly NOT resolved

The independent canonical EFT1 has not been *algebraically* identified with the saved original fermion-first EOMEFT. The previous Y restricted scalar-bubble cancellation does **not** certify complete EFT1/EFT2 counterterms or a scheme-correct finite threshold; the inherited LLXX terms, X matching and intermediate running are also unresolved. The script never outputs a physical C5 matrix. In particular, it does not promote formal 32-term Y constants to a physical coefficient.
