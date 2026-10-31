# Workstream C v5 — independent scalar-bubble UV counterterm (experimental)

## Validated input and precise scope

The user's latest **v4** local log states:

- `B_-1: pole matches=True; formal finite matches XY=True; YY/XY=18/32; full-mass bubble F_XY=0.131434492`
- `C_-1: pole matches=True; formal finite matches XY=True; YY/XY=18/32; full-mass bubble F_XY=0.132107353`

The v3 native SU(2) check previously found exactly 50 coefficient ratios of `1` for each benchmark. This v5 audit is new: derive a **UV counterterm from the two-scalar bubble loop integral and the original EFT1 tree-vertex product**, then verify cancellation of the native pole. Do not construct the CT as `-P_native` and mistake the identity for an independent renormalization derivation.

There are only two source scalar pair classes that contain a Y propagator: `YY` and `XY`. A pure `XX` loop belongs to the lower-energy EFT side of matching. The input and output normalization is the Matchete ordered **neutral** `LLHH`, `GammaCC/Proj` coefficient, with common `hbar=1/(16 pi^2)` and `1/MF` factored out by the pre-existing v3 parser. The result is **not** the symmetrized two-component Weinberg coefficient.

## Independent UV calculation

In Euclidean signature (suppressing overall vertex/sign conventions already kept in the native W-Q contraction), the massive-bubble radial integral is

\[
J(L;a,b)=\int_0^L\frac{t\,dt}{(t+a)(t+b)}, \qquad a=m_a^2,\ b=m_b^2,\ L=\Lambda^2.
\]

The coefficient of the UV logarithm follows directly:

\[
\lim_{L\to\infty}\frac{dJ}{d\ln L}
=\lim_{L\to\infty}\frac{L^2}{(L+a)(L+b)}=1.
\]

The UV residue is **independent of the two masses**. The infrared endpoints are finite for both relevant cases:

- YY: `t/(t+M_Y^2)^2 -> 0` as `t -> 0`;
- XY, keeping leading dimension-five term `M_X -> 0`: `t/((t+M_Y^2)*t) -> 1/M_Y^2`.

Consequently, **the one-loop YY/XY bubble itself has a UV pole but no IR pole** at zero external momenta. This does not classify every divergence in the full Matchete/EFT1/EFT2 matching calculation.

In conventional dimensional regularization and MS-bar, the scalar integral has, after extracting `hbar`,

\[
I_{ab}^{\overline{\mathrm{MS}}}\sim
\frac{1}{\bar\epsilon}+ \int_0^1dx\log\frac{\mu^2}{xm_a^2+(1-x)m_b^2}
\quad (+\,O(\epsilon)),
\]

in the same overall sign and normalization as the independently checked native Wick contraction. At `mu=M_Y`, the **strict dimension-five/light-X** limit gives

\[
I_{YY}^{\rm fin}=0,\qquad I_{XY}^{\rm fin}|_{M_X\to0}=1.
\]

If instead `M_X` is kept **fully massive**, the complete zero-momentum bubble is

\[
I_{XY}^{\rm fin}|_{\mu=M_Y}=1+\frac{r\ln r}{1-r},\quad r=M_X^2/M_Y^2.
\]

The latter is `~0.13` for B/C, as the v4 diagnostic found. It is **not** interchangeable with the strict dimension-five hard matching piece: its terms in powers of `M_X^2/M_Y^2` involve soft-region and higher-dimensional EFT information. A sharp Euclidean cutoff fixes the universal UV-log residue; it **does not** fix the MS-bar finite constant (the finite scheme is established separately by the dimensional-regularization Feynman-parameter expression above).

## Derived restricted counterterm

Compute the contraction **without using the native pole**:

\[
K_Y=K_{YY}^{LLSS\times HHSS}+K_{XY}^{LLSS\times HHSS},
\]

where the native SU(2) CG tensors, Weyl/Dirac projector structure, and Wick pairing multiplicities are provided by the already verified v3/v4 code. In the **native ordered neutral-component operator convention**, an independently derived MS-bar counterterm for this scalar-bubble contribution is

\[
\boxed{\delta\mathcal L_{LLHH}\big|_{\lambda\cdot LLSS,Y}
=-\frac{hbar}{M_F\bar\epsilon}K_Y\,\mathcal O_{LLHH}^{\rm native}}
\]

subject to the usual interpretation that `K_Y` is a set of flavour/ordered-component coefficient monomials, not a single scalar number. The script computes `-K_Y` **from source W-Q vertices only**, then independently checks it against the archived Matchete Y pole. It does not infer a CT by taking the negative of the saved hard pole. The native pole is only an external cross-check.

Within this **bubble-only scalar-quartic sector**, the subtracted finite expression at `mu=M_Y` and the `m_X -> 0` dim-five projection is `hbar*K_XY/M_F`, because YY contributes zero and XY contributes one. We call this the **conditional scalar-bubble-sector finite operator**, not a physical `C5`.

Remaining unresolved before a physical Y-threshold coefficient can be certified:

1. Full EFT1 and EFT2 dimension-five operator renormalization, counterterm **difference** and any additional field/operator mixing;
2. EFT2 soft-overlap and its higher-dimensional `LLHHXX` analogues, if finite `M_X` effects are to be kept;
3. Precise dimensional-regulator/MS-bar conventions of the installed Matchete version beyond this checked bubble subset;
4. `GammaCC`/Majorana to flavour-symmetric Weinberg `C5` normalization;
5. EFT1 running to `muY`, X-threshold inheritance, then final `C5`.

## Project-relative installation

Extract the ZIP into repository root; it installs **two complete new files** (no existing script replaced):

- `scripts/DeriveT3BCYBubbleMSbarCT.experimental.py`
- `docs/WorkstreamC_YBubbleMSbarCounterterm_README.md`

Requires pre-existing v2/v3/v4 scripts and their **original local native data**. It does not write or alter GitHub or production matching code.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive `
  -LiteralPath "$HOME\Downloads\T3_WorkstreamC_YBubbleMSbarCT_v5.zip" `
  -DestinationPath . -Force

python scripts\DeriveT3BCYBubbleMSbarCT.experimental.py --self-test
python scripts\DeriveT3BCYBubbleMSbarCT.experimental.py `
  --output output\audits\t3_bc_y_scalar_bubble_msbar_v1\summary.json
```

Expect for both models **18 YY / 32 XY**, independent bubble-CT cancellation against the 50 native pole coefficient monomials, and exact XY finite agreement. This is **not** proof of full EFT1/EFT2 renormalization or physical `C5`.

Diagnostic status and anti-overclaim flags are included with all 50 CT coefficients and 32 conditional finite coefficients in the generated JSON. A changed source hash, changed tensor, altered prior v4 result or coefficient mismatch causes failure.
