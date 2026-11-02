# T3 simultaneous scalar-threshold matching to SMEFT C5

## Scope and status

**Exact threshold sequence implemented by this isolated module:**

`UV -> [all F at mu_F] -> EFT1(H,S1,S2; C_LLSS and C5) -> [all real components of S1 and S2 together at mu_S] -> SMEFT C5`.

There is **one scalar threshold, not two**, for **all 16 T3 models**. The module does not modify or launch fermion matching, EFT1 RGE flow or SMEFT RG flow. It reads the already-transported EFT1 tensors at **the same mu_S**. Real mass-matrix diagonalization treats normal, self-bilinear/split-real and mixed B(-1)/C(-1) scalar sectors in the same contraction.

**Verified implementation:** the complete one-insertion `C_LLSS x quartic_HHSS` renormalized zero-momentum scalar-pair bubble in the dimension-five, unbroken-EW and leading `1/mu_F` truncation, including *all* scalar mass eigenstate pairs, finite MSbar B0, full complex 3x3 flavour, independent quartic contributions supplied in the tensor, and the inherited direct `C5_EFT1(mu_S)` added exactly once. No dominant-contribution filtering occurs. Generic quadratic mass-matrix diagonalization and tensor rotations are implemented. Numerical unit tests are run without Wolfram.

**Not independently verified:** the complete physical scalar threshold of an arbitrary full EFT1 Lagrangian, including external-field threshold renormalization, independently computed finite counterterms and possible additional EFT1 operators, and the Matchete/Chat-1 CG and flavour normalization. Those steps were already marked unresolved in the source repository's B/C scalar-bubble diagnostic. No physical neutrino-mass prediction or certified full one-loop SMEFT matching is claimed by the fixture results. `physical_full_one_loop_certified` remains `false`. The source files do not automatically provide a Chat-1 normalized `(3,3,2,2,n,n)` tensor for every model.

### Audited reference code (GitHub `main` commit `c4788e1e6e22341385b155ae06cec664e740204a`, read only)

- `Lagrangian/model/T3ModelCatalog.wl`: `(dS1,dS2,dF)` mapping for classes A-E.
- `Lagrangian/interactions/ScalarPotential.wl`, `Lagrangian/interactions/GeneralScalarOperators.wl`: baseline and extended independent quartics; extra quadratic terms require real mass treatment.
- `RGE/group_factors/core/QuarticTensorAlgebra.py`, `RGE/running/intermediate/ScalarOnlyTensorAdapters.py`: fully symmetric real-scalar potential convention and EFT1 Wilson representation.
- `RGE/running/intermediate/ScalarOnlyWilsonTensorRGE.py`: scalar quartic contribution to the EFT1 Wilson beta.
- `Numerical/experimental/T3BCGroupedFiniteScalarBubble.py`, `Numerical/experimental/T3BCScalarBubbleThreshold.py`: existing, **experimental**, bubble formula and explicit finite-matching caveats.
- `Numerical/experimental/GeneralT3MixedScalarMassBasis.py`: B/C doublet mixed-basis convention. The new matcher uses the equivalent all-real orthogonal mass basis supplied by Chat 1, *without* replacing that module.
- `docs/T3_COMPLETE_GENERAL_SCALAR_SYMBOLIC.md`: six self-bilinear cases and B/C mixed cases; physical full-scalar integration was not certified.

## Analytic matching

Take canonical, real heavy scalars \(s_a\) (all components of both multiplets) and real Higgs coordinates \(h_A\), with

\[
 V_4=\frac1{4!}\lambda_{ABab}h_Ah_Bs_as_b+\cdots,\quad
 16\pi^2\,\beta(C_{LLhh})\big|_{\lambda}=\lambda_{hhss}C_{LLss}.
\]

Here the displayed \(\lambda_{ABab}\) is the **component of the fully symmetric four-real-scalar potential tensor**, not the coefficient of the quartic monomial before combinatoric symmetrization. Let \(R^T M_S^2 R=\mathrm{diag}(m_a^2)\) with \(s_{\rm old}=R s_{\rm eigen}\). Transform the *entire* EFT1 Wilson tensor and the quartic simultaneously:

\[
 \widehat C^{fgij}_{ab}=C^{fgij}_{cd}R_{ca}R_{db},\qquad
 \widehat\lambda_{ABab}=\lambda_{ABcd}R_{ca}R_{db}.
\]

The one-loop **joint** scalar threshold in MSbar at \(\mu_S\) is

\[
 \Delta\mathcal C^{fgij}_{AB}(\mu_S)=-\frac1{2(16\pi^2)}
 \sum_{a,b=1}^{2d_{S1}+2d_{S2}}
 \widehat C^{fgij}_{ab}\widehat\lambda_{ABab}
 \overline B_0(0;m_a^2,m_b^2;\mu_S^2),
\]

\[
 \overline B_0(0;x,y;\mu^2)=1-\frac{x\ln(x/\mu^2)-y\ln(y/\mu^2)}{x-y},
 \quad\overline B_0(0;x,x;\mu^2)=-\ln(x/\mu^2).
\]

The only scalar-potential quartics contributing directly in this one-vertex dimension-five bubble are those with two external Higgs legs and two heavy scalar legs; quartics with no Higgs legs enter through upstream masses or RG running, not by inserting an extra diagram at this order. The code accepts independently labelled quartic sub-tensors and sums them **without filtering small entries**; the output reports both each quartic channel and each eigenstate pair. Each individual projected channel can have a non-Weinberg residual even when the complete sum is Weinberg-invariant. The code sums **unordered eigenpairs** with multiplicity \(1\) (equal indices) or \(2\) (distinct indices), reproducing the ordered sum exactly. The finite constant `+1`, unequal-mass logarithms and smooth degeneracies are retained. Splittings and mass mixing enter through \(M_S^2\) and transformed tensors, not through model-specific corrections applied after projection.

We use \(H=(h_0+ih_1,h_2+ih_3)/\sqrt2\), \(\epsilon_{01}=+1\), and define the Weinberg external tensor \(T_{ijAB}\) by

\[
 (L^i\epsilon_{ij}H^j)(L^k\epsilon_{k\ell}H^\ell)
 = (L^i L^k)\,T_{ikAB}h_Ah_B.
\]

The projection \(\Delta C_5^{fg}=\langle T,\Delta\mathcal C^{fg}\rangle/\langle T,T\rangle\) is checked by reconstructing the full LLHH external tensor. The relative norm of the orthogonal residual must be below `1e-8` unless the user explicitly selects the **diagnostic-only** override. With input \(C_5^{\rm EFT1}(\mu_S)\) already renormalized, the returned total is

\[
 \boxed{C_5^{\rm SMEFT}(\mu_S)=C_5^{\rm EFT1}(\mu_S)+\Delta C_5^{\rm joint\ scalar}(\mu_S).}
\]

No complex conjugation of either flavour Wilson coefficient or the real orthogonal mixing matrix is applied in this holomorphic contraction. This formula **assumes the input Lorentz bilinear, lepton flavour factors, quartic Feynman-rule normalization and Weinberg operator convention match the declared adapter**. The RGE sign fixes the logarithmic piece relative to the source project, not all possible source convention ambiguities.

### Matching logs vs EFT1 RG logs

For any *common* fixed reference scale \(\mu_*\):

\[
 \overline B_0(\mu_S)=\overline B_0(\mu_*)+2\ln(\mu_S/\mu_*).
\]

The output labels the two terms `scalar_bubble_finite_at_reference` and `scalar_bubble_explicit_matching_scale_log`. **They are a partition of the same single threshold term, not additional EFT running contributions.** \(C_5^{\rm EFT1}(\mu_S)\) is supplied by the upstream EFT1 flow and contains its own \(\ln(\mu_F/\mu_S)\) effects. The joint matcher adds **no** such log again. At fixed couplings/masses, the test verifies

\[
 \frac{d\Delta\mathcal C_{LLhh}}{d\ln\mu_S}
 =-\frac{\lambda_{hhss}C_{LLss}}{16\pi^2},
\]

which cancels the dropped heavy-scalar contribution to the EFT1 anomalous dimension. The complete running-scale invariance of a physical observable requires independent verification of all EFT1 and SMEFT sectors.

### EFT truncation and mass splitting

- The calculation is at zero external momentum, O(C_LLSS lambda) and dimension five, leading in \(1/M_F\). Dimension-seven terms from the fermion expansion, non-scalar-bubble loops or omitted scalar-dependent wavefunction effects are not calculated.
- `scalar_mass_origin=unbroken_electroweak` is the clean dimension-five interpretation. If a mass splitting comes only from EWSB, input `EWSB_resummed`: its exact eigenmass dependence can resum \(v^2/M_S^2\) higher-dimension terms. A nonzero Weinberg-projection residual is an explicit warning; suppressing it does not certify SMEFT matching.
- The full *unbroken-EW* quadratic mass matrix may contain gauge-invariant self-bilinear splitting or B/C mixing. The canonical, real mass matrix must be constructed by Chat 1 with correct epsilon and SU(2) component ordering. Supplying arbitrary numerical mass matrices does not by itself establish gauge invariance.
- Gauge-dependent finite parts, external-field normalization, counterterm scheme conversions, and possible other independent operators must be derived separately if present in the chosen EFT1 truncation. For the restricted LLSS x quartic sector, the bubble is complete at the specified coupling order.

## 16 models and scalar cases

| Class | (dS1,dS2,dF) | alpha values |
|---|---|---|
| A | (1,3,2) | -4, -2, 0 |
| B | (2,2,1) | -3, -1, +1 |
| C | (2,2,3) | -3, -1, +1 |
| D | (3,1,2) | -2, 0, +2 |
| E | (3,3,2) | -4, -2, 0, +2 |

- Self-bilinear/split-real: `A_+0`, `D_+0`, `E_+0` (S1) and `A_-2`, `D_-2`, `E_-2` (S2).
- Two-scalar mixed quadratic sector: `B_-1`, `C_-1`.
- Other eight model points: ordinary (no additional quadratic mixing).

## Required Chat 1 normalized adapter contract

A pair of files is required **for each model**, e.g. `A_m2.json` and `A_m2.npz`. The CLI scans metadata JSON by its *model value*, not filename, and rejects missing/duplicate models. Use only the following axes, all in the **canonical real scalar component basis**, with original order **all real S1 components, then all real S2 components**. Every tensor must be at the same \(\mu_S\), in the same `MSbar` scheme:

| NPZ key | dtype | array shape | meaning |
|---|---|---|---|
| `mass_squared_gev2` | real float | `(n,n)` | symmetric, positive real mass-squared matrix, `n=2*dS1+2*dS2`; *includes every scalar bilinear* |
| `C_LLSS` | complex | `(3,3,2,2,n,n)` | `(flavour f,flavour g,lepton SU2 i,j,real scalar a,b)`; includes **all independent Wilson contractions**, has symmetry `(f,i)<->(g,j)` and `a<->b` |
| `lambda_HHSS` | real float | `(4,4,n,n)` | full symmetric potential tensor `lambda[A,B,a,b]` with **all independent Higgs/scalar quartics and their relative signs/CGs already summed**, not just `lambdaT3` |
| `C5_EFT1` | complex | `(3,3)` | inherited **already-run** SMEFT-holomorphic Weinberg coefficient at \(\mu_S\); set explicit all-zero array if absent, do not copy first-threshold value at \(\mu_F\) |
| `lambda_HHSS_channels` | real float | `(k,4,4,n,n)` | **optional but recommended** per-independent-quartic source-normalized tensor; must sum to `lambda_HHSS`; no channel is discarded |

Metadata example (JSON; replace the values by actual matched ones):

```json
{
  "schema_version": 1,
  "model": "B_-1",
  "tensor_convention": "T3_REAL_LLSS_HHSS_WEINBERG_v1",
  "renormalization_scheme": "MSbar",
  "mu_f_gev": 5000.0,
  "mu_s_gev": 1000.0,
  "reference_scale_gev": 1100.0,
  "scalar_mass_origin": "unbroken_electroweak",
  "tensor_npz": "B_m1.npz",
  "quartic_channel_names": ["lambdaT3", "lambdaExtra1"],
  "synthetic_fixture": false,
  "provenance": {"source": "Chat 1 normalized tensor export (specify snapshot and conventions)"}
}
```

The illustrated `quartic_channel_names` is an example: its length must **exactly equal** the first dimension of `lambda_HHSS_channels` in the NPZ, and its entries must be unique. If per-channel tensors are unavailable, the whole `lambda_HHSS` is reported as `aggregated_quartics_independence_not_certified`.

**No automatic import** from raw Matchete coupling names, RGBeta symbolic exports, benchmark Yukawas or existing hard-mass formulae is attempted. Chat 1 must verify and supply the above *normalized component tensors*. Without those the generated 16-model results are **synthetic stress-test fixtures only**, not physical values. The adapter fails on shape/asymmetry/nonfinite data, missing provenance, unspecified renormalization convention and nonpositive eigenmasses. Mixed B/C complex scalar rotations must be fully realified by the upstream exporter; a complex `mass_squared_gev2` in a real-component basis is rejected.

## Full installation and exact PowerShell commands

From PowerShell (the ZIP preserves project-relative paths):

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_joint_scalar_matching_v1.zip' -DestinationPath . -Force

# Offline algebra checks (no Wolfram, no existing pipeline re-run):
python -m unittest discover -s tests -p 'test_t3_joint_scalar_matching.py' -v

# Regenerate all 16 SYNTHETIC-only JSON matrices and NPZ input fixtures:
python -m scripts.RunT3JointScalarMatching --fixtures

# Check the single-model adapter on one synthetic normalized fixture:
python -m scripts.RunT3JointScalarMatching `
    --input-dir 'output\general_observable\lo_joint_scalar_matching\fixtures_demo\fixtures' `
    --model 'B_-1' `
    --output-dir 'output\general_observable\lo_joint_scalar_matching\adapter_demo'

# After Chat 1 exports correctly normalized *physical* EFT1 inputs for all 16:
python -m scripts.RunT3JointScalarMatching --input-dir 'C:\path\to\chat1_joint_scalar_inputs'

# Examine results:
Get-Content .\output\general_observable\lo_joint_scalar_matching\fixtures_demo\summary.json
Get-Content .\output\general_observable\lo_joint_scalar_matching\fixtures_demo\B_m1_c5.json
```

Requires Python 3.10+ and NumPy. Importing the package uses no `Matchete`, `RGBeta` or `Wolfram` dependency. All **created outputs** are under `output/general_observable/lo_joint_scalar_matching/` only. Real-data results use per-model `*_c5.json` and one `summary.json` directly in the allowed output root. The `--fixtures` flag defaults to `output/general_observable/lo_joint_scalar_matching/fixtures_demo/`, with synthetic tensor inputs in its `fixtures/` child, so synthetic results **cannot overwrite real-data outputs by default**. The generated archive additionally includes the precomputed synthetic outputs for inspection. Synthetic fixtures are segregated from real outputs by default. If providing `--output-dir` manually, choose separate non-overlapping subdirectories for experimental and physical runs.

**Files in this delivery:**

```
Numerical/joint_scalar_matching/__init__.py
Numerical/joint_scalar_matching/core.py
Numerical/joint_scalar_matching/adapter.py
Numerical/joint_scalar_matching/fixtures.py
Numerical/joint_scalar_matching/README.md
scripts/RunT3JointScalarMatching.py
tests/test_t3_joint_scalar_matching.py
output/general_observable/lo_joint_scalar_matching/fixtures_demo/  # synthetic demonstration only
```

No replacement or modification of the repository's existing EFT1-mixing, first fermion-threshold, scalar matching, fitter, or any prior numerical output code is required. There are no GitHub writes.
