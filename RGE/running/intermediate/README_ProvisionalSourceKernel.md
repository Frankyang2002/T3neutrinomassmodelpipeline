# All-16 provisional F-tree LLSS → quartic mixing kernel

**Paths owned by this workstream:**

- `RGE/running/intermediate/ProvisionalSourceKernel.py`
- `RGE/running/intermediate/README_ProvisionalSourceKernel.md`
- `scripts/RunProvisionalSourceKernel.py`
- `tests/test_provisional_source_kernel.py`

These are all **new** files. No Matchete, F matching, scalar threshold, existing RGE or previously generated output is modified.

## Why this intermediate stage exists

As verified from public `main` commit `e5e4a32f6f5a96eb92ac4fc2c9884127e56a2ce4` (2026-11-02), Chat 1 has now published source-resolved *numerical* LLSS boundary coordinates for 16/16 models under

`output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/`.

All 16 are still explicitly marked `NumericalLLSSBySourceVerified_SpinorNormalizationOpen`. The `GammaCC` → canonical left-Weyl sign/factor remains `null`; the recorded local Matchete probe did **not** load Matchete. Likewise, a complete **physical** EFT1 scalar-quartic tensor with certified Matchete/RGBeta conventions has not been handed to this workstream. The existing certified-only mixer must continue to reject these inputs.

This bridge reads the **actual source-resolved** LLSS coefficients to generate the *linear response* of the Weinberg **pivot** to every independent real-scalar quartic tensor component. It uses the correct source component coordinates and their 1/2/4 distinct permutation orbits without guessing a normalization factor. Unlike synthetic fixtures, these coefficients come from Chat 1's saved source calculations, but the output remains *conditional*, not physical.

## Derivation and precision

With conventions from `RGE/general/WilsonTensorRGE.py` and `RGE/matching/WeinbergTensorAdapter.py`:

\[
 V_4 = \frac{1}{4!}\lambda_{abcd}\phi_a\phi_b\phi_c\phi_d,\quad
 \mathcal L_5=C_{ijab}\psi_i\psi_j\phi_a\phi_b+\mathrm{h.c.},
 \quad C_{\nu_p\nu_q,r,r}=\frac14 C_5^{pq}.
\]

The scalar-pair quartic part of the master RGE, for \(r=\operatorname{Re}H^0\), gives

\[
 16\pi^2\frac{d C_5^{pq}}{d\ln\mu}\bigg|_{\lambda,\mathrm{pivot}}
   = 4\sum_{a,b\in\mathrm{heavy}}\lambda_{rrab}C_{\nu_p\nu_q ab}
   = \sum_{a\le b} K_{ab}^{pq}\lambda_{rrab},\qquad
 K_{ab}^{pq}=4(2-\delta_{ab})C^{\rm source}_{\nu_p\nu_q ab}.
\]

The implicit **unknown** common `GammaCC`/Weyl normalization \(N_\Gamma\) must multiply the entire RHS. This is **not** a full Weinberg-subspace projection: the other external Higgs components and SU(2) Ward residual need a *complete, certified* quartic tensor. No quartic coefficient is set to a guessed value; no numerical physical \(\beta_{C_5}\) or \(\Delta C_5\) is claimed. The fixed-order scale dependence, once normalized couplings are supplied, is

\[
 \Delta C_5^{\rm EFT1,log}=\frac{\ln(\mu_S/\mu_F)}{16\pi^2}\,\beta_{C_5}^{(1)}(\mu_F).
\]

Every original PL Matchete source and every heavy-scalar pair contributes independently, even when a source vanishes at a selected benchmark. Mixed B/C and split-real A/D/E sources are preserved; no dominant-channel filter is applied. The program checks each component's orbit, summed source tensor, flavour symmetry, and source-by-quartic recomposition before writing anything.

The other six master beta terms are retained in `GeneralEFT1C5Mixing.py`; this small bridge **only** derives the scalar-pair contribution's coupling-independent response. Complete leading-one-loop operator-class closure beyond matched `psi^2 phi^2` has not been established. The joint-scalar threshold belongs to the other workstream.

## Commands — PowerShell from repository root

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
python -m unittest discover -s tests -p 'test_provisional_source_kernel.py' -v

# This is a read-only scan of the immutable Chat 1 source release.
python scripts\RunProvisionalSourceKernel.py `
  --release-root 'output\general_observable\lo_fermion_boundary\releases\source_021dcc1e8fe8ba0f' `
  --output-root 'output\general_observable\lo_eft1_mixing\source_kernel_v1'
```

Expected outputs:

```text
output/general_observable/lo_eft1_mixing/source_kernel_v1/
  summary.json
  A_m4/quartic_pivot_kernel.json
  A_m2/quartic_pivot_kernel.json
  ... (all 16 tags)
```

Outputs are **new** and non-overwriting. The script loads **all** inputs and validates source sums before writing any files. Run to a distinct output directory if `source_kernel_v1` already exists. If Chat 1 publishes a newer source release with different flags, the program deliberately fails closed until the new format is reviewed.

**Next blocking physics task:** Independently certify Matchete `GammaCC` spinor mapping and supply the source-consistent full EFT1 quartic tensor. Then connect the normalised F-tree boundary to the already tested 7-term EFT1 mixer and its full SU(2) Ward check. **Do not** promote the conditional kernel to a physical C5 by inserting arbitrary quartic values or assuming `GammaCC=+1`.
