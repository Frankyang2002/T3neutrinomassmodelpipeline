# T3 Batch 3 — simultaneous heavy-scalar one-loop Higgs-potential threshold (16 models)

## Scope and upstream data

**Required:** all 16 models completed in
`output/parallel_t3/strict_one_loop_all16_batch2_v1/`.

This batch never changes Batch 1/2 or the saved Workstream B inputs. It evaluates
the full source-conditional real-scalar Hessian and Higgs-portal quartic tensors
already constructed in Batch 2, in the **unbroken electroweak/Z2 phase**.
Both scalar multiplets, including the B/C mixed-doublet and all self-split
neutral scalar modes, are integrated out at ONE common `MS = 1e3 GeV`.

**Actual new one-loop matching:** the Coleman–Weinberg contribution to the
renormalisable SM Higgs quadratic and quartic potential. This is not another
one-loop Weinberg matching graph; the existing strict Weinberg C5 at MS remains
bitwise unchanged. Source normalization is conditional (as in Batch 2).

## Mathematics and normalization

For heavy **real** components `s_a` and four Higgs real components `h_i`:

\[
V \supset \frac12 s_a (M_S^2)_{ab}s_b
 +\frac14 (\lambda_{HHSS})_{ijab}h_i h_j s_a s_b.
\]

The scalar mass operator is

\[
\mathcal M^2(h) = M_S^2 +\tfrac12\lambda_{ijab}h_i h_j.
\]

The one-loop heavy-scalar contribution to the effective potential is

\[
\Delta V(H)=\frac{1}{64\pi^2}\operatorname{Tr}\left[
\mathcal M^4(H)\left(\ln\frac{\mathcal M^2(H)}{M_S^2}
-\frac32\right)\right].
\]

With `m_a^2` the complete heavy-scalar mass eigenvalues, and `Q^r_ijab`
individual named HHSS quartic channels in that eigenbasis,

\[
\Delta(m_H^2)^{(r)}_{ij}
=\frac{1}{64\pi^2}\sum_a 2m_a^2
\left(\ln\frac{m_a^2}{M_S^2}-1\right) Q^{(r)}_{ijaa}.
\]

The complete fourth-derivative tensor is

\[
\Delta\lambda_{ijkl}^{(r,s)}
=\frac{24}{512\pi^2}\,\mathrm{Sym}_{ijkl}
\sum_{a,b}F_{ab}Q^{(r)}_{ijab}Q^{(s)}_{klba},
\]

\[
F_{ab}=\begin{cases}
\displaystyle\frac{2m_a^2(\ln(m_a^2/M_S^2)-1)-2m_b^2(\ln(m_b^2/M_S^2)-1)}{m_a^2-m_b^2},&a\ne b,\\
2\ln(m_a^2/M_S^2),&a=b.
\end{cases}
\]

The implementation uses a stable divided difference for nearly degenerate
masses and keeps all ordered `(r,s)` interference combinations. The projected
`delta_lambdaH_source_convention` uses the already selected real potential
normalization `V4 = lambdaH * (Hdag H)^2 / 2`. SU2 invariance is checked via
the `delta(m_H^2) \propto I_4` and `delta lambda4 \propto` Higgs quartic
singlet projection residuals, without inventing missing gauge thresholds.

### One-loop counting limitation

The strict Batch 2 result is an **O(hbar)** calculation for *included* C5
topologies. The CW quartic threshold itself is also **O(hbar)**, but feeding
`delta_lambdaH` into a **one-loop** SMEFT beta function for an existing
`C5 = O(hbar)` generates an **O(hbar^2)** effect in neutrino masses. The
runner therefore exports BOTH:

- `masses_eV_Batch2_strict_baseline`: unchanged strict Batch 2 output;
- `masses_eV_CW_feedback_higher_order`: optional subsequent one-loop SM +
  Weinberg evolution starting from `lambdaH + delta_lambdaH` at MS, clearly
  an exploratory higher-order feedback rather than a new strict O(hbar) mass.

The default conditional `--lambdaH-transfer-factor 1` uses the same Higgs
quartic-coordinate identification as the original reduced benchmark; alternate
normalizations can be tested with the explicit parameter. No claim is made
that this is a certified native Matchete/RGBeta matching convention.

## Matchete native F-only one-loop inventory

For EACH model, the runner reads the archived
`output/general_observable/lo_certification_spinor/pole_log_e515d82d118284ba_bridge01/<slug>/trios.json`
file, checks the kernel SHA256 and Laurent triple structure, inventories
holomorphic PL versus HC PR native terms and reports which extra quartic
symbols have a recorded *conditional* RGBeta dictionary. When the latter map
has non-diagonal blocks, every potentially contributing RGBeta symbol is
listed; a zipped one-to-one alias is NOT asserted.

**No raw native pole/constant/log coefficient is inserted into C5.**
Those source records still need EFT1 soft-region subtraction, counterterms,
and an absolute spinor normalization. Simply retaining the raw finite piece
would not constitute valid hard matching and could double-count Batch 1's
triangle. This inventory is a prepared numerical source ledger for the next
first-threshold matching implementation, not a one-loop Weinberg result.

## Project paths (full files)

1. `Numerical/parallel_t3/threshold/scalar_cw_batch3.py`
2. `Numerical/parallel_t3/threshold/native_source_inventory_batch3.py`
3. `scripts/RunT3ScalarCWBatch3.py`
4. `scripts/RunT3ScalarCWBatch3.ps1`
5. `tests/test_t3_scalar_cw_batch3.py`
6. `README_T3_SCALAR_CW_BATCH3.md`

Install into `C:\thesistemp\T3neutrinomassmodelpipeline`, preserving the
paths above.

## PowerShell commands

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_scalar_cw_batch3_all16.zip" `
  -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_scalar_cw_batch3.py
.\.venv\Scripts\python.exe scripts\RunT3ScalarCWBatch3.py
.\.venv\Scripts\python.exe scripts\RunT3ScalarCWBatch3.py --show
```

Resume existing run while preserving success:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ScalarCWBatch3.py --resume
```

Recompute in a new isolated output using the same saved Batch 2 inputs:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ScalarCWBatch3.py --overwrite
```

Only compute the *strict one-loop renormalisable threshold* without a
higher-order low-energy feedback ODE:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ScalarCWBatch3.py `
  --output output\parallel_t3\scalar_cw_all16_batch3_no_feedback_v1 `
  --skip-low-feedback
```

**No code writes to GitHub.** No previously generated numerical result is
modified. The output directory is
`output/parallel_t3/scalar_cw_all16_batch3_v1/`.

Per model the package writes:

- `results/<slug>/scalar_CW_one_loop_MS.npz` — complete real-matrix
  one-loop Higgs mass/quartic shifts; **all ordered quartic-channel pairs**.
- `results/<slug>/native_F_Laurent_not_injected.json` — exact native first
  threshold source coverage without invented MSbar matching.
- `results/<slug>/CW_feedback_SM_Weinberg_low.npz` — optional exploratory
  higher-order feedback result.
- `results/<slug>/CW_feedback_observables.json` — corresponding masses/PMNS.
- `results/<slug>/model_manifest.json` — all hashes, status and numerical shifts.

Across all 16 models, `summary.json` and `scalar_CW_all16.csv` contain the
model-by-model status and primary threshold shifts.

## Next batch

The outstanding one-loop **Weinberg** matching work is the native first-threshold
hard/soft/counterterm subtraction, including the extra native LLHH quartic
channels, and any non-bubble scalar-threshold dimension-five diagrams allowed
by the original Lagrangian. None is asserted to vanish without checking its
operator topology. The source-conditional six-term EFT1 ADM and spinor
normalisation also remain provisional.
