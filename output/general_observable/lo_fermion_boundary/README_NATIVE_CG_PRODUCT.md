# A_-4 native SU(2) contraction check

**Intended project directory:** `output/general_observable/lo_fermion_boundary/`

Input: output from `native_cg_a_m4_probe.wl`, already generated locally under
`output/general_observable/lo_fermion_boundary/normalization_probes/native_cg_a_m4_<uuid>.json`.

Run from the repository root in PowerShell:

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$dir = 'output\general_observable\lo_fermion_boundary'
python -X utf8 "$dir\test_native_cg_a_m4.py"
$cg = Get-ChildItem "$dir\normalization_probes\native_cg_a_m4_*.json" |
    Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName
python -X utf8 "$dir\check_native_cg_a_m4.py" --report $cg
```

The diagnostic reads, but does not change, the source CG report. It writes a **new** JSON
file `normalization_probes/native_cg_product_a_m4_<sha16>.json`, refusing overwrites.

## What is actually calculated

The report is a numerical component-level test of the expected native holomorphic A_-4
SU(2) contraction. The two Yukawa contractions in the installed UV Lagrangian are

- `Bar[T3Y1CG]` with external lepton `i` and internal heavy doublet `k`;
- `Bar[T3Y2CG]` with external lepton `j`, the conjugate heavy doublet `k`, and scalar-triplet `u`.

On the nominal heavy-fermion `delta[k,k']` contraction, the diagnostic computes

```
P[i,j,u] = sum_k conj(Y1[i,k]) * conj(Y2[j,k,u])
```

and compares it with the native Matchete matched `T3MixCG[i,j,u]` using

```
k = sum(conj(T)*P) / sum(abs(T)**2)
residual = norm(P-k*T) / norm(P)
```

It also reports the swapped-lepton-index comparison and the unbarred-CG diagnostic,
but **does not automatically select the one with the smallest residual**. The primary
orientation is taken from the saved A_-4 native UV source and could require an
explicit Matchete index-convention clarification if the result is not a unit match.

Passing requires `k=+1` and a relative residual <= `1e-10`. Other outcomes exit
with code 2, preserving the diagnostic JSON for review; such an outcome does not
by itself prove the physical EFT coefficient is wrong.

## Certification limits

This checks SU(2) component contractions **only**, separately from the already-run
signed p=0 Grassmann spinor check. It does not establish the unprovided absolute
`GammaCC` phase relative to the master Weyl bilinear. It does not alter, rescale,
or overwrite any of the 16-model boundary releases, or claim physical certification.
