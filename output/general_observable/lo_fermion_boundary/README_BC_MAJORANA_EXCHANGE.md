# B_-1 / C_-1 Majorana signed tree LLSS audit (conditional)

**Installation paths** (all three files in this archive):

- `output/general_observable/lo_fermion_boundary/bc_majorana_exchange_audit.py`
- `output/general_observable/lo_fermion_boundary/test_bc_majorana_exchange_audit.py`
- `output/general_observable/lo_fermion_boundary/README_BC_MAJORANA_EXCHANGE.md`

No code is changed outside this folder. The existing source-resolved 16-model release is **read only**.

## Verified operations

1. Use the existing exact rational exterior algebra (`signed_zero_momentum_check.G`) with left Grassmann differentiation to solve the heavy Majorana field equations at zero momentum for **four independent source currents**. In the explicitly chosen convention

   \(\epsilon_{12}=+1,\quad \mathcal L_{\rm UV}= -\frac M2 h^T\epsilon h+\sum_{a=1}^4 y_a h^T\epsilon \ell_a\),

   the exact result is

   \(\mathcal L_{\rm EFT}^{(0)} =\sum_a\frac{y_a^2}{2M}(\ell_a\ell_a)+\sum_{a<b}\frac{y_a y_b}{M}(\ell_a\ell_b)\).

   This is a **conditional spinor identity**, not a determination of Matchete's charge-conjugation phase.

2. Recompute the B_-1 and C_-1 source provenance against the immutable release `releases/source_021dcc1e8fe8ba0f`, checking all **10 B_-1 PL** and **12 C_-1 PL** source terms, their barred couplings, mixing spurions, X/Y scalar legs, explicit signs, literal mass denominators, and native CG families; exclude the 10/12 PR Hermitian conjugates.

3. Keep the **C_-1 bridge contractions at source indices 23 and 24** separate. They have \(\pm \sqrt3/(2M)\) *times their bridge CG structures*, not generic \(\pm 1/(2M)\) coefficients. Their net effect can only be decided by an independently verified CG contraction.

## Not verified here

This does **not** extract B/C mixed-eigenfield Feynman rules, contract a full native UV four-point amplitude, check native SU(2) CG identities, or fix `GammaCC` relative to the master's Weyl spinor convention. Do **not** promote the 16-model LLSS boundary to a physical Wilson normalization on the basis of this diagnostic.

## Run in PowerShell at repository root

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$dir='output\general_observable\lo_fermion_boundary'
python -X utf8 "$dir\test_bc_majorana_exchange_audit.py"
python -X utf8 "$dir\bc_majorana_exchange_audit.py" --offline
python -X utf8 "$dir\bc_majorana_exchange_audit.py" --check-only
python -X utf8 "$dir\bc_majorana_exchange_audit.py"
```

The final command writes a deterministic, idempotent report to `output/general_observable/lo_fermion_boundary/normalization_probes/bc_majorana_exchange_<SHA16>.json`. If an identical file already exists, it reports `UNCHANGED`, avoiding the previous probe's unnecessary `FileExistsError`. To prevent any output, use `--check-only`.

Next certification step: calculate B/C native mass-eigenfield Feynman-rule UV amplitudes and evaluate their SU(2) CG contractions before fixing the absolute master-Weyl phase.
