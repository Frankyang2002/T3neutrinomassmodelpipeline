# B_-1 / C_-1 independent native SU(2) Yukawa-CG product audit

Project directory: `output/general_observable/lo_fermion_boundary/`.

## Scope and provenance

`native_bc_yukawa_cg_probe.wl` evaluates the same `PhysicalSU2InvariantBasis`
arguments and original-S2 epsilon transport as the current
`Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl`; it **does not** rerun the
heavy-fermion matching, modify any pipeline code, or touch certified releases.
The saved native Matchete matched-source CG registry is independently loaded
by `check_bc_native_cg_products.py` from
`output/general_observable/ordered_wilson_tensor_v4/{B_m1,C_m1}/mixed_native_tree_llss_cg_registry.json`.
The existing `bc_majorana_exchange_audit.source_audit` authenticates all
22 PL sources and their original release checksums.

The B_-1 test compares independently extracted native `T3MENativeY1CG` and
`T3MENativeY2CG` matrices, component by component, to the matrices saved in
the matched registry. The C_-1 test constructs the actual heavy-triplet
index contraction of the *barred* Yukawa CGs,

    K[i,j,a,b] = sum_m conjugate(Y[i,m,a]) conjugate(Y[j,m,b])

and compares it against the independent matched invariant tensors,

    Q[i,j,a,b] = Bar[T3ExtraCGQ020020I1][i,j,a,b]
    E[i,j,a,b] = Bar[eps[SU2L]][i,j] T3MEBridgeBB[b,a].

For each Yukawa pair 11, 12, 22, the **test target** is

    K = Q / sqrt(3) + sqrt(3) E / 2.

The tensors and their relative phase are not silently reoriented to force
agreement. The minus sign in the C_-1 UV **y2 interaction** is separately
verified from its `-alpha` scalar epsilon-transport factor (alpha=+1).
B_-1 has alpha=-1 and a positive explicit y2 factor.

The separate source 23/24 *full Majorana exchange* coefficients (external
leg antisymmetrization, mixing spurions and scalar field order) are **not**
independently contracted by this package. Neither the physical Matchete
`GammaCC` to master Weyl phase nor absolute 16-model physical normalization
is certified by this diagnostic. If a predicted native identity fails,
there is a nonzero exit code and a report must **not** be labelled verified.

## Run from the project root (PowerShell)

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$dir = 'output\general_observable\lo_fermion_boundary'
$zip = "$HOME\Downloads\T3_BC_Native_CG_Products.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

python -X utf8 "$dir\test_bc_native_cg_products.py"
& wolframscript -file "$dir\native_bc_yukawa_cg_probe.wl" *> "$dir\normalization_probes\bc_native_yukawa_cg_run.log"
Get-Content "$dir\normalization_probes\bc_native_yukawa_cg_run.log" -Tail 15

$report = Get-ChildItem "$dir\normalization_probes\bc_native_yukawa_cg_*.json" |
    Sort-Object LastWriteTime -Descending |
    Select-Object -First 1 -ExpandProperty FullName
python -X utf8 "$dir\check_bc_native_cg_products.py" --report $report --check-only
python -X utf8 "$dir\check_bc_native_cg_products.py" --report $report
```

The last command is idempotent: an identical report is shown as `UNCHANGED`.
No `T3_EXPERIMENTAL_*` flags are required by the self-contained native
invariant extractor. Matchete v0.6.0 and a functioning Wolfram front end are
required for the Wolfram stage; NumPy is required for Python.

**Next:** if this succeeds, evaluate *full* C_-1 sources 23/24 and native
Majorana Feynman-rule contractions. The 16-model formal release remains
unaltered.
