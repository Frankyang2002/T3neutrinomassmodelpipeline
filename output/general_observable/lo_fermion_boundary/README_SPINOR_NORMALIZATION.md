# GammaCC -> Weyl spinor normalization: next verification stage

**All files belong in:** `output/general_observable/lo_fermion_boundary/`.

This package adds **new** diagnostic files. It does not replace `build_boundary.py`, `test_boundary.py`, or the certified-as-source-resolved 16-model release. No files outside this worker's output directory are written.

## What is independently checked

In the Weyl representation with `metric = diag(+,-,-,-)` and

```
  gamma^0 = ((0,I),(I,0)), gamma^k = ((0,sigma^k),(-sigma^k,0)),
  gamma^5 = diag(-I,+I), P_L = (1-gamma^5)/2,
  epsilon = i*sigma^2 = ((0,1),(-1,0)),
  C = i*gamma^2*gamma^0 = diag(epsilon,-epsilon), C^-1 = -C,
```

exact four-component matrix multiplication gives

```
  (chi,0)^T C P_L (eta,0)    = chi^T epsilon eta
  (chi,0)^T C^-1 P_L (eta,0) = -chi^T epsilon eta
```

The identities are exact. The master `WeinbergTensorAdapter.py` extracts a coefficient with `(1/4)*d_i*d_j*d_a*d_b` from the commuting-placeholder polynomial; the 1/2/4 distinct-orbit rules are checked exactly. This demonstrates no extra *formal polynomial* factor.

**It does NOT establish whether the Matchete installation defines `GammaCC` as `C`, `C^-1`, `-C`, or uses an additional field/Grassmann ordering convention.** A numerical factor of +1 or -1 is therefore *conditional*, not physically certified. The running/scalar workers must continue to reject the current source-resolved release as a fully normalized physical input.

## PowerShell commands

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$dir = 'output\general_observable\lo_fermion_boundary'

# Unit tests: no write, no Matchete required
python -X utf8 "$dir\test_spinor_normalization.py"

# Matchete evidence: read-only local Mathematica package inspection
& wolframscript -file "$dir\probe_matchete_gamma_cc.wl" *> "$dir\matchete_gamma_cc_evidence.txt"

# Build a new source-stamped, non-certified diagnostic report
python -X utf8 "$dir\probe_spinor_normalization.py" --matchete-evidence "$dir\matchete_gamma_cc_evidence.txt"
```

Run the final Python command **without** `--matchete-evidence` if `wolframscript` or Matchete is not installed; the algebra is still checked but does not resolve the sign:

```powershell
python -X utf8 "$dir\probe_spinor_normalization.py"
```

To preview without writing any file:

```powershell
python -X utf8 "$dir\probe_spinor_normalization.py" --check-only
```

New reports are written only into `output/general_observable/lo_fermion_boundary/normalization_probes/` and cannot overwrite an existing report with the same source/evidence hashes. The Mathematica console evidence is captured only inside `lo_fermion_boundary`.

## How to complete physical certification

Inspect or share `matchete_gamma_cc_evidence.txt`. Specifically, locate a Matchete package **definition** for `GammaCC`, `CC`, `CConj`, and `Proj` (not merely a name, usage message, or our own `RunMatching.wl` label). Verify which charge-conjugation matrix is used and how `Transp[l] GammaCC P_L l` is ordered. Confirm the master Weyl epsilon phase and a signed tree-level heavy-fermion exchange amplitude in the same field conventions. Only after these independent checks should a new *normalized* release be derived; do not edit the existing `source_021dcc1e8fe8ba0f` release.

**Not yet covered:** one-loop hard matching of heavy fermions and any suppressed leading-order operators beyond the saved F-tree LLSS source set. These must be investigated separately before claiming complete leading one-loop physical matching.
