# A(-4) source CG basis bridge — v8

This package fixes one **identified numerical-basis mismatch** in the rough
fermion-first pipeline, while leaving the original UV, Stream 1, Stream 2, and
Stream 3 source and output files untouched.

## Verified source and map

Source file (tracked in the project):

`output/general_observable/fermion_first_hermitian_portal_v1/A_m4/data/eft1_after_F_wilson_seed.json`

The Matchete PL tree term is `L L NewScalar1† NewScalar2`; `T3MixCG` is stored as
`SparseArray[Automatic,{2,2,3},0,...]` with six nonzero coefficients. The
matrix `G[i,j,k]` for the **Matchete native triplet slot** k is

- `G[0,0,0]=3^(-1/4)`, `G[0,0,2]=i 3^(-1/4)`;
- `G[0,1,1]=G[1,0,1]=-3^(-1/4)`;
- `G[1,1,0]=-3^(-1/4)`, `G[1,1,2]=i 3^(-1/4)`.

Its unique unitary equivalence to the normalised canonical spherical singlet
coefficient `T[i,j,m]` is

```
U_native_from_canonical = [[-1/sqrt(2), 0,  1/sqrt(2)],
                           [ 0,          1,  0        ],
                           [-i/sqrt(2),  0, -i/sqrt(2)]]
G_native @ U_native_from_canonical = sqrt(6)*3^(-1/4)*T_spherical
```

**This factor is a CG-normalisation comparison, NOT a multiplication of the
Wilson coefficient.** The algorithm applies only the *unitary* U. A global
signed GammaCC-to-left-Weyl factor remains a separate upstream assumption.

With `S2_native = U S2_canonical`, real interleaved coordinates obey

```
R2[2i:2i+2, 2j:2j+2] = [[Re(U[i,j]), -Im(U[i,j])],
                          [Im(U[i,j]),  Re(U[i,j])]]
R_heavy_source_from_canonical = diag(identity_2_for_S1, R2)
```

Every real scalar source index is transformed by this orthogonal 8x8 map:

`C_ijAB = sum_{ab} C_ijab(source) R[a,A] R[b,B]`

Thus the tensor norm and original flavour information are conserved. This is
NOT the v7 Ward-singlet projection, which retained only ~25.5% of the LLSS
norm. The raw B/Stream1 scalar coordinate contract is explicitly checked,
along with the exact original source SHA256 and native `T3MixCG` text.

An independent arithmetic example using the archived Matchete CG exactly
reproduces the measured identity-map SU2 Ward norm `1.5086889589691435`;
applying R reduces it to floating point noise, retaining 100% of the norm.

## Files

- `Numerical/parallel_t3/wilson/matchete_a4_map_v8.py` (source CG/read-only transformation)
- `scripts/RunT3MatcheteMappedA4V8.py` (new isolated runner)
- `tests/test_t3_matchete_a4_map_v8.py` (new tests)
- `Numerical/parallel_t3/wilson/README_SOURCE_CG_A4_V8.md` (this file)

No existing source files are replaced. The original rough Stream1/2/3 outputs
are not edited. The new runner reads the same `A_m4` Stream1 real boundary,
replaces ONLY the in-memory identity scalar map by the source-derived R, and
uses the EXISTING Stream2 EFT1 DOP853, beta, matching and SMEFT RGEs.

## Commands (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

$zip = Join-Path $env:USERPROFILE 'Downloads\t3_rough_A4_source_CG_map_v8.zip'
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_matchete_a4_map_v8.py

# Verify source tensor, no norm loss, Ward identities and MF beta; no ODE.
.\.venv\Scripts\python.exe scripts\RunT3MatcheteMappedA4V8.py --boundary-only --overwrite

$r = Get-Content output\parallel_t3\source_cg_mapped_a4_v8\source_cg_map_summary.json -Raw | ConvertFrom-Json
$r.MF_source_boundary_map | Select-Object LLSS_norm_fraction_retained,source_cg_mapped_ward
$r.MF_mapped_boundary_beta_Ward

# If successful, run existing EFT1 -> one common scalar MS -> SMEFT -> masses.
.\.venv\Scripts\python.exe scripts\RunT3MatcheteMappedA4V8.py --overwrite

$r = Get-Content output\parallel_t3\source_cg_mapped_a4_v8\source_cg_map_summary.json -Raw | ConvertFrom-Json
$r | Select-Object stage,masses_eV_CONDITIONAL,MS_full_Ward
```

Outputs are isolated in
`output/parallel_t3/source_cg_mapped_a4_v8/`, with all original
boundary and transformed endpoint evidence recorded in
`source_cg_map_summary.json`. A `*.blocked` sibling is written if a check
fails; the program exits nonzero. A case that fails the exact identity-source
reconstruction gate is never run through the EFT1 ODE.

## EFT and certification

The order is **UV → integrate out ALL fermions together at MF=10^5 GeV → EFT1
with BOTH scalar multiplets → integrate out BOTH scalars together at MS=10^3
GeV → SMEFT → low-energy neutrino matrix**.

Results retain the status `SourceCGMappedA4NumericalCandidateNotCertified`.
The triplet *coordinate* error has been corrected from the Matchete source,
but GammaCC/Weyl overall normalisation, other UV one-loop hard matching,
the complete native scalar/fermion operator basis, and the candidate EFT1 ADM
have not been independently certified. Do NOT use these numbers as a physics
exclusion or validated neutrino prediction.

## Verification from this delivery

Local isolated tests pass, using repository-compatible minimal basis/test
stubs; the actual Windows repository test and genuine numerical ODE must be
executed with the supplied PowerShell commands. This ZIP does not contain
fake Matchete or B trajectories.
