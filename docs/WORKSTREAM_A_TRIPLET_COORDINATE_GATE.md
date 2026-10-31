# Workstream A: native triplet-coordinate gate

This opt-in stage addresses the unresolved physical component-basis bridge between
Matchete and RGBeta. It is **not** a Matchete--RGBeta physical normalization certificate.
It uses no source mutations or imported, guessed RGBeta generator matrices.

## Inputs (existing local files)

- `Lagrangian/RunT3UVNativeBridge.py`: the **unchanged previous-stage** reader.
- `Reports/output/workstream_a/native_cg_bridge/native_cg_bridge.json` and
  all 16 `native/*.json` reports created by the preceding **48/48** native-CG audit.
- `RGE/running/rgbeta/T3RGBetaModel.wl`, evaluated using locally installed RGBeta.

This ZIP adds only:

1. `Lagrangian/RunT3RGBetaTripletCoordinateProbe.wl`
2. `Lagrangian/RunT3UVTripletBridgeGate.py`
3. `docs/WORKSTREAM_A_TRIPLET_COORDINATE_GATE.md`

The Wolfram driver evaluates **installed RGBeta** `delS2[SU2L,k,i,j]` and
`tGen[SU2L[S2],a,i,j]` for literal integer indices, after building the
A(-2) model. This is a **probe**: symbolic results are an expected possible
outcome and do *not* imply an invalid representation. It records metadata and,
only if returned values are exact numeric expressions, component arrays.

The Python verifier takes the actual registered `T3MixCG` of A and D
(dS1,dS2)=(1,3) and (3,1) for all three hypercharges per family. Write the
native symmetric topology tensor as \(Q_{(ij),a}=C_{ij a}\), \(Q\in\mathbb C^{4\times3}\),
using its saved **live Matchete** components. Test

\[
Q^\dagger Q=g\mathbf1_3,\quad g>0,\qquad
E=\begin{pmatrix}1&0&0\\0&1/\sqrt2&0\\0&1/\sqrt2&0\\0&0&1\end{pmatrix},
\quad U=E^\dagger Q/\sqrt g.
\]

The code verifies \(U^\dagger U=1\) and \(EU=Q/\sqrt g\). For
\(t_a=\sigma_a/2\), define the covariant Higgs-pair generator
\(D_a=t_a^T\otimes1+1\otimes t_a^T\), and the **derived** Matchete triplet
representation

\[
R_a=-\frac{1}{g}Q^\dagger D_aQ,\qquad
D_aQ+QR_a=0.
\]

It checks exact Hermiticity, \([R_1,R_2]=iR_3\) and cyclic permutations,
\(\sum_a R_a^2=2\mathbf1\), and
\(\operatorname{tr}(R_aR_b)=2\delta_{ab}\). For the standard symmetric-pair
singlet metric

\[
J_{\rm ref}=\begin{pmatrix}0&0&1\\0&-1&0\\1&0&0\end{pmatrix},
\qquad J_M=U^TJ_{\rm ref}U,
\]

it verifies \(J_M^T=J_M\), \(J_M^\dagger J_M=1\) and
\(R_a^TJ_M+J_MR_a=0\). This is a **derived Matchete-side** intertwiner;
Matchete's actual independent native triplet generator matrices have not
been exported by this stage.

For E(-4,-2,0,+2), the stage additionally probes whether the native
\(HH S_1S_2^\dagger\) tensor is covariant using the separately derived
A-side \(S_2^\dagger\) and D-side \(S_1\) generators. A mismatch is
preserved in the report (no phase/basis adjustment is silently inserted).

If RGBeta exposes *both* exact native `delS2` and `tGen` integer-index
arrays, the Python gate tests their Lie algebra, the `delS2` covariance,
and a candidate Matchete/RGBeta unitary matrix. This is still a *candidate*
until the full field, conjugation and signed action/coupling normalization
are independently verified. **No candidate is promoted to a physical map**.

## Execute in PowerShell (from the project root)

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = "$env:USERPROFILE\Downloads\WorkstreamA_TripletCoordinateGate.zip"
Set-Location $repo

foreach ($p in @(
    'Lagrangian\RunT3RGBetaTripletCoordinateProbe.wl',
    'Lagrangian\RunT3UVTripletBridgeGate.py',
    'docs\WORKSTREAM_A_TRIPLET_COORDINATE_GATE.md'
)) {
    if (Test-Path $p) { throw "Target exists; do not overwrite: $p" }
}

Expand-Archive -LiteralPath $zip -DestinationPath $repo
python Lagrangian\RunT3UVTripletBridgeGate.py --self-test
if ($LASTEXITCODE -ne 0) { throw 'Triplet self-test failed' }
python Lagrangian\RunT3UVTripletBridgeGate.py --all
if ($LASTEXITCODE -ne 0) { throw 'Triplet native-coordinate gate failed; inspect the log' }
```

The runner performs one RGBeta Wolfram process and reads already saved
Matchete data. It does **not** rerun any of the 16 Matchete builds.

Outputs (new, never auto-overwritten):

- `Reports/output/workstream_a/triplet_coordinate_gate/rgbeta_native_probe.json`
- `Reports/output/workstream_a/triplet_coordinate_gate/rgbeta_native_probe.log`
- `Reports/output/workstream_a/triplet_coordinate_gate/triplet_coordinate_gate.json`
- `Reports/output/workstream_a/triplet_coordinate_gate/triplet_coordinate_gate.md`

To revalidate existing outputs without another Wolfram call, run
`python Lagrangian\RunT3UVTripletBridgeGate.py --all --use-existing`.

If both RGBeta arrays remain symbolic, the **diagnostic** (not physical
normalization) result is
`LiveMatcheteTripletRepresentationVerifiedRGBetaCoordinatesUnavailable`.
A successful candidate result would be
`CrossPackageTripletCoordinateCandidateVerifiedNotPhysicalDictionary`.
This stage never reports full UV physical normalization certified.

The existing 46 Stage23 coupling-map blocks, the 50 additional quartics,
and all previous 16-model Lagrangian completeness reports remain unchanged.
The physical signed Matchete/RGBeta coupling dictionary, independent native
RGBeta coordinate basis, explicit \(Z_2\) registration, and EFT matching are
not certified by the diagnostic.
