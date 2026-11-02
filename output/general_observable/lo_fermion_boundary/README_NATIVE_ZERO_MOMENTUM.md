# A_-4 independent signed Dirac exchange: native Matchete source

This diagnostic is intentionally confined to `output/general_observable/lo_fermion_boundary/`.
It neither overwrites the released 16-model boundary nor modifies the existing
native Feynman-rule report. The only output is a new JSON in `normalization_probes/`.

## What it computes

The checked A_-4 native UV Lagrangian has the heavy mass and selected Yukawa orientations

\[
\mathcal L_{\rm UV}\supset -M_F\bar F F
+\bar y_{1i}\,L_i^T C P_L F\,S_1^*\,\mathcal G_1
+\bar y_{2j}\,\bar F P_L L_j\,S_2\,\mathcal G_2+\mathrm{h.c.}
\]

Here the independent SU(2) CG tensors \(\mathcal G_{1,2}\) are **not** silently
set equal to the matched `T3MixCG`; that identity must be verified separately.

The four-component Dirac action is contracted by exact Grassmann algebra in the
specified convention \(C=\mathrm{diag}(\epsilon,-\epsilon)\),
\(\epsilon_{12}=+1\), and \(P_L=\mathrm{diag}(1,1,0,0)\).
The zero-momentum propagator is \(S_F(0)=-i/M_F\). The calculation includes
both terms in the second-order interaction expansion and both external-lepton
assignments. In this convention the resulting local spinor kernel is

\[
+\frac{\bar y_{1i}\bar y_{2j}}{M_F}L_i^T C P_L L_j
+\frac{\bar y_{1j}\bar y_{2i}}{M_F}L_j^T C P_L L_i,
\]

with respective SU(2) factors suppressed. The native Matchete EFT1 Feynman
rule is checked **independently** from the previously generated JSON for its
`+i` contact sign and two lepton orderings. It does not try to multiply
scalar vertex prefactors without the Grassmann and fermion-flow signs.

**Limits.** This is a signed zero-momentum tree exchange computation, not an
all-momentum scattering calculation. The full SU(2) CG tensor identity is
not certified by this Python script. `C -> -C` changes the UV result and the
native contact together, so this comparison cannot determine the *absolute*
`GammaCC` to master-Weyl phase. Majorana B_-1 and C_-1 are not part of this check.

## Run in PowerShell

From the project root:

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
Set-Location $repo
$dir = 'output\general_observable\lo_fermion_boundary'
python -X utf8 "$dir\test_native_zero_momentum_a_m4.py"
python -X utf8 "$dir\native_zero_momentum_a_m4.py" --check-only
python -X utf8 "$dir\native_zero_momentum_a_m4.py"
```

The default source is the exact existing on-disk report
`normalization_probes/native_feynman_a_m4_45f6bbfe96e54fc79eb25362cb3e9d92.json`.
Pass `--report PATH` if the report was relocated. Source InputForm digests,
fermion-flow structures and prior non-certification flags are validated.
For a given source file the output name is deterministic and an existing
output is **never overwritten**. Use `--check-only` to rerun validation.

### Optional next CG evidence

For the native Yukawa versus `T3MixCG` identity, the included Wolfram script
exports all three exact installed-Matchete CG tensors with flattened real/imag
components. It does not declare their contraction complete. With the four
experimental flags already used by the earlier native FR script set to `1`:

```powershell
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS  = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS = '1'
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS = '1'
$env:T3_EXPERIMENTAL_FLAVOR_MASSES = '1'
& wolframscript -file "$dir\native_cg_a_m4_probe.wl" *> "$dir\normalization_probes\native_cg_a_m4_run.log"
Get-Content "$dir\normalization_probes\native_cg_a_m4_run.log" -Tail 8
```

The CG report is written under the same `normalization_probes/` folder.
No GitHub push or repository modification is needed for these commands.
