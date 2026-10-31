# Workstream C — v19: restricted one-loop hard/soft MSbar subtraction

## Why

The v18 three-fermion reference topology passed, but the difference between
raw pole cancellation and an actual high-minus-low amplitude subtraction
remained explicit. This calculation identifies the pole and finite pieces
of the **scalar-quartic/LLSS triangle topology only**, at leading order in
`1/MF`, retaining the full `MX/MY` ratio in the renormalized intermediate-EFT scalar bubble.

## Exact identities

For `A=MF^2`, `B=Ma^2`, `C=Mb^2`, and `mu2=mu^2`:

```
Phi_strict = (B*log(B/A)-C*log(C/A))/(B-C)
B0fin      = 1-(B*log(B/mu2)-C*log(C/mu2))/(B-C)
H(mu)      = 1+log(mu2/A)
S(mu)      = -B0fin
Phi_strict = H(mu)+S(mu)
```

At `B=C`, the continuous result is `Phi_strict=1+log(B/A)`.
`H` and `S` have opposite scale slopes `+1` and `-1`.
The corresponding restricted UV pole residues are `+1` for the hard
region and `-1` for the EFT soft loop. This is a method-of-regions
identity in dimensional regularization, with the scalar EFT bubble
MSbar-subtracted at the common scale.

From source-reconstructed couplings `WXX`, `WXY`, `WYY` and the existing
F/Y/X pole polynomials, the code **requires**, for each of the three
Yukawa structures, the exact identities

```
WXX + WXY + WYY == PF
-(WXX + WXY + WYY) == PY + PX
```

It then reconstructs the B/C three-fermion **frozen-input reference**
`C5-like` matrix from `H` and `S` separately, matching the v18 **finite_XY**
matrix at `mu=MY`. The v18 `strict` matrix uses an extra `MX/MY -> 0`
expansion in the mixed channel, which is NOT part of the leading `1/MF`
expansion itself. The difference quantifies higher-dimensional operators
needed at the Y threshold in a strict sequential truncation. Three other scales are checked analytically to prevent
an accidental chosen-scale cancellation.

This is *not* a comprehensive matching theorem for all operators. It does
not perform the finite coupling/field redefinitions, prove the UV full
amplitude vertex sign, run EFT1/SM couplings or implement the tower of higher-dimensional
`MX/MY` corrections in the sequential Y-threshold EFT (these corrections
are retained only in the full intermediate-EFT bubble here). Consequently **full physical C5 = null**.

## Files

- `scripts/VerifyT3BCHardSoftSubtraction.experimental.py`
- `scripts/RunT3BCHardSoftSubtraction.ps1`
- `docs/WorkstreamC_HardSoftSubtraction_v19_README.md`

## Run from repository root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_HardSoftSubtraction_v19.zip" -DestinationPath . -Force
& .\scripts\RunT3BCHardSoftSubtraction.ps1
```

Output:
`output/audits/t3_bc_hard_soft_subtraction_v19/summary.json`.

Requires existing local v14-v18 summaries and original Matchete source
artifacts, together with the project numerical frozen-reference files.
Does not edit upstream scripts, production matching code or GitHub.

## Remaining physics gate

The scalar triangle's **restricted hard-minus-EFT soft** coefficient is
independently defined after MSbar subtraction, but the complete physical
T3 matching still requires establishing finite parameter/EOM shifts,
other operator contributions and running. No additional frozen-input
phenomenology is needed until this is addressed.
