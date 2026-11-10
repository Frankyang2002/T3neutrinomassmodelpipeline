# T3 A(−4): source-grounded direct F-hard / joint-scalar-soft Laurent bridge

## Scope

The source files in `thesis_a4_fresh_f_v3/A_-4` show successful *fresh* matching
of all three F flavours at one common threshold, followed by a native one-loop
joint match of scalar multiplets S1 and S2. This patch inspects the six
first-threshold **direct LLHH** source terms (both chiralities × constant,
1/epsilon, log) against the fresh one-loop **InputForm**, and symbolically
verifies the complete A(-4) joint-scalar **soft Weinberg** coefficient.

The comparison is deliberately restricted to the **ordered direct LLHH
channel**. It is not a projection of the full 135-operator F source onto a
verified physical EFT basis, nor a local counterterm calculation. It may not
be labelled `C5` or used as an input to a NuFIT fit.

## Exact algebra checked in source

For a real diagonal heavy mass matrix and **fixed** F-scale couplings, write

```
Wbar[p,q] = sum_r conj(y1[p,r])*conj(y2[q,r])/M[r]
K = 2*lambdaT3/(sqrt(3)*16*pi^2)
```

The first-threshold direct holomorphic ordered coefficient is

```
A_hard[p,q] = K * sum_r conj(y1[p,r])*conj(y2[q,r])/M[r]
                   * (1 + 1/eps + log(mu_bar^2/M[r]^2))
```

The actual scalar-soft source gives

```
A_soft[p,q] = -K * Wbar[p,q] *
    (1 + 1/eps + (a log(mu_bar^2/a)-b log(mu_bar^2/b))/(a-b))
```

with `a = MS1^2`, `b = MS2^2` in the Matchete symbolic expression.
The source-level `1/eps` poles and constants cancel in the same *stripped
ordered operator* convention, leaving the scale-independent formal term

```
A_direct,formal[p,q] = K * sum_r conj(y1[p,r])*conj(y2[q,r])/M[r] *
       (a log(a/M[r]^2)-b log(b/M[r]^2))/(a-b)
```

A stable `a≈b` implementation approaches `log(b/M[r]^2)+1` and avoids
catastrophic subtraction for degenerate scalar masses.

**Warning:** The source `ExtractWeinbergCoefficient` helper removes SU(2)
Clebsches and `GammaCC` operators to isolate an ordered prefactor. Equality
of those *stripped* prefactors is not a proof of physical operator tensor
identity. The raw Matchete `hbar` loop normalization was treated formally as
`1/(16*pi^2)` per the existing hard-source diagnostic; physical sign and
Weyl/SMEFT normalization require a separate source-level derivation. The
RGBeta S-endpoint scalar masses are used only as an explicitly uncertified
numerical coordinate mapping. The F-scale coupling is held fixed for the
formal one-loop pairing; the actual EFT1 Wilson RG is still absent.

## Replacement/new files

All paths relative to `C:\thesistemp\T3neutrinomassmodelpipeline`:

- `Numerical/parallel_t3/complete/a4_direct_hard_soft.py` — new, complete source algebra and numerical direct-channel implementation
- `scripts/RunT3A4DirectHardSoft.py` — new, independent runner consuming saved sources only
- `tests/test_t3_a4_direct_hard_soft.py` — new unit tests
- `docs/T3_A4_DIRECT_HARD_SOFT_SOURCE_BRIDGE_v1.md` — new documentation

No existing files are overwritten; no GitHub changes, Wolfram re-matching,
new Batch1/B2 runs, or approximated surrogate scalar match.

## PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_A4_DirectHardSoft_SourceBridge_v1.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p test_t3_a4_direct_hard_soft.py -v
python -X utf8 -m scripts.RunT3A4DirectHardSoft
Get-Content output\parallel_t3\a4_direct_hard_soft_bridge_v1\A_-4\A4_direct_hard_soft_formal_audit.json
```

The JSON declares `full_F_finite_MSbar_Wilson_boundary_computed=false`,
`complete_eft1_wilson_anomalous_dimension_running_performed=false`, and
`physical_C5_computed=false`. It contains exact file SHA256 inputs, the
unambiguous A-direct Laurent audit and **formal** 3x3 numerical remainder.

## Subsequent required physics

1. Prove the field/CG/GammaCC-to-Weyl operator identity with indices and
   normalization and include **all** F-threshold operators that mix into LLHH.
2. Derive/insert the complete finite UV-hard minus EFT1-soft plus local
   MSbar counterterm matching, preserving the remainder of the full raw
   matched operator ledger, not only the six direct terms.
3. Evolve the entire relevant EFT1 Wilson operator multiplet from F to S
   with the existing running renormalisable couplings.
4. At S, combine inherited hard and jointly matched S1/S2 contributions,
   verify pole cancellation in the physical Weinberg basis for A and D, and
   produce the finite 3x3 symmetric C5 matrix needed by the existing SMEFT
   and Yukawa-fitting machinery.
