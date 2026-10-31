# Workstream A: ordinary T3 UV-operator completeness — static source audit

**Status:** source-inspected and independently hypercharge-enumerated; proposed
Wolfram replacements are **not executed in Matchete** in this environment.
No full physical completeness certificate is issued.

## Source provenance

- GitHub `main` at `305035e9df8b90e34fc74c7a3f57453ea6c6aae3` (read only).
- `Lagrangian/model/T3ModelCatalog.wl` blob `9a6c47e10fc83d9480f08353682151180c9fbb4a`.
- `Lagrangian/model/T3Fields.wl` blob `cde2b4a85445de9f4336cee509de2081a1e74437`.
- `Lagrangian/model/LagrangianBuilder.wl` blob `048534bee5e1c360139761895abf8d5d44f2dd3a`.
- `Lagrangian/interactions/T3Topology.wl` blob `f05883c3eb07ce52c821c66cd898b49e37949053`.
- `docs/general_scalar_trial_operator_manifest.json` blob `d48ff0ec1932751aef44ce3498d209c5de633273`.

## Field content and symmetries

Each row adds `S1` and `S2` complex scalars and one three-flavour fermion multiplet
`F` to the SM. They are colour singlets and **assumed odd** under exact
`Z2`; all SM fields are **even**. `T3Fields.wl` does *not* register `Z2`
as a Matchete charge. Gauge symmetry is `SU(3)c × SU(2)L × U(1)Y`,
with `Q = T3 + Y`, `Y(L) = -1/2`, `Y(H) = +1/2`.

`Y(S1)=alpha/2`, `Y(S2)=(alpha+2)/2`, `Y(F)=(alpha+1)/2`.

All dimension-three scalar monomials are absent: `Z2` would require
exactly two odd scalars and one SM Higgs, but the two odd SU(2) reps
always have the same integer/half-integer isospin parity, so their
product contains only integer-isospin irreps and cannot combine with H
into a singlet. An odd scalar singlet tadpole is also forbidden by Z2.

## Sixteen-model operator inventory

`Ytot` is the number of allowed complex Yukawa flavour matrices (three
by three for the stated multiplicity), **not** the number of real
physical parameters. `QR/QC` are BSM scalar quartic real/complex
couplings, including the baseline T3 topology quartic and the static
extra quartic manifest; SM-only Higgs quartic is excluded. `B` is the
number of additional complex scalar bilinear couplings besides the
canonically diagonal `S1†S1`, `S2†S2` mass terms.

| Model | (d1,d2,dF) | (Y1,Y2,YF) | F mass | Ytot | Missing Yukawa | QR | QC | B |
|---|---|---|---|---:|---|---:|---:|---:|
| A_-4 | (1,3,2) | (-2,-1,-3/2) | Dirac | 2 | — | 7 | 1 | 0 |
| A_-2 | (1,3,2) | (-1,0,-1/2) | Dirac | 3 | Lbar F S2 | 7 | 6 | 1 |
| A_0 | (1,3,2) | (0,1,1/2) | Dirac | 3 | Lbar F^c S1^dag | 7 | 6 | 1 |
| B_-3 | (2,2,1) | (-3/2,-1/2,-1) | Dirac | 2 | — | 8 | 2 | 0 |
| B_-1 | (2,2,1) | (-1/2,1/2,0) | Majorana | 2 | — | 8 | 8 | 1 |
| B_1 | (2,2,1) | (1/2,3/2,1) | Dirac | 2 | — | 8 | 2 | 0 |
| C_-3 | (2,2,3) | (-3/2,-1/2,-1) | Dirac | 2 | — | 8 | 2 | 0 |
| C_-1 | (2,2,3) | (-1/2,1/2,0) | Majorana | 2 | — | 8 | 8 | 1 |
| C_1 | (2,2,3) | (1/2,3/2,1) | Dirac | 2 | — | 8 | 2 | 0 |
| D_-2 | (3,1,2) | (-1,0,-1/2) | Dirac | 3 | Lbar F S2 | 7 | 6 | 1 |
| D_0 | (3,1,2) | (0,1,1/2) | Dirac | 3 | Lbar F^c S1^dag | 7 | 6 | 1 |
| D_2 | (3,1,2) | (1,2,3/2) | Dirac | 2 | — | 7 | 1 | 0 |
| E_-4 | (3,3,2) | (-2,-1,-3/2) | Dirac | 2 | — | 11 | 1 | 0 |
| E_-2 | (3,3,2) | (-1,0,-1/2) | Dirac | 3 | Lbar F S2 | 11 | 7 | 1 |
| E_0 | (3,3,2) | (0,1,1/2) | Dirac | 3 | Lbar F^c S1^dag | 11 | 7 | 1 |
| E_2 | (3,3,2) | (1,2,3/2) | Dirac | 2 | — | 11 | 1 | 0 |

**Totals:** six missing independent complex Yukawa matrices, 50 extra
complex scalar quartics, eight extra complex scalar quadratic couplings.

### Yukawa operator derivation

Lorentz-scalar chirality choices are represented by `Lbar F S`,
`Lbar F^c S` and their scalar conjugates. For each candidate test
`Y(Lbar) + Y(F or F^c) + Y(S or S^dag) = 0`, `Z2 = +1`, and
`2 ⊗ dF` contains `dS`. The last condition has singlet multiplicity
one for each allowed scalar since SU(2) tensor products are multiplicity
free. Identify `F = F^c` only for neutral Majorana B(-1)/C(-1).

The existing topology has `Lbar F^c S1` and `Lbar F S2^dag`, each + HC.
For alpha=0, `S1` is complex and Y=0: independently add
`Lbar F^c S1^dag` + HC. For alpha=-2, `S2` is complex and Y=0:
independently add `Lbar F S2` + HC. None are HC duplicates because
the same Lbar/F orientation occurs in both partner operators.

Each extra flavour matrix has nine complex entries before flavour-basis
quotienting. The added operators have different fermion species and need
no 1/2 identical-fermion symmetry factor. The replacement preserves
`PlusHc` and the existing CG-generation conventions; CG normalization
and `CheckLagrangian` still require live testing.

### Fermion mass structures

- B(-1) and C(-1): real integer-isospin and Y=0 support
  `-1/2 F_r^T C M_rs F_s + HC` with symmetric `M_rs`; Takagi
  diagonalization yields three real nonnegative masses.
- Remaining fourteen models: the current four-component non-self-conjugate
  F field is vector-like, supporting `-Fbar_L,r M_rs F_R,s + HC`;
  biunitary rotation yields three nonnegative singular values.
- `T3Fields.wl` currently defines NFlavor multiplicity 3 but mass
  `Mass -> {Heavy, MF}` (one unindexed symbol); this does not
  visibly represent three independent eigenvalues. Confirm Matchete
  flavour mass semantics before claiming an unrestricted mass basis.
- `Z2` forbids odd-fermion/SM-fermion mass mixing.

### Scalar quadratic, quartic and trilinear basis

The current `GeneralScalarOperators.wl` and
`GeneralScalarVerifiedQuartics.wl` implement a symmetry-projected scalar
quartic basis, but the generated interactions are **opt-in**.
Additional scalar bilinears for the eight models above are holomorphic
neutral `S_i S_i` or, for B(-1)/C(-1), `S1 S2` plus HC.
`B(-1)/C(-1)` has mixing between charged scalar components; the
non-diagonal heavy scalar mass basis fails the documented Matchete
`HeavyMassBasis` check and requires an explicit diagonalization.
An operator allowed by gauge/Z2 is not removed when the matching basis
is inconvenient. No scalar trilinear interactions are allowed.

### Implementation status and validation separation

1. **Inspected:** builder always generates only two topology Yukawas in its
   original version; conjugate-neutral-scalar Yukawas were absent.
2. **Modified, opt-in:** new `T3_EXPERIMENTAL_GENERAL_FERMIONS=1`
   creates the additional CG tensor and independent Yukawa coupling,
   adds one extra candidate in each of the six models, and rejects a
   failed required extra Yukawa instead of silently dropping it.
3. **Builder interface:** added expected/included/individually-validated
   extra Yukawa counts and an `FermionYukawaInventoryAssembled` flag.
   `FullSymmetryCompletenessCertified` and
   `FermionFlavorMassGeneralityCertified` intentionally remain `False`.
4. **Tested locally here:** independent Python enumeration passes all
   16 static model checks, including 6/50/8 totals.
5. **Not tested here:** no Wolfram kernel/Matchete installation in this
   environment; no real `CheckLagrangian` or compiled CG checks.
6. **Other workstreams:** RGE, matching and numerical code are untouched.
   The new Yukawa coefficients must not be silently fed into old RGEs.

### Canonical output for downstream workstreams

`BuildT3Lagrangian[model]` retains the original `LUV`, `LInt`, `LFree`,
`AllowedInteractions` and `Status` keys and now exposes explicit
`AdditionalYukawasExpected`, `AdditionalYukawasIncluded`,
`AdditionalYukawasIndividuallyValidated`,
`FermionYukawaInventoryAssembled`, plus scalar inventory counts.
The standalone `RunGeneralFermionAudit.wl` writes both an
`*.interaction.inputform.txt` and an `*.uv.inputform.txt` file alongside
each JSON report, without running Match or RGBeta.
Downstream code must require the requested sector inventory and
full gauge checks **separately**. In a full general scalar run,
`B(-1)/C(-1)` remain blocked for mixed heavy scalar masses.

### How to validate on Windows (do not push)

From repository root after confirming that local files match the
inspected GitHub revision and installing full replacement files:

```powershell
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS = "1"
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS = "0"
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS = "0"
python Lagrangian/RunGeneralFermionStaticAudit.py --outdir docs
foreach ($case in @("A -2", "A 0", "B -1", "C -1", "D -2", "D 0", "E -2", "E 0")) {
  $fields = $case.Split(" ")
  $cls = $fields[0]; $alpha = $fields[1]
  $tag = "$($cls)_$($alpha.Replace("-","m"))"
  wolframscript -file Lagrangian/RunGeneralFermionAudit.wl "output/uv_fermion_audit/$tag.json" $cls $alpha
  if ($LASTEXITCODE -ne 0) { throw "Matchete fermion test failed: $case" }
}
```

After these smoke tests, use both `T3_EXPERIMENTAL_GENERAL_FERMIONS=1`
and `T3_EXPERIMENTAL_GENERAL_QUADRATICS=1` for a separate fully assembled
UV operator test. Do not interpret mixed scalar heavy-mass failures as
absence of a gauge-invariant operator. Full flavour masses and vacuum
conditions remain unresolved.
