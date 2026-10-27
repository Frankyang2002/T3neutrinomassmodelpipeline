# Fermion-first EFT1: raw hard LLHH pole/log audit

This is a **diagnostic**, not a renormalized or physical matched Weinberg
coefficient. It analyses the successful local Matchete exports produced by
`RunGeneralT3FermionFirstEFT1OneLoop` for `B_-1` and `C_-1`. It reads those
files verbatim; it does **not** redo one-loop Matchete matching, evaluate
SU(2) Clebsch–Gordan tensors, alter quartic normalisations or perform EFT1
running.

## Verified directly against the uploaded user exports

| Item | B_-1 | C_-1 |
|---|---:|---:|
| Full one-loop EFT1 correction terms | 9,438 | 9,380 |
| Raw LLHH terms, signed subset of full correction | 204 | 204 |
| Pole terms containing Unicode `ϵ` | 68 | 68 |
| Log terms containing `Log[μbar2/Coupling[MF, {}, 0]^2]` | 68 | 68 |
| Finite, log-free terms | 68 | 68 |
| Independent signed tensor/flavour monomials in each sector | 68 | 68 |
| Holomorphic `Proj[-1]` terms in each sector | 34 | 34 |
| Conjugate `Proj[1]` terms in each sector | 34 | 34 |

The 68 monomials in each sector contain 28 contributions with `lambdaT3`,
20 with `lambdaT3ExtraQ020020I1`, and 20 with
`lambdaT3ExtraQ020200I1`. Their **full tensor/operator coefficients**
are identical in the pole, log and constant pieces after removal of *only*
their explicit regulator or logarithm factor. The computed raw-operator identity is

\[
\Delta\mathcal L_{LLHH}^{(1)} =
\mathcal A_{LLHH}\left(1+\frac1\epsilon+
\log\frac{\bar\mu^2}{M_F^2}\right),
\]

where \(\mathcal A_{LLHH}\) already includes Matchete's `hbar` and the
**uncontracted** operator tensors. This is a termwise identity, not a scalar
number or a physical Wilson coefficient. At the diagnostic matching point
`μbar2 = MF^2`, its **pole-discarded finite piece** is
`A_LLHH`, subject to counterterm and renormalisation consistency.

A careful distinction is essential: the identical pole and log coefficients
show that \(\partial_{\ln\mu}\log(\bar\mu^2/M_F^2)=2\) multiplies the
same tensor residue. They do **not** independently establish the EFT1
`LLSS → C5` beta function, cancellation of \(1/\epsilon\) by the
appropriate counterterm, or the physical Majorana operator normalisation.

## Source and outputs

Source directory:

`output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop/`

The following existing source files must be present **for each model**:

- `<stem>.json` (successful Matchete matching certificate)
- `<stem>_eft1_one_loop_correction.inputform.txt` (full 9,438/9,380-term correction)
- `<stem>_raw_hard_llhh_sector.inputform.txt` (204 extracted terms)
- `<stem>_llss_tree.inputform.txt` (tree EFT1 source sector)

Run from the repository root in PowerShell:

```powershell
python -m pytest -q tests\python\test_general_t3_eft1_hard_pole_audit.py
python -m scripts.RunGeneralT3EFT1HardPoleAudit --all
```

Outputs:

`output/general_t3_mass_eigen_matchete/eft1_hard_pole_audit/`

Each model has a provenance JSON with SHA-256 hashes of the four input
artifacts, the counts and *explicit false* physics-certification fields, plus
InputForm operator expressions for (i) the raw pole/log/finite sectors,
(ii) the full-operator pole residue, log coefficient, and finite piece at
\(\bar\mu=M_F\), and (iii) the separate `Proj[-1]` and `Proj[1]` finite
branches. Each text output is an uncontracted LLHH **operator expression**,
not a flavour-space \(C_5\) matrix.

The Python parser splits only at top-level Mathematica `+`/`-` symbols and
checks exact reconstructibility. It does not use regular expressions to
rewrite CG tensors, masses, flavour indices or nested algebra; factor
removal is restricted to the two exact exported regulator/log tokens. It
asserts that every raw sector term is an identical signed term of the full
Matchete one-loop correction.

## Next genuine physics step

Use the **independent EFT1 Wilson tensor and direct Weinberg beta function**
to check whether the hard pole residual is cancelled, with the correct CG
and flavour projection. Subsequently construct the renormalized hard
threshold coefficient and direct `LLSS → C5` contribution and combine at the
scalar threshold. These are not yet performed. In particular the one-scale
neutrino masses cannot be relabelled as sequential EFT1 predictions.
