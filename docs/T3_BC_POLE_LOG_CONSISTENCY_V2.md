# T3 mixed B(-1)/C(-1): exact formal pole/log coefficient audit v2

## Scope and verified statement

Input is the **v1 finite-term audit output** of the native Matchete Y-only and X-only LLSS scalar one-loop amplitudes. This stage checks, algebraically and term-by-term in a free commutative ring of *opaque Matchete factors*, that

\[
\mathcal A^{\rm raw}_{1\ell} = \frac{\mathcal P}{\epsilon}
  + \mathcal L\ln\frac{\bar\mu^2}{M_{Y,X}^2} + \mathcal K,
\qquad \mathcal P=\mathcal L.
\]

The identity was independently replayed from the B/C uploaded Matchete evidence, including hashes of the **original raw source files** and of all three v1 split products. The mathematical equality holds in the opaque-factor algebra; it is **not** a proof of a physical RGE decoupling identity or a computed MS-bar matching counterterm. The derivative of the pole-free expression with respect to `ln(mu)` at *fixed masses/couplings* is formally `2 P`, but the true EFT RGEs, anomalous dimensions and normalizations must be checked separately.

A full `NCM[...]` (including fermions/chirality matrices) is deliberately one indivisible atom. Likewise `CG[...]`, `Field[...]`, `Coupling[...]`, `Bar[...]` and all other Mathematica calls are opaque. The arithmetic parser supports exact rational multiplication, integer powers, products/divisions and overall addition only. It does **not** use string substitutions that assume SU(2)/Fierz identities. Unrecognized structures abort before any model outputs are written. The parser does not run Mathematica or execute source expressions.

## Audited results

| Model | Sector | Pole terms | Log terms | Constant terms | Identity |
|---|---|---:|---:|---:|---|
| B(-1) | Y/LLHH | 50 | 50 | 32 | P = L |
| B(-1) | Y/LLHX | 0 | 0 | 0 | Trivial |
| B(-1) | Y/LLXX | 240 | 240 | 186 | P = L |
| B(-1) | X/LLHH | 18 | 18 | 0 | P = L |
| C(-1) | Y/LLHH | 50 | 50 | 32 | P = L |
| C(-1) | Y/LLHX | 0 | 0 | 0 | Trivial |
| C(-1) | Y/LLXX | 198 | 198 | 186 | P = L |
| C(-1) | X/LLHH | 18 | 18 | 0 | P = L |

### Formal finite constants by source quartic

Identical raw additive-term counts in B and C; this is a partition, **not** a statement of independent physical operators. These are *raw Matchete scale-equal-mass* terms, not renormalized coefficients:

| Y threshold operator sector | Quartic coupling | Additive terms |
|---|---|---:|
| LLHH | `lambdaT3` | 16 |
| LLHH | `lambdaT3ExtraQ020020I1` | 8 |
| LLHH | `lambdaT3ExtraQ020200I1` | 8 |
| LLXX | `lambda12Inv2` | 30 |
| LLXX | `lambdaS1Inv1` | 16 |
| LLXX | `lambdaS2Inv1` | 16 |
| LLXX | `lambdaT3ExtraQ000112I1` | 46 |
| LLXX | `lambdaT3ExtraQ000202I1` | 32 |
| LLXX | `lambdaT3ExtraQ001201I1` | 46 |

Individual exact, unhashed-source-dependent InputForm sums are written into `{model}/Y_{LLHH|LLXX}_{quartic}_finite_at_mass.inputform.txt`. Existing Matchete field/CG conventions are left untouched. The X/LLHH local constants at `mubar=MSX` are formally zero in the saved local amplitude, which **does not prove** full finite physical X matching zero. No `MSX`/`MSY` dependence appears in the explicit constant terms of these archived Y sources, so this is a strict hard-region expansion and not a finite-mass-resummed bubble result.

## Installation / execution

Put the complete ZIP into the repository root, retaining all relative project paths.

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_bc_pole_log_consistency_v2.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_t3_bc_pole_log_consistency_v2.py' -v
python -u -X utf8 -m scripts.AuditT3BCPoleLogConsistency --all
```

Requires the successful earlier `python -u -X utf8 -m scripts.ExtractT3BCFiniteScalarTerms --all` run and original saved Y/X raw Matchete inputs. If only **v1 split files** are archived, an explicit `--skip-original-source-validation` option can audit their reference hashes; it never promotes physical certification.

Outputs go under `output/audits/t3_bc_pole_log_consistency_v2/`:

- `summary.csv`, `summary.json`: eight sector results, explicit negative certification flags.
- `{B_m1|C_m1}/{Y_LLHH|Y_LLHX|Y_LLXX|X_LLHH}_pole_log_audit.json`: exact coefficient checks.
- `{B_m1|C_m1}/Y_{LLHH|LLXX}_{quartic}_finite_at_mass.inputform.txt`: 9 grouped finite-constant source expressions per model.

No existing files are overwritten outside this dedicated output directory.

## Remaining physics, not claimed

1. Verify dimensional-regularization pole cancellation with Matchete/EFT counterterms and light-field loop subtraction, including precise `hbar`, operator, and scheme conventions.
2. Prove equivalence between the original F-first EOMEFT and the independent canonical EFT1 construction (or derive their EFT field redefinition).
3. Translate finite `LLHH` and `LLXX` into the numerical real-doublet Wilson basis, retaining flavor, SU(2), conjugation, chirality and Majorana normalization.
4. Add inherited finite Y matching, substitute the **RG-evolved** masses/couplings at `mu_Y` and `mu_X`, and compare the threshold beta discontinuity.

**Do not** add these formal constants to a certified neutrino mass. No numerical masses are altered by v2.
