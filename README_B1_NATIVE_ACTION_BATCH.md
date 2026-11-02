# B1 next big batch: source-level factors, chirality, and installed abstract C

**Intended repository:** `C:\thesistemp\T3neutrinomassmodelpipeline`

**Source checked:** GitHub `main` commit `5e88158147bfce94e34da534fa0be3c905561e1c`, one commit after the user's previous local `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4`. That commit adds workers' evidence but does not change the seven existing source files guarded by the B1 audit. The user's Windows local source guard must pass again; GitHub was not modified. The current user's 16-model numerical residuals are supplied in the conversation (max mass 2.182344e-16, Y 1.664054e-16, Z zero), not rerun in this container.

## Scope and measured/pending levels

The independent `b1_native_action_contract.py` parser audits **fully emitted Matchete `InputForm` terms** in all three previously executed native models `A_p0`, `B_m1`, `C_m1`, including the exact rational numerical prefactor, the `CG`, `Coupling`, scalar orientation, Yukawa Hermitian conjugate, left/heavy field ordering, `Proj[1]`/`Proj[-1]`, and `GammaCC` placement. No user-authored Mathematica is evaluated by Python. A wrong value of `1/2`, doubling of a Yukawa, an HC scalar flip, a missing charge conjugation, a wrong projector, an extra operator, an absent native term, or a false physical certificate blocks execution.

These *native* operators have (from inspected, committed **executed** native JSON files, not from the old component builder):

| Model | Mass monomial | Coefficient | y1 original | y2 original |
|---|---|---|---|---|
| A_+0 Dirac | `Bar[F] F` | `-1` | `Bar[l] C Proj[1] Transp[Bar[F]]` | `Bar[l] Proj[1] F` |
| B_-1 Majorana | `Transp[F] C F` | `-1/2` | `Bar[l] Proj[1] F` | `Bar[l] Proj[1] F` |
| C_-1 Majorana | `Transp[F] C F` | `-1/2` | `Bar[l] Proj[1] F` | `Bar[l] Proj[1] F` |

Each original Yukawa and its explicit HC carries coefficient `+1`. The parser verifies the HC has the correct scalar/bar, CG/bar, coupling/bar, projector, and fermion order. This establishes the symbolic native action coefficients, **not** the final physically normalized component coefficients without the concrete C/Weyl and SU(2) map.

The independent 4→2 Grassmann derivation in existing `b1_matchete_weyl.py` remains:

\[
C=+i\gamma^2\gamma^0,\quad
\overline L F^c=(\ell\chi)^\dagger,\quad
\overline L F=(\ell\xi)^\dagger,
\quad \overline F F=(\chi\xi)+\mathrm{h.c.}
\]

For the numerical convention \(\mathcal L_Y=-\nu\psi(Y\phi+Z\phi^*)+\mathrm{h.c.}\), the **chosen-basis** candidate is `Z_y1=-conj(y1*C_y1)` and `Y_y2=-conj(y2*C_y2)` (with the pre-existing Majorana weight-basis metric). The previous 16/16 frozen candidate match uses one simultaneous external neutrino rephasing; it does **not** prove the native `GammaCC` sign relative to the chosen concrete matrix or B2's independent RGBeta epsilon/del axes.

The new Wolfram probe checks the actual installed Matchete's **abstract** GammaCC identities and symbol-value metadata in an independently restricted fresh kernel. It intentionally does not promote those identities to a measured numeric matrix. Its execution is optional to the Python source-action audit, but encouraged.

## Install and execute on Windows

Download `B1_Native_Action_Coefficients_Batch.zip` into Downloads. From PowerShell:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Get-ChildItem "$env:USERPROFILE\Downloads" -Filter 'B1_Native_Action_Coefficients_Batch*.zip' -File |
    Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $zip) { throw 'B1 ZIP absent from Downloads' }
Expand-Archive -LiteralPath $zip.FullName -DestinationPath (Get-Location).Path -ErrorAction Stop

.\.venv\Scripts\python.exe -m unittest discover -s tests -p 'test_physical_b1_native_action_contract.py' -v

$native = 'output\general_observable\physical_normalisation\b1_02\native'
$oldroot = 'output\general_observable\physical_normalisation'
$prior = Get-ChildItem $oldroot -Directory -Filter 'b1_retry_*' |
    Where-Object { Test-Path (Join-Path $_.FullName 'all16_component_check.json') } |
    Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $prior) { throw 'Previous 16/16 b1_retry audit missing' }

# Fresh-kernel abstract CC measurement, independent of the three already executed native models
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$cc = "output\general_observable\physical_normalisation\b1_cc_probe_$stamp\A_p0.json"
wolframscript -file scripts\RunPhysicalB1CCConventionProbe.wl $cc
if ($LASTEXITCODE -ne 0) { throw 'Native abstract CC probe failed; inspect the CC JSON' }

$out = "output\general_observable\physical_normalisation\b1_action_contract_$stamp"
.\.venv\Scripts\python.exe scripts\RunPhysicalB1NativeCoefficientAudit.py `
    --native-dir $native `
    --conditional-audit (Join-Path $prior.FullName 'all16_component_check.json') `
    --cc-abstract-json $cc `
    --output-dir $out
if ($LASTEXITCODE -ne 0) { throw 'B1 native coefficient audit blocked' }

$evidence = Join-Path $out 'native_action_coefficients.json'
$r = Get-Content $evidence -Raw | ConvertFrom-Json
$r | Select-Object status,native_actions_structurally_verified_count,physical_source_certified
$r.previous_conditional_16 | Select-Object conditional_models,max_mass_relative_residual,max_Y_relative_residual,max_Z_relative_residual
```

If the optional `wolframscript` probe fails, run the Python action audit without the `--cc-abstract-json $cc` argument; the output explicitly states `NOT_EXECUTED`. Do not re-use the failed probe as an approved witness.

### Expected (conditional) status

`ThreeNativeActionsCoefficientAndChiralityVerified_AbsoluteCCAndB2Pending`; `native_actions_structurally_verified_count = 3`; `conditional_models = 16`; `physical_source_certified = false`. **These are expectations, not execution claims**. The fresh native probe is **untested in this delivery environment**. On a failure, leave production files untouched and provide the complete terminal output and generated JSON file(s).

## File manifest (all new paths)

- `Numerical/lo_certification/physical_normalisation/b1_native_action_contract.py`
- `scripts/RunPhysicalB1NativeCoefficientAudit.py`
- `scripts/RunPhysicalB1CCConventionProbe.wl`
- `tests/test_physical_b1_native_action_contract.py`
- `output/general_observable/physical_normalisation/b1_native_action_batch_01/HANDOFF.md`
- `output/general_observable/physical_normalisation/b1_native_action_batch_01/evidence.json`
- `README_B1_NATIVE_ACTION_BATCH.md`

Local validation: **44/44** B1 Python tests passed in the assembled test workspace (20 earlier + 24 new). These are synthetic/code tests, not a fresh native Mathematica run.

The Python runner creates separate unique runtime `b1_action_contract_<timestamp>` folders and refuses to overwrite any previous reports. The new static `b1_native_action_batch_01` is explicitly an *unexecuted delivery handoff*.

## Remaining certification steps

A physical certificate requires transport of the source's abstract C convention into one explicit consistent 4×4 Weyl representation (up to allowed simultaneous field rephasings) and B2's **actual** Matchete/RGBeta relative SU(2) component factors and phases. The three new action contracts cannot substitute for the missing B2 measurement or an independent physical one-loop amplitude. The user's fixed-scale UV pipeline and frozen benchmarks remain unchanged. No EFT work: if undertaken later, the required sequence is UV → integrate out all heavy fermions → EFT with both scalars → integrate out both scalars → SMEFT.
