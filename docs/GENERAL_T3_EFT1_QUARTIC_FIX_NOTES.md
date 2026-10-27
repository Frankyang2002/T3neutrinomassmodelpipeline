# General T3 scalar EFT1 — first narrow generalisation correction

## Verified source boundary

Source inspected: GitHub `main` at the time of packaging (read-only).
Files inspected: `RGE/general/generated/general_rgbeta_quartic_registry.json`,
`Lagrangian/interactions/{GeneralScalarOperators,GeneralScalarVerifiedQuartics,ScalarPotential}.wl`,
`RGE/running/rgbeta/{T3RGBetaModel,RunGeneralT3RGBeta}.wl`,
`Numerical/experimental/GeneralT3ScalarRGState.py`,
`scripts/{Stage33FermionThresholdHandoff,Stage34IntermediateEFTRunning,ExportGeneralT3Stage34EFT1OneLoop}.py/.wl`,
and documentation for Stages 22–29 and 34–36.

The public source contains all 16 generated RGBeta extension files, the
38 new coupling declarations, and all 12 B/C reuse records.
The ordinary Stage34 exporter **explicitly omits** the new Stage29 quartics.

## First correction

Install **only** `scripts/RunGeneralT3EFT1Quartics.py`. This separate
experimental entrypoint uses the current `Stage33` EFT1 boundary and its
38 new quartic sidecars, converts them to the existing
`GeneralT3ScalarRGState`, regenerates **EFT1** full-projected RGBeta beta
functions using the present local `T3RGBetaModel.wl`, and numerically evolves
all baseline, reused, and new quartics from 1e5 to 1e3 GeV.

Pre-existing `Stage34` reports and `Stage35–36` neutrino-mass processing remain
unchanged. This is **renormalisable EFT1 only** (not full EFT including
higher-dimensional Wilson-coefficient running).

## Essential limitations

- Native RGBeta/Matchete physical component/phase/potential conventions are
  **not certified**. The initial new couplings are used in the explicitly
  convention-labelled RGBeta basis already present in Stage33.
- Eight allowed scalar quadratic structures remain absent from scalar mass
  matrices. Six are neutral singlet/triplet holomorphic bilinears; two B/C
  structures are inter-doublet mixing terms. In a mass-eigenstate basis some
  off-diagonal mass parameters are removable via field rotation, but this
  changes Yukawa/quartic tensors. This does **not** justify simply dropping
  interactions from a fully general fixed-basis Lagrangian.
- Fermion-threshold finite matching and the resulting higher-dimensional
  Wilson coefficients are **not** recalculated.
- The new quartic flow is **not** inserted into scalar threshold matching,
  `C5`, or neutrino mass observables.
- Do not feed the output into Stage35 or Stage36 and describe it as a fully
  general T3 result.
- Older `output/general_t3_stage20/*_eft1_rgbeta.json` must **not** be reused
  blindly; the runner regenerates expressions after the later B/C covariance fixes.

## Installation, PowerShell (project root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath C:\path\to\T3_Generalisation_EFT1_Fix.zip -DestinationPath . -Force
python -m py_compile .\scripts\RunGeneralT3EFT1Quartics.py
python -m scripts.RunGeneralT3EFT1Quartics --model B_+1
```

For all 16 once the single-model check succeeds:

```powershell
python -m scripts.RunGeneralT3EFT1Quartics --all
```

`wolframscript` must be on `PATH` with RGBeta installed. The Stage33 files
and hashed Stage30/31 source artifacts must exist at their current standard
location. A missing dependency or mismatched RGE-key set causes a failure,
not substitution of zero beta functions.

## Output project locations

- `output/general_t3_eft1_quartics/<MODEL>_eft1_live_full_rgbeta.json`
- `output/general_t3_eft1_quartics/<MODEL>_eft1_full_quartic_running.json`

The `<MODEL>` in filenames uses `B_p1`, `A_p0`, `E_m2`, etc.

`docs/general_t3_operator_gap_inventory.csv` enumerates 38 new, 12 reused,
and eight blocked quadratic model-specific slots. "Additional quartic slots"
are model-specific: same named symbolic contraction in two different models
is counted twice. The 12 reuse records are six in B_-1 and six in C_-1,
not 12 distinct globally named RGBeta couplings.
