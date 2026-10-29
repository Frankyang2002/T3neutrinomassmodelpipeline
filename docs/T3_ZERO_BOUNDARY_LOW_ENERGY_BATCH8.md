# Batch 8 — zero-additional-coupling B/C low-energy neutrino continuation

## Purpose

Bridge the successful **Batch 7** zero-extra UV+EFT1 renormalisable RGE benchmark to the **existing Stage 36** model-specific hierarchical one-loop \(C_5\) matching and SM+Weinberg evolution. Deliver concrete B(-1)/C(-1) low-energy mass matrices, Takagi masses and splittings for the zero-extra-coupling boundary. Do not re-fit.

The same preexisting scalar/Yukawa/mass representation and the same frozen fitted UV benchmark are used. All seven extra quartics plus the allowed extra quadratic mixing start at **exactly zero at \(10^7\) GeV**. The solver **replays** Batch 7 from the hashed candidate UV/EFT1 RGBeta inputs and requires agreement with all four saved Batch 7 snapshots, rather than infer a complete RG numerical state from compact scalar snapshots.

1. Native candidate one-loop RGBeta, \(10^7 \to 10^5\) GeV with the complete UV state; cross-check the zero-mixing UV beta vector against the independently saved Stage31 UV export.
2. Tree-project at the fermion threshold; native candidate scalar EFT1 one-loop RGBeta, \(10^5 \to 10^3\) GeV, including the allowed quadratic coefficient.
3. Fail if saved Batch 7 input hashes or endpoint snapshots disagree, if scalar squared masses are nonpositive, or if the extra quadratic mixing is appreciable: the existing Stage36 matching presumes diagonal complex scalars.
4. Reuse the model's **own existing validated hierarchical** fixed-order \(C_5\) artifact; evaluate it with the newly RG-evolved Yukawas, masses and baseline \(\lambda_{T3}\) using `Numerical/running/FinalC5TrajectoryAdapter.py` and `SMWeinbergStage.py`.
5. Reuse SM+Weinberg one-loop evolution \(10^3 \to 100\) GeV and existing neutrino observables code, following `scripts/FinalWeinbergRunning.py`.
6. Compare \(\|C_5^{\rm zero-boundary}-C_5^{\rm Stage36}\|_F/\|C_5^{\rm Stage36}\|_F\) at \(10^3\) and \(100\) GeV, and output both neutrino spectra. **This is a boundary-scenario comparison**: original Stage36 starts with nonzero Stage29 additional quartics, not the zero boundary. It is NOT \(\epsilon_{C_5}\) from inserting a generated extra quartic alone.

## Delivered project paths

- `scripts/RunT3ZeroBoundaryLowEnergyC5.py` — new complete batch runner
- `tests/python/test_t3_zero_boundary_low_energy_batch8.py` — offline tests
- `docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md` — permanent restricted-potential conclusion
- `docs/T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md` — this implementation and interpretation guide

Existing production files are not replaced. The batch **depends on** the previously installed Batch 7 script (`scripts/RunT3ZeroExtraBoundaryBenchmark.py`) and the existing Stage31/Stage36 numerical modules.

## PowerShell (project root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Scalar_Conclusion_and_LowEnergy_Batch8.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_zero_boundary_low_energy_batch8.py
python -m scripts.RunT3ZeroBoundaryLowEnergyC5 --all
Get-Content -Encoding UTF8 output\scalar_relevance\zero_boundary_low_energy\summary.md
```

If `output/scalar_relevance/zero_boundary_benchmark/B_m1/result.json` and its C counterpart do not exist, first run:

```powershell
python -m scripts.RunT3ZeroExtraBoundaryBenchmark --all
```

`--model B_-1` or `--model C_-1` can be used for focused runs.

## New numerical output paths

- `output/scalar_relevance/zero_boundary_low_energy/B_m1/result.json`
- `output/scalar_relevance/zero_boundary_low_energy/C_m1/result.json`
- `output/scalar_relevance/zero_boundary_low_energy/summary.json`
- `output/scalar_relevance/zero_boundary_low_energy/summary.md`

Per-model reports contain C5 at both matching scales, the low-scale neutrino mass matrix, masses, mass-squared splittings, absolute PMNS entries and source SHA-256 fingerprints. They also record the preexisting Stage36 C5 matrix differences and the fact that they compare **different additional-quartic boundary scenarios**.

## Limitations (must accompany quoted numbers)

- This uses an **existing** hierarchical, model-specific fixed-order one-loop matched Weinberg artifact. No extra-scalar-quartic one-loop Matchete matching is performed.
- Native candidate RGBeta additional-quartic normalization is not physically certified against Matchete. The run uses the evolved baseline \(\lambda_{T3}\), not a fabricated mapping of the new quartics into the matching coefficient.
- Extra quadratic mixing is checked to be numerically negligible; generic nonzero mass mixing is not evaluated.
- No EFT1 \(C_{12}\) or additional Wilson-operator running, finite threshold corrections, or new fit. The result is a **conditional benchmark estimate** rather than a new complete physical prediction.
- Original Stage36 uses Stage29 *nonzero* extra UV quartics. Its difference from this zero-extra scenario includes backreaction on retained coupling running, and must not be called the error from neglecting the direct extra C5 quartic.
- The existing Stage36 common grouped heavy/scalar threshold approximations remain in force.

The Batch 8 test suite runs offline with fixtures and mock metadata; the full B/C RG + matched-C5 + SM evolution must execute in the user's actual local repository before numerical values are claimed.
