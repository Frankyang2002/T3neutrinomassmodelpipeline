# Batch 16 — UV-running B(-1)/C(-1) hierarchy-safe neutrino normalization probes

## Purpose and verified source interfaces

Batch 15's smaller-Yukawa scan found seven *mass-screen-passing* candidate points, not neutrino oscillation fits:

- B(-1): common Yukawa factors `0.005,0.01,0.015,0.02`.
- C(-1): common Yukawa factors `0.005,0.01,0.015`.

With the original `lambdaT3=5e-9`, scaling both Yukawa matrices can strongly suppress the one-loop neutrino mass. Test an **illustrative normalization restoration** by varying

\[
(y_1,y_2)_{UV}=f(y_1,y_2)_{\rm seed},\qquad
\lambda_{T3,UV}=\frac{5\times10^{-9}}{f^2}\,k,
\]

where `k=0.1,1,10` by default. The candidate seven extra quartics at the UV boundary are still `0.1*abs(lambdaT3)`, so they scale with the topology coupling. The UV input scalar masses and heavy-fermion masses are unchanged.

**This is a probe, not a refit.** Holding each Yukawa matrix texture fixed means its neutrino mixing pattern and mass-ratio pattern may still be wrong. Scaling lambda principally changes neutrino mass normalization; optimizing 18 Yukawas and lambda is the subsequent step if needed. No first-principles relationship `lambda ∝ 1/f²` is asserted: it is a seed heuristic, and the RG equations are re-integrated for every point.

## Inputs and provenance

Reads original historical `small_uv_running_5e-9` files and source-hash-verified prior `general_fresh_v1_B_m1_cont`, `general_fresh_v1_C_m1_cont` outputs via Batch 14 `fit_source`. Two optional UV Yukawa texture sources are evaluated:

- `saved_poor_fit`: the saved general-scalar exploratory best Yukawa matrices;
- `historical_fit`: the original restricted-potential UV-running fit Yukawas.

Each trial changes copies of only the 18 real Yukawa entries and the `lambdaT3` UV boundary. `source['lambdaT3']` **and** `source['cfg']['base_state']['ordinary']['lambdaT3']` are changed together; this is essential because the existing Batch 9 boundary builder reads the former and Stage 31 state initialization reads the latter.

The established Batch 11 logic supplies candidate one-loop RGBeta UV/EFT1, the heavy mass-basis rotation, the saved one-scale ordered Matchete C5 with the **conditional** Stage23/26 minus-sign map, and SM+Weinberg running to 100 GeV. All five NuFIT observables and their approximate chi-square are recomputed for every eligible candidate.

Before evaluating C5, the runner requires positive 2×2 running scalar mass eigenvalues and both:

\[
\rho_F=\frac{m_{S,\max}(10^5\,\mathrm{GeV})}{M_{F,\min}(10^5\,\mathrm{GeV})}\le 0.5,
\qquad m_S(10^3\,\mathrm{GeV})\in[0.5,5]\,\mathrm{TeV}.
\]

Failures at the fermion threshold skip EFT1; failures at the scalar threshold skip the neutrino calculation. The one-scale C5 is never presented as a certified physical or sequential matched result. Running masses are not pole masses, and the scalar window is a numerical benchmark convention, not a validated decoupling prescription.

A candidate with chi-square <=5.86 is a conditional seed-level fit only. Because no optimization occurs, neither `chi2` nor any parameter is promised to improve.

## Files (new; no existing file replacements)

- `scripts/ProbeT3BCUVNeutrinoNormalization.py` — full standalone numerical probe runner.
- `tests/python/test_t3_bc_uv_neutrino_normalization.py` — 25 offline contract and mock-interface tests.
- `docs/T3_BC_UV_NEUTRINO_NORMALIZATION_BATCH16.md` — this file.

All inputs are read-only. The results are new-only under:

```
output/studies/general_scalar_extended/diagnostics/<new-run-id>/
  probe.json
  probe.csv
  summary.md
```

`probe.json` records source hashes, full trial predictions (including C5 matrices where computed), the threshold mass screen, and explicit noncertification flags. `probe.csv` is a compact flat table.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_BC_UV_Neutrino_Normalization_Batch16.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_bc_uv_neutrino_normalization.py

# First: inexpensive 14-trial fixed-normalization probe (two Yukawa textures).
python -m scripts.ProbeT3BCUVNeutrinoNormalization `
  --run-id uv_neutrino_norm_k1_v1 `
  --lambda-multipliers "1"

Get-Content -Encoding UTF8 output\studies\general_scalar_extended\diagnostics\uv_neutrino_norm_k1_v1\summary.md

# If needed, a separate 42-trial lambda sensitivity run; never reuses an ID.
python -m scripts.ProbeT3BCUVNeutrinoNormalization `
  --run-id uv_neutrino_norm_grid_v1 `
  --lambda-multipliers "0.1,1,10"

Get-Content -Encoding UTF8 output\studies\general_scalar_extended\diagnostics\uv_neutrino_norm_grid_v1\summary.md
```

To run fewer trials, select `--model B_-1 --seed-source historical_fit`, or supply `--yukawa-factors "0.005,0.01"` with `--model`. Use a unique run ID for every new set.

## Validation boundary

Offline tests cover boundary consistency, seed ordering, mass-screen gating, skip-C5 on invalid hierarchy, successful UV/EFT1/C5/low-energy mock-interface sequence, input parsing and nonoverwrite guards. No live RGE/C5 calculation is claimed to have run in the offline staging tree: the full numerical services and saved parameter inputs require the user's Windows repository.

If the new candidate chi-square remains large, the next step is a **hierarchy-constrained 18-Yukawa plus lambda refinement** starting from the best mass-safe point, not another unconstrained historical-UV fit.
