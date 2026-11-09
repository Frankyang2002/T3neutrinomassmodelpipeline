# Workstream C v2 — real A/B handoff, preflight, counting comparison

**Purpose:** consume the *actual* Workstream A full evolved EFT1 Wilson endpoint and Workstream B full scalar/SM couplings, then perform the **one common** scalar threshold and call the project's one-loop SM Weinberg ODE. This is an additive C-only batch on top of `README_C.md`/C v1. No A/B-owned modules, source files, native scalar potential, or GitHub contents are modified.

**Inspected public source:** `Frankyang2002/T3neutrinomassmodelpipeline` `main` at `a944e0030b8ef24b23342d99e0d8f33cef3cffd5` on 9 Nov 2026. Neither production A endpoint modules under `Numerical/parallel_t3/wilson/` nor B endpoint modules under `Numerical/parallel_t3/background/` were present in that public snapshot. The user's local repository may differ. This batch deliberately **cannot synthesize real A/B data**.

## Exact EFT order

```
UV --[integrate ALL fermions at MF]--> EFT1 (BOTH scalars, full one-loop dim-5 evolution A)
    --[at ONE MS, diagonalize and integrate ALL heavy scalar eigenstates together]-->
    SMEFT C5(MS) --> existing one-loop SM+Weinberg ODE --> C5(mu_low), m_nu.
```

At matching, the RGE-improved mode uses

```
C5_SMEFT(MS) = C5_EFT1_evolved(MS)
              + Delta_C5_LLSS-HHSS_scalar_bubble(MS)
              + independent_finitematching_C5(MS) [ONLY if supplied].
```

The C5 already induced by EFT1 running is not added a second time as an old LLSS→C5 leading log. Independent missing diagrams remain labelled **uncomputed**. The scalar bubble uses the full B-supplied real quartic tensor's HHSS sector and transforms the complete supplied LLSS tensor by the explicit A→B real rotation and B mass-eigenvectors. Both scalars are always integrated at the same MS; no intermediate scalar EFT is created.

### Critical v2 strict-loop improvement

The previous C v1 `strict_one_loop` mode used `C_psi2phi2_tree_mf` in the scalar bubble but took the evolved (possibly resummed) C5 from the generic A endpoint. It therefore only made the **bubble** strict. **C v2 corrects this:**

- `rge_improved`: A NPZ `C_psi2phi2_ms` provides the already evolved full tensor, including C5 and LLSS.
- `strict_one_loop`: additionally needs A NPZ `C_psi2phi2_ms_strict_1loop` (the **full** 1-loop-truncated C at MS), A NPZ `C_psi2phi2_tree_mf` (the **full** order-zero tensor at MF), and A JSON `"strict_1loop_endpoint_at_ms": true`. C5 is projected from `C_psi2phi2_ms_strict_1loop` and the one-loop scalar bubble is contracted with `C_psi2phi2_tree_mf`. No synthetic or guessed A one-loop coefficient substitutes for it.
- A optional `C5_EFT1_ms` remains a **redundancy check against the RGE-improved full tensor**, not an additive term.
- Both modes use **B's supplied mass/quartic parameters at MS**, without a perturbative expansion of their own one-loop RG trajectories. Therefore `strict_one_loop` is *strict in A Wilson/bubble counting* but not a proof of a fully expanded renormalisable-coupling perturbation series. C itself remains an incomplete independent one-loop scalar threshold.

All optional A/B quality flags default to **uncertified**. The strict C5 contribution's origin, both MS EFT1 C5 projections and the bubble's order-zero/evolved origin are recorded in C's manifest and `C5_smeft_at_MS.npz`.

## A/B endpoint contract (no guessed conversion)

All arrays are real `float64` or complex `complex128` NumPy arrays in sibling NPZ with `allow_pickle=False`. They must be already projected into the master real-field tensor conventions and use **the same matching state**. Paths below are **examples of output layout**, not evidence that they already exist locally:

```
output/parallel_t3/wilson/<A_run_id>/A_m4/wilson.json
output/parallel_t3/wilson/<A_run_id>/A_m4/wilson_at_MS.npz
output/parallel_t3/background/<B_run_id>/A_m4/background.json
output/parallel_t3/background/<B_run_id>/A_m4/background.npz
```

Each manifest's `tensor_npz` names its adjacent NPZ. For all16, use slugs `A_m4`, `A_m2`, `A_p0`, etc. A/B matching fields must agree **exactly** on `model`, `matching_context_id`, `matched_seed_id`, `lagrangian_source_id`, `mu_f_gev`, `mu_s_gev`, scheme `MSbar`, and scalar-basis bridge ID. `tensor_npz_sha256` is checked if present. An old pre-portal A(-4) source must not be combined with its corrected source merely because the numerical scales match.

A JSON required fields and NPZ keys are specified in `README_C.md` and `contract_schema_v1.json`. In addition, C v2 strict requires:

```json
{
  "strict_1loop_endpoint_at_ms": true,
  "strict_1loop_origin": "A_EFT1_one_loop_truncated_at_MS"
}
```

And `wilson_at_MS.npz` must contain these three independent source tensors:

| NPZ key | Shape | Meaning |
|---|---|---|
| `C_psi2phi2_ms` | `(Nweyl,Nweyl,4+n,4+n)` | Full **RGE-improved** evolved tensor, includes EFT1 C5, all scalar components |
| `C_psi2phi2_ms_strict_1loop` | same | One-loop **truncated** full A endpoint at MS, for strict C5 projection |
| `C_psi2phi2_tree_mf` | same | Full zeroth-order A Wilson at MF for strict bubble |

B NPZ must contain `mass_squared_gev2_ms`, `R_wilson_from_background` with `s_A = R s_B`, `quartic_real_full_ms` of shape `(4+n)^4`, `ye`, `yu`, `yd`, and if available `quartic_real_channels_ms` summing to the full B tensor. B JSON must explicitly indicate `quartic_inventory_status="full_tensor"` and provide the four SM couplings at MS. **No reduced portal is accepted**. The A(-4) anti-Hermitian Inv1 correction is preserved upstream by accepting the **complete** B potential; C never reconstructs coupling tensors from names.

If A/B produce *different JSON/NPZ formats*, do **not** rename or pretend the formats coincide. Have A/B emit a canonical C input pair or provide the actual current files for a C-owned conversion adapter. C preflight intentionally fails on unrecognized formats.

## Installation (PowerShell from the repository root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$env:USERPROFILE\Downloads\T3_Workstream_C_RealAB_Integration_v2.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
.\.venv\Scripts\python.exe -m pytest -q `
  tests\test_t3_parallel_c_threshold.py `
  tests\test_t3_parallel_c_real_integration.py
```

### Read-only availability check

```powershell
Get-ChildItem output\parallel_t3\wilson -Recurse -Filter wilson.json -ErrorAction SilentlyContinue | Select-Object FullName
Get-ChildItem output\parallel_t3\background -Recurse -Filter background.json -ErrorAction SilentlyContinue | Select-Object FullName
```

Use the actual A/B **run roots** from these results, not arbitrary old benchmark outputs. Probe **one real A(-4)** before full 16:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ParallelCRealIntegration.py `
  --wilson-root "output\parallel_t3\wilson\<A_run_id>" `
  --background-root "output\parallel_t3\background\<B_run_id>" `
  --models A_-4 --mode both `
  --output output\parallel_t3\threshold\real_ab_v2
```

The **default is PREVIEW**. This only reads A/B and produces the C-owned `real_integration_summary.json` readiness report. It checks all tensor dimensions, full Wilson pair symmetry, consistent source hashes/seed IDs, canonical C5 projection, both scalar mass spectra, full quartic/channel inventory, physical positivity, the single-scale scalar-bubble projection and valid SM couplings. It does not perform low-scale evolution in preview mode. Unavailable or stale A/B data are marked BLOCKED rather than replaced.

Once `ready_conditional` is returned, add `--run` to the same command to write C results and call the actual existing SM Weinberg runner:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ParallelCRealIntegration.py `
  --wilson-root "output\parallel_t3\wilson\<A_run_id>" `
  --background-root "output\parallel_t3\background\<B_run_id>" `
  --models A_-4 --mode both --run `
  --output output\parallel_t3\threshold\real_ab_v2
```

To process all 16, **omit `--models A_-4`**. If an A strict endpoint is unavailable, RGE-improved can be preflight-ready while strict is BLOCKED; v2 reports both independently. For a C-only diagnostic without the SMEFT runner, add `--match-only` (this does **not** produce a low-energy evolved prediction). Use `--overwrite` to replace existing C outputs.

### Output layout and counts

With `--mode both --run`, the individual modes are separated:

```
output/parallel_t3/threshold/real_ab_v2/
    real_integration_summary.json
    A_m4/
        strict_one_loop/
            manifest.json
            scalar_mass_basis.npz
            C5_smeft_at_MS.npz
            C5_smeft_low.npz
        rge_improved/
            manifest.json
            scalar_mass_basis.npz
            C5_smeft_at_MS.npz
            C5_smeft_low.npz
        counting_comparison.json
```

`counting_comparison.json` gives Frobenius norms of the strict-to-RGE-improved difference in the scalar bubble and combined C5(MS). The LLSS bubble difference is **higher order in the fixed-order counting** (assuming A supplied a genuine one-loop-truncated endpoint). Differences in EFT1 C5 can also reflect A's high-order resummation. These are conditional numerical comparisons, not a matching-scale independence test.

### Status and unresolved physics

- **Verified in C v2**: all16 synthetic A/B contract preflight; positive eigenmass/rotation and quartic checks; strict vs RGE-improved output separation; no duplicate direct leading log; explicit preflight failures for missing/stale/incorrect A/B input. These are implementation tests, **not native Matchete/RGBeta validation**.
- **Verified from public repo source**: scalar `core.py`, native mixed-doublet and self-bilinear conventions, existing SM+Weinberg ODE, A(-4) portal phase fix.
- **Waiting on actual A/B outputs**: canonical full evolved Wilson `C_psi2phi2_ms`, optional one-loop-truncated strict C, complete B full scalar Hessian/quartic and mapping, aligned SHA/basis IDs. The public HEAD inspected does not provide these production A/B endpoints.
- **Still uncomputed by C**: other one-loop finite scalar diagrams, wavefunction/EOM terms, any additional dimension-five-class threshold matching; the full one-loop scalar matching claim remains **false**. An EWSB-resummed split spectrum need not project exactly onto the unbroken SMEFT Weinberg operator (higher-dimensional contamination).
- **Not certified by C**: signed GammaCC→Weyl factor, Matchete↔RGBeta full scalar operator normalization, one-loop EFT1 dimension-five ADM/closure, physical neutrino masses/pole masses.
