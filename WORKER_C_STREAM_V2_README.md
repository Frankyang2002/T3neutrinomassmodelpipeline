# Worker C — established-scale joint scalar matching handoff v2

**Install all files at the paths encoded in the ZIP. No existing production source is overwritten.** This release does not touch GitHub, the production coupled DOP853 solver, archived rough16 runs, or the old 3 TeV toy D(-2) benchmark.

## Physics scope and bookkeeping

Mandatory sequence:

    UV(10^7 GeV) -> integrate ALL fermions at muF=10^5 GeV
    -> EFT1(SM+S1+S2) -> integrate BOTH scalar multiplets at muS=10^3 GeV
    -> SMEFT -> 10^2 GeV

Once *actual evolved* Wilson tensors at `muS` are provided, the endpoint mode uses the already-implemented `worker_c_joint.match_joint_scalars` and the existing `core.match_c5`:

    C5_SMEFT(muS) = C5_EFT1_evolved(muS) + DeltaC5_LLSS_HHSS_bubble(muS)

There is **no second EFT1 leading-log term** in endpoint mode; inherited C5 is supplied already evolved. The matching includes the finite MSbar zero-momentum scalar bubble (all eigenpairs, 1/2 symmetry factor, both S1 and S2, complete HHSS tensor). This is not the complete finite scalar one-loop matching for arbitrary EFT1 operators. The signed Matchete GammaCC/Weyl conversion and full non-bubble scalar matching remain separate prerequisites.

## Files (project-relative)

- `Numerical/joint_scalar_matching/worker_c_stream_v2.py` — new; canonical endpoint validation, actual endpoint matching, all16 source/background survey, SHA-bound archive regression, missing-work inventory.
- `scripts/RunT3WorkerCStreamV2.py` — new; CLI.
- `tests/test_t3_worker_c_stream_v2.py` — new; 17 offline unit tests including integration mocks and rejection cases.
- `WORKER_C_STREAM_V2_README.md` — new.

Dependencies are existing repository modules, not replacements: `Numerical.joint_scalar_matching.adapter`, `.core`, `.large_batch`, `.fermion_boundary_bridge`, `.worker_c_joint` (the latter from earlier Worker C package). Runtime also needs NumPy.

## Exact PowerShell commands

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_worker_C_stream_V2_20261112.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p "test_t3_worker_c_stream_v2.py" -v
python -m scripts.RunT3WorkerCStreamV2 --remaining

# First: focused A(-4) and D(-2) historical-source reevaluation, using the
# established 1e5 -> 1e3 GeV scalar threshold, NOT the earlier 4/3 TeV toy.
python -m scripts.RunT3WorkerCStreamV2 --survey-all16 --models A_-4 D_-2

# Next: evaluate every available 16-model source/background; any failure
# produces an individual BLOCKED report and nonzero exit code.
python -m scripts.RunT3WorkerCStreamV2 --survey-all16
```

Default source/boundary inputs are **read-only**:

- `output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/{slug}/boundary.json`
- `output/general_observable/lo_joint_scalar_matching/exploratory_all16_v6_native_fix/numeric_backgrounds/{slug}.json` plus sibling `{slug}.npz`
- Historical comparison: `output/general_observable/lo_joint_scalar_matching/exploratory_all16_v6_native_fix/C5_per_unknown_NGamma_exploratory/{slug}.json`

Where `slug(D_-2)=D_m2`, `slug(A_-4)=A_m4`. Override using `--source-dir`, `--background-dir`, `--archive-dir`, or `--no-archive`. The archive is compared only when all three SHA fields and both scales are identical. Every run creates a unique fresh folder under `output/audits/worker_c_all16_*`. It saves per-model reports and `summary.json` with conditional C5 matrices and verification metrics; it never writes back to the source folders. `PASS_CONDITIONAL` means the formula/inputs computed, **not** that a physical neutrino mass or restricted Yukawa completion is certified.

### Evolved-Wilson endpoint handoff

Generate an **incomplete** contract (no fictional matrices) for Worker B / D:

```powershell
python -m scripts.RunT3WorkerCStreamV2 --contract D_-2
```

The generated template requires the existing canonical adapter's NPZ arrays at `muS=1000 GeV`:

| NPZ key | Shape | Source |
|---|---|---|
| `mass_squared_gev2` | `(n,n)` | Worker B full real scalar mass matrix at muS |
| `lambda_HHSS` | `(4,4,n,n)` | Worker B **source-retained** complete HHSS tensor at muS |
| `C_LLSS` | `(3,3,2,2,n,n)` | Worker D evolved 7-term EFT1 Wilson sector at muS |
| `C5_EFT1` | `(3,3)` | Worker D inherited/evolved Weinberg coefficient at muS |
| `lambda_HHSS_channels` | optional `(k,4,4,n,n)` | Worker B additive quartic channel decomposition |

`n=2(dS1+dS2)`, with both scalar multiplets in one mass matrix. The endpoint JSON is adapter schema v1 plus explicit `yukawa_sector = Y1_Y2_only_no_independent_conjugate`, ownership and scale declarations. The endpoint fails closed if a stage-one `muF` Wilson coefficient is mislabeled as a `muS` endpoint. B(-1) / C(-1) require the additional `mixed_scalar_mass_basis_verified` flag. These flags must be set **from verified source evidence**, not just toggled to satisfy validation. Use the actual local evolved Wilson array exporter; this package does not reconstruct arbitrary source-to-Weyl normalization.

Then run (once the real metadata and NPZ are exported):

```powershell
python -m scripts.RunT3WorkerCStreamV2 --endpoint "path\to\actual_D_m2.json"
```

Endpoint output: finite complex-symmetric `c5_inherited_gev_inv`, `c5_finite_scalar_LLSS_HHSS_gev_inv`, `c5_total_one_bubble_only_gev_inv`, eigenpair and quartic-channel contributions, mass eigenvalues, Ward residual and all source SHA values. `complete_scalar_matching_certified` remains `false` even if all provided sources are validated: no missing non-bubble hard matching is invented.

## Verified and not verified

**Verified upstream source inspection:** `core.py`, `adapter.py`, `large_batch.py`, `fermion_boundary_bridge.py`, the archived all16 conditional background/release filenames, D(-2) source and tensor metadata. These files already implement the scalar-bubble loop and independent conditional source/flavour decompositions.

**Verified in this package:** 17 offline tests for the new code, including scale exclusions, complex 3x3 encoding, nonoverwriting outputs, source/basis/yukawa guardrails, inherited C5 once, eigenpair summation, mock survey end-to-end, hash-equivalent archive regression.

**Not yet verified in this environment:** Execution against the **user's latest local** 16-model NPZs or evolving seven-term solver endpoint. No `muS` actual evolved Wilson endpoint was supplied in the conversation. The survey uses an existing conditional *muF source* with a fixed-order generated leading log, not the coupled DOP853-evolved Wilson state. In particular it cannot replace the user's successful 16/16 / 256/256 provisional runs.

## What is left

1. **Worker A:** final **finite renormalized** hard one-loop fermion-boundary coefficients, including UV-minus-EFT1 soft subtraction, not the raw Matchete Laurent regular part.
2. **Worker B:** same-scale all16 scalar tensor/mass production in canonical real basis; retain established HHSS source tensor; include exceptional quadratic and B/C mixed terms. Explicit source hashes and basis map.
3. **Worker D:** export actual evolved `C_LLSS(muS)` and `C5_EFT1(muS)` in matching convention from the coupled seven-term solver. Ensure all16 strictly restricted `Y1,Y2` source (legacy `y2Conj` evidence cannot be ignored).
4. **Worker C:** identify/add *independent finite scalar threshold* corrections other than the already implemented bubble (if supported), with counterterms, scheme, and overlap subtraction. The old D(-2) **3 TeV** finite-part diagnostic is a different illustrative source and must NOT be combined with the 1 TeV pipeline by numerical subtraction.
5. **Worker A/B/C:** signed Matchete `GammaCC` / Weyl / CG normalization bridge. For rough conditional results the unknown signed factor may remain explicitly symbolic/conditional.
6. **Worker D:** only after handoffs, couple to production solver, rerun A(-4), D(-2), then 16 models and 256 stage-freeze cases; regenerate comparisons without altering archived outputs.

The next numerical handoff to prioritize is **the actual evolved EFT1 Wilson endpoint at 1000 GeV**. A source-only conditional survey is useful for regression, but cannot substitute for that endpoint.
