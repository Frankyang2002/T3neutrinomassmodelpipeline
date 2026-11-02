# Workstream 2 / V5: exact source Matchete↔RGBeta quartic dictionary audit

This is a **source-level** all-16 extension of V4. It does **not** constitute a
physical 24-channel coupling dictionary or an EFT1 Wilson-coefficient solver.

## Scope and source provenance

Consumes the actual original Matchete `ScalarQuarticTerms` and `CGRegistry`
from `output/full16_symbolic/T3_<family>_alpha_<tag>/data/eft1_after_F_scalar_quartic_seed.json`,
the local candidate RGBeta registry and `AddQuartic` declarations, and the
successful V4 reports. V4 file hashes, the Matchete quartic seed, and local
RGBeta Wolfram source hashes are checked against the previous audit. A changed
source is not silently treated as physically normalized.

The source inventory includes **all 50 additional RG-closed quartic entries**
(24 HHSS and 26 SSSS), rather than only the 24 HHSS named vertices. The native
Matchete tensor and term census is independently checked for all 16 models.

### Newly implemented

- `source_cg_v5.py`: exact Matchete CSR `SparseArray` decoder, full tensor
  coordinate, representation, sign, radical and complex-phase preservation;
  strict source field and named CG usage checks. No approximate floats.
- `matchete_rgbeta_dictionary_v5.py`: exact full-vector tensor-basis solver,
  source-sector identity matching, and candidate CG index/dual-map search.
  Candidate orientation searches are **not** physical equivalences.
- `named_dictionary_batch_v5.py`: all-16 source crosswalk and fail-closed
  provenance; checks previous V4 and preserves no-overwrite output policy.
- `scripts/RunLOMatcheteRGBetaDictionary.py`: CLI.

The exact tensor equation solved when both bases exist is

\[
 T_{\rm target}=T_{\rm source}\,M,
\qquad T_{\rm source}M-T_{\rm target}=0.
\]

No least-squares approximations, smoothing, assumed normalization factors or
unjustified signs are accepted. A *formal* CG tensor equation is not sufficient
to give a physical coupling dictionary: the original Matchete `CG[...]
Index[...]` slot association, `PlusHc` sign/factor, and RGBeta's *internal*
`AddQuartic` factor must still be verified.

### Critical source result

The read-only November 2026 repository audit of the 16 `full16_symbolic`
Matchete seeds found **no direct native source term with the field monomial
of any of the 24 additional HHSS operators**. Hence V5 cannot certify 24
Matchete↔RGBeta extra-channel coupling conversions solely from these seeds.
The native seed needs to be enlarged with normalized extra operators or an
independent fully normalized exact Matchete-to-RGBeta tensor export supplied.

The generated RGBeta registry itself is labeled
`GeneratedCandidateBasisNotWolframVerified`. Verified source parsing is *not*
Wolfram verification, RG closure of that dictionary, physical running, or a
calculation of neutrino mass.

## Owned paths and run instructions

All new files are exclusively in:
`Numerical/lo_certification/couplings/`, `scripts/`, `tests/`, and
`output/general_observable/lo_certification_couplings/`.

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WS2_Matchete_RGBeta_Dictionary_V5.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_lo_couplings_*.py' -v
python -m scripts.RunLOMatcheteRGBetaDictionary --repo-root . --run-name ws2_matchete_rgbeta_dictionary_v1 --v4-run-name ws2_named_hhss_v1
Get-Content 'output\general_observable\lo_certification_couplings\ws2_matchete_rgbeta_dictionary_v1\summary.json' -Raw
```

The runner does not overwrite previous outputs. If a model fails CG decoding
or source verification, it emits a local `failure.json` and exits nonzero for
an all-16 run. Do **not** bypass the failure by zero-filling a coupling.

### Output layout

```
output/general_observable/lo_certification_couplings/
    ws2_matchete_rgbeta_dictionary_v1/
        summary.json
        A_m4/matchete_rgbeta_exact_source_crosswalk.json
        A_m2/matchete_rgbeta_exact_source_crosswalk.json
        ... (16 models)
```

Per-model reports include complete extra-operator inventories, all Matchete
CG exact digests and norms, native field counts and couplings, matching-sector
flags, candidate orientations (if any), running coefficient presence at
both scales, explicit blockers, and source SHA-256s.

No physical seven-term or simultaneous scalar-matching solver input is
exported. The existing `WilsonTensorRGE`,
`GeneralEFT1C5Mixing`, and `IntegratedEFT1JointScalarBridge` remain untouched.

**Thresholds**: UV → integrate out **all fermions** at μF → EFT1 with **both
scalar multiplets** → integrate out **both scalars together** at μS → SMEFT.
