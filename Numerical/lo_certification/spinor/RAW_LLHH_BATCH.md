# Raw LLHH sector inventory after all-16 simultaneous fermion matching

**File destination:** `Numerical/lo_certification/spinor/RAW_LLHH_BATCH.md`

The preceding `RunLOSpinorHardThreshold.py --run-id native01 --phase-probe` returned 16/16 native F-only one-loop results at one common `mu_F=100000 GeV`, with both BSM scalar multiplets still active and SM fermions retained. This step **does not repeat Matchete**; it classifies every additive term in those saved raw expressions.

`ExtractLOSpinorRawLLHH.py` verifies the original native summary, the 16 per-model inspection reports, native logs, raw-result SHA256 hashes, and the original heavy/light threshold labels. The output is a new immutable directory within `output/general_observable/lo_certification_spinor/`.

For each model the new `terms.json` keeps every matched **candidate** with its exact Matchete InputForm, source additive-term index, SHA256, all referenced coupling names, chiral projector (PL or PR), explicit epsilon pole, and logarithm marker. A separate record retains any matching LLHH structures without `hbar`, so they are not silently counted in the loop sector. The `PL_raw_llhh.inputform.txt` and `PR_raw_llhh.inputform.txt` exports keep the Hermitian-conjugate sectors separate. They are not added to one another. All small-coupling terms are retained.

Exact candidate rule: exactly two `Field[l,Fermion,...]` and two `Field[H,Scalar,...]` *including field powers*, no other external fields, one `DiracProduct[GammaCC,Proj[±1]]`, and explicit `hbar` for the one-loop inventory. `Field[NewScalar*]`, split neutral scalar labels and additional SM external fields do not enter this **LLHH** sub-sector. This is **not** a complete operator basis projection of all EFT1 loop terms; the original full one-loop file remains authoritative and is never modified.

## Install and run (PowerShell)

At `C:\thesistemp\T3neutrinomassmodelpipeline` after extracting the ZIP:

```powershell
$env:PYTHONPATH = (Get-Location).Path
python -X utf8 -m unittest discover -s tests -p 'test_lo_spinor_*.py' -v
python -X utf8 scripts\ExtractLOSpinorRawLLHH.py --check-only
python -X utf8 scripts\ExtractLOSpinorRawLLHH.py --run-id inventory01
```

For an alternative native batch, pass `--native-run output\general_observable\lo_certification_spinor\<your-native-run>`. The source must be an existing complete all-16 native run in the exclusive spinor output directory. Never overwrite `inventory01`: choose a distinct `--run-id` if making a new snapshot.

Outputs:

```text
output/general_observable/lo_certification_spinor/raw_llhh_e515d82d118284ba_inventory01/
    summary.json
    A_m4/terms.json
    A_m4/PL_raw_llhh.inputform.txt
    A_m4/PR_raw_llhh.inputform.txt
    ... one folder per model (16)
```

## Why this is **not** the physical C5

The original general-14 `eft1_after_F_one_loop_inputform.txt` files include terms without `hbar`; the raw B/C LLHH candidate sectors have explicit `1/ϵ` terms. Simply dropping pole terms, setting `ϵ -> ∞`, retaining only a finite marker, or adding PL/PR would be wrong. Properly identify the threshold's hard region and EFT1 soft subtraction, choose the renormalization convention, handle counterterms and wave-function factors, and match to the master Weinberg operator with the independent Matchete-to-Weyl spinor phase. All physical normalization keys remain JSON `null`/`false` until those checks are independently demonstrated.

The required EFT ordering is unchanged: **UV → integrate all BSM fermions together at `mu_F` → EFT1 with both BSM scalar multiplets → integrate both scalars together at `mu_S` → SMEFT**.

## Verified upstream evidence and limitations

- The Windows log supplied by the user reports 16/16 native runs and 6/6 abstract GammaCC identities, **not** a numerical physical GammaCC phase.
- A read-only GitHub audit of 9 available raw native expressions on 2026-11-03 finds 6 LLHH candidates in each of A_-4, B_-3, B_+1, C_-3, C_+1; 12 in each of A_-2, A_+0; and 204 in each of B_-1, C_-1. The read-only repository snapshot did not expose the D/E raw expressions for independent inspection during this batch; their counts must come from the local Windows extraction run.
- The Python tests exercise full 16-model provenance and source hashing against synthetic fixtures; **local Windows execution against all 16 real files is required** to confirm this new extractor.
- No physical or renormalized first-threshold Wilson coefficient is certified by this package.
