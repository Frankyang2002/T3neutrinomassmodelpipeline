# Batch 1: T3 all-16 general source/basis contract

Baseline read-only GitHub `main`: `ed8f6ba13c9a803e698ac155d24d8c025c53bc79`, 9 November 2026.
All files in the ZIP are **new paths**, not changes to the existing `ScalarPotential.wl`.

## Scope

- Audits the current local registry for exactly the 16 ordinary T3 model points and its 50 extra quartic slots (38 new RGBeta vertices, 12 reused declarations). Checks representation/hypercharge and Z2 for each registered quartic, and checks new RGBeta declarations are in the local generated source.
- Records SHA-256 of the actual local source and `git status`/HEAD; compares selected Git blobs to the independently read GitHub baseline, without overwriting local work.
- Optional execution of the existing native one-loop general UV/EFT1 quadratic RGBeta generators for **all 16** and native **all-fermions-first** Matchete matching for the **14 nonmixed** models. Both scalars remain light at MF. `B_-1/C_-1` are explicitly held for their separate mass-eigen-doublet backend; no fake ordinary-source substitutions.
- The A(-4) correction is regression audited against the **actual matched** EFT1 `tree_inputform` and native validation log, not mere source symbols. It accepts a specified corrected pre-existing local A(-4) result directory.
- Optional all-model native portal Hermiticity census; optional Stage 22 extra-quartic Matchete component and RGBeta *declaration* exports.
- Generic real-parameter scalar quartic dictionary solver `T3Batch1TensorMap.py` using explicit component polynomials, random complex-field contractions, and holdout residuals. The existing RGBeta declaration exporter is *not* sufficient: actual native RGBeta **components** must be provided to activate the mapping. No guessed dictionary is output.
- Records **physics blockers** on GammaCC-to-Weyl signed normalization, seven-vs-six ADM, and mixed/phase mass bases. Does not declare physical `C5`.

## Install and run (PowerShell from project root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Batch1_General_Contract.zip' -DestinationPath . -Force

# Verify that your local edits are not accidentally overwritten
 git status --short
 git rev-parse HEAD

# Fast unit tests (NumPy and existing project environment)
.\.venv\Scripts\python.exe -X utf8 -m unittest discover -s tests -p 'test_t3_batch1_contract.py' -v

# First run: read-only audit of local sources and any existing batch output
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3Batch1GeneralContract --run-id batch1_inventory_v1

# Native 16-model renormalisable one-loop betas + 14 F-only tree thresholds
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3Batch1GeneralContract `
  --run-id batch1_native_v1 --run-native `
  --a4-source 'output\general_observable\fermion_first_hermitian_portal_v1\A_m4'

# Optional: exact Matchete portal-CG Hermiticity census across all 16
wolframscript -file scripts\RunT3Batch1PortalAudit.wl `
  output\batch1_general_contract\batch1_native_v1\native\portal_census.json

# Optional Stage22 matched extra-quartic tensor/declaration evidence
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3Batch1GeneralContract `
  --run-id batch1_stage22_v1 --export-stage22
```

Use a **new run-id** for each experiment; `--resume` allows explicit reruns. The all-model native run may take substantial local computation. If one model fails, the runner records its failure and continues through the model list.

## Outputs

`output/batch1_general_contract/RUN_ID/`:

- `summary.json` — 16 model statuses and precise blockers
- `provenance.json` — local SHA-256 hashes, Git HEAD and dirty files
- `registry_audit.json` — all-16 symmetry/extra-quartic inventory
- `models/A_m4.json` ... `models/E_p2.json` — one detailed status per model
- `rgbeta/{modeltag}/{uv,eft1}.json` — native 1-loop renormalisable betas if executed
- `matchete/{modeltag}/native_status.json` etc. — native fermion-only matching for the nonmixed cases
- `logs/{modeltag}/...` — UTF-8 logs, no PowerShell `Tee-Object` encoding issue
- `native/portal_census.json` — optional native Higgs-portal phase census

## Claim boundary and next implementation

The runner creates an *operational evidence and orchestration layer*, not an all-16 physical one-loop neutrino mass calculation. In particular, no numerical Matchete/RGBeta map is invented from RGBeta declaration strings, and no signed GammaCC-to-Weyl conversion is assumed. These remain machine-readable blockages pending explicit native components. The separate seven-term Wilson ODE is not silently promoted to the correct anomalous dimension. UV→all F out→EFT1 with S1 and S2→both scalars out together→SMEFT remains the enforced target for later batches.
