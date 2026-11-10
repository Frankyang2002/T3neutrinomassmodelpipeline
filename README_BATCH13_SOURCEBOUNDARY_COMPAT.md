# Batch 13 / Batch 15 compatibility runner

Verified from GitHub `main` on 2026-11-11: 
`Numerical.joint_scalar_matching.fermion_boundary_bridge.SourceBoundary`
has a `path` field, while the CG-enabled `two_model_master_wilson_batch13.py`
reads `seed_path`. This package avoids editing the CG-enabled source file.

New project path:

`scripts/RunT3MasterWilsonBatch13CGFix.py`

From the project root in PowerShell:

```powershell
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3MasterWilsonBatch13CGFix --all --output output/parallel_t3/thesis_master_wilson_batch13_cgfix_v3 --debug
if ($LASTEXITCODE -ne 0) { throw 'Batch 13 failed' }
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3PhysicalClosureBatch14 --all --generator-details --output output/parallel_t3/thesis_physical_closure_batch14_cgfix_v2 --debug
```

If the local repository has already corrected the attribute name, the runner
simply delegates to the original command. Both computed physics results and
certification status remain conditional; this package addresses only the
Python API mismatch.
