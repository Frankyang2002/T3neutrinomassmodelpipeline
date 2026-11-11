# Batch 3B source-path correction

Full replacement file: `scripts/BuildT3ScalarQuarticResponseA4D2.py`.

The existing `Numerical.joint_scalar_matching.fermion_boundary_bridge.SourceBoundary` dataclass names the original source path `path`, not `seed_path`. The runner previously accessed a nonexistent attribute twice. This version uses `Path(boundary.path)` for both JSON loading and SHA-256 provenance and verifies the file exists before reading it.

The numerical response evaluator (`Numerical/running/T3ScalarQuarticLinearResponse.py`) is unchanged, as are the scalar physics conventions and all prior result files. This is a compatibility fix, not physical scalar normalization certification.
