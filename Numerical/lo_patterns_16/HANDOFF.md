# Chat A pre-run handoff (deliverable source package)

## Provenance

- Remote `main` verified on 2026-11-03: `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4`.
- Repository source inspected via read-only GitHub connector; the actual Windows local working tree is **not mounted** in this environment. The runner will record the actual local `git rev-parse HEAD`, dirty status, and SHA256 of the frozen inputs and source audits.
- `Numerical/direct_one_loop/core.py` observed remote Git blob: `abbf261fbf2ca2561482d3e53b84452d1654610e`.
- `Numerical/lo_certification/direct_uv/native_uv_loop_reference.py` observed remote Git blob: `f48128b1af6e48038f1658f98f7d3fa4e7c03092`.
- `Numerical/lo_certification/direct_uv/channels.py` observed remote Git blob: `930e9516d537f49b1a8438cefda13314e0bcbe15`.
- `Numerical/lo_certification/direct_uv/native_amplitude_contract.py` observed remote Git blob: `6aebf081f6a9a8a89796356b2141f3972420afc7`.
- Source `charge_resolved_loop_weights.json` remote blob: `eb9cf25563b541ff92660ddee839f93aa1703ed2`.
- Git **blob IDs are not SHA256**; actual run computes SHA256 from local bytes.

## Computed vs pending

- Actual local Windows frozen benchmark run: **NOT EXECUTED HERE**; those files cannot be accessed directly in this sandbox. Do not interpret the enclosed synthetic checks or algebraic channel table as a 16/16 real run.
- Source convention: candidate Weyl `y1` uses conjugate scalar channel (`Z`), candidate `y2` uses field scalar channel (`Y`). Source CG audit contraction: `Y1 * Y2 * conjugate(T3MixCG)` in common weight indices; E has two channels, A–D one.
- Conditional exact algorithm: production `evaluate` (Takagi/real Hessian) against independent Hermitian spectral reference, UV-pole checks and original-charge partitions.
- Analytical lead: cross-S1/S2 scalar mixing enters through divided difference of finite B0; one-Yukawa rescalings linear if same-vertex terms vanish; lambda mixing linear only at first order with odd higher corrections.
- Missing for physical certificate: B1 four-spinor→left-Weyl phases, mass/sign/symmetry factors; B2 RGBeta/Matchete SU2 and scalar coupling normalisation; a coordinator-verified physical dictionary. The installed CLI rejects unverified certificate files.
- Claims safe to integrate: **source-conditional diagnostic pipeline only**; not physically normalized neutrino masses.

## Output after user executes on Windows

The runner writes new `output/general_observable/lo_pattern_study16/provisional_01/` including per-model JSON, all16 CSV/JSON, scans CSV, native CG weights JSON, common-kernel comparison, plots when available, and a run-specific `HANDOFF.md` + `HANDOFF.json`. An existing nonempty output directory is never overwritten.
