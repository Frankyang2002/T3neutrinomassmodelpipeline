# Workstream 4: executable delivery and evidence report

**Source checked**: public GitHub `main`, commit
`2c29e5ba2c0350b7e9be292e02ceb8bc7c52c0d3`, read only. No branch,
commit, push, or modifications to frozen/local source were made.

**Physical conclusion**: 16/16 existing restricted exact `p=0` UV loop
numerical results remain **provisional**; 0/16 independently
Matchete-to-Weyl source-certified. The native portal geometry has been
validated in 32/32 scalar sectors, but this is not a physical coupling map.

**All-16 existing matrices**: `github_provisional_matrix_reference.json`
contains each complete symmetric real+imag `3x3` matrix in GeV, the three
Takagi masses in eV, the cancellation ratio, the old result path and the
source Git commit. All 16 matrices were independently checked with direct
NumPy SVD: maximum relative mass-spectrum discrepancy `8.59e-14`.
These are previous provisional outputs, **not** 16 new physical runs.

**Source status per model**: `source_coverage_all16.json`,
`existing_provisional_mass_summary.csv`. Exactly 6 assignments A, D, E
at alpha=-2,0 allow a new neutral-scalar conjugate Yukawa in the opt-in
source but no corresponding frozen coupling matrix exists. Of the 16
ordinary restricted models, 14 use Dirac doubled-Weyl fermions, 2 use
Majorana charge-pair metric (B_-1 and C_-1). The common physical-charge
intersection contains 1 channel per A-D model and 2 per E model,
20 source-allowed charge-channel labels in total. The source's
conditional native CG contraction `Y1*Y2*conjugate(T3MixCG)` gives total
group coefficients `-2*sqrt(3)/3` (A, B, D), `2/3` (C), and
`2*sqrt(6)/3` (E, two charge channels). These are **not** proven
physical spinor-normalised multiplicity factors. Relative phases remain
unfixed.

**New implementation**:
`Numerical/lo_certification/direct_uv/{certification.py,channels.py,runner.py}`
and `scripts/RunLODirectUVCertification.py`. No existing loop engine,
matching, RGE, or EFT stage was rewritten.

**Actual test results in artifact construction environment**:
17 `unittest` methods passed, including a 16-case synthetic subtest, an
entire 16-case isolated synthetic runner integration, complex quadratic
Hessian test, CP-phase single-scalar loop normalisation, full contributions,
channel reconstruction, zero-mass cancellations, EFT equal-input guards and
16 retrieved-source-matrix SVD checks. Because the user's Windows working tree
is not mounted, integration tests were run with a scratch transcription of
retrieved `core.py` and `comparison.py` that is **not** included in the ZIP.
Run `python -m unittest discover -s tests -p 'test_lo_direct_uv_*.py' -v`
in the real repo to validate against the local engine.

**Not yet executed on the user's live project**: newly generated per-channel,
per-heavy-fermion, per-real-scalar breakdown for the 16 frozen source-derived
component-input files. `scripts/RunLODirectUVCertification.py` performs that
calculation in the local repo without overwriting previous results.

**Blocking physical conventions**: (1) four-spinor `GammaCC` to Weyl phases,
(2) Dirac/Majorana mass and charge-conjugation normalisation from the native
free Lagrangian, (3) RGBeta-to-Matchete signed and normalised portal plus
lambdaT3 dictionary, (4) extra symmetry-allowed interaction numerical values
and native component CG when used, (5) robust complex-cancellation check for
near-zero numerical 3x3 masses.

**EFT comparison limitation**: UV frozen reference values are specified at
`mu=1e7 GeV`, while the two physical EFT thresholds are
`mu_F=1e5 GeV` (all heavy fermions together) and
`mu_S=1e3 GeV` (both heavy scalars together). The comparison gate never
pretends that an EFT result at another scale is the same renormalized point.
