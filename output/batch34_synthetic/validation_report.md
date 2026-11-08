# T3 A(-4) Batches 3–4: seven-term Wilson-tensor flow — validation

**Status:** synthetic numerical verification only. **No physical SMEFT C5 certification.**

## Revision and source audit

- GitHub `main` HEAD inspected read-only: `12128de602ac6bf109dbad50dbac8446ed629434` (identical to previously inspected SHA).
- The numerical `RGE/running/intermediate/GeneralEFT1C5Mixing.py::calculate_master_beta` and the symbolic `RGE/general/WilsonTensorRGE.py` implement the same seven structures. They share a crossed-Yukawa **two simultaneous exchange** choice; the four independent exchanges denoted by `sigma({a,b} x {i,j})` in Aebischer–Bresciani–Selimović arXiv:2512.16890v1 Eq. (4.85) differ. New evaluator implements four and offers `corrected_crossed=False` *for audit only*; the existing upstream modules remain unchanged. Other six numerical contractions agree with the expressions in the inspected symbolic source, with synthetic component-loop tests.
- The old `GeneralEFT1C5Mixing.solve_boundary` is **fixed-order one-loop logarithmic mixing with frozen inputs**, not the complete coupled Wilson ODE. New ODE uses `16*pi^2 dC/dln(mu)=sum beta_n(C,C*,g(mu),y(mu),lambda(mu))`.
- Source paths inspected: `RGE/general/{WilsonTensorRGE,AnomalousDimensions,FermionBasis,ScalarBasis,GaugeGenerators}.py`, `RGE/running/intermediate/GeneralEFT1C5Mixing.py`, `Numerical/running/{GeneralT3RealLinearWilsonFlow,GeneralT3WilsonEvolution,GeneralT3OrderedTensor}.py`, `RGE/matching/WeinbergTensorAdapter.py`, `scripts/RunGeneralEFT1C5Mixing.py`, and related unittest modules. No README/status documentation was used.

## Operator basis and normalization

`L5=C_ijab psi_i psi_j phi_a phi_b + h.c.` (no factorial); `V4=lambda_abcd phi_a phi_b phi_c phi_d / 4!`. Each fermion pair `i<->j` and scalar pair `a<->b` is separately symmetric. Scalars use interleaved `Re, Im` components; Higgs components `[Re H+, Im H+, Re H0, Im H0]`. The full three-generation SM Weyl inventory contains `n_f=45`, plus two heavy scalar multiplets with `n_s=12` real scalar coordinates for A(-4). Unrestricted pair-symmetric envelope size `(45*46/2)*(12*13/2)=80,730` complex coordinates. It includes all initially zero positions rather than treating the initial LLSS support as an invariant subspace. These are tensor envelope coordinates; **not** a claim that 80,730 are all gauge-independent physical operators.

The induced LLHH 3x3 projection uses `4 C[nu_i,nu_j,Re H0,Re H0]`, under the boundary convention `L5^W=(1/2)C5_ij(L_i.H)(L_j.H)+h.c.`. This is **not** the fully matched SMEFT `C5`: the simultaneous scalar threshold matching is not included.

## Synthetic numerical runs (generated data)

All inputs below are explicitly synthetic, including gauge/Yukawa/quartic values and canonical boundaries; they are not gauge-invariant A(-4) EFT1 physical benchmarks. Thresholds are 100,000 GeV → 2,000 GeV. Only the dimension and tensor architecture of the large run mimic the A(-4) model.

| Quantity | Small source-decomposed fixture | Full-size scaling fixture |
| --- | ---: | ---: |
| Fermion, real-scalar components | 6, 12 | 45, 12 |
| Independent separately symmetric coordinates | 1,638 | 80,730 |
| Initially nonzero coordinates | 8 | 3 |
| Newly nonzero coordinates at M_S | 161 | 36 |
| Coupling knots | 3 | 2 |
| RHS evaluations | 91 (2 segments) | 44 (1 segment) |
| Max per-sample beta pair-symmetry residual | 3.01e-16 | 1.06e-17 |
| Conjugate-term full beta norm at M_F | 1.746e-12 | 2.029e-12 |
| Conjugate-term full beta norm at M_S | 1.196e-12 | 1.349e-12 |
| Relative \|C_ODE(M_S)-C_fixed_order(M_S)\| / \|C(M_F)\| | 6.610e-4 | 1.119e-3 |
| Source-channel recomposition relative residual | 1.287e-16 | not calculated |

A separate small-fixture tolerance-convergence check at `rtol = 1e-6, atol = 1e-18` (61 RHS evaluations) versus `rtol = 1e-9, atol = 1e-21` (85 RHS evaluations) gave a relative endpoint difference of **1.217e-16**.

The complete symmetric coordinate trajectory, initial/final fully expanded complex tensor and seven separate full beta tensors at both thresholds are in each `seven_term_trajectory.npz`. JSON `seven_term_summary.json` contains scales, seven per-term 3x3 LLHH beta projections at each output time, source decomposition, all 7 Frobenius norms, anti-linear term norm, zero-start generated coordinates, normalization and solver metadata.

The leading **direct** LLSS→LLHH projection at the high threshold in these fixtures comes solely from `scalar_pair`; each of the other six projected contributions is zero there (verified numerically). After evolving the **full** C tensor, the other terms can contribute to lower-threshold LLHH derivatives via the generated coordinates. This is an iteration of a one-loop anomalous dimension, i.e., higher-order logarithmic resummation; it is **not** an extra fixed-order one-loop matching source.

## Tests and limits

Local `python -m unittest discover -s tests -p test_general_eft1_seven_term_flow.py -v` outcome: **16 tests, 14 passed, 2 skipped**. Passed: seven independent component-sum contractions; crossed term discrepancy; pair symmetries; anti-linearity and real-linearity; real orthogonal coordinate covariance; gauge-only SU2/U1 Weinberg Ward covariance; full-size coordinate inventory; generated non-lepton coordinates; all-scale source superposition; fixed-order infinitesimal-log limit; threshold direction; RK convergence; missing input/unsymmetric boundary rejection; scalar-pair direct LLHH projection. Skipped: symbolic upstream import and installed repository fixed-order integration because the exact repository source is not mounted in this execution container. **Run the same tests under the Windows project root to exercise those two integrations.**

Coupling interpolation in log(mu) is a numerical handoff convention, not a replacement for the renormalizable Chat B RGEs. The solver integrates separately across sample knots. It does not generate RGEs for the couplings itself. For a certified run Chat B must supply: full 45-Weyl, 12-real-scalar basis maps; complete running 3x3-family Yukawa, scalar quartic, and `SU3C`, `SU2L`, `U1Y` gauge couplings; all 8+3+1 constant representation generators, consistent basis/order; scale samples spanning both thresholds; explicit verification flag. Chat A must provide the certified `t3-eft1-certified-master-weyl-v1` boundary with canonical complex-flavour components and its full certificate dictionary. The engine retains exact zeros outside the source's canonical components.

### Still unverified for physical A(-4)

1. Chat A's external Matchete→Weyl normalization, physical operators and their certified coefficients.
2. Chat B's fully verified general EFT1 coupling tensors and their own one-loop running (including quark/gauge sectors), indexed in the *same* basis.
3. Full gauge-covariant A(-4) direct-projection zero checks using **physical** Yukawas/quartics and a Ward-closed matched boundary. The current zero checks are only synthetic.
4. Absolute normalization/UV log partition, finite hard-scalar simultaneous matching, and physical SMEFT C5. These are outside Batches 3–4.
5. Full integration with the installed repository modules on the user's Windows checkout: two tests are skipped locally by design.

**Installation scope:** Only additive files inside `RGE/running/intermediate`, `scripts`, `tests` and `output/batch34_synthetic`; no existing Chat A/B code and no GitHub source is altered.
