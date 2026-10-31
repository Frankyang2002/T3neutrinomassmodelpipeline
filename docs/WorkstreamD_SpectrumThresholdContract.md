# Workstream D: candidate scalar spectrum and threshold-handoff contract

## Scope and source status

`Numerical/experimental/GeneralT3SpectrumThresholdAdapter.py` is a **new additive adapter**. It imports, rather than replaces:

- `GeneralT3ScalarQuadraticMasses`: the 16-model trial quadratic manifest, analytic spectra, and Cartesian-real Hessians.
- `GeneralT3MixedScalarMassBasis`: the existing B(-1)/C(-1) complex 2×2 mass basis and `U` convention.

No existing runner, source matching file, symbolic Lagrangian, coupling map, or benchmark input is changed. This component does **not** perform EFT matching or evolve a state. It is a spectrum/threshold scheduling input to a future stage driver. The optional read-only `scripts/AuditGeneralT3MassCrossings.py` consumes existing `continuous_trajectory.json` reports and writes a **new** JSON audit, without altering the legacy trajectory or threshold sequence.

## Input contract

A (model inventory) supplies a **canonical model key** like `B_-1`, its declared SU(2) dimensions and allowed quadratic invariant, and the source manifest `docs/general_scalar_trial_operator_manifest.json`. All 16 keys are verified by `load_manifest`. A must supply an authoritative replacement map before candidate Hessians are described as physical mass matrices.

B (RGE) supplies a **single-stage continuous or sampled trajectory**. A `RunningMassSample` contains:

| Field | Units | Interpretation |
|---|---|---|
| `mu_gev` | GeV | Renormalization scale, positive; samples in strictly descending order |
| `m_s1_sq_gev2`, `m_s2_sq_gev2` | GeV² | Real scalar quadratic coefficients in the declared spherical/RGBeta basis |
| `mu2_gev2` | GeV² | Complex holomorphic bilinear, `0j` for models without an allowed operator |
| `fermion_eigenmasses_gev` | GeV | Optional *positive* heavy eigenmasses, consistently labelled across the same stage |

All remaining UV state components—SM gauge/Yukawa couplings, complex `y1/y2`, complex symmetric `MF`, scalar quartics, and extra per-model RGBeta quartics—remain stored in the existing `T3UVState`/`T3IntermediateScalarState` and `GeneralT3ScalarRGState` state vector, not duplicated here. The existing packing uses real/imag blocks. `mu2` lives in an **extra pair of real ODE entries** appended by `RunT3FullScalarContinuous`; its mass dimension is 2, unlike dimensionless quartics. All beta exports must carry the existing loop-factor convention `16*pi^2*dX/dln(mu)` and candidate basis status. **No unverified RGBeta↔Matchete identification is performed.**

The adapter's `scalar_spectrum(...)` returns distinct mass branches with their real degree-of-freedom multiplicities, the full Hessian and its real orthogonal eigensystem, and for B/C mixed models **the exact legacy mixed-doublet `U` object**. It checks eigenvalues against the legacy Hessian and analytic spectrum. The Hessian is *unbroken-phase quadratic only*, so it does not include EWSB shifts or pole-mass corrections. Degenerate real eigenvectors should not be interpreted as uniquely labelled SU(2) component fields.

`ordered_mass_events(spectrum, fermions)` orders positive masses at one scale. `find_mass_crossings(model, samples)` finds **mu = M(mu)**, re-diagonalizing scalar quadratic matrices during interpolation in `ln(mu)`. It reports roots **within the supplied stage only**; unbracketed roots are not invented. All results carry `matching_required_not_applied` status, and a model-specific EFT stage cannot yet be selected automatically.

## Handoff contract (future stage driver; NOT implemented by this adapter)

C (matching) must provide a callable tied to `(model, EFT stage, incoming active-field list, outgoing active-field list, threshold masses, field basis and U/CG maps, operator basis, scale, scheme)`. It returns a *complete* new Wilson tensor dictionary and a finite complex correction **with explicit operator normalization**; a missing coefficient must produce a hard blocked/missing status, **never** `0`. For the low-energy Weinberg operator, every stage must preserve the full complex symmetric 3×3 coefficient (`C5 = C5.T`, **not** Hermitian). Its physical conversion convention `m_nu = N * v² * C5` must be supplied and certified by C; D must not guess `N`. After charged-lepton basis rotation, the Takagi factorization obeys `U.T @ m_nu @ U = diag(m_i >= 0)` and is distinct from an ordinary Hermitian eigendecomposition. The existing SM+C5 running and neutrino observable modules remain the intended consumers.

D must not integrate beyond a threshold until its outgoing EFT's verified one-loop beta functions and matching boundary are available. The stage transition must report the difference `C_out(mu_t) - C_in(mu_t)` separately from numerical integration error; a finite jump is not by itself an error. Physical threshold continuity requires checking the **matched amplitude**, not asserting that Wilson coefficients are unchanged.

## Regression and limitations

Run `python -m unittest discover -s tests -p test_general_t3_spectrum_threshold_adapter.py -v` from the project root after copying the ZIP contents. The unit suite validates all 16 quadratic branch spectra against the existing Hessian, positive/tachyonic behavior, the B/C `U` and eigenvalues, non-hardcoded mass ordering, and running threshold roots. Synthetic mass values in the tests are **not benchmark UV inputs**.

To additionally verify saved B/C conditional reports (if they exist locally), run:

```powershell
python -m unittest discover -s tests -p test_general_t3_bc_conditional_regression.py -v
```

The regression checks existing report masses without overwriting or recomputing the conditional matching. If reports are not present, it skips with an explicit message; missing reports are **not** counted as reproduced results. The numerical values tested are the unchanged public main-branch saved conditional masses, with a 1e-10 eV absolute tolerance. Full end-to-end benchmark replay using the user's *current local files* remains a separate, unverified integration test.

**Not certified:** full physical RGBeta/Matchete map, physical EWSB scalar spectrum and component basis, finite Y/X one-loop matching, arbitrary-order EFT beta coverage, and three-family Wilson anomalous dimensions. None is replaced with a numerical guess.

## Run on saved all-model continuous-v1 trajectory

```powershell
python -m scripts.AuditGeneralT3MassCrossings --all
```

The audit reads the original UV/EFT1 **scalar** trajectories, locates only bracketed scalar mass crossings, and records the heavy fermion masses at the **saved grouped F boundary**. The continuous report does not store a sampled UV heavy-fermion Takagi trajectory: **fermion crossing scales cannot be recovered from those saved masses**. An uncrossed or tachyonic state is reported, not silently dropped. Output is a new `workstream_d_threshold_audit_v1/audit_summary.json` under the input run directory; the command refuses to overwrite it.
