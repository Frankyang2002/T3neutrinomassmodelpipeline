# T3 Batch 1 — A_-4 fermion-boundary normalization audit

## Current result

**Physical certification: BLOCKED.** The declared Weyl-basis identities establish
`C PL / epsilon = +1` and `C^-1 PL / epsilon = -1`, but the Matchete native
`GammaCC` phase and the signed master Wilson-operator prefactor are unmeasured.
Both candidate signs satisfy the same abstract conjugation identities. Neither
sign has been assigned to the physical boundary.

## Source-checked provisional data

- Native model: A_-4; SU(2) reps (S1,S2,F)=(1,3,2), Dirac F.
- V6 additional HHSS channels: 0.
- PL source terms: 1; independent LLSS coordinates: 60.
- Formal sample C[1,1,1,1,5,7] = 2.1535596891394695e-13 GeV^-1.
- Already-applied symmetric orbit factors are NOT reapplied. PR conjugates are NOT counted.
- This is a numeric **formal** tensor; no physical overall sign/prefactor is certified.

## Frozen source benchmark (UV values, not RG-evolved to mu_F)

- mu_UV=1e+07 GeV.
- mu_F=100000 GeV; mu_S=1000 GeV.
- Dirac F masses: [100000.0, 110000.0, 120000.0] GeV.
- Unbroken scalar mass parameters: [1142.857142857143, 1000.0] GeV.
- min(MF)/max(mS)=87.5; hence a perturbative fermion-heavy expansion is plausible.
- lambda_T3={'real': 0.1, 'imag': 0.0}. See JSON input for full Yukawa matrices and baseline quartics.
- Additional independent couplings fixed to zero in this **restricted benchmark at mu_UV**, not claimed RG-stable.
- Broken-phase scalar masses, physical source portal normalisation and running between mu_UV/mu_F remain uncertified.

## GitHub snapshot drift audit

Read-only main at `e581a434c4114ea1b100c50378c085529fb85be7`. Comparison to local checkout:
- `Lagrangian/model/T3ModelCatalog.wl`: SOURCE_DRIFT
- `output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/A_m4/boundary.json`: SOURCE_DRIFT
- `output/general_observable/lo_certification_couplings/ws2_extra_source_v1/A_m4/V6_source_evidence.json`: SOURCE_DRIFT

## Next missing amplitude and prescription

Compute one signed, connected native Matchete UV four-point amplitude for
`l l S1 S2` in the Dirac-doublet A_-4 model at zero external momentum, including
both exchange orientations, native `T3MixCG`, Grassmann ordering and the operator
prefactor. Project it independently onto a **declared** two-component Weyl
`epsilon_12=+1` amplitude. Do not choose `+1` from `C=i gamma2 gamma0` alone.
Then check whether already-Matchete-matched one-loop F terms contain hard/soft
subtraction before adding any EFT1 subtraction or counterterm; avoid double counting.

No C5 or m_nu calculation is claimed by this batch.
