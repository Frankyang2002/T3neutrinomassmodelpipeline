# Independent full-UV leading one-loop neutrino mass

## Scope and present certification status

This directory is **new**; it imports no EFT matching, EFT1 mixing, scalar threshold, or RGE module. `core.py` evaluates a **physical mass-eigenstate** Weyl-fermion/real-scalar loop exactly at external momentum `p=0`, for arbitrary complex symmetric Weyl mass matrices, general complex scalar quadratic forms, and any number of component vertices. It does **not** approximate scalar-to-fermion mass ratios.

The core is tested. The automatic conversion of the project's Matchete four-component UV Lagrangian into correctly phased **two-component** UV vertices is **not yet certified**. Accordingly, `BuildT3DirectOneLoopInputs.py` builds explicitly **uncertified/provisional** source-dependent component files; `RunT3DirectOneLoop.py` refuses to label them certified without an independent UV spinor, CG, scalar-potential-sign, and normalization validation. **Do not describe provisional output as a verified T3 prediction.** No 16-model physical mass matrices are supplied in this package.

The 16 frozen source benchmarks are in `Numerical/benchmarks/fitted_reference/T3_*.json`. They contain three heavy-flavour masses, two flavour Yukawa matrices and ordinary scalar couplings, but they do **not** fix every symmetry-allowed general UV coupling. A model-wide general-Lagrangian calculation therefore also needs the values of any additional Yukawas, quartics and scalar bilinears. The builder explicitly restricts itself to the frozen ordinary benchmark and rejects nonzero unsupported extra portal coefficients; the general kernel itself handles any additional vertices when supplied.

## Direct physical mass formula

In a left-handed two-component Weyl basis,

```
L_mass = -1/2 psi^T M_F psi + h.c.
L_Y    = - nu_alpha psi_i [Y_{alpha i s} phi_s + Z_{alpha i s} phi_s^*] + h.c.
V_2    = phi^dagger A phi + [phi^T B phi + h.c.]/2
```

`M_F` is complex symmetric; `A` is Hermitian; `B` is symmetric. Their Takagi and real-Hessian diagonalizations give

```
U_F.T M_F U_F = diag(M_k),        M_k >= 0
phi_s = sum_a C_{s a} R_a,       R_a canonically normalized real scalars
Gamma_{alpha k a} = sum_{i,s} U_F[i,k]
                      * (Y[alpha,i,s] C[s,a] + Z[alpha,i,s] conj(C[s,a]))
```

The complete leading broken-phase mass at zero external momentum is

```
(m_nu)_{alpha beta} = [1/(16*pi^2)] * sum_{k,a}
                      M_k Gamma_{alpha k a} Gamma_{beta k a}
                      * B0_MSbar(0;M_k^2,m_a^2;mu)
B0_MSbar = 1 - [M_k^2*ln(M_k^2/mu^2) - m_a^2*ln(m_a^2/mu^2)]/(M_k^2-m_a^2)
```

`B0` is evaluated using `log1p`/Taylor limits at coincident masses. The UV-pole coefficient is `sum_{k,a} M_k Gamma_{alpha k a} Gamma_{beta k a}`. Its cancellation is essential for a finite radiatively generated mass in the absence of an allowed counterterm; the script tests it. The result is a **complex symmetric** matrix, never made Hermitian. Takagi gives nonnegative singular values, with a separate residual. For the convention `L_Y=-h nu N phi`, `phi=(R+iI)/sqrt(2)`, the single-neutral-fermion check gives `h^2 M [B0_R-B0_I]/(32*pi^2)`: compare to formulas using a differently normalized complex-scalar Yukawa only after translating conventions.

Every output contains a matrix for each **vertex pair**, plus individual matrices for each heavy-fermion and real-scalar eigenstate, the complex total, three neutrino masses (eV), UV-pole diagnostic, Takagi residual, exact mass hierarchy ratios and a cancellation factor. Individual contributions depend on the subtraction convention; the finite sum is physical after cancellations.

## Physical input JSON contract

One JSON object per model under the **new** project path `input/lo_direct_one_loop/A_-4.json`, etc. All arrays are complex-encoded as `{"real": [...], "imag": [...]}`. The main required keys are:

- `model`: `A_-4`, `B_-1`, etc.
- `fermion_mass_gev`: `nF x nF` symmetric Weyl mass matrix; include both chiralities of vectorlike fermions, and every SU(2) component.
- `scalar_A_gev2`, `scalar_B_gev2`: complex scalar quadratic coefficients for **all** relevant charge components. Tachyons and massless scalar IR cases are rejected.
- `vertices`: list of independent source-labeled entries containing `field` and/or `conjugate`, each shaped `3 x nF x nS`. All *allowed* source interactions must be included, even if suppressed. Do **not** collapse independent Yukawa structures.
- `provenance.benchmark_sha256`: SHA-256 of the frozen input file, checked automatically.
- `provenance.validated_uv_component_conventions`: explicit certification flag; the provisional builder sets it to `false`.

This input format allows all 16 representations and their charged/neutral internal states, but it is **not** proof that a particular Matchete-to-Weyl input mapping is correct. A failed UV-pole cancellation also blocks that case.

## Source-checked model catalog

The repository's `Lagrangian/model/T3ModelCatalog.wl` defines A=(1,3,2), B=(2,2,1), C=(2,2,3), D=(3,1,2), E=(3,3,2), with the 16 alpha assignments in the runner. The Yukawa operator structures verified from `Lagrangian/interactions/T3Topology.wl` are `Bar[L] F^c S1` and `Bar[L] F S2^dagger`; its mixing operator is `H H S1 S2^dagger + h.c.`. The repository has an **experimental** general-scalar/general-fermion sector. Source CG arrays must be used, not guessed from the representation dimension or a proposed multiplicity factor.

## Install and run (PowerShell, from repository root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_direct_one_loop_partial.zip" -DestinationPath . -Force
python -m pip install numpy sympy
python -m unittest discover -s tests/python -p "test_t3_direct_one_loop*.py" -v
```

If the local source has already generated `Reports/output/full/multiplicity/component_support_audit.json`, `matchete_generator_analysis.json` and `cg_covariance_audit.json`, generate *provisional* inputs without changing existing outputs:

```powershell
python scripts/BuildT3DirectOneLoopInputs.py --project-root . `
  --audit-dir Reports/output/full/multiplicity `
  --output input/lo_direct_one_loop
python scripts/RunT3DirectOneLoop.py --project-root . `
  --physical-input input/lo_direct_one_loop `
  --output output/general_observable/lo_direct_one_loop `
  --allow-uncertified
```

`--allow-uncertified` is **diagnostics only**, not verification. Omit that option when the component UV convention is independently certified and `validated_uv_component_conventions` is sourced from that verification. The all-16 command exits with code `2` when any benchmark is incomplete; it still writes `all_16_status.json` and the successful individual model outputs. It does **not** alter prior numerical outputs.

If any of the three source CG audit files are missing, the provisional component builder stops rather than reconstructing or overwriting the repository's original audit files. Supply copies of the completed, source-verified audit artifacts in a dedicated new directory and pass that directory with `--audit-dir`. The registered CG exports and SU(2) generator probes are prerequisites. No existing EFT or audit output needs to be regenerated or modified by this package.

## EFT comparison

An optional `--eft-reference-dir <directory>` accepts `A_-4.json` etc with the **physical** `3x3` neutrino matrix in GeV, `physical_dimension="mass_matrix_gev"`, `benchmark_sha256`, `renormalization_scale_gev`, `scheme` and `same_renormalized_uv_parameters=true`. Without the exact same renormalized input and scale, the comparison is **skipped**, not misreported. Existing EFT results at `mu_S`/the low scale cannot be compared directly to the frozen **UV `base_state`** at `1e7 GeV` without consistent running. Differences can be of order `m_S^2/M_F^2` (and logs), and exact equality outside the expansion's validity is not a requirement.

## Coverage and current limitations

- Verified in synthetic tests: exact loop function including degenerate masses; scalar-mixing and LNV zero limits; arbitrary complex symmetric/degenerate fermion Takagi; Dirac-fermion paired-chirality loop; complex scalar real Hessian; contribution bookkeeping; EFT same-point comparison guards.
- Verified from live GitHub source: the 16-model catalog and frozen fitted-reference input format; UV two-Yukawa and mixing vertices; the presence of additional general-scalar and general-fermion interaction options.
- **Not yet verified:** specific source component phase/normalization mapping, physical CG multiplicities, scalar-source-vs-quadratic potential coefficient signs, full source component data for optional interactions. The source adapter is explicitly provisional on these grounds. This package cannot truthfully report 16 source-certified physical matrices until those checks are completed.
