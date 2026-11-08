# A(-4): conditional full Wilson evolution after the verified fresh-index RGBeta fix

## Scope

The fixed native RGBeta scalar/Yukawa/gauge running background is reused from
`output/a_minus4_next_batch/fresh_index_production/`. The F-only tree-level
fermion threshold is taken from the **explicitly declared**, but not native
spinor-normalized, master-Weyl convention in
`output/general_observable/a_minus4_boundary_batch1/declared_convention_boundary.json`.

Numerical field basis: **45 light left-Weyl fermions**, **12 real scalar
components** (H, S1, S2). The new code maps the declared local six leptons
to the matching 45-fermion slots by field name. It never maps a lepton local
index 1 to the first quark component.

Full-sum normalization uses exactly

\[\mathcal L_5=\sum_{ijab}C_{ijab}(\psi_i^T\epsilon\psi_j)\phi_a\phi_b+\mathrm{h.c.}\]

without inserting a `1/4` or `1/2`. Only the two separate pair exchange
symmetries are applied. There are
\((45\cdot46/2)(12\cdot13/2)=80,730\)
independent **complex** coordinates. The saved trajectory contains all of them,
including coordinates that vanish at the initial threshold.

## Formula and authority

Use the existing `GeneralEFT1SevenTermFlow.evolve` implementation in mode
`fonseca_eq22_six_term`, which omits the disputed crossed-Yukawa channel and
retains the candidate transposed gauge-Casimir correction. Its gaugeless
limit has been cross-checked against the independent Eq. (22), and the gauge
representation has undergone prior Ward audits. The **complete physical
gauge/Yukawa one-loop Wilson RGE remains unverified**: numerical closure is
not equivalent to loop-level certification.

The conditional branches are `rho=+1` and `rho=-1` for a fixed set of native
yukawa inputs. Under **global sign reversal**, the homogeneous real-linear
Wilson ODE satisfies `beta(-C)=-beta(C)` and therefore `C_-(mu)=-C_+(mu)`.
Do not extrapolate this statement to general complex phases: `beta(C)`
contains terms in `C*`, and a generic complex prefactor need not commute with
the evolution.

The code does **not** identify the native Matchete `GammaCC` matrix/sign, certify
the physical Matchete-to-RGBeta quartic dictionary, include the hard one-loop
F threshold, integrate out the two heavy scalars, produce the matched SMEFT
Weinberg coefficient, or calculate neutrino masses.

## Files installed (full files, no diffs)

- `Numerical/running/Aminus4ConditionalFullWilson.py`
- `scripts/RunAminus4ConditionalFullWilson.py`
- `tests/test_a_minus4_conditional_full_wilson.py`
- `docs/A_MINUS4_CONDITIONAL_FULL_WILSON_EFT1.md`

No previous production files are replaced.

## PowerShell from the project root

```powershell
Set-Location "C:\thesistemp\T3neutrinomassmodelpipeline"
Expand-Archive `
    -Path "C:\Users\popok\Downloads\Aminus4_ConditionalFullWilson_EFT1_Benchmark.zip" `
    -DestinationPath "C:\thesistemp\T3neutrinomassmodelpipeline" `
    -Force
python -m unittest discover -s tests -p "test_a_minus4_conditional_full_wilson.py" -v
python scripts\RunAminus4ConditionalFullWilson.py `
    --output-root "output\a_minus4_next_batch\conditional_full_wilson_v1"
```

If the first integration is too slow for an exploratory benchmark, do **not**
truncate the 45-Weyl basis. The default uses RK23 at `rtol=5e-5` with 3
reported output scales, but it still integrates across every corrected
coupling interpolation knot. A more accurate reference run (in a NEW
output folder) can use `--method RK45 --rtol 1e-7 --points 7`.

## Output and physical meaning

`output/a_minus4_next_batch/conditional_full_wilson_v1/summary.json` records
input hashes, signed branch relation, global gauge Ward residuals, solver
statistics, full-tensor change, and the directly generated LLHH projection.

`output/a_minus4_next_batch/conditional_full_wilson_v1/full_conditional_C_trajectory.npz`
contains: `indices_1based`, `scales_gev`, the complete independent symmetric
trajectory, the complete 45x45x12x12 boundary/final rho=+1 tensor, and the
negative rho=-1 endpoint.

**The directly generated LLHH is NOT the matched SMEFT C5.** To complete the
EFT sequence, integrate out *all* scalar mass eigenstates simultaneously,
including both hard matching and running-to-threshold effects, then match and
run SMEFT. There is no numerical neutrino-mass claim in this deliverable.
