# T3 self-bilinear normal-propagator matching continuation

## Project destination

This ZIP contains **new files only**; extract into `C:\thesistemp\T3neutrinomassmodelpipeline`, preserving paths:

- `Numerical/experimental/T3SelfBilinearMassBasis.py`
- `scripts/RunT3SelfBilinearNeutrino.py`
- `tests/test_t3_self_bilinear_mass_basis.py`
- `tests/test_t3_self_bilinear_neutrino.py`
- `docs/T3_SELF_BILINEAR_NORMAL_C5_BATCH.md`

No root README, prior benchmark, source matching artifact, existing result, or core production file is overwritten.

## Run (PowerShell, from repository root)

```powershell
python -X utf8 -m pytest -q `
    tests/test_t3_self_bilinear_mass_basis.py `
    tests/test_t3_self_bilinear_neutrino.py

python -m scripts.RunT3SelfBilinearNeutrino `
    --all `
    --run-id continuous_v1
```

To resume an interrupted run:

```powershell
python -m scripts.RunT3SelfBilinearNeutrino `
    --all `
    --run-id continuous_v1 `
    --resume
```

The runner requires the existing `continuous_v1` scalar trajectory and the previously completed `neutrino_v1` 16-model output. It also reads the already-saved Matchete extra-quartic tensor source `output/general_t3_stage22/matchete_tensors.json` and existing model-specific matching files. It **does not** run Wolfram, regenerate betas, rerun UV/EFT1 ODEs, or refit Yukawas.

## Physical computation performed

For a neutral real SU(2) representation of dimension 1 or 3, the declared quadratic potential is

\[
V_2=m_S^2 S^\dagger S+\mu_S^2 (S^\dagger)^T C (S^\dagger)+\text{h.c.}
\]

with **no** additional `1/2` on the holomorphic term. The real scalar mass multiplets are

\[
m_{R,I}^2=m_S^2\pm 2|\mu_S^2|.
\]

For the original one-loop topology with one **normal** propagator of that species,

\[
G_{SS^\dagger}(p^2)=\frac12\left(\frac{1}{p^2+m_R^2}+\frac{1}{p^2+m_I^2}\right),
\]

so the matching contribution from **that same topology** is evaluated as

\[
C_5^{\rm topology,normal}(\mu_S^2)=
\frac12 C_5^{\rm original}(m_R^2)+\frac12 C_5^{\rm original}(m_I^2).
\]

Here `C5_original` is the exact **existing model-specific hierarchical or legacy repository kernel** and retains the previously evolved UV Yukawa and heavy mass basis. The scalar-only EFT state is evaluated at the two **actual split** masses, not at their mean. The selected original kernel is a declared approximation if other vertices/anomalous propagators also contribute.

At the scalar endpoint the script computes the resulting 3×3 symmetric candidate `C5`. It then evolves SM+Weinberg running from 1 TeV to 100 GeV using the existing implementation, and calculates the three Takagi masses and oscillation observables. Every self-model case records the zero-self-mass control and the relative normal-channel difference.

## Additional quartic vertex transport

The script also uses the *saved native Matchete tensor components* for the 5 or 6 extra quartics in each of A(+0), A(-2), D(+0), D(-2), E(+0), E(-2). It rotates each supplied tensor factor into the real `R/I` scalar basis and writes every rotated ordered component tensor. A numerical polynomial identity is checked for each vertex.

**Important convention limitation:** This vertex transport applies the declared spherical-component intertwiner to Matchete tensor axes. An independent Matchete-to-spherical CG-component identification is **not certified**. Rotated vertices are generated but are **not** inserted into the neutrino-mass matching amplitude. They prepare the next actual vertex-inclusive matching implementation; they are not a physical C5 correction.

## Output files

All outputs are under:

`output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/self_bilinear_v1/`

- `A_p0_self_bilinear_c5.json`, `A_m2_self_bilinear_c5.json`, etc. (6): normal-channel resummed C5, low-energy masses, source hashes and limitations.
- `A_p0_rotated_extra_vertices.json`, etc. (6): exact algebraic R/I rotations of all 5/6 saved extra-quartic tensor components under the declared candidate basis.
- `summary.json`, `summary.csv`: six self-bilinear models.
- `combined_16_summary.json`, `combined_16_summary.csv`: ten existing candidates + six new restricted normal-channel calculations. These are **not** 16 certified full-general-scalar C5 predictions.

## Explicitly NOT implemented in this batch

- One-loop diagrams with **anomalous** `S-S` propagators and their required additional independent `HH` scalar quartic vertices.
- Certified relative CG normalization/signs between native Matchete and RGBeta tensors.
- Full fermion-first Wilson coefficient matching, intermediate Wilson running, and scalar threshold matching.
- Refit of neutrino oscillation parameters.

The numerical status is `CandidateSelfBilinearNormalPropagatorResummedC5AnomalousNotMatched`. Never represent this as the complete physical general-scalar neutrino mass.
