# EFT1 Stage A: source mapping and first numerical master-tensor transport

**Read-only repository inspection:** `Frankyang2002/T3neutrinomassmodelpipeline`, branch `main`, commit `ea93a730812cc6a288ae7eba363874a2166ea6c6` (9 November 2026). **Target:** user's local Windows repository. No GitHub writes. All files in this ZIP are additive.

## Verified pre-existing source facts

- `Numerical/running/GeneralT3DynamicFermionWilson.py::_dense_operator` gives `source.raw[p,q,i,j,a,b]` with shape `(3,3,2,2,n,n)`, real **heavy-only** scalar axes and exact **combined** lepton and scalar exchange symmetry. It explicitly flags `absolute_matchete_to_weyl_normalization_certified: False`.
- `Numerical/running/GeneralT3OrderedTensor.py::scalar_basis` places four Higgs real components at indices 1..4; heavy scalar real components start at 5. For A(-4), `n=8` (complex singlet + complex triplet). **No implicit scalar/Higgs exchange factors can be inserted.**
- `Numerical/joint_scalar_matching/fermion_boundary_bridge.py::write_provisional_source` writes `C_LLSS_formal_muF` and sidecar JSON, with GammaCC normalization intentionally unresolved.
- `Numerical/running/Aminus4EFT1Couplings.py` defines the 45-left-Weyl/12-real-scalar RGBeta reference basis, full SM Yukawa and gauge generators, quartic tensor, and native EFT1 trajectory snapshots. It explicitly flags the **Matchete-to-RGBeta scalar/invariant basis map as uncertified**.
- `RGE/running/intermediate/GeneralEFT1C5Mixing.py::calculate_master_beta` implements the seven self-class contributions to `16 pi^2 dC_{ijab}/dln(mu)` (the `psi^2 phi^2 <- psi^2 phi^2` part of Eq. (4.85), arXiv:2512.16890). Eqs. (4.86)-(4.91) also contain mixings from *other* dimension-five classes; those have **not** been ruled out or included by this pilot.
- `scripts/RunGeneralDMinus2RealEFT1RGE.py` provides a physical real-basis **one-generation LLHH beta** from native D(-2) seed, but explicitly does not perform full-flavour Wilson transport, threshold matching or running renormalizable coupling certification.

## Added functionality

1. `Numerical/running/T3EFT1MasterSourceBridge.py`: reads the **actual** formal F-threshold export, checks original axis orbit symmetries, implements the explicit `L[p,i]` to 45-Weyl-index embedding, transforms the heavy real-scalar coordinates with an explicit real orthogonal matrix if supplied, and applies a supplied signed GammaCC factor **without pretending it has been proven**. No guessed 1/2 or 1/4 coefficients are inserted.
2. `Numerical/running/T3EFT1WilsonODETransport.py`: generic complex DOP853 solver for all `C_ijab` with a scale-dependent seven-term master beta, plus an A(-4) adapter consuming the **running** 45-Weyl/12-real-scalar `A_m4_eft1_tensors.npz` and `A_m4_eft1_tensor_basis.json` snapshots. Within `psi^2 phi^2`, both LLSS self-running and induced LLHH are advanced simultaneously. This is provisional, not a physical EFT1 prediction.
3. `scripts/ExportAminus4FermionMasterSource.py`: obtains the A(-4) source using the existing native full-flavour F-matching function at a caller-supplied `mu_F`, and writes the formal NPZ plus audit.
4. `scripts/RunAminus4EFT1WilsonODEPilot.py`: CLI for a **provisional** full master ODE run, requiring the user to explicitly request an identity basis-map diagnostic if no external scalar map is provided. Refuses to overwrite output directories.
5. Two tests: coordinate embedding/rotations/checks and synthetic, nonphysical ODE transport.

## Field and coefficient coordinates

Let `T[p,q,i,j,a,b]` be the formal output of F matching, with `p,q=0..2`, `i,j=0..1`, and `a,b=0..n-1` labelling only heavy real scalar components. The target fermion basis explicitly names `L[p+1,i+1]`; Higgs coordinates are `[Re H+, Im H+, Re H0, Im H0]`.

For the real-scalar change of basis `phi_source[a] = R[a,A] phi_target[A]`, the candidate embedded tensor is

\[
C^{\mathrm{target}}_{L_{pi} L_{qj},\,4+A,\,4+B}
  = k_{\Gamma\mathrm{CC}\to\mathrm{Weyl}}
    \sum_{a,b} T_{pqijab}R_{aA}R_{bB}.
\]

This expression defines a coordinate transport. The **signed factor** `k` and the **Matchete-to-RGBeta map** `R` are *not* established by the files inspected. Therefore, a formal run with `k=1`, `R=I` is explicitly **not a physically normalized Wilson evolution**. The mapping of real scalar coordinates must also respect the SU(2), U(1) generators and the interaction tensors, not merely orthogonality.

The code's status is always `ProvisionalFullMasterTensorODE_NotPhysicalPrediction` even if the user supplies candidate matrices.

## Verified locally in this deliverable

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
python -X utf8 -m unittest discover -s tests -p 'test_t3_eft1_*.py' -v
```

Seven **synthetic implementation** tests passed in the isolated execution environment; `compileall` passed. Neither a native full-flavour F-tree numerical boundary nor an A(-4) RGBeta coupling trajectory was available in the uploaded materials, so **no physical benchmark was executed** here.

## How to assemble the two physical-input files locally

The actual F-matching export requires an existing JSON with `MF,y1,y2` at the **fermion threshold**, and the repository's Matchete source files under the existing `F_ROOTS` search paths. Example:

```powershell
$uv = 'PATH_TO_A_m4_MF_y1_y2_AT_muF.json'
$muF = 100000.0  # REPLACE with your actual benchmark muF
python -X utf8 -m scripts.ExportAminus4FermionMasterSource `
  --uv-state $uv --mu-f $muF `
  --output 'output\eft1_stageA\A_m4_formal_F_source.npz'
```

For the 17-component renormalizable EFT1 state and beta export, reuse the repository's existing exporter, rather than inventing RGEs:

```powershell
$eft1 = 'PATH_TO_A_m4_running_couplings_muF_muS.json'
$beta = 'PATH_TO_A_m4_native_RGBeta_17_beta_export.json'
$muS = 1000.0  # REPLACE with the actual common scalar threshold
python -X utf8 -m scripts.ExportAminus4EFT1Couplings `
  --initial-file $eft1 --beta-file $beta --mu-s $muS `
  --samples 7 --output-dir 'output\eft1_stageA\A_m4_running'
```

Then run a **formal-only numerical diagnostic** (do not use its `C5` as a physical neutrino result):

```powershell
python -X utf8 -m scripts.RunAminus4EFT1WilsonODEPilot `
  --source 'output\eft1_stageA\A_m4_formal_F_source.npz' `
  --couplings 'output\eft1_stageA\A_m4_running\A_m4_eft1_tensors.npz' `
  --formal-identity-map `
  --output 'output\eft1_stageA\A_m4_provisional_ode'
```

With the actual scales and state files, the script either produces a clearly provisional `pilot_summary.json` and `formal_master_C_at_muS.npz`, or fails with explicit basis/schema/solver diagnostics. **These commands were not run against the user's private local source/results.**

## Remaining gates before calling this *physical* EFT1 running

- Derive and independently verify the Matchete PL/GammaCC to master two-component Weyl **signed** normalization from the exact F-matching Lagrangian.
- Derive the exact transformation between Matchete scalar CG/component coordinates and the full RGBeta renormalizable quartic/gauge/Yukawa basis, including complex phases and charge-conjugation conventions.
- Prove or enlarge dimension-five operator-class closure beyond the current `psi^2 phi^2 <- psi^2 phi^2` seven-term subsystem (see the other sources in Eqs. (4.86)-(4.91)).
- Reconcile the explicit EFT1-generated logarithm with the old fixed-order scalar bubble: subtract the shared infrared log once; retain the valid finite common-scalar threshold contribution, not an extra duplicated log.
- Add D(-2) complete running renormalizable-coupling tensor in the correct real-field basis; test both physical benchmark neutrino matrices only after the above gates.

The local preliminary `GeneralEFT1WilsonODE.py` and `TestEFT1WilsonODE.py` mentioned in the request do **not** appear on the inspected public branch and were **not uploaded**. Accordingly these additions do not overwrite or assume their local interfaces.
