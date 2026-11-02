# T3 fermion-first native hard-one-loop evidence batch (provisional)

**Install root:** `C:\thesistemp\T3neutrinomassmodelpipeline\`.

This batch runs the repository's *existing* read-only Matchete matching drivers in fresh Wolfram kernels, writing ALL outputs only beneath `output/general_observable/lo_certification_spinor/`.

- The 14 ordinary/nonmixed models use `Lagrangian/RunT3GeneralFThreshold.wl <outputdir> MODEL 1`. This driver calls `RunSequentialT3Matching[..., {{"F"}}, 5, 1, ...]`. It retains both scalar multiplets. Full-general experimental flags are `1`.
- Mixed `B_-1` and `C_-1` use `Lagrangian/RunGeneralT3FermionFirstEFT1OneLoop.wl <existing-mass-basis-input> <output.json>`. That driver rebuilds its own seven additional quartics and rejects the general-quartics opt-in. A separate kernel with experimental flags `0` is used. The exact numeric mass-basis input is checked, never synthesized.
- Optionally `--phase-probe` invokes `Lagrangian/RunT3YukawaSpinorPhaseNative.wl` in a separate kernel, preserving its dedicated flag requirements. It establishes abstract CConj identities, **not** the numerical charge-conjugation phase.

The Python runner validates the upstream 16/16 source-resolved F-tree release and SHA256 references *before* touching any output. It strictly prevents output outside this worker's directory. It never changes the immutable provisional boundary. It refuses to overwrite any prior run ID; use a fresh ID for each batch.

## PowerShell

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$env:PYTHONPATH = (Get-Location).Path
python -X utf8 -m unittest discover -s tests -p 'test_lo_spinor_*.py' -v
python -X utf8 scripts\RunLOSpinorHardThreshold.py --preflight --phase-probe
python -X utf8 scripts\RunLOSpinorHardThreshold.py --run-id native01 --phase-probe
```

The full native execution runs 16 fresh Matchete matching kernels plus one native abstract CConj kernel, and may require substantial computing resources. It writes:

```text
output/general_observable/lo_certification_spinor/
  native_hard_<source-summary-sha-prefix>_native01/
    A_m4/  ...  E_p2/
      wolfram.log
      inspection.json
      native_status.json   (ordinary models only)
      data/...             (ordinary Matchete outputs)
      <tag>_native_one_loop.json  (mixed models only)
      ...raw matching InputForm files...
    phase_probe/
      native_A_p0.json
      wolfram.log
      inspection.json
    summary.json
```

If a kernel fails, its exit code/log are retained and `summary.json` reports a **blocked** model. The command exits code `2` unless all requested models plus optional phase probe return readable raw evidence; it never replaces missing results with zeros. `--allow-partial` changes only the shell exit status, not the reported physics status. `--models A_-4 B_-1 C_-1` is available when investigating blockers, but the default runs all 16.

## Interpretation

Even if 17 kernels succeed:

1. Ordinary models have a saved, SHA-checked raw one-loop correction, **not** a renormalized LLHH Wilson coefficient. The parser does not guess hard LLHH operators by splitting nested InputForm strings.
2. B/C models have a saved raw hard LLHH candidate sector and number of candidate terms, **not** the renormalized scheme-matched coefficient. A zero candidate is **not** proof of zero hard matching.
3. Matchete's abstract `GammaCC`, `CConj`, and `Proj` identities do not alone assign the native `GammaCC` phase relative to the master Weyl epsilon.
4. No native connected UV four-point signed Dirac/Majorana+CG amplitude is independently evaluated by this batch. This remains a separate physical blocker.
5. No scalar matching or EFT1 running is performed. There is no intermediate EFT with a subset of heavy fermions or only one heavy scalar integrated out.

**Therefore:** `physical_master_wilson_boundary=null`, `native_GammaCC_to_master_Weyl_phase=null`, `hard_one_loop_MSbar_certified=false`, and `complete_leading_one_loop_neutrino_mass_certified=false` in every batch summary. Do not feed raw evidence to the physical EFT1/scalar workers. Preserve individual LLSS source terms from the existing spinor handoff.

Existing native-driver interfaces were inspected from read-only GitHub snapshot `2c29e5ba2c0350b7e9be292e02ceb8bc7c52c0d3`. Current local versions are checked by SHA256 in the preflight output, not assumed identical. Matchete is unavailable in the packaging environment, so actual native matching must execute in the Windows installation.
