# Rough UV -> EFT1 -> SMEFT -> neutrino mass benchmark

**Files:**

- `output/general_observable/lo_fermion_boundary/run_rough_smeft_neutrinos.py`
- `output/general_observable/lo_fermion_boundary/test_run_rough_smeft_neutrinos.py`
- `output/general_observable/lo_fermion_boundary/README_ROUGH_SMEFT_PILOT.md`

## What it actually does

The repository already contains a successful `16/16` approximate LLSS-to-Weinberg scalar matching run under
`output/general_observable/lo_joint_scalar_matching/exploratory_all16_v6_native_fix/`.
That calculation is produced by `scripts/RunT3JointScalarAll16V6.py` from the 16-model fermion-first boundary, renormalisable EFT1 running snapshots, the full original native HHSS quartic tensor, and real heavy scalar mass matrices. It integrates **both scalars together** at one scalar matching scale and records

`C5_LLSS_total_one_loop_per_N_gev_inverse = C5_EFT1_induced_leading_log_per_N_gev_inverse + C5_joint_scalar_finite_bubble_per_N_gev_inverse`.

This new command **consumes that calculated numerical scalar result** (does not rematch or double count either term), sets an explicitly *assumed* signed spinor normalization `N=+1` by default, then:

1. Reads `output/general_observable/all16_v1/{model}/running.json` for the SM couplings at the same `muS`.
2. Runs `Numerical.running.SMWeinbergEvolution.evolve_sm_weinberg` from `muS=1000 GeV` to `mu_low=91.1876 GeV`.
3. Computes `m_nu = -(v^2/2) C5`, uses `physics.NeutrinoTrajectory.charged_lepton_mass_basis_matrix`, and obtains Takagi neutrino masses, mass splittings and `|PMNS|` using `physics.NeutrinoObservables`.
4. Writes explicit complex `3x3` matrices and observables into new local JSON files and a CSV.

This is a **rough numerical implementation**, not a final UV prediction. In particular, `N=+1`, the native Matchete/RGBeta quartic mapping, benchmark scalar basis alignment, and the limited EFT1 leading-log truncation are modeling choices, **not proven physical inputs**. Some supplied benchmark Yukawas are diagonal and therefore do not produce a realistic PMNS matrix. No fitting, new Matchete evaluation or full RGE/threshold campaign is attempted. Additional hard-fermion one-loop matching and full EFT1 Wilson evolution are not included.

## PowerShell

From the repository root:

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$dir = 'output\general_observable\lo_fermion_boundary'
$zip = "$HOME\Downloads\T3_Rough_SMEFT_Neutrinos.zip"
if (-not (Test-Path -LiteralPath $zip)) { throw "ZIP not found: $zip" }
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

python -X utf8 "$dir\test_run_rough_smeft_neutrinos.py"
if ($LASTEXITCODE -ne 0) { throw 'Smoke tests failed' }

python -X utf8 "$dir\run_rough_smeft_neutrinos.py"
if ($LASTEXITCODE -ne 0) { throw 'Rough SMEFT/low-energy calculation failed' }

Import-Csv "$dir\rough_smeft_pilot_v1\rough_neutrino_masses.csv" | Format-Table -AutoSize
```

Default models: `A_-4`, `B_-1`, `C_-1`.

To run all sixteen models into a **fresh** directory:

```powershell
python -X utf8 "$dir\run_rough_smeft_neutrinos.py" --all16 `
  --output-dir "$dir\rough_smeft_all16_v1"
```

To select an alternative preliminary signed normalization or low scale:

```powershell
python -X utf8 "$dir\run_rough_smeft_neutrinos.py" --model B_-1 `
  --normalization -1 --mu-low 91.1876 `
  --output-dir "$dir\rough_smeft_B_m1_negativeN_v1"
```

The default output paths are:

- `output/general_observable/lo_fermion_boundary/rough_smeft_pilot_v1/A_m4.json`
- `output/general_observable/lo_fermion_boundary/rough_smeft_pilot_v1/B_m1.json`
- `output/general_observable/lo_fermion_boundary/rough_smeft_pilot_v1/C_m1.json`
- `output/general_observable/lo_fermion_boundary/rough_smeft_pilot_v1/rough_neutrino_masses.csv`
- `output/general_observable/lo_fermion_boundary/rough_smeft_pilot_v1/summary.json`

The script never writes into the existing 16-model fermion release or the existing joint-scalar outputs. To repeat a run, select a fresh `--output-dir`.

## If the existing scalar C5 inputs are absent

The source solver is already implemented; run it directly **once** to regenerate 16 benchmark scalar outputs into a new folder in its permitted area:

```powershell
python -X utf8 -m scripts.RunT3JointScalarAll16V6 `
  --allow-exploratory-map --allow-exploratory-special-bases `
  --output-dir 'output/general_observable/lo_joint_scalar_matching/quick_recomputed_v6'

python -X utf8 "$dir\run_rough_smeft_neutrinos.py" `
  --source-dir 'output/general_observable/lo_joint_scalar_matching/quick_recomputed_v6' `
  --output-dir "$dir\rough_smeft_recomputed_v1"
```

The scalar solver's output must be numerically successful for the models selected. Install `numpy`, `scipy`, and `sympy` in the active `.venv` if absent. Wolfram is not required to consume the existing cached numeric results.
