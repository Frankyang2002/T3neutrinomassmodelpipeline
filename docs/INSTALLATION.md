# Installation and reproducible execution

This guide describes the **Python software baseline** and the separate requirements for live Wolfram/Matchete/RGBeta calculations. It is based on read-only inspection of the source tree at Git commit `806cdf08608dbe3f5e0dcd8af785cfad16f26af6` (8 November 2026). A clean installation and a fresh all-model symbolic rebuild have **not** been verified from this document alone.

## 1. Python environment (Windows PowerShell)

Python **3.12** is a practical starting point because the project's last reported passing core tests ran in a Windows Python 3.12 environment. Source syntax uses Python 3.10+ features, but Python 3.10/3.11 are **not** independently tested here.

From the repository root:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
if (-not (Test-Path .\.venv\Scripts\python.exe)) {
    py -3.12 -m venv .venv
}
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip check
```

If `.venv` already exists and is working, **do not recreate it**. Activate it, then run the `pip` commands. If activation is blocked by local PowerShell policy, invoke `.\.venv\Scripts\python.exe` directly in place of `python` rather than changing machine-wide execution policy.

`requirements.txt` lists directly observed Python imports: `numpy`, `scipy`, `sympy`, `pytest`, `matplotlib`, and `mpmath`. NumPy/SciPy cover numerical evolution and fitting; SymPy covers symbolic/group-theory adapters and tests; pytest runs the regression suite; Matplotlib and mpmath support additional plots and research numerics. The file deliberately has **no version pins**: it is a bootstrap dependency list, **not** a proven lockfile for every historical script or a guarantee of cross-platform reproduction.

## 2. Confirm the basic entry points

These commands do **not** perform live Matchete threshold matching:

```powershell
python -X utf8 pipeline.py --help
python -X utf8 full16_symbolic.py --help
python -X utf8 -m scripts.RunT3GeneralSymbolicStudy --self-test
.\scripts\RunCoreTests.ps1
```

The core test selection is defined in `scripts/RunCoreTests.ps1`. The most recent user-reported local execution passed **189 tests** on 8 November 2026. That is a recorded baseline, **not** a fresh verification of every future checkout.

For the broader retained test suite:

```powershell
.\scripts\RunCoreTests.ps1 -Full
```

The extended suite has previously reported missing generated artifacts, superseded test interfaces and scientific/numerical mismatches. A failure there is **not** evidence that the `pip` installation alone is wrong. See [testing strategy](TEST_STRATEGY.md) and [test-failure triage](TEST_FAILURE_TRIAGE_2026_11_08.md).

## 3. Live symbolic calculations: separate external dependencies

The following tools are **not installed by pip**:

| Tool | Where the source uses it | Required for |
| --- | --- | --- |
| Wolfram Language / `wolframscript` | `Lagrangian/WolframRunner.py` launches `wolframscript -file Lagrangian/RunModel.wl` | New symbolic model and threshold calculations |
| Matchete | `Lagrangian/RunModel.wl` evaluates ``Needs["Matchete`"]`` | Native invariant Lagrangian and matching |
| RGBeta | `RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl` evaluates ``Needs["RGBeta`"]`` | Fresh beta-function generation |
| LaTeX engine, where relevant | `Reports/` produces LaTeX/PDF materials | PDF compilation only, not basic Python tests |

Check that the Wolfram executable is discoverable:

```powershell
Get-Command wolframscript -ErrorAction SilentlyContinue
```

Successful executable discovery **does not** prove that Matchete or RGBeta loads, nor that the installed versions reproduce saved source conventions. Configure these packages in the local Wolfram environment according to their own installation instructions, then validate a small model before attempting the 16-model study. No unverified installation paths or package versions are prescribed here.

The general symbolic entry point is:

```powershell
python -X utf8 full16_symbolic.py --help
```

Running it without `--help` may start substantial Wolfram work. The source-level self-test above does not require live Matchete/RGBeta and does **not** certify their physics output.

## 4. Reproducing one numerical calculation

The numerical scan entry point is `Numerical/fitting/ScanCLI.py`. Its JSON configuration specifies paths to UV and intermediate-EFT RGBeta exports, the final Weinberg JSON, and an oscillation-fit target. It **does not generate those files** during the scan.

For an existing, compatible local configuration:

```powershell
python -X utf8 -m Numerical.fitting.ScanCLI path\to\existing_scan.json --dry-run
```

Replace the placeholder with a **real** existing configuration; do not assume a generic file is included. After reviewing its referenced paths, omit `--dry-run` to run the scan. A saved-input one-point smoke utility is also present:

```powershell
python -X utf8 scripts/maintenance/SmokeTestFullPipeline.py
```

That utility inspects existing configurations and writes under `output/audit/end_to_end_smoke/`. It may fail if generated JSON exports have been moved to an external archive. Restore or regenerate the *correct* source-matched artifacts rather than substituting unrelated JSON files. The last reported operational smoke test covered one ordinary T3-B(-1) benchmark; it was not a clean-room physics certification.

## 5. Capture the actual environment for reproducibility

When a local setup works, save **its installed versions** alongside the intended benchmark inputs and source commit:

```powershell
New-Item -ItemType Directory -Force output\maintenance | Out-Null
python --version | Out-File -Encoding utf8 output\maintenance\python_version.txt
python -m pip freeze | Out-File -Encoding utf8 output\maintenance\python_packages.txt
git rev-parse HEAD | Out-File -Encoding utf8 output\maintenance\source_commit.txt
```

Retain the selected JSON configuration, exported source files, source fingerprints, fitting seeds and output manifests. Avoid presenting a run from cached artifacts as a fresh Matchete/RGBeta regeneration. `output/` includes scientific inputs and is **not** safe to empty indiscriminately.

For scientific capabilities and certification boundaries, see [implementation status](IMPLEMENTATION_STATUS.md).
