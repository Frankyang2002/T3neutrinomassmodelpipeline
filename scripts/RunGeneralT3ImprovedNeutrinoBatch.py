"""Evaluate the 14 existing matched T3 C5 expressions using improved running outputs.

Run from the repository root:
 python -X utf8 -m scripts.RunGeneralT3ImprovedNeutrinoBatch

This is exploratory direct matching; it does not certify physical normalization
or higher-dimensional EFT1 Wilson-coefficient evolution.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from Numerical.matching.GeneralT3MultiModelNeutrinos import parser as numerical_parser, run

MODELS = (
    "A_-4", "A_-2", "A_+0", "B_-3", "B_+1", "C_-3", "C_+1",
    "D_-2", "D_+0", "D_+2", "E_-4", "E_-2", "E_+0", "E_+2",
)


def tag(model: str) -> str:
    return model.replace("-", "m").replace("+", "p")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model", action="append", choices=MODELS)
    p.add_argument("--running-root", type=Path,
                   default=Path("output/general_numerical_all14_improved/running"))
    p.add_argument("--matching-root", type=Path,
                   default=Path("output/general_physical_matching_batch"))
    p.add_argument("--output-root", type=Path,
                   default=Path("output/general_numerical_all14_improved/neutrinos"))
    p.add_argument("--mu-low", type=float, default=91.1876)
    p.add_argument("--quartic-real", type=float, default=0.02)
    p.add_argument("--quartic-imag", type=float, default=0.0)
    p.add_argument("--kappa-real", type=float)
    p.add_argument("--kappa-imag", type=float, default=0.0)
    p.add_argument("--strict-mapping", action="store_true")
    p.add_argument("--wolframscript", default="wolframscript")
    args = p.parse_args(argv)

    selected = tuple(dict.fromkeys(args.model)) if args.model else MODELS
    args.output_root.mkdir(parents=True, exist_ok=True)
    rows = []

    for model in selected:
        t = tag(model)
        running = args.running_root / f"{t}.json"
        coefficient = args.matching_root / t / "c5_coefficient.txt"
        output = args.output_root / t / "neutrinos.json"
        row = {"model": model, "running_json": str(running),
               "coefficient": str(coefficient), "output_json": str(output)}
        try:
            if not running.is_file():
                raise FileNotFoundError(f"Improved running missing: {running}")
            if not coefficient.is_file():
                raise FileNotFoundError(f"Symbolic C5 missing: {coefficient}")

            running_data = json.loads(running.read_text(encoding="utf-8-sig"))
            if running_data.get("model") != model:
                raise ValueError(f"Running model mismatch for {model}")
            if running_data.get("status") != "NumericalRenormalisableRunningOnly":
                raise ValueError("Running JSON is not a completed UV/EFT1 result")
            if not (running_data.get("uv_solver", {}).get("success") and
                    running_data.get("eft1_solver", {}).get("success")):
                raise ValueError("Running solver success flags are not both true")

            num = numerical_parser().parse_args([
                "--running-json", str(running),
                "--coefficient", str(coefficient),
                "--output-json", str(output),
                "--mu-low", str(args.mu_low),
                "--wolframscript", args.wolframscript,
                "--quartic-real", str(args.quartic_real),
                "--quartic-imag", str(args.quartic_imag),
            ])
            num.use_running_bilinear = True
            num.strict_mapping = args.strict_mapping
            num.kappa_real = args.kappa_real
            num.kappa_imag = args.kappa_imag

            result = run(num)
            status = result.get("low_energy_rge_status")
            row["low_energy_rge_status"] = status
            row["neutrino_masses_ev"] = result.get("neutrino_masses_ev")
            row["observables_scale_gev"] = result.get("observables_scale_gev")
            row["status"] = "Success" if status == "Success" else "LowEnergyRGEFailed"
            if status != "Success":
                row["error"] = result.get("low_energy_rge_errors", [])
        except Exception as exc:
            row["status"] = "Failed"
            row["error"] = f"{type(exc).__name__}: {exc}"

        rows.append(row)
        summary = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "scope": "exploratory_direct_C5_matching_and_SM_Weinberg_running",
            "physical_normalization_certified": False,
            "eft1_higher_dimensional_wilson_running_performed": False,
            "attempted": len(rows),
            "successful_to_mz": sum(r["status"] == "Success" for r in rows),
            "models": rows,
        }
        (args.output_root / "summary.json").write_text(
            json.dumps(summary, indent=2, allow_nan=False) + "\\n", encoding="utf-8"
        )
        print(f"{model}: {row['status']}", flush=True)
        if row["status"] != "Success":
            print(f"  {row.get('error', '')}", flush=True)

    success = sum(r["status"] == "Success" for r in rows)
    print(f"SM+Weinberg to MZ: {success}/{len(rows)} successful")
    print(f"Summary: {args.output_root / 'summary.json'}")
    return 0 if success == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
