"""Stage 12: extra-only B/C(-1) mass impact, phase scan and two-model CSV."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Numerical.experimental.GeneralT3NeutrinoImpact import load_and_analyze


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=Path("output/general_t3_numerical_pilot"))
    parser.add_argument("--output-root", type=Path, default=Path("output/general_t3_neutrino_impact"))
    parser.add_argument("--vev-gev", type=float, default=246.22)
    options = parser.parse_args()
    options.output_root.mkdir(parents=True, exist_ok=True)
    rows = []
    errors = []
    for model in ("B_m1", "C_m1"):
        source = options.input_root / f"{model}_matchete_hard_c5.json"
        target = options.output_root / f"{model}_extra_mass_impact.json"
        try:
            report = load_and_analyze(source, vev_gev=options.vev_gev)
            report["input_path"] = str(source)
            target.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
            data = report["scenarios"]["both_as_input"]
            vals = data["extra_only_takagi_singular_values_ev"]
            row = {
                "model": model,
                "extra_c5_frobenius_gev_inverse": data["c5"]["frobenius_norm"],
                "extra_mnu_frobenius_ev": data["extra_neutrino_mass_matrix_ev"]["frobenius_norm"],
                "extra_only_singular_value_1_ev": vals[0],
                "extra_only_singular_value_2_ev": vals[1],
                "extra_only_singular_value_3_ev": vals[2],
                "phase_scan_min_mnu_norm_ev": report["minimum_mass_norm_phase"]["mnu_norm_ev"],
                "phase_scan_max_mnu_norm_ev": report["maximum_mass_norm_phase"]["mnu_norm_ev"],
            }
            rows.append(row)
            print(f"{model} |extra m_nu|_F = {row['extra_mnu_frobenius_ev']:.6e} eV; "
                  f"singular values={vals}")
            print("Report:", target)
        except (ValueError, OSError, TypeError, KeyError) as exc:
            errors.append({"model": model, "reason": str(exc)})
            print(f"BLOCKED {model}: {exc}", file=sys.stderr)
    output = options.output_root / "extra_mass_impact_summary.json"
    result = {"status": "Success" if not errors else "Incomplete", "rows": rows, "errors": errors}
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if rows:
        csv_path = options.output_root / "extra_mass_impact_summary.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print("CSV:", csv_path)
    return 0 if not errors and len(rows) == 2 else 2


if __name__ == "__main__":
    raise SystemExit(main())
