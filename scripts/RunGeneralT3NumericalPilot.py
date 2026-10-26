"""Evaluate the extra B/C(-1) Matchete hard-C5 response, not production RGEs."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Numerical.experimental.GeneralT3C5Pilot import load_kernels, evaluate_response


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("output/general_t3_trial"))
    parser.add_argument("--model", choices=("B_m1", "C_m1"), required=True)
    parser.add_argument("--input", type=Path, required=True,
                        help="Numerical JSON with MF, y1_real/imag, y2_real/imag and Matchete extra couplings")
    parser.add_argument("--output", type=Path, required=True)
    ns = parser.parse_args()
    try:
        atoms = ns.root / ns.model / "general_c5_algebraic" / "delta_c5_atoms.json"
        kernels = load_kernels(atoms, ns.model)
        result = evaluate_response(kernels, json.loads(ns.input.read_text(encoding="utf-8-sig")))
        result["model"] = ns.model
        result["atom_source"] = str(atoms)
        result["input_file"] = str(ns.input)
        result["kernels"] = [
            {"coupling": k.coupling, "conjugated": k.conjugated,
             "yukawa": k.yukawa, "finite_hard_kernel": str(k.expression)} for k in kernels
        ]
        ns.output.parent.mkdir(parents=True, exist_ok=True)
        ns.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(ns.model, "hard C5 numerical response =", result["delta_c5"]["frobenius_norm"], "GeV^-1")
        print("Report:", ns.output)
        return 0
    except (ValueError, KeyError, TypeError, OSError) as exc:
        print(f"BLOCKED: {exc}", file=sys.stderr)
        # Preserve machine-readable diagnostic even when a kernel cannot be
        # evaluated.  Keep full original expressions for an exact follow-up.
        try:
            atoms = ns.root / ns.model / "general_c5_algebraic" / "delta_c5_atoms.json"
            payload = json.loads(atoms.read_text(encoding="utf-8-sig")) if atoms.is_file() else {}
            diagnostic = {
                "status": "Blocked", "model": ns.model, "reason": str(exc),
                "source": str(atoms),
                "raw_response_atoms": payload.get("atoms", []),
                "input_file": str(ns.input),
            }
            ns.output.parent.mkdir(parents=True, exist_ok=True)
            ns.output.write_text(json.dumps(diagnostic, indent=2) + "\n", encoding="utf-8")
            print(f"Diagnostic report: {ns.output}", file=sys.stderr)
        except (OSError, ValueError):
            pass
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
