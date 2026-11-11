"""Run and validate candidate UV RGBeta mass/quartic exports for all 16 T3 models.

Project path: scripts/RunBatch13UVExportsAll16.py
Run at repository root with: python scripts/RunBatch13UVExportsAll16.py

This is a symbolic *candidate* export audit, not numerical UV running or
Matchete/RGBeta normalisation certification. Does not modify source files.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "RGE/general/generated/general_rgbeta_quartic_registry.json"
RUNNER = ROOT / "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl"
DEFAULT_OUTPUT = ROOT / "output/batch13/all16_uv"

CLASS_DIMENSIONS = {
    "A": (1, 3, 2), "B": (2, 2, 1), "C": (2, 2, 3),
    "D": (3, 1, 2), "E": (3, 3, 2),
}
CLASS_ALPHAS = {
    "A": (0, -2, -4), "B": (1, -1, -3),
    "C": (1, -1, -3), "D": (0, 2, -2), "E": (0, 2, -2, -4),
}
EXPECTED_QUADRATIC = {
    "A_+0": "muT3GeneralS1Sq", "A_-2": "muT3GeneralS2Sq",
    "B_-1": "muT3GeneralS1S2", "C_-1": "muT3GeneralS1S2",
    "D_+0": "muT3GeneralS1Sq", "D_-2": "muT3GeneralS2Sq",
    "E_+0": "muT3GeneralS1Sq", "E_-2": "muT3GeneralS2Sq",
}


def model_name(label: str, alpha: int) -> str:
    return f"{label}_{alpha:+d}"


def integer_token(number: int) -> str:
    return f"m{-number}" if number < 0 else str(number)


def expected_models() -> list[tuple[str, tuple[int, int, int], int]]:
    return [(model_name(label, alpha), CLASS_DIMENSIONS[label], alpha)
            for label, values in CLASS_ALPHAS.items() for alpha in values]


def load_registry(path: Path = REGISTRY) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    models = manifest.get("models")
    expected = expected_models()
    if not isinstance(models, dict) or set(models) != {name for name, _, _ in expected}:
        raise ValueError("Registry does not match the complete 16-model T3 catalogue")
    for name, dims, alpha in expected:
        entry = models[name]
        if entry.get("dimensions") != list(dims) or entry.get("alpha") != alpha:
            raise ValueError(f"Registry representation mismatch: {name}")
        mass_count = entry.get("quadratics_not_implemented")
        if mass_count != int(name in EXPECTED_QUADRATIC):
            raise ValueError(f"Registry additional-quadratic mismatch: {name}")
    return models


def validate_export(name: str, dims: tuple[int, int, int], alpha: int,
                    entry: dict[str, Any], path: Path) -> list[str]:
    errors: list[str] = []
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return [f"JSON read/parse failed: {exc}"]
    if not isinstance(data, dict):
        return ["Export is not a JSON object"]
    if data.get("status") != "SuccessCandidateMassBeta":
        errors.append(f"Unexpected export status: {data.get('status')!r}")
    meta = data.get("metadata")
    if not isinstance(meta, dict):
        return errors + ["Missing metadata object"]
    if meta.get("GeneralScalarModel") != name or meta.get("GeneralScalarStage") != "UV":
        errors.append("Model name or UV stage mismatch")
    for key, val in zip(("dS1", "dS2", "dF", "alpha"), (*dims, alpha)):
        if meta.get(key) != val:
            errors.append(f"Metadata mismatch: {key}")
    if meta.get("QuadraticMassSymbol") != EXPECTED_QUADRATIC.get(name):
        errors.append("Additional quadratic mass symbol mismatch")
    if meta.get("GeneralScalarAdditionalQuartics") != entry["total_additional_quartics"]:
        errors.append("Additional quartic count mismatch")
    if meta.get("PhysicalMatcheteRGBetaMassBasisCertified") is not False:
        errors.append("Unexpected physical mass-basis certification")
    if data.get("physical_mapping_certified") is not False:
        errors.append("Unexpected physical mapping certification")
    symbols = data.get("beta_symbols")
    betas = data.get("report_betas")
    if not isinstance(symbols, list) or not all(isinstance(s, str) for s in symbols):
        return errors + ["Invalid beta_symbols list"]
    if len(symbols) != len(set(symbols)):
        errors.append("Repeated beta symbols")
    if not isinstance(betas, dict):
        return errors + ["Invalid report_betas object"]
    if set(betas) != set(symbols):
        errors.append("Beta symbol/report key mismatch")
    required = {"gY", "g2", "g3", "yu", "yd", "ye", "y1", "y2",
                "MF", "mS1Sq", "mS2Sq", "lambdaH", "lambdaT3"}
    required.update(q["name"] for q in entry["new_rgbeta_symbols"])
    extra = EXPECTED_QUADRATIC.get(name)
    if extra:
        required.add(extra)
    missing = sorted(required - set(symbols))
    if missing:
        errors.append(f"Missing beta symbols: {missing}")
    forbidden = sorted({"y1Conj", "y2Conj"} & set(symbols))
    if forbidden:
        errors.append(f"Restricted Yukawa violation: {forbidden}")
    for symbol, expr in betas.items():
        if not isinstance(expr, str) or not expr.strip():
            errors.append(f"Empty/non-string beta expression: {symbol}")
        elif any(marker in expr for marker in ("$Failed", "$Aborted", "Indeterminate", "BetaTerm[")):
            errors.append(f"Unresolved beta expression: {symbol}")
    return errors


def execute(*, executable: str, output_dir: Path, registry_path: Path,
            runner_path: Path, dry_run: bool = False) -> tuple[int, dict[str, Any]]:
    models = load_registry(registry_path)
    if not runner_path.is_file():
        raise FileNotFoundError(f"Missing Wolfram runner: {runner_path}")
    output_dir.mkdir(parents=True, exist_ok=True)
    records: dict[str, Any] = {}
    for name, dims, alpha in expected_models():
        output = output_dir / f"{name}_UV.json"
        logfile = output_dir / f"{name}_UV.log"
        diagfile = Path(str(output) + ".diagnostic.json")
        cmd = [executable, "-file", str(runner_path), "UV",
               *(str(x) for x in dims), integer_token(alpha), name, str(output)]
        record: dict[str, Any] = {"command": cmd, "output": str(output),
                                  "log": str(logfile), "status": "pending", "errors": []}
        records[name] = record
        print(f"[{name}] {'DRY RUN' if dry_run else 'RUN'}: {' '.join(cmd)}", flush=True)
        if dry_run:
            record["status"] = "dry_run"
            continue
        # Stale outputs must never be counted as current successful results.
        for stale in (output, diagfile):
            stale.unlink(missing_ok=True)
        env = os.environ.copy()
        # UV restriction is explicit; do not inherit an experimental setting.
        env["T3_EXPERIMENTAL_GENERAL_FERMIONS"] = "0"
        try:
            with logfile.open("w", encoding="utf-8", errors="replace") as log:
                result = subprocess.run(cmd, cwd=ROOT, env=env, stdout=log,
                                        stderr=subprocess.STDOUT, check=False)
            record["exit_code"] = result.returncode
            if result.returncode != 0:
                record["errors"].append(f"wolframscript exit code {result.returncode}")
            if output.exists():
                record["errors"].extend(validate_export(name, dims, alpha, models[name], output))
            else:
                record["errors"].append("No export JSON created")
        except OSError as exc:
            record["errors"].append(f"Process launch failed: {exc}")
        record["status"] = "failed" if record["errors"] else "passed"
        print(f"  {record['status'].upper()}" +
              (f": {'; '.join(record['errors'])}" if record["errors"] else ""), flush=True)
    passed = sum(r["status"] == "passed" for r in records.values())
    failed = sum(r["status"] == "failed" for r in records.values())
    summary = {"format": "batch13-uv-export-audit-v1", "stage": "UV",
               "restricted_yukawa": True, "model_count": len(records),
               "passed": passed, "failed": failed, "dry_run": dry_run,
               "physical_normalisation_certified": False,
               "numerical_running_performed": False, "models": records}
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Summary: {summary_path} | passed={passed} failed={failed}", flush=True)
    return (1 if failed else 0), summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--wolframscript", default="wolframscript")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    try:
        rc, _ = execute(executable=args.wolframscript,
                        output_dir=args.output_dir.resolve(),
                        registry_path=REGISTRY, runner_path=RUNNER,
                        dry_run=args.dry_run)
        return rc
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"Batch 13 configuration error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
