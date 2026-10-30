"""Opt-in, fail-closed general-scalar T3 symbolic assembly + UV/EFT1 RG export.

Matchete: complete scalar-operator Lagrangian and detailed operator check.
RGBeta: all extra quartics + quadratic mass beta functions (UV and EFT1).
Does NOT claim mass-eigenbasis matching, physical Matchete/RGBeta coupling
normalisation, or EFT Wilson-operator RG closure.

PowerShell:
 python -X utf8 -m scripts.RunT3CompleteGeneralScalarSymbolic
 python -X utf8 -m scripts.RunT3CompleteGeneralScalarSymbolic --model A_-2 --run --run-id a_m2_v1
 python -X utf8 -m scripts.RunT3CompleteGeneralScalarSymbolic --all --run --run-id all16_v1
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import subprocess
from typing import Any

from scripts.AuditGeneralT3QuadraticBasis import audit, ROOT
from Numerical.experimental.GeneralT3ScalarQuadraticMasses import REP, load_manifest, quadratic_kind

REGISTRY = ROOT / "RGE/general/generated/general_rgbeta_quartic_registry.json"
MANIFEST = ROOT / "docs/general_scalar_trial_operator_manifest.json"
LAGRANGIAN = ROOT / "Lagrangian/RunFullGeneralScalarBuild.wl"
MAT_AUDIT = ROOT / "Lagrangian/RunGeneralQuadraticMatcheteAudit.wl"
RGBETA = ROOT / "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl"
OUT = ROOT / "output/studies/general_scalar_extended/complete_scalar_symbolic"
MASS_SYMBOLS = {"self_s1": "muT3GeneralS1Sq",
                "self_s2": "muT3GeneralS2Sq",
                "mixed": "muT3GeneralS1S2"}
MIXED = frozenset(("B_-1", "C_-1"))
FORMAT = "t3-full-symmetry-allowed-scalar-symbolic-integration-v1"


def read(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return obj


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def tag(model: str) -> str:
    return model.replace("+", "p").replace("-", "m")


def arguments(model: str) -> tuple[str, str, list[str]]:
    family, signed_alpha = model.split("_")
    alpha = int(signed_alpha)
    d1, d2, df = REP[family]
    return family, (f"m{-alpha}" if alpha < 0 else str(alpha)), [
        str(d1), str(d2), str(df),
    ]


def execute(name: str, command: list[str], path: Path, log: Path,
            *, env: dict[str, str] | None = None) -> dict:
    if path.exists():
        raise FileExistsError(f"Refusing to reuse stale output: {path}")
    try:
        p = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, check=False, env=env)
    except OSError as exc:
        raise RuntimeError(f"{name}: could not launch {command[0]}: {exc}") from exc
    log.write_text(p.stdout, encoding="utf-8")
    if p.returncode or not path.is_file():
        raise RuntimeError(f"{name}: Wolfram failed (exit={p.returncode}); "
                           f"see {log}; tail={p.stdout[-1200:]!r}")
    return read(path)


def validate_build(model: str, payload: dict, count: int) -> None:
    status = ("RequiresMixedScalarMassDiagonalization" if model in MIXED
              else "Success")
    if (payload.get("format") != "t3-full-general-scalar-lagrangian-candidate-v1"
        or payload.get("model") != model
        or payload.get("status") != status
        or payload.get("quadratic_expected") != count
        or payload.get("quadratic_included") != count
        or payload.get("quadratic_mass_basis_blocked") != (count if model in MIXED else 0)
        or payload.get("matching_mass_basis_ready") != (model not in MIXED)
        or payload.get("matching_performed") is not False
        or payload.get("physical_matchete_rgbeta_mapping_certified") is not False):
        raise ValueError(f"{model}: general Lagrangian build not validated: {payload}")


def validate_betas(model: str, stage: str, payload: dict,
                   entry: dict, kind: str) -> dict:
    meta = payload.get("metadata")
    names = payload.get("beta_symbols")
    betas = payload.get("report_betas")
    added = {q["name"] for q in entry["new_rgbeta_symbols"]}
    masses = ({MASS_SYMBOLS[kind]} if kind != "none" else set())
    if not isinstance(meta, dict) or not isinstance(betas, dict) or not isinstance(names, list):
        raise ValueError(f"{model}/{stage}: missing beta data")
    expected = {
        "format": "general-t3-quadratic-rgbeta-experimental-v1",
        "status": "SuccessCandidateMassBeta",
        "physical_mapping_certified": False,
        "mass_eigenstate_threshold_matching_performed": False,
    }
    for key, val in expected.items():
        if payload.get(key) != val or type(payload.get(key)) is not type(val):
            raise ValueError(f"{model}/{stage}: bad export field {key}")
    expected_meta = {
        "GeneralScalarModel": model,
        "GeneralScalarStage": stage,
        "GeneralScalarQuadraticMassMixing": "IncludedCandidate",
        "QuadraticKind": kind,
        "QuadraticMassSymbol": MASS_SYMBOLS.get(kind),
        "GeneralScalarAdditionalQuartics": entry["total_additional_quartics"],
        "GeneralScalarQuarticBetaBasis": "RGBetaFullQuarticProjection",
        "PhysicalMatcheteRGBetaMassBasisCertified": False,
    }
    for key, val in expected_meta.items():
        if meta.get(key) != val or type(meta.get(key)) is not type(val):
            raise ValueError(f"{model}/{stage}: bad beta metadata {key} ({meta.get(key)!r})")
    if (set(meta.get("GeneralScalarNewBetaKeys", [])) != added
            or set(names) != set(betas)
            or not (added | masses | {"mS1Sq", "mS2Sq"}) <= set(betas)
            or len(names) != len(set(names))):
        raise ValueError(f"{model}/{stage}: missing/duplicated beta function")
    for name, expr in betas.items():
        if (not isinstance(expr, str) or not expr.strip() or
            any(forbidden in expr for forbidden in
                ("BetaTerm[", "Missing[", "$Failed", "$Aborted", "Indeterminate"))):
            raise ValueError(f"{model}/{stage}: unresolved beta {name}")
    return {"beta_count": len(betas), "new_quartic_betas": len(added),
            "quadratic_mass_beta": MASS_SYMBOLS.get(kind),
            "diagonal_mass_beta_keys": ["mS1Sq", "mS2Sq"]}


def check_sources(models: list[str], registry: dict, manifest: dict,
                  basis: dict) -> list[dict]:
    for required in [LAGRANGIAN, MAT_AUDIT, RGBETA, MANIFEST, REGISTRY,
                     ROOT / "Lagrangian/model/LagrangianBuilder.wl",
                     ROOT / "Lagrangian/interactions/T3Topology.wl"]:
        if not required.is_file():
            raise FileNotFoundError(f"Required symbolic source missing: {required}")
    result = []
    for model in models:
        family, alpha, _ = arguments(model)
        entry = registry["models"][model]
        if (entry["alpha"] != int(model.split("_")[1])
            or entry["dimensions"] != list(REP[family])):
            raise ValueError(f"Registry representation mismatch: {model}")
        kind = quadratic_kind(manifest, model)
        expected = basis["models"][model]["independent_complex_mass_parameters_in_declared_basis"]
        if expected != int(kind != "none") or expected != entry["quadratics_not_implemented"]:
            raise ValueError(f"Quadratic inventory mismatch: {model}")
        result.append({"model": model, "quadratic_kind": kind,
                       "extra_bilinears": expected,
                       "extra_quartics": entry["total_additional_quartics"],
                       "matchete_heavy_mass_basis": "Blocked" if model in MIXED else "Permitted"})
    return result


def run_one(model: str, run_dir: Path, executable: str,
            entry: dict, kind: str) -> dict:
    family, alpha, dims = arguments(model)
    folder = run_dir / tag(model)
    folder.mkdir(parents=True, exist_ok=False)
    base = {"model": model, "quadratic_kind": kind, "status": "Failed",
            "stages": {}, "errors": [], "physical_mapping_certified": False,
            "mixed_scalar_matching_performed": False,
            "eft_wilson_running_performed": False}
    try:
        environment = os.environ.copy()
        environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] = "1"
        build_path = folder / "lagrangian.json"
        lag = execute("Lagrangian", [executable, "-file", str(LAGRANGIAN),
                    str(build_path), family, alpha], build_path,
                    folder / "lagrangian.log", env=environment)
        validate_build(model, lag, int(kind != "none"))
        base["stages"]["lagrangian"] = lag["status"]
        base["sources"] = {"lagrangian_sha256": digest(build_path),
                           "lagrangian_interaction_sha256": digest(
                               Path(lag["interaction_inputform_file"])),
                           "lagrangian_uv_sha256": digest(
                               Path(lag["uv_lagrangian_inputform_file"]))}
        # Uses the existing Matchete baseline-calibrated detailed checker.
        # B/C(-1) pass non-mass checks but intentionally fail HeavyMassBasis.
        from scripts.RunGeneralT3QuadraticMatcheteAudit import run as audit_matchete
        mat = audit_matchete(model, folder / "matchete", executable)
        base["stages"]["matchete"] = mat["status"]
        for stage in ("UV", "EFT1"):
            path = folder / f"{stage.lower()}_rgbeta.json"
            command = [executable, "-file", str(RGBETA), stage, *dims,
                       alpha, model, str(path)]
            data = execute(stage, command, path,
                           folder / f"{stage.lower()}_rgbeta.log")
            detail = validate_betas(model, stage, data, entry, kind)
            base["stages"][stage] = detail
            base["sources"][f"{stage.lower()}_rgbeta_sha256"] = digest(path)
        base["status"] = ("SymbolicRGExportedMassBasisBlocked" if model in MIXED
                          else "SymbolicRGExportedMatchingUnverified")
    except (OSError, RuntimeError, ValueError, KeyError, TypeError) as exc:
        base["errors"].append(f"{type(exc).__name__}: {exc}")
    (folder / "status.json").write_text(json.dumps(base, indent=2) + "\n",
                                        encoding="utf-8")
    return base


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    models = sorted(audit()["models"])
    group = ap.add_mutually_exclusive_group()
    group.add_argument("--model", choices=models)
    group.add_argument("--all", action="store_true")
    ap.add_argument("--run", action="store_true", help="Execute live Matchete/RGBeta")
    ap.add_argument("--run-id", default="symbolic_v1")
    ap.add_argument("--output-root", type=Path, default=OUT)
    ap.add_argument("--wolframscript", default="wolframscript")
    args = ap.parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", args.run_id):
        ap.error("Invalid --run-id")
    selected = models if args.all or args.model is None else [args.model]
    registry, manifest, basis = read(REGISTRY), load_manifest(MANIFEST), audit()
    inv = check_sources(selected, registry, manifest, basis)
    print("Symmetry-allowed T3 scalar inventory (live source verified):", flush=True)
    for line in inv:
        print(f"{line['model']}: bilinears={line['extra_bilinears']} "
              f"quartics={line['extra_quartics']} "
              f"MatcheteHeavyMassBasis={line['matchete_heavy_mass_basis']}", flush=True)
    if not args.run:
        print("Source checks only; add --run to execute symbolic UV/EFT1 exports.")
        return 0
    run_dir = args.output_root / args.run_id
    if run_dir.exists():
        raise FileExistsError(f"Refusing to overwrite previous symbolic campaign: {run_dir}")
    run_dir.mkdir(parents=True)
    rows = []
    for row in inv:
        model = row["model"]
        result = run_one(model, run_dir, args.wolframscript,
                         registry["models"][model], row["quadratic_kind"])
        rows.append(result)
        print(f"{model}: {result['status']} {result['errors'][:1]}", flush=True)
        (run_dir / "summary.json").write_text(json.dumps({
            "format": FORMAT, "models": rows,
            "matching_mass_basis_certified": False,
            "physical_matchete_rgbeta_mapping_certified": False,
            "full_sequential_wilson_matching_performed": False,
            "numerical_benchmarks_modified": False,
        }, indent=2) + "\n", encoding="utf-8")
    with (run_dir / "summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["model", "quadratic_kind", "status",
                                           "matchete", "uv_betas", "eft1_betas", "error"])
        w.writeheader()
        for r in rows:
            w.writerow({"model": r["model"], "quadratic_kind": r["quadratic_kind"],
                        "status": r["status"], "matchete": r["stages"].get("matchete"),
                        "uv_betas": r["stages"].get("UV", {}).get("beta_count"),
                        "eft1_betas": r["stages"].get("EFT1", {}).get("beta_count"),
                        "error": "; ".join(r["errors"])})
    return 0 if all(not r["errors"] for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
