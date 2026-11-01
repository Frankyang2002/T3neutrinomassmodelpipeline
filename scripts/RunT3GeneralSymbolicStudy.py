"""Reproducible, resumable, general-interaction symbolic T3 study (16 models).

This is the canonical entry point for the *general* UV BSM Lagrangian and
UV/EFT1 **renormalisable** one-loop RGEs. It reuses the repository's existing
Matchete and RGBeta implementations; no numerical values are introduced.

It deliberately does not perform extended-theory threshold matching, EFT1
higher-dimensional Wilson running, or final C5 matching. Those remain the
restricted pipeline's results until independently validated in the general
interaction basis. The two mixed-scalar-mass models are explicitly blocked.

Run from project root:
    python -X utf8 full16_symbolic.py
    python -X utf8 full16_symbolic.py --reports-only
    python -X utf8 full16_symbolic.py --model A_-2 --study general_a_m2
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import re
from typing import Any

from scripts import RunT3CompleteGeneralScalarSymbolic as scalar
from scripts.AuditGeneralT3QuadraticBasis import audit
from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
    REP,
    load_manifest,
    quadratic_kind,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_ROOT = ROOT / "output/studies/general_scalar_extended/complete_scalar_symbolic"
DEFAULT_REPORT_ROOT = ROOT / "Reports/output"
FORMAT = "t3-general-interactions-16-symbolic-study-v1"
EXTRA_Y_MODELS = {"A_-2", "A_+0", "D_-2", "D_+0", "E_-2", "E_+0"}
MIXED_MASS_MODELS = {"B_-1", "C_-1"}
REQUIRED_SOURCES = (
    "Lagrangian/RunFullGeneralScalarBuild.wl",
    "Lagrangian/model/LagrangianBuilder.wl",
    "Lagrangian/model/T3ModelCatalog.wl",
    "Lagrangian/model/T3Fields.wl",
    "Lagrangian/interactions/T3Topology.wl",
    "Lagrangian/interactions/ScalarPotential.wl",
    "Lagrangian/interactions/GeneralScalarOperators.wl",
    "Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl",
    "Lagrangian/model/SU2Invariants.wl",
    "Lagrangian/RunGeneralQuadraticMatcheteAudit.wl",
    "RGE/running/rgbeta/T3RGBetaModel.wl",
    "RGE/running/rgbeta/T3RGBetaRunnerCommon.wl",
    "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl",
    "RGE/general/generated/general_rgbeta_quartic_registry.json",
    "docs/general_scalar_trial_operator_manifest.json",
    "scripts/RunT3CompleteGeneralScalarSymbolic.py",
    "scripts/RunGeneralT3QuadraticMatcheteAudit.py",
)


def load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(data, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return data


def model_tag(model: str) -> str:
    return model.replace("+", "p").replace("-", "m")


def fingerprint() -> dict[str, str]:
    """Invalidate cached calculations when their principal source changes."""
    result: dict[str, str] = {}
    for relative in REQUIRED_SOURCES:
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(f"Missing symbolic source: {path}")
        result[relative] = sha256(path.read_bytes()).hexdigest()
    # The registered quartic component tensors are source code, not merely
    # generated output. A changed tensor basis must invalidate prior betas.
    generated = sorted((ROOT / "RGE/general/generated").glob("*_extra_quartics.wl"))
    if not generated:
        raise FileNotFoundError("No generated general-quartic RGBeta definitions found")
    for path in generated:
        result[path.relative_to(ROOT).as_posix()] = sha256(path.read_bytes()).hexdigest()
    return result


def expected_yukawa(model: str) -> str | None:
    if model not in EXTRA_Y_MODELS:
        return None
    return "y1Conj" if model.endswith("+0") else "y2Conj"


def validate_cached(
    model: str,
    folder: Path,
    entry: dict[str, Any],
    kind: str,
    current_sources: dict[str, str],
) -> dict[str, Any]:
    """Reject partial, stale, or restricted-basis files; never trust a filename."""
    stamp = load_json(folder / "provenance.json")
    if stamp.get("model") != model or stamp.get("source_sha256") != current_sources:
        raise ValueError("source fingerprint differs (or missing provenance)")
    lag = load_json(folder / "lagrangian.json")
    scalar.validate_build(model, lag, int(kind != "none"))
    expected_quartics = int(entry["total_additional_quartics"])
    if (lag.get("quartic_expected") != expected_quartics or
        lag.get("quartic_included") != expected_quartics or
        lag.get("extra_yukawas_requested") is not True or
        lag.get("general_scalar_quadratics_requested") is not True):
        raise ValueError("general scalar/Yukawa mode or quartic inventory missing")
    ykey = expected_yukawa(model)
    yn = int(ykey is not None)
    if lag.get("extra_yukawa_expected") != yn or lag.get("extra_yukawa_included") != yn:
        raise ValueError("extra Yukawa inventory incomplete")
    interactions = lag.get("interactions", [])
    if not isinstance(interactions, list) or len(interactions) != len(set(interactions)):
        raise ValueError("interaction inventory invalid or duplicated")
    quartic_names = [x for x in interactions if x.startswith("T3ExtraQ_")]
    if len(quartic_names) != expected_quartics:
        raise ValueError("general quartics not present in BSM Lagrangian")
    y_names = [x for x in interactions if x.startswith("Yukawa1ConjugateScalar_")
               or x.startswith("Yukawa2ConjugateScalar_")]
    if len(y_names) != yn:
        raise ValueError("independent extra Yukawa vertex missing or spurious")
    if y_names and ("Yukawa1" if ykey == "y1Conj" else "Yukawa2") not in y_names[0]:
        raise ValueError("incorrect extra Yukawa species")
    for key in ("interaction_inputform_file", "uv_lagrangian_inputform_file",
                "bsm_lagrangian_inputform_file"):
        raw = lag.get(key)
        if not isinstance(raw, str) or not Path(raw).is_file():
            raise ValueError(f"Missing full symbolic expression: {key}")
        if not Path(raw).read_text(encoding="utf-8").strip():
            raise ValueError(f"Empty symbolic expression: {key}")
    stages: dict[str, Any] = {}
    expected_keys = {ykey} if ykey else set()
    for stage in ("UV", "EFT1"):
        native = load_json(folder / f"{stage.lower()}_rgbeta.json")
        result = scalar.validate_betas(model, stage, native, entry, kind)
        rge_keys = set(native["report_betas"])
        observed = rge_keys & {"y1Conj", "y2Conj"}
        if observed != (expected_keys if stage == "UV" else set()):
            raise ValueError(f"{stage}: missing/spurious additional Yukawa beta")
        stages[stage] = {**result, "path": str(folder / f"{stage.lower()}_rgbeta.json")}
    status = load_json(folder / "status.json")
    if status.get("model") != model or status.get("errors") != []:
        raise ValueError("model-specific symbolic stage failed")
    matchete = status.get("stages", {}).get("matchete")
    expected_matchete = (
        "MatcheteNonMassBasisDetailedChecksPassedHeavyMassBasisBlocked"
        if model in MIXED_MASS_MODELS else "MatcheteFullDetailedChecksPassed"
    )
    if matchete != expected_matchete:
        raise ValueError(f"Matchete diagnostic is not expected status: {matchete}")
    return {
        "model": model,
        "status": status.get("status"),
        "extra_quartics": expected_quartics,
        "extra_bilinears": int(kind != "none"),
        "extra_yukawa": ykey,
        "mass_basis": "Blocked" if model in MIXED_MASS_MODELS else "CheckedNotMatched",
        "matching_performed": False,
        "physical_rgbeta_matchete_mapping_certified": False,
        "eft1_wilson_rge_performed": False,
        "lagrangian": lag,
        "rges": stages,
        "folder": str(folder),
    }


def archive_stale_folder(folder: Path) -> Path:
    """Keep earlier results; do not silently delete or overwrite a run."""
    for n in range(1, 10000):
        candidate = folder.with_name(f"{folder.name}.superseded_{n:03d}")
        if not candidate.exists():
            folder.rename(candidate)
            return candidate
    raise RuntimeError("Too many archived revisions for " + str(folder))


def build_model(
    model: str,
    run_dir: Path,
    registry: dict[str, Any],
    manifest: dict[str, Any],
    source_hashes: dict[str, str],
    wolframscript: str,
    reports_only: bool,
) -> dict[str, Any]:
    entry = registry["models"][model]
    kind = quadratic_kind(manifest, model)
    folder = run_dir / model_tag(model)
    if folder.exists():
        try:
            cached = validate_cached(model, folder, entry, kind, source_hashes)
            cached["cache"] = "Reused"
            print(f"{model}: verified source-matched symbolic files reused", flush=True)
            return cached
        except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
            if reports_only:
                return {"model": model, "status": "Unavailable", "error": str(exc)}
            previous = archive_stale_folder(folder)
            print(f"{model}: archived stale/incomplete results -> {previous.name} ({exc})", flush=True)
    elif reports_only:
        return {"model": model, "status": "Unavailable", "error": "No saved extended symbolic output"}

    # Existing scalar driver handles the live Matchete detailed audit and both
    # UV/EFT1 RGBeta passes. Extend only this subprocess environment to include
    # the gauge-allowed conjugate-scalar Yukawas in the six applicable models.
    old = os.environ.get("T3_EXPERIMENTAL_GENERAL_FERMIONS")
    os.environ["T3_EXPERIMENTAL_GENERAL_FERMIONS"] = "1"
    try:
        result = scalar.run_one(model, run_dir, wolframscript, entry, kind)
    finally:
        if old is None:
            os.environ.pop("T3_EXPERIMENTAL_GENERAL_FERMIONS", None)
        else:
            os.environ["T3_EXPERIMENTAL_GENERAL_FERMIONS"] = old
    if result.get("errors"):
        return {"model": model, "status": "Failed", "error": "; ".join(result["errors"])}
    (folder / "provenance.json").write_text(
        json.dumps({"model": model, "source_sha256": source_hashes}, indent=2) + "\n",
        encoding="utf-8",
    )
    try:
        checked = validate_cached(model, folder, entry, kind, source_hashes)
        checked["cache"] = "Generated"
        return checked
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
        return {"model": model, "status": "ValidationFailed", "error": str(exc)}


def su2_invariants(dimension: int) -> tuple[str, str]:
    # Exact quadratic Casimir and Dynkin index; T(fund)=1/2.
    return (
        str(Fraction(dimension * dimension - 1, 4)),
        str(Fraction(dimension * (dimension * dimension - 1), 12)),
    )


def link_to(target: Path, parent: Path) -> str:
    return Path(os.path.relpath(target, parent)).as_posix()


def report_one(row: dict[str, Any], report_dir: Path) -> None:
    model = row["model"]
    file = report_dir / "models" / (model_tag(model) + ".md")
    file.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# T3 {model}: general-interaction symbolic results", ""]
    if "error" in row:
        lines += [f"**Unavailable:** {row['error']}", ""]
        file.write_text("\n".join(lines), encoding="utf-8")
        return
    d1, d2, df = REP[model.split("_")[0]]
    a = int(model.split("_")[1])
    lines += [
        f"Fields: S1=({d1},{Fraction(a, 2)}), S2=({d2},{Fraction(a + 2, 2)}), "
        f"F=({df},{Fraction(a + 1, 2)}) [SU(2) dimension, hypercharge].",
        "",
        f"Additional independent quartics: **{row['extra_quartics']}**; "
        f"additional complex bilinears: **{row['extra_bilinears']}**; "
        f"additional Yukawa: **{row['extra_yukawa'] or 'none'}**.",
        "",
        f"Heavy scalar mass basis: **{row['mass_basis']}**. "
        "General sequential matching: **not performed**. "
        "Matchete/RGBeta physical coupling dictionary: **not certified**.",
        "",
        "## Independent interaction inventory",
        "",
    ]
    lag = row["lagrangian"]
    lines += [f"- `{name}`" for name in lag["interactions"]]
    lines += [
        "", "## Complete BSM UV Lagrangian (native Matchete InputForm)", "",
        "The expression below includes BSM free fields and interactions; "
        "it excludes the pure Standard Model Lagrangian.", "",
        "```text", Path(lag["bsm_lagrangian_inputform_file"]).read_text(encoding="utf-8").strip(),
        "```", "", "## One-loop renormalisable beta functions", "",
        "RGBeta's report convention is `16*pi^2*dX/dln(mu)`. "
        "These are candidate native RGBeta expressions, not independently "
        "mapped/certified Matchete coefficients.", "",
    ]
    for stage in ("UV", "EFT1"):
        native = load_json(Path(row["rges"][stage]["path"]))
        lines += [f"### {stage}: {len(native['report_betas'])} beta functions", ""]
        for name, expr in native["report_betas"].items():
            lines += [f"**beta_{name}**", "", "```text", expr, "```", ""]
    lines += [
        "## Interpretation and limitations", "",
        "- UV RGEs include symmetry-allowed scalar quartics, bilinears, and the "
        "additional Yukawa where gauge invariance permits it.",
        "- EFT1 RGEs describe the renormalisable SM+S1+S2 theory; the heavy "
        "fermion and its Yukawas have been removed.",
        "- EFT1 higher-dimensional Wilson RGEs and final C5 have **not** been "
        "rederived with the full general interactions.",
        "- Mixed scalar mass terms require mass-basis diagonalisation for B(-1) "
        "and C(-1) before general matching.",
        "- The existing heavy-fermion flavour-mass parameterisation is not "
        "certified as the most general non-diagonal flavour mass basis.",
        "",
    ]
    file.write_text("\n".join(lines), encoding="utf-8")


def write_reports(
    selected: list[str],
    rows: list[dict[str, Any]],
    report_dir: Path,
    run_dir: Path,
) -> None:
    report_dir.mkdir(parents=True, exist_ok=True)
    for row in rows:
        report_one(row, report_dir)
    with (report_dir / "group_invariants.csv").open("w", newline="", encoding="utf-8") as fp:
        writer = csv.writer(fp)
        writer.writerow(["model", "field", "d_SU2", "Y", "C2_SU2", "T_SU2"])
        for model in selected:
            cls, alpha_s = model.split("_")
            a = int(alpha_s)
            d1, d2, df = REP[cls]
            for label, d, y in (
                ("S1", d1, Fraction(a, 2)),
                ("S2", d2, Fraction(a + 2, 2)),
                ("F", df, Fraction(a + 1, 2)),
            ):
                c2, t2 = su2_invariants(d)
                writer.writerow([model, label, d, str(y), c2, t2])
    with (report_dir / "summary.csv").open("w", newline="", encoding="utf-8") as fp:
        writer = csv.writer(fp)
        writer.writerow(["model", "status", "extra_quartics", "extra_bilinears",
                         "extra_yukawa", "uv_betas", "eft1_betas", "mass_basis", "cache", "error"])
        for r in rows:
            writer.writerow([r["model"], r["status"], r.get("extra_quartics"),
                             r.get("extra_bilinears"), r.get("extra_yukawa"),
                             r.get("rges", {}).get("UV", {}).get("beta_count"),
                             r.get("rges", {}).get("EFT1", {}).get("beta_count"),
                             r.get("mass_basis"), r.get("cache"), r.get("error")])
    summary = {
        "format": FORMAT,
        "model_count_requested": len(selected),
        "model_count_validated": sum("error" not in r for r in rows),
        "all_requested_models_validated": len(rows) == len(selected) and all("error" not in r for r in rows),
        "data_root": str(run_dir),
        "model_statuses": [{k: v for k, v in r.items() if k not in {"lagrangian", "rges"}}
                           for r in rows],
        "all_symmetry_allowed_interaction_inventory_required": True,
        "mass_eigenstate_sequential_matching_performed": False,
        "general_wilson_operator_rge_performed": False,
        "physical_matchete_rgbeta_mapping_certified": False,
        "standard_pipeline_reports_promoted_to_general": False,
    }
    (report_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    index = [
        "# Full general-interaction T3 symbolic study", "",
        f"Validated: **{summary['model_count_validated']}/{len(selected)}** model(s).", "",
        "The BSM UV Lagrangian and **UV/EFT1 renormalisable beta functions** "
        "come from the general-interaction Matchete/RGBeta workflows, not from "
        "the restricted `pipeline.py` matching outputs.", "",
        "The beta expressions contain their group-dependent coefficients, but "
        "no independent Matchete-to-RGBeta physical coupling-map or term-by-term "
        "group-factor identity is claimed. `group_invariants.csv` gives the "
        "exact SU(2) Casimirs and Dynkin indices as a separate representation check.", "",
        "**Not computed:** general-interaction threshold Wilson matching, "
        "EFT1 higher-dimensional Wilson RGEs, and final Weinberg C5. "
        "For B(-1)/C(-1), mixed scalar masses require diagonalisation.", "",
        "| Model | Extra quartics | Extra bilinears | Extra Yukawa | Status | Report |",
        "|---|---:|---:|---|---|---|",
    ]
    for r in rows:
        model = r["model"]
        index.append(
            f"| {model} | {r.get('extra_quartics', '—')} | "
            f"{r.get('extra_bilinears', '—')} | {r.get('extra_yukawa') or '—'} | "
            f"{r['status']} | [Detail](models/{model_tag(model)}.md) |"
        )
    index += ["", "See also [exact group invariants](group_invariants.csv) and "
              "[machine-readable status](summary.json).", ""]
    (report_dir / "README.md").write_text("\n".join(index), encoding="utf-8")


def self_test() -> None:
    assert len(EXTRA_Y_MODELS) == 6 and len(MIXED_MASS_MODELS) == 2
    assert expected_yukawa("A_+0") == "y1Conj"
    assert expected_yukawa("D_-2") == "y2Conj"
    assert expected_yukawa("C_-1") is None
    assert su2_invariants(1) == ("0", "0")
    assert su2_invariants(2) == ("3/4", "1/2")
    assert su2_invariants(3) == ("2", "2")
    assert len(audit()["models"]) == 16
    print("PASS: 16 models; six extra Yukawas; two mixed scalar masses; SU(2) factors")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", help="Single model, e.g. A_-2; default all 16")
    parser.add_argument("--study", default="full16_general_symbolic", help="Report study folder")
    parser.add_argument("--run-id", default="full16_general_v1", help="Persistent data-run ID")
    parser.add_argument("--reports-only", action="store_true", help="Never run Wolfram; fail closed on missing/stale outputs")
    parser.add_argument("--wolframscript", default="wolframscript")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args(argv)
    if args.self_test:
        self_test()
        return 0
    for key in ("study", "run_id"):
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", getattr(args, key)):
            parser.error(f"--{key.replace('_', '-')} must be a simple directory name")
    basis = audit()
    models = sorted(basis["models"])
    if args.model is not None and args.model not in models:
        parser.error(f"Unknown model {args.model!r}; choose one of {', '.join(models)}")
    selected = [args.model] if args.model else models
    source_hashes = fingerprint()
    registry = load_json(scalar.REGISTRY)
    manifest = load_manifest(scalar.MANIFEST)
    scalar.check_sources(selected, registry, manifest, basis)
    run_dir = DEFAULT_DATA_ROOT / args.run_id
    report_dir = DEFAULT_REPORT_ROOT / args.study / "GeneralSymbolic"
    if not args.reports_only:
        run_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for n, model in enumerate(selected, 1):
        try:
            row = build_model(model, run_dir, registry, manifest, source_hashes,
                              args.wolframscript, args.reports_only)
        except (OSError, RuntimeError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
            row = {"model": model, "status": "Failed", "error": str(exc)}
        rows.append(row)
        print(f"[{n}/{len(selected)}] {model}: {row['status']}"
              + (" — " + row["error"] if "error" in row else ""), flush=True)
        write_reports(selected[:len(rows)], rows, report_dir, run_dir)
    success = all("error" not in r for r in rows)
    print(f"Validated {sum('error' not in r for r in rows)}/{len(selected)}. "
          f"Reports: {report_dir / 'README.md'}", flush=True)
    print("Extended threshold matching / physical RGBeta dictionary: NOT CERTIFIED", flush=True)
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
