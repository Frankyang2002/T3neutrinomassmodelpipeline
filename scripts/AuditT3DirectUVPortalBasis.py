"""Audit native Matchete Higgs-portal mass contributions; DO NOT assume a UV/RGBeta map.

For each of the 16 frozen benchmarks this script reads independently exported,
registered native HiggsPortal{1,2}Inv{k} CG tensors, rotates into the verified
candidate charge basis, and projects H0^*H0. The native Matchete coupling
lambdaH{which}Inv{k} is NOT the numerical RGBeta lambdaH{which} without a
separately established coupling dictionary. In particular, the current
provisional adapter assumes a universal -lambdaH{which}v^2/2 scalar shift.

No EFT code, matching, running, benchmark, or existing provisional result is
changed. This is a source-tensor audit, not physical normalization certification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.RunT3DirectOneLoop import models, CLASSES  # noqa: E402
from scripts.BuildT3DirectOneLoopInputs import _rotations, _tensor, _matrix  # noqa: E402

FLAGS = (True, False, True, False)  # Bar[H], H, Bar[Si], Si
FORMAT = "t3-direct-uv-native-portal-projector-audit-v1"


def require(condition: bool, explanation: str) -> None:
    if not condition:
        raise ValueError(explanation)


def read_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"Required source not found: {path}")
    out = json.loads(path.read_text(encoding="utf-8-sig"))
    require(isinstance(out, dict), f"Not a JSON object: {path}")
    return out


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_contract(root: Path) -> dict[str, str]:
    """Fail if source field ordering or PlusHc potential prefactors change."""
    name = root / "Lagrangian/interactions/ScalarPotential.wl"
    topology = root / "Lagrangian/interactions/T3Topology.wl"
    source = name.read_text(encoding="utf-8-sig")
    tsource = topology.read_text(encoding="utf-8-sig")
    start = source.find("BuildHiggsPortalInvariantCandidates[")
    end = source.find("BuildCrossScalarInvariantCandidates[", start)
    require(start >= 0 and end > start, "Portal source function not found")
    fn = source[start:end]
    expected = (
        '"T3HiggsPortal" <> ToString[which] <> "CG"',
        '"lambdaH" <> ToString[which]',
        "{2, 2, d, d}",
        "{True, False, True, False}",
        "{1, 2}",
        '"Expression" ->',
        "1/2 PlusHc[",
        "coupling[] fieldProduct",
    )
    require(all(token in fn for token in expected),
            "Native Matchete portal source changed (orientation, factor, or fields)")
    require("lambdaT3[] H[i] H[j] s1 s2 CG[T3MixCG, labels]" in tsource,
            "Native T3 scalar-mixing source changed")
    return {str(name.relative_to(root)): digest(name),
            str(topology.relative_to(root)): digest(topology)}


def projected_neutral_portal(tensor: np.ndarray, dim: int) -> list[float]:
    """Extract native invariant at H0*H0 and enforce charge-diagonal Si*Si.

    Tensor is ordered Bar[H], H, Bar[Si], Si in the candidate T3 weight basis.
    The neutral H0 is weight index 1 (T3=-1/2, Y=+1/2).
    """
    arr = np.asarray(tensor, dtype=np.complex128)
    require(arr.shape == (2, 2, dim, dim), "Wrong Higgs-portal CG tensor shape")
    require(np.isfinite(arr).all(), "Nonfinite portal component")
    sub = arr[1, 1, :, :]
    scale = max(1., float(np.max(np.abs(sub))))
    residual = float(np.linalg.norm(sub-np.diag(np.diag(sub))))
    require(residual <= 1e-11 * scale,
            f"H0*H0 portal has off-diagonal charge components: {residual:g}")
    # Source is 1/2 PlusHc[lambda*...], with lambdaH{which}Inv{k}
    # declared real. The physically relevant quadratic form is the
    # Hermitian part of the native tensor, not the unprojected complex CG.
    hermitian_projector = (sub + sub.conj().T) / 2
    diag = np.diag(hermitian_projector)
    require(float(np.max(np.abs(diag.imag))) <= 1e-11 * scale,
            "Projected native portal diagonal is not real")
    return [float(v.real) for v in diag]


def _native_portal_components(seed: dict[str, Any], dim: int, which: int,
                              us: dict[int, sp.Matrix]) -> list[dict[str, Any]]:
    from RGE.running.intermediate.IntermediateMatcheteParsing import load_cg_registry
    from Numerical.diagnostics.AuditT3ComponentMultiplicity import expand_tensor
    registry = load_cg_registry(seed)
    found: list[dict[str, Any]] = []
    for record in seed.get("CGRegistry", []):
        require(record.get("WhichScalar") in (1, 2), "Unexpected scalar label")
        if record["WhichScalar"] != which:
            continue
        k = int(record["Invariant"])
        name = f"T3HiggsPortal{which}CGInv{k}"
        coupling = f"lambdaH{which}Inv{k}"
        require(record["Name"] == name and record["Coupling"] == coupling
                and record["SourceInteraction"] == f"HiggsPortal{which}Inv{k}",
                f"Unexpected Matchete portal name/coupling: {record}")
        dims = (2, 2, dim, dim)
        require(tuple(record["Dimensions"]) == dims,
                f"Wrong native portal dimensions for {name}")
        require(name in registry, f"Missing registered tensor {name}")
        components = [
            {"indices": list(idx), "cg": str(value)}
            for idx, value in expand_tensor(registry[name].tensor, dims)
        ]
        weight = _tensor({"dimensions": list(dims), "components": components}, FLAGS, us)
        diag = projected_neutral_portal(weight, dim)
        found.append({
            "native_coupling": coupling,
            "source_interaction": record["SourceInteraction"],
            "invariant_number": k,
            "H0bar_H0_scalar_diagonal": diag,
            "conditional_identification_with_universal_identity": bool(
                np.allclose(diag, np.ones(dim), rtol=0, atol=1e-11)
            ),
            "projector_spread": float(max(diag)-min(diag)),
        })
    found.sort(key=lambda row: row["invariant_number"])
    require(1 <= len(found) <= 2 and
            [row["invariant_number"] for row in found] == list(range(1, len(found)+1)),
            f"Missing or nonconsecutive native portal invariants for scalar {which}")
    return found


def project_model(seed: dict[str, Any], model: str, expected: tuple[int, int, int, int],
                  us: dict[int, sp.Matrix], ordinary: dict[str, Any], v: float,
                  source_path: Path) -> dict[str, Any]:
    require(seed.get("status") == "Success" and
            seed.get("format") == "t3-native-direct-uv-portal-cg-v1",
            f"Invalid native portal registry for {model}")
    got = tuple(int(seed["metadata"][k]) for k in ("dS1", "dS2", "dF", "alpha"))
    require(got == expected, f"Native portal model metadata mismatch for {model}: {got}")
    d1, d2, _, _ = expected
    rows = []
    for which, dim in ((1, d1), (2, d2)):
        portal = _native_portal_components(seed, dim, which, us)
        # This is a WHAT-IF comparison, not a proposed source/RGBeta dictionary.
        legacy = float(ordinary[f"lambdaH{which}"])
        assumed_shift = -legacy * v*v/2
        conditional_native_shift = [
            -legacy * v*v/2 * z for z in portal[0]["H0bar_H0_scalar_diagonal"]
        ]
        rows.append({
            "scalar": f"S{which}",
            "baseline_rgbeta_parameter": f"lambdaH{which}",
            "baseline_rgbeta_value": legacy,
            "native_matchete_invariants": portal,
            "provisional_adapter_universal_shift_gev2": assumed_shift,
            "if_lambdaHInv1_equals_lambdaH_and_other_invariants_zero": {
                "source_invariant_1_scalar_shifts_gev2": conditional_native_shift,
                "max_abs_difference_from_adapter_gev2": float(max(
                    abs(shift-assumed_shift) for shift in conditional_native_shift
                )),
                "is_assumed_universal_identity": portal[0]["conditional_identification_with_universal_identity"],
            },
        })
    require(len(seed["CGRegistry"]) == sum(len(x["native_matchete_invariants"]) for x in rows),
            "Unexpected extra Matchete portal tensor")
    return {"model": model, "dimensions_alpha": list(expected),
            "native_export_path": str(source_path),
            "native_export_sha256": digest(source_path),
            "scalars": rows,
            "physical_rgbeta_to_matchete_portal_dictionary_certified": False}


def run(project: Path, portal_root: Path, audit_dir: Path,
        output_file: Path) -> dict[str, Any]:
    source_hashes = source_contract(project)
    analysis = read_json(audit_dir / "matchete_generator_analysis.json")
    covariance = read_json(audit_dir / "cg_covariance_audit.json")
    require(covariance.get("status") == "cg_generator_covariance_verified" and
            len(covariance.get("models", [])) == 16 and
            all(row.get("all_covariant") is True for row in covariance["models"]),
            "Prior all-16 source CG covariance check is missing")
    us = _rotations(analysis)
    results = []
    for short, frozen_name in models():
        frozen_path = project / "Numerical/benchmarks/fitted_reference" / (frozen_name + ".json")
        frozen = read_json(frozen_path)
        r = frozen["representation"]
        dims = tuple(int(r[k]) for k in ("d_s1", "d_s2", "d_f", "alpha"))
        expected_dims = (*CLASSES[short[0]][:3], int(short[2:]))
        require(dims == expected_dims,
                f"Frozen representation mismatch for {short}: {dims} != {expected_dims}")
        source = portal_root / frozen_name / "native_higgs_portal_registry.json"
        obj = project_model(read_json(source), short, dims, us,
                            frozen["base_state"]["ordinary"],
                            float(frozen["numerical"]["vev_gev"]),source)
        obj["frozen_benchmark_sha256"] = digest(frozen_path)
        results.append(obj)
        print(f"{short:6s} native portal invariants: " +
              ", ".join(f"{s['scalar']}={len(s['native_matchete_invariants'])}" for s in obj["scalars"]),
              flush=True)
    require(len(results) == 16 and len({x["model"] for x in results}) == 16,
            "Expected 16 distinct audited models")
    report = {
        "format": FORMAT,
        "status": "native_matchete_portal_CG_projectors_extracted_physical_coupling_map_pending",
        "models_checked": 16,
        "source_sha256": source_hashes,
        "source_higgs_vev_component_zero_based": 1,
        "source_scalar_potential_convention":
            "+1/2 PlusHc[lambdaHInvk Hdag H Sdag S CG] in L; V=-L",
        "limitation": "Native portal couplings lambdaH1Invk/lambdaH2Invk are NOT yet mapped "
                      "to RGBeta/frozen lambdaH1/lambdaH2. Invariant-1-only comparison "
                      "is explicitly conditional. Does not establish physical UV masses.",
        "physical_portal_dictionary_certified": False,
        "models": results,
    }
    output_file.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if output_file.exists():
        require(output_file.read_text(encoding="utf-8-sig") == text,
                f"Refusing to overwrite changed existing output: {output_file}")
    else:
        with output_file.open("x", encoding="utf-8") as handle:
            handle.write(text)
    print("Native portal projectors evaluated: 16/16")
    print("Physical Matchete/RGBeta portal coupling dictionary: NOT CERTIFIED")
    print("Output:", output_file)
    return report


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, default=ROOT)
    p.add_argument("--portal-root", type=Path, required=True)
    p.add_argument("--audit-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    try:
        run(args.project_root.resolve(), args.portal_root.resolve(),
            args.audit_dir.resolve(), args.output.resolve())
    except (ValueError, OSError, KeyError, TypeError, IndexError, sp.SympifyError) as exc:
        p.exit(2, "BLOCKED: " + str(exc) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
