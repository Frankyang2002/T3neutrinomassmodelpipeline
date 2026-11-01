"""Run real general-T3 UV -> EFT1 RGBeta exports through the existing numerical runners.

This is an orchestration entry point, not a second ODE implementation.  It
uses Numerical.running.UVRunner and IntermediateScalarRunner unchanged.
The scalar/Weinberg matching calculation is NOT performed here.

Example (PowerShell):
 python -X utf8 -m Numerical.orchestration.GeneralRunningBenchmark `
   --uv-json path/to/A_m2_uv.json --eft1-json path/to/A_m2_eft1.json `
   --output-json output/general_A_m2_running.json
"""
from __future__ import annotations

import argparse
from dataclasses import fields
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from Numerical.core.GeneralCouplingLayout import spec_for
from Numerical.core.State import SMNumericalState, T3Representation, T3UVState
from Numerical.core.StateVector import pack_uv_state
from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state
from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
from Numerical.running.UVRunner import run_uv_segment


REP = {"A": (1, 3, 2), "B": (2, 2, 1), "C": (2, 2, 3),
       "D": (3, 1, 2), "E": (3, 3, 2)}
_UV_FIELDS = frozenset(f.name for f in fields(T3UVState))
_SM_NAMES = frozenset(("gY", "g2", "g3", "yu", "yd", "ye", "lambdaH"))
_MATRIX_NAMES = frozenset(("yu", "yd", "ye", "y1", "y2", "MF", "y1Conj", "y2Conj"))


def _load(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected an object in {path}")
    return value


def _model_representation(payload: Mapping[str, Any]) -> tuple[str, T3Representation]:
    meta = payload.get("metadata")
    if not isinstance(meta, Mapping):
        raise ValueError("UV payload missing metadata")
    try:
        dims = (int(meta["dS1"]), int(meta["dS2"]), int(meta["dF"]))
        alpha = int(meta["alpha"])
    except (KeyError, ValueError, TypeError) as exc:
        raise ValueError("UV payload missing valid T3 representation metadata") from exc
    matches = [name for name, d in REP.items() if d == dims]
    if len(matches) != 1:
        raise ValueError(f"Not one of the ordinary T3 representations: {dims}")
    rep = T3Representation(*dims, alpha, shared_scalar=False).validated()
    model = f"{matches[0]}_{alpha:+d}"
    declared = meta.get("GeneralScalarModel")
    if declared is not None and declared != model:
        raise ValueError(f"Model mismatch: {declared!r} != {model!r}")
    return model, rep


def _baseline(name: str, *, mu_f: float, mu_s: float) -> Any:
    if name in ("gY", "g2", "g3"):
        return {"gY": 0.36, "g2": 0.65, "g3": 1.06}[name]
    if name == "lambdaH":
        return 0.13
    if name == "yu":
        return np.diag([1e-5, 0.005, 0.8]).astype(complex)
    if name == "yd":
        return np.diag([2e-5, 5e-4, 0.02]).astype(complex)
    if name == "ye":
        return np.diag([3e-6, 6e-4, 0.01]).astype(complex)
    if name == "MF":
        return np.diag([mu_f * 0.9, mu_f, mu_f * 1.1]).astype(complex)
    if name in ("y1", "y2", "y1Conj", "y2Conj"):
        return np.eye(3, dtype=complex) * 0.02
    if name in ("mS1Sq", "mS2Sq"):
        return float((1.3 * mu_s) ** 2 if name == "mS1Sq" else (1.5 * mu_s) ** 2)
    if name.startswith("muT3General"):
        return complex(0.02 * mu_s**2, 0.01 * mu_s**2)
    if spec_for(name).kind == "real":
        return 0.02
    return 0.01 + 0.002j


def _read_complex(value: Any, *, name: str) -> complex:
    if isinstance(value, Mapping) and set(value) == {"real", "imag"}:
        return complex(float(value["real"]), float(value["imag"]))
    if isinstance(value, (int, float)):
        return complex(value)
    raise ValueError(f"{name}: complex values must be numbers or {{real, imag}}")


def _override_value(name: str, value: Any) -> Any:
    spec = spec_for(name)
    if spec.kind == "matrix":
        if not isinstance(value, list):
            raise ValueError(f"{name}: expected nested 3x3 matrix")
        return spec.validate(np.array([[_read_complex(x, name=name) for x in row] for row in value]))
    if spec.kind == "complex":
        return spec.validate(_read_complex(value, name=name))
    return spec.validate(value)


def initial_state(uv_payload: Mapping[str, Any], *, mu_uv: float,
                  mu_f: float, mu_s: float,
                  overrides: Mapping[str, Any] | None = None) -> tuple[str, T3UVState]:
    """Construct one explicit benchmark using *all* UV export coupling keys."""
    model, rep = _model_representation(uv_payload)
    beta = uv_payload.get("report_betas")
    if not isinstance(beta, Mapping) or not beta:
        raise ValueError("UV report_betas is empty or malformed")
    values = {name: _baseline(name, mu_f=mu_f, mu_s=mu_s) for name in beta}
    if overrides:
        unknown = set(overrides).difference(values)
        if unknown:
            raise ValueError(f"Benchmark overrides not in UV RGBeta export: {sorted(unknown)}")
        values.update({name: _override_value(name, value) for name, value in overrides.items()})
    for name, value in values.items():
        spec_for(name).validate(value, majorana=rep.majorana_fermion)
    missing_sm = _SM_NAMES.difference(values)
    if missing_sm:
        raise ValueError(f"UV export missing required SM couplings: {sorted(missing_sm)}")
    sm = SMNumericalState(**{name: values[name] for name in _SM_NAMES})
    mandatory = ("y1", "y2", "MF", "mS1Sq", "mS2Sq", "lambdaS1", "lambdaS2",
                 "lambdaH1", "lambdaH2", "lambda12", "lambdaT3")
    absent = set(mandatory).difference(values)
    if absent:
        raise ValueError(f"UV export missing required ordinary T3 couplings: {sorted(absent)}")
    args = {name: values[name] for name in mandatory}
    # The standard state requires existing representation-specific structural quartics.
    # They MUST have beta functions in the UV payload: never create unrun state values.
    special = rep.d_s1 == 2 and rep.d_s2 == 2 and rep.alpha == -1
    req = {"lambdaH1Adj": rep.d_s1 > 1, "lambdaH2Adj": rep.d_s2 > 1,
           "lambdaS1Adj": rep.d_s1 == 3, "lambdaS2Adj": rep.d_s2 == 3,
           "lambda12Adj": rep.d_s1 > 1 and rep.d_s2 > 1,
           "lambda12Cross": rep.d_s1 == 3 and rep.d_s2 == 3}
    complex_old = ("lambdaHHdagS2S2", "lambdaHHdagS1barS1bar",
                   "lambdaS1bar2S2bar2", "lambdaS1barS2S2bar2",
                   "lambdaS1S1bar2S2bar", "lambdaHHdagS1barS2barCross")
    req.update({name: special for name in complex_old})
    for name, required in req.items():
        if required and name not in values:
            raise ValueError(f"UV export missing required structural quartic {name}")
        if not required and name in values:
            raise ValueError(f"UV export has unexpected structural quartic {name}")
        args[name] = values.get(name)
    standard = set(mandatory) | set(req) | _SM_NAMES
    extra = {name: value for name, value in values.items() if name not in standard}
    state = T3UVState(mu_gev=mu_uv, representation=rep, sm=sm,
                      general_couplings=extra, **args).validated()
    _, layout = pack_uv_state(state)
    if set(layout.names) != set(beta):
        raise ValueError(f"UV state/export mismatch: missing={sorted(set(beta)-set(layout.names))}, "
                         f"extra={sorted(set(layout.names)-set(beta))}")
    return model, state


def _encode(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return _encode(value.tolist())
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.generic):
        return _encode(value.item())
    if isinstance(value, (list, tuple)):
        return [_encode(x) for x in value]
    if isinstance(value, Mapping):
        return {str(k): _encode(v) for k, v in value.items()}
    return value


def _state_values(state: Any) -> dict[str, Any]:
    if isinstance(state, T3UVState):
        _, layout = pack_uv_state(state)
    else:
        _, layout = pack_intermediate_scalar_state(state)
    sm = {name: getattr(state.sm, name) for name in _SM_NAMES}
    rest = {name: getattr(state, name) for name in layout.names
            if name not in sm and name not in state.general_couplings}
    return _encode({**sm, **rest, **state.general_couplings})


def load_canonical_general_exports(model: str, run_id: str, *,
                                   generate: bool = False,
                                   wolframscript: str = "wolframscript",
                                   study: str = "general_numerical") -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Load only source-fingerprint-validated *full-general* symbolic exports.

    The canonical builder validates the Matchete interaction inventory, extra
    conjugate Yukawa, exact RGBeta keys, and source fingerprints. It does not
    establish the physical Matchete <-> RGBeta parameter dictionary.
    """
    from scripts import RunT3GeneralSymbolicStudy as canonical
    from scripts import RunT3CompleteGeneralScalarSymbolic as scalar
    from scripts.AuditGeneralT3QuadraticBasis import audit
    from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
        load_manifest, quadratic_kind,
    )

    if not isinstance(run_id, str) or not run_id or not all(
        c.isalnum() or c in "_-" for c in run_id
    ):
        raise ValueError("run_id must contain only letters, digits, underscore, hyphen")
    models = audit()["models"]
    if model not in models:
        raise ValueError(f"Unknown ordinary T3 model {model!r}")
    if generate:
        args = ["--model", model, "--run-id", run_id,
                "--study", study, "--wolframscript", wolframscript]
        if canonical.main(args) != 0:
            raise RuntimeError("Canonical full-general symbolic generation failed")

    folder = canonical.DEFAULT_DATA_ROOT / run_id / canonical.model_tag(model)
    registry = canonical.load_json(scalar.REGISTRY)
    manifest = load_manifest(scalar.MANIFEST)
    source_hashes = canonical.fingerprint()
    entry = registry["models"][model]
    kind = quadratic_kind(manifest, model)
    validation = canonical.validate_cached(model, folder, entry, kind, source_hashes)
    if validation["matching_performed"] is not False:
        raise ValueError("Unexpected claim of general threshold matching")
    if validation["physical_rgbeta_matchete_mapping_certified"] is not False:
        raise ValueError("Unexpected physical coupling-basis certification")
    uv_path = folder / "uv_rgbeta.json"
    eft_path = folder / "eft1_rgbeta.json"
    uv = _load(uv_path)
    eft = _load(eft_path)
    expected_y = canonical.expected_yukawa(model)
    actual_y = set(uv["report_betas"]) & {"y1Conj", "y2Conj"}
    if actual_y != ({expected_y} if expected_y else set()):
        raise ValueError("Canonical UV extra-Yukawa inventory differs from full Lagrangian")
    if set(eft["report_betas"]) & {"y1", "y2", "MF", "y1Conj", "y2Conj"}:
        raise ValueError("EFT1 export still contains integrated-out fermion couplings")
    info = {
        "source": "canonical_full_general_symbolic",
        "model": model,
        "run_id": run_id,
        "extra_yukawa": expected_y,
        "additional_quartics": validation["extra_quartics"],
        "additional_bilinears": validation["extra_bilinears"],
        "mass_basis": validation["mass_basis"],
        "lagrangian_inventory_checked": True,
        "rge_inventory_checked": True,
        "physical_matchete_rgbeta_mapping_certified": False,
        "physical_neutrino_mass_prediction": False,
        "uv_json": str(uv_path),
        "eft1_json": str(eft_path),
        "uv_sha256": hashlib.sha256(uv_path.read_bytes()).hexdigest(),
        "eft1_sha256": hashlib.sha256(eft_path.read_bytes()).hexdigest(),
    }
    return uv, eft, info


def run(uv_payload: Mapping[str, Any], eft1_payload: Mapping[str, Any], *,
        mu_uv: float = 1e7, mu_f: float = 1e5, mu_s: float = 1e3,
        overrides: Mapping[str, Any] | None = None,
        rtol: float = 1e-7, atol: float = 1e-10) -> dict[str, Any]:
    """Use the existing UV and EFT1 ODE runners; return running parameter snapshots."""
    if not (np.isfinite(mu_uv) and np.isfinite(mu_f) and np.isfinite(mu_s)
            and mu_uv > mu_f > mu_s > 0):
        raise ValueError("Require mu_uv > mu_f > mu_s > 0")
    model, state = initial_state(uv_payload, mu_uv=mu_uv, mu_f=mu_f,
                                 mu_s=mu_s, overrides=overrides)
    uv = run_uv_segment(state, uv_payload, mu_f, rtol=rtol, atol=atol)
    threshold = build_intermediate_scalar_boundary(uv.final_state)
    _, eftlayout = pack_intermediate_scalar_state(threshold)
    eftkeys = set(eft1_payload.get("report_betas", {}))
    if set(eftlayout.names) != eftkeys:
        raise ValueError(f"EFT1 state/export mismatch: missing={sorted(eftkeys-set(eftlayout.names))}, "
                         f"extra={sorted(set(eftlayout.names)-eftkeys)}")
    eft = run_intermediate_scalar_segment(threshold, eft1_payload, mu_s,
                                          rtol=rtol, atol=atol)
    return {"model": model, "status": "NumericalRenormalisableRunningOnly",
            "general_matching_performed": False,
            "scales_gev": {"UV": mu_uv, "F": mu_f, "S": mu_s},
            "uv_solver": {"success": uv.solver_success, "nfev": uv.nfev},
            "eft1_solver": {"success": eft.solver_success, "nfev": eft.nfev},
            "uv_initial": _state_values(uv.initial_state),
            "uv_at_fermion_threshold": _state_values(uv.final_state),
            "eft1_at_fermion_threshold": _state_values(threshold),
            "eft1_at_scalar_threshold": _state_values(eft.final_state)}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--canonical-run-id", help="Validated full-general symbolic run ID")
    source.add_argument("--uv-json", type=Path, help="Direct UV candidate JSON (not canonically validated)")
    p.add_argument("--eft1-json", type=Path, help="Required with --uv-json")
    p.add_argument("--model", help="Model name, required with --canonical-run-id (e.g. A_-2)")
    p.add_argument("--generate-symbolic", action="store_true",
                   help="Generate full-general symbolic export before numerical running")
    p.add_argument("--wolframscript", default="wolframscript")
    p.add_argument("--study", default="general_numerical")
    p.add_argument("--output-json", type=Path, required=True)
    p.add_argument("--initial-json", type=Path, help="Optional overrides for initial couplings")
    p.add_argument("--mu-uv", type=float, default=1e7)
    p.add_argument("--mu-f", type=float, default=1e5)
    p.add_argument("--mu-s", type=float, default=1e3)
    a = p.parse_args(argv)
    overrides = _load(a.initial_json) if a.initial_json else None
    if a.canonical_run_id:
        if not a.model:
            p.error("--model is required with --canonical-run-id")
        uv_payload, eft_payload, provenance = load_canonical_general_exports(
            a.model, a.canonical_run_id, generate=a.generate_symbolic,
            wolframscript=a.wolframscript, study=a.study)
    else:
        if a.eft1_json is None:
            p.error("--eft1-json is required with --uv-json")
        if a.generate_symbolic:
            p.error("--generate-symbolic requires --canonical-run-id")
        uv_payload, eft_payload = _load(a.uv_json), _load(a.eft1_json)
        provenance = {"source": "unvalidated_direct_rgbeta_exports",
                      "lagrangian_inventory_checked": False,
                      "rge_inventory_checked": False,
                      "physical_matchete_rgbeta_mapping_certified": False,
                      "physical_neutrino_mass_prediction": False}
    result = run(uv_payload, eft_payload, mu_uv=a.mu_uv,
                 mu_f=a.mu_f, mu_s=a.mu_s, overrides=overrides)
    result["symbolic_provenance"] = provenance
    result["physical_mass_parameters_available"] = False
    result["basis_mapping_required_for_physical_interpretation"] = True
    a.output_json.parent.mkdir(parents=True, exist_ok=True)
    a.output_json.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n",
                             encoding="utf-8")
    print(f"{result['model']}: UV and EFT1 numerical RGEs completed; {a.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
