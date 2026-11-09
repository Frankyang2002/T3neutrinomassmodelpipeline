"""Compare saved B_-1/C_-1 one-scale and dynamic sequential C5 at muS.

Read-only physics diagnostic. This script DOES NOT perform matching, RG running,
optimisation, or change any existing result. The best complex rescaling is only
an algebraic flavour-shape diagnostic, NOT a physical renormalisation.

Default inputs (run from project root):
  output/studies/general_scalar_extended/thesis_general_all16/
    thesis16_general_c5_full_v1/{B_m1,C_m1}/result.json
  output/general_sequential_pilot/all16_dynamic_scalar_threshold/
    batch16_v1/{B_m1,C_m1}/result.json

Usage (PowerShell):
  python -X utf8 -m scripts.CompareT3MixedC5 --run-id batch16_v1
  python -X utf8 -m scripts.CompareT3MixedC5 --run-id batch16_v1 `
    --output-json Reports/output/mixed_c5_comparison/batch16_v1.json `
    --print-matrices

Requires NumPy. Does not write or contact GitHub.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODELS = ("B_-1", "C_-1")
FIT_ROOT = (ROOT / "output/studies/general_scalar_extended/thesis_general_all16/"
            "thesis16_general_c5_full_v1")
SEQUENTIAL_ROOT = ROOT / "output/general_sequential_pilot/all16_dynamic_scalar_threshold"


def _slug(model: str) -> str:
    return model.replace("-", "m").replace("+", "p")


def _load(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing input: {path}")
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return obj


def _complex_value(raw: Any, label: str) -> complex:
    if isinstance(raw, dict) and set(raw) == {"real", "imag"}:
        a = complex(float(raw["real"]), float(raw["imag"]))
    elif isinstance(raw, (float, int, complex)):
        a = complex(raw)
    else:
        raise ValueError(f"{label}: malformed scalar encoding")
    if not math.isfinite(a.real) or not math.isfinite(a.imag):
        raise ValueError(f"{label}: nonfinite scalar")
    return a


def _matrix(raw: Any, label: str) -> np.ndarray:
    if isinstance(raw, dict) and set(raw) == {"real", "imag"}:
        re = np.asarray(raw["real"], dtype=float)
        im = np.asarray(raw["imag"], dtype=float)
        if re.shape != (3, 3) or im.shape != (3, 3):
            raise ValueError(f"{label}: incorrect real/imag matrix shape")
        a = re + 1j * im
    else:
        a = np.asarray(raw, dtype=complex)
    if a.shape != (3, 3) or not np.isfinite(a).all():
        raise ValueError(f"{label}: expected finite complex 3x3 matrix")
    symmetry_residual = float(np.linalg.norm(a-a.T) / max(float(np.linalg.norm(a)), 1e-300))
    if symmetry_residual > 1e-8:
        raise ValueError(f"{label}: C5 is not symmetric (relative={symmetry_residual:.3e})")
    return a


def _encode(a: complex | np.ndarray) -> dict[str, Any]:
    z = np.asarray(a, dtype=complex)
    if z.ndim == 0:
        return {"real": float(z.real), "imag": float(z.imag)}
    return {"real": z.real.tolist(), "imag": z.imag.tolist()}


def _reldiff(x: complex, y: complex) -> float:
    return float(abs(x-y)/max(abs(x),abs(y),1e-300))


def _stage_comparisons(old: dict[str, Any], new: dict[str, Any]) -> dict[str, Any]:
    """Check overlapping muS scalar couplings; no equality is presumed."""
    old_stage = old.get("details", {}).get("trajectory", {}).get("eft1_end")
    new_stage = new.get("rge", {}).get("renormalisable_parameters_at_S")
    if not isinstance(old_stage, dict) or not isinstance(new_stage, dict):
        raise ValueError("Missing old/new scalar EFT1 endpoints")
    pairs: dict[str, Any] = {
        "mS1Sq": (old_stage["mS1Sq_gev2"], new_stage["mS1Sq"]),
        "mS2Sq": (old_stage["mS2Sq_gev2"], new_stage["mS2Sq"]),
        "lambdaT3": (old_stage["lambdaT3"], new_stage["lambdaT3"]),
        "muT3GeneralS1S2": (
            old_stage["quadratic_mu2_gev2"],new_stage["muT3GeneralS1S2"]),
    }
    for sector in ("new_quartics", "reused_quartics"):
        for k,v in old_stage.get(sector, {}).items():
            if k not in new_stage:
                raise ValueError(f"Sequential scalar EFT1 endpoint missing old coupling {k}")
            pairs[k] = (v, new_stage[k])
    rows = {}
    for name,(v_old,v_new) in sorted(pairs.items()):
        x=_complex_value(v_old, f"old muS {name}")
        y=_complex_value(v_new, f"new muS {name}")
        rows[name] = {"old": _encode(x), "sequential": _encode(y),
                      "relative_difference": _reldiff(x,y)}
    return {"couplings": rows,
            "maximum_relative_difference": max(r["relative_difference"] for r in rows.values()),
            "same_EFT1_endpoint_to_1percent": all(r["relative_difference"] <= .01
                                                  for r in rows.values())}


def _comparison(model: str, fit_root: Path, sequence_root: Path) -> dict[str, Any]:
    tag = _slug(model)
    old_folder = fit_root / tag
    old_path = old_folder / "result.json"
    new_path = sequence_root / tag / "result.json"
    old = _load(old_path)
    new = _load(new_path)
    if old.get("model") != model or new.get("model") != model:
        raise ValueError(f"{model}: model mismatch in source results")
    if old.get("status") != "ExploratoryGeneralScalarOneScaleConditionalFit":
        raise ValueError(f"{model}: original is not the expected saved one-scale fit")
    if new.get("status") != "ConditionalDynamicFToBothScalarsToSMEFTCompleted":
        raise ValueError(f"{model}: sequential result has unexpected status")

    reference = new.get("reference", {})
    # Same saved fitted boundary; does not assert identical downstream RG running.
    if Path(str(reference.get("fit_root", ""))).resolve() != old_folder.resolve():
        raise ValueError(f"{model}: sequential runner used a DIFFERENT saved fit folder: "
                         f"{reference.get('fit_root')!r}")
    old_chi2 = float(old["chi2_final"])
    new_source_chi2 = float(reference["original_one_scale_fit_chi2"])
    if not math.isclose(old_chi2, new_source_chi2, rel_tol=1e-10, abs_tol=1e-12):
        raise ValueError(f"{model}: saved-fit chi2 does not match sequential provenance")
    old_config = _load(old_folder / "fitted_config.json")
    if old_config.get("representation", {}).get("alpha") != -1:
        raise ValueError(f"{model}: old fit representation is not the mixed model")
    muS = float(old_config["scales"]["mu_scalar_threshold_gev"])
    if not math.isclose(muS, float(new["scales_gev"]["S"]), rel_tol=1e-12):
        raise ValueError(f"{model}: comparison scales differ")

    # As in the existing physics runners, both matrices are in GeV^-1 at muS.
    c_old = _matrix(old["details"]["C5_muS_gev_inv"], "old C5(muS)")
    scalar = new["scalar_matching"]
    c_seq = _matrix(scalar["C5_used_at_S_gev_inverse"], "sequential C5(muS)")
    c_log = _matrix(scalar["direct_log_gev_inverse"], "sequential direct log")
    c_bubble = _matrix(scalar["both_scalar_finite_bubble_gev_inverse"], "sequential bubble")
    n_old = float(np.linalg.norm(c_old))
    n_seq = float(np.linalg.norm(c_seq))
    if n_old == 0 or n_seq == 0:
        raise ValueError(f"{model}: C5 norm is zero; complex rescaling undefined")
    if not np.allclose(c_seq, c_log+c_bubble, rtol=1e-9, atol=1e-28):
        raise ValueError(f"{model}: sequential C5 does not equal log + bubble")

    # Minimise ||C_seq-a*C_old||^2 over complex a.
    a = complex(np.vdot(c_old, c_seq)/np.vdot(c_old,c_old))
    residual = float(np.linalg.norm(c_seq-a*c_old)/n_seq)
    cosine = float(abs(np.vdot(c_old,c_seq))/(n_old*n_seq))
    cosine = min(1.0,max(0.0,cosine))
    stage = _stage_comparisons(old, new)
    result = {
        "model": model,
        "comparison_scope": "saved one-scale conditional C5 versus sequential soft-only C5 at same saved UV fit",
        "old_result": str(old_path.resolve()),
        "sequential_result": str(new_path.resolve()),
        "common_saved_fit_folder_confirmed": True,
        "same_muS_gev": muS,
        "old_fit_chi2": old_chi2,
        "sequential_chi2_at_old_fit_point": float(new["neutrinos"]["oscillation_chi2"]),
        "old_C5_norm_gev_inverse": n_old,
        "sequential_C5_norm_gev_inverse": n_seq,
        "sequential_direct_log_norm_gev_inverse": float(np.linalg.norm(c_log)),
        "sequential_bubble_norm_gev_inverse": float(np.linalg.norm(c_bubble)),
        "norm_ratio_sequential_over_old": n_seq/n_old,
        "optimal_complex_a_in_Cseq_equals_a_times_Cold": _encode(a),
        "optimal_scale_magnitude": abs(a),
        "optimal_scale_phase_degrees": float(np.angle(a, deg=True)),
        "shape_alignment_cosine_magnitude": cosine,
        "shape_residual_after_optimal_complex_rescaling": residual,
        "unscaled_relative_difference_to_old": float(np.linalg.norm(c_seq-c_old)/n_old),
        "muS_scalar_endpoint_comparison": stage,
        "old_C5_muS_gev_inverse": _encode(c_old),
        "sequential_C5_muS_gev_inverse": _encode(c_seq),
        "sequential_direct_log_muS_gev_inverse": _encode(c_log),
        "sequential_bubble_muS_gev_inverse": _encode(c_bubble),
        "interpretation_limits": [
            "Matching prescriptions differ; results do not isolate an RG correction.",
            "Independent hard one-loop matching is unavailable for these sequential points.",
            "Absolute GammaCC/Weyl normalisation and mixed muF/muS scalar-basis identification remain unverified.",
            "Same saved UV fit is confirmed by file path/chi2; overlapping muS endpoint couplings are compared independently.",
        ],
    }
    return result


def _print_matrix(label: str, raw: dict[str, Any]) -> None:
    a = _matrix(raw,label)
    print(f"  {label} [GeV^-1]:")
    for row in a:
        print("    " + "   ".join(f"{x.real:+.6e}{x.imag:+.6e}j" for x in row))


def _print_record(r: dict[str, Any], *, matrices: bool) -> None:
    print(f"\n{r['model']}: old one-scale vs sequential C5 at muS={r['same_muS_gev']:g} GeV")
    print(f"  saved fit chi2                 = {r['old_fit_chi2']:.9g}")
    print(f"  sequential chi2 (same fit)     = {r['sequential_chi2_at_old_fit_point']:.9g}")
    print(f"  ||C5_old||_F                  = {r['old_C5_norm_gev_inverse']:.9e} GeV^-1")
    print(f"  ||C5_sequential||_F           = {r['sequential_C5_norm_gev_inverse']:.9e} GeV^-1")
    print(f"  ||sequential|| / ||old||       = {r['norm_ratio_sequential_over_old']:.9g}")
    a = r['optimal_complex_a_in_Cseq_equals_a_times_Cold']
    print(f"  best complex a                = {a['real']:+.9g}{a['imag']:+.9g}j")
    print(f"  |a|, arg(a)                   = {r['optimal_scale_magnitude']:.9g}, "
          f"{r['optimal_scale_phase_degrees']:.6g} deg")
    print(f"  shape alignment |cos|         = {r['shape_alignment_cosine_magnitude']:.9g}")
    print(f"  residual r after complex a    = {r['shape_residual_after_optimal_complex_rescaling']:.9g}")
    print(f"  raw relative difference       = {r['unscaled_relative_difference_to_old']:.9g}")
    stage = r['muS_scalar_endpoint_comparison']
    print(f"  muS overlapping couplings max relative difference = "
          f"{stage['maximum_relative_difference']:.9g}")
    if not stage['same_EFT1_endpoint_to_1percent']:
        print("  CAUTION: muS EFT1 endpoints differ >1%; this comparison does not "
              "isolate the matching kernel.")
    print("  muS overlapping coupling differences:")
    for name,entry in stage['couplings'].items():
        print(f"    {name:36s} {entry['relative_difference']:.6g}")
    if matrices:
        _print_matrix("old C5",r['old_C5_muS_gev_inverse'])
        _print_matrix("sequential C5",r['sequential_C5_muS_gev_inverse'])
        _print_matrix("sequential log",r['sequential_direct_log_muS_gev_inverse'])
        _print_matrix("sequential bubble",r['sequential_bubble_muS_gev_inverse'])


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-id",default="batch16_v1", help="Existing sequential batch run ID")
    p.add_argument("--fit-root",type=Path,default=FIT_ROOT)
    p.add_argument("--sequential-root",type=Path,default=SEQUENTIAL_ROOT)
    p.add_argument("--output-json",type=Path, help="Optional new diagnostic JSON path")
    p.add_argument("--print-matrices",action="store_true")
    args = p.parse_args(argv)
    if not args.run_id or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in args.run_id):
        p.error("--run-id must contain only letters, digits, underscores or hyphens")
    sequence = args.sequential_root / args.run_id
    records = [_comparison(m,args.fit_root,sequence) for m in MODELS]
    print("C5 diagnostic: identical stored one-scale fit source, both coefficients at muS")
    for record in records:
        _print_record(record,matrices=args.print_matrices)
    print("\nNOTE: small residual r means predominantly overall complex rescaling. "
          "Large r means the two flavour matrices are not proportional.")
    print("Matching completeness/basis certification are NOT established by this test.")
    if args.output_json is not None:
        path=args.output_json
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite report: {path}")
        path.parent.mkdir(parents=True,exist_ok=True)
        doc={"format":"t3-mixed-c5-matrix-comparison-v1","run_id":args.run_id,
             "models":{r['model']:r for r in records}}
        path.write_text(json.dumps(doc,indent=2,allow_nan=False)+"\n",encoding="utf-8")
        print(f"Report: {path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
