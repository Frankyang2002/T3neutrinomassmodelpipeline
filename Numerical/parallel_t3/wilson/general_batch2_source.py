"""Source-preserving Batch-2 first-threshold adapter (not a physical Wilson solver).

This module deliberately refuses to call a bare Matchete one-loop EFT a finite
numerical Wilson coefficient.  Native Matchete, RGBeta, and the flavour basis
have not been connected by a certified coupling/CG dictionary in this project.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from Numerical.parallel_t3.background.general_batch1_api import (
    decode_native, encode_native,
)

SCHEMA = "t3-thesis-full-general-batch2-raw-fermion-v1"
MODELS = ("A_-4", "D_-2")
HEAVY = frozenset(("MF", "y1", "y2", "y1Conj", "y2Conj"))
D_QUARTICS = frozenset((
    "lambdaT3GeneralQ000004I1", "lambdaT3GeneralQ000013I1",
    "lambdaT3GeneralQ001102I1", "lambdaT3GeneralQ020101I1",
    "lambdaT3GeneralQ110002I1",
))
BLOCKERS = (
    "Source-normalised RGBeta to Matchete coupling/CG and heavy-flavour-mass dictionary",
    "Complete one-loop EFT1 operator projection, including all generated operators",
    "MSbar UV-hard minus EFT1-soft subtraction/counterterm with pole cancellation",
    "Full EFT1 one-loop Wilson anomalous dimensions and operator mixing",
    "Full scalar-mass-eigenfield tensor mapping for subsequent scalar threshold",
)


class MatchingIncomplete(RuntimeError):
    """A required physical step cannot be inferred from a symbolic match."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise MatchingIncomplete("Batch2: " + message)


def sha(path: Path) -> str:
    path = Path(path)
    require(path.is_file(), "missing source file: " + str(path))
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    require(isinstance(data, dict), "expected JSON object: " + str(path))
    return data


def write_json(path: Path, data: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    partial.write_text(json.dumps(data, indent=2, sort_keys=True,
                                  allow_nan=False) + "\n", encoding="utf-8")
    partial.replace(path)


def native_numbers(values: Mapping[str, Any], model: str, *, at_f: bool = True) -> dict:
    """Check every native beta coordinate and retain all permitted interactions.

    The 3x3 MF rotations are exposed without guessing their chiral assignment.
    They are *not* a Matchete/RGBeta physical dictionary.
    """
    require(model in MODELS, "unsupported model")
    names = set(values)
    must = {"MF", "y1", "y2", "mS1Sq", "mS2Sq", "lambdaT3"}
    require(must <= names, "missing native input(s): " + str(sorted(must - names)))
    expected_extra = {"y2Conj", "muT3GeneralS2Sq"} | D_QUARTICS if model == "D_-2" else set()
    found_extra = {n for n in names if n in {"y1Conj", "y2Conj", "muT3GeneralS1Sq",
                                               "muT3GeneralS2Sq", "muT3GeneralS1S2"}
                   or n.startswith("lambdaT3GeneralQ")}
    require(found_extra == expected_extra,
            "incorrect model's symmetry-allowed live inventory: missing="
            + str(sorted(expected_extra - found_extra)) + ", extra="
            + str(sorted(found_extra - expected_extra)))
    native = {key: np.asarray(decode_native(val), dtype=complex)
              for key, val in values.items()}
    require(all(np.isfinite(x).all() for x in native.values()), "nonfinite native input")
    require(all(native[n].shape == (3, 3) for n in ("MF", "y1", "y2")),
            "heavy flavour matrices must all be 3x3")
    if "y2Conj" in native:
        require(native["y2Conj"].shape == (3, 3), "D_-2 y2Conj not 3x3")
    masses = np.linalg.svd(native["MF"], compute_uv=False)
    require(np.min(masses) > 0, "MF must be invertible")
    u, singular, vh = np.linalg.svd(native["MF"])
    residual = np.max(abs(u.conj().T @ native["MF"] @ vh.conj().T - np.diag(singular)))
    require(residual <= 1e-9 * max(1., float(max(singular))),
            "heavy-mass singular-value decomposition did not reconstruct")
    mu2 = complex(native["muT3GeneralS2Sq"]) if model == "D_-2" else 0j
    # The existing full-general D_-2 Matchete source supports a *real*
    # holomorphic self-bilinear, but not an independent complex phase.
    require(abs(mu2.imag) <= 1e-11 * max(1., abs(mu2.real)),
            "D_-2 native holomorphic mass has a physical unimplemented complex phase")
    m1 = complex(native["mS1Sq"]); m2 = complex(native["mS2Sq"])
    require(abs(m1.imag) < 1e-7 and abs(m2.imag) < 1e-7,
            "diagonal scalar masses must be real")
    scalar_eigenvalues = ([m1.real] * 2 if model == "A_-4" else [m1.real] * 6)
    if model == "A_-4":
        scalar_eigenvalues += [m2.real] * 6
    else:
        scalar_eigenvalues += [m2.real - 2 * abs(mu2), m2.real + 2 * abs(mu2)]
    require(min(scalar_eigenvalues) > 0, "tachyonic scalar input at fermion threshold")
    require(len(scalar_eigenvalues) == 8, "wrong number of real scalar components")
    return {
        "model": model,
        "scale_GeV": 1e5 if at_f else None,
        "native_parameter_values": {k: encode_native(v) for k, v in native.items()},
        "native_parameter_names": sorted(native),
        "heavy_fermion_singular_masses_GeV": sorted(float(x) for x in singular),
        "fermion_left_svd_matrix": encode_native(u),
        "fermion_right_svd_matrix": encode_native(vh.conj().T),
        "heavy_svd_is_only_a_diagnostic_not_Matchete_rotation": True,
        "candidate_scalar_quadratic_mass_eigenvalues_GeV2": sorted(scalar_eigenvalues),
        "neutral_scalar_self_bilinear_GeV2": encode_native(mu2),
        "all_f_fermion_modes_still_present": True,
        "scalar_multiplets_light_at_f": ["S1", "S2"],
        "matchete_to_rgbeta_dictionary_applied": False,
        "source_exact_one_loop_wilson_numbers_produced": False,
    }


def verify_native_wolfram(folder: Path, model: str) -> dict:
    """Reuse existing RunT3GeneralFThresholdBatch's source/CG guard."""
    from scripts.RunT3GeneralFThresholdBatch import verify_native_tree
    result = verify_native_tree(folder, model, 1)
    report = load_json(folder / "native_status.json")
    require(report.get("raw_one_loop_F_matching_included") is True,
            "Matchete did not request the complete one-loop F threshold")
    data_dir = folder / "data"
    prefix = "eft1_after_F"
    files = {
        "tree": data_dir / f"{prefix}_tree_inputform.txt",
        "one_loop": data_dir / f"{prefix}_one_loop_inputform.txt",
        "compact": data_dir / f"{prefix}_compact_inputform.txt",
        "rge_export": data_dir / f"{prefix}_rge_export.json",
        "tree_wilson_seed": next(iter(data_dir.glob("*wilson_seed.json")), None),
    }
    require(all(p is not None and p.is_file() and p.stat().st_size > 0
                for p in files.values()), "incomplete native Matchete matching outputs")
    rge = load_json(files["rge_export"])
    require(rge.get("IntegratedFields") == ["F"]
            and rge.get("ActiveHeavyFields") == ["S1", "S2"],
            "one-loop matched EFT is not fermion-first with both scalars light")
    text = files["one_loop"].read_text(encoding="utf-8")
    require("hbar" in text, "one-loop matched EFT contains no loop-bookkeeping factor")
    require("Field[NewFermion" not in text,
            "one-loop matched EFT retains integrated fermion field")
    return {
        **result,
        "native_source_files_sha256": {n: sha(p) for n, p in files.items()},
        "native_rge_export": rge,
        "raw_matchete_loop_different_from_zero": text.strip() not in ("0", "0.")
    }


def verify_ledger(report: Mapping[str, Any], raw_path: Path) -> dict:
    """Full operator ledger must account for the exact unmodified bare source.

    The ledger records ALL terms. It does NOT project to only LLSS/LLHH.
    """
    require(report.get("format") == "t3-batch2-complete-loop-operator-ledger-v1",
            "wrong Matchete operator ledger format")
    require(report.get("status") == "UnrenormalisedOneLoopTermsEnumerated",
            "Matchete operator ledger did not finish successfully")
    require(report.get("source_sha256") == sha(raw_path),
            "operator ledger source SHA differs from native one-loop output")
    n = report.get("term_count")
    terms = report.get("terms")
    require(type(n) is int and n > 0 and isinstance(terms, list) and len(terms) == n,
            "one-loop terms are missing or incomplete")
    require(all(type(v) is dict and type(v.get("raw_inputform")) is str
                and v["raw_inputform"] for v in terms),
            "invalid Matchete full operator term entry")
    require(report.get("exact_reconstruction_by_expansion") is True,
            "Matchete full one-loop expression was not reconstructed from every term")
    require(report.get("only_fermions_integrated_out") is True,
            "F-only threshold field criterion failed")
    require(report.get("finite_msbar_wilson_coefficients_computed") is False,
            "bare one-loop ledger incorrectly promoted to finite Wilson coefficients")
    return {
        "term_count": n,
        "one_loop_input_sha256": sha(raw_path),
        "operator_ledger_classification_complete": True,
        "finite_msbar_wilson_coefficients_computed": False,
        "raw_one_loop_operator_sector_counts": report.get("sector_counts", {}),
    }


def batch2_scope() -> dict:
    return {
        "source_interactions": "general symmetry-allowed Matchete native UV",
        "matched_scope": "all heavy fermions simultaneously, all other fields light, through one loop",
        "matching_scale_GeV": 1e5,
        "subtraction_scheme_target": "MSbar: UV one-loop hard - EFT1 soft + counterterms",
        "subtraction_scheme_implemented": False,
        "numerical_complete_Wilson_projection_implemented": False,
        "eft1_all_operator_one_loop_ADM_implemented": False,
        "physical_matching_succeeded": False,
        "missing_physics": list(BLOCKERS),
        "restricted_analytic_triangle_used": False,
    }
