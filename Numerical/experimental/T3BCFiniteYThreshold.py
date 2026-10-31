"""B(-1)/C(-1) Y-threshold subtraction CONTRACT (not a physical certificate).

This module does not parse raw Matchete InputForm or manufacture counterterms.
It evaluates, only after a caller supplies *separately projected* ordered
Weinberg amplitudes with explicit Laurent pole and finite pieces,

    Delta A_Y = [A_EFT1^1PI + A_EFT1^CT
                 - A_EFT2^1PI - A_EFT2^CT
                 + A_ext/EOM/parameter]_{finite},
    Delta C5_Y = Delta A_Y + Delta A_Y^T.

The EFT2 loop contribution must EXCLUDE the unknown local threshold insertion
being solved for, while including existing lower-order EFT operators.  Inputs
must already contain hbar = 1/(16*pi^2), the same flavor basis, the same IR
regulator and the same MSbar convention.  Amplitude projection, MSbar
counterterm construction and EFT1/EFT2 source equivalence are NOT implemented.

Conventions: repository Lagrangian/RunMatching.wl, WeinbergOperatorConvention[],
uses L5 = +1/2 C5^{ij} O_ij + h.c. and C5 = A5 + A5^T.  No Majorana mass
conversion/sign is inferred from an unverified GammaCC spinor projection.

A returned candidate is a CONDITIONAL NUMERICAL SUBTRACTION, not a renormalized
physical Wilson coefficient. The 'physical_c5' output is always null until
independent B/C evidence and matching normalization checks are implemented.

Intended repository path: Numerical/experimental/T3BCFiniteYThreshold.py
Run: python -m Numerical.experimental.T3BCFiniteYThreshold --self-test
     python -m Numerical.experimental.T3BCFiniteYThreshold --input your.json
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

MODELS = frozenset({"B_-1", "C_-1"})
CHECKS = (
    "original_eomeft_equivalent_to_scalar_source",
    "ordered_gamma_cc_chirality_and_sign_verified",
    "full_and_eft_1pi_amplitudes_use_same_ir_prescription",
    "msbar_counterterms_and_evanescents_verified",
    "external_legs_eom_and_parameter_shifts_included",
    "fermion_boundary_and_eft1_running_accounted_for",
    "scalar_mixing_and_quartic_normalizations_verified",
)


class MatchingInputError(ValueError):
    """Invalid or incomplete subtraction inputs: no physical result is returned."""


def _matrix(value: Any, label: str) -> np.ndarray:
    """A 3x3 matrix encoded as real nested arrays or {real,imag} arrays."""
    try:
        if isinstance(value, dict):
            if "real" not in value or set(value) - {"real", "imag"}:
                raise ValueError("matrix object must contain real and optional imag")
            re = np.asarray(value["real"], dtype=float)
            im = np.asarray(value.get("imag", np.zeros((3, 3))), dtype=float)
            if re.shape != im.shape:
                raise ValueError("real/imag shapes differ")
            arr = re + 1j * im
        else:
            arr = np.asarray(value, dtype=float).astype(complex)
    except (TypeError, ValueError) as exc:
        raise MatchingInputError(f"{label}: invalid numeric matrix: {exc}") from exc
    if arr.shape != (3, 3) or not np.isfinite(arr).all():
        raise MatchingInputError(f"{label}: expected finite complex 3x3 matrix")
    return arr


@dataclass(frozen=True)
class Laurent:
    pole: np.ndarray
    finite: np.ndarray

    @classmethod
    def from_dict(cls, data: Any, name: str) -> "Laurent":
        if not isinstance(data, dict) or set(data) != {"pole", "finite"}:
            raise MatchingInputError(f"{name}: require exactly 'pole' and 'finite'")
        return cls(_matrix(data["pole"], f"{name}.pole"),
                   _matrix(data["finite"], f"{name}.finite"))


def _symmetric(m: np.ndarray, label: str, rtol: float) -> None:
    scale = max(np.max(np.abs(m)), 1e-30)
    if np.max(np.abs(m - m.T)) > rtol * scale:
        raise MatchingInputError(f"{label} is not complex SYMMETRIC (transpose, not adjoint)")


def _encode(m: np.ndarray) -> dict[str, list[list[float]]]:
    return {"real": m.real.tolist(), "imag": m.imag.tolist()}


def compute_finite_y_subtraction(data: dict, *, rtol: float = 1e-10) -> dict:
    """Conditionally subtract projected Laurent matrices without false certification.

    Required JSON keys: model, mu_GeV, scheme='MSbar', above_c5,
    tree_shift_c5, full_1pi, low_1pi, full_counterterm, low_counterterm,
    external_eom_parameter_shift, and evidence_checks. Every Laurent input
    has pole and finite matrices. `above_c5` and tree_shift_c5 are already
    normalized symmetric C5; others are ORDERED A5 contributions.

    Pole cancellation is necessary but insufficient for physical matching.
    Even if evidence_checks are all true, this interface does not authenticate
    those claims, so physical_matching_certified remains false.
    """
    if not isinstance(data, dict):
        raise MatchingInputError("input must be a JSON object")
    required = {"model", "mu_GeV", "scheme", "above_c5", "tree_shift_c5",
                "full_1pi", "low_1pi", "full_counterterm", "low_counterterm",
                "external_eom_parameter_shift", "evidence_checks"}
    missing = required - data.keys()
    if missing:
        raise MatchingInputError("required input absent: " + ", ".join(sorted(missing)))
    model = data["model"]
    if model not in MODELS:
        raise MatchingInputError("only B_-1 and C_-1 are supported")
    mu = data["mu_GeV"]
    if not isinstance(mu, (float, int)) or not np.isfinite(mu) or mu <= 0:
        raise MatchingInputError("mu_GeV must be positive finite GeV")
    if data["scheme"] != "MSbar":
        raise MatchingInputError("only an explicitly MSbar-labelled input is supported")
    if not np.isfinite(rtol) or rtol <= 0:
        raise MatchingInputError("rtol must be positive finite")

    c_above = _matrix(data["above_c5"], "above_c5")
    c_tree = _matrix(data["tree_shift_c5"], "tree_shift_c5")
    _symmetric(c_above, "above_c5", rtol)
    _symmetric(c_tree, "tree_shift_c5", rtol)
    items = {name: Laurent.from_dict(data[name], name) for name in (
        "full_1pi", "low_1pi", "full_counterterm", "low_counterterm",
        "external_eom_parameter_shift")}
    u, l = items["full_1pi"], items["low_1pi"]
    cu, cl = items["full_counterterm"], items["low_counterterm"]
    e = items["external_eom_parameter_shift"]

    pole = u.pole + cu.pole - l.pole - cl.pole + e.pole
    finite = u.finite + cu.finite - l.finite - cl.finite + e.finite
    # Normalize against input Laurent residues, NOT against an arbitrary GeV=1.
    norm = max(np.max(np.abs(z)) for z in (
        u.pole, cu.pole, l.pole, cl.pole, e.pole, finite, np.zeros((3, 3))))
    scale = max(norm, 1e-30)
    max_residual = float(np.max(np.abs(pole)))
    if max_residual > rtol * scale:
        raise MatchingInputError(
            f"uncancelled 1/epsilon residue {max_residual:.12g} "
            f"> {rtol * scale:.12g}: cannot extract a finite threshold")

    proof = data["evidence_checks"]
    if not isinstance(proof, dict) or any(
        not isinstance(proof.get(k), bool) for k in CHECKS
    ):
        raise MatchingInputError("evidence_checks must contain every named boolean check")
    delta_c5 = finite + finite.T  # Explicit repository ordered-A5 convention.
    formal_below = c_above + c_tree + delta_c5
    _symmetric(formal_below, "formal_below_c5", rtol)
    blockers = [k for k in CHECKS if not proof[k]]
    return {
        "format": "t3-bc-y-finite-subtraction-contract-v1",
        "model": model, "mu_GeV": float(mu), "scheme_claimed": "MSbar",
        "operator_convention": "L5=+1/2 C5^(ij) O_ij+h.c.; C5=A5+A5^T",
        "source_projection": "ordered_A5_supplied_by_caller_not_reconstructed_from_raw_Matchete",
        "pole_residue_max_abs": max_residual,
        "numerical_pole_cancellation": True,
        "formal_delta_c5": _encode(delta_c5),
        "formal_c5_below_y": _encode(formal_below),
        "unverified_prerequisites": blockers,
        "evidence_checks_are_caller_assertions_only": True,
        "physical_c5": None,
        "physical_matching_certified": False,
        "y_llxx_inherited_to_x_certified": False,
        "note": "No physical Y matching is certified from supplied data alone.",
    }


def _self_test() -> None:
    z = np.zeros((3, 3))
    eye = np.eye(3)
    def m(v: np.ndarray) -> dict:
        return _encode(v.astype(complex))
    def lp(p: np.ndarray, f: np.ndarray) -> dict:
        return {"pole": m(p), "finite": m(f)}
    above = np.diag([1., 2., 3.]).astype(complex)
    ordered = np.array([[1 + 1j, 2, 0], [3, 0, 0], [0, 0, -1]], complex)
    d = {"model": "B_-1", "mu_GeV": 1187.806, "scheme": "MSbar",
         "above_c5": m(above), "tree_shift_c5": m(z),
         "full_1pi": lp(eye, ordered), "low_1pi": lp(z, z),
         "full_counterterm": lp(-eye, z), "low_counterterm": lp(z, z),
         "external_eom_parameter_shift": lp(z, z),
         "evidence_checks": {k: True for k in CHECKS}}
    r = compute_finite_y_subtraction(d)
    expected = above + ordered + ordered.T
    got = np.asarray(r["formal_c5_below_y"]["real"]) + 1j * np.asarray(r["formal_c5_below_y"]["imag"])
    assert np.allclose(got, expected)
    assert r["physical_c5"] is None and r["physical_matching_certified"] is False
    assert not np.allclose(got, np.conjugate(got.T))  # No accidental Hermitian conjugation.
    bad = dict(d, full_counterterm=lp(z, z))
    try:
        compute_finite_y_subtraction(bad)
    except MatchingInputError as err:
        assert "uncancelled" in str(err)
    else:
        raise AssertionError("Residual pole must be rejected")
    for mod in ("B_-1", "C_-1"):
        trial = dict(d, model=mod)
        assert compute_finite_y_subtraction(trial)["model"] == mod
    print("PASS: B/C conditional subtraction, ordered symmetrization, pole guard, uncertified physical flag")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--input", type=Path, help="projected B/C Laurent-amplitude JSON")
    p.add_argument("--output", type=Path, help="optional JSON report")
    args = p.parse_args()
    if args.self_test:
        _self_test()
        return 0
    if not args.input:
        p.error("provide --self-test or --input")
    try:
        data = json.loads(args.input.read_text(encoding="utf-8"))
        report = compute_finite_y_subtraction(data)
    except (ValueError, OSError) as exc:
        p.exit(2, f"NOT CERTIFIED: {exc}\n")
    out = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(out, encoding="utf-8")
    else:
        print(out, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
