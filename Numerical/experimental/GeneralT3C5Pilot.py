"""Experimental *hard-boundary* B/C(-1) C5 response to Matchete quartics.

This is NOT the hierarchical renormalized C5 or an RGBeta-to-Matchete
coupling map. It lifts the coefficient of a single symmetric y_i*y_i
flavour structure in the verified Matchete formal derivative export.
No UV poles may be discarded: the pilot refuses an unresolved pole.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any

import numpy as np
import sympy as sp

from RGE.matching.MatcheteC5Parsing import parse_matchete_c5

ACTIVE = {"lambdaT3ExtraQ020020I1": "y2", "lambdaT3ExtraQ020200I1": "y1"}
MODELS = ("B_m1", "C_m1")
_NAME = re.compile(r"lambdaT3ExtraQ\d+I\d+")


@dataclass(frozen=True)
class ResponseKernel:
    coupling: str
    conjugated: bool
    yukawa: str
    expression: sp.Expr
    raw_expression: str


def _atom_info(atom: str) -> tuple[str, bool]:
    names = _NAME.findall(atom)
    if len(names) != 1:
        raise ValueError(f"Expected one extra coupling in atom: {atom}")
    conjugated = atom.startswith(("Bar[", "Conjugate["))
    if not conjugated and not atom.startswith("Coupling["):
        raise ValueError(f"Unrecognised Mathematica coupling atom: {atom}")
    return names[0], conjugated


# Matchete uses a renormalisation mass squared.  It is set to MF^2 for
# each heavy mass eigenstate; finite logarithms of physical mass RATIOS remain.
_MATCHING_SCALE_SQUARED = {"μbar2", "\\[Mu]bar2", "mubar2"}
_POLE_SYMBOLS = {"ϵ", "\\[Epsilon]", "Epsilon", "epsilon", "eps", "ε"}
_SCALAR_MASS_SYMBOLS = {"MS1", "MS2", "mS1Sq", "mS2Sq", "MS"}
_ALLOWED_KERNEL_SYMBOLS = {"MF", "hbar"} | _SCALAR_MASS_SYMBOLS


def _safe_kernel(expr_text: str, yukawa: str) -> sp.Expr:
    """Extract a finite one-generation response, preserving physical logs.

    Matchete's μbar^2 is a renormalisation scale, not a physical scalar mass.
    At μbar^2=MF^2, Log[MF^2/MS1^2] is a genuine finite threshold
    logarithm and MUST NOT be replaced by zero.  Any remaining symbols,
    unvalidated poles or unknown renormalisation scales cause a hard failure.
    """
    expr = parse_matchete_c5(expr_text)
    own = sp.Symbol(yukawa)
    other = sp.Symbol("y1" if yukawa == "y2" else "y2")
    product = sp.conjugate(own) ** 2
    if expr.has(other, sp.conjugate(other)):
        raise ValueError("Coefficient contains the other Yukawa species")
    kernel = sp.cancel(expr / product)
    if sp.simplify(kernel * product - expr) != 0:
        raise ValueError(f"Cannot factor conjugate({yukawa})^2 from coefficient")
    if kernel.has(own, sp.conjugate(own)):
        raise ValueError("Yukawa still present after response factorisation")

    matching_scales = [s for s in kernel.free_symbols
                       if str(s) in _MATCHING_SCALE_SQUARED]
    if len(matching_scales) > 1:
        raise ValueError("Ambiguous matching-scale symbols: "
                         + ", ".join(sorted(map(str, matching_scales))))
    if matching_scales:
        kernel = sp.simplify(kernel.subs(matching_scales[0], sp.Symbol("MF") ** 2))

    pole_symbols = {str(s) for s in kernel.free_symbols} & _POLE_SYMBOLS
    if pole_symbols:
        raise ValueError(f"Unresolved regulator pole ({sorted(pole_symbols)}); "
                         "no subtraction assumed")

    unknown = {str(s) for s in kernel.free_symbols} - _ALLOWED_KERNEL_SYMBOLS
    if unknown:
        remaining_logs = sorted(str(node) for node in sp.preorder_traversal(kernel)
                                if node.func == sp.log)
        raise ValueError(
            "Unsupported symbols in matching kernel " + str(sorted(unknown)) +
            "; remaining logarithms=" + str(remaining_logs) +
            ". Supply actual scalar masses only for supported physical "
            "mass-ratio logs; never delete logs or UV poles arbitrarily."
        )

    # Reject dimensional logarithms like log(MF) even if their symbol is
    # otherwise supported.  A log argument must have mass dimension zero.
    for node in sp.preorder_traversal(kernel):
        if node.func == sp.log:
            d = _mass_dimension(node.args[0])
            if d != 0:
                raise ValueError("Logarithm argument is not demonstrably "
                                 f"dimensionless: {node.args[0]} (dimension={d})")
    if kernel == 0:
        raise ValueError("Expected nonzero hard response kernel")
    return sp.factor(kernel)


def _mass_dimension(expr: sp.Expr) -> sp.Expr | None:
    """Conservative dimensional analysis of a logarithm argument."""
    if expr.is_Number:
        return sp.S.Zero
    if expr.is_Symbol:
        name = str(expr)
        return ({"MF": 1, "MS1": 1, "MS2": 1, "MS": 1,
                 "mS1Sq": 2, "mS2Sq": 2, "hbar": 0}).get(name)
    if expr.is_Mul:
        dims = [_mass_dimension(arg) for arg in expr.args]
        return None if None in dims else sp.Add(*dims)
    if expr.is_Pow:
        dimension = _mass_dimension(expr.base)
        return None if dimension is None or not expr.exp.is_number else sp.simplify(dimension * expr.exp)
    if expr.is_Add:
        dims = [_mass_dimension(arg) for arg in expr.args]
        return dims[0] if dims and None not in dims and len(set(dims)) == 1 else None
    return None


def _mass_substitutions(parameters: dict[str, Any]) -> dict[sp.Symbol, float]:
    """Positive physical masses in GeV, with squared forms derived safely."""
    result: dict[sp.Symbol, float] = {}
    for label in ("MS1", "MS2", "MS"):
        if label not in parameters:
            continue
        value = float(parameters[label])
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{label} must be a finite positive mass in GeV")
        result[sp.Symbol(label)] = value
    for label in ("mS1Sq", "mS2Sq"):
        number = label[2]  # 1 or 2
        mass = result.get(sp.Symbol(f"MS{number}"))
        supplied = parameters.get(label)
        if supplied is not None:
            val = float(supplied)
            if not np.isfinite(val) or val <= 0.0:
                raise ValueError(f"{label} must be positive in GeV^2")
            if mass is not None and not np.isclose(val, mass * mass, rtol=1e-9):
                raise ValueError(f"Inconsistent {label} and MS{number}")
            result[sp.Symbol(label)] = val
        elif mass is not None:
            result[sp.Symbol(label)] = mass * mass
    return result


def load_kernels(atom_json: Path, model: str) -> tuple[ResponseKernel, ...]:
    if model not in MODELS:
        raise ValueError("Numerical pilot is only defined for ordinary B/C alpha=-1")
    obj = json.loads(Path(atom_json).read_text(encoding="utf-8-sig"))
    if obj.get("model") != model:
        raise ValueError("C5 response artifact model mismatch")
    for key in (
        "response_reconstruction_residual_proven_zero",
        "response_at_zero_proven_zero",
        "polynomial_in_extra_atoms",
    ):
        if obj.get(key) is not True:
            raise ValueError(f"Unproven formal response property: {key}")
    atoms = obj.get("atoms")
    if not isinstance(atoms, list) or len(atoms) != 2:
        raise ValueError("Expected exactly two C5-active derivative atoms")
    out = []
    seen = set()
    for item in atoms:
        coupling, conjugated = _atom_info(item["coupling_atom"])
        if coupling not in ACTIVE or coupling in seen:
            raise ValueError("Unexpected or duplicated C5-active coupling")
        if item.get("coefficient_formally_zero") is True or item.get("depends_on_other_extra_atoms") is True:
            raise ValueError("Invalid non-independent derivative coefficient")
        yukawa = ACTIVE[coupling]
        raw = item["linear_coefficient_inputform"]
        out.append(ResponseKernel(coupling, conjugated, yukawa, _safe_kernel(raw, yukawa), raw))
        seen.add(coupling)
    if seen != set(ACTIVE):
        raise ValueError("Missing an active C5 response coupling")
    return tuple(sorted(out, key=lambda k: k.coupling))


def _decode_matrix(real: Any, imag: Any, n_heavy: int, name: str) -> np.ndarray:
    r = np.asarray(real, dtype=float)
    i = np.asarray(imag, dtype=float)
    if r.shape != (3, n_heavy) or i.shape != r.shape:
        raise ValueError(f"{name} must be a 3 x {n_heavy} complex matrix")
    m = r + 1j * i
    if not np.isfinite(m).all():
        raise ValueError(f"{name} contains nonfinite values")
    return m


def _decode_complex(value: Any, name: str) -> complex:
    if not isinstance(value, dict) or set(value) != {"real", "imag"}:
        raise ValueError(f"{name} must contain real and imag keys")
    v = complex(float(value["real"]), float(value["imag"]))
    if not np.isfinite([v.real, v.imag]).all():
        raise ValueError(f"{name} not finite")
    return v


def evaluate_response(kernels: tuple[ResponseKernel, ...], parameters: dict[str, Any]) -> dict[str, Any]:
    masses = np.asarray(parameters["MF"], dtype=float)
    if masses.ndim != 1 or len(masses) == 0 or not np.isfinite(masses).all() or np.any(masses <= 0):
        raise ValueError("MF must contain finite positive masses in GeV")
    y1 = _decode_matrix(parameters["y1_real"], parameters["y1_imag"], len(masses), "y1")
    y2 = _decode_matrix(parameters["y2_real"], parameters["y2_imag"], len(masses), "y2")
    extras = parameters.get("matchete_extra_couplings")
    if not isinstance(extras, dict) or set(extras) != set(ACTIVE):
        raise ValueError("Provide exactly the two Matchete extra couplings (complex real/imag)")
    hbar = float(parameters.get("hbar", 1.0 / (16.0 * np.pi ** 2)))
    if hbar <= 0 or not np.isfinite(hbar):
        raise ValueError("Invalid hbar")
    masses_subs = _mass_substitutions(parameters)
    contributions = {}
    total = np.zeros((3, 3), dtype=complex)
    for item in kernels:
        lam = _decode_complex(extras[item.coupling], item.coupling)
        if item.conjugated:
            lam = np.conjugate(lam)
        ys = y1 if item.yukawa == "y1" else y2
        mat = np.zeros((3, 3), dtype=complex)
        for r, mass in enumerate(masses):
            evaluated = item.expression.subs({sp.Symbol("MF"): mass, sp.Symbol("hbar"): hbar, **masses_subs})
            if evaluated.free_symbols:
                raise ValueError(f"Unresolved kernel symbols for {item.coupling}: "
                                 f"{sorted(map(str, evaluated.free_symbols))}; "
                                 "provide physical MS1/MS2 masses in the input JSON")
            pref = complex(sp.N(evaluated, 18)) * lam
            # Identical Yukawa species: symmetric matrix by construction.
            v = np.conjugate(ys[:, r])
            mat += pref * np.outer(v, v)
        if not np.isfinite(mat).all():
            raise ValueError("Nonfinite C5 response")
        total += mat
        contributions[item.coupling] = _matrix_json(mat)
    if not np.allclose(total, total.T, rtol=1e-12, atol=1e-20):
        raise ValueError("Result not symmetric")
    return {
        "status": "ExperimentalHardBoundaryC5",
        "matrix_dimension": 3,
        "units": "GeV^-1",
        "contributions": contributions,
        "delta_c5": _matrix_json(total),
        "matching_scale_convention": "μbar^2 = MF_r^2 separately for each heavy mass eigenstate",
        "finite_mass_ratio_logs_preserved": True,
        "scalar_masses_gev": {key: float(parameters[key]) for key in ("MS1", "MS2", "MS") if key in parameters},
        "warning": "Experimental fixed-order Matchete response at the selected MF matching scale, retaining physical mass-ratio logs. NOT renormalised hierarchical final C5, not RG evolved, and the Matchete-to-RGBeta map has NOT been verified.",
    }


def _matrix_json(m: np.ndarray) -> dict[str, Any]:
    return {"real": m.real.tolist(), "imag": m.imag.tolist(),
            "frobenius_norm": float(np.linalg.norm(m, ord="fro"))}
