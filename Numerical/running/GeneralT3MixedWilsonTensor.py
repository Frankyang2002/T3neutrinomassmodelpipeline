"""Source-faithful B_-1/C_-1 mixed-doublet tree LLSS -> ordered real tensors.

Requires the actual Matchete native CG registry exported during the same
mass-basis F-only matching run. Never infers CG phases or operator symmetry.
This is a first-boundary numeric tensor, NOT canonical C_ijab / EFT1 running.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any, Mapping

import numpy as np

from Numerical.running.GeneralT3OrderedTensor import project_ordered_term
from scripts.RunT3GeneralWilsonBoundary import _factors, _MF, flavour_matrix

MIXED = ("B_-1", "C_-1")
FORMAT = "t3-mixed-native-matchete-ordered-tensor-v1"
_SPURION = re.compile(r"(?P<bar>Bar\[)?Coupling\[(?P<name>T3MEU[12][12]),\s*\{\},\s*0\](?P<close>\])?")
_ALL_COUPLING = re.compile(r"Coupling\[([^,\]]+)")


def sha(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    obj = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return obj


def save_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _is_matchete_sum_separator(c: str, levels: tuple[int, int, int]) -> bool:
    return c in "+-" and levels == (0, 0, 0)


def split_matchete_sum(raw: str) -> list[str]:
    """Expand a flattened InputForm Plus into signed terms without SymPy guesses.

    Matchete prints complete multiplicative expressions surrounded by (), so
    minus/plus signs *inside* those factors are never regarded as Plus heads.
    """
    text = raw.strip()
    if not text:
        raise ValueError("Empty LLSS source")
    starts = []
    sq = curly = paren = 0
    start = 0
    for i, c in enumerate(text):
        if c == "[": sq += 1
        elif c == "]": sq -= 1
        elif c == "{": curly += 1
        elif c == "}": curly -= 1
        elif c == "(": paren += 1
        elif c == ")": paren -= 1
        elif _is_matchete_sum_separator(c, (sq, curly, paren)) and i != start:
            # '+'/'-' may also occur in a Mathematica scientific exponent.
            if i and text[i-1] in "eE" and i + 1 < len(text) and text[i+1].isdigit():
                continue
            starts.append(text[start:i].strip())
            start = i
        if min(sq, curly, paren) < 0:
            raise ValueError("Unbalanced LLSS expression")
    if (sq, curly, paren) != (0, 0, 0):
        raise ValueError("Unbalanced LLSS expression")
    starts.append(text[start:].strip())
    result = []
    for term in starts:
        if term.startswith("+"):
            term = term[1:].strip()
        if not term or term in ("-", "+") or "GammaCC" not in term:
            raise ValueError("Non-LLSS/empty term in matched source")
        result.append(term)
    return result


def _complex_matrix(obj: Any, shape: tuple[int, int]) -> np.ndarray:
    if not isinstance(obj, Mapping) or set(obj) != {"real", "imag"}:
        raise ValueError("Invalid complex matrix encoding")
    array = np.asarray(obj["real"], dtype=float) + 1j * np.asarray(obj["imag"], dtype=float)
    if array.shape != shape or not np.isfinite(array).all():
        raise ValueError("Invalid finite complex matrix shape")
    return array


def _decode_complex(z: Any) -> complex:
    if isinstance(z, Mapping) and set(z) == {"real", "imag"}:
        out = complex(float(z["real"]), float(z["imag"]))
    else:
        raise ValueError("Malformed complex value")
    if not np.isfinite(out.real) or not np.isfinite(out.imag):
        raise ValueError("Nonfinite complex value")
    return out


def _spurion_map(basis: Mapping[str, Any]) -> dict[str, complex]:
    mat = _complex_matrix(basis["basis"]["unitary_flavour_rotation_U"], (2, 2))
    if not np.allclose(mat.conj().T @ mat, np.eye(2), atol=1e-10, rtol=1e-10):
        raise ValueError("Uncertified scalar rotation U")
    return {f"T3MEU{i+1}{j+1}": complex(mat[i,j]) for i in range(2) for j in range(2)}


def factorise_matched_term(raw: str, index: int, spurions: Mapping[str, complex]) -> dict:
    """Remove exactly 2 ordered y and 1 inverse MF; retain ALL other terms.

    Spurion couplings are numerically substituted only after their four values
    have been validated against the synchronized scalar mass-basis U.
    """
    factors = _factors(raw)
    if len(factors) != 2 or len({f.heavy_index for f in factors}) != 1:
        raise ValueError(f"Term {index}: requires two Yukawas with common internal flavour")
    if len({f.external_index for f in factors}) != 2:
        raise ValueError(f"Term {index}: external flavour labels must be distinct")
    if any(f.name not in ("y1", "y2") for f in factors):
        raise ValueError(f"Term {index}: unsupported source Yukawa convention")
    mf = list(_MF.finditer(raw))
    if len(mf) != 1:
        raise ValueError(f"Term {index}: expected one heavy MF denominator")
    mass = mf[0]
    if mass.group("bar") or mass.group("close") or (
        mass.group("internal") not in (None, factors[0].heavy_index)
    ):
        raise ValueError(f"Term {index}: malformed MF factor")
    # The source has exactly one inverse mass. In all supplied Matchete terms
    # it occurs after '/', followed by a closing ')' possibly with products.
    before = raw[:mass.start()]
    if not (re.search(r"/\s*$", before) or
            re.search(r"/\s*\([^()]*\*\s*$", before)):
        raise ValueError(f"Term {index}: heavy mass must occur in denominator")
    original = raw
    if ("Proj[1]" in original) == ("Proj[-1]" in original):
        raise ValueError(f"Term {index}: missing or ambiguous chirality")
    # This MUST match RunMatching.wl's native Matchete convention:
    # Proj[-1] is the holomorphic L^T GammaCC PL L operator; Proj[1]
    # accompanies the Bar[L] ... Transp[Bar[L]] Hermitian conjugate.
    # The earlier mixed parser accidentally inverted the two labels.
    chirality = "PL" if "Proj[-1]" in original else "PR"
    if chirality == "PL":
        if "Transp[Field[l," not in original or "Transp[Bar[Field[l," in original:
            raise ValueError(f"Term {index}: PL projector does not have L^T GammaCC L")
    else:
        if "Transp[Bar[Field[l," not in original or "Transp[Field[l," in original:
            raise ValueError(f"Term {index}: PR projector does not have barred-lepton conjugate structure")
    matches = sorted([f.span for f in factors] + [mass.span()], reverse=True)
    for a, b in matches:
        raw = raw[:a] + "1" + raw[b:]
    if re.search(r"Coupling\[(?:MF|y1|y2)\b", raw):
        raise ValueError(f"Term {index}: residual source Yukawa/MF factor")
    seen = {m.group("name") for m in _SPURION.finditer(raw)}
    if seen - set(spurions):
        raise ValueError(f"Term {index}: undefined scalar U spurion")
    def number_text(z: complex) -> str:
        return f"({z.real:.17g}+({z.imag:.17g})*I)"
    def replace_spurion(match: re.Match[str]) -> str:
        barred = match.group("bar") is not None
        if barred != (match.group("close") is not None):
            raise ValueError(f"Term {index}: malformed barred mixing spurion")
        value = spurions[match.group("name")]
        return number_text(value.conjugate() if barred else value)
    template = _SPURION.sub(replace_spurion, raw)
    remaining = _ALL_COUPLING.findall(template)
    if remaining:
        raise ValueError(f"Term {index}: unsupported couplings remain: {remaining}")
    return {
        "source_index": index,
        "chirality": chirality,
        "operator_template_unit_flavor_inputform": template,
        "source_term_inputform": original,
        "factors": [{"name": f.name, "conjugated": f.conjugated,
                     "external_index": f.external_index, "heavy_index": f.heavy_index}
                    for f in factors],
        "flavour_matrix_axes_external_indices": [f.external_index for f in factors],
        "spurion_names": sorted(seen),
        "symmetric_flavour_combination_assumed": False,
    }


def build_mixed_source_terms(model: str, llss: str, basis: Mapping[str, Any],
                             masses: np.ndarray, yukawas: Mapping[str, np.ndarray],
                             *, expected_count: int | None = None) -> list[dict]:
    if model not in MIXED or basis.get("model") != model:
        raise ValueError("Model/mass-basis mismatch")
    terms = split_matchete_sum(llss)
    if expected_count is not None and len(terms) != expected_count:
        raise ValueError(f"{model}: {len(terms)} LLSS terms, expected {expected_count}")
    spurions = _spurion_map(basis)
    out = []
    for i, raw in enumerate(terms, 1):
        item = factorise_matched_term(raw, i, spurions)
        mat = flavour_matrix(item, masses, yukawas)
        item["ordered_flavour_matrix_gev_inv"] = {
            "real": mat.real.tolist(), "imag": mat.imag.tolist(),
            "frobenius_norm": float(np.linalg.norm(mat)),
        }
        out.append(item)
    return out


def _cg_used(terms: list[dict]) -> set[str]:
    from RGE.running.intermediate.ScalarOnlyTensorAdapters import _find_wilson_cg_calls
    used = set()
    for term in terms:
        for call in _find_wilson_cg_calls(term["operator_template_unit_flavor_inputform"]):
            used.add(call.name)
    return used


def project_mixed_terms(model: str, terms: list[dict], registry_payload: Mapping[str, Any]) -> list[dict]:
    if model not in MIXED or registry_payload.get("model") != model or (
        registry_payload.get("format") != "t3-mixed-matchete-native-cg-registry-v1"
    ):
        raise ValueError("CG registry provenance mismatch")
    registry = registry_payload.get("CGRegistry")
    if not isinstance(registry, list):
        raise ValueError("Matchete CG registry missing")
    names = [entry["Name"] for entry in registry]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate CG definitions")
    selected = [t for t in terms if t["chirality"] == "PL"]
    if len(selected) != len(terms)//2 or not selected:
        raise ValueError("No complete PL/PR source inventory")
    required = _cg_used(terms)
    if required != set(names):
        raise ValueError(f"Missing/unexpected CG entries: needed={required}, actual={set(names)}")
    return [project_ordered_term(t, registry, model) for t in selected]
