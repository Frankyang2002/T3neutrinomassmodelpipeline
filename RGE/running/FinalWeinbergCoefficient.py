from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any


WEINBERG_KEYS = (
    "Weinberg",
    "weinberg",
    "weinberg_coefficient",
    "weinberg_running_coefficient",
    "direct_weinberg",
    "generated_weinberg",
)


# ---------------------------------------------------------------------------
# Pole/RGE consistency payload normalization
# Consolidated from the former standalone helper module.
# ---------------------------------------------------------------------------

def normalize_pole_rge_consistency(path: Path) -> dict[str, Any]:
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))

    if payload.get("Status") != "Success":
        raise ValueError("Pole/RGE diagnostic itself is not successful.")

    validated = bool(
        payload.get("DirectLogEqualsTwicePoleResidueOneGeneration", False)
    )

    payload["ConsistencyConvention"] = "strict_one_generation_direct_vs_hard"
    payload["ConsistencyValidated"] = validated
    payload["ConsistencyReason"] = (
        "Validated: direct running log coefficient = 2 x hard pole residue"
        if validated
        else
        "Failed strict regression: direct running log coefficient != "
        "2 x hard pole residue. No normalization fallback was applied."
    )

    # Remove fields left by the temporary orientation-resolved workaround so
    # old cached JSON cannot make a failed strict comparison look successful.
    for key in (
        "OriginalDirectLogEqualsTwicePoleResidueOneGeneration",
        "DirectFlavorOrientationCount",
        "HardPoleFlavorOrientationCount",
        "FlavorOrientationMultiplicity",
        "OrientationResolvedDirectLogToPoleRatio",
        "OrientationResolvedDirectLogEqualsTwicePoleResidue",
    ):
        payload.pop(key, None)

    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload

def _read_text(path: Path) -> str:
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError(f"Empty coefficient file: {path}")
    return text


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return payload


def _normalise_expression(value: Any) -> str | None:
    """Return a compact symbolic expression when a JSON value contains one."""
    if isinstance(value, str):
        value = value.strip()
        return value or None

    if isinstance(value, (int, float)):
        return str(value)

    if isinstance(value, dict):
        # Common expression-bearing field names used throughout the RGE files.
        for key in (
            "expression",
            "coefficient",
            "value",
            "running_correction",
            "delta",
            "beta",
            "transported",
        ):
            if key in value:
                result = _normalise_expression(value[key])
                if result is not None:
                    return result

    return None


def _find_direct_weinberg(payload: Any, path: tuple[str, ...] = ()) -> tuple[str, str] | None:
    """Find the direct running-generated Weinberg coefficient.

    Current EFT1WilsonFlavorRunning.py writes the result as

        running_corrections["Weinberg"]["running_tensor_text"]

    where ``running_tensor_text`` is already the complete full-flavor
    leading-log correction, e.g. a sum of flavor-blind and
    ``He.C + C.He^T`` tensor structures.  Prefer that exact schema.

    Older/development schemas are handled by the recursive fallback below.
    """
    if isinstance(payload, dict):
        corrections = payload.get("running_corrections")
        if isinstance(corrections, dict):
            weinberg = corrections.get("Weinberg")
            if isinstance(weinberg, dict):
                for key in (
                    "running_tensor_text",
                    "running_tensor_wolfram",
                    "expression",
                    "coefficient",
                ):
                    expr = _normalise_expression(weinberg.get(key))
                    if expr is not None:
                        return (
                            expr,
                            f"running_corrections.Weinberg.{key}",
                        )

        # Fallback for older/development layouts.  When a Weinberg-labelled
        # object is a dictionary, inspect expression-bearing child fields
        # rather than requiring the dictionary itself to normalise.
        for key, value in payload.items():
            key_text = str(key)
            lower = key_text.lower()
            if "weinberg" in lower:
                expr = _normalise_expression(value)
                if expr is not None:
                    return expr, ".".join((*path, key_text))
                if isinstance(value, dict):
                    for child_key in (
                        "running_tensor_text",
                        "running_tensor_wolfram",
                        "expression",
                        "coefficient",
                        "value",
                    ):
                        expr = _normalise_expression(value.get(child_key))
                        if expr is not None:
                            return (
                                expr,
                                ".".join((*path, key_text, child_key)),
                            )

        for key, value in payload.items():
            result = _find_direct_weinberg(value, (*path, str(key)))
            if result is not None:
                return result

    elif isinstance(payload, list):
        for index, value in enumerate(payload):
            result = _find_direct_weinberg(value, (*path, str(index)))
            if result is not None:
                return result

    return None


def _find_scalar_metadata(payload: Any, key_name: str) -> Any | None:
    if isinstance(payload, dict):
        if key_name in payload:
            return payload[key_name]
        for value in payload.values():
            found = _find_scalar_metadata(value, key_name)
            if found is not None:
                return found
    elif isinstance(payload, list):
        for value in payload:
            found = _find_scalar_metadata(value, key_name)
            if found is not None:
                return found
    return None


def _strip_outer_quotes(text: str) -> str:
    if len(text) >= 2 and text[0] == text[-1] and text[0] in {"'", '"'}:
        return text[1:-1].strip()
    return text


def _looks_flavor_indexed(expr: str) -> bool:
    """Conservative indicator that an expression explicitly retains p,q."""
    patterns = (
        r"\bp\b",
        r"\bq\b",
        r"\[[^\]]*p[^\]]*,[^\]]*q[^\]]*\]",
        r"\([^\)]*p[^\)]*,[^\)]*q[^\)]*\)",
    )
    return any(re.search(pattern, expr) for pattern in patterns)


def _contains_uv_pole(expr: str) -> bool:
    """Detect common textual forms of an explicit dimensional-regulator pole.

    Matchete/Wolfram InputForm commonly serialises epsilon as ``\\[Epsilon]``.
    The raw expression is always preserved.  A pole is removed from the
    separate renormalized-hard object only after the independent pole/RGE
    consistency diagnostic has validated the subtraction.
    """
    compact = expr.replace(" ", "")
    lower = compact.lower()

    # Wolfram escaped symbol form, e.g. 1/(\\[Epsilon]*MF).
    if "\\[epsilon]" in lower:
        return True

    patterns = (
        r"1/[εϵ]",
        r"1/epsilon\b",
        r"[εϵ]\*\*-1",
        r"epsilon\*\*-1",
        r"Power\[[εϵ],-1\]",
        r"Power\[epsilon,-1\]",
    )
    return any(
        re.search(pattern, compact, flags=re.IGNORECASE)
        for pattern in patterns
    )



def _load_pole_rge_consistency(
    threshold_c5_path: Path,
    explicit_path: Path | None = None,
) -> tuple[dict[str, Any] | None, Path | None]:
    """Load the independent hard-pole/RGE consistency diagnostic.

    RunThresholdStage.wl writes ``c5_pole_rge_consistency.json`` next to the
    authoritative threshold C5.  Keep this optional for backward
    compatibility, but never call the pole subtraction validated unless this
    file explicitly says that the direct running log coefficient equals twice
    the hard 1/epsilon residue.
    """
    candidate = (
        Path(explicit_path)
        if explicit_path is not None
        else threshold_c5_path.parent / "c5_pole_rge_consistency.json"
    )

    if not candidate.exists():
        return None, candidate

    payload = _read_json(candidate)
    return payload, candidate


def _pole_rge_validation_status(
    payload: dict[str, Any] | None,
) -> tuple[bool, str]:
    if payload is None:
        return False, "Pole/RGE consistency diagnostic is missing"

    status = str(payload.get("Status", payload.get("status", ""))).lower()
    if status and status != "success":
        return False, "Pole/RGE consistency diagnostic did not succeed"

    relation = payload.get(
        "DirectLogEqualsTwicePoleResidueOneGeneration",
        payload.get("direct_log_equals_twice_pole_residue_one_generation"),
    )
    ratio = payload.get(
        "DirectLogToPoleRatioOneGenerationInputForm",
        payload.get("direct_log_to_pole_ratio_one_generation"),
    )

    if relation is not True:
        return (
            False,
            "Direct-running log coefficient was not validated as twice the "
            "hard UV-pole residue",
        )

    if ratio is not None and str(ratio).strip() not in {"2", "2.0"}:
        return (
            False,
            "Pole/RGE diagnostic reported a non-2 direct-log/pole ratio",
        )

    return True, "Validated: direct running log coefficient = 2 x hard pole residue"


def _drop_wolfram_epsilon_terms(expr: str) -> tuple[str, int]:
    """Drop additive terms proportional to the explicit Wolfram UV pole.

    This is intentionally narrow: it recognizes the Matchete InputForm shape

        +(numerator)/(\\[Epsilon]*denominator)

    seen in the authoritative threshold coefficient.  The raw Matchete
    expression is always retained separately.  Call this only after the
    pole/RGE relation has been independently validated.
    """
    pattern = re.compile(
        r"\s*\+\s*"
        r"\((?P<num>[^()]*)\)"
        r"/\("
        r"\\\[Epsilon\]\*"
        r"(?P<den>[^()]*)"
        r"\)"
    )
    renormalized, count = pattern.subn("", expr)
    return renormalized, count


def _set_matchete_matching_scale_to_mf(expr: str) -> tuple[str, int]:
    """Set the explicit Matchete MS-bar scale to the first threshold MF.

    At mu_bar^2 = MF^2 the hard matching logarithm vanishes.  This turns the
    finite hard threshold into the boundary coefficient at the matching scale;
    subsequent MF -> MS evolution is supplied by the separate EFT1 running
    contribution.
    """
    patterns = (
        re.compile(
            r"Log\["
            r"\\\[Mu\]bar2/"
            r"Coupling\[MF,\s*\{\},\s*0\]\^2"
            r"\]"
        ),
        re.compile(r"Log\[\\\[Mu\]bar2/MF\^2\]"),
    )

    out = expr
    replacements = 0
    for pattern in patterns:
        out, count = pattern.subn("0", out)
        replacements += count
    return out, replacements


def _build_msbar_hard_at_mf(
    threshold_expr: str,
    *,
    validation_ok: bool,
) -> dict[str, Any]:
    """Construct the finite hard threshold coefficient at mu = MF.

    The subtraction is performed only after the independent pole/RGE
    consistency check succeeds.  Otherwise the object is returned as blocked
    and the raw threshold expression remains the only authoritative hard
    expression.
    """
    if not validation_ok:
        return {
            "status": "Blocked",
            "scheme": "MSbar",
            "matching_scale": "MF",
            "expression": None,
            "pole_terms_removed": 0,
            "matching_logs_set_to_zero": 0,
            "reason": (
                "UV-pole subtraction requires a validated pole/RGE "
                "consistency relation"
            ),
        }

    without_pole, pole_count = _drop_wolfram_epsilon_terms(threshold_expr)
    at_mf, log_count = _set_matchete_matching_scale_to_mf(without_pole)
    still_has_pole = _contains_uv_pole(at_mf)

    status = "Success" if pole_count > 0 and not still_has_pole else "Failed"
    reason = None
    if pole_count == 0:
        reason = (
            "Validated subtraction requested, but no recognized Matchete "
            "1/epsilon additive term was removed"
        )
    elif still_has_pole:
        reason = "Renormalized hard expression still contains a UV pole"

    return {
        "status": status,
        "scheme": "MSbar",
        "matching_scale": "MF",
        "expression": at_mf if status == "Success" else None,
        "pole_terms_removed": pole_count,
        "matching_logs_set_to_zero": log_count,
        "contains_uv_pole_after_subtraction": still_has_pole,
        "reason": reason,
    }



def _swap_matchete_external_flavor_indices(
    expr: str,
    p_index: str,
    q_index: str,
) -> str:
    """Swap only the two external Matchete Flavor indices.

    The same printed dummy name may also occur in NFlavor, so replacements
    must be restricted to ``Index[..., Flavor]`` and must not touch the
    internal fermion-generation index.
    """
    p_token = f"Index[{p_index}, Flavor]"
    q_token = f"Index[{q_index}, Flavor]"
    marker = "__T3_EXTERNAL_FLAVOR_SWAP__"

    if p_token not in expr or q_token not in expr:
        return expr

    return (
        expr.replace(p_token, marker)
        .replace(q_token, p_token)
        .replace(marker, q_token)
    )


def _build_physical_majorana_hard(
    renormalized_hard_expr: str | None,
    threshold_info: dict[str, Any],
) -> dict[str, Any]:
    """Construct the physical symmetric Majorana coefficient.

    Matchete returns an ordered p,q representative.  The Weinberg coefficient
    is symmetric in its two lepton-flavor indices, so construct

        C_phys[p,q] = C_ordered[p,q] + C_ordered[q,p].

    In one generation the physical coefficient is twice the ordered Matchete coefficient.
    """
    if not isinstance(renormalized_hard_expr, str):
        return {
            "status": "Blocked",
            "expression": None,
            "reason": "Renormalized hard threshold is unavailable",
        }

    external = list(threshold_info.get("external_flavor_indices", []))
    if len(external) < 2:
        return {
            "status": "Blocked",
            "expression": None,
            "reason": "Two explicit external Flavor indices were not found",
        }

    p_index, q_index = external[:2]
    swapped = _swap_matchete_external_flavor_indices(
        renormalized_hard_expr,
        p_index,
        q_index,
    )

    if swapped == renormalized_hard_expr:
        return {
            "status": "Blocked",
            "expression": None,
            "reason": "Could not construct a distinct p<->q swapped hard term",
        }

    symmetric = f"(({renormalized_hard_expr}) + ({swapped}))"

    return {
        "status": "Success",
        "symmetrization": "C[p,q] + C[q,p]",
        "expression": symmetric,
        "ordered_expression": renormalized_hard_expr,
        "swapped_expression": swapped,
        "external_flavor_indices": [p_index, q_index],
        "one_generation_ordered_to_physical_factor": 2,
        "reason": None,
    }


def _swap_generic_pq(expr: str) -> str:
    """Swap generic p,q labels without touching substrings in other names."""
    marker = "__T3_GENERIC_PQ_SWAP__"
    out = expr
    replacements = (
        ("[p,", f"[{marker},"),
        ("[q,", "[p,"),
        (f"[{marker},", "[q,"),
        (",p]", f",{marker}]"),
        (",q]", ",p]"),
        (f",{marker}]", ",q]"),
    )
    for old, new in replacements:
        out = out.replace(old, new)
    return out


def _direct_running_is_manifestly_symmetric(expr: str) -> bool:
    """Conservative text-level p<->q symmetry check for expanded C12 kernels."""
    compact = re.sub(r"\s+", "", expr)
    swapped = re.sub(r"\s+", "", _swap_generic_pq(expr))

    if compact == swapped:
        return True

    # The current C12 expansion is a symmetric sum of the two ordered kernels.
    required = (
        "conjugate(y1[p,r])conjugate(y2[q,r])",
        "conjugate(y2[p,r])conjugate(y1[q,r])",
    )
    return all(piece in compact for piece in required)


def _contains_explicit_hbar(expr: str) -> bool:
    compact = expr.replace(" ", "")
    return bool(
        re.search(r"\bhbar\b", compact, flags=re.IGNORECASE)
        or "1/(16*pi**2)" in compact.lower()
        or "1/(16*pi^2)" in compact.lower()
    )


def _contains_unexpanded_c12(expr: str) -> bool:
    """True when the direct running still refers to the symbolic C12 tensor."""
    return bool(
        re.search(r"\bC12\s*[\[(]", expr)
        or re.search(r"\bEFT1C12\s*[\[(]", expr)
    )


def _expand_boundary_tensors(
    expr: str,
    payload: dict[str, Any],
) -> tuple[str, dict[str, str]]:
    """Expand symbolic Cij[p,q] tensors using transport boundary kernels."""
    boundaries = payload.get("boundary_tensors")
    if not isinstance(boundaries, dict):
        return expr, {}

    expanded = expr
    replacements: dict[str, str] = {}

    for name in ("C11", "C12", "C22"):
        info = boundaries.get(name)
        if not isinstance(info, dict):
            continue

        definition = _normalise_expression(info.get("definition_text"))
        if definition is None:
            continue

        replacements[name] = definition
        for pattern in (
            rf"\b{re.escape(name)}\s*\[\s*p\s*,\s*q\s*\]",
            rf"\b{re.escape(name)}\s*\(\s*p\s*,\s*q\s*\)",
        ):
            expanded = re.sub(
                pattern,
                lambda _m, d=definition: f"({d})",
                expanded,
            )

    return expanded, replacements


def _matchete_threshold_flavor_info(expr: str) -> dict[str, Any]:
    """Inspect a Matchete C5 coefficient for explicit flavor structure."""
    flavor_tokens = re.findall(r"Index\[([^,\]]+),\s*Flavor\]", expr)
    nflavor_tokens = re.findall(r"Index\[([^,\]]+),\s*NFlavor\]", expr)

    unique_flavor = list(dict.fromkeys(flavor_tokens))
    unique_nflavor = list(dict.fromkeys(nflavor_tokens))

    if len(unique_flavor) < 2:
        return {
            "explicit_full_flavor": False,
            "external_flavor_indices": unique_flavor,
            "internal_nflavor_indices": unique_nflavor,
            "normalised_expression": None,
        }

    normalised = expr
    for symbol, replacement in zip(unique_flavor[:2], ("p", "q")):
        normalised = re.sub(
            rf"Index\[{re.escape(symbol)},\s*Flavor\]",
            replacement,
            normalised,
        )

    internal_labels = ("r", "s", "t", "u")
    for symbol, replacement in zip(unique_nflavor, internal_labels):
        normalised = re.sub(
            rf"Index\[{re.escape(symbol)},\s*NFlavor\]",
            replacement,
            normalised,
        )

    normalised = re.sub(
        r"Bar\[Coupling\[([A-Za-z0-9_]+),\s*\{([^}]*)\},\s*0\]\]",
        lambda m: f"conjugate({m.group(1)}[{m.group(2).replace(' ', '')}])",
        normalised,
    )
    normalised = re.sub(
        r"Coupling\[([A-Za-z0-9_]+),\s*\{\},\s*0\]",
        r"\1",
        normalised,
    )
    normalised = normalised.replace("Sqrt[3]", "sqrt(3)")
    normalised = re.sub(r"Log\[([^\]]+)\]", r"log(\1)", normalised)

    return {
        "explicit_full_flavor": True,
        "external_flavor_indices": unique_flavor[:2],
        "internal_nflavor_indices": unique_nflavor,
        "normalised_expression": normalised,
    }


def build_final_weinberg_coefficient(
    *,
    threshold_c5_path: Path,
    flavor_transport_path: Path,
    output_path: Path | None = None,
    pole_rge_consistency_path: Path | None = None,
) -> dict[str, Any]:
    """Build the authoritative final-C5 bookkeeping object.

    Renormalized hierarchical fixed-order structure:
        C5_final^{pq}(MS)
          = C5_hard,MSbar^{pq}(MF)
          + hbar * Delta C5_run,direct^(1),pq(MF -> MS).

    The raw Matchete threshold expression is preserved unchanged for
    diagnostics.  A separate finite hard threshold is constructed only when
    ``c5_pole_rge_consistency.json`` validates that the direct-running
    logarithmic coefficient equals twice the hard 1/epsilon residue.

    Here the EFT1 flavor-running JSON stores the order-one one-loop object
    Delta C5_run,direct^(1), not the physical hbar-weighted correction.

    The current Matchete `c5_coefficient.txt` is a stripped coefficient: its
    operator flavor indices are no longer explicit.  We therefore must NOT
    silently pretend that scalar text is already a full-flavor matrix.

    Until the threshold contribution itself has a full-flavor lift, the
    combined full-flavor result is represented formally as

        C5Threshold[p,q] + DeltaC5Run[p,q]

    while preserving the exact scalar-stripped Matchete coefficient alongside
    it.  This prevents loss of the p,q dependence generated by the EFT1 RGE.
    """
    threshold_c5_path = Path(threshold_c5_path)
    flavor_transport_path = Path(flavor_transport_path)

    threshold_expr = _strip_outer_quotes(_read_text(threshold_c5_path))
    transport = _read_json(flavor_transport_path)

    pole_rge_payload, resolved_pole_rge_path = _load_pole_rge_consistency(
        threshold_c5_path,
        pole_rge_consistency_path,
    )
    pole_rge_validated, pole_rge_reason = _pole_rge_validation_status(
        pole_rge_payload
    )

    found = _find_direct_weinberg(transport)
    if found is None:
        raise ValueError(
            "Could not find the direct Weinberg running correction in "
            f"{flavor_transport_path}. Expected "
            "running_corrections['Weinberg']['running_tensor_text']."
        )

    direct_expr_raw, direct_source = found
    direct_expr, boundary_replacements = _expand_boundary_tensors(
        direct_expr_raw,
        transport,
    )

    log_ratio = _find_scalar_metadata(transport, "log_ratio")
    equal_scale = _find_scalar_metadata(
        transport, "equal_scale_running_vanishes"
    )
    one_gen_matches = _find_scalar_metadata(
        transport, "one_generation_reduction_matches"
    )

    direct_is_flavor_indexed = _looks_flavor_indexed(direct_expr)
    threshold_info = _matchete_threshold_flavor_info(threshold_expr)

    threshold_has_uv_pole = _contains_uv_pole(threshold_expr)
    threshold_has_explicit_hbar = _contains_explicit_hbar(threshold_expr)
    direct_has_explicit_hbar = _contains_explicit_hbar(direct_expr)
    direct_requires_hbar = not direct_has_explicit_hbar
    c12_boundary_expanded = not _contains_unexpanded_c12(direct_expr)

    renormalized_hard = _build_msbar_hard_at_mf(
        threshold_expr,
        validation_ok=pole_rge_validated,
    )
    renormalized_hard_expr = renormalized_hard.get("expression")
    renormalized_hard_ok = (
        renormalized_hard.get("status") == "Success"
        and isinstance(renormalized_hard_expr, str)
        and not _contains_uv_pole(renormalized_hard_expr)
    )

    physical_hard = _build_physical_majorana_hard(
        renormalized_hard_expr,
        threshold_info,
    )
    physical_hard_expr = physical_hard.get("expression")
    physical_hard_ok = (
        physical_hard.get("status") == "Success"
        and isinstance(physical_hard_expr, str)
    )
    direct_running_symmetric = _direct_running_is_manifestly_symmetric(
        direct_expr
    )
    # The direct-running object was built in the previous averaged convention.
    # Convert it only at the final physical-C5 boundary.
    physical_direct_expr = (
        f"2*({direct_expr})" if direct_running_symmetric else None
    )

    threshold_flavor_status = (
        "explicit_full_flavor_matchete"
        if threshold_info["explicit_full_flavor"]
        else "stripped_scalar_coefficient"
    )
    direct_flavor_status = (
        "explicit_full_flavor"
        if direct_is_flavor_indexed
        else "flavor_kernel_or_symbolic_expression"
    )

    formal_threshold = "C5HardMSbarAtMF[p,q]"
    formal_direct = "DeltaC5RunDirect1[p,q]"
    formal_combined = f"{formal_threshold} + hbar*{formal_direct}"

    actual_combined_ordered = (
        f"({renormalized_hard_expr}) + hbar*({direct_expr})"
        if renormalized_hard_ok
        else None
    )
    actual_combined_physical = (
        f"({physical_hard_expr}) + hbar*({physical_direct_expr})"
        if physical_hard_ok and physical_direct_expr is not None
        else None
    )

    structurally_full_flavor = bool(
        threshold_info["explicit_full_flavor"]
        and direct_is_flavor_indexed
    )
    ready_for_ordered_full_flavor_numerics = bool(
        structurally_full_flavor
        and c12_boundary_expanded
        and pole_rge_validated
        and renormalized_hard_ok
    )
    ready_for_physical_majorana_numerics = bool(
        ready_for_ordered_full_flavor_numerics
        and physical_hard_ok
        and direct_running_symmetric
    )

    result: dict[str, Any] = {
        "status": "Success",
        "scheme": "fixed_order_one_loop_MSbar_hierarchical",
        "formula": (
            "C5_final[p,q](MS) = C5_hard_MSbar[p,q](MF) "
            "+ hbar*DeltaC5_run_direct^(1)[p,q](MF->MS)"
        ),
        "bookkeeping": {
            "threshold_expression_already_contains_hbar": (
                threshold_has_explicit_hbar
            ),
            "direct_running_expression_already_contains_hbar": (
                direct_has_explicit_hbar
            ),
            "direct_running_requires_hbar": direct_requires_hbar,
            "direct_running_object_order": "one_loop_order_one_coefficient",
        },
        "threshold_contribution_raw": {
            "source": str(threshold_c5_path),
            "expression": threshold_expr,
            "contains_uv_pole": threshold_has_uv_pole,
            "role": "unrenormalized Matchete hard-region diagnostic",
            "flavor_status": threshold_flavor_status,
            "external_flavor_indices": threshold_info["external_flavor_indices"],
            "internal_nflavor_indices": threshold_info["internal_nflavor_indices"],
            "normalised_full_flavor_expression": threshold_info["normalised_expression"],
        },
        "renormalized_hard_threshold": {
            **renormalized_hard,
            "source": str(threshold_c5_path),
            "full_flavor_symbol": formal_threshold,
            "flavor_status": threshold_flavor_status,
        },
        "physical_majorana_hard_threshold": {
            **physical_hard,
            "source": str(threshold_c5_path),
            "flavor_status": (
                "explicit_full_flavor_symmetric"
                if physical_hard_ok
                else "blocked"
            ),
        },
        "direct_running_contribution": {
            "source": str(flavor_transport_path),
            "source_location": direct_source,
            "raw_expression": direct_expr_raw,
            "expression": direct_expr,
            "physical_majorana_expression": physical_direct_expr,
            "physical_majorana_conversion_factor": 2,
            "boundary_tensor_replacements": boundary_replacements,
            "flavor_status": direct_flavor_status,
            "full_flavor_symbol": formal_direct,
            "manifestly_symmetric_under_pq": direct_running_symmetric,
        },
        "combined": {
            "full_flavor_expression": formal_combined,
            "ordered_expression": actual_combined_ordered,
            "physical_majorana_expression": actual_combined_physical,
            "structurally_full_flavor": structurally_full_flavor,
            "ready_for_ordered_full_flavor_numerics": (
                ready_for_ordered_full_flavor_numerics
            ),
            "ready_for_full_flavor_numerics": (
                ready_for_physical_majorana_numerics
            ),
            "ready_for_physical_majorana_numerics": (
                ready_for_physical_majorana_numerics
            ),
            "physical_majorana_reason": (
                None
                if ready_for_physical_majorana_numerics
                else "; ".join(
                    reason
                    for reason in (
                        (
                            physical_hard.get("reason") or
                            "Physical Majorana hard coefficient unavailable"
                            if not physical_hard_ok
                            else ""
                        ),
                        (
                            "Direct running contribution is not manifestly "
                            "symmetric under p<->q"
                            if not direct_running_symmetric
                            else ""
                        ),
                    )
                    if reason
                )
            ),
            "reason": (
                None
                if ready_for_ordered_full_flavor_numerics
                else "; ".join(
                    reason
                    for reason in (
                        (
                            "At least one contribution does not expose explicit "
                            "p,q flavor structure"
                            if not structurally_full_flavor
                            else ""
                        ),
                        (
                            pole_rge_reason
                            if not pole_rge_validated
                            else ""
                        ),
                        (
                            renormalized_hard.get("reason") or
                            "Renormalized hard threshold was not constructed"
                            if not renormalized_hard_ok
                            else ""
                        ),
                        (
                            "Direct running still contains symbolic C12[p,q] "
                            "instead of the expanded y1/y2 flavor kernel"
                            if not c12_boundary_expanded
                            else ""
                        ),
                    )
                    if reason
                )
            ),
        },
        "diagnostics": {
            "threshold_has_uv_pole": threshold_has_uv_pole,
            "raw_threshold_preserved": True,
            "pole_rge_consistency_source": (
                str(resolved_pole_rge_path)
                if resolved_pole_rge_path is not None
                else None
            ),
            "pole_rge_consistency_validated": pole_rge_validated,
            "pole_rge_consistency_reason": pole_rge_reason,
            "pole_rge_consistency_payload": pole_rge_payload,
            "renormalized_hard_threshold_constructed": renormalized_hard_ok,
            "physical_majorana_hard_constructed": physical_hard_ok,
            "direct_running_manifestly_symmetric": direct_running_symmetric,
            "direct_running_requires_hbar": direct_requires_hbar,
            "c12_boundary_expanded": c12_boundary_expanded,
            "c12_symmetry_convention_resolved": True,
        },
        "checks": {
            "equal_scale_running_vanishes": equal_scale,
            "one_generation_reduction_matches": one_gen_matches,
            "log_ratio": log_ratio,
        },
    }

    # A one-generation expression is still useful and is unambiguous: in one
    # generation the p,q structure collapses to the scalar coefficient.
    result["one_generation_combined_expression"] = (
        f"({physical_hard_expr}) + hbar*({physical_direct_expr})"
        if physical_hard_ok and physical_direct_expr is not None
        else None
    )

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(result, indent=2),
            encoding="utf-8",
        )

    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Combine the threshold-generated Weinberg coefficient with the "
            "direct EFT1 running-generated Weinberg contribution."
        )
    )
    parser.add_argument(
        "threshold_c5",
        type=Path,
        help="Final threshold Matchete c5_coefficient.txt",
    )
    parser.add_argument(
        "flavor_transport",
        type=Path,
        help="eft1_wilson_flavor_at_S_threshold.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("final_weinberg_coefficient.json"),
        help="Output JSON path",
    )
    parser.add_argument(
        "--pole-rge-consistency",
        type=Path,
        default=None,
        help=(
            "Optional c5_pole_rge_consistency.json. If omitted, the file is "
            "auto-discovered next to the threshold C5."
        ),
    )
    return parser


def main() -> int:
    args = _parser().parse_args()

    try:
        result = build_final_weinberg_coefficient(
            threshold_c5_path=args.threshold_c5,
            flavor_transport_path=args.flavor_transport,
            output_path=args.output,
            pole_rge_consistency_path=args.pole_rge_consistency,
        )
    except Exception as exc:
        print(
            json.dumps(
                {
                    "status": "Failed",
                    "error": str(exc),
                },
                indent=2,
            )
        )
        return 1

    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
