"""Extract representation-dependent C12 group factors from saved EFT1 tensors.

This does not derive a new RGE. It reads the actual Matchete Wilson CG seed
and the saved scalar-only master-tensor beta, isolates the mixed y1* y2*
channel, and verifies coefficient collinearity on the C12 seed support.

It deliberately reports no unique coefficient when proportionality fails.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

import sympy as sp


_SECTORS = (
    "yukawa_wavefunction", "scalar_pair", "mixed_yukawa_gauge",
    "crossed_yukawa", "conjugate_coefficient_yukawa",
    "scalar_anomalous_dimension", "fermion_anomalous_dimension", "total",
)


def _expression(value: object) -> sp.Expr:
    return sp.sympify(str(value), locals={"sqrt": sp.sqrt, "conjugate": sp.conjugate, "I": sp.I})


def _canonical(indices: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    i, j, a, b = indices
    return min(i, j), max(i, j), min(a, b), max(a, b)


def _mixed_coefficient(expression: sp.Expr) -> sp.Expr:
    """The independent mixed tree-level y1* y2* matching channel."""
    u = sp.conjugate(sp.Symbol("y1"))
    v = sp.conjugate(sp.Symbol("y2"))
    return sp.expand(expression).coeff(u).coeff(v)


def project_c12_group_factors(
    seed_entries: Mapping[tuple[int, int, int, int], sp.Expr],
    betas: Mapping[str, Mapping[str, object]],
) -> dict:
    """Return exact seed-support group factors or explicitly fail collinearity.

    `seed_entries` must use the matching adapter's separately symmetrised
    real-scalar tensor convention. `betas` uses independent component keys
    from the saved master RGE (e.g. `1,1,5,11`).
    """
    seed: dict[tuple[int, int, int, int], sp.Expr] = {}
    for indices, expression in seed_entries.items():
        key = _canonical(tuple(indices))
        value = sp.simplify(_mixed_coefficient(_expression(expression)))
        if value == 0:
            continue
        if key in seed and sp.simplify(seed[key] - value) != 0:
            raise ValueError(f"Inconsistent symmetric C12 seed component: {key}")
        seed[key] = value

    if not seed:
        raise ValueError("The matched Wilson seed has no mixed conjugate(y1)*conjugate(y2) C12 channel")

    saved: dict[tuple[int, int, int, int], Mapping[str, object]] = {}
    for raw_key, sectors in betas.items():
        key = _canonical(tuple(int(x.strip()) for x in raw_key.split(",")))
        if len(key) != 4:
            raise ValueError(f"Invalid saved Wilson component: {raw_key}")
        if key in saved:
            raise ValueError(f"Duplicate independent beta component: {key}")
        saved[key] = sectors

    results: dict[str, dict] = {}
    support = sorted(seed)
    for name in _SECTORS:
        ratios = []
        for component in support:
            raw = saved.get(component, {}).get(name, "0")
            numerator = _mixed_coefficient(_expression(raw))
            ratios.append(sp.simplify(sp.cancel(numerator / seed[component])))
        reference = ratios[0]
        failures = [
            {"component": list(key), "residual": sp.sstr(sp.factor(sp.simplify(value - reference)))}
            for key, value in zip(support[1:], ratios[1:])
            if sp.simplify(value - reference) != 0
        ]
        results[name] = {
            "collinear_on_c12_seed_support": not failures,
            "coefficient": sp.sstr(sp.factor(reference)) if not failures else None,
            "coefficient_latex": sp.latex(sp.factor(reference)) if not failures else None,
            "noncollinear_component_count": len(failures),
            "failures": failures[:12],
        }

    generated_mixed_components = []
    for key, sectors in sorted(saved.items()):
        if key in seed:
            continue
        if _mixed_coefficient(_expression(sectors.get("total", "0"))) != 0:
            generated_mixed_components.append(list(key))

    return {
        "projection_method": "exact mixed-yukawa coefficient; per-component seed-support collinearity",
        "seed_support_component_count": len(support),
        "reference_component": list(support[0]),
        "sectors": results,
        "c12_seed_support_is_invariant_under_projected_beta": results["total"]["collinear_on_c12_seed_support"],
        "generated_mixed_components_outside_c12_seed_support": generated_mixed_components,
        "full_operator_closure_verified": False,
    }


def derive_from_saved_data(data_dir: Path) -> dict:
    """Compute exact factors from existing saved output; never hard-code beta values."""
    data_dir = Path(data_dir)
    seed_file = data_dir / "eft1_after_F_wilson_seed.json"
    rgbeta_file = data_dir / "eft1_rgbeta_rge.json"
    beta_file = data_dir / "eft1_wilson_rge.json"
    for path in (seed_file, rgbeta_file, beta_file):
        if not path.is_file():
            raise FileNotFoundError(f"Required representation-specific RGE input missing: {path}")

    # Import only when actually deriving, so the projection and HTML tests can
    # run independently from the full Matchete/CG adapter.
    from RGE.running.intermediate.ScalarOnlyTensorAdapters import load_and_build_eft1_wilson_tensor

    tensor = load_and_build_eft1_wilson_tensor(
        seed_file, rgbeta_file, chirality="PL", project_operator_symmetry=True,
    )
    payload = json.loads(beta_file.read_text(encoding="utf-8"))
    rgbeta = json.loads(rgbeta_file.read_text(encoding="utf-8"))
    if payload.get("status") != "Success":
        raise ValueError("Saved EFT1 master Wilson beta does not have status Success")
    if payload.get("scan_mode") == "seed_support_only":
        raise ValueError("Complete independent-component beta scan is required for this audit")
    meta = rgbeta["metadata"]
    representation = {
        "dS1": int(meta["dS1"]), "dS2": int(meta["dS2"]),
        "dF": int(meta["dF"]), "alpha": int(meta["alpha"]),
        "shared_scalar": bool(meta.get("SharedScalar", False)),
    }
    if bool(meta.get("SharedScalar", False)):
        raise NotImplementedError("Shared-scalar C12 projector requires a separate channel definition")
    saved_meta = payload.get("metadata", {})
    for key in ("dS1", "dS2"):
        if int(saved_meta[key]) != representation[key]:
            raise ValueError(f"Matching and Wilson RGE disagree on {key}")

    result = project_c12_group_factors(tensor.entries, payload["betas"])
    result["representation"] = representation
    result["sources"] = {
        "matched_wilson_seed": str(seed_file),
        "rgbeta_representation_metadata": str(rgbeta_file),
        "master_wilson_beta": str(beta_file),
    }

    if representation == {"dS1": 2, "dS2": 2, "dF": 1, "alpha": -1, "shared_scalar": False}:
        total = result["sectors"]["total"]
        if total["coefficient"] is not None:
            raw = _expression(total["coefficient"])
            converted = sp.simplify(raw.subs({
                sp.Symbol("lambda12Inv2"): (
                    -sp.sqrt(3) * sp.Symbol("lambda12") / 2
                    + sp.sqrt(3) * sp.Symbol("lambda12Adj") / 8
                ),
                sp.Symbol("lambda12Inv1"): (
                    -sp.Symbol("lambda12") / 2
                    - 3 * sp.Symbol("lambda12Adj") / 8
                ),
            }))
            expected = (
                -3 * sp.Symbol("g2") ** 2 + 2 * sp.Symbol("lambda12")
                - sp.Symbol("lambda12Adj") / 2
                + sp.Symbol("ye") * sp.conjugate(sp.Symbol("ye"))
            )
            result["t3b_paper_basis"] = {
                "coefficient": sp.sstr(sp.factor(converted)),
                "coefficient_latex": sp.latex(sp.factor(converted)),
                "matches_previous_two_operator_diagonal_coefficient": (
                    sp.simplify(converted - expected) == 0
                ),
                "basis_mapping": {
                    "lambda12Inv1": "-lambda12/2 - 3*lambda12Adj/8",
                    "lambda12Inv2": "-sqrt(3)*lambda12/2 + sqrt(3)*lambda12Adj/8",
                },
            }
    return result
