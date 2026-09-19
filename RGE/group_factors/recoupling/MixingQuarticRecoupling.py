from __future__ import annotations

"""Exact scalar-quartic recoupling factors in beta_lambdaT3.

This module uses the existing Matchete -> real-scalar quartic tensor adapter.
No A--E coefficient table is encoded.

Convention
----------
The adapter constructs the fully symmetric real tensor lambda_abcd from

    V4 = (1/4!) lambda_abcd phi_a phi_b phi_c phi_d.

For this convention the pure-scalar one-loop contribution is

    beta_abcd|_{lambda^2}
      = sum_{e,f} [
          lambda_abef lambda_efcd
        + lambda_acef lambda_efbd
        + lambda_adef lambda_efbc
        ].

This is equivalent to (1/8) times the full 24-permutation sum.

For a target real quartic X, write

    lambda_abcd = lambdaT3 * T_abcd + X * X_abcd + ...

The coefficient linear in lambdaT3 * X is then built exactly from the
cross-contractions T*X + X*T in the three channels and projected back onto T.

The output reports
- the exact projection coefficient R_X,
- the projection residual,
- whether the full generated tensor is proportional to the lambdaT3 tensor.

The lambdaT3 coupling is kept holomorphic: conjugate(lambdaT3) is treated as
an independent direction.  The adjoint/cross quartics tested here are real, so
their conjugates are identified with themselves before extracting their basis
tensors.
"""

import argparse
from dataclasses import dataclass, asdict
import json
from pathlib import Path
import sys
from typing import Iterable

import sympy as sp


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.group_factors.core.MixingQuarticTensorAlgebra import (
    basis_tensor,
    build_cross_tensor,
    project_onto_direction,
    scalar_dimension,
    tensor_inner_product,
)
from RGE.running.EFT1TensorAdapters import (
    SparseQuarticTensor,
    load_and_build_eft1_quartic_tensor,
)


REAL_TARGETS = (
    "lambdaH1Adj",
    "lambdaH2Adj",
    "lambda12Adj",
    "lambda12Cross",
)


@dataclass(frozen=True)
class ProjectionResult:
    coupling: str
    present_in_seed: bool
    t3_nonzero_components: int
    generated_nonzero_components: int
    projection_coefficient: str | None
    norm_t3: str
    max_residual: str
    residual_nonzero_components: int
    proportional_to_t3: bool


def _text(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.simplify(expr)))


def _max_residual_text(residual: SparseQuarticTensor) -> str:
    if not residual.entries:
        return "0"

    # Exact expressions do not have a useful universal ordering.  Return a
    # representative simplified nonzero residual rather than a floating norm.
    key = sorted(residual.entries)[0]
    return _text(residual.entries[key])


def recoupling_for_seed(
    quartic_seed_path: Path,
    rgbeta_path: Path,
    *,
    targets: Iterable[str] = REAL_TARGETS,
) -> dict:
    full = load_and_build_eft1_quartic_tensor(
        quartic_seed_path,
        rgbeta_path,
    )
    scalar_dimension_value = scalar_dimension(full)

    # lambdaT3 is complex.  Keep only its holomorphic tensor direction.
    t3 = basis_tensor(
        full,
        "lambdaT3",
        identify_conjugate=False,
    )
    if not t3.entries:
        raise ValueError(
            f"No holomorphic lambdaT3 tensor found in {quartic_seed_path}."
        )

    results: list[ProjectionResult] = []

    for coupling_name in targets:
        other = basis_tensor(
            full,
            coupling_name,
            identify_conjugate=True,
        )

        if not other.entries:
            results.append(
                ProjectionResult(
                    coupling=coupling_name,
                    present_in_seed=False,
                    t3_nonzero_components=len(t3.entries),
                    generated_nonzero_components=0,
                    projection_coefficient=None,
                    norm_t3=_text(tensor_inner_product(t3, t3)),
                    max_residual="0",
                    residual_nonzero_components=0,
                    proportional_to_t3=True,
                )
            )
            continue

        generated = build_cross_tensor(
            t3,
            other,
            scalar_dimension_value,
        )
        coefficient, residual, norm = project_onto_direction(
            generated,
            t3,
        )

        results.append(
            ProjectionResult(
                coupling=coupling_name,
                present_in_seed=True,
                t3_nonzero_components=len(t3.entries),
                generated_nonzero_components=len(generated.entries),
                projection_coefficient=_text(coefficient),
                norm_t3=_text(norm),
                max_residual=_max_residual_text(residual),
                residual_nonzero_components=len(residual.entries),
                proportional_to_t3=not residual.entries,
            )
        )

    return {
        "status": (
            "Success"
            if all(item.proportional_to_t3 for item in results)
            else "Failed"
        ),
        "quartic_seed": str(quartic_seed_path),
        "rgbeta_file": str(rgbeta_path),
        "scalar_dimension": scalar_dimension_value,
        "convention": (
            "V4=lambda_abcd phi_a phi_b phi_c phi_d/4!; "
            "beta_lambda|scalar^2 has the three pair-contraction channels"
        ),
        "lambdaT3_holomorphic_nonzero_components": len(t3.entries),
        "results": [asdict(item) for item in results],
    }


def _default_rgbeta(seed: Path) -> Path:
    candidate = seed.with_name("eft1_rgbeta_rge.json")
    if not candidate.is_file():
        raise FileNotFoundError(
            f"Could not infer EFT1 RGBeta metadata next to seed: {candidate}"
        )
    return candidate


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Derive lambdaT3 non-singlet scalar-quartic recoupling factors "
            "from exact Matchete real-scalar tensors."
        )
    )
    parser.add_argument(
        "quartic_seed_json",
        type=Path,
        help="eft1_after_F_scalar_quartic_seed.json",
    )
    parser.add_argument(
        "--rgbeta",
        type=Path,
        default=None,
        help="eft1_rgbeta_rge.json; inferred from the seed directory by default.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional JSON output path.",
    )
    args = parser.parse_args()

    rgbeta = args.rgbeta or _default_rgbeta(args.quartic_seed_json)
    payload = recoupling_for_seed(
        args.quartic_seed_json,
        rgbeta,
    )

    print(f"seed: {args.quartic_seed_json}")
    print(f"real scalar dimension: {payload['scalar_dimension']}")
    print(
        "lambdaT3 holomorphic tensor components: "
        f"{payload['lambdaT3_holomorphic_nonzero_components']}"
    )
    print()
    print(
        f"{'coupling':<20} {'present':<8} {'R_X':>12} "
        f"{'residual nz':>12} {'closure':>9}"
    )
    print("-" * 68)

    for item in payload["results"]:
        coefficient = item["projection_coefficient"]
        print(
            f"{item['coupling']:<20} "
            f"{str(item['present_in_seed']):<8} "
            f"{str(coefficient):>12} "
            f"{item['residual_nonzero_components']:>12} "
            f"{('PASS' if item['proportional_to_t3'] else 'FAIL'):>9}"
        )

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(payload, indent=2),
            encoding="utf-8",
        )
        print()
        print(f"JSON summary: {args.output}")

    print(f"OVERALL: {payload['status']}")
    return 0 if payload["status"] == "Success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
