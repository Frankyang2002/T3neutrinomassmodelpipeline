from __future__ import annotations

"""Exact SU(2) group-factor diagnostics for the T3 Yukawa vertices.

Input
-----
A UV CG registry JSON produced by ``tools/export_uv_cg_registry.wl``.

Physics
-------
The Matchete invariant tensor is first contracted exactly.  Matchete's
``InvariantTensors`` normalization is not, in general, the orthonormal
Clebsch-Gordan normalization used to expose the representation-theory group
factor.  We therefore report both:

1. ``raw_cg_norm_factor``:
       N_raw = Tr(M_raw) / d_S

2. ``canonical_cg_norm_factor``:
   the same contraction after rescaling the complete CG tensor so that

       sum |C|^2 = d_> = max(d_F, d_S),

   which is the standard orthonormal SU(2) Clebsch-Gordan normalization for
   2 x d_< -> d_> when d_S = d_F +/- 1.

In that canonical normalization Schur's lemma gives

    sum_(a,A) C_(a A B) C*_(a A B')
        = [max(d_F,d_S)/d_S] delta_(B B').

Hence

    G_S(d_S,d_F) = max(d_F,d_S)/d_S
                 = 1                       if d_S = d_F + 1,
                   d_F/d_S                 if d_S = d_F - 1.

No fitting to the A--E beta functions is used to obtain this expression.
"""

import argparse
import itertools
import json
from dataclasses import asdict, dataclass
from pathlib import Path
import sys

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.running.MatcheteParsing import load_cg_registry


@dataclass(frozen=True)
class YukawaScalarGroupFactor:
    cg_name: str
    dF: int
    dS: int
    original_dimensions: tuple[int, int, int]
    retained_original_positions: tuple[int, ...]
    scalar_tensor_axis: int | None
    raw_total_cg_norm: str
    raw_cg_norm_factor: str
    canonical_total_cg_norm: str
    canonical_rescaling_squared: str
    canonical_cg_norm_factor: str
    closed_form_group_factor: str
    canonical_minus_closed_form: str
    raw_contraction_matrix: list[list[str]]
    canonical_contraction_matrix: list[list[str]]
    raw_schur_identity_check: bool
    canonical_schur_identity_check: bool


def _matrix_strings(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [str(sp.simplify(matrix[i, j])) for j in range(matrix.cols)]
        for i in range(matrix.rows)
    ]


def _load_seed(path: Path) -> dict:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError("UV CG registry JSON must have a dictionary top level.")
    if payload.get("status") not in {None, "Success"}:
        raise ValueError(
            f"UV CG registry export status is {payload.get('status')!r}."
        )
    return payload


def _scalar_leg_raw_matrix(
    tensor: sp.MutableDenseNDimArray,
    *,
    original_dimensions: tuple[int, int, int],
    scalar_original_position: int = 2,
) -> tuple[sp.Matrix, tuple[int, ...], int | None]:
    retained = tuple(
        i for i, dimension in enumerate(original_dimensions)
        if int(dimension) > 1
    )
    dS = int(original_dimensions[scalar_original_position])

    if len(tensor.shape) != len(retained):
        raise ValueError(
            f"Parsed tensor rank {len(tensor.shape)} does not match the "
            f"{len(retained)} retained non-singlet legs {retained}. "
            f"Tensor shape={tensor.shape}; original dimensions="
            f"{original_dimensions}."
        )

    if dS == 1:
        total = sp.S.Zero
        for index in itertools.product(
            *[range(int(dimension)) for dimension in tensor.shape]
        ):
            total += tensor[index] * sp.conjugate(tensor[index])
        return sp.Matrix([[sp.simplify(total)]]), retained, None

    scalar_axis = retained.index(scalar_original_position)
    matrix = sp.zeros(dS, dS)
    other_axes = [
        axis for axis in range(len(tensor.shape))
        if axis != scalar_axis
    ]
    other_ranges = [
        range(int(tensor.shape[axis]))
        for axis in other_axes
    ]

    for b in range(dS):
        for bp in range(dS):
            total = sp.S.Zero
            for other_values in itertools.product(*other_ranges):
                idx = [0] * len(tensor.shape)
                idxp = [0] * len(tensor.shape)

                for axis, value in zip(other_axes, other_values):
                    idx[axis] = value
                    idxp[axis] = value

                idx[scalar_axis] = b
                idxp[scalar_axis] = bp

                total += (
                    tensor[tuple(idx)]
                    * sp.conjugate(tensor[tuple(idxp)])
                )

            matrix[b, bp] = sp.simplify(total)

    return matrix, retained, scalar_axis


def _is_schur_identity(matrix: sp.Matrix) -> bool:
    n = sp.simplify(sp.trace(matrix) / matrix.rows)
    residual = (matrix - n * sp.eye(matrix.rows)).applyfunc(sp.simplify)
    return all(entry == 0 for entry in residual)


def yukawa_scalar_group_factor(
    seed_path: Path,
    *,
    cg_name: str,
    scalar_key: str,
) -> YukawaScalarGroupFactor:
    seed = _load_seed(seed_path)
    meta = seed["metadata"]

    dF = int(meta["dF"])
    dS = int(meta[scalar_key])

    registry = load_cg_registry(seed)
    if cg_name not in registry:
        available = ", ".join(sorted(registry)) or "(none)"
        raise KeyError(
            f"CG tensor {cg_name!r} not found. Available: {available}"
        )

    tensor = registry[cg_name].tensor
    original_dimensions = (2, dF, dS)

    raw_matrix, retained, scalar_axis = _scalar_leg_raw_matrix(
        tensor,
        original_dimensions=original_dimensions,
        scalar_original_position=2,
    )

    raw_total = sp.simplify(sp.trace(raw_matrix))
    raw_factor = sp.simplify(raw_total / dS)

    # For the allowed T3 Yukawa, 2 x d_< contains d_> with
    # d_> = max(dF,dS).  Orthonormal CG coefficients obey
    # sum |C|^2 = d_>.
    canonical_total = sp.Integer(max(dF, dS))
    rescale_squared = sp.simplify(canonical_total / raw_total)
    canonical_matrix = raw_matrix.applyfunc(
        lambda entry: sp.simplify(rescale_squared * entry)
    )
    canonical_factor = sp.simplify(
        sp.trace(canonical_matrix) / dS
    )

    closed_form = sp.Rational(max(dF, dS), dS)
    residual = sp.simplify(canonical_factor - closed_form)

    return YukawaScalarGroupFactor(
        cg_name=cg_name,
        dF=dF,
        dS=dS,
        original_dimensions=original_dimensions,
        retained_original_positions=retained,
        scalar_tensor_axis=scalar_axis,
        raw_total_cg_norm=str(raw_total),
        raw_cg_norm_factor=str(raw_factor),
        canonical_total_cg_norm=str(canonical_total),
        canonical_rescaling_squared=str(rescale_squared),
        canonical_cg_norm_factor=str(canonical_factor),
        closed_form_group_factor=str(closed_form),
        canonical_minus_closed_form=str(residual),
        raw_contraction_matrix=_matrix_strings(raw_matrix),
        canonical_contraction_matrix=_matrix_strings(canonical_matrix),
        raw_schur_identity_check=_is_schur_identity(raw_matrix),
        canonical_schur_identity_check=_is_schur_identity(canonical_matrix),
    )


def y1_scalar_wavefunction_group_factor(
    seed_path: Path,
) -> YukawaScalarGroupFactor:
    return yukawa_scalar_group_factor(
        seed_path,
        cg_name="T3Y1CG",
        scalar_key="dS1",
    )


def y2_scalar_wavefunction_group_factor(
    seed_path: Path,
) -> YukawaScalarGroupFactor:
    return yukawa_scalar_group_factor(
        seed_path,
        cg_name="T3Y2CG",
        scalar_key="dS2",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Contract exact UV T3 Yukawa CGs and convert them to the "
            "orthonormal SU(2) CG convention."
        )
    )
    parser.add_argument("uv_cg_registry_json", type=Path)
    parser.add_argument(
        "--yukawa",
        choices=("y1", "y2"),
        default="y1",
    )
    args = parser.parse_args()

    result = (
        y1_scalar_wavefunction_group_factor(args.uv_cg_registry_json)
        if args.yukawa == "y1"
        else y2_scalar_wavefunction_group_factor(args.uv_cg_registry_json)
    )

    print(json.dumps(asdict(result), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
