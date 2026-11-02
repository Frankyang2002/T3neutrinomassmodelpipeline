"""Source-resolved, normalization-*free* LLSS -> LLHH quartic kernel.

Read the actual Chat-1 16-model provisional F-tree LLSS release, without
assigning GammaCC/Weyl normalization or constructing an uncertified quartic.

For master convention V=lambda_abcd phi_a phi_b phi_c phi_d /4!, and
L5=C_ijab psi_i psi_j phi_a phi_b + hc, the scalar-pair contribution is
  (16 pi^2) d C[nu_p,nu_q,r,r]/d ln(mu) = sum_ab lambda[r,r,a,b] C[nu_p,nu_q,a,b].
Consequently its *conditional* C5 pivot is
  beta_C5[p,q] = sum_{a<=b} K[p,q,a,b] lambda[r,r,a,b],
  K[p,q,a,b] = 4 (2-delta_ab) C[nu_p,nu_q,a,b].

All K entries use *provisional* source coefficients: their common unknown
Matchete GammaCC -> Weyl factor/sign multiplies the entire kernel.  The
result is NOT an independently Weinberg-projected physical anomalous dimension:
there is no quartic tensor supplied and the complete SU(2) Ward orbit cannot
be checked.  No scalar matching, additional operator classes or physical C5
are inferred.  This is a diagnostic/handoff, not a certified numerical C5.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
from math import isfinite
from pathlib import Path
from typing import Any, Mapping
import json

import numpy as np

MODELS = tuple(f"{family}_{a:+d}" for family, alphas in {
    "A": (-4, -2, 0), "B": (-3, -1, 1), "C": (-3, -1, 1),
    "D": (-2, 0, 2), "E": (-4, -2, 0, 2),
}.items() for a in alphas)
SOURCE_FORMAT = "t3-lo-fermion-boundary-source-resolved-provisional-v1"
RESULT_FORMAT = "t3-eft1-conditional-quartic-pivot-kernel-v1"
EXPECTED_INDEX_ORDER = ["lepton_su2_i", "lepton_flavour_p", "lepton_su2_j",
                        "lepton_flavour_q", "real_scalar_a", "real_scalar_b"]


class SourceKernelError(ValueError):
    pass


def tag(model: str) -> str:
    if model not in MODELS:
        raise SourceKernelError(f"Unknown T3 model {model!r}")
    return model.replace("+", "p").replace("-", "m")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise SourceKernelError(message)


def _orbit(k: tuple[int, ...]) -> set[tuple[int, ...]]:
    i, p, j, q, a, b = k
    return {(i, p, j, q, a, b), (j, q, i, p, a, b),
            (i, p, j, q, b, a), (j, q, i, p, b, a)}


def _parse_coordinates(rows: Any, ns: int, label: str) -> dict[tuple[int, ...], complex]:
    _require(isinstance(rows, list), f"{label}: expected coordinate list")
    assigned: dict[tuple[int, ...], complex] = {}
    independent = set()
    for n, row in enumerate(rows):
        _require(isinstance(row, dict), f"{label}[{n}]: expected object")
        raw = row.get("coordinate")
        _require(isinstance(raw, list) and len(raw) == 6 and
                 all(type(v) is int for v in raw), f"{label}[{n}]: invalid coordinate")
        k = tuple(raw)
        i, p, j, q, a, b = k
        _require(i in (1, 2) and j in (1, 2) and
                 1 <= p <= 3 and 1 <= q <= 3 and 5 <= a <= ns and 5 <= b <= ns,
                 f"{label}[{n}]: out-of-bounds or non-LLSS coordinate {k}")
        orbit = _orbit(k)
        _require(k == min(orbit), f"{label}[{n}]: noncanonical orbit representative")
        _require(row.get("orbit_size") == len(orbit),
                 f"{label}[{n}]: wrong orbit multiplicity")
        _require(k not in independent, f"{label}[{n}]: duplicate independent representative")
        independent.add(k)
        try:
            z = complex(float(row["coefficient_real"]), float(row["coefficient_imag"]))
        except (TypeError, ValueError, KeyError, OverflowError) as exc:
            raise SourceKernelError(f"{label}[{n}]: invalid complex coefficient") from exc
        _require(isfinite(z.real) and isfinite(z.imag), f"{label}[{n}]: nonfinite coefficient")
        for member in orbit:
            _require(member not in assigned, f"{label}[{n}]: overlapping representative orbit")
            assigned[member] = z
    return assigned


def _matrix(z: np.ndarray) -> dict[str, list]:
    return {"real": z.real.tolist(), "imag": z.imag.tolist()}


def _combine(sources: list[dict[tuple[int, ...], complex]]) -> dict[tuple[int, ...], complex]:
    total = defaultdict(complex)
    for source in sources:
        for k, v in source.items():
            total[k] += v
    return dict(total)


def _compare(actual: Mapping, expected: Mapping, label: str) -> None:
    keys = set(actual) | set(expected)
    for k in keys:
        x, y = actual.get(k, 0j), expected.get(k, 0j)
        if not np.isclose(x, y, rtol=1e-9, atol=5e-26):
            raise SourceKernelError(f"{label}: independent sum mismatch at {k}: {x} != {y}")


def _make_kernels(source: Mapping[tuple[int, ...], complex], ns: int) -> list[dict]:
    """Each entry is a separate 3x3 coefficient multiplying lambda[H0R,H0R,a,b]."""
    out = []
    for a in range(5, ns + 1):
        for b in range(a, ns + 1):
            z = np.zeros((3, 3), dtype=complex)
            for p in range(1, 4):
                for q in range(1, 4):
                    z[p - 1, q - 1] = 4 * (1 if a == b else 2) * source.get((1,p,1,q,a,b), 0j)
            _require(np.allclose(z, z.T, rtol=1e-12, atol=1e-27),
                     f"Kernel for real-scalar pair ({a},{b}) is not flavour symmetric")
            if np.any(z):
                out.append({"heavy_scalar_pair_one_based": [a, b],
                            "multiplicity_in_ordered_ab_sum": 1 if a == b else 2,
                            "indeterminate_quartic": f"lambda_[H0R,H0R,{a},{b}]",
                            "conditional_beta_C5_per_unit_quartic": _matrix(z)})
    return out


def kernel_for_boundary(payload: Mapping[str, Any], *, source_sha256: str | None = None) -> dict:
    """No normalization or scalar couplings are introduced. Never labels result certified."""
    _require(isinstance(payload, Mapping) and payload.get("format") == SOURCE_FORMAT,
             "Requires Chat 1's source-resolved *provisional* boundary format")
    model = payload.get("model")
    tag(model)  # validates known model
    _require(payload.get("status") == "NumericalLLSSBySourceVerified_SpinorNormalizationOpen",
             f"{model}: unexpected source verification status")
    _require(payload.get("units") == "GeV^-1" and payload.get("index_base") == 1 and
             payload.get("index_order") == EXPECTED_INDEX_ORDER,
             f"{model}: changed index/units convention")
    spinor = payload.get("spinor_normalization", {})
    _require(isinstance(spinor, dict) and spinor.get("certified") is False and
             spinor.get("gamma_cc_to_master_weyl_factor") is None and
             spinor.get("gamma_cc_to_master_weyl_sign") is None and
             spinor.get("do_not_use_as_physical_master_wilson_boundary") is True,
             f"{model}: source normalization flags changed; inspect new release")
    completeness = payload.get("completeness", {})
    _require(completeness.get("F_tree_LLSS_included") is True and
             completeness.get("distinct_PL_sources_retained") is True and
             completeness.get("PR_Hermitian_conjugate_not_double_counted") is True and
             completeness.get("complete_physical_leading_one_loop_neutrino_mass_certified") is False,
             f"{model}: unexpected F-first completeness metadata")
    ns = payload.get("real_scalar_dimension")
    _require(type(ns) is int and 5 <= ns <= 64, f"{model}: invalid scalar dimension")
    _require(isinstance(payload.get("scalar_field_offsets_one_based"), dict) and
             isinstance(payload.get("scalar_field_complex_dimensions"), dict),
             f"{model}: missing scalar component layout")
    mu_f = payload.get("mu_F_gev")
    _require(isinstance(mu_f, (float, int)) and isfinite(mu_f) and mu_f > 0,
             f"{model}: invalid mu_F")
    raw_sources = payload.get("sources")
    _require(isinstance(raw_sources, list) and len(raw_sources) > 0,
             f"{model}: no F-tree sources")
    ids = set()
    source_maps: list[dict[tuple[int, ...], complex]] = []
    per_source = []
    for index, item in enumerate(raw_sources):
        _require(isinstance(item, dict) and
                 type(item.get("matchete_source_index")) is int and
                 item["matchete_source_index"] > 0 and
                 item["matchete_source_index"] not in ids and
                 item.get("native_chirality") == "PL",
                 f"{model}: invalid/duplicate PL source {index}")
        source_id = item["matchete_source_index"]
        ids.add(source_id)
        coords = _parse_coordinates(item.get("independent_symmetric_coordinates"), ns,
                                    f"{model}/source_{source_id}")
        if item.get("zero_at_this_benchmark") is True:
            _require(not coords, f"{model}/source_{source_id}: zero source has nonzero tensor")
        source_maps.append(coords)
        per_source.append({"matchete_source_index": source_id,
                           "zero_at_this_benchmark": not bool(coords),
                           "component_coordinate_count": len(coords),
                           "kernel_channels": _make_kernels(coords, ns)})
    total = _combine(source_maps)
    aggregate = _parse_coordinates(payload.get("summed_independent_symmetric_coordinates"),
                                   ns, f"{model}/source_sum")
    _compare(total, aggregate, f"{model}: all-F-PL source sum")
    by_pair = _make_kernels(total, ns)
    # Verify linearity using the actual complete flavour matrices for each pair.
    for pair in by_pair:
        key = pair["heavy_scalar_pair_one_based"]
        subtotal = np.zeros((3, 3), dtype=complex)
        for source in per_source:
            match = next((r for r in source["kernel_channels"]
                          if r["heavy_scalar_pair_one_based"] == key), None)
            if match is not None:
                subtotal += np.array(match["conditional_beta_C5_per_unit_quartic"]["real"]) + 1j * np.array(match["conditional_beta_C5_per_unit_quartic"]["imag"])
        target = np.array(pair["conditional_beta_C5_per_unit_quartic"]["real"]) + 1j * np.array(pair["conditional_beta_C5_per_unit_quartic"]["imag"])
        if not np.allclose(subtotal, target, rtol=1e-9, atol=5e-26):
            raise SourceKernelError(f"{model}: source-resolved scalar-pair mixing fails recomposition for pair {key}")
    return {
        "format": RESULT_FORMAT, "model": model,
        "status": "ProvisionalLinearQuarticKernel_NotPhysicalC5",
        "units": "GeV^-1 per unit dimensionless quartic, before loop factor",
        "mu_F_gev": float(mu_f), "source_boundary_sha256": source_sha256,
        "real_scalar_dimension": ns,
        "H0R_real_scalar_index_one_based": 3,
        "per_source": per_source,
        "summed_kernel_channels": by_pair,
        "source_count": len(per_source),
        "distinct_nonzero_quartic_pairs": len(by_pair),
        "conditional_beta_formula": "beta_C5_16pi2 = SUM_(a<=b) K[a,b] * lambda_[H0R,H0R,a,b] * N_GammaCC_to_Weyl",
        "conditional_delta_formula": "delta_C5_log = beta_C5_16pi2 * ln(mu_S/mu_F)/(16*pi^2)",
        "physical_normalization": {"certified": False, "N_GammaCC_to_Weyl": None},
        "quartic_component_couplings": "UNSUPPLIED - no numerical C5 is computed",
        "full_weinberg_SU2_Ward_check": "NOT POSSIBLE without full lambda_[H_A,H_B,a,b] tensor",
        "other_dim5_operator_classes": "NOT AUDITED - no physical completeness assertion",
        "finite_joint_scalar_matching": "NOT INCLUDED",
        "provenance": {
            "source_format": SOURCE_FORMAT,
            "chiral_sources": "PL only, no PR Hermitian conjugate added",
            "orbit_symmetrization": "distinct members only; no post-symmetry extra factor",
            "normalization_status": "the mathematical 4*(2-delta_ab) is canonical; unknown Matchete-to-Weyl factor remains external",
        },
    }


def load_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def file_digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()
