#!/usr/bin/env python3
"""Workstream A: exact component Ward test of the 38 T3 Yukawa matrices.

Source data:
  - 16 existing native Matchete baseline CG exports (32 Yukawas), and
  - 6 NEW live Matchete extra-CG exports (independent scalar orientation).

Checks SU(2), U(1)_Y, and the trivial SU(3)_c singlet action. Compares each
additional CG with its corresponding native baseline *component by component*.
The SU(2) triplet matrices are the exact reference-coordinate matrices from
Stage27, NOT independently extracted Matchete/RGBeta native generators.

This IS a concrete component Ward test in that representation basis. It is NOT
an action-level Matchete/RGBeta coupling normalization certificate, and it
makes no changes to production RGEs, matching, or existing source results.

Commands (repository root):
  python Lagrangian/RunT3YukawaComponentWard.py --self-test
  python Lagrangian/RunT3YukawaComponentWard.py --all
  python Lagrangian/RunT3YukawaComponentWard.py --all --no-wolfram
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys

import sympy as sp
from sympy.parsing.mathematica import parse_mathematica

ROOT = Path(__file__).resolve().parent.parent
EXTRAS = {"A_m2", "A_p0", "D_m2", "D_p0", "E_m2", "E_p0"}
MODELS = {
    "A": ((1, 3, 2), (-4, -2, 0)),
    "B": ((2, 2, 1), (-3, -1, 1)),
    "C": ((2, 2, 3), (-3, -1, 1)),
    "D": ((3, 1, 2), (-2, 0, 2)),
    "E": ((3, 3, 2), (-4, -2, 0, 2)),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def model_labels() -> dict[str, tuple[int, int, int, int]]:
    result = {}
    for family, (dims, alphas) in MODELS.items():
        for alpha in alphas:
            tag = f"m{-alpha}" if alpha < 0 else f"p{alpha}"
            result[f"{family}_{tag}"] = (*dims, alpha)
    require(len(result) == 16, "16 model definitions expected")
    return result


def exact_number(raw: str) -> sp.Expr:
    """Parse exact numeric Mathematica InputForm; reject symbolic leftovers."""
    require(isinstance(raw, str) and len(raw) <= 200, "Invalid numeric component")
    allowed = set("0123456789-+*/^().[] I,Sqrt")
    require(set(raw) <= allowed, f"Unsupported component expression: {raw}")
    val = sp.simplify(parse_mathematica(raw))
    require(not val.free_symbols and val.is_number is True and val.is_finite is True,
            f"Component is not finite/exact numeric: {raw}")
    return val


def tensor_from_record(record: dict, dims: tuple[int, int, int], name: str) -> sp.MutableDenseNDimArray:
    active_dims = [d for d in dims if d != 1]
    require(record.get("dimensions") == active_dims,
            f"{name}: CG shape differs from physical dimensions")
    require(record.get("tensor_accessor") == "Matchete`GetCGTensor"
            and record.get("registered_in_assembled_interaction") is True,
            f"{name}: tensor is not proven native and registered")
    vals = [exact_number(x) for x in record.get("components_inputform", [])]
    require(len(vals) == sp.prod(dims) and any(x != 0 for x in vals),
            f"{name}: empty/wrong-length tensor")
    return sp.MutableDenseNDimArray(vals, dims)


T2 = (
    sp.Matrix([[0, 1], [1, 0]]) / 2,
    sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2,
    sp.diag(1, -1) / 2,
)
# Exact Matchete-coordinate spin-one reference candidate from Stage27.
# Basis is determined by the native symmetric doublet-triplet CG, not by
# independent numerical RGBeta components. Never label it a native extraction.
T3 = (
    sp.Matrix([[0, 0, 0], [0, 0, sp.I], [0, -sp.I, 0]]),
    sp.Matrix([[0, sp.I, 0], [-sp.I, 0, 0], [0, 0, 0]]),
    sp.Matrix([[0, 0, -sp.I], [0, 0, 0], [sp.I, 0, 0]]),
)


def native_generators(dim: int) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    require(dim in (1, 2, 3), "Unsupported SU2 dimension")
    return (sp.zeros(1),) * 3 if dim == 1 else T2 if dim == 2 else T3


def ward_residuals(c: sp.MutableDenseNDimArray, dims: tuple[int, int, int],
                   conjugated: tuple[bool, bool, bool]) -> tuple[list[sp.Expr], int]:
    """TENSOR-index Ward equation: sum_r C[...,r,...] T[r,index] = 0.

    A conjugated field transforms with -T*, not T. The transpose location
    is essential: confusing T[r,i] with T[i,r] breaks non-real generators.
    """
    gs = [tuple(-g.conjugate() if conj else g for g in native_generators(dim))
          for dim, conj in zip(dims, conjugated)]
    out = []
    for a in range(3):
        for i in range(dims[0]):
            for j in range(dims[1]):
                for k in range(dims[2]):
                    val = sum(c[r, j, k] * gs[0][a][r, i] for r in range(dims[0]))
                    val += sum(c[i, r, k] * gs[1][a][r, j] for r in range(dims[1]))
                    val += sum(c[i, j, r] * gs[2][a][r, k] for r in range(dims[2]))
                    val = sp.simplify(val)
                    if val != 0:
                        out.append(val)
    return out, 3 * int(sp.prod(dims))


def u1_charges(alpha: int, which: int, extra: bool) -> tuple[Fraction, Fraction, Fraction]:
    y1, yf, y2 = Fraction(alpha, 2), Fraction(alpha + 1, 2), Fraction(alpha + 2, 2)
    if which == 1:
        return (Fraction(1, 2), -yf, -y1 if extra else y1)
    return (Fraction(1, 2), yf, y2 if extra else -y2)


def norm_squared(c: sp.MutableDenseNDimArray) -> sp.Expr:
    return sp.simplify(sum(sp.conjugate(v) * v for v in c._array))


def check_yukawa(model: str, dims: tuple[int, int, int], alpha: int,
                 which: int, extra: bool, record: dict) -> dict:
    # Always field axes = (Bar[L], F^c/F, S1/S2). Matchete's CG record
    # omits only dimension-1 axes; tensor_from_record restores them.
    cdims = (2, dims[2], dims[which - 1])
    conjugated = (True, which == 1, (which == 1) == extra)
    # Baseline: Y1 scalar S1 (not conjugate), Y2 scalar Bar[S2].
    # Extra:    Y1 scalar Bar[S1], Y2 scalar S2.
    require(conjugated[2] == (extra if which == 1 else not extra),
            "Scalar orientation mismatch")
    name = f"{model}/y{which}{'Conj' if extra else ''}"
    c = tensor_from_record(record, cdims, name)
    residuals, equations = ward_residuals(c, cdims, conjugated)
    charges = u1_charges(alpha, which, extra)
    total_charge = sum(charges, Fraction(0))
    require(not residuals, f"{name}: nonzero SU2 Ward residual: {residuals[:4]}")
    require(total_charge == 0, f"{name}: nonzero hypercharge sum: {total_charge}")
    return {
        "model": model,
        "coupling": f"y{which}{'Conj' if extra else ''}",
        "field_order": ["Bar[L]", "F^c" if which == 1 else "F",
                        ("Bar[S]" if conjugated[2] else "S")],
        "tensor_dimensions_with_singlets": list(cdims),
        "native_cg_norm_squared_exact": str(norm_squared(c)),
        "conjugated_axes": list(conjugated),
        "u1_hypercharges": [str(q) for q in charges],
        "u1_ward_residual": str(total_charge),
        "su2_ward_equations_checked": equations,
        "su2_ward_nonzero_residuals": len(residuals),
        "su3_color_singlet_ward": "Trivial (all three project fields are SU3 singlets)",
        "passed": True,
    }


def source_hash(root: Path, relative: str) -> str:
    p = root / relative
    require(p.is_file(), f"Missing source: {p}")
    return sha256(p.read_bytes()).hexdigest()


def check_stage27(root: Path) -> None:
    """Guard the fixed triplet reference basis against a source change."""
    p = root / "output/general_t3_stage22/stage27_triplet_generator_witness.json"
    source = json.loads(p.read_text(encoding="utf-8"))
    require(source.get("status") ==
            "ExactSingleMatcheteCGSixTripletsGlobalRGBetaBasisUncertified",
            "Stage27 triplet witness unavailable/changed")
    mats = source.get("matchete_triplet_candidate_generators", [])
    require(len(mats) == 3, "No Stage27 triplet-generator candidate")
    for a in range(3):
        g = sp.Matrix([[exact_number(x) for x in row] for row in mats[a]])
        require(g == T3[a], "Triplet reference basis changed since Stage27")
    require(source.get("native_rgbeta_s2_generator_components_exported") is False and
            source.get("global_matchete_rgbeta_cartesian_basis_certified") is False,
            "This audit must not silently upgrade Stage27's status")


def build_report(root: Path) -> dict:
    check_stage27(root)
    baseline_dir = root / "Reports/output/workstream_a/native_cg_bridge/native"
    extras_dir = root / "Reports/output/workstream_a/yukawa_component_ward/native_extras"
    rows, extra_equal = [], []
    for name, (s1, s2, f, alpha) in model_labels().items():
        p = baseline_dir / f"{name}.json"
        require(p.is_file(), f"Missing registered Matchete CG: {p}")
        native = json.loads(p.read_text(encoding="utf-8"))
        require(native.get("status") == "LiveRegisteredBaselineCGsCaptured" and
                native.get("model") == name and native.get("full_matchete_check_passed") is True,
                f"{name}: previous live Matchete baseline record rejected")
        original = native.get("native_tensors", {})
        for which in (1, 2):
            key = f"Yukawa{which}"
            rows.append(check_yukawa(name, (s1, s2, f), alpha, which, False, original[key]))
        if name in EXTRAS:
            p_extra = extras_dir / f"{name}.json"
            require(p_extra.is_file(), f"Missing NEW live extra-CG export: {p_extra}")
            live = json.loads(p_extra.read_text(encoding="utf-8"))
            which = 1 if alpha == 0 else 2
            base = original[f"Yukawa{which}"]
            require(live.get("model") == name and live.get("status") ==
                    "LiveExtraCGEqualityVerified" and
                    live.get("extra_yukawa") == f"y{which}Conj" and
                    live.get("full_matchete_check_passed") is True and
                    live.get("extra_equals_baseline_components_exactly") is True,
                    f"{name}: extra CG is not a verified native export")
            new = live["extra_tensor"]
            ref = live["baseline_tensor"]
            active = tuple(x for x in (2, f, s1 if which == 1 else s2) if x != 1)
            for rec in (new, ref, base):
                require(rec.get("dimensions") == list(active), f"{name}: CG dims changed")
            for x, y in zip((new, ref), (base, base)):
                require(len(x["components_inputform"]) == len(y["components_inputform"]) and
                        all(sp.simplify(exact_number(a) - exact_number(b)) == 0
                            for a, b in zip(x["components_inputform"], y["components_inputform"])),
                        f"{name}: live extra/baseline native CG mismatch")
            rows.append(check_yukawa(name, (s1, s2, f), alpha, which, True, new))
            extra_equal.append(name)
    require(len(rows) == 38 and len(extra_equal) == 6, "Missing Yukawa checks")
    return {
        "format": "t3-workstream-a-exact-yukawa-component-ward-v1",
        "status": "ExactSU2U1ComponentWardPassedPhysicalCrossPackageDictionaryPending",
        "models_checked": 16,
        "baseline_yukawas_checked": 32,
        "new_native_extra_yukawas_checked": 6,
        "extra_native_cgs_equal_to_baseline_exact": extra_equal,
        "su2_ward_component_equations_checked": sum(x["su2_ward_equations_checked"] for x in rows),
        "u1_hypercharge_checks": len(rows),
        "su3_c": "Trivial on the model-source color-singlet fields; no colored external indices",
        "triplet_representation_basis": "Stage27 exact CG-derived Matchete reference candidate",
        "independent_native_triplet_generator_coordinates_verified": False,
        "rgbeta_native_ward_check_passed": False,
        "physical_matchete_rgbeta_yukawa_normalization_certified": False,
        "physical_action_equivalence_verified": False,
        "production_rgbeta_modified": False,
        "next_step": "Reconstruct both packages' physical Yukawa actions using the Weyl/Dirac and complex-real scalar conventions; compare exact signed coefficients.",
        "source_sha256": {
            p: source_hash(root, p) for p in (
                "Lagrangian/model/SU2Invariants.wl",
                "Lagrangian/interactions/T3Topology.wl",
                "RGE/running/rgbeta/T3RGBetaModel.wl",
                "output/general_t3_stage22/stage27_triplet_generator_witness.json",
            )
        },
        "yukawas": rows,
    }


def live_export(root: Path, timeout: int) -> None:
    out_dir = root / "Reports/output/workstream_a/yukawa_component_ward/native_extras"
    out_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update({
        "T3_EXPERIMENTAL_GENERAL_FERMIONS": "1",
        "T3_EXPERIMENTAL_FLAVOR_MASSES": "1",
        "T3_EXPERIMENTAL_GENERAL_QUARTICS": "1",
        "T3_EXPERIMENTAL_GENERAL_QUADRATICS": "0",
    })
    wl = root / "Lagrangian/RunT3YukawaExtraCGExport.wl"
    require(wl.is_file(), f"Missing export script: {wl}")
    for name in sorted(EXTRAS):
        cls, tag = name.split("_")
        dest = out_dir / f"{name}.json"
        proc = subprocess.run(
            ["wolframscript", "-file", str(wl), str(dest), cls, tag],
            env=env, cwd=root, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, timeout=timeout, check=False,
        )
        (out_dir / f"{name}.log").write_text(proc.stdout, encoding="utf-8")
        require(proc.returncode == 0, f"{name} live Matchete CG export failed (exit={proc.returncode}); read {name}.log")
        print(f"Live extra CG: {name} PASS", flush=True)


def self_test() -> None:
    e = {"tensor_accessor": "Matchete`GetCGTensor",
         "registered_in_assembled_interaction": True, "dimensions": [2, 2],
         "components_inputform": ["0", "1", "-1", "0"]}
    ident = {**e, "components_inputform": ["1", "0", "0", "1"]}
    r1 = check_yukawa("A_p0", (1, 3, 2), 0, 1, True, e)
    r2 = check_yukawa("D_m2", (3, 1, 2), -2, 2, True, ident)
    require(r1["passed"] and r2["passed"], "Independent singlet fixture failed")
    tampered = {**e, "components_inputform": ["1", "1", "-1", "0"]}
    try:
        check_yukawa("A_p0", (1, 3, 2), 0, 1, True, tampered)
    except ValueError as exc:
        require("SU2 Ward residual" in str(exc), "Wrong failure for corrupted CG")
    else:
        raise ValueError("Negative control incorrectly accepted corrupted CG")
    print("PASS: two exact singlet component Ward fixtures and one nonzero-residual negative control")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--all", action="store_true", help="Export six native Matchete CGs, then audit 38 Yukawas")
    parser.add_argument("--no-wolfram", action="store_true", help="Use previously exported native extra CG files")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--timeout", type=int, default=600, help="Seconds per Wolfram export")
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            return 0
        require(args.all, "Use --all or --self-test")
        root = args.repo_root.resolve()
        if not args.no_wolfram:
            live_export(root, args.timeout)
        report = build_report(root)
        dest = root / "Reports/output/workstream_a/yukawa_component_ward/yukawa_component_ward.json"
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"{report['status']}: {len(report['yukawas'])}/38 Yukawas; "
              f"{report['su2_ward_component_equations_checked']} exact SU(2) equations")
        print("Cross-package physical coupling normalization: NOT CERTIFIED")
        print(f"Report: {dest}")
        return 0
    except (ValueError, OSError, json.JSONDecodeError, subprocess.TimeoutExpired) as exc:
        print("Yukawa component Ward audit BLOCKED: " + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
