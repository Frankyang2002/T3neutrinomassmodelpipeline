#!/usr/bin/env python3
"""Read-only A(0) Yukawa action bridge, including independent y1Conj.

Uses existing *live* Matchete CG exports, the existing RGBeta AddYukawa
registry probe, and the 38-Yukawa normalization and Ward audit reports.
Does NOT execute Wolfram or edit any production code.

The two physical Yukawa actions are expanded in the common formal monomial
basis (S1 or S1†) x (l F or its Hermitian conjugate) x SU2 indices.
The Matchete four-component bilinear is converted with a declared *unknown*
unit-modulus Weyl phase eta. The absolute orientation of RGBeta's symbolic
SU2 epsilon is also not independently numerically exported; rho=+/-1 is
checked as a convention choice. Thus action equality and coupling map are
CONDITIONAL, not a native cross-package physical certificate.

Usage (repo root):
  python Lagrangian/RunT3YukawaA0ActionBridge.py --self-test
  python Lagrangian/RunT3YukawaA0ActionBridge.py --all

Output:
  Reports/output/workstream_a/yukawa_a0_action_bridge/a0_action_bridge.json
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT_REL = Path("Reports/output/workstream_a/yukawa_a0_action_bridge/a0_action_bridge.json")
BASE_REL = Path("Reports/output/workstream_a/native_cg_bridge/native/A_p0.json")
EXTRA_REL = Path("Reports/output/workstream_a/yukawa_component_ward/native_extras/A_p0.json")
EXTRA_RGBETA_REL = Path("Reports/output/workstream_a/rgbeta_extra_yukawa_bridge/native_rgbeta/A_p0.json")
PROBE_REL = Path("Reports/output/workstream_a/yukawa_normalization_bridge/native_rgbeta/A_p0.json")
NORM_REL = Path("Reports/output/workstream_a/yukawa_normalization_bridge/yukawa_normalization_bridge.json")
WARD_REL = Path("Reports/output/workstream_a/yukawa_component_ward/yukawa_component_ward.json")
MAT_SOURCE = Path("Lagrangian/interactions/T3Topology.wl")
RGB_SOURCE = Path("RGE/running/rgbeta/T3RGBetaModel.wl")


class BridgeError(ValueError):
    pass


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise BridgeError(reason)


def load(root: Path, path: Path) -> dict:
    filename = root / path
    require(filename.is_file(), f"Missing prior native evidence: {filename}")
    obj = json.loads(filename.read_text(encoding="utf-8"))
    require(isinstance(obj, dict), f"Malformed report: {filename}")
    return obj


def exact_cg(record: dict) -> sp.Matrix:
    require(record.get("tensor_accessor") == "Matchete`GetCGTensor"
            and record.get("registered_in_assembled_interaction") is True
            and record.get("dimensions") == [2, 2],
            "A(0) tensor is not an actual 2x2 registered native Matchete CG")
    components = record.get("components_inputform")
    require(isinstance(components, list) and len(components) == 4,
            "A(0): wrong number of CG components")
    # A(0) native singlet tensor is real integer. Refuse arbitrary Wolfram
    # code, complex numerical approximations, or nontrivial expressions here.
    require(all(value in ("0", "1", "-1") for value in components),
            "A(0): native CG has unexpected noninteger entries")
    return sp.Matrix(2, 2, [sp.Integer(v) for v in components])


def exact_equal(x: sp.Expr, y: sp.Expr) -> bool:
    return sp.simplify(x - y) == 0


def source_guards(root: Path) -> dict:
    mat = root / MAT_SOURCE
    rgb = root / RGB_SOURCE
    for source in (mat, rgb):
        require(source.is_file(), f"Missing source: {source}")
    ms = mat.read_text(encoding="utf-8")
    rs = rgb.read_text(encoding="utf-8")
    # Read-only source guards: inspect exactly which field orientations and
    # chirality were installed; do not treat old reports as timeless proof.
    for needle in (
        "scalarBar = which === 2;", "useCConj = which === 1;",
        '"Expression" -> PlusHc[',
        "y[p, r] NCM[Bar[l[i, p]], f] scalar CG[cg, labels]",
        "conjugated = which === 1;",
        "coupling[p, r] NCM[Bar[l[i, p]], f] scalar CG[cg, labels]",
    ):
        require(needle in ms, f"Matchete source convention changed: {needle}")
    for needle in (
        "T3RGBetaYukawaInvariant[1, 2]", "eps[SU2L @ fund, #2, #3]",
        "{Bar @ S1, l, F}", "Chirality -> Right",
    ):
        require(needle in rs, f"RGBeta source convention changed: {needle}")
    return {str(path): sha256((root / path).read_bytes()).hexdigest()
            for path in (MAT_SOURCE, RGB_SOURCE)}


def verify_evidence(root: Path) -> dict:
    source = source_guards(root)
    base = load(root, BASE_REL)
    extra = load(root, EXTRA_REL)
    reg = load(root, EXTRA_RGBETA_REL)
    probe = load(root, PROBE_REL)
    norm = load(root, NORM_REL)
    ward = load(root, WARD_REL)

    require(base.get("model") == "A_p0"
            and base.get("status") == "LiveRegisteredBaselineCGsCaptured"
            and base.get("dimensions") == [1, 3, 2]
            and base.get("alpha") == 0
            and base.get("full_matchete_check_passed") is True
            and base.get("extra_yukawas_included") == 1,
            "Missing accepted A(0) Matchete action/CG source evidence")
    require(extra.get("model") == "A_p0"
            and extra.get("status") == "LiveExtraCGEqualityVerified"
            and extra.get("extra_yukawa") == "y1Conj"
            and extra.get("extra_equals_baseline_components_exactly") is True
            and extra.get("full_matchete_check_passed") is True,
            "A(0) live Matchete y1Conj CG report is missing/uncertified")
    baseline = exact_cg(base["native_tensors"]["Yukawa1"])
    extra_base = exact_cg(extra["baseline_tensor"])
    extra_actual = exact_cg(extra["extra_tensor"])
    require(baseline == extra_base == extra_actual,
            "A(0) actual baseline/extra Matchete CGs do not agree")
    epsilon = sp.Matrix([[0, 1], [-1, 0]])
    require(baseline == epsilon,
            "A(0) native Matchete CG not in expected epsilon coordinates")

    require(reg.get("status") == "ExtraYukawaRegisteredNativeWardUncertified"
            and reg.get("model") == "A_p0"
            and reg.get("extra_yukawa") == "y1Conj"
            and reg.get("extra_fields_inputform") == "{S1, l, F}"
            and reg.get("baseline_fields_inputform") == "{Bar[S1], l, F}"
            and reg.get("extra_invariant_inputform") ==
                "eps[SU2L[fund], t3ry2, t3ry3]"
            and reg.get("baseline_invariant_inputform") ==
                "eps[SU2L[fund], t3ry2, t3ry3]"
            and reg.get("same_native_su2_invariant_as_baseline") is True,
            "A(0) native RGBeta invariant or scalar orientation changed")

    require(probe.get("status") == "NativeRGBetaYukawaInternalPrefactorsVerified"
            and probe.get("model") == "A_p0"
            and probe.get("internal_real_scalar_tensor_prefactor") == "Sqrt[2]"
            and probe.get("physical_action_equivalence_verified") is False,
            "A(0) native RGBeta prefactor probe unavailable/changed")
    rows = probe.get("baseline_yukawas", [])
    extra_row = probe.get("extra_yukawa")
    require(isinstance(rows, list) and len(rows) == 2 and
            isinstance(extra_row, dict), "A(0) native RGBeta Yukawa rows missing")
    checked = [r for r in rows if r.get("coupling") == "y1"] + [extra_row]
    require([r.get("coupling") for r in checked] == ["y1", "y1Conj"],
            "A(0) y1/y1Conj not independently registered")
    for r in checked:
        require(r.get("chirality") == "Right"
                and r.get("native_complex_scalar_yukawa_tensor_prefactor_verified") is True
                and r.get("coefficient_of_Matrix_Bar_y") in
                    ("Sqrt[2]", "2^Rational[1, 2]")
                and r.get("coefficient_of_Matrix_y_in_CouplingBar") in
                    ("Sqrt[2]", "2^Rational[1, 2]"),
                "A(0) chirality or native RGBeta coefficient changed")
    require(checked[0].get("fields_inputform") == "{Bar[S1], l, F}"
            and checked[1].get("fields_inputform") == "{S1, l, F}",
            "A(0) probe native scalar orientation changed")
    require(norm.get("status") ==
                "ExactCGMagnitudesAndNativeRGBetaInternalPrefactorsVerifiedActionDictionaryPending"
            and norm.get("native_matchete_yukawas_checked") == 38
            and norm.get("rgbeta_native_probe_models") == 8
            and norm.get("matchete_rgbeta_absolute_coupling_dictionary_certified") is False,
            "Complete 38-Yukawa normalization audit missing/changed")
    require(ward.get("status") ==
                "ExactSU2U1ComponentWardPassedPhysicalCrossPackageDictionaryPending"
            and ward.get("su2_ward_component_equations_checked") == 1032
            and ward.get("new_native_extra_yukawas_checked") == 6,
            "Complete component Ward audit missing/changed")
    return {"cg_matchete": baseline, "source_sha256": source,
            "prior_ward_su2_equations": 1032,
            "prior_38_yukawa_norms": 38}


# A coordinate monomial is (scalar= S or Sbar, bilinear= W or Wbar,
# lepton_SU2_index, F_SU2_index). W = (l_i F_j) Weyl bilinear.
# The Q_M = NCM[Bar[L_i], CConj[F_j]] -> eta * Wbar mapping is an explicit
# hypothesis; eta's complex conjugate is eta^-1 because |eta|=1.
Monomial = tuple[str, str, int, int]


def expand_matchete(c: sp.Matrix, *, eta: sp.Expr) -> dict[Monomial, sp.Expr]:
    y, yc, yb, ycb = sp.symbols("yM ycM yMbar ycMbar")
    coeff: dict[Monomial, sp.Expr] = {}
    for i in range(2):
        for j in range(2):
            coeff[("S", "Wbar", i, j)] = eta * y * c[i, j]
            coeff[("Sbar", "W", i, j)] = eta**-1 * yb * sp.conjugate(c[i, j])
            coeff[("Sbar", "Wbar", i, j)] = eta * yc * c[i, j]
            coeff[("S", "W", i, j)] = eta**-1 * ycb * sp.conjugate(c[i, j])
    return coeff


def expand_rgbeta(eps: sp.Matrix, *, ry: sp.Expr, rybar: sp.Expr,
                  ryc: sp.Expr, rycbar: sp.Expr) -> dict[Monomial, sp.Expr]:
    """RGBeta MBL L=-y phi psi psi+HC, with Chirality=Right.

    sqrt2 native Yuk tensor x 1/sqrt2 complex-real scalar projection,
    x two (l,F) permutations x (-1/2) gives -1 per physical monomial.
    """
    coeff: dict[Monomial, sp.Expr] = {}
    factor = -sp.Rational(1, 2) * 2 * (sp.sqrt(2) / sp.sqrt(2))
    require(factor == -1, "RGBeta real/complex scalar symmetry factor changed")
    for i in range(2):
        for j in range(2):
            coeff[("Sbar", "W", i, j)] = factor * rybar * eps[i, j]
            coeff[("S", "Wbar", i, j)] = factor * ry * sp.conjugate(eps[i, j])
            coeff[("S", "W", i, j)] = factor * rycbar * eps[i, j]
            coeff[("Sbar", "Wbar", i, j)] = factor * ryc * sp.conjugate(eps[i, j])
    return coeff


def compare_actions(m: dict[Monomial, sp.Expr],
                    r: dict[Monomial, sp.Expr]) -> int:
    require(m.keys() == r.keys() and len(m) == 16,
            "Missing scalar/bilinear/SU2 component monomials")
    nonzero = 0
    for key in sorted(m):
        require(exact_equal(m[key], r[key]),
                f"Action mismatch in {key}: Matchete={m[key]}, RGBeta={r[key]}")
        nonzero += int(m[key] != 0)
    return nonzero


def solve_action(c: sp.Matrix, rho: int) -> dict:
    require(rho in (-1, 1), "Epsilon orientation must be +/-1")
    eps = rho * sp.Matrix([[0, 1], [-1, 0]])
    # A common SU(2) field basis is DECLARED, not proved from RGBeta numeric
    # indices; rho represents the remaining epsilon orientation choice.
    ratios = []
    for i in range(2):
        for j in range(2):
            if eps[i, j] != 0:
                ratios.append(sp.simplify(c[i, j] / sp.conjugate(eps[i, j])))
            else:
                require(c[i, j] == 0, "Matchete CG has support outside RGBeta epsilon")
    require(len(ratios) == 2 and all(exact_equal(x, ratios[0]) for x in ratios),
            "No single epsilon-proportional action coefficient")
    q = ratios[0]
    eta = sp.Symbol("eta", nonzero=True)
    y, yc, yb, ycb = sp.symbols("yM ycM yMbar ycMbar")
    # Both independent Yukawas must share the same spinor conversion and SU2
    # reference basis. Their physical RGBeta couplings are NOT H.c. of each
    # other; only ybar is the conjugate of y, and ycbar of yc.
    ry = -eta * q * y
    ryc = -eta * q * yc
    rybar = -eta**-1 * sp.conjugate(q) * yb
    rycbar = -eta**-1 * sp.conjugate(q) * ycb
    m = expand_matchete(c, eta=eta)
    r = expand_rgbeta(eps, ry=ry, rybar=rybar,
                      ryc=ryc, rycbar=rycbar)
    nonzero = compare_actions(m, r)
    require(nonzero == 8, "Expected 8 nonzero SU2-spinor-scalar monomials")

    # Cross-check independent y1/y1Conj: they must appear with distinct
    # scalar orientations and are not Hermitian conjugates of one another.
    require(sp.expand(m[("S", "Wbar", 0, 1)]).has(y)
            and sp.expand(m[("Sbar", "Wbar", 0, 1)]).has(yc),
            "Independent Yukawa/source scalar orientation mixed")
    return {
        "epsilon_relative_orientation_rho": rho,
        "matchete_to_rgbeta_cg_ratio_q": str(q),
        "physical_candidate_mapping_y1_R": str(sp.simplify(ry)),
        "physical_candidate_mapping_y1Conj_R": str(sp.simplify(ryc)),
        "conjugate_candidate_mapping_y1_Rbar": str(sp.simplify(rybar)),
        "conjugate_candidate_mapping_y1Conj_Rbar": str(sp.simplify(rycbar)),
        "all_monomial_channels_compared": len(m),
        "nonzero_exact_monomial_coefficients_equal": nonzero,
        "all_residuals_zero_exact": True,
        "conditional_on_common_su2_coordinates": True,
        "conditional_on_undetermined_weyl_bilinear_phase_eta": True,
        "physical_absolute_mapping_certified": False,
    }


def self_test() -> None:
    c = sp.Matrix([[0, 1], [-1, 0]])
    for rho in (-1, 1):
        result = solve_action(c, rho)
        require(result["nonzero_exact_monomial_coefficients_equal"] == 8,
                "Formal A(0) H.c. test failed")
    # Wrong CG support and mismatched antisymmetric orientation must fail.
    for corrupt in (sp.eye(2), sp.Matrix([[0, 1], [1, 0]])):
        try:
            solve_action(corrupt, 1)
        except BridgeError:
            pass
        else:
            raise BridgeError("Negative control accepted an invalid CG")
    # Check that the opposite-scalar coupling cannot silently be identified
    # with the baseline coupling in the independent symbolic M action.
    eta = sp.Symbol("eta", nonzero=True)
    terms = expand_matchete(c, eta=eta)
    require(terms[("Sbar", "Wbar", 0, 1)].has(sp.Symbol("ycM")) and
            not terms[("Sbar", "Wbar", 0, 1)].has(sp.Symbol("yM")),
            "Independent extra Yukawa collapsed into baseline")
    # Negative controls catch a missing H.c. channel and an incorrect
    # extra-Yukawa coefficient even with the same SU2 CG tensor.
    rgb = dict(terms)
    rgb.pop(("Sbar", "W", 0, 1))
    try:
        compare_actions(terms, rgb)
    except BridgeError:
        pass
    else:
        raise BridgeError("Negative control accepted missing Hermitian conjugate")
    rgb = dict(terms)
    rgb[("Sbar", "Wbar", 0, 1)] *= -1
    try:
        compare_actions(terms, rgb)
    except BridgeError:
        pass
    else:
        raise BridgeError("Negative control accepted wrong extra-Yukawa sign")
    print("PASS: both epsilon orientations, 16 formal channels, Hermitian conjugates, negative controls")


def run(root: Path) -> dict:
    evidence = verify_evidence(root)
    results = [solve_action(evidence["cg_matchete"], rho) for rho in (1, -1)]
    return {
        "format": "t3-workstream-a-a0-conditional-yukawa-action-v1",
        "model": "A_p0",
        "status": "A0ExactFormalWeylActionIdentityBothIndependentYukawasPhysicalPhasePending",
        "source": "Native Matchete y1/y1Conj CG, live RGBeta AddYukawa registry and prefactor, RGBeta MBL Yukawa sign",
        "previous_38_yukawa_cg_and_prefactor_audit": "Passed",
        "previous_1032_component_ward_equations": "Passed",
        "matchete_cg_dimensions": [2, 2],
        "matchete_cg_components": ["0", "1", "-1", "0"],
        "matchete_nonhc_terms": [
            "y1M S1 NCM[Bar[L],CConj[F]] epsilon",
            "y1ConjM Bar[S1] NCM[Bar[L],CConj[F]] epsilon",
        ],
        "rgbeta_native_nonhc_terms": [
            "-Bar[y1R] Bar[S1] (l F) epsilon",
            "-Bar[y1ConjR] S1 (l F) epsilon",
        ],
        "rgbeta_generic_yukawa_tensor_factor": "2/sqrt(2)",
        "rgbeta_complex_scalar_projection": "1/sqrt(2)",
        "rgbeta_distinct_fermion_permutations": 2,
        "rgbeta_generic_lagrangian_yukawa_prefactor": "-1/2",
        "rgbeta_physical_nonhc_coefficient": "-1",
        "physical_fermion_mapping_assumption":
            "NCM[Bar[L_i],CConj[F_j]] = eta * (l_i F_j)^dagger, |eta|=1, in a declared common spinor convention",
        "native_matchete_spinor_conversion_eta_measured": False,
        "native_rgbeta_absolute_epsilon_coordinates_measured": False,
        "conditional_candidate_maps": results,
        "relative_y1Conj_to_y1_mapping_coherent": True,
        "formal_action_identity_under_declared_assumptions": True,
        "absolute_matchete_rgbeta_yukawa_coupling_mapping_certified": False,
        "physical_cross_package_action_certified": False,
        "production_files_modified": False,
        "next_required_measurement":
            "Anchor the four-component Matchete CConj bilinear to the RGBeta left-Weyl bilinear in an explicit canonical spinor convention, then independently fix the RGBeta fundamental epsilon orientation in the same physical component basis.",
        "source_sha256": evidence["source_sha256"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--all", action="store_true", help="Audit A(0) using the existing local live reports")
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            return 0
        require(args.all, "Specify --self-test or --all")
        root = args.repo_root.resolve()
        result = run(root)
        dest = root / OUT_REL
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(result["status"])
        print("A(0) y1 and y1Conj: 8 nonzero coefficients per epsilon orientation, including H.c.")
        print("Spinor phase eta / RGBeta epsilon orientation: NOT NATIVE-CERTIFIED")
        print("Physical cross-package coupling dictionary: NOT CERTIFIED")
        print("Report: " + str(dest))
        return 0
    except (BridgeError, OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        print("A(0) physical Yukawa action bridge BLOCKED: " + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
