"""Run and audit general D(-2) F matching in the diagonal real-S2 basis.

The symmetry-allowed holomorphic singlet bilinear must be absorbed into the
physical scalar mass basis BEFORE Matchete matches F.  Neither the original
complex-S2 EFT1 tensor adapter nor the scalar threshold is invoked here.

From the repository root:
    python -X utf8 -m scripts.RunGeneralDMinus2Stage1
    python -X utf8 -m scripts.RunGeneralDMinus2Stage1 --audit-only

No restricted-production fallback and no GitHub modifications.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "output/general_sequential_pilot/D_m2_stage1"
WOLFRAM_SCRIPT = ROOT / "Lagrangian/RunGeneralDMinus2Stage1.wl"
MODEL = "D_-2"
EXTRA_Q_RE = re.compile(r"lambdaT3ExtraQ\d+I\d+")
EXTRA_M = "lambdaT3ExtraM000002I1"


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Required output missing: {path}")
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(data, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return data


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def extra_quartic_coupling(operator_name: str) -> str:
    """T3ExtraQ_0_0_0_0_0_4_1 -> lambdaT3ExtraQ000004I1."""
    parts = operator_name.split("_")
    if len(parts) != 8 or parts[0] != "T3ExtraQ":
        raise ValueError(f"Unexpected general quartic operator: {operator_name!r}")
    if not all(part.isdecimal() for part in parts[1:]):
        raise ValueError(f"Invalid general quartic operator: {operator_name!r}")
    return "lambdaT3ExtraQ" + "".join(parts[1:-1]) + "I" + parts[-1]


def audit_stage1(output: Path) -> dict[str, Any]:
    """Validate actual real-basis stage-1 outputs; fail on stale prior files."""
    output = Path(output)
    inventory = read_json(output / "general_uv_inventory.json")
    basis = read_json(output / "scalar_mass_basis.json")
    stage = read_json(output / "first_threshold_status.json")
    wilson = read_json(output / "data/eft1_after_F_realbasis_wilson_inventory.json")
    physical_path = output / "data/physical_uv_realbasis_inputform.txt"
    tree_path = output / "data/eft1_after_F_realbasis_tree_inputform.txt"
    loop_path = output / "data/eft1_after_F_realbasis_loop_inputform.txt"
    require(physical_path.is_file() and tree_path.is_file() and loop_path.is_file(),
            "Physical-basis UV/EFT1 expression export missing")
    physical = physical_path.read_text(encoding="utf-8-sig")
    tree = tree_path.read_text(encoding="utf-8-sig")

    require(inventory.get("Model") == MODEL and inventory.get("BuildStatus") == "Success",
            "Matched theory is not the validated general D(-2) model")
    require(inventory.get("HeavyFieldsAtFirstThreshold") == ["F"],
            "The first threshold must integrate out F only")
    accepted = inventory.get("AllowedInteractions", [])
    require(isinstance(accepted, list) and len(accepted) == len(set(accepted)),
            "Invalid UV interaction inventory")
    q = inventory.get("AdditionalQuartics", [])
    m = inventory.get("AdditionalQuadratics", [])
    y = inventory.get("AdditionalYukawas", [])
    require(isinstance(q, list) and len(q) == 5
            and inventory.get("AdditionalQuarticsCount") == 5
            and all(isinstance(name, str) and name in accepted
                    and name.startswith("T3ExtraQ_") for name in q),
            "Full general quartic inventory is incomplete")
    require(isinstance(m, list) and len(m) == 1
            and inventory.get("AdditionalQuadraticsCount") == 1
            and m[0] in accepted and m[0] == "T3ExtraM_0_0_0_0_0_2_1",
            "D(-2) holomorphic mass bilinear missing")
    require(isinstance(y, list) and len(y) == 1
            and inventory.get("AdditionalYukawasCount") == 1
            and y[0] in accepted and y[0].startswith("Yukawa2ConjugateScalar"),
            "D(-2) additional Yukawa missing")

    require(basis.get("Model") == MODEL
            and basis.get("RealFields") == ["T3S2R", "T3S2I"]
            and basis.get("OriginalComplexField") == "NewScalar2"
            and basis.get("RemovedConjugateBilinearTermCount") == 2
            and basis.get("HolomorphicCoupling") == EXTRA_M,
            "Neutral complex scalar was not fully rotated to two real fields")
    require(basis.get("PhysicalMassRelationsImposedInMatchete") is False
            and basis.get("ScalarThresholdPerformed") is False,
            "Physical mass-constraint or scalar-threshold claim is unsupported")
    require(basis.get("RealMassSquared") == "MS2^2-2*Abs[kappa]"
            and basis.get("ImagMassSquared") == "MS2^2+2*Abs[kappa]",
            "Mass eigenvalue convention does not match the UV bilinear sign")
    required_q = {extra_quartic_coupling(name) for name in q}
    require(set(basis.get("GeneralQuarticCouplingNames", [])) == required_q,
            "General quartic names disagree with original UV inventory")
    for name in required_q:
        require(name in physical and name in tree,
                f"General quartic lost during real-basis rotation or first threshold: {name}")
    require("y2Conj" in physical and "y2Conj" in tree,
            "Independent conjugate-scalar Yukawa lost in matching")
    require(EXTRA_M not in physical and EXTRA_M not in tree,
            "Off-diagonal holomorphic bilinear was left in physical matching input")
    require("Field[NewScalar2" not in physical and "Field[NewScalar2" not in tree,
            "Original complex scalar incorrectly retained after rotation")

    require(stage.get("Status") == "Success" and stage.get("Model") == MODEL
            and stage.get("ThresholdPlan") == [["F"]]
            and stage.get("IntegratedFields") == ["F"]
            and stage.get("ActiveFields") == ["S1", "S2R", "S2I"],
            "Incorrect fermion-first matching stage or active field basis")
    for key in ("FermionAbsentFromEFT1Tree", "Scalar1PresentInEFT1Tree",
                "RealScalarRPresentInEFT1Tree", "RealScalarIPresentInEFT1Tree",
                "OriginalComplexS2AbsentFromEFT1Tree", "GeneralExtraQuarticsRetained",
                "ExtraYukawaRetainedInEFT1", "ScalarMassBasisDiagonalized"):
        require(stage.get(key) is True, f"EFT1 real-basis audit failed: {key}")
    for key in ("ScalarThresholdPerformed", "IntermediateWilsonRunningPerformed",
                "PhysicalMassRelationsImposed", "StandardComplexScalarTensorAdapterCompatible"):
        require(stage.get(key) is False, f"Unsupported completion claim: {key}")
    require("Field[NewFermion" not in tree
            and "Field[NewScalar1" in tree
            and "Field[T3S2R" in tree and "Field[T3S2I" in tree,
            "Incorrect EFT1 dynamical field content")

    terms = wilson.get("TreeWilsonTermsInputForm", [])
    require(wilson.get("Status") == "Success"
            and wilson.get("StandardComplexScalarTensorAdapterCompatible") is False
            and wilson.get("PhysicalScalarNames") ==
                ["NewScalar1", "T3S2R", "T3S2I"]
            and isinstance(terms, list) and bool(terms)
            and len(terms) == wilson.get("TreeWilsonTermCount")
            and len(terms) == stage.get("TreeWilsonTermCount"),
            "Real-basis LLSS seed missing or falsely marked as legacy compatible")
    require(all("GammaCC" in term for term in terms),
            "Non-LLSS term misidentified as real-basis Wilson seed")
    for filename in (
        "eft1_after_F_realbasis_tree.wxf",
        "eft1_after_F_realbasis_loop.wxf",
        "eft1_after_F_realbasis_transition_tree.wxf",
        "eft1_after_F_realbasis_transition_full.wxf",
    ):
        require((output / "data" / filename).is_file(),
                f"Native Wolfram continuation missing: {filename}")

    return {
        "status": "GeneralDMinus2RealMassBasisFermionThresholdVerified",
        "model": MODEL,
        "general_uv": {
            "interaction_count": len(accepted),
            "extra_quartics": len(q),
            "neutral_scalar_bilinears": len(m),
            "extra_yukawas": len(y),
            "symmetry_completeness_certified": bool(
                inventory.get("SymmetryCompletenessCertified", False)),
            "fermion_flavor_mass_generality_certified": bool(
                inventory.get("GeneralFermionFlavorMassCertified", False)),
        },
        "first_matching": {
            "integrated_out": ["F"],
            "remaining_physical_BSM_scalars": ["NewScalar1", "T3S2R", "T3S2I"],
            "tree_wilson_terms": len(terms),
            "additional_quartics_retained": len(required_q),
            "additional_yukawa_retained": True,
        },
        "mass_basis": {
            "status": basis["Status"],
            "field_decomposition": basis["PhysicalFieldDefinition"],
            "real_mass_squared": basis["RealMassSquared"],
            "imag_mass_squared": basis["ImagMassSquared"],
            "phase_constraint": basis["PhaseConstraint"],
            "physical_mass_relations_imposed": False,
        },
        "eft1_wilson_running_performed": False,
        "second_scalar_threshold_performed": False,
        "standard_complex_scalar_tensor_adapter_compatible": False,
        "full_general_sequential_C5_certified": False,
    }


def run_wolfram(output: Path, executable: str) -> int:
    if not WOLFRAM_SCRIPT.is_file():
        raise FileNotFoundError(f"Wolfram script missing: {WOLFRAM_SCRIPT}")
    output.mkdir(parents=True, exist_ok=True)
    command = [executable, "-file", str(WOLFRAM_SCRIPT), str(output)]
    env = dict(os.environ)
    env["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] = "1"
    env["T3_EXPERIMENTAL_GENERAL_FERMIONS"] = "1"
    log_path = output / "wolfram_stage1.log"
    print("Matching D(-2) F in the physical real-S2 mass basis...", flush=True)
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command, cwd=ROOT, env=env, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            encoding="utf-8", errors="replace",
        )
        assert process.stdout is not None
        for line in process.stdout:
            log.write(line)
            log.flush()
            print(line, end="", flush=True)
        code = process.wait()
    print(f"Wolfram exit code: {code}", flush=True)
    return code


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--wolframscript", default="wolframscript")
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args(argv)
    output = args.output.resolve()
    summary = output / "data/general_d_m2_stage1_summary.json"
    try:
        if not args.audit_only:
            code = run_wolfram(output, args.wolframscript)
            require(code == 0,
                    f"Wolfram stage1 returned code {code}; inspect {output/'wolfram_stage1.log'}")
        result = audit_stage1(output)
        summary.parent.mkdir(parents=True, exist_ok=True)
        summary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(f"{result['status']} -> {summary}")
        return 0
    except (OSError, RuntimeError, ValueError, KeyError, TypeError,
            json.JSONDecodeError) as exc:
        failed = {"status": "Failed", "model": MODEL,
                  "error": f"{type(exc).__name__}: {exc}"}
        summary.parent.mkdir(parents=True, exist_ok=True)
        summary.write_text(json.dumps(failed, indent=2) + "\n", encoding="utf-8")
        print(f"ERROR: {failed['error']}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
