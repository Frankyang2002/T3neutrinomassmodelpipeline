"""Calculate the one-generation direct Weinberg beta in general D(-2) EFT1.

Required first-stage producer:
    python -X utf8 -m scripts.RunGeneralDMinus2Stage1

This first-stage script must include the real-basis tensor-seed exporters in
Lagrangian/RunGeneralDMinus2Stage1.wl (both stage-1 matching and exports run
in the same kernel; no restricted theory is ever substituted).

Run from repository root:
    python -X utf8 -m scripts.RunGeneralDMinus2RealEFT1RGE

This calculates O(1/16pi^2) LLSS -> LLHH mixing, NOT a full-flavor
insertion or S1/R/I threshold. The neutral mass constraints are symbolic.
"""
from __future__ import annotations

import argparse
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sys
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DIR = ROOT / "output/general_sequential_pilot/D_m2_stage1"
EXPECTED_QUARTICS = {
    "lambdaT3ExtraQ000004I1", "lambdaT3ExtraQ000013I1",
    "lambdaT3ExtraQ001102I1", "lambdaT3ExtraQ020101I1",
    "lambdaT3ExtraQ110002I1",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing exact physical-basis input: {path}")
    doc = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(doc, dict):
        raise ValueError(f"Invalid JSON object: {path}")
    return doc


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str)+"\n", encoding="utf-8")


def validate_seeds(output: Path, wilson: dict, quartic: dict) -> dict:
    from scripts.RunGeneralDMinus2Stage1 import audit_stage1
    audit = audit_stage1(output)
    require(audit["status"] == "GeneralDMinus2RealMassBasisFermionThresholdVerified",
            "General D(-2) physical stage-1 audit failed")
    require(wilson.get("Status") == quartic.get("Status") == "Success",
            "A tensor seed is not successful")
    require(wilson.get("Model") == quartic.get("Model") == "D_-2",
            "A tensor seed belongs to another model")
    require(wilson.get("PhysicalScalarNames") == ["NewScalar1","T3S2R","T3S2I"],
            "Tree seed is not physical D(-2) field basis")
    require(quartic.get("PhysicalScalarNames") == ["H","NewScalar1","T3S2R","T3S2I"],
            "Quartic seed is not physical D(-2) field basis")
    terms = wilson.get("TreeWilsonTerms", [])
    require(len(terms) == wilson.get("TreeWilsonTermCount") == 8,
            "Expected eight tree-level LLSS terms from verified stage-1")
    qterms = quartic.get("ScalarQuarticTerms", [])
    require(bool(qterms) and len(qterms) == quartic.get("ScalarQuarticTermCount"),
            "Invalid real-basis quartic term inventory")
    qtext = "\n".join(x.get("TermInputForm", "") for x in qterms)
    missing = sorted(q for q in EXPECTED_QUARTICS if q not in qtext)
    require(not missing, "General D(-2) quartic missing from real basis: " + ",".join(missing))
    for coll in (terms, qterms):
        for record in coll:
            text = record.get("TermInputForm", "")
            require("Field[NewScalar2" not in text,
                    "An original complex S2 field entered the real-basis tensor seed")
            require("lambdaT3ExtraM000002I1" not in text,
                    "Off-diagonal holomorphic bilinear was not absorbed into masses")
    require(all(x.get("Chirality") in ("PL","PR") for x in terms),
            "Tree-level Wilson chirality not classified")
    for seed, term_key in ((wilson, "TreeWilsonTerms"),
                           (quartic, "ScalarQuarticTerms")):
        declared = {x.get("Name") for x in seed.get("CGRegistry", [])}
        references = set()
        from RGE.running.intermediate.ScalarOnlyTensorAdapters import (
            _find_cg_calls,
        )
        for rec in seed[term_key]:
            references.update(call.name for call in _find_cg_calls(rec["TermInputForm"]))
        require(references <= declared,
                f"Same-kernel SU(2) CG registry incomplete: {sorted(references-declared)}")
    return {
        "real_scalar_field_content_verified": True,
        "fermion_threshold_verified": True,
        "tree_wilson_terms": len(terms),
        "left_chiral_terms": sum(x["Chirality"] == "PL" for x in terms),
        "right_chiral_terms": sum(x["Chirality"] == "PR" for x in terms),
        "scalar_quartic_terms": len(qterms),
        "all_five_additional_quartics_present": True,
        "same_kernel_cg_registries_complete": True,
    }


def compute_rge(output: Path) -> dict[str, Any]:
    work = output / "data/general_d_m2_eft1_rge"
    summary = work / "pilot_summary.json"
    result: dict[str, Any] = {
        "status": "Started", "model": "D_-2",
        "scope": "general real-basis D(-2) EFT1 one-generation one-loop LLSS-to-Weinberg mixing",
        "output_directory": str(output),
        "physical_mass_relations_imposed": False,
        "full_flavor_running_performed": False,
        "scalar_threshold_performed": False,
        "renormalisable_rgbeta_mapping_certified": False,
    }
    write_json(summary, result)
    try:
        wpath = output / "data/eft1_after_F_realbasis_wilson_tensor_seed.json"
        qpath = output / "data/eft1_after_F_realbasis_quartic_tensor_seed.json"
        wilson, quartic = read_json(wpath), read_json(qpath)
        result["seed_audit"] = validate_seeds(output, wilson, quartic)
        print("[D(-2)] General mass-basis EFT1 seeds and CG registries verified.",flush=True)

        from RGE.general.ScalarBasis import ScalarBasis
        from RGE.general.FermionBasis import build_gauge_sectors
        from RGE.general.WilsonTensorRGE import (
            WilsonRGEInputs, calculate_wilson_tensor_rge,
        )
        from RGE.general.AnomalousDimensions import with_collinear_anomalous_dimensions
        from RGE.matching.WeinbergTensorAdapter import build_sm_eft, build_sm_yukawa
        from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import (
            weinberg_subspace_diagnostics,
        )
        from RGE.running.intermediate.GeneralDMinus2RealTensor import (
            build_real_wilson, build_real_quartic,
        )
        # Gauge structure is equivalent to a neutral complex singlet *as a
        # real vector space*, even though R,I are physically distinct real
        # mass eigenfields. Hypercharge=0 => their gauge generators are zero.
        model = ScalarBasis.t3(3, sp.Integer(-1), 1, sp.Integer(0))
        sm_model, fb = build_sm_eft()
        l1, l2 = fb.global_index("L",1),fb.global_index("L",2)
        wc, used_y = build_real_wilson(wilson,(l1,l2))
        require(bool(wc.entries),"Real-basis PL Wilson tensor is identically zero")
        require(not any(wc.symmetry_residuals().values()),
                "Real-basis LLSS tensor violates index symmetries")
        require("y2Conj" in used_y and any(
            sp.Symbol("y2Conj") in v.free_symbols
            for v in wc.entries.values()),
            "Independent conjugate-scalar Yukawa absent from the nonzero PL Wilson tensor")
        qt = build_real_quartic(quartic)
        require(bool(qt.entries),"General scalar quartic tensor is zero")
        free = set().union(*(v.free_symbols for v in qt.entries.values()))
        seen_extra = {s.name for s in free if s.name in EXPECTED_QUARTICS}
        require(seen_extra == EXPECTED_QUARTICS,
                "General quartic lost in real tensor conversion: " +
                str(sorted(EXPECTED_QUARTICS-seen_extra)))
        print(f"[D(-2)] C: {len(wc.entries)} sparse components; "
              f"lambda: {len(qt.entries)} symmetric components.",flush=True)

        inp = WilsonRGEInputs(
            fermion_dimension=fb.dimension,
            yukawa=build_sm_yukawa(model,fb),
            quartic=qt,
            gauge_sectors=build_gauge_sectors(model,fb),
        )
        inp = with_collinear_anomalous_dimensions(model=model,inputs=inp)
        h = list(sm_model.block("H").indices)
        candidates = [
            (i,j,a,b)
            for i,j in combinations_with_replacement((l1,l2),2)
            for a,b in combinations_with_replacement(h,2)
        ]
        betas = {}
        for k,component in enumerate(candidates,1):
            entry = calculate_wilson_tensor_rge(
                model=model, inputs=inp, output_component=component,
                coefficient=wc, simplify_each=False,
            )
            if entry["total"] != 0:
                betas[component] = entry
            if k % 10 == 0:
                print(f"[D(-2)] Evaluated {k}/{len(candidates)} LLHH components.",flush=True)
        # The shared Weinberg template is evaluated on the SAME real H 1..4
        # coordinates; hence no extra normalisation is introduced here.
        class Context:
            pass
        context = Context()
        context.model = model
        context.fermion_basis = fb
        validation = weinberg_subspace_diagnostics(context,betas)
        rge_path = work/"weinberg_mixing_beta.json"
        payload = {
            "status": "Success" if validation["matches_weinberg_subspace"] else "WeinbergProjectionFailed",
            "operator": "one-generation LLSS -> LLHH",
            "loop": "16*pi^2*dC/dln(mu) = beta[C_tree]",
            "scalar_order": ["H_R1","H_I1","H_R2","H_I2",
                             "S1_R1","S1_I1","S1_R2","S1_I2","S1_R3","S1_I3",
                             "T3S2R","T3S2I"],
            "tree_wilson_component_count": len(wc.entries),
            "quartic_symmetric_component_count": len(qt.entries),
            "quartic_general_couplings_in_tensor": sorted(seen_extra),
            "yukawa_symbols_in_PL_boundary": sorted(used_y),
            "checked_LLHH_components": len(candidates),
            "nonzero_LLHH_beta_components": len(betas),
            "weinberg_subspace_validation": validation,
            "components": {
                ",".join(map(str,key)): {
                    name: sp.sstr(sp.factor(v))
                    for name,v in parts.items()
                }
                for key,parts in sorted(betas.items())
            },
            "full_flavor_transport_completed": False,
            "scalar_threshold_completed": False,
            "mass_eigenvalue_relations_imposed": False,
        }
        write_json(rge_path,payload)
        result["beta"] = {
            "path": str(rge_path),
            "LLHH_component_count": len(candidates),
            "nonzero_component_count": len(betas),
            "weinberg_subspace_validation": validation,
            "general_quartics_included": True,
            "extra_Yukawa_included": True,
        }
        require(validation["matches_weinberg_subspace"],
                "One-loop LLHH beta is not proportional to the Weinberg template")
        require(len(betas)>0,"No one-loop LLSS-to-LLHH mixing generated")
        result["status"] = "GeneralDMinus2RealBasisWeinbergBetaVerifiedOneGeneration"
        print("[D(-2)] One-generation general EFT1 Weinberg beta validated.",flush=True)
    except Exception as exc:
        result["status"]="Failed"
        result["error"]=f"{type(exc).__name__}: {exc}"
        print(f"ERROR: {result['error']}",file=sys.stderr)
    finally:
        write_json(summary,result)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=DEFAULT_DIR)
    args = parser.parse_args(argv)
    result = compute_rge(args.output.resolve())
    print(f"{result['status']} -> {args.output.resolve()/'data/general_d_m2_eft1_rge/pilot_summary.json'}")
    return 0 if result["status"] == "GeneralDMinus2RealBasisWeinbergBetaVerifiedOneGeneration" else 1


if __name__ == "__main__":
    raise SystemExit(main())
