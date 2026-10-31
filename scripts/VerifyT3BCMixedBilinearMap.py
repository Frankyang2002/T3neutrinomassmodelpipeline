"""Workstream B: fail-closed B(-1)/C(-1) mixed-bilinear vertex audit.

Requires LIVE outputs from scripts/ProbeT3BCMixedBilinear{Matchete,RGBeta}.wl.
No production file modifications. A candidate common-component map is computed,
but no global field-basis or physical matching certificate is issued.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
MODELS = ("B_-1", "C_-1")
EPS = sp.Matrix([[0, 1], [-1, 0]])
EPS_NUM = np.array(EPS.tolist(), dtype=complex)
MAT_FORMAT = "t3-bc-live-mixed-matchete-bilinear-v1"
RGB_FORMAT = "t3-bc-live-mixed-rgbeta-bilinear-v1"
COUPLING_M = "lambdaT3ExtraM000101I1"
COUPLING_R = "muT3GeneralS1S2"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError("B/C live mixed bilinear map: " + message)


def parse_component(text: str) -> sp.Expr:
    require(isinstance(text, str) and len(text) < 100,
            "Invalid native CG component")
    require(re.fullmatch(r"[0-9\s+*/().,\[\]^\-A-Za-z]+", text) is not None,
            "Unsupported native CG expression")
    tokens = set(re.findall(r"[A-Za-z][A-Za-z0-9]*", text))
    require(tokens <= {"Sqrt", "sqrt", "Rational", "Power", "I", "Pi", "pi"},
            "Unrecognized native CG symbols: " + str(tokens))
    converted = text.replace("[", "(").replace("]", ")").replace("^", "**")
    value = sp.sympify(converted, locals={"Sqrt": sp.sqrt, "sqrt": sp.sqrt, "Rational": sp.Rational,
                  "Power": sp.Pow, "I": sp.I, "Pi": sp.pi})
    require(not value.free_symbols and value.is_finite is True,
            "Native CG component is not a finite constant")
    return sp.simplify(value)


def read_json(path: Path) -> dict:
    require(path.is_file(), "Missing live result: " + str(path))
    document = json.loads(path.read_text(encoding="utf-8-sig"))
    require(isinstance(document, dict), "Malformed live report: " + str(path))
    return document


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_source(root: Path) -> dict:
    paths = (
        "Lagrangian/interactions/GeneralScalarOperators.wl",
        "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl",
        "RGE/running/rgbeta/T3RGBetaModel.wl",
    )
    contents = {}
    for p in paths:
        q = root / p
        require(q.is_file(), "Missing source: " + p)
        contents[p] = q.read_text(encoding="utf-8-sig")
    a, b, c = (contents[p] for p in paths)
    require('"Expression" -> PlusHc[contribution]' in a and
            'contribution = coupling[] product' in a and
            'DefineCoupling[coupling, SelfConjugate -> False]' in a,
            "Matchete source bilinear/HC convention changed")
    require(re.search(r'AddScalarMass\[muT3GeneralS1S2,\s*\{Bar\s*@\s*S1,\s*Bar\s*@\s*S2\}', b),
            "RGBeta mixed scalar mass fields changed")
    require('GroupInvariant -> (eps[SU2L @ fund,#1,#2] &)' in b and
            'SelfConjugate -> False' in b,
            "RGBeta mixed scalar mass invariant/H.c. changed")
    require('T3RGBetaBuildEFT1[' in c and 'T3RGBetaAddScalarMasses[dS1, dS2]' in c,
            "RGBeta EFT1 scalar mass builder changed")
    return {p: sha(root / p) for p in paths}


def check_matchete(payload: dict, model: str) -> sp.Matrix:
    require(payload.get("format") == MAT_FORMAT and
            payload.get("physical_mapping_certified") is False,
            "Wrong Matchete report version or certification state")
    rows = payload.get("models", {})
    require(set(rows) == set(MODELS), "Incomplete Matchete model collection")
    item = rows[model]
    require(item.get("model") == model and
            item.get("status") == "LiveMatcheteMixedBilinearCGCaptured" and
            item.get("counts") == [0,0,0,1,0,1] and
            item.get("matchete_coupling") == COUPLING_M and
            item.get("matchete_cg") == "T3ExtraCGM000101I1" and
            item.get("tensor_accessor") == "Matchete`GetCGTensor" and
            item.get("cg_present_in_actual_plus_hc_expression") is True and
            item.get("physical_mapping_certified") is False,
            model + ": Matchete native source contract incomplete")
    raw = item.get("cg_components_inputform")
    require(isinstance(raw, list) and len(raw) == 2 and
            all(isinstance(x,list) and len(x)==2 for x in raw),
            model + ": Matchete CG is not 2x2")
    return sp.Matrix([[parse_component(x) for x in row] for row in raw])


def check_rgbeta(payload: dict, model: str) -> dict:
    require(payload.get("format") == RGB_FORMAT and
            payload.get("physical_mapping_certified") is False,
            "Wrong RGBeta report version or certification state")
    rows = payload.get("models", {})
    require(set(rows) == set(MODELS), "Incomplete RGBeta model collection")
    item = rows[model]
    require(item.get("model") == model and
            item.get("status") == "LiveRGBetaMixedBilinearRegisteredAndLamDetected" and
            item.get("rgbeta_mass_symbol") == COUPLING_R and
            item.get("rgbeta_fields") == ["Bar[S1]", "Bar[S2]"] and
            item.get("rgbeta_self_conjugate") is False and
            item.get("rgbeta_unique_permutations") == 12 and
            item.get("native_lam_symbol") == "RGBeta`PackageScope`Lam" and
            item.get("native_massive_lam_evaluated") is True and
            item.get("native_lam_mass_sparse_array_rules_count",0) > 0 and
            item.get("native_rgbeta_epsilon_component_phase_verified") is False and
            item.get("physical_mapping_certified") is False,
            model + ": RGBeta native source contract incomplete")
    coefficient = item.get("rgbeta_mass_coupling_function_inputform", "")
    require(''.join(coefficient.split()) == "12*Matrix[muT3GeneralS1S2][]",
            model + ": internal real mass-tensor normalization not 12")
    invariant = ''.join(item.get("rgbeta_invariant_inputform", "").split())
    require(invariant in ("eps[SU2L[fund],i,j]", "eps[SU2L@fund,i,j]"),
            model + ": native epsilon invariant changed")
    require(item.get("rgbeta_potential_convention") ==
            "+muR * eps * Bar[S1] Bar[S2] + h.c.; - sign in Lagrangian",
            model + ": RGBeta source sign convention changed")
    return {"internal_real_tensor_factor": 12,
            "distinct_massive_Lam_arrangements": 12,
            "massive_Lam_witness_count": item["native_lam_mass_sparse_array_rules_count"],
            "physical_complex_potential_prefactor_under_package_convention": 1,
            "native_epsilon_component_phase_verified": False}


def bilinear(z1: np.ndarray, z2: np.ndarray, tensor: np.ndarray) -> complex:
    return complex(np.einsum("i,ij,j->", np.conjugate(z1), tensor,
                                   np.conjugate(z2)))


def hessian_numeric(func, step: float = 1.e-3) -> np.ndarray:
    result = np.empty((8,8), dtype=float)
    e = np.eye(8)*step
    for i in range(8):
        for j in range(8):
            result[i,j] = (func(e[i]+e[j])-func(e[i]-e[j])-
                           func(-e[i]+e[j])+func(-e[i]-e[j])) / (4*step*step)
    return result


def check_equivalence(mat_tensor: sp.Matrix, k: sp.Expr) -> dict:
    mat = np.array(mat_tensor.evalf().tolist(), dtype=complex)
    factor = complex(k.evalf())
    rng = np.random.default_rng(20261031)
    max_poly = max_hessian = max_mass_rotation = 0.0
    negative_conj = negative_factor = 0.0
    for _ in range(24):
        mu_m = complex(*rng.uniform(-1.0,1.0,size=2))
        mu_r = -factor * mu_m  # SAME holomorphic S1dag S2dag orientation
        m1 = float(rng.uniform(2,3))
        m2 = float(rng.uniform(4,5))
        v1 = rng.normal(size=2)+1j*rng.normal(size=2)
        v2 = rng.normal(size=2)+1j*rng.normal(size=2)
        poly = bilinear(v1,v2,mat)
        rgb_poly = bilinear(v1,v2,EPS_NUM)
        lmat = 2*np.real(mu_m*poly)
        vrgb = 2*np.real(mu_r*rgb_poly)
        scale = max(1.0,abs(lmat),abs(vrgb))
        max_poly = max(max_poly,abs(lmat+vrgb)/scale)
        wrong_conj = 2*np.real(-factor*mu_m.conjugate()*rgb_poly)
        wrong_factor = 2*np.real(-2*factor*mu_m*rgb_poly)
        negative_conj = max(negative_conj,abs(lmat+wrong_conj)/scale)
        negative_factor = max(negative_factor,abs(lmat+wrong_factor)/scale)
        def potential(q, mode):
            z=(q[:4]+1j*q[4:])/np.sqrt(2.0)
            s1,s2=z[:2],z[2:]
            diag = m1*np.vdot(s1,s1).real + m2*np.vdot(s2,s2).real
            if mode == "mat":
                return float(diag-2*np.real(mu_m*bilinear(s1,s2,mat)))
            return float(diag+2*np.real(mu_r*bilinear(s1,s2,EPS_NUM)))
        hm = hessian_numeric(lambda q:potential(q,"mat"))
        hr = hessian_numeric(lambda q:potential(q,"rgb"))
        max_hessian=max(max_hessian,float(np.max(np.abs(hm-hr))))
        eig=np.linalg.eigvalsh(hr)
        avg=(m1+m2)/2
        split=np.sqrt(((m1-m2)/2)**2+abs(mu_r)**2)
        expected=np.sort(np.array([avg-split]*4+[avg+split]*4))
        max_mass_rotation=max(max_mass_rotation,
                              float(np.max(abs(eig-expected))))
    require(max_poly < 1e-12, "Matchete + RGBeta potential sign/prefactor disagrees")
    require(max_hessian < 1e-8, "Independent real mass Hessians disagree")
    require(max_mass_rotation < 1e-8, "Mixed scalar spectrum/basis discrepancy")
    require(negative_conj > 1.e-3 and negative_factor > 1.e-3,
            "Conjugation/factor-two negative controls are not discriminating")
    return {"max_action_relative_residual":max_poly,
            "max_real_hessian_difference":max_hessian,
            "max_mass_eigenvalue_residual":max_mass_rotation,
            "wrong_complex_conjugation_negative_control":negative_conj,
            "wrong_factor_two_negative_control":negative_factor,
            "numerical_trials":24}


def audit(mat_payload:dict, rgb_payload:dict, root:Path=ROOT) -> dict:
    provenance=check_source(root)
    for report,name in ((mat_payload,"Matchete"),(rgb_payload,"RGBeta")):
        actual = report.get("source_sha256")
        expected_paths = (("Lagrangian/interactions/GeneralScalarOperators.wl",)
                          if name=="Matchete" else
                          ("RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl",
                           "RGE/running/rgbeta/T3RGBetaModel.wl"))
        require(isinstance(actual,dict) and set(actual)==set(expected_paths),
                name+": missing native source digests")
        for filename in expected_paths:
            require(actual[filename]==provenance[filename],
                    name+": source changed since native probe: "+filename)
    rows=[]
    for model in MODELS:
        tensor=check_matchete(mat_payload,model)
        rgbmeta=check_rgbeta(rgb_payload,model)
        ratio=sp.simplify(tensor[0,1])
        require(ratio != 0 and ratio.is_zero is False,
                model+": Matchete CG epsilon factor vanished")
        require(all(sp.simplify(x)==0 for x in (tensor-ratio*EPS)),
                model+": live CG tensor not proportional to declared epsilon")
        tests=check_equivalence(tensor,ratio)
        rows.append({"model":model,
            "native_matchete_registered_cg":[[str(tensor[i,j]) for j in range(2)] for i in range(2)],
            "matchete_over_declared_epsilon_exact":str(ratio),
            "conditional_mu_rgbeta_over_mu_matchete_exact":str(-ratio),
            "conditional_mu_matchete_over_mu_rgbeta_exact":str(-1/ratio),
            "requires_conjugation_under_same_holomorphic_field_word":False,
            "holomorphic_field_order":["S1dag","S2dag"],
            "basis_rotation":"(S1, epsilon * conjugate(S2)) = U (X,Y)",
            "rgbeta_mass_meta":rgbmeta,
            "tests":tests,
            "physical_mapping_certified":False})
    require(sp.simplify(parse_component(rows[0]["matchete_over_declared_epsilon_exact"])-
                            parse_component(rows[1]["matchete_over_declared_epsilon_exact"]))==0,
            "B/C live Matchete mixed CG phase is unexpectedly different")
    return {"format":"t3-bc-independent-mixed-bilinear-convention-audit-v1",
        "status":"NativeSourceAndSameComponentActionMapCheckedConditionalNotPhysicalCertificate",
        "models":rows,
        "source_sha256":provenance,
        "assumptions":[
           "Canonical RGBeta component epsilon_12=+1 aligned with Matchete labeled components",
           "RGBeta -V sign and massive Lam/4! internal tensor reconstruction",
           "Matchete extra operator enters +PlusHc with no factor 1/2"],
        "limitations":[
           "No independently live-verified RGBeta epsilon component-phase map",
           "No full Matchete-RGBeta kinetic/generator intertwiner certificate",
           "No native full mass-basis Matchete matching for these running values",
           "No Workstream A final operator-basis certificate supplied"],
        "live_matchete_cg_tensor_checked":True,
        "live_rgbeta_massive_lam_symbol_detected":True,
        "cross_tool_same_component_action_identity_checked":True,
        "physical_mapping_certified":False,
        "production_modified":False}


def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matchete",type=Path,required=True)
    p.add_argument("--rgbeta",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--root",type=Path,default=ROOT)
    args=p.parse_args(argv)
    require(args.output!=args.matchete and args.output!=args.rgbeta,
            "Output must not replace native evidence")
    require(not args.output.exists(), "Refusing to overwrite existing audit report")
    result=audit(read_json(args.matchete),read_json(args.rgbeta),args.root)
    result["live_evidence_sha256"]={str(args.matchete):sha(args.matchete),
                                    str(args.rgbeta):sha(args.rgbeta)}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    for x in result["models"]:
        print(x["model"],"conditional mu_RGB/mu_M =",
              x["conditional_mu_rgbeta_over_mu_matchete_exact"],
              "max action residual =",x["tests"]["max_action_relative_residual"])
    print(result["status"])
    return 0

if __name__=="__main__":
    raise SystemExit(main())
