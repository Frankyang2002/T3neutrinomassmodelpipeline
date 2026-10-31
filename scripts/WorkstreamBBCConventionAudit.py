"""Read-only B/C(-1) scalar RGE convention audit. No physical-map certification.

Run from repository root:
  python -X utf8 -m scripts.WorkstreamBBCConventionAudit \
      --output output/audits/workstream_b_bc_conventions.json
  python -X utf8 -m scripts.WorkstreamBBCConventionAudit --self-test

The audit deliberately uses independently derived algebraic identities rather
than treating a previous report's `status` as proof. Source observations and
mathematical tests are kept distinct from physical Matchete/RGBeta equivalence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
MODELS = ("B_-1", "C_-1")
MASS = "muT3GeneralS1S2"
M1 = "mS1Sq"
M2 = "mS2Sq"
L1 = "lambdaS1S1bar2S2bar"
L2 = "lambdaS1barS2S2bar2"
L3 = "lambdaS1bar2S2bar2"
MATS = {
    "000112": ("lambdaT3ExtraQ000112I1",),
    "000202": ("lambdaT3ExtraQ000202I1",),
    "001201": ("lambdaT3ExtraQ001201I1",),
    "020020": ("lambdaT3ExtraQ020020I1",),
    "020200": ("lambdaT3ExtraQ020200I1",),
    "110101": ("lambdaT3ExtraQ110101I1", "lambdaT3ExtraQ110101I2"),
}
RGBS = {
    "000112": (L2,), "000202": (L3,), "001201": (L1,),
    "020020": ("lambdaHHdagS2S2",),
    "020200": ("lambdaHHdagS1barS1bar",),
    "110101": ("lambdaT3GeneralQ110101I1", "lambdaHHdagS1barS2barCross"),
}
A110101 = sp.Matrix([[sp.Rational(1, 2), 1], [-sp.sqrt(3)/2, 0]])
EXPECTED_110101 = -A110101.inv()
EPS = np.array([[0, 1], [-1, 0]], dtype=complex)
POLY_SYMBOLS = set((MASS, M1, M2, L1, L2, L3, "lambdaS1", "lambdaS2",
                    "lambda12", "lambda12Adj", "g2", "gY"))
REAL_SYMBOLS = POLY_SYMBOLS - {MASS, L1, L2, L3}

def symbol(name: str):
    return sp.Symbol(name, real=True) if name in REAL_SYMBOLS else sp.Symbol(name, complex=True)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def source(root: Path, path: str) -> str:
    return (root / path).read_text(encoding="utf-8-sig")


def read_json(root: Path, path: str) -> dict:
    obj = json.loads(source(root, path))
    require(isinstance(obj, dict), f"Invalid JSON object: {path}")
    return obj


def source_digest(root: Path, path: str) -> str:
    return hashlib.sha256((root/path).read_bytes()).hexdigest()


def exact(value: str):
    # Only expressions from previously generated symbolic matrices; no eval.
    return sp.sympify(value, locals={"sqrt": sp.sqrt}, evaluate=True)


def matrix(values) -> sp.Matrix:
    require(isinstance(values, list) and bool(values), "Missing exact matrix")
    return sp.Matrix([[exact(x) for x in row] for row in values])


def same_matrix(x: sp.Matrix, y: sp.Matrix) -> bool:
    return x.shape == y.shape and all(sp.simplify(z) == 0 for z in x-y)


def audit_operator_map(root: Path) -> dict:
    registry = read_json(root, "RGE/general/generated/general_rgbeta_quartic_registry.json")
    manifest = read_json(root, "docs/general_scalar_trial_operator_manifest.json")
    handoff = read_json(root, "output/general_t3_stage22/stage29_convention_uv_handoff.json")
    stage23 = read_json(root, "output/general_t3_stage22/stage23_quartic_prefactor_candidates.json")
    require(manifest.get("model_count") == 16 and len(manifest.get("models", {})) == 16,
            "Existing Workstream A precursor manifest is incomplete")
    require(registry["totals"]["matchete_extra_quartics"] == 50 and
            registry["totals"]["quadratic_couplings_blocked"] == 8,
            "General coupling registry inventory changed")
    require(handoff.get("physical_coupling_mapping_certified") is False and
            stage23.get("signed_physical_coupling_map_certified") is False,
            "Candidate signed mapping status changed: inspect independently")
    require(handoff.get("models_checked") == 16 and handoff.get("slots_checked") == 50,
            "Incomplete Stage29 handoff")
    mat_src = source(root, "Lagrangian/interactions/GeneralScalarOperators.wl")
    rgb_src = source(root, "RGE/running/rgbeta/T3RGBetaModel.wl")
    require('"Expression" -> PlusHc[contribution]' in mat_src,
            "Matchete extra interaction source convention changed")
    for pattern in ("{Bar @ H, Bar @ H, S2, S2}",
                    "{Bar @ H, Bar @ H, Bar @ S1, Bar @ S1}",
                    "{H, Bar @ H, Bar @ S1, Bar @ S2}"):
        require(pattern in rgb_src, "Missing RGBeta registered field ordering: "+pattern)
    rows = []
    for model in MODELS:
        m = manifest["models"][model]
        r = registry["models"][model]
        h = handoff["models"][model]
        require(m["extra_quartic_complex_couplings"] == 7 and
                m["quadratic_complex_couplings_not_enabled"] == 1,
                model+": operator inventory changed")
        require([x["counts"] for x in m["quadratic_structures_not_enabled"]] ==
                [[0,0,0,1,0,1]], model+": mixed quadratic changed")
        new = {x["name"] for x in r["new_rgbeta_symbols"]}
        reused = {x["existing_name"] for x in r["preexisting_rgbeta_symbols"]}
        require(new == {"lambdaT3GeneralQ110101I1"} and len(reused) == 6 and
                reused == (set(sum((list(v) for v in RGBS.values()), [])) - new),
                model+": RGBeta registry does not cover expected 7 quartics")
        blocks = {x["block"]: x for x in h["blocks"]}
        s23 = {x["block"]: x for x in stage23["blocks"] if x.get("model")==model}
        require(set(blocks)==set(MATS)==set(s23), model+": 6 independent monomial blocks required")
        for key in MATS:
            b, p = blocks[key], s23[key]
            require(tuple(b["matchete_symbols"]) == MATS[key] and
                    tuple(b["rgbeta_symbols"]) == RGBS[key] and
                    tuple(p["matchete_couplings"]) == MATS[key] and
                    tuple(p["rgbeta_couplings"]) == RGBS[key],
                    model+"/"+key+": source symbol orientation changed")
            C = matrix(b["signed_map_exact"])
            P = matrix(p["lambda_rgbeta_from_lambda_matchete_if_assumptions_hold"])
            R = matrix(p["r_in_m"])
            require(same_matrix(C,P) and same_matrix(R*C,-sp.eye(len(MATS[key]))),
                    model+"/"+key+": exact source geometry/sign identity failed")
            require(p.get("physical_mapping_certified") is False and
                    p.get("assumed_relative_potential_sign")=="-1",
                    model+"/"+key+": conditional map contract changed")
            require(same_matrix(C,EXPECTED_110101 if key=="110101" else
                                sp.Matrix([[-1]])),
                    model+"/"+key+": expected B/C signed convention changed")
            # Numeric inverse check on arbitrary COMPLEX independent couplings.
            v = sp.Matrix([sp.Rational(2+i,7) + sp.I*sp.Rational(i+1,11)
                           for i in range(len(MATS[key]))])
            require(same_matrix(R*(C*v),-v), model+"/"+key+": complex map test failed")
        rows.append({"model": model, "quartic_slots": 7, "monomial_blocks": 6,
                     "new_rgbeta_quartics": 1, "reused_rgbeta_quartics": 6,
                     "quadratic_slots": 1,
                     "conditional_110101_rgbeta_from_matchete":
                        [[str(EXPECTED_110101[i,j]) for j in range(2)] for i in range(2)],
                     "physical_mapping_certified": False})
    return {"models": rows, "conditional_exact_matrix_checked": True,
            "independent_runtime_matchete_rgbeta_action_subtraction": False,
            "workstream_a_dependency": "Existing 16-model manifest checked; final A basis certificate not supplied"}


def eft1_poly(text: str):
    require("Tr[" not in text and "BetaTerm" not in text, "Non-polynomial EFT1 mass beta")
    # The EFT1 mass sector uses only simple Bar[name] and rational arithmetic.
    require(re.fullmatch(r"[A-Za-z0-9_+*/^()., \t\n\[\]-]+", text) is not None,
            "Unsupported beta syntax")
    value = re.sub(r"Bar\[([A-Za-z][A-Za-z0-9]*)\]", r"conjugate(\1)", text)
    require("[" not in value and "]" not in value, "Unexpected unparsed beta tensor")
    syms = {name:symbol(name) for name in POLY_SYMBOLS}
    expr = sp.sympify(value.replace("^", "**"),
                      locals={**syms,"conjugate":sp.conjugate},evaluate=True)
    require(expr.free_symbols <= set(syms.values()), "Unexpected EFT1 mass variable")
    return sp.expand(expr)


def evaluate(expr, p: dict[str,complex]) -> complex:
    subs = {symbol(k): v for k,v in p.items()}
    return complex(expr.evalf(16, subs=subs))


def mass_potential(q: np.ndarray, m1:float, m2:float, mu:complex) -> float:
    # q=[Re S1_0,Re S1_1,Re S2_0,Re S2_1,Im ...]; canonical /sqrt2.
    z=(q[:4]+1j*q[4:])/np.sqrt(2)
    s1,s2=z[:2],z[2:]
    return float(m1*np.vdot(s1,s1).real+m2*np.vdot(s2,s2).real+
                 2*np.real(mu*np.vdot(s1,EPS@s2.conjugate())))


def finite_difference_hessian(m1,m2,mu,step=0.001):
    d=8
    H=np.zeros((d,d))
    f=lambda q:mass_potential(q,m1,m2,mu)
    for i in range(d):
        for j in range(d):
            ei=np.eye(d)[i]*step
            ej=np.eye(d)[j]*step
            H[i,j]=(f(ei+ej)-f(ei-ej)-f(-ei+ej)+f(-ei-ej))/(4*step**2)
    return H


def expected_hessian(m1,m2,mu):
    A=np.diag([m1,m1,m2,m2])
    B=np.zeros((4,4),complex)
    B[:2,2:]=np.conjugate(mu)*EPS
    B[2:,:2]=B[:2,2:].T
    return np.block([[A+B.real,-B.imag],[-B.imag,A-B.real]])


def mathematical_tests():
    rng=np.random.default_rng(19031)
    errors={"hessian":0.,"mass_spectrum":0.,"rephasing":0.,
            "mass_beta_eigenvalue":0.}
    for _ in range(24):
        m1,m2=2.+float(rng.random()),4.+float(rng.random())
        mu=complex(rng.normal(scale=.2),rng.normal(scale=.2))
        H=expected_hessian(m1,m2,mu)
        D=finite_difference_hessian(m1,m2,mu)
        errors["hessian"]=max(errors["hessian"],float(np.max(np.abs(H-D))))
        evals=np.linalg.eigvalsh(H)
        split=math.sqrt(((m1-m2)/2)**2+abs(mu)**2)
        analytic=np.array([.5*(m1+m2)-split]*4+[.5*(m1+m2)+split]*4)
        errors["mass_spectrum"]=max(errors["mass_spectrum"],
                                    float(np.max(np.abs(evals-analytic))))
        a,b=rng.uniform(-math.pi,math.pi,2)
        sr=np.exp(1j*(a+b))
        tr=np.linalg.eigvalsh(expected_hessian(m1,m2,mu*sr))
        errors["rephasing"]=max(errors["rephasing"],float(np.max(np.abs(tr-evals))))
    require(errors["hessian"]<5e-9,"Independent Cartesian Hessian mismatch")
    require(errors["mass_spectrum"]<5e-12,"Mixed scalar analytic spectrum failed")
    require(errors["rephasing"]<5e-12,"Mixed scalar spectrum violates rephasing")
    return errors


def audit_mass_betas(root:Path) -> dict:
    results=[]
    # Different B/C heavy SU(2) dimensions generate different UV Yukawa traces;
    # compare EFT1 only, where F has been removed.
    saved_eft1=None
    for model in MODELS:
        tag=model.replace("-","m")
        base=f"output/general_t3_quadratic_eft1/{tag}_"
        uv=read_json(root,base+"uv_quadratic_live_rgbeta.json")
        eft=read_json(root,base+"eft1_quadratic_live_rgbeta.json")
        for stage,payload in (("UV",uv),("EFT1",eft)):
            require(payload.get("status")=="SuccessCandidateMassBeta" and
                    payload["metadata"]["QuadraticKind"]=="mixed" and
                    payload["metadata"]["QuadraticMassSymbol"]==MASS and
                    payload["metadata"]["GeneralScalarStage"]==stage and
                    payload.get("physical_mapping_certified") is False,
                    model+"/"+stage+": live quadratic beta contract changed")
            require({M1,M2,MASS} <= set(payload["beta_symbols"]),
                    model+"/"+stage+": incomplete coupled mass beta equations")
        b={k:eft["report_betas"][k] for k in (M1,M2,MASS)}
        if saved_eft1 is not None:
            require(saved_eft1==b,"B/C EFT1 mass betas differ in identical scalar sector")
        saved_eft1=b
        expr={k:eft1_poly(t) for k,t in b.items()}
        masses={symbol(MASS):0,sp.conjugate(symbol(MASS)):0}
        z=sp.simplify(expr[MASS].subs(masses))
        target=6*(symbol(L1)*symbol(M1)+symbol(L2)*symbol(M2))
        require(sp.simplify(z-target)==0,model+": zero mixed mass source changed")
        require(sp.simplify(expr[M1]-sp.conjugate(expr[M1]))==0 and
                sp.simplify(expr[M2]-sp.conjugate(expr[M2]))==0,
                model+": diagonal mass beta not manifestly real")
        rng=np.random.default_rng(3101 if model=="B_-1" else 3102)
        max_cov,max_eigen=0.,0.
        for _ in range(40):
            p={M1:3.+rng.random(),M2:5.+rng.random(),MASS:complex(.2+rng.random(),.1+rng.random()),
               L1:complex(.02*rng.normal(),.02*rng.normal()),
               L2:complex(.02*rng.normal(),.02*rng.normal()),
               L3:complex(.02*rng.normal(),.02*rng.normal()),
               "lambdaS1":float(rng.random()),"lambdaS2":float(rng.random()),
               "lambda12":float(rng.random()),"lambda12Adj":float(rng.random()),
               "g2":.6,"gY":.35}
            bm1,bm2,bmu=[evaluate(expr[k],p) for k in (M1,M2,MASS)]
            require(abs(bm1.imag)<1e-11 and abs(bm2.imag)<1e-11,
                    model+": real mass beta acquired imaginary part")
            phase=np.exp(1j*rng.uniform(-math.pi,math.pi))
            rotated=dict(p)
            for key in (L1,L2,MASS):rotated[key]*=phase
            rotated[L3]*=phase**2
            rb=[evaluate(expr[k],rotated) for k in (M1,M2,MASS)]
            max_cov=max(max_cov,abs(rb[0]-bm1),abs(rb[1]-bm2),
                        abs(rb[2]-phase*bmu))
            a=np.array([[p[M1],p[MASS]],[p[MASS].conjugate(),p[M2]]],complex)
            da=np.array([[bm1.real,bmu],[np.conjugate(bmu),bm2.real]],complex)
            eigen, U=np.linalg.eigh(a)
            projected=np.real(np.diag(U.conj().T@da@U))
            step=1e-6
            slope=(np.linalg.eigvalsh(a+step*da)-np.linalg.eigvalsh(a-step*da))/(2*step)
            max_eigen=max(max_eigen,float(np.max(abs(projected-slope))))
        require(max_cov<2e-10, model+": EFT1 mass beta rephasing covariance fails")
        require(max_eigen<2e-8,model+": mass beta eigenvalue transport fails")
        results.append({"model":model,"eft1_mass_beta_keys":list(b),
                        "uv_mass_beta_keys":[M1,M2,MASS],
                        "mu_beta_at_mu_zero":"6*lambdaS1S1bar2S2bar*mS1Sq + 6*lambdaS1barS2S2bar2*mS2Sq",
                        "zero_mu_alone_one_loop_closed":False,
                        "max_rephasing_covariance_residual":max_cov,
                        "max_eigenvalue_beta_residual":max_eigen,
                        "full_uv_yukawa_flavour_beta_checked":False})
    return {"results":results,"physical_mass_parameter_map_certified":False}


def audit_source_coverage(root:Path) -> dict:
    rgb=source(root,"RGE/running/rgbeta/T3RGBetaModel.wl")
    runner=source(root,"RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl")
    source_m=source(root,"Lagrangian/interactions/GeneralScalarOperators.wl")
    require("AddScalarMass[muT3GeneralS1S2,{Bar @ S1,Bar @ S2}" in runner and
            "SelfConjugate -> False" in runner and "eps[SU2L @ fund,#1,#2]" in runner,
            "Mixed RGBeta quadratic declaration changed")
    require('"Expression" -> PlusHc[contribution]' in source_m,
            "Matchete extra quadratic insertion changed")
    # Detect an actually registered Higgs mass, not occurrences in comments.
    addmass=re.findall(r"AddScalarMass\s*\[\s*([A-Za-z][A-Za-z0-9]*)\s*,\s*\{([^}]*)\}",rgb)
    higgs=[x for x,fields in addmass if re.search(r"\bH\b",fields)]
    # Main base beta tables are the full renormalizable-export contract.
    base_keys=set(re.findall(r'"([A-Za-z][A-Za-z0-9]*)"\s*->\s*Quiet\[BetaTerm\[',rgb))
    return {"rgbeta_registered_higgs_quadratic_mass_symbols":higgs,
            "rgbeta_beta_keys_include_higgs_quadratic":any(k in base_keys for k in higgs),
            "higgs_quadratic_mass_gap":not bool(higgs),
            "rgbeta_mixed_bilinear_source_declared":True,
            "matchete_plus_hc_source_declared":True,
            "matchete_mixed_mass_relative_sign_normalisation_certified":False}


def run(root:Path) -> dict:
    map_result=audit_operator_map(root)
    mass_math=mathematical_tests()
    mass_beta=audit_mass_betas(root)
    coverage=audit_source_coverage(root)
    inputs=["RGE/general/generated/general_rgbeta_quartic_registry.json",
            "docs/general_scalar_trial_operator_manifest.json",
            "output/general_t3_stage22/stage29_convention_uv_handoff.json",
            "output/general_t3_stage22/stage23_quartic_prefactor_candidates.json",
            "Lagrangian/interactions/GeneralScalarOperators.wl",
            "RGE/running/rgbeta/T3RGBetaModel.wl",
            "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl"]
    inputs += [f"output/general_t3_quadratic_eft1/{tag}_{stage}_quadratic_live_rgbeta.json"
               for tag in ("B_m1","C_m1") for stage in ("uv","eft1")]
    return {"status":"IndependentBCConventionAndMassRGAlgebraChecksPassedNotPhysicalCertificate",
            "map":map_result,"hessian_math":mass_math,"mass_beta":mass_beta,
            "coverage":coverage,"basis_certificate_dependency":"FINAL_WORKSTREAM_A_CERTIFICATE_REQUIRED",
            "source_sha256":{x:source_digest(root,x) for x in inputs},
            "physical_matchete_rgbeta_coupling_map_certified":False,
            "physical_mass_threshold_matching_certified":False,
            "renormalizable_rge_complete":False,
            "three_generation_wilson_evolution_certified":False}


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--output",type=Path,default=None)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args(argv)
    if a.self_test:
        print(json.dumps({"self_test":"PASS","errors":mathematical_tests()},indent=2))
        return 0
    result=run(a.root.resolve())
    if a.output is not None:
        out=a.output if a.output.is_absolute() else a.root/a.output
        out.parent.mkdir(parents=True,exist_ok=True)
        out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],
                      "models":[x["model"] for x in result["map"]["models"]],
                      "higgs_quadratic_mass_gap":result["coverage"]["higgs_quadratic_mass_gap"],
                      "physical_mapping_certified":False},indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
