"""Six-model A/D/E self-bilinear tensor, RGE and one-scale C5 batch audit.

A(-2), A(0), D(-2), D(0), E(-2), E(0).
Production source and fitted benchmark files are read only. Uses the exact
published one-scale Matchete kernel already saved in this repository; it is
NOT a new mass-eigenstate Matchete Match[] and NOT sequential EFT matching.

Run from repository root:
 python -u -X utf8 -m scripts.VerifySixT3SelfBilinearBatch --all --run-id six_self_v1
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import re
import subprocess
from typing import Any

import mpmath as mp
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
MODELS = ("A_-2", "A_+0", "D_-2", "D_+0", "E_-2", "E_+0")
REP = {"A": (1, 3, 2), "D": (3, 1, 2), "E": (3, 3, 2)}
MASS_NAMES = {-2: "muT3GeneralS2Sq", 0: "muT3GeneralS1Sq"}
EXTRA_NAMES = {-2: ("lambdaT3ExtraQ020101I1", "lambdaT3GeneralQ020101I1"),
                0: ("lambdaT3ExtraQ021010I1", "lambdaT3GeneralQ021010I1")}
SYMBOLIC = ROOT / "output/studies/general_scalar_extended/complete_scalar_symbolic"
OUTPUT = ROOT / "output/studies/general_scalar_extended/six_self_bilinear_batch"
NATIVE = ROOT / "Lagrangian/RunSixSelfBilinearNativeTensor.wl"
STAGE22 = ROOT / "output/general_t3_stage22/matchete_tensors.json"
REGISTRY = ROOT / "RGE/general/generated/general_rgbeta_quartic_registry.json"
FORMAT = "t3-six-self-bilinear-rge-tensor-c5-candidate-v1"


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError("T3 six-self audit: " + reason)


def load_json(path: Path) -> dict:
    require(path.is_file(), "missing file: " + str(path))
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    require(isinstance(obj, dict), "JSON object required: " + str(path))
    return obj


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def cjson(z: complex) -> dict:
    z = complex(z)
    require(np.isfinite([z.real, z.imag]).all(), "nonfinite complex number")
    return {"real": float(z.real), "imag": float(z.imag)}


def mjson(a: np.ndarray) -> dict:
    a = np.asarray(a, dtype=complex)
    require(a.shape == (3, 3) and np.isfinite(a).all(), "bad C5 matrix")
    return {"real": a.real.tolist(), "imag": a.imag.tolist(),
            "frobenius_norm": float(np.linalg.norm(a))}


def decode_matrix(payload: dict, shape: tuple[int, int]) -> np.ndarray:
    require(isinstance(payload, dict) and set(payload) == {"real", "imag"},
            "complex matrix schema")
    a = np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)
    require(a.shape == shape and np.isfinite(a).all(), "complex matrix dimensions")
    return a


def tag(model: str) -> str:
    return model.replace("-", "m").replace("+", "p")


def verify_summary(folder: Path) -> dict:
    o = load_json(folder / "summary.json")
    rows = o.get("models")
    require(isinstance(rows, list) and len(rows) == 16, "expected 16-model symbolic summary")
    statuses = {r.get("model"): r.get("status") for r in rows}
    require(len(statuses) == 16 and all(not r.get("errors") for r in rows),
            "symbolic batch incomplete or contains errors")
    require(all(statuses[m] == "SymbolicRGExportedMatchingUnverified" for m in MODELS),
            "self-model symbolic RG exports failed")
    require(all(statuses[m] == "SymbolicRGExportedMassBasisBlocked"
                for m in ("B_-1", "C_-1")), "unexpected mixed-model mass status")
    require(o.get("physical_matchete_rgbeta_mapping_certified") is False and
            o.get("full_sequential_wilson_matching_performed") is False,
            "symbolic summary certification flags changed")
    return {"models": 16, "statuses": statuses, "sha256": sha(folder / "summary.json")}


def verify_rg_exports(folder: Path, model: str) -> dict:
    alpha = int(model.split("_")[1])
    sym = MASS_NAMES[alpha]
    row = folder / tag(model)
    status, lag = load_json(row / "status.json"), load_json(row / "lagrangian.json")
    require(status.get("status") == "SymbolicRGExportedMatchingUnverified" and not status.get("errors") and
            lag.get("status") == "Success" and lag.get("quadratic_included") == 1 and
            lag.get("quadratic_expected") == 1 and lag.get("matching_performed") is False,
            model + ": incomplete UV Lagrangian")
    result = {"mass_symbol": sym}
    for stage in ("uv", "eft1"):
        p = row / (stage + "_rgbeta.json")
        o = load_json(p)
        meta = o.get("metadata", {})
        betas = o.get("report_betas", {})
        require(o.get("format") == "general-t3-quadratic-rgbeta-experimental-v1" and
                o.get("status") == "SuccessCandidateMassBeta" and
                o.get("physical_mapping_certified") is False and
                o.get("mass_eigenstate_threshold_matching_performed") is False and
                meta.get("GeneralScalarModel") == model and
                meta.get("GeneralScalarStage") == stage.upper() and
                meta.get("QuadraticKind") == ("self_s1" if alpha == 0 else "self_s2") and
                meta.get("QuadraticMassSymbol") == sym and
                meta.get("GeneralScalarQuarticBetaBasis") == "RGBetaFullQuarticProjection" and
                isinstance(betas, dict) and {"mS1Sq", "mS2Sq", sym} <= set(betas) and
                set(o.get("beta_symbols", [])) == set(betas),
                model + ": malformed " + stage + " RGE export")
        require(all(isinstance(x, str) and x.strip() and
                    not any(bad in x for bad in ("BetaTerm[", "Missing[", "$Failed", "$Aborted", "Indeterminate"))
                    for x in betas.values()), model + ": unresolved " + stage + " beta")
        require(sym in betas[sym], model + ": mass beta lacks its own mass parameter")
        result[stage] = {"beta_count": len(betas), "mass_beta": betas[sym],
                         "source_sha256": sha(p)}
    return result


def verify_convention_sources(model: str) -> dict:
    alpha = int(model.split("_")[1])
    mat = ROOT / "Lagrangian/interactions/GeneralScalarOperators.wl"
    rgb = ROOT / "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl"
    mt, rt = mat.read_text(encoding="utf-8"), rgb.read_text(encoding="utf-8")
    require('"Expression" -> PlusHc[contribution]' in mt and
            'contribution = coupling[] product' in mt and
            'V2=mu2*(Bar[S])^2+h.c.' in rt,
            "Matchete Lagrangian/RGBeta scalar potential source convention changed")
    symbol = MASS_NAMES[alpha]
    field = "S1" if alpha == 0 else "S2"
    pattern = r"AddScalarMass\[" + symbol + r",\{Bar @ " + field + r",Bar @ " + field + r"\}"
    require(re.search(pattern, rt) is not None, "RGBeta holomorphic mass registration changed")
    source = ROOT / "Numerical/experimental/GeneralT3LagrangianConvention.py"
    require("RGBETA_MBL_EQUATION" in source.read_text(encoding="utf-8"),
            "missing RGBeta Lagrangian sign source")
    registry = load_json(REGISTRY)["models"][model]
    mat_name, rgb_name = EXTRA_NAMES[alpha]
    x = [q for q in registry["new_rgbeta_symbols"] if q["name"] == rgb_name]
    expected = (["Bar @ H", "Bar @ H", "Bar @ S1", "Bar @ S2"] if alpha == -2
                else ["Bar @ H", "Bar @ H", "S1", "S2"])
    require(len(x) == 1 and x[0]["fields"] == expected and
            x[0]["matchete_basis_name"] == mat_name[:-1] + "#" and
            x[0]["mapping_status"] == "NotNormalised", model + ": unexpected quartic registry")
    return {"matchete_mass_operator": "+muMat Q+h.c. in Lagrangian",
            "rgbeta_mass_operator": "+muRGB Q+h.c. in scalar potential",
            "conditional_mass_map": "muRGB = -kNative*muMat after component-basis alignment",
            "conditional_extra_quartic_map": "lambdaMat = -lambdaRGB after component-basis alignment",
            "matchete_extra": mat_name, "rgbeta_extra": rgb_name,
            "matchete_source_sha256": sha(mat), "rgbeta_source_sha256": sha(rgb),
            "registry_sha256": sha(REGISTRY),
            "triplet_component_basis_mapping_certified": False}


def parse_native_number(raw: str) -> complex:
    """Restricted numeric Matchete InputForm; no symbolic evaluation or code."""
    require(isinstance(raw,str) and re.fullmatch(r"[0-9ISqrt\s+\-*/^()\[\]]+",raw) is not None,
            "unsafe or nonnumeric native CG entry")
    text=raw.replace("Sqrt[","sqrt(").replace("]",")").replace("^","**")
    expr=sp.sympify(text,locals={"I":sp.I,"sqrt":sp.sqrt})
    require(not expr.free_symbols and expr.is_finite is True,"native CG entry not finite number")
    v=complex(sp.N(expr,24))
    require(np.isfinite([v.real,v.imag]).all(),"native CG entry overflow")
    return v


def validate_native_result(obj: dict, model: str) -> complex:
    alpha = int(model.split("_")[1])
    d = REP[model[0]][0 if alpha == 0 else 1]
    species = ("S1dag" if alpha == 0 else "S2dag")
    require(obj.get("format") == "t3-six-self-native-bilinear-tensor-v1" and
            obj.get("status") == "NativeMassTensorExtracted" and
            obj.get("model") == model and obj.get("dimension") == d and
            obj.get("species") == [species, species] and
            obj.get("native_metric_proportional_identity") is True and
            obj.get("mass_operator_in_full_lagrangian") is True and
            len(obj.get("operator_names", [])) == 1 and
            obj.get("rgbeta_matchete_vertex_map_certified") is False,
            model + ": invalid native tensor certificate")
    raw = obj.get("native_metric_inputform")
    require(isinstance(raw, list) and len(raw) == d and
            all(isinstance(a, list) and len(a) == d for a in raw),
            "native CG dimension mismatch")
    tensor = np.asarray([[parse_native_number(x) for x in a] for a in raw], dtype=complex)
    factor = complex(parse_native_number(obj["native_metric_factor_inputform"]))
    require(abs(factor) > 1e-14 and np.allclose(tensor, factor*np.eye(d),rtol=2e-12,atol=2e-13),
            "native CG not proportional to unit matrix")
    return factor


def triangle_integral(mf: float, ms1_sq: float, ms2_sq: float, *, dps: int = 52) -> mp.mpf:
    """Independent Euclidean loop: ∫_0∞ t dt/[(t+MF²)(t+m1²)(t+m2²)]."""
    require(min(mf,ms1_sq,ms2_sq) > 0 and np.isfinite([mf,ms1_sq,ms2_sq]).all(),
            "triangle mass arguments invalid")
    with mp.workdps(dps):
        mass = mp.mpf(str(mf))
        a = mp.mpf(str(ms1_sq))/mass**2
        b = mp.mpf(str(ms2_sq))/mass**2
        return +(mp.quad(lambda t: t/((t+1)*(t+a)*(t+b)),
                         [0,mp.mpf("0.00001"),mp.mpf("0.001"),mp.mpf("0.1"),1,10,mp.inf])/mass**2)


def kernel_weights(source: Path, mf: np.ndarray, ms1sq: float, ms2sq: float,
                   split_species: str, split_sq: tuple[float,float],
                   *, dps: int = 52) -> tuple[np.ndarray,np.ndarray,dict]:
    """Exact one-generation source C5 divided by Yukawas, with lambdaT3=1.

    Compare each absolute mass-dependent value against an independent integral.
    The SU2 factor is determined by the *source* and checked to be identical
    at six separate mass points, not advertised as independently certified.
    """
    from RGE.matching.MatcheteC5Parsing import parse_matchete_c5
    from RGE.matching.WeinbergFlavorMatching import extract_one_generation_loop_kernel
    raw=source.read_text(encoding="utf-8-sig")
    require(all(t in raw for t in ("Coupling[MF,", "Coupling[MS1,", "Coupling[MS2,",
                                     "Coupling[lambdaT3,", "Coupling[y1,", "Coupling[y2,")),
            "mass-sensitive one-scale Matchete source not complete")
    expr=extract_one_generation_loop_kernel(parse_matchete_c5(raw))
    symbols = tuple(sp.Symbol(x) for x in ("MF","MS1","MS2","lambdaT3","hbar"))
    require(expr.free_symbols <= set(symbols) and
            set(symbols[:3]) | {symbols[3]} <= expr.free_symbols,
            "Matchete one-scale loop kernel has unexpected symbols")
    require(sp.simplify(expr.subs(symbols[3], 0)) == 0,
            "source C5 not zero at lambdaT3=0")
    require(sp.simplify(expr - symbols[3]*expr.subs(symbols[3],1)) == 0,
            "source C5 not lambdaT3-linear")
    fn=sp.lambdify(symbols,expr,modules="mpmath")
    source_weights=[[],[]]
    factors=[]; comparisons=[]
    with mp.workdps(dps):
        hbar=1/(16*mp.pi**2)
        for idx,m0 in enumerate(mf):
            m=mp.mpf(str(float(m0)))
            for branch,mass_sq in enumerate(split_sq):
                a=split_sq[branch] if split_species == "S1" else ms1sq
                b=split_sq[branch] if split_species == "S2" else ms2sq
                exact=mp.mpc(fn(m,mp.sqrt(mp.mpf(str(a))),mp.sqrt(mp.mpf(str(b))),1,hbar))
                integral=triangle_integral(float(m0),a,b,dps=dps)
                reference=hbar*m*integral
                require(abs(reference)>mp.mpf("1e-80") and abs(exact)>mp.mpf("1e-80"),
                        "one-scale loop is unexpectedly zero")
                fac=exact/reference
                require(abs(fac.imag)<mp.mpf("1e-10")*max(abs(fac),mp.mpf("1e-50")),
                        "unexpected complex source SU2 coefficient")
                factors.append(fac)
                comparisons.append({"heavy_generation":idx+1,"mass_gev":float(m0),
                                    "branch":"plus" if branch==0 else "minus",
                                    "source_kernel_gev_inverse":float(exact.real),
                                    "integral_hbar_MF_gev_inverse":float(reference),
                                    "source_over_integral":cjson(complex(fac))})
                source_weights[branch].append(float(exact.real))
        k=factors[0]
        scatter=max(float(abs(x-k)/max(abs(x),abs(k))) for x in factors)
        require(scatter < 2e-9, "Matchete kernel fails independent triangle mass-dependence consistency")
        if str(source).replace('\\','/').endswith('T3_A_alpha_m2/data/c5_coefficient.txt'):
            require(abs(k+2/mp.sqrt(3))<mp.mpf('2e-9'),
                    "A(-2) absolute source group coefficient changed")
    result={"kernel_mass_dependence_verified":True,
            "source_su2_coefficient_from_matchete":cjson(complex(k)),
            "coefficient_constancy_max_relative_error":scatter,
            "independent_triangle_integral_points":6,
            "comparison":comparisons,
            "absolute_su2_coefficient_independently_derived":False,
            "source_sha256":sha(source)}
    return np.asarray(source_weights[0]),np.asarray(source_weights[1]),result


def lift_flavour(weights: np.ndarray, y1: np.ndarray, y2: np.ndarray) -> np.ndarray:
    w=np.asarray(weights,dtype=complex)
    require(w.shape==(3,) and y1.shape == y2.shape == (3,3),"flavour matrix shapes")
    a,b=y1.conjugate(),y2.conjugate()
    c=(a*w)@b.T+(b*w)@a.T
    require(np.isfinite(c).all() and np.allclose(c,c.T,rtol=1e-11,atol=1e-27),
            "C5 not complex symmetric")
    return c


def candidate_channels(plus: np.ndarray, minus: np.ndarray, *, lam_t3: complex,
                       extra_rgb: complex, mu_rgb: complex, alpha: int,
                       cg_ratio: complex = 1+0j, matchete_sign: int = -1) -> dict:
    """Mass-split normal/anomalous spectral continuation; conditional map.

    alpha=-2 uses <S2 S2>, phase mu/abs(mu);
    alpha= 0 uses <S1dag S1dag>, phase conjugate(mu)/abs(mu).
    """
    require(alpha in (-2,0) and matchete_sign in (-1,1),"invalid spectral sign/alpha")
    p,m=np.asarray(plus,dtype=complex),np.asarray(minus,dtype=complex)
    require(p.shape==m.shape==(3,3) and np.isfinite(p).all() and np.isfinite(m).all(),
            "invalid normal/anomalous C5 unit kernels")
    z=complex(mu_rgb)
    phase = 0j if z==0 else (z/abs(z) if alpha==-2 else z.conjugate()/abs(z))
    normal=complex(lam_t3)*.5*(p+m)
    abnormal=matchete_sign*complex(extra_rgb).conjugate()*phase*complex(cg_ratio)*.5*(p-m)
    combined=normal+abnormal
    require(all(np.isfinite(t).all() and np.allclose(t,t.T,rtol=1e-11,atol=1e-26)
                for t in (normal,abnormal,combined)), "bad symmetric C5 channels")
    den=float(np.linalg.norm(normal))
    return {"normal":normal,"anomalous":abnormal,"total":combined,
            "relative_anomalous_norm":float(np.linalg.norm(abnormal)/den) if den>0 else None,
            "phase":phase}


def verify_geometry_and_mass_map(model: str, k: complex,
                                 masssq: float, mu: complex) -> dict:
    from Numerical.experimental.T3SelfBilinearMassBasis import (
        SelfBilinearMassBasis, invariant_metric)
    alpha=int(model.split('_')[1]); d=REP[model[0]][0 if alpha==0 else 1]
    field="S1" if alpha==0 else "S2"
    basis=SelfBilinearMassBasis.build(field,d,masssq,mu)
    require(not basis.tachyonic,"tachyonic reference split")
    # Compare native Matchete Cartesian k*delta and RGBeta spherical C
    # only as an SU(2) metric statement; this is not all vertex intertwiners.
    C=invariant_metric(d)
    W=basis.intertwiner
    require(np.allclose(W.T@C@W,np.eye(d),atol=2e-13),"spherical/real metric mismatch")
    mat_mu=-mu/k
    max_err=0.
    rng=np.random.default_rng(790+d+alpha)
    for _ in range(20):
        r,i=rng.normal(size=(2,d))
        s=basis.fields(r,i)
        native=W.conjugate().T@s
        vmat=masssq*np.vdot(native,native).real-2*np.real(mat_mu*k*np.dot(native.conjugate(),native.conjugate()))
        vrgb=basis.potential_old(s)
        vdiag=basis.potential_mass_basis(r,i)
        err=max(abs(vmat-vrgb),abs(vdiag-vrgb))/max(abs(vdiag),1.)
        max_err=max(max_err,float(err))
        require(err<2e-12,"Matchete/RGBeta self mass metric convention mismatch")
    return {"native_bilinear_coefficient":cjson(k),"matchete_mass_coefficient":cjson(mat_mu),
            "rgbeta_mass_coefficient":cjson(mu),
            "eigen_mass_squared_gev2":list(basis.eigen_masses_squared),
            "mass_potential_max_relative_error":max_err,
            "basis_metric_intertwiner_verified":True,
            "all_yukawa_and_quartic_vertices_mapped":False}


def reference_data(model: str) -> dict:
    family,alpha_s=model.split('_')
    alpha=int(alpha_s); d1,d2,df=REP[family]
    token=('m' if alpha<0 else 'p')+str(abs(alpha))
    path=ROOT/f'Numerical/benchmarks/fitted_reference/T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{token}.json'
    obj=load_json(path)
    rep=obj.get('representation',{})
    require((rep.get('d_s1'),rep.get('d_s2'),rep.get('d_f'),rep.get('alpha'))==(d1,d2,df,alpha),
            "reference benchmark representation mismatch")
    o=obj['base_state']['ordinary']
    m=decode_matrix(o['MF'],(3,3))
    require(np.max(abs(m-np.diag(np.diag(m))))<1e-10 and
            np.max(abs(np.diag(m).imag))<1e-10 and np.min(m.diagonal().real)>0,
            "benchmark heavy MF not positive real diagonal")
    y1=decode_matrix(o['y1'],(3,3));y2=decode_matrix(o['y2'],(3,3))
    return {"config":obj,"ordinary":o,"file":str(path),"sha256":sha(path),
            "mf":m.diagonal().real,"y1":y1,"y2":y2}


def native_tensor_extract(model: str, out: Path, *, exe: str) -> dict:
    require(NATIVE.is_file(), "missing native tensor extraction Wolfram source")
    family,alpha=model.split('_')
    env=os.environ.copy()
    env['T3_EXPERIMENTAL_GENERAL_QUADRATICS']='1'
    env['T3_EXPERIMENTAL_GENERAL_QUARTICS']='1'
    p=out/'native_bilinear.json'
    cmd=[exe,'-file',str(NATIVE),str(p),family,str(int(alpha))]
    res=subprocess.run(cmd,cwd=ROOT,env=env,check=False,capture_output=True,text=True)
    (out/'native_matchete.log').write_text(res.stdout+'\n'+res.stderr,encoding='utf-8')
    require(res.returncode==0 and p.is_file(),
            model + ": live Matchete native tensor failed; see native_matchete.log: "
            +res.stdout[-600:]+res.stderr[-600:])
    return load_json(p)


def run_one(model: str, dest: Path, symbolic_run: Path, *, exe: str,
            mu_fraction: float, mu_phase: float, extra: complex,
            continuous_run: Path | None = None) -> dict:
    from Numerical.experimental.T3SelfBilinearAnomalousKernel import tensor_ratio
    from Numerical.experimental.T3SelfBilinearMassBasis import SelfBilinearMassBasis
    from Numerical.experimental.T3SelfBilinearMassDependentKernel import source_path
    model_dir=dest/tag(model)
    model_dir.mkdir()
    alpha=int(model.split('_')[1]);d1,d2,_=REP[model[0]]
    rge=verify_rg_exports(symbolic_run,model)
    convention=verify_convention_sources(model)
    data=reference_data(model)
    o=data['ordinary']
    species='S1' if alpha==0 else 'S2'
    chosen_m=float(o['mS1Sq' if alpha==0 else 'mS2Sq'])
    mu=chosen_m*mu_fraction*complex(math.cos(mu_phase),math.sin(mu_phase))
    native=native_tensor_extract(model,model_dir,exe=exe)
    k=validate_native_result(native,model)
    mass=verify_geometry_and_mass_map(model,k,chosen_m,mu)
    b=SelfBilinearMassBasis.build(species,d1 if alpha==0 else d2,chosen_m,mu)
    stage=load_json(STAGE22)
    cg=ROOT/f'output/group_factors/T3_{model[0]}_uv_cg_registry.json'
    checked=tensor_ratio(stage,load_json(cg),model,d1,d2)
    require(checked['matchete_extra_symbol']==EXTRA_NAMES[alpha][0] and
            checked['rgbeta_extra_symbol']==EXTRA_NAMES[alpha][1] and
            abs(complex(checked['ratio'])-1)<2e-12 and
            checked['relative_tensor_residual']<2e-12 and
            checked['native_self_metric_residual']<2e-12,
            model+": native Matchete anomalous HH tensor not aligned with T3 topology")
    source=source_path(ROOT,model)
    require(source.is_file(),"missing native one-scale C5 source: "+str(source))
    wp,wm,loop=kernel_weights(source,data['mf'],float(o['mS1Sq']),float(o['mS2Sq']),
                               species,tuple(b.eigen_masses_squared))
    plus=lift_flavour(wp,data['y1'],data['y2'])
    minus=lift_flavour(wm,data['y1'],data['y2'])
    lam=complex(float(o['lambdaT3']['real']),float(o['lambdaT3']['imag']))
    chan=candidate_channels(plus,minus,lam_t3=lam,extra_rgb=extra,mu_rgb=mu,
                            alpha=alpha,cg_ratio=complex(checked['ratio']),matchete_sign=-1)
    opposite=candidate_channels(plus,minus,lam_t3=lam,extra_rgb=extra,mu_rgb=mu,
                            alpha=alpha,cg_ratio=complex(checked['ratio']),matchete_sign=1)
    require(np.allclose(chan['normal'],opposite['normal'],rtol=0,atol=0) and
            np.allclose(chan['anomalous'],-opposite['anomalous'],rtol=1e-13,atol=1e-30),
            "sign-sensitivity scenario inconsistent")
    zero=candidate_channels(plus,minus,lam_t3=lam,extra_rgb=0,mu_rgb=mu,alpha=alpha)
    require(np.linalg.norm(zero['anomalous'])==0,"zero quartic yields nonzero anomalous contribution")
    if abs(mu)>0 and abs(extra)>0:
        require(np.linalg.norm(chan['anomalous'])>0,"nonzero mass splitting/quartic anomalous channel vanished")
    output={
        "format":FORMAT,"status":"SixSelfCandidateOneScaleTensorAndRGVerified",
        "model":model,"inputs_are_illustrative_not_fitted":True,
        "reference_benchmark_sha256":data['sha256'],
        "reference_heavy_masses_gev":data['mf'].tolist(),
        "reference_scalar_masses_squared_gev2":{"S1":float(o['mS1Sq']),"S2":float(o['mS2Sq'])},
        "reference_lambda_t3":cjson(lam),"extra_quartic_rgbeta":cjson(extra),
        "mu2_fraction":mu_fraction,"mu2_phase_rad":mu_phase,
        "mass":mass,"symbolic_rg":rge,"source_conventions":convention,
        "native_tensors":{"anomalous_to_topology_factor":cjson(complex(checked['ratio'])),
                          "relative_tensor_residual":checked['relative_tensor_residual'],
                          "native_self_metric_residual":checked['native_self_metric_residual'],
                          "stage22_sha256":sha(STAGE22),"group_cg_sha256":sha(cg),
                          "native_bilinear_sha256":sha(model_dir/'native_bilinear.json')},
        "one_scale_kernel":loop,
        "normal_c5":mjson(chan['normal']),"anomalous_c5":mjson(chan['anomalous']),
        "total_c5":mjson(chan['total']),
        "relative_anomalous_to_normal":chan['relative_anomalous_norm'],
        "anomalous_scalar_phase":cjson(chan['phase']),
        "opposite_sign_scenario_total_c5":mjson(opposite['total']),
        "opposite_sign_scenario_not_claimed_physical":True,
        "rgbeta_matchete_full_component_basis_map_certified":False,
        "quartic_source_potential_relative_sign_verified_given_same_component_fields":True,
        "independent_one_loop_mass_eigen_matchete_matching_performed":False,
        "sequential_EFT_Wilson_matching_performed":False,
        "UV_and_EFT1_running_evaluated_on_this_benchmark":False,
        "neutrino_masses_fitted_or_predicted":False,
        "scope":"Conditional source-convention spectral one-scale continuation with native Matchete SU2 tensors; RGEs checked for completeness but not evaluated as a trajectory."
    }
    if continuous_run is not None:
        # Reuse the EXISTING conditional mass-dependent one-scale matching
        # and low-energy SM+Weinberg evolution on a previously generated
        # UV/EFT1 continuous trajectory. This is independent of the illustrative
        # static benchmark above; the two parameter sets MUST NOT be conflated.
        from scripts.RunT3SelfBilinearMassDependentC5 import calculate
        evolved = calculate(model,source_root=continuous_run,
                            output_root=model_dir,sign=-1,resume=False)
        require(evolved.get('status')=='CandidateMassDependentOneScaleNormalPlusAnomalousHH' and
                evolved.get('physical_full_general_scalar_matching_certified') is False and
                evolved.get('full_sequential_wilson_matching_performed') is False,
                'unexpected continuous candidate physics-certification metadata')
        output['continuous_RG_candidate']={
          'source_run_dir':str(continuous_run),
          'source_trajectory_sha256':sha(continuous_run/tag(model)/'continuous_trajectory.json'),
          'one_scale_matching_status':evolved['status'],
          'result_file':str(model_dir/(tag(model)+'_massdependent_c5.json')),
          'relative_low_energy_mnu_change_vs_normal':evolved.get('relative_low_energy_mnu_change_vs_normal'),
          'anomalous_over_normal_c5':evolved.get('anomalous_over_normal_c5'),
          'low_energy_takagi_masses_ev_descending':evolved.get('takagi_masses_ev_descending'),
          'is_conditionally_mapped_not_certified':True,
        }
        output['UV_and_EFT1_running_evaluated_on_this_benchmark']=False
        output['preexisting_UV_and_EFT1_trajectory_reused_for_separate_candidate']=True
    (model_dir/'certificate.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return output


def main(argv: list[str] | None = None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--all',action='store_true')
    group.add_argument('--model',choices=MODELS)
    ap.add_argument('--symbolic-run-id',default='all16_v1')
    ap.add_argument('--run-id',required=True)
    ap.add_argument('--output-root',type=Path,default=OUTPUT)
    ap.add_argument('--wolframscript',default='wolframscript')
    ap.add_argument('--mu2-fraction',type=float,default=0.01)
    ap.add_argument('--mu2-phase',type=float,default=0.35)
    ap.add_argument('--extra-rgbeta-real',type=float,default=0.02)
    ap.add_argument('--extra-rgbeta-imag',type=float,default=0.01)
    ap.add_argument('--continuous-run-id',default=None,
                    help='Also evaluate existing evolved trajectory via existing conditional C5 runner')
    ns=ap.parse_args(argv)
    require(re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_-]{0,63}',ns.run_id) is not None,
            'unsafe --run-id')
    require(re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_-]{0,63}',ns.symbolic_run_id) is not None,
            'unsafe --symbolic-run-id')
    require(0<ns.mu2_fraction<0.45 and np.isfinite([ns.mu2_fraction,ns.mu2_phase,
          ns.extra_rgbeta_real,ns.extra_rgbeta_imag]).all(),"invalid illustrative parameters")
    sym=SYMBOLIC/ns.symbolic_run_id
    summary=verify_summary(sym)
    chosen=MODELS if ns.all else (ns.model,)
    continuous_run=None
    if ns.continuous_run_id is not None:
        require(re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_-]{0,63}',ns.continuous_run_id) is not None,
                'unsafe --continuous-run-id')
        continuous_run=ROOT/'output/studies/general_scalar_extended/continuous_full_scalar'/ns.continuous_run_id
        for m in chosen:
            require((continuous_run/tag(m)/'continuous_trajectory.json').is_file(),
                    'missing preexisting continuous trajectory for '+m)
    out=ns.output_root/ns.run_id
    require(not out.exists(),"output already exists; use new --run-id: "+str(out))
    out.mkdir(parents=True)
    results=[]
    errors=[]
    for model in chosen:
        try:
            record=run_one(model,out,sym,exe=ns.wolframscript,
                           mu_fraction=ns.mu2_fraction,mu_phase=ns.mu2_phase,
                           extra=complex(ns.extra_rgbeta_real,ns.extra_rgbeta_imag),
                           continuous_run=continuous_run)
            results.append({"model":model,"status":record['status'],
                            "native_metric_factor":record['mass']['native_bilinear_coefficient'],
                            "source_SU2_prefactor":record['one_scale_kernel']['source_su2_coefficient_from_matchete'],
                            "loop_relative_scatter":record['one_scale_kernel']['coefficient_constancy_max_relative_error'],
                            "anomalous_over_normal":record['relative_anomalous_to_normal'],
                            "certificate":str(out/tag(model)/'certificate.json'),
                            "continuous_low_energy_mnu_relative_change":
                               (record.get('continuous_RG_candidate') or {}).get('relative_low_energy_mnu_change_vs_normal')})
            print(model,': checked, |anomalous|/|normal| =',
                  format(record['relative_anomalous_to_normal'],'.9g'),flush=True)
        except (ValueError,TypeError,KeyError,RuntimeError,OSError,ArithmeticError,
                np.linalg.LinAlgError) as exc:
            errors.append({"model":model,"error":f'{type(exc).__name__}: {exc}'})
            print(model,': FAILED:',str(exc),flush=True)
    report={"format":"t3-six-self-audit-summary-v1",
            "status":"CandidateSixSelfOneScaleBatchComplete" if not errors and len(results)==len(chosen)
                     else "IncompleteOrFailed",
            "models_requested":list(chosen),"models_passed":len(results),
            "models_failed":len(errors),"results":results,"errors":errors,
            "symbolic_16_model_inventory":summary,
            "triplet_component_basis_fully_mapped":False,
            "new_mass_eigenstate_Matchete_match_performed":False,
            "sequential_wilson_matching_performed":False,
            "RGE_trajectory_integrated_in_this_batch":False,
            "existing_RGE_trajectory_reused":continuous_run is not None,
            "conditional_low_energy_C5_evaluated":continuous_run is not None,
            "fitted_neutrino_mass_prediction":False}
    (out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    with (out/'summary.csv').open('w',newline='',encoding='utf-8') as f:
        wr=csv.DictWriter(f,fieldnames=['model','status','anomalous_over_normal',
                      'loop_relative_scatter','certificate','continuous_low_energy_mnu_relative_change'])
        wr.writeheader()
        wr.writerows([{k:r.get(k) for k in wr.fieldnames} for r in results])
    print('Summary:',out/'summary.json',flush=True)
    return 0 if not errors else 2

if __name__=='__main__':
    raise SystemExit(main())
