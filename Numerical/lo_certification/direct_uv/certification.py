"""Additive, fail-closed source audit for direct T3 broken-phase UV inputs.

This module does not modify the existing direct_one_loop engine or EFT pipeline.
It distinguishes numerical/algebraic verification from physical UV certification.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.direct_one_loop.core import _complex, encoded

CLASSES = {
    "A": (1, 3, 2, (-4, -2, 0)),
    "B": (2, 2, 1, (-3, -1, 1)),
    "C": (2, 2, 3, (-3, -1, 1)),
    "D": (3, 1, 2, (-2, 0, 2)),
    "E": (3, 3, 2, (-4, -2, 0, 2)),
}

SOURCE_GUARDS = (
    "Lagrangian/interactions/T3Topology.wl",
    "Lagrangian/interactions/ScalarPotential.wl",
    "Lagrangian/model/T3Fields.wl",
    "Lagrangian/model/SU2Invariants.wl",
    "Numerical/direct_one_loop/core.py",
)

RESTRICTED_EXTRA_COUPLINGS = (
    "lambdaHHdagS2S2", "lambdaHHdagS1barS1bar",
    "lambdaHHdagS1barS2barCross", "lambdaH1Adj", "lambdaH2Adj",
    "lambdaS1Adj", "lambdaS2Adj",
)


def all_models():
    for kind, (s1, s2, f, alphas) in CLASSES.items():
        for alpha in alphas:
            suffix = f'{"m" if alpha < 0 else "p"}{abs(alpha)}'
            yield f"{kind}_{alpha:+d}", f"T3_dS1_{s1}_dS2_{s2}_dF_{f}_alpha_{suffix}"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _charge(dim: int, hypercharge: float, weight: int) -> float:
    return (dim - 1) / 2 - weight + hypercharge


def charge_bases(kind: str, alpha: int) -> dict[str, Any]:
    ds1, ds2, df, _ = CLASSES[kind]
    qs1 = np.array([_charge(ds1, alpha / 2, s) for s in range(ds1)])
    qs2 = np.array([_charge(ds2, (alpha + 2) / 2, s) for s in range(ds2)])
    qf = np.array([_charge(df, (alpha + 1) / 2, s) for s in range(df)])
    self_conjugate = (df % 2 == 1 and alpha == -1)
    # Provisional adapter's ordering: [F, F^c] for Dirac;
    # [F charge-paired components] for self-conjugate Majorana.
    qpsi = (np.tile(qf, 3) if self_conjugate else
            np.r_[np.tile(-qf, 3), np.tile(qf, 3)])
    return {"scalar": np.r_[qs1, qs2], "weyl": qpsi,
            "source_fermion_weights": qf, "self_conjugate": self_conjugate,
            "dimensions": (ds1, ds2, df)}


def witness_quadratic(A: np.ndarray, B: np.ndarray) -> dict[str, Any]:
    """Independently test the entire complex quadratic action against the Hessian.

    V(phi) = phi^dag A phi + Re(phi^T B phi), phi=(x+iy)/sqrt(2).
    No assumption about whether the scalar bilinear B is zero.
    """
    A = _complex(A, "A"); B = _complex(B, "B", A.shape)
    n = len(A)
    H = np.block([[A.real + B.real, -A.imag - B.imag],
                  [A.imag - B.imag, A.real - B.real]])
    xs = [np.arange(2*n, dtype=float)+1.,
          np.array([(-1.)**i/(i+1) for i in range(2*n)]),
          np.array([np.sin(i + .37) for i in range(2*n)])]
    discrepancies = []
    for x in xs:
        phi = (x[:n] + 1j*x[n:]) / np.sqrt(2.)
        v = (np.vdot(phi, A@phi) + np.real(phi.T@B@phi)).real
        reconstructed = .5*x@H@x
        discrepancies.append(float(abs(v-reconstructed)/max(abs(v),abs(reconstructed),1.)))
    return {"quadratic_action_max_relative_residual": max(discrepancies),
            "real_hessian_symmetry_residual": float(
                np.linalg.norm(H-H.T)/max(1.,np.linalg.norm(H))),
            "real_hessian_gev2": encoded(H)}


def _bad_entries(tensor: np.ndarray, allowed: np.ndarray, *, tolerance: float) -> list:
    amplitude = np.abs(tensor)
    scale = max(float(np.max(amplitude,initial=0)), 1e-100)
    return [tuple(int(j) for j in ij) for ij in
            np.argwhere((amplitude>tolerance*scale)&(~allowed))[:24]]


def validate_components(data: dict[str,Any], kind: str, alpha: int) -> dict[str, Any]:
    """Exact U(1)_em checks for ALL nonzero mass, bilinear and Yukawa entries."""
    base = charge_bases(kind, alpha)
    qf, qs = base["weyl"], base["scalar"]
    M = _complex(data["fermion_mass_gev"], "fermion mass", (len(qf),len(qf)))
    A = _complex(data["scalar_A_gev2"], "scalar A", (len(qs),len(qs)))
    B = _complex(data["scalar_B_gev2"], "scalar B", (len(qs),len(qs)))
    checks = {
        "fermion_mass_symmetric": bool(np.allclose(M,M.T,rtol=1e-11,atol=1e-12)),
        "scalar_A_hermitian": bool(np.allclose(A,A.conj().T,rtol=1e-11,atol=1e-12)),
        "scalar_B_symmetric": bool(np.allclose(B,B.T,rtol=1e-11,atol=1e-12)),
    }
    forbidden = {
        "fermion_mass": _bad_entries(M,np.isclose(qf[:,None]+qf[None,:],0),tolerance=1e-11),
        "scalar_A": _bad_entries(A,np.isclose(qs[:,None]-qs[None,:],0),tolerance=1e-11),
        "scalar_B": _bad_entries(B,np.isclose(qs[:,None]+qs[None,:],0),tolerance=1e-11),
    }
    vertex_details=[]
    channels=set()
    names=[]
    for row in data.get('vertices',[]):
        name = str(row['name']); names.append(name)
        record = {"name":name,"active_component_vertices":0,"field":0,
                  "conjugate":0,"charge_channels":[],"forbidden":[]}
        for field, sign in (("field",+1),("conjugate",-1)):
            if field not in row: continue
            v = _complex(row[field],name+'.'+field,(3,len(qf),len(qs)))
            allowed = np.isclose(qf[:,None]+sign*qs[None,:],0)[None,:,:]
            bad = _bad_entries(v,np.broadcast_to(allowed,v.shape),tolerance=1e-11)
            record['forbidden'] += [{"field":field,"indices":ij} for ij in bad]
            present = np.argwhere(np.abs(v)>1e-14*max(float(np.max(np.abs(v))),1e-100))
            record[field] += len(present)
            record['active_component_vertices'] += len(present)
            for _a, _f, s in present:
                q = float(qs[s]); channels.add(q)
                if q not in record['charge_channels']:record['charge_channels'].append(q)
        record['charge_channels'].sort()
        vertex_details.append(record)
    checks['distinct_vertex_names'] = (len(names)==len(set(names)))
    checks['vertices_present'] = (len(names)>0)
    checks['ordinary_source_topology_vertices_present'] = {'y1_HC','y2_HC'}.issubset(set(names))
    checks['all_mass_and_vertex_charge_conservation'] = (
        all(not x for x in forbidden.values()) and all(not r['forbidden'] for r in vertex_details))
    checks['real_hessian_action'] = witness_quadratic(A,B)['quadratic_action_max_relative_residual'] < 1e-12
    if not all(checks.values()) or not checks['all_mass_and_vertex_charge_conservation']:
        raise ValueError('Source candidate violates structural/charge/Hessian checks: '
                         + str({"checks":checks,"forbidden":forbidden,"vertices":vertex_details}))
    return {"status":"candidate_charge_and_quadratic_action_verified",
            "checks":checks,"forbidden_components":forbidden,
            "charge_bases": {k:(v.tolist() if isinstance(v,np.ndarray) else v)
                             for k,v in base.items()},
            "observed_scalar_charge_channels": sorted(channels),
            "vertex_dictionary": vertex_details,
            "hessian": {k:v for k,v in witness_quadratic(A,B).items() if k!="real_hessian_gev2"}}


def source_evidence(project:Path, short:str, key:str,
                    inp:dict, benchmark:dict, audit:dict|None,
                    portal:dict|None, paths:dict[str,Path]) -> dict:
    """Source consistency with exact file bytes, without a declaration of physical proof."""
    bfile=paths['benchmark']
    flags={"source_benchmark_matches":
           bool(inp.get('provenance',{}).get('benchmark_sha256')==digest(bfile)),
           "source_representation_matches": False,
           "input_at_source_base_scale":False,
           "source_guard_files_current":False,
           "prior_candidate_evidence_matches":False,
           "native_portal_tensor_geometry_representable":False,
           "cross_package_portal_dictionary_certified":False,
           "native_weyl_spinor_phases_certified":False,
           "complete_general_interactions_numerically_fixed":False}
    kind,alpha=short[0],int(short[2:]); d1,d2,df,_=CLASSES[kind]
    r=benchmark['representation']
    flags['source_representation_matches'] = (r['d_s1'],r['d_s2'],r['d_f'],r['alpha'])==(d1,d2,df,alpha)
    flags['input_at_source_base_scale'] = bool(np.isclose(
        float(inp['renormalization_scale_gev']),float(benchmark['base_state']['mu_gev']),rtol=1e-12))
    mismatches={}
    if audit:
        expected=audit.get('source_guard_sha256',{})
        for path in SOURCE_GUARDS:
            f=project/path
            if not f.is_file() or expected.get(path)!=digest(f):
                mismatches[path]={"audit":expected.get(path),"current":digest(f) if f.is_file() else None}
        flags['source_guard_files_current'] = not mismatches
        entry=next((e for e in audit.get('models',[]) if e.get('model')==short),None)
        flags['prior_candidate_evidence_matches'] = bool(entry and
            entry.get('benchmark_sha256')==digest(bfile) and
            entry.get('physical_uv_mass_certified') is False and
            audit.get('physical_models_certified')==0)
    pm=next((e for e in (portal or {}).get('models',[]) if e.get('model')==short),None)
    if pm:
        flags['native_portal_tensor_geometry_representable']=bool(
            len(pm.get('scalars',[]))==2 and
            all(s.get('full_tensor_identity_in_real_native_span') is True for s in pm['scalars']))
    # Intentional fail-closed independent physical certification:
    # the supplied source reports do NOT anchor four-spinor phases or the
    # RGBeta vs Matchete sign/normalisation dictionary.  They are not inferred.
    ordinary=benchmark['base_state']['ordinary']
    nonzero_extras=[]
    for name in RESTRICTED_EXTRA_COUPLINGS:
        val=ordinary.get(name,0)
        if isinstance(val,dict):val=complex(val['real'],val['imag'])
        if abs(complex(val))>0:nonzero_extras.append(name)
    extended_yukawa_allowed= alpha in (-2,0) and kind in ('A','D','E')
    return {"checks":flags,"source_guard_mismatches":mismatches,
            "native_general_extra_yukawa_allowed":extended_yukawa_allowed,
            "extra_yukawa_input":"not_present_in_restricted_frozen_benchmark" if extended_yukawa_allowed else "not_allowed_by_source_hypercharges",
            "additional_scalar_couplings_nonzero_without_verified_component_map":nonzero_extras,
            "physically_certified":False,
            "blocking_physical_conventions":[
                "independent native four-spinor to Weyl Feynman-rule phase/normalisation witness",
                "RGBeta to Matchete Higgs-portal and lambdaT3 coupling sign/normalisation map",
                "source-derived kinetic/mass charge-conjugation phase and Majorana normalisation",
                "numerical input and CG for allowed extra Yukawas/bilinears when nonzero",
            ]}
