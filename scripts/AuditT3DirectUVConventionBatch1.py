#!/usr/bin/env python3
"""Read-only, fail-closed all-16 direct-UV source-to-component convention audit.

This audits the EXISTING restricted frozen point. It does not derive a native
Matchete 4-spinor -> Weyl phase or identify RGBeta/Matchete couplings.
Successful algebraic checks MUST NOT set validated_uv_component_conventions.

Run from project root:
  python scripts/AuditT3DirectUVConventionBatch1.py --project-root .

Output: output/general_observable/lo_direct_one_loop/convention_batch1/
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
from itertools import product
import sys
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Independent consumer of the existing CG reports; no import from the
# provisional builder so the array contraction has an independent code path.
CLASSES = {'A':(1,3,2,(-4,-2,0)), 'B':(2,2,1,(-3,-1,1)),
           'C':(2,2,3,(-3,-1,1)), 'D':(3,1,2,(-2,0,2)),
           'E':(3,3,2,(-4,-2,0,2))}
FLAGS = {'T3Y1CG':(True,True,False),
         'T3Y2CG':(True,False,True),
         'T3MixCG':(False,False,False,True)}


def models():
    for label,(d1,d2,df,alphas) in CLASSES.items():
        for alpha in alphas:
            m='m' if alpha<0 else 'p'
            yield (f'{label}_{alpha:+d}',
                   f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{m}{abs(alpha)}')


def _number(value):
    if isinstance(value,dict):
        return complex(value['real'],value['imag'])
    return complex(value)


def _matrix(value):
    if isinstance(value,dict) and set(value)=={'real','imag'}:
        return np.asarray(value['real'],float)+1j*np.asarray(value['imag'],float)
    return np.asarray(value,complex)


def _spnum(s):
    return sp.sympify(str(s),locals={'I':sp.I,'Sqrt':sp.sqrt})


def _rotations(analysis):
    us={1:sp.eye(1)}
    for rec in analysis['representations']:
        require(rec.get('status')=='generator_algebra_verified',
                'SU2 candidate generator report did not pass')
        d=int(rec['dimension'])
        U=sp.Matrix([[_spnum(e) for e in row]
                     for row in rec['U_columns_are_weight_eigenvectors']])
        require(U.shape==(d,d) and all(sp.simplify(z)==0 for z in
                U.H*U-sp.eye(d)), 'SU2 generator eigenvector rotation not unitary')
        us[d]=U
    require(set(us)=={1,2,3},'Expected SU2 dimensions 1,2,3')
    return us


def _tensor(vertex, flags, us):
    dims=tuple(int(d) for d in vertex['dimensions'])
    require(len(dims)==len(flags),'Native tensor index rank mismatch')
    maps=[us[d].conjugate() if conj else us[d]
          for d,conj in zip(dims,flags)]
    out=np.zeros(dims,complex)
    for rec in vertex['components']:
        old=tuple(int(i)-1 for i in rec['indices'])
        require(all(0<=i<d for i,d in zip(old,dims)),
                'Native tensor CG index out of range')
        val=_spnum(rec['cg'])
        for new in product(*(range(d) for d in dims)):
            term=val
            for axis in range(len(dims)):
                term*=maps[axis][old[axis],new[axis]]
            out[new]+=complex(sp.N(term,16))
    return out

OUT = Path('output/general_observable/lo_direct_one_loop/convention_batch1')
CG_DIR = Path('output/general_observable/lo_direct_one_loop/source_audits')
CONDITIONAL = 'restricted_UV_action_algebra_checked_physical_coupling_dictionary_PENDING'
EPS = 2e-10


def require(cond: bool, why: str) -> None:
    if not cond:
        raise ValueError(why)


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f'Missing required evidence: {path}')
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), f'Expected JSON object: {path}')
    return obj


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root)).replace('\\', '/')
    except ValueError:
        return str(path)


def complex_array(obj: Any, label: str) -> np.ndarray:
    z = _matrix(obj)
    require(np.isfinite(z).all(), f'{label}: nonfinite number')
    return np.asarray(z, dtype=np.complex128)


def near(actual: np.ndarray, expected: np.ndarray, label: str,
         tol: float=EPS) -> float:
    a = np.asarray(actual, dtype=complex)
    b = np.asarray(expected, dtype=complex)
    require(a.shape == b.shape, f'{label}: shapes {a.shape} != {b.shape}')
    residual = float(np.linalg.norm(a-b) / max(np.linalg.norm(b), 1.0))
    require(residual < tol, f'{label}: residual {residual:.3e} > {tol:.3e}')
    return residual


def check_source_contracts(root: Path) -> dict[str, str]:
    """Reject changed native operator orientations; not a spinor proof."""
    required = {
        'Lagrangian/interactions/T3Topology.wl': (
            'useCConj = which === 1;',
            'scalarBar = which === 2;',
            'y[p, r] NCM[Bar[l[i, p]], f] scalar CG[cg, labels]',
            'lambdaT3[] H[i] H[j] s1 s2 CG[T3MixCG, labels]',
            '"Expression" -> PlusHc[',
        ),
        'Lagrangian/model/SU2Invariants.wl': (
            'T3Y1CG, {2, dF, d1}, {True, True, False}',
            'T3Y2CG, {2, dF, d2}, {True, False, True}',
            '{2, 2, d1, d2}', '{False, False, False, True}',
        ),
        'Lagrangian/model/T3Fields.wl': (
            'selfConj = TrueQ[PossibleZeroQ[f["Y"]]] && OddQ[f["SU2"]]',
            'DefineField[', 'NewFermion, Fermion,',
            'Mass -> T3MassSpec["F", Lookup[f, "MassSymbol", MF]]',
        ),
        'Lagrangian/interactions/ScalarPotential.wl': (
            'BuildHiggsPortalInvariantCandidates',
            '1/2 PlusHc[',
        ),
        'Numerical/direct_one_loop/core.py': (
            'def scalar_real_basis(', 'def evaluate(',
            '1.0/(16.0*math.pi**2)',
        ),
    }
    result = {}
    for name, phrases in required.items():
        path = root / name
        require(path.is_file(), f'Missing source for version guard: {path}')
        code = path.read_text(encoding='utf-8-sig')
        for phrase in phrases:
            require(phrase in code, f'Native source contract changed: {name}: {phrase}')
        result[name] = digest(path)
    # Ensure no novel topology normalization silently changes the operator.
    topology = (root/'Lagrangian/interactions/T3Topology.wl').read_text(encoding='utf-8-sig')
    require(re.search(r'"Name"\s*->\s*"T3ScalarMix_GenericCG".*?"Expression"\s*->\s*PlusHc\[\s*lambdaT3\[\]\s*H\[i\]\s*H\[j\]\s*s1\s*s2\s*CG\[T3MixCG,\s*labels\]', topology, re.S) is not None,
            'Mixing operator source sign/prefactor not recognized')
    return result


def mass_tensor(mf: np.ndarray, rep_dim: int, majorana: bool,
                pairing: np.ndarray) -> np.ndarray:
    """Independent all-left-handed matrix for -psi^T M psi/2 + h.c.

    The Dirac representation uses (F_L,F_R^c), and M has M_D and M_D^T
    off-diagonal. Self-conjugate real SU2 representations use the symmetric
    weight-basis charge-pairing metric. This is an algebraic candidate, not a
    derivation from Matchete's native FreeLag spinor mass term.
    """
    require(mf.shape == (3, 3), 'Expected three-flavour benchmark MF')
    n = 3*rep_dim
    if majorana:
        require(pairing.shape == (rep_dim,rep_dim), 'Bad Majorana pairing shape')
        require(np.linalg.norm(pairing-pairing.T) < EPS, 'Nonsymmetric Majorana charge metric')
        require(np.linalg.norm(pairing.conj().T@pairing-np.eye(rep_dim)) < EPS,
                'Majorana charge metric not unitary')
        mass = np.kron(mf, pairing)
        require(np.linalg.norm(mass-mass.T) < EPS*max(np.linalg.norm(mass),1),
                'MF is not symmetric for a Majorana multiplet')
    else:
        md = np.kron(mf, np.eye(rep_dim))
        mass = np.block([[np.zeros((n,n),complex), md],
                         [md.T, np.zeros((n,n),complex)]])
    return np.asarray(mass,complex)


def scalar_tensor(m1: float, m2: float, h1: float, h2: float, lam: complex,
                  vev: float, mix: np.ndarray) -> np.ndarray:
    """Candidate V2 quadratic from L_mix=+lambda H H S1 S2*+h.c.

    Uses H0=v/sqrt(2), L=-V, and native H0,H0 mixing CG. Because this is
    quadratic in complex scalars with canonical kinetic terms, the mixing
    occupies the off-diagonal Hermitian A (not the symmetric B).
    """
    require(mix.ndim == 2, 'Topology CG VEV slice must be a matrix')
    d1,d2=mix.shape
    a=np.diag([m1-h1*vev*vev/2]*d1+[m2-h2*vev*vev/2]*d2).astype(complex)
    # V_mix = -lambda*v^2/2 * sum_{s,t} CG_{00st} S1_s S2_t^* + h.c.
    a[d1:, :d1] += -(lam*vev*vev/2)*mix.T
    a[:d1, d1:] = a[d1:, :d1].conj().T
    return a


def scalar_quadratic_action(a: np.ndarray, b: np.ndarray, s: np.ndarray) -> complex:
    return np.vdot(s, a@s) + 0.5*(s.T@b@s + (s.T@b@s).conjugate())


def scalar_action_source(s: np.ndarray, d1: int, m1: float, m2: float,
                         h1: float, h2: float, lam: complex,
                         vev: float, mix: np.ndarray) -> float:
    s1,s2=s[:d1],s[d1:]
    diag=(m1-h1*vev*vev/2)*np.vdot(s1,s1)
    diag+=(m2-h2*vev*vev/2)*np.vdot(s2,s2)
    z=lam*vev*vev/2*(s1@mix@s2.conj())
    return float(np.real(diag-z-z.conjugate()))


def scalar_action_witness(a: np.ndarray, b: np.ndarray, d1: int,
                          o: dict, v: float, mix: np.ndarray) -> dict:
    # Two deterministic genuinely complex field configurations, not just a
    # diagonal or real-only probe. These detect index-transpose and HC errors.
    n=len(a)
    s=np.asarray([(i+2)/11+1j*((i*i+3)%7+1)/9 for i in range(n)])
    s2=np.asarray([(-1)**i*(i+1)/7+1j*((3*i+1)%5+2)/13 for i in range(n)])
    vals=[]
    for t in (s,s2):
        lhs=scalar_quadratic_action(a,b,t)
        rhs=scalar_action_source(t,d1,float(o['mS1Sq']),float(o['mS2Sq']),
                                 float(o['lambdaH1']),float(o['lambdaH2']),
                                 _number(o['lambdaT3']),v,mix)
        require(abs(lhs-rhs) < EPS*max(1.0,abs(rhs)),
                'Source quadratic action is not represented by A/B')
        vals.append(float(abs(lhs-rhs)))
    require(np.linalg.norm(a-a.conj().T)<EPS*max(1,np.linalg.norm(a)),
            'Scalar A is not Hermitian')
    require(np.linalg.norm(b-b.T)<EPS*max(1,np.linalg.norm(b)),
            'Scalar B is not symmetric')
    return {'independent_complex_field_witnesses':2,
            'maximum_absolute_source_potential_discrepancy_gev4':max(vals)}


def charge_channel_witness(d1: int, d2: int, df: int, alpha: int,
                           pairing: np.ndarray, c1: np.ndarray,
                           c2: np.ndarray, mix: np.ndarray,
                           majorana: bool) -> dict:
    """Independent explicit Q=T3+Y/2 check on all NONZERO CG channels."""
    q=lambda d,y,j: (d-1)/2-j+y
    qs1=[q(d1,alpha/2,j) for j in range(d1)]
    qs2=[q(d2,(alpha+2)/2,j) for j in range(d2)]
    qf=[q(df,(alpha+1)/2,j) for j in range(df)]
    c=0
    for f in range(df):
        for j in range(d1):
            if abs(c1[0,f,j])>1e-12:
                require(abs(qf[f]-qs1[j])<1e-12,
                        'Yukawa1 native CG violates electric charge')
                c+=1
        for j in range(d2):
            if abs(c2[0,f,j])>1e-12:
                require(abs(qf[f]-qs2[j])<1e-12,
                        'Yukawa2 native CG violates electric charge')
                c+=1
    for i in range(d1):
        for j in range(d2):
            if abs(mix[i,j])>1e-12:
                require(abs(qs1[i]-qs2[j])<1e-12,
                        'Topology scalar quadratic mixes different charges')
                c+=1
    if majorana:
        for f in range(df):
            for j in range(df):
                if abs(pairing[f,j])>1e-12:
                    require(abs(qf[f]+qf[j])<1e-12,
                            'Majorana charge-pair metric couples nonopposite Q')
                    c+=1
    require(c>0,'No nonzero native charge channels')
    return {'native_component_charge_channels_checked':c,
            'majorana_charge_pairing_checked':majorana}


def couplings_candidate(bench: dict, c1: np.ndarray,c2: np.ndarray,
                        pairing: np.ndarray, rep_dim: int,
                        majorana: bool, d1: int, d2: int) -> tuple[np.ndarray,np.ndarray]:
    """Restricted candidate Y/Z from native CG charge-basis coefficients.

    Complex conjugation is that of the EXISTING provisional builder and is
    checked, not certified. Overall source spinor signs and RGBeta/Matchete
    coupling identification are deliberately LEFT UNRESOLVED.
    """
    o=bench['base_state']['ordinary']
    y1,y2=complex_array(o['y1'],'y1'),complex_array(o['y2'],'y2')
    require(y1.shape==y2.shape==(3,3),'Yukawa flavour sizes changed')
    n=3*rep_dim
    nf=n if majorana else 2*n
    field=np.zeros((3,nf,d1+d2),complex)
    conj=np.zeros_like(field)
    for r in range(3):
        for f in range(rep_dim):
            for s in range(d1):
                k=r*rep_dim+f if majorana else n+r*rep_dim+f
                conj[:,k,s] += y1[:,r].conj()*c1[0,f,s].conj()
            for t in range(d2):
                if majorana:
                    for j in range(rep_dim):
                        field[:,r*rep_dim+j,d1+t] += (y2[:,r].conj() *
                            c2[0,f,t].conj()*pairing[f,j])
                else:
                    field[:,r*rep_dim+f,d1+t] += y2[:,r].conj()*c2[0,f,t].conj()
    return conj,field


def candidate_negative_controls(a: np.ndarray, mix: np.ndarray,
                                vev: float, lam: complex, h1: float,
                                h2: float, d1: int) -> dict:
    """Check that wrong offdiagonal signs and missing 1/2 are detectable.

    A zero coefficient is 'not applicable', never a claimed passed control.
    """
    delta=-(lam*vev*vev/2)*mix.T
    against=[]
    if np.linalg.norm(delta) > 1e-16:
        bad=a.copy()
        bad[d1:,:d1] -= 2*delta
        bad[:d1,d1:]=bad[d1:,:d1].conj().T
        against.append(('mixing_sign_flip',bool(np.linalg.norm(bad-a)>1e-16)))
        bad=a.copy()
        bad[d1:,:d1] += delta
        bad[:d1,d1:]=bad[d1:,:d1].conj().T
        against.append(('missing_vev_one_half',bool(np.linalg.norm(bad-a)>1e-16)))
    if max(abs(h1),abs(h2))>1e-16:
        bad=a.copy()
        bad[:d1,:d1] += np.eye(d1)*h1*vev*vev
        bad[d1:,d1:] += np.eye(len(a)-d1)*h2*vev*vev
        against.append(('portal_sign_flip',bool(np.linalg.norm(bad-a)>1e-16)))
    require(all(ok for _,ok in against), 'Negative-control mutation went undetected')
    return {'controls':dict(against), 'nonzero_controls_checked':len(against)}


def audit_case(root: Path, short: str, key: str,
               cases: dict, us: dict, physical: Path, audits: Path) -> dict:
    ref=root/'Numerical/benchmarks/fitted_reference'/f'{key}.json'
    inp=physical/f'{short}.json'
    bench=load_json(ref)
    raw=load_json(inp)
    dim=(*CLASSES[short[0]][:3],int(short[2:]))
    info=bench['representation']
    require(tuple(info[k] for k in ('d_s1','d_s2','d_f','alpha'))==dim,
            f'{short}: benchmark representation differs from model inventory')
    require(raw.get('model')==short,f'{short}: component file model mismatch')
    prov=raw.get('provenance') or {}
    require(prov.get('benchmark_sha256')==digest(ref),
            f'{short}: component input not matched to original benchmark SHA')
    require(prov.get('cg_support_sha256')==digest(audits/'component_support_audit.json')
            and prov.get('generator_analysis_sha256')==digest(audits/'matchete_generator_analysis.json')
            and prov.get('cg_generator_covariance_checked') is True,
            f'{short}: component input native-CG/generator provenance mismatch')
    require(prov.get('validated_uv_component_conventions') is False,
            f'{short}: provisional tensor was unexpectedly promoted to certified')
    case=cases[dim]
    require(tuple(int(case['model'][k]) for k in ('dS1','dS2','dF','alpha'))==dim,
            f'{short}: CG source dimensions disagree')
    d1,d2,df,alpha=dim
    cg={name:_tensor(case['vertices'][name],flags,us)
        for name,flags in FLAGS.items()}
    mf=complex_array(bench['base_state']['ordinary']['MF'],'benchmark MF')
    majorana=(df%2==1 and alpha==-1)
    pairing=np.asarray(us[df].T*us[df],dtype=complex)
    charge=charge_channel_witness(d1,d2,df,alpha,pairing,
                  cg['T3Y1CG'],cg['T3Y2CG'],cg['T3MixCG'][1,1,:,:],majorana)
    expected_m=mass_tensor(mf,df,majorana,pairing)
    input_m=complex_array(raw['fermion_mass_gev'],'input MF')
    mres=near(input_m,expected_m,f'{short}: double-Weyl mass candidate')
    o=bench['base_state']['ordinary']
    vev=float(bench['numerical']['vev_gev'])
    require(vev>0, f'{short}: nonpositive VEV')
    mixing=cg['T3MixCG'][1,1,:,:]
    expect_a=scalar_tensor(float(o['mS1Sq']),float(o['mS2Sq']),
                           float(o['lambdaH1']),float(o['lambdaH2']),
                           _number(o['lambdaT3']),vev,mixing)
    inp_a=complex_array(raw['scalar_A_gev2'],'input scalar A')
    inp_b=complex_array(raw['scalar_B_gev2'],'input scalar B')
    ares=near(inp_a,expect_a,f'{short}: scalar A candidate')
    bres=near(inp_b,np.zeros_like(expect_a),f'{short}: restricted scalar B')
    source_witness=scalar_action_witness(inp_a,inp_b,d1,o,vev,mixing)
    minus=candidate_negative_controls(inp_a,mixing,vev,_number(o['lambdaT3']),
                                    float(o['lambdaH1']),float(o['lambdaH2']),d1)
    cy1,cy2=couplings_candidate(bench,cg['T3Y1CG'],cg['T3Y2CG'],
                               pairing,df,majorana,d1,d2)
    verts=raw.get('vertices')
    require(isinstance(verts,list) and len(verts)==2 and
            [p.get('name') for p in verts]==['y1_HC','y2_HC'],
            f'{short}: expected precisely two restricted Yukawa vertices')
    v1=complex_array(verts[0]['conjugate'],'y1_HC')
    v2=complex_array(verts[1]['field'],'y2_HC')
    y1res=near(v1,cy1,f'{short}: conditional y1 spinor map')
    y2res=near(v2,cy2,f'{short}: conditional y2 spinor map')
    require(np.linalg.norm(v1)>0 and np.linalg.norm(v2)>0,
            f'{short}: one Yukawa operator vanishes')
    # Only existing candidate orientations are compared here. A missing
    # native concrete C, epsilon, or benchmark coupling phase is never
    # inferred from cancellation or these numerical equalities.
    return {
        'model':short,'source_benchmark':relative(ref,root),
        'benchmark_sha256':digest(ref),
        'input_file':relative(inp,root),'input_sha256':digest(inp),
        'fermion_type':'Majorana' if majorana else 'Dirac',
        'electric_charge_channel_witness':charge,
        'fermion_double_Weyl_candidate_residual':mres,
        'fermion_charge_pair_metric_symmetric':bool(np.linalg.norm(pairing-pairing.T)<EPS),
        'scalar_A_candidate_residual':ares, 'scalar_B_candidate_residual':bres,
        'scalar_source_action':source_witness,'negative_controls':minus,
        'conditional_yukawa1_candidate_residual':y1res,
        'conditional_yukawa2_candidate_residual':y2res,
        'source_group_charge_basis_is_candidate_generator_not_live_export':True,
        'native_spinor_phase_certified':False,
        'native_fermion_mass_normalization_certified':False,
        'absolute_portal_and_topology_coupling_dictionary_certified':False,
        'physical_uv_mass_certified':False,
        'status':CONDITIONAL,
    }


def _previous_evidence(root: Path, audits: Path) -> dict:
    cg=load_json(audits/'component_support_audit.json')
    analysis=load_json(audits/'matchete_generator_analysis.json')
    cov=load_json(audits/'cg_covariance_audit.json')
    require(not cg.get('errors') and len(cg.get('cases',[]))==16,
            'Native CG component support report incomplete')
    require(cov.get('status')=='cg_generator_covariance_verified'
            and len(cov.get('models',[]))==16
            and all(r.get('all_covariant') is True for r in cov['models']),
            'All-16 SU2 covariance evidence missing')
    rot=_rotations(analysis)
    ids={tuple(int(v['model'][k]) for k in ('dS1','dS2','dF','alpha')):v
         for v in cg['cases']}
    require(len(ids)==16,'16 distinct native CG model cases required')
    portal=load_json(audits/'native_portal_identity_decomposition.json')
    require(portal.get('models_checked')==16 and
            portal.get('scalar_portal_sectors_checked')==32 and
            portal.get('full_tensor_identities_representable')==32 and
            len(portal.get('models',[]))==16,
            'Prior 32-sector portal identity proof missing')
    require(portal.get('status')==
            'all_32_geometric_identities_representable_physical_dictionary_pending',
            'Portal identity certificate status changed')
    return {'cases':ids,'rotations':rot,'source_evidence_hashes':{
        relative(audits/p,root):digest(audits/p) for p in (
            'component_support_audit.json', 'matchete_generator_analysis.json',
            'cg_covariance_audit.json', 'native_portal_identity_decomposition.json')}}


def prior_spinor_evidence(root: Path) -> dict:
    """Reuse A(0) exact/natively checked spinor phase witness, if present.

    Missing prior evidence is explicitly reported as a blocker, and a previous
    conditional witness is never counted as a full 16-model source mapping.
    """
    rel=Path('Reports/output/workstream_a/yukawa_spinor_phase/a0_spinor_phase.json')
    path=root/rel
    if not path.is_file():
        return {'present':False,'reason':'A(0) native spinor-phase witness not found',
                'cross_package_physical_dictionary_certified':False}
    obj=load_json(path)
    require(obj.get('status')==
            'NativeChargeConjugationCheckedCanonicalWeylPhaseFixedByChoicePhysicalDictionaryPending'
            and obj.get('model')=='A_p0'
            and obj.get('matchete_native_abstract_cc_algebra_verified') is True
            and obj.get('matchete_installed_concrete_CC_orientation_proved') is False
            and obj.get('rgbeta_installed_numeric_epsilon_orientation_proved') is False
            and obj.get('global_dirac_mass_and_y2_phase_transport_proved') is False
            and obj.get('physical_cross_package_yukawa_dictionary_certified') is False,
            'Prior spinor-phase report incompatible with its documented limited scope')
    return {'present':True, 'path':str(rel).replace('\\','/'),
            'sha256':digest(path),
            'A0_native_abstract_CC_and_chosen_basis_weyl_identity_passed':True,
            'global_dirac_mass_and_y2_phase_transport_proved':False,
            'cross_package_physical_dictionary_certified':False}


def audit_all(root: Path, physical: Path, audits: Path) -> dict:
    """Never fail-open; keep per-model errors and no physically certified flag."""
    contracts=check_source_contracts(root)
    saved=_previous_evidence(root,audits)
    spinor=prior_spinor_evidence(root)
    results=[]
    for short,key in models():
        try:
            row=audit_case(root,short,key,saved['cases'],saved['rotations'],physical,audits)
        except (ValueError,KeyError,IndexError,TypeError,OSError) as exc:
            row={'model':short,'status':'BLOCKED','reason':str(exc),
                 'physical_uv_mass_certified':False}
        results.append(row)
    passed=sum(x['status']==CONDITIONAL for x in results)
    require(len(results)==16,'Model inventory changed unexpectedly')
    return {
        'format':'t3-direct-uv-integrated-conventions-v1',
        'status':'physical_uv_normalization_NOT_CERTIFIED',
        'scope':'restricted frozen benchmark, 16 direct UV models; no EFT/RGE',
        'algebraic_candidate_models_passed':passed,
        'algebraic_candidate_models_blocked':16-passed,
        'physical_models_certified':0,
        'source_guard_sha256':contracts,
        'existing_evidence_sha256':saved['source_evidence_hashes'],
        'prior_A0_spinor_phase_evidence':spinor,
        'checks_covered':[
            'both native topology Yukawa field orientations and prior source CG rotations',
            'all nonzero component Q-conservation channels and Majorana Q pairing',
            'candidate Dirac two-Weyl or Majorana charge-metric fermion mass tensor',
            'H0=v/sqrt(2) topology mixing in Hermitian scalar A with -lambda*v^2/2',
            'portal mass shifts and restricted B=0 quadratic action',
            'two independent complex-field quadratic-action witnesses',
            'wrong topology mixing sign/half and portal-sign negative controls',
            'candidate complex-conjugated Yukawa tensors and flavour indexing',
        ],
        'limitations_not_silently_certified':[
            'Matchete abstract GammaCC has no native numerical 4-spinor sign anchor',
            'candidate SU2 generator basis was validated by covariance, not exported live',
            'RGBeta fundamental epsilon orientation and benchmark y1/y2 dictionary pending',
            'source FreeLag Weyl mass/charge conjugation and Majorana 1/2 not independently projected',
            'the full native portal invariant combination and physical benchmark coupling dictionary pending',
            'extra symmetry-allowed interactions absent from restricted fitted benchmarks',
        ],
        'never_modify_frozen_benchmarks_or_provisional_results':True,
        'input_provenance_never_promoted_to_validated':True,
        'models':results,
    }


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=ROOT)
    ap.add_argument('--physical-input',type=Path,default=None)
    ap.add_argument('--audit-dir',type=Path,default=None)
    ap.add_argument('--output-dir',type=Path,default=None)
    args=ap.parse_args()
    root=args.project_root.resolve()
    physical=(args.physical_input or root/'input/lo_direct_one_loop').resolve()
    audits=(args.audit_dir or root/CG_DIR).resolve()
    out=(args.output_dir or root/OUT).resolve()
    try:
        report=audit_all(root,physical,audits)
    except (ValueError,KeyError,TypeError,IndexError,OSError) as exc:
        report={'format':'t3-direct-uv-integrated-conventions-v1',
                'status':'BLOCKED_missing_or_changed_source_evidence',
                'reason':str(exc),'physical_models_certified':0,
                'source_component_inputs_not_modified':True}
    out.mkdir(parents=True,exist_ok=True)
    path=out/'integrated_convention_audit.json'
    path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(f"Batch1: {report['status']}")
    print('Physically certified: 0/16')
    if 'algebraic_candidate_models_passed' in report:
        print('Candidate source algebra checks:',
              f"{report['algebraic_candidate_models_passed']}/16")
        for rec in report['models']:
            print(f" {rec['model']:6s}: {rec['status']}",rec.get('reason',''))
    elif report.get('reason'):
        print('BLOCKED:',report['reason'])
    print('Report:',path)
    # A source-algebra witness alone is never an admissible physical certificate.
    return 2


if __name__=='__main__':
    raise SystemExit(main())
