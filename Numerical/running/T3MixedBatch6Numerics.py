"""Restricted, source-driven B/C mixed scalar hard C5 lift to three heavy masses.

Only a native scalar Match/ExtractWeinbergCoefficient export may be used. This
module deliberately does NOT infer a hard coefficient from the full-theory loop,
ordinary-model SU(2) prefactors, or a synthetic scalar tensor. A pole cancellation
certificate needs an independent first-threshold inherited source, which is not
produced by the archived mixed F-first TREE and remains blocked here.
"""
from __future__ import annotations

import ast
from dataclasses import dataclass
from hashlib import sha256
import math
from pathlib import Path
import re
from typing import Any

import mpmath as mp
import numpy as np

# Standalone identical I/O and Majorana mass conventions to Batch 5; no
# dependence on a previously generated package version for these primitives.
def typed_complex_matrix(raw: Any, label: str) -> np.ndarray:
    if not isinstance(raw, list) or len(raw) != 3:
        raise ValueError(label+': 3x3 complex-list matrix expected')
    a = np.empty((3, 3), complex)
    for i, row in enumerate(raw):
        if not isinstance(row, list) or len(row) != 3:
            raise ValueError(label+': invalid complex matrix row')
        for j, z in enumerate(row):
            if not isinstance(z, dict) or not {'real','imag'} <= set(z):
                raise ValueError(label+': complex entry requires real/imag')
            a[i,j] = complex(float(z['real']), float(z['imag']))
    if not np.isfinite(a).all():
        raise ValueError(label+': nonfinite complex matrix')
    return a


def scalar_complex(z: Any, label: str) -> complex:
    if isinstance(z, dict):
        z = complex(float(z['real']), float(z['imag']))
    else:
        z = complex(z)
    if not (math.isfinite(z.real) and math.isfinite(z.imag)):
        raise ValueError(label+': nonfinite number')
    return z


def encode_complex(a: np.ndarray) -> dict:
    a = np.asarray(a, dtype=complex)
    if a.shape != (3,3) or not np.isfinite(a).all():
        raise ValueError('Expected finite 3x3 symmetric C5 matrix')
    return {'real':a.real.tolist(), 'imag':a.imag.tolist()}


def metrics(a: np.ndarray, vev_gev: float) -> dict:
    a = np.asarray(a, dtype=complex)
    if a.shape != (3,3) or not np.isfinite(a).all():
        raise ValueError('Bad C5 matrix')
    err = np.linalg.norm(a-a.T)/max(1e-40, np.linalg.norm(a))
    if err > 1e-11:
        raise ValueError('Majorana C5 not symmetric')
    return {'conditional_takagi_masses_ev_ascending':
        np.sort(np.linalg.svd(a,compute_uv=False)*vev_gev**2/2*1e9).tolist(),
        'c5_symmetry_relative_residual':float(err)}

ALLOWED = {'MF','MSX','MSY','T3MEU11','T3MEU12','T3MEU21','T3MEU22',
           'lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1',
           'y1a','y1b','y2a','y2b','hbar','EPS','MU2'}
BASIC = {'MF','MSX','MSY','T3MEU11','T3MEU12','T3MEU21','T3MEU22','lambdaT3',
         'y1a','y1b','y2a','y2b','hbar','EPS','MU2'}
QUARTICS = {'lambdaT3ExtraQ020020I1':'lambdaHHdagS2S2',
            'lambdaT3ExtraQ020200I1':'lambdaHHdagS1barS1bar'}
LOOP = 1/(16*math.pi**2)


class NativeSourceUnsupported(ValueError):
    pass


def require(ok: bool, message: str) -> None:
    if not ok:
        raise NativeSourceUnsupported(message)


def safe_compiler(raw: str) -> tuple[Any,set[str]]:
    """Parse only the source's allowed one-generation Matchete scalar kernel.

    An unknown native CG, field, spinor or quartic is a BLOCKER; it must not be
    ignored or silently mapped to the ordinary-model kernel.
    """
    require(isinstance(raw,str) and 100<len(raw)<2_000_000,'Native hard C5 missing/oversized')
    src=raw.strip()
    # Only a simple epsilon pole is supported: the finite/pole separation
    # E(EPS=+1) +/- E(EPS=-1) is exact within this verified source language.
    require('ϵ' in src and 'hbar' in src,'Native source has no explicit pole/hbar')
    require('ϵ^' not in src,'Higher epsilon-pole order not supported')
    for name in ('y1','y2'):
        for idx,label in (('d$$1','a'),('d$$2','b')):
            expr=f'Coupling[{name}, {{Index[{idx}, Flavor], Index[d$$1, NFlavor]}}, 0]'
            src=src.replace(f'Bar[{expr}]',f'conj({name}{label})')
            src=src.replace(expr,name+label)
    present=set(re.findall(r'Coupling\[([A-Za-z][A-Za-z0-9]*),',src))
    unknown=present-(ALLOWED-{'y1a','y1b','y2a','y2b','EPS','MU2','hbar'})
    require(not unknown,'Unsupported native coupling inventory: '+str(sorted(unknown)))
    for name in sorted(present,key=len,reverse=True):
        sym=f'Coupling[{name}, {{}}, 0]'
        src=src.replace(f'Bar[{sym}]',f'conj({name})').replace(sym,name)
    src=src.replace('μbar2','MU2').replace('ϵ','EPS')
    for head in ('Log','Sqrt'):
        src=src.replace(head+'[',head.lower()+'(')
    src=src.replace(']',')').replace('^','**')
    require(not any(x in src for x in ('Coupling','Index','Field','Gamma','CG[','NCM','Bar[','d$$')),
            'Native hard result contains unprojected fields/CG/flavour indices')
    require('[' not in src and ']' not in src,'Unsupported Mathematica head in native result')
    require(not re.search(r'[^A-Za-z0-9_.,+\-*/()\s]',src),
            'Unrecognized symbol/function in native expression')
    try: tree=ast.parse(src,mode='eval')
    except (SyntaxError,RecursionError) as exc:
        raise NativeSourceUnsupported('Native hard C5 not parseable') from exc
    for n in ast.walk(tree):
        require(isinstance(n,(ast.Expression,ast.Load,ast.Add,ast.Sub,ast.Mult,ast.Div,
                 ast.Pow,ast.UAdd,ast.USub,ast.BinOp,ast.UnaryOp,ast.Constant,
                 ast.Call,ast.Name)),f'Unsafe expression node: {type(n).__name__}')
        if isinstance(n,ast.Constant):
            require(type(n.value) in (int,float),'Unexpected source literal')
        if isinstance(n,ast.Name):
            require(n.id in ALLOWED|{'conj','log','sqrt'},'Unexpected native symbol: '+n.id)
        if isinstance(n,ast.Call):
            require(isinstance(n.func,ast.Name) and n.func.id in ('conj','log','sqrt')
                    and len(n.args)==1 and not n.keywords,'Unsupported function call')
    return compile(tree,'<checked-native-mixed-C5>',mode='eval'),present


def _finite_complex(x: Any, label: str) -> complex:
    z=scalar_complex(x,label)
    require(np.isfinite(z.real) and np.isfinite(z.imag),label+' nonfinite')
    return z


def cval(obj: Any) -> complex:
    if isinstance(obj,dict):
        return complex(float(obj['real']),float(obj['imag']))
    return complex(obj)


def _matrix(raw: Any) -> np.ndarray:
    a=np.asarray(raw['real'],float)+1j*np.asarray(raw['imag'],float)
    require(a.shape==(2,2) and np.isfinite(a).all(),'Invalid mixed-doublet U')
    return a


def numerical_inputs(model: str, running: dict, params: dict,
                     *, conditional_quartic_map: bool) -> dict:
    require(model in ('B_-1','C_-1'),'Only mixed B/C(-1)')
    require(running.get('model')==model and
            running.get('status')=='NumericalRenormalisableRunningOnly',
            'Unexpected renormalisable running report')
    uv=running['uv_at_fermion_threshold']; eft=running['eft1_at_fermion_threshold']
    y1=typed_complex_matrix(uv['y1'],'y1(muF)')
    y2=typed_complex_matrix(uv['y2'],'y2(muF)')
    fmat=typed_complex_matrix(uv['MF'],'MF(muF)')
    masses=np.diag(fmat)
    require(np.linalg.norm(fmat-np.diag(masses))<1e-12*max(np.linalg.norm(fmat),1.)
            and np.max(abs(masses.imag))<1e-12*max(np.linalg.norm(fmat),1.)
            and np.min(masses.real)>0,'MF must be diagonal, positive, real: no invented Takagi map')
    require(np.max(masses.real)/np.min(masses.real)>1.00001,'Need three nondegenerate MF')
    U=_matrix(params['unitary_U'])
    masses_scalar=np.asarray(params['heavy_mass_squared_gev2'],float)
    require(np.all(masses_scalar>0) and masses_scalar.shape==(2,) and
            np.allclose(U.conj().T@U,np.eye(2),atol=1e-10),
            'Invalid native X/Y mass eigenbasis')
    lam=_finite_complex(eft['lambdaT3'],'lambdaT3')
    muF=float(running['scales_gev']['F'])
    require(muF>max(np.sqrt(masses_scalar)),'F threshold must exceed both X and Y masses')
    couplings={'lambdaT3':lam}
    if conditional_quartic_map:
        couplings.update({m:-_finite_complex(eft[r],r) for m,r in QUARTICS.items()})
    return {'y1':y1,'y2':y2,'mf':masses.real,'U':U,
            'ms':np.sqrt(masses_scalar),'couplings':couplings,'muF':muF,
            'quartic_map_certified':False}


def _matrix_at_scale(code: Any, names: set[str], inputs: dict,
                     mu: float, eps: int) -> np.ndarray:
    U=inputs['U']; ms=inputs['ms'];
    env={f'T3MEU{i+1}{j+1}':mp.mpc(complex(U[i,j]))
         for i in range(2) for j in range(2)}
    env.update({k:mp.mpc(v) for k,v in inputs['couplings'].items()})
    env.update(MSX=mp.mpf(str(ms[0])),MSY=mp.mpf(str(ms[1])),
               MU2=mp.mpf(str(mu))**2, EPS=mp.mpf(eps),
               hbar=mp.mpf(1)/(16*mp.pi**2),conj=mp.conj,log=mp.log,sqrt=mp.sqrt)
    for name in names:
        require(name in env or name in ('y1','y2','MF'),
                f'Numerical evaluation needs independent source coupling {name}; no automatic guessing')
    result=np.zeros((3,3),complex)
    for i in range(3):
        for j in range(3):
            z=mp.mpc(0)
            for r,m in enumerate(inputs['mf']):
                local=dict(env,MF=mp.mpf(str(m)))
                for f in ('y1','y2'):
                    local[f+'a']=mp.mpc(complex(inputs[f][i,r]))
                    local[f+'b']=mp.mpc(complex(inputs[f][j,r]))
                val=eval(code,{'__builtins__':{}},local)
                require(mp.isfinite(val.real) and mp.isfinite(val.imag),'Nonfinite native hard C5')
                z+=val
            result[i,j]=complex(z)
    return (result+result.T)/2  # symmetric Majorana C5, NOT Hermitian projection


def evaluate_hard(model: str, source_file: Path, running: dict,
                  params: dict, *, conditional_quartic_map: bool,
                  muS: float|None=None, vev_gev:float=246.22) -> dict:
    raw=source_file.read_text(encoding='utf-8-sig')
    code,names=safe_compiler(raw)
    inp=numerical_inputs(model,running,params,
                         conditional_quartic_map=conditional_quartic_map)
    missing=(names-set(inp['couplings']))-{'MF','MSX','MSY','T3MEU11','T3MEU12',
      'T3MEU21','T3MEU22'}
    require(not missing,'Native coupling map not provided/certified: '+str(sorted(missing)))
    matching_scale=float(muS or max(inp['ms']))
    require(0<matching_scale<inp['muF'] and math.isfinite(matching_scale),
            'Invalid scalar matching scale')
    require(math.isfinite(vev_gev) and vev_gev>0,'Invalid Higgs vev')
    with mp.workdps(45):
        plus=_matrix_at_scale(code,names,inp,matching_scale,1)
        minus=_matrix_at_scale(code,names,inp,matching_scale,-1)
        finite=(plus+minus)/2
        pole=(plus-minus)/2
        # source has at most a simple epsilon pole; test E(EPS=2) consistency
        twice=_matrix_at_scale(code,names,inp,matching_scale,2)
        predicted=finite+pole/2
        scale=max(np.linalg.norm(finite),np.linalg.norm(pole),1e-40)
        eps_error=float(np.linalg.norm(twice-predicted)/scale)
        require(eps_error<2e-10,'Non-simple epsilon structure in native hard C5')
        shifted_mu=matching_scale*math.exp(0.037)
        plus2=_matrix_at_scale(code,names,inp,shifted_mu,1)
        minus2=_matrix_at_scale(code,names,inp,shifted_mu,-1)
        finite2=(plus2+minus2)/2
        pole2=(plus2-minus2)/2
    # This is an explicit mass/log closure of the scalar-only hard expression.
    # It does not establish the independent EFT1 anomalous dimension or CT.
    log_residual=float(np.linalg.norm((finite2-finite)/(2*.037)-pole)/scale)
    pole_invariance=float(np.linalg.norm(pole-pole2)/scale)
    require(pole_invariance<5e-10,'Hard pole unexpectedly depends on scalar matching scale')
    # Do not subtract either the archived one-scale baseline or raw F pole.
    return {'format':'t3-mixed-batch6-source-native-hard-three-flavour-v1',
        'model':model,'status':'NativeHardThreeFlavorPoleSeparated_NotInheritedPoleCancelled',
        'eft_order':'UV -> ALL F -> EFT1(SM+X+Y) -> ALL X,Y -> SMEFT',
        'hard_scalar_source_sha256':sha256(source_file.read_bytes()).hexdigest(),
        'source_couplings':sorted(names),'fermion_mass_eigenvalues_gev':inp['mf'].tolist(),
        'scalar_mass_eigenvalues_gev':inp['ms'].tolist(),
        'scalar_U':{'real':inp['U'].real.tolist(),'imag':inp['U'].imag.tolist()},
        'common_couplings_at':'muF','muF_gev':inp['muF'],'muS_gev':matching_scale,
        'hard_finite_c5_gev_inverse':encode_complex(finite),
        'hard_pole_residue_c5_gev_inverse':encode_complex(pole),
        'hard_finite_only_conditional_takagi_masses_ev':metrics(finite,vev_gev)['conditional_takagi_masses_ev_ascending'],
        'checks':{'three_nondegenerate_fermions_summed':True,
           'two_mixed_scalar_mass_eigenstates_used':True,
           'native_one_eigenvalue_template_evaluated':True,
           'three_distinct_MF_lift_inferred_from_common_MF_symbol':True,
           'independent_three_MF_native_matching_performed':False,
           'simple_pole_decomposition_relative_error':eps_error,
           'pole_matching_scale_invariance_error':pole_invariance,
           'native_hard_log_vs_pole_relative_residual':log_residual,
           'direct_EFT1_operator_rge_verified':False,
           'first_threshold_inherited_pole_cancelled':False,
           'matched_MSbar_counterterms_verified':False,
           'full_fixed_order_c5_assembled':False},
        'physical_C5_certified':False,'quartic_map_is_conditional':conditional_quartic_map,
        'source_warning':'Saved source uses one unindexed common MF; numerically evaluating it at three eigenmasses is a provisional per-mass lift. Masses use hard finite piece ONLY, not complete physical neutrino masses',
        'pole_warning':'A hard pole is not a pole-cancellation check. Independent inherited/soft CT unavailable.'}
