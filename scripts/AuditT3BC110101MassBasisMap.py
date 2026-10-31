"""Source-native B/C(-1) two-contraction 110101 tensor map in mass basis.

Uses live Wolfram RGBeta SU(2) covariance evidence, saved Matchete native
components and the existing B/C scalar basis. Resolves the *tensor geometry*
only. No RGBeta potential prefactor, Lagrangian sign, Yukawa map, or Matchete
one-loop physical matching is certified.

python -u -X utf8 -m scripts.AuditT3BC110101MassBasisMap --all --run-id bc_110101_map_v1
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from pathlib import Path
import json
import re

import numpy as np
import sympy as sp

from scripts.AuditT3BCNativeQuarticMassBasis import (
    native_tensors, load_basis, rotate_tensor,
)

ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT / 'output/general_t3_stage22/matchete_tensors.json'
RGB_NATIVE = ROOT / 'output/general_t3_stage22/rgbeta_native_invariants.json'
RGB_LIVE = ROOT / 'output/general_t3_stage22/stage22_110101_wolfram_native_covariance.json'
MASS_BASIS = ROOT / 'output/general_t3_mass_basis'
OUTPUT = ROOT / 'output/audits/t3_bc_110101_mass_basis_map'
MODELS = ('B_-1', 'C_-1')
MAT = ('lambdaT3ExtraQ110101I1', 'lambdaT3ExtraQ110101I2')
RGB = ('lambdaT3GeneralQ110101I1', 'lambdaHHdagS1barS2barCross')
SPECIES = ('H', 'Hdag', 'S1dag', 'S2dag')
EXPR = (
    'del[SU2L[fund],t3i1,t3i4]*eps[SU2L[fund],t3i2,t3i3]',
    'del[SU2L[fund],t3i1,t3i3]*eps[SU2L[fund],t3i2,t3i4]',
)


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise ValueError('B/C live RGBeta-Matchete 110101 mass-basis map: ' + msg)


def load(path: Path) -> dict:
    require(path.is_file(), 'missing source ' + str(path))
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data, dict), 'JSON object required: '+str(path))
    return data


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def exact_number(raw: str) -> sp.Expr:
    require(isinstance(raw,str) and re.fullmatch(r'[0-9+\-*/^()\[\]Sqrt\s]+',raw) is not None,
            'not a restricted native exact component: '+repr(raw))
    s=raw.replace('Sqrt[','sqrt(').replace(']',')').replace('^','**')
    x=sp.sympify(s,locals={'sqrt':sp.sqrt})
    require(not x.free_symbols and x.is_finite is True,
            'not an exact finite numeric component')
    return sp.simplify(x)


def source_tensor(i: int) -> np.ndarray:
    """Independently evaluate the two exact live del-epsilon expressions.

    Component convention: delta_ij, eps_01=+1, eps_10=-1.
    The expressions and field order must also match the live Wolfram export.
    """
    delta=np.eye(2,dtype=int)
    eps=np.array([[0,1],[-1,0]],dtype=int)
    t=np.zeros((2,)*4,dtype=int)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    t[a,b,c,d]=(delta[a,d]*eps[b,c] if i==0
                                  else delta[a,c]*eps[b,d])
    return t


def check_live(live: dict, invariants: dict, model: str) -> tuple[np.ndarray,np.ndarray]:
    require(live.get('format')=='stage22-live-rgbeta-su2-covariance-v1' and
            live.get('status')=='NativeWolframSourceCovarianceDiagnosticNotPhysicalCouplingCertificate' and
            live.get('potential_prefactors_verified') is False and
            live.get('physical_coupling_mapping_certified') is False and
            set(live.get('models',{}))==set(MODELS),
            'live native RGBeta covariance evidence has changed or claims physical certification')
    require(invariants.get('format')=='stage22-native-rgbeta-invariants-v1' and
            invariants.get('source')=='RGBetaRegisteredQuarticInvariantFunctions' and
            invariants.get('native_potential_prefactors_verified') is False and
            invariants.get('native_invariant_component_alignment_verified') is False,
            'saved RGBeta registered-invariant provenance or certification changed')
    row=live['models'][model]
    require(row.get('native_registered_symbols')==list(RGB),
            model+': unexpected native RGBeta couplings')
    observed=row.get('native_expression_inputform')
    require(isinstance(observed,list) and len(observed)==2 and
            [''.join(x.split()) for x in observed]==list(EXPR),
            model+': live registered invariant expression changed')
    nat=invariants.get('models',{}).get(model,{}).get('invariants',[])
    require(isinstance(nat,list),'missing native registered RGBeta invariants')
    pair={r['symbol']:r for r in nat if r.get('symbol') in RGB}
    require(set(pair)==set(RGB) and len(pair)==2,model+': missing or duplicate native entries')
    result=[]
    for i,name in enumerate(RGB):
        ref=pair[name]
        require(ref.get('native_fields_inputform')=='{H, Bar[H], Bar[S1], Bar[S2]}' and
                ''.join(ref.get('native_invariant_expression_inputform','').split())==EXPR[i] and
                ''.join(ref.get('native_coupling_function_inputform','').split())==
                    f'6*Matrix[{name}][##1]&' and
                ref.get('component_tensor_exported') is False and
                ref.get('verified_potential_prefactor') is None,
                model+'/'+name+': changed native expression/potential metadata')
        audit=row.get('native_covariance_audits',[])
        require(len(audit)==2 and audit[i].get('su2_covariant') is True and
                audit[i].get('nonzero_residual_components_by_generator')==[0,0,0],
                model+'/'+name+': native SU2 covariance failed')
        tensor=source_tensor(i)
        require(audit[i].get('component_tensor_flat_inputform')==
                [str(int(z)) for z in tensor.ravel()],
                model+'/'+name+': exported live RGBeta components differ from del-epsilon identity')
        result.append(tensor.astype(complex))
    return result[0],result[1]


def original_matchete_records(native: dict, model: str) -> tuple[list[str],list[str]]:
    """The exact unrounded InputForm values from native Matchete export."""
    records=native.get('models',{}).get(model,{}).get('tensors',[])
    require(isinstance(records,list),'missing Matchete native tensor records')
    rec={r['symbol']:r for r in records if r.get('symbol') in MAT}
    require(set(rec)==set(MAT) and len(rec)==2,'missing Matchete 110101 doublet pair')
    out=[]
    for name in MAT:
        row=rec[name]
        require(tuple(row.get('species',[]))==SPECIES and
                row.get('dimensions')==[2]*4 and
                row.get('potential_prefactor')=='1' and
                len(row.get('components',[]))==16,
                model+'/'+name+': invalid native Matchete tensor contract')
        for value in row['components']:
            exact_number(value)
        out.append(list(row['components']))
    return out[0],out[1]


def exact_map(native:dict, model:str, rgb_tensors:tuple[np.ndarray,np.ndarray]) -> sp.Matrix:
    """Find A such that R_i = sum_j M_j A[j,i] with unrounded native inputs.

    Exact SymPy expressions are parsed DIRECTLY from saved Matchete InputForm
    strings. No numerical inference, nsimplify, or tolerance establishes the
    tensor-basis map. Numeric rotation tests happen only after this proof.
    """
    records=original_matchete_records(native,model)
    M=sp.Matrix.hstack(*[sp.Matrix([exact_number(x) for x in row])
                           for row in records])
    require(all(np.max(np.abs(np.imag(t))) < 1e-14 and
                np.max(np.abs(np.real(t)-np.rint(np.real(t)))) < 1e-14
                for t in rgb_tensors),
            'native RGBeta components must be exact integers')
    R=sp.Matrix.hstack(*[sp.Matrix([sp.Integer(int(np.real(x))) for x in np.ravel(t)])
                           for t in rgb_tensors])
    require(M.rank()==2 and R.rank()==2,'expected two independent tensors in each basis')
    A=(M.T*M).inv()*M.T*R
    require(all(sp.simplify(z)==0 for z in M*A-R),
            'RGBeta native invariant not exactly in Matchete tensor span')
    require(sp.simplify(A.det())!=0,'noninvertible native contraction map')
    return A.applyfunc(sp.simplify)


def original_matchete_exact(native: dict, model: str) -> tuple[np.ndarray,np.ndarray]:
    records=original_matchete_records(native,model)
    return tuple(np.asarray([complex(sp.N(exact_number(x),35)) for x in row],
                            dtype=complex).reshape((2,)*4) for row in records)


def compare_rotations(mat:tuple[np.ndarray,np.ndarray], rgb:tuple[np.ndarray,np.ndarray],
                      A:sp.Matrix, U:np.ndarray, trials:int=12) -> dict:
    """Independently rotate M and R; compare every transformed component."""
    mm=[rotate_tensor(t,SPECIES,U) for t in mat]
    rr=[rotate_tensor(t,SPECIES,U) for t in rgb]
    require(all(set(t)==set(mm[0]) for t in mm+rr),
            'native rotations produce inconsistent mass-eigenstate field terms')
    max_component=0.
    numeric=np.asarray([[complex(sp.N(A[i,j],35)) for j in range(2)] for i in range(2)])
    for j in range(2):
        for fields,actual in rr[j].items():
            expected=sum(numeric[i,j]*mm[i][fields] for i in range(2))
            error=np.max(abs(actual-expected))/max(1.,np.max(abs(actual)),np.max(abs(expected)))
            max_component=max(max_component,float(error))
    require(max_component < 2e-12,'rotated native RGBeta and Matchete tensors fail exact geometry')
    rng=np.random.default_rng(2201)
    max_polynomial=0.
    for _ in range(trials):
        random=lambda: rng.normal(size=2)+1j*rng.normal(size=2)
        h,x,y=random(),random(),random()
        s1=U[0,0]*x+U[0,1]*y
        tilde=U[1,0]*x+U[1,1]*y
        eps=np.array([[0,1],[-1,0]])
        s2=-eps@tilde.conjugate()
        old={'H':h,'Hdag':h.conjugate(),'S1dag':s1.conjugate(),'S2dag':s2.conjugate()}
        new={'H':h,'Hdag':h.conjugate(),'X':x,'Xdag':x.conjugate(),
             'Y':y,'Ydag':y.conjugate()}
        for j in range(2):
            pre=np.einsum('abcd,a,b,c,d->',rgb[j],*(old[f] for f in SPECIES))
            post=sum(np.einsum('abcd,a,b,c,d->',t,*(new[f] for f in fields))
                     for fields,t in rr[j].items())
            err=abs(pre-post)/max(1.,abs(pre),abs(post))
            max_polynomial=max(max_polynomial,float(err))
    require(max_polynomial<2e-12,'native RGBeta polynomial not conserved by mass rotation')
    return {'mass_basis_rotated_component_max_relative_error':max_component,
            'rgbeta_original_vs_rotated_polynomial_max_relative_error':max_polynomial,
            'mass_basis_monomial_terms_per_rgbeta_invariant':len(rr[0]),
            'exact_mass_basis_tensor_map_verified':True}


def run_one(model: str, *, mat_source:dict, live:dict, inv:dict,
            mass_root:Path=MASS_BASIS) -> dict:
    require(model in MODELS,'unsupported model')
    native=native_tensors(mat_source,model)
    require(set(MAT)<=set(native),'Matchete source missing independent 110101 basis')
    mat=original_matchete_exact(mat_source,model)
    rgb=check_live(live,inv,model)
    A=exact_map(mat_source,model,rgb)
    expected=sp.Matrix([[sp.Rational(1,2),sp.Integer(1)],[-sp.sqrt(3)/2,sp.Integer(0)]])
    require(all(sp.simplify(z)==0 for z in A-expected),
            model+': source-native exact 110101 map changed; inspect field convention')
    tag=model.replace('-','m')
    basis=load_basis(mass_root/f'{tag}_final_mass_basis.json',model)
    result=compare_rotations(mat,rgb,A,basis.U)
    return {'model':model,'status':'Native110101TensorMapMassBasisChecked',
            'native_matchete_symbols':list(MAT),
            'native_rgbeta_symbols':list(RGB),
            'rgbeta_in_matchete_exact_columns':[[str(A[i,j]) for j in range(2)] for i in range(2)],
            'basis_convention':'H U, Bar[H/S1/S2] conjugate(U), epsilon_12=+1',
            'scalar_mass_squared_eigenvalues_gev2':basis.eigenvalues_sq.tolist(),
            **result,
            'verified_native_rgbeta_coupling_prefactors':False,
            'verified_native_rgbeta_matchete_lagrangian_sign':False,
            'physical_coupling_basis_map_certified':False,
            'full_mass_eigenstate_Matchete_matching_performed':False,
            'sequential_EFT_matching_performed':False}


def main(argv:list[str]|None=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    s=p.add_mutually_exclusive_group(required=True)
    s.add_argument('--all',action='store_true')
    s.add_argument('--model',choices=MODELS)
    p.add_argument('--run-id',required=True)
    p.add_argument('--output-root',type=Path,default=OUTPUT)
    args=p.parse_args(argv)
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id) is not None,
            'unsafe --run-id')
    src=(NATIVE,RGB_NATIVE,RGB_LIVE)
    for path in src:require(path.is_file(),'missing source '+str(path))
    data=[load(path) for path in src]
    folder=args.output_root/args.run_id
    require(not folder.exists(),'output already exists: choose a new --run-id')
    folder.mkdir(parents=True)
    rows=[];errors=[]
    for model in MODELS if args.all else (args.model,):
        try:
            result=run_one(model,mat_source=data[0],inv=data[1],live=data[2])
            tag=model.replace('-','m')
            result['source_sha256']={str(f.relative_to(ROOT)):digest(f) for f in src}
            result['mass_basis_report_sha256']=digest(MASS_BASIS/f'{tag}_final_mass_basis.json')
            (folder/f'{tag}.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            rows.append({'model':model,'status':result['status'],
                        'map':result['rgbeta_in_matchete_exact_columns'],
                        'rotated_component_residual':result['mass_basis_rotated_component_max_relative_error'],
                        'rotated_polynomial_residual':result['rgbeta_original_vs_rotated_polynomial_max_relative_error']})
            print(model+': exact 110101 rank-2 map checked; mass-basis residual=',
                  f'{result["mass_basis_rotated_component_max_relative_error"]:.3e}',flush=True)
        except (ValueError,KeyError,TypeError,OSError,RuntimeError) as e:
            errors.append({'model':model,'error':f'{type(e).__name__}: {e}'})
            print(model+': BLOCKED '+str(e),flush=True)
    summary={'format':'t3-bc-native-110101-mass-basis-tensor-map-v1',
             'status':'NativeTensorGeometryChecked' if not errors else 'IncompleteOrFailed',
             'models_requested':2 if args.all else 1,
             'models_passed':len(rows),'results':rows,'errors':errors,
             'source_files_sha256':{str(f.relative_to(ROOT)):digest(f) for f in src},
             'physical_coupling_mapping_certified':False,
             'sequential_EFT_matching_performed':False}
    (folder/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Summary:',folder/'summary.json',flush=True)
    return 2 if errors else 0


if __name__=='__main__':
    raise SystemExit(main())
