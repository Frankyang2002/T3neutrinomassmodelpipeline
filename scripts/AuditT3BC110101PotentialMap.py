"""B/C(-1) 110101 quartic potential normalization and sign (source reconstruction).

Read the saved live RGBeta invariant/coupling records, exact native Matchete
component tensors, and B/C scalar mass-basis reports. Reconstruct the RGBeta
real-to-complex quartic normalization following AddQuartic -> Lam/4! and
compare the RGBeta scalar POTENTIAL against minus the Matchete LAGRANGIAN.

The sign assumes the independently documented RGBeta generic Lagrangian
L_4 = -Lambda_abcd Phi_a Phi_b Phi_c Phi_d/24 (Thomsen, EPJC 81, 408,
2021, Eq. 2.1). The source-convention reconstruction is NOT an independent
runtime evaluation of installed RGBeta Lam or native Matchete Match[].

    python -u -X utf8 -m scripts.AuditT3BC110101PotentialMap --all --run-id bc_110101_potential_v1
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import re
from pathlib import Path

import numpy as np
import sympy as sp

from scripts.AuditT3BC110101MassBasisMap import (
    MODELS, MAT, RGB, SPECIES, NATIVE, RGB_NATIVE, RGB_LIVE, MASS_BASIS,
    load, check_live, exact_map, original_matchete_exact, source_tensor,
)
from scripts.AuditT3BCNativeQuarticMassBasis import load_basis, rotate_tensor

ROOT = Path(__file__).resolve().parents[1]
MAT_SOURCE = ROOT / 'Lagrangian/interactions/GeneralScalarOperators.wl'
OUTPUT = ROOT / 'output/audits/t3_bc_110101_potential_map'
NORMALIZATION = sp.Integer(6)
N_REAL_LEGS = 4
PERMUTATIONS = tuple(itertools.permutations(range(4)))
EXPECTED_A = sp.Matrix([[sp.Rational(1, 2), sp.Integer(1)],
                        [-sp.sqrt(3)/2, sp.Integer(0)]])


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('B/C 110101 potential normalization: '+message)


def source_guard(wolfram: str) -> None:
    """Check the Matchete-side source +h.c. convention is still in place."""
    compact = ''.join(wolfram.split())
    require('"Class"->If[degree===4,"ExtraQuartic","ExtraQuadratic"]' in compact,
            'Matchete extra quartic construction changed')
    require('"Expression"->PlusHc[contribution]' in compact,
            'Matchete scalar Lagrangian +PlusHc convention changed')
    require('contribution=coupling[]productIf[cg===None,1,CG[cg,labels]]' in compact,
            'Matchete unit extra quartic coefficient source changed')


def metadata_guard(native_rgbeta: dict, model: str) -> dict:
    records = native_rgbeta['models'][model]['invariants']
    chosen = [next((row for row in records if row.get('symbol') == name), None)
              for name in RGB]
    require(all(isinstance(row, dict) for row in chosen), 'missing native RGBeta couplings')
    for name, row in zip(RGB, chosen):
        factor=''.join(row['native_coupling_function_inputform'].split())
        require(factor==f'6*Matrix[{name}][##1]&',
                model+'/'+name+': RGBeta internal real-tensor factor is not 6')
        # Four distinguishable field species: all 24 permutations are unique.
        perms = row['native_unique_arrangements_inputform']
        require(isinstance(perms,str), 'missing RGBeta permutation witness')
        nums = re.findall(r'\{([1-4](?:,\s*[1-4]){3})\}', perms)
        parsed = [tuple(int(i) for i in s.split(',')) for s in nums]
        require(len(parsed)==24 and set(parsed)==set(itertools.permutations((1,2,3,4))),
                model+'/'+name+': expected exactly 24 native arrangements')
        require(row.get('self_conjugate') is False,
                model+'/'+name+': RGBeta must include Hermitian conjugate')
    return {'rgbeta_internal_coupling_factor':6, 'unique_arrangements':24,
            'hermitian_conjugate_enabled':True}


def normalization_identity() -> sp.Expr:
    """Lambda/4! with 6 and Phi_complex=(Phi_R+iPhi_I)/sqrt(2)."""
    return sp.simplify(NORMALIZATION/sp.factorial(N_REAL_LEGS) * sp.sqrt(2)**4)


def contract(t: np.ndarray, vectors: list[np.ndarray] | tuple[np.ndarray,...]) -> complex:
    require(t.shape==(2,)*4 and len(vectors)==4, 'incorrect native tensor dimensions')
    return complex(np.einsum('abcd,a,b,c,d->',t,*vectors,optimize=True))


def contraction_in_mass_basis(tensors: dict, new: dict) -> complex:
    return sum((contract(t,[new[f] for f in fields]) for fields,t in tensors.items()),0j)


def compute_potentials(mat:tuple[np.ndarray,np.ndarray],
                       rgb:tuple[np.ndarray,np.ndarray],
                       A:sp.Matrix, U:np.ndarray,
                       g:np.ndarray, hfield:np.ndarray,
                       x:np.ndarray, y:np.ndarray) -> dict[str,float | complex]:
    """Source-normalized V_RGB and L_M for arbitrary complex field values.

    For the four complex fields the RGBeta real tensor has factor 6. Its
    fully symmetrized potential is 6/4! times the product of the unnormalized
    real-combination fields Z=R+iI=sqrt(2) Phi. The two effects give 1.
    Both scalar interactions include h.c. The generic RGBeta L4=-V4 sign
    is then compared with Matchete's +PlusHc Lagrangian.
    """
    g=np.asarray(g,dtype=complex)
    require(g.shape==(2,), 'two independent complex RGBeta couplings required')
    coeff=np.asarray([[complex(sp.N(-A[i,j],35)) for j in range(2)]
                      for i in range(2)],dtype=complex)@g
    eps=np.array([[0.,1.],[-1.,0.]])
    s1=U[0,0]*x+U[0,1]*y
    s2=-eps@(U[1,0]*x+U[1,1]*y).conjugate()
    old={'H':hfield,'Hdag':hfield.conjugate(),
         'S1dag':s1.conjugate(),'S2dag':s2.conjugate()}
    new={'H':hfield,'Hdag':hfield.conjugate(),
         'X':x,'Xdag':x.conjugate(),'Y':y,'Ydag':y.conjugate()}
    z=[np.sqrt(2)*old[f] for f in SPECIES]
    # RGBeta's 1/4! real-scalar action with factor 6 and all 24 permutations.
    # Scalar factors commute: averaging over all permutations preserves each
    # invariant polynomial. We explicitly test that average below.
    r_polys=np.array([contract(t,z) for t in rgb])
    rgb_ordered=complex(np.dot(g,r_polys))
    rgb_fourth=complex(float(NORMALIZATION)/24.0 * rgb_ordered)
    rgb_potential=float(2*np.real(rgb_fourth))
    rgb_complex=complex(np.dot(g,[contract(t,[old[f] for f in SPECIES]) for t in rgb]))
    mat_lagrangian=float(2*np.real(np.dot(coeff,[contract(t,[old[f] for f in SPECIES]) for t in mat])))
    rot_rgb=[rotate_tensor(t,SPECIES,U) for t in rgb]
    rot_mat=[rotate_tensor(t,SPECIES,U) for t in mat]
    rgb_rotated=float(2*np.real(np.dot(g,[contraction_in_mass_basis(t,new) for t in rot_rgb])))
    mat_rotated=float(2*np.real(np.dot(coeff,[contraction_in_mass_basis(t,new) for t in rot_mat])))
    return {'rgbeta_v_from_real_lam':rgb_potential,
            'rgbeta_v_from_complex_invariants':float(2*np.real(rgb_complex)),
            'matchete_l_from_native_hc':mat_lagrangian,
            'rgbeta_v_in_mass_basis':rgb_rotated,
            'matchete_l_in_mass_basis':mat_rotated,
            'complex_rgbeta_polynomial':rgb_complex,
            'complex_matchete_polynomial':complex(np.dot(coeff,[contract(t,[old[f] for f in SPECIES]) for t in mat]))}


def check_potential_ensemble(mat:tuple[np.ndarray,np.ndarray],
                             rgb:tuple[np.ndarray,np.ndarray],
                             A:sp.Matrix, U:np.ndarray,
                             *, trials:int=32, seed:int=331) -> dict:
    require(trials>=4,'at least four trials required')
    rng=np.random.default_rng(seed)
    samples=[np.array([1+0j,0+0j]),np.array([0+0j,1+0j]),
             np.array([0.12+0.84j, 0+0j]),np.array([0+0j,-0.45+0.33j])]
    while len(samples)<trials:
        samples.append(rng.normal(size=2)+1j*rng.normal(size=2))
    max_norm=max_sign=max_rot=max_complex=0.
    positive_cases=0
    for g in samples:
        rand=lambda: rng.normal(size=2)+1j*rng.normal(size=2)
        vals=compute_potentials(mat,rgb,A,U,g,rand(),rand(),rand())
        v=vals['rgbeta_v_from_real_lam']; c=vals['rgbeta_v_from_complex_invariants']
        l=vals['matchete_l_from_native_hc']; rv=vals['rgbeta_v_in_mass_basis']
        rl=vals['matchete_l_in_mass_basis']
        scale=max(1.,abs(v),abs(c),abs(l),abs(rv),abs(rl))
        max_norm=max(max_norm,abs(v-c)/scale)
        max_sign=max(max_sign,abs(v+l)/scale)
        max_rot=max(max_rot,abs(v-rv)/scale,abs(l-rl)/scale)
        complex_scale=max(1.,abs(vals['complex_rgbeta_polynomial']),
                          abs(vals['complex_matchete_polynomial']))
        max_complex=max(max_complex,abs(vals['complex_rgbeta_polynomial']+
                                        vals['complex_matchete_polynomial'])/complex_scale)
        if abs(v)>1e-8:positive_cases+=1
    require(positive_cases>=trials-2, 'insufficient nonzero potential samples')
    require(max_norm<5e-13 and max_sign<5e-13 and
            max_rot<5e-13 and max_complex<5e-13,
            f'potential conversion failed: norm={max_norm}, sign={max_sign}, mass={max_rot}, complex={max_complex}')
    return {'trials':trials,'max_rgbeta_real_to_complex_potential_relative_error':max_norm,
            'max_rgbeta_potential_plus_matchete_lagrangian_relative_error':max_sign,
            'max_mass_basis_potential_relative_error':max_rot,
            'max_complex_polynomial_relative_error':max_complex}


def run_one(model:str, native:dict, inv:dict, live:dict,
            *, mass_root:Path=MASS_BASIS,
            mat_source_path:Path=MAT_SOURCE, trials:int=32) -> dict:
    require(model in MODELS,'unsupported model')
    source_guard(mat_source_path.read_text(encoding='utf-8-sig'))
    meta=metadata_guard(inv,model)
    rgb=check_live(live,inv,model)
    A=exact_map(native,model,rgb)
    require(all(sp.simplify(x)==0 for x in A-EXPECTED_A),
            model+': exact 110101 tensor map changed')
    m=original_matchete_exact(native,model)
    basis=load_basis(mass_root/f'{model.replace("-","m")}_final_mass_basis.json',model)
    require(not basis.tachyonic and not basis.degenerate,'invalid scalar mass basis')
    identity=normalization_identity()
    require(identity==1,'RGBeta real-complex normalization no longer equals one')
    result=check_potential_ensemble(m,rgb,A,basis.U,trials=trials)
    return {'model':model,'status':'SourceReconstructedQuarticPotentialMapChecked',
            'sector':'110101','rgbeta_symbols':list(RGB),'matchete_symbols':list(MAT),
            'rgbeta_in_matchete_exact_columns':[[str(A[i,j]) for j in range(2)] for i in range(2)],
            'rgbeta_potential_normalization_from_6_over_24_and_four_sqrt2':str(identity),
            'matchete_lagrangian_couplings_in_terms_of_rgbeta_potential':[
                ['-'+str(A[i,j]) for j in range(2)] for i in range(2)],
            'coupling_map_exact_h1':'-g1/2 - g2',
            'coupling_map_exact_h2':'sqrt(3)*g1/2',
            'source_argument':'RGBeta generic L4=-Lambda Phi^4/24; AddQuartic factor=6; '
                              'Matchete PlusHc coefficient=+1; scalar potential V=-L',
            'rgbeta_package_source_version_verified_locally':False,
            'installed_rgbeta_lam_runtime_reconstruction_performed':False,
            'installed_matchete_lagrangian_runtime_evaluated':False,
            'physical_coupling_mapping_certified':False,
            'one_loop_matching_performed':False,
            'sequential_eft_matching_performed':False,
            **meta,**result}


def digest(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv:list[str]|None=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    choose=parser.add_mutually_exclusive_group(required=True)
    choose.add_argument('--all',action='store_true')
    choose.add_argument('--model',choices=MODELS)
    parser.add_argument('--run-id',required=True)
    parser.add_argument('--output-root',type=Path,default=OUTPUT)
    parser.add_argument('--trials',type=int,default=32)
    args=parser.parse_args(argv)
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id) is not None,
            'unsafe run ID')
    sources=[NATIVE,RGB_NATIVE,RGB_LIVE,MAT_SOURCE]
    for source in sources:require(source.is_file(),'missing input '+str(source))
    native,inv,live=[load(f) for f in sources[:3]]
    folder=args.output_root/args.run_id
    require(not folder.exists(), 'output already exists, choose a new --run-id')
    rows=[];errors=[]
    for model in MODELS if args.all else (args.model,):
        try:
            row=run_one(model,native,inv,live,trials=args.trials)
            report={**row,'source_sha256':{f.relative_to(ROOT).as_posix():digest(f) for f in sources},
                    'mass_basis_sha256':digest(MASS_BASIS/f'{model.replace("-","m")}_final_mass_basis.json')}
            tag=model.replace('-','m')
            rows.append({'model':model,'status':row['status'],
                         'potential_sign_residual':row['max_rgbeta_potential_plus_matchete_lagrangian_relative_error'],
                         'mass_basis_residual':row['max_mass_basis_potential_relative_error']})
            folder.mkdir(parents=True,exist_ok=True)
            (folder/f'{tag}.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            print(f'{model}: source reconstructed 110101 potential; sign residual='
                  f'{row["max_rgbeta_potential_plus_matchete_lagrangian_relative_error"]:.3e}; '
                  f'mass-basis residual={row["max_mass_basis_potential_relative_error"]:.3e}',flush=True)
        except (ValueError,KeyError,TypeError,RuntimeError,OSError) as exc:
            errors.append({'model':model,'error':str(exc)})
            print(f'{model}: BLOCKED {exc}',flush=True)
    folder.mkdir(parents=True,exist_ok=True)
    summary={'format':'t3-bc-110101-potential-map-reconstruction-v1',
             'status':'SourcePotentialReconstructionChecked' if not errors else 'IncompleteOrFailed',
             'models_requested':2 if args.all else 1,'models_passed':len(rows),
             'results':rows,'errors':errors,
             'physical_coupling_mapping_certified':False,
             'native_runtime_action_checked':False,
             'sequential_eft_matching_performed':False}
    (folder/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Summary:',folder/'summary.json',flush=True)
    return 2 if errors else 0

if __name__=='__main__':
    raise SystemExit(main())
