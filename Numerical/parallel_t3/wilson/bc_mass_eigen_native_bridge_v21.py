"""A v21: source-guarded B/C(-1) exact Matchete complex quartic mass-basis bridge.

REUSES the repository's already implemented MixedDoubletMassBasis and field_map;
never invents another diagonalizer or declares physical Matchete/RGBeta matching.
Only source CGs from v18 and the conditional six-of-seven v19 map are used.
"""
from __future__ import annotations

from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from typing import Any
import re
import numpy as np

from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
    EPS, MixedDoubletMassBasis, field_map, quadratic_mass_tensor_transport,
    verify_su2_covariance,
)
from Numerical.parallel_t3.wilson.complex_scalar_bridge_v19 import (
    BC, MAT_FORMAT, MAT_STATUS, SPECIAL, V19Blocked,
    native_document, sha as digest, matchete_number,
    FORMAT as V19_FORMAT, STATUS as V19_STATUS,
)
from Numerical.parallel_t3.wilson.full_scalar_census_v20 import (
    FORMAT as V20_FORMAT, STATUS as V20_STATUS,
)

FORMAT = 't3-parallel-a-bc-complex-cg-mass-eigen-transport-v21'
STATUS = 'BCMinus1MatcheteNativeEightComplexCGsPerModelRotated_MassEigenActionVerifiedConditional'
MASS_FORMAT = 't3-mixed-doublet-mass-basis-transport-v1'
MASS_STATUS = 'ExperimentalMassBasisTensorTransportNoPhysicalMatching'
STATES = ('H','Hdag','S1','S1dag','S2','S2dag')
EIGEN = ('H','Hdag','X','Xdag','Y','Ydag')

class BCBridgeBlocked(ValueError):
    pass

def require(ok: Any, msg: str) -> None:
    if not bool(ok):
        raise BCBridgeBlocked('Workstream A v21: '+msg)

def read_json(path: Path) -> dict:
    path=Path(path)
    require(path.is_file(),'missing witness '+str(path))
    try:
        doc=json.loads(path.read_text(encoding='utf-8-sig'),
                       parse_constant=lambda x: (_ for _ in ()).throw(BCBridgeBlocked('nonfinite JSON '+x)))
    except (ValueError,TypeError) as exc:
        raise BCBridgeBlocked('invalid witness JSON '+str(path)) from exc
    require(isinstance(doc,dict),'witness not an object: '+str(path))
    return doc

def finite_matrix(raw: dict, dims: tuple[int,int]) -> np.ndarray:
    require(isinstance(raw,dict) and set(raw)=={'real','imag'},'split complex matrix schema')
    try:
        real=np.asarray(raw['real'],dtype=float)
        imag=np.asarray(raw['imag'],dtype=float)
    except (TypeError,ValueError,OverflowError) as exc:
        raise BCBridgeBlocked('invalid split complex matrix') from exc
    require(real.shape==dims and imag.shape==dims and
            np.isfinite(real).all() and np.isfinite(imag).all(),
            'invalid complex matrix values/shape')
    return real+1j*imag


def finite_v19_kernel(raw: dict) -> np.ndarray:
    """Read the v19 SVD null vector, which v19 writes as a flat 7-vector.

    v19 complex_scalar_bridge_v19.py computes `null = vh[-1].conj()`;
    `.tolist()` produces real/imag arrays of shape (7,), NOT (1,7).
    Accept the older synthetic fixture's single row too, but never allow
    NumPy broadcasting or a column/transposed/multidimensional kernel.
    These are representation adapters only: rank and B@n checks follow.
    """
    require(isinstance(raw,dict) and set(raw)=={'real','imag'},
            'v19 kernel: split complex vector schema')
    try:
        real=np.asarray(raw['real'],dtype=float)
        imag=np.asarray(raw['imag'],dtype=float)
    except (TypeError,ValueError,OverflowError) as exc:
        raise BCBridgeBlocked('v19 kernel: nonnumeric split complex vector') from exc
    require(real.shape==imag.shape and real.shape in ((7,), (1,7)) and
            np.isfinite(real).all() and np.isfinite(imag).all(),
            'v19 kernel: expected seven entries (or one legacy synthetic row), finite real/imag shapes equal')
    return (real+1j*imag).reshape(7)

def require_false_flags(doc:dict, keys:tuple[str,...], label:str)->None:
    for key in keys:
        require(doc.get(key) is False,label+': missing or false certification guard '+key)

def read_mass_basis(root:Path,model:str,path:Path)->tuple[MixedDoubletMassBasis,dict]:
    data=read_json(path)
    require(data.get('format')==MASS_FORMAT and data.get('status')==MASS_STATUS and
            data.get('model')==model and data.get('endpoint')=='final',
            model+': incompatible saved final scalar mass basis')
    require_false_flags(data,('physical_mass_matching_performed',
           'MatcheteRGBetaPhysicalCouplingMapVerified',
           'MatcheteHeavyMassBasisCheckOnRedefinedLagrangian'),model+' mass source')
    require(data.get('quadratic_hessian_agreement') is True and
            data.get('source_quartic_count')==17 and
            data.get('max_source_tensor_transport_relative_error',float('inf'))<3e-12,
            model+': incomplete existing mass-basis source audit')
    hashes=data.get('source_sha256')
    require(isinstance(hashes,dict) and hashes,model+': mass source hash inventory missing')
    for rel,h in hashes.items():
        require(isinstance(rel,str) and not Path(rel).is_absolute() and '..' not in Path(rel).parts and
                (root/rel).is_file() and digest(root/rel)==h,
                model+': mass source hash mismatch '+str(rel))
    raw=data.get('basis',{})
    require(raw.get('model')==model and raw.get('tachyonic') is False and
            raw.get('degenerate_basis_nonunique') is False and
            raw.get('matchete_field_redefinition_installed') is False and
            raw.get('rgbeta_matchete_component_map_certified') is False and
            raw.get('threshold_matching_performed') is False and
            raw.get('epsilon')==EPS.real.tolist() and
            raw.get('eigenstate_hypercharge')==-0.5,
            model+': unsafe or incompatible existing eigenbasis')
    A=finite_matrix(raw.get('mass_squared_flavour_matrix_gev2'),(2,2))
    U=finite_matrix(raw.get('unitary_flavour_rotation_U'),(2,2))
    require(np.allclose(A,A.conj().T,rtol=0,atol=2e-7),
            model+': mass source not Hermitian')
    basis=MixedDoubletMassBasis.build(model,float(A[0,0].real),float(A[1,1].real),complex(A[0,1]))
    require(not basis.tachyonic and not basis.degenerate and
            np.all(basis.eigenvalues_sq>0) and
            np.allclose(basis.A,A,rtol=1e-12,atol=2e-7) and
            np.allclose(basis.U,U,rtol=1e-12,atol=2e-10) and
            np.allclose(basis.eigenvalues_sq,raw.get('masses_squared_gev2'),rtol=1e-11,atol=2e-6),
            model+': saved scalar mixing rotation/spectrum differs from source quadratic')
    return basis,data

def convert_fields(fields:list[str]|tuple[str,...])->tuple[str,...]:
    mapping={'H':'H','Bar[H]':'Hdag','Hdag':'Hdag',
             'S1':'S1','Bar[S1]':'S1dag','S1dag':'S1dag',
             'S2':'S2','Bar[S2]':'S2dag','S2dag':'S2dag'}
    require(len(fields)==4 and all(f in mapping for f in fields),'unsupported native CG field list')
    return tuple(mapping[f] for f in fields)

def tensor_record(record:dict)->tuple[str,tuple[str,...],np.ndarray]:
    sym=record.get('symbol')
    dims=record.get('dimensions');vals=record.get('components_inputform')
    require(isinstance(sym,str) and sym and dims==[2]*4 and
            isinstance(vals,list) and len(vals)==16 and
            record.get('assembled_lagrangian_contains_operator') is True,
            'missing/invalid native doublet CG record')
    t=np.asarray([matchete_number(v) for v in vals],complex).reshape((2,)*4)
    require(np.isfinite(t).all() and np.linalg.norm(t)>1e-12,'native Matchete CG zero/invalid')
    return sym,convert_fields(record['fields']),t

def rotate_tensor(t:np.ndarray,fields:tuple[str,...],U:np.ndarray)->dict[tuple[str,...],np.ndarray]:
    require(t.shape==(2,2,2,2) and len(fields)==4 and all(f in STATES for f in fields),
            'only four explicit SU2 doublet axes permitted')
    out:dict[tuple[str,...],np.ndarray]={}
    for choice in product(*(field_map(f,U) for f in fields)):
        labels=tuple(k for k,_ in choice)
        vals=np.einsum('abcd,ai,bj,ck,dl->ijkl',t,*(p for _,p in choice),optimize=True)
        out[labels]=out.get(labels,np.zeros((2,)*4,dtype=complex))+vals
    require(out and all(k in EIGEN for fields in out for k in fields),
            'rotated tensor contains unsupported mass eigenfield')
    return out

def contract(t:np.ndarray, fields:tuple[str,...],vectors:dict[str,np.ndarray])->complex:
    return complex(np.einsum('abcd,a,b,c,d->',t,*(vectors[f] for f in fields),optimize=True))

def sample_vectors(basis:MixedDoubletMassBasis,rng:np.random.Generator):
    h=rng.normal(size=2)+1j*rng.normal(size=2)
    x=rng.normal(size=2)+1j*rng.normal(size=2)
    y=rng.normal(size=2)+1j*rng.normal(size=2)
    s1,s2=basis.original_fields(x,y)
    old={'H':h,'Hdag':h.conj(),'S1':s1,'S1dag':s1.conj(),
         'S2':s2,'S2dag':s2.conj()}
    new={'H':h,'Hdag':h.conj(),'X':x,'Xdag':x.conj(),
         'Y':y,'Ydag':y.conj()}
    return old,new

def action_residual(t:np.ndarray,fields:tuple[str,...],rotated:dict,
                    basis:MixedDoubletMassBasis,*,seed:int=123,trials:int=8)->float:
    require(1<=trials<=32 and isinstance(seed,int),'invalid independent action-sampling request')
    rng=np.random.default_rng(seed);worst=0.
    for _ in range(trials):
        old,new=sample_vectors(basis,rng)
        coeff=complex(*rng.normal(size=2))
        direct=2*np.real(coeff*contract(t,fields,old))
        transformed=2*np.real(coeff*sum(contract(z,k,new) for k,z in rotated.items()))
        # PlusHc coefficient 1: Matchete complex quartic convention.
        err=abs(direct-transformed)/max(1.,abs(direct),abs(transformed))
        worst=max(worst,float(err))
        require(err<2e-11,'source native complex tensor fails mass-basis PlusHc identity')
    return worst

def special_source_permutation(v19: dict, model: str,
                               source_symbols: list[str]) -> tuple[list[str], list[int]]:
    """Return v19 historical coupling order and column indices in NEW native order.

    v19 used the historical Stage22 Matchete order; v18's independent native
    exporter uses the sorted verified-structure order. These span the SAME
    seven-dimensional source space but need not list its axes in the same order.
    Match by exact coupling symbol, never by position or tensor norm.
    """
    d = v19.get('BC_baseline_six_of_seven_conditional_source_maps', {}).get(model, {})
    historical = d.get('matchete_seven_source_symbols')
    require(isinstance(historical, list) and len(historical) == 7 and
            all(isinstance(s, str) and s.startswith('lambdaT3ExtraQ') for s in historical) and
            len(set(historical)) == 7,
            model+': invalid/duplicate historical Matchete coupling list in v19')
    require(isinstance(source_symbols, list) and len(source_symbols) == 7 and
            all(isinstance(s, str) for s in source_symbols) and
            len(set(source_symbols)) == 7 and set(source_symbols) == set(historical),
            model+': native v18 and historical v19 Matchete coupling inventories differ')
    return historical, [historical.index(s) for s in source_symbols]


def validate_special_matrix(v19:dict,model:str,source_symbols:list[str])->tuple[np.ndarray,np.ndarray]:
    d=v19.get('BC_baseline_six_of_seven_conditional_source_maps',{}).get(model,{})
    historical, indices=special_source_permutation(v19,model,source_symbols)
    require(d.get('baseline_rgbeta_six_symbols')==list(SPECIAL) and
            isinstance(d.get('optional_extra_rgbeta_symbol'),str) and
            d['optional_extra_rgbeta_symbol'].startswith('lambdaT3GeneralQ') and
            d.get('baseline_map_rank')==6 and d.get('optional_map_rank')==7 and
            d.get('baseline_kernel_dimension_in_Matchete_seven_space')==1 and
            d.get('source_monomial_not_native_projector_normalization') is True,
            model+': invalid v19 baseline six-of-seven conditional map')
    historical_B=finite_matrix(d.get('baseline_rgbeta_from_Matchete_seven_complex_matrix'),(6,7))
    historical_N=finite_v19_kernel(d.get('kernel_direction_normalized_complex'))
    require(np.linalg.matrix_rank(historical_B,tol=1e-10)==6 and
            np.linalg.norm(historical_B@historical_N)<1e-9 and
            abs(np.linalg.norm(historical_N)-1)<1e-10,
            model+': v19 conditional complex subspace is rank defective')
    # Same source action in the new v18 order: r = M_old q_old =
    # M_old[:, indices] q_new; the null vector is relabelled identically.
    B=historical_B[:,indices]
    N=historical_N[indices]
    require(np.linalg.matrix_rank(B,tol=1e-10)==6 and
            np.linalg.norm(B@N)<1e-9 and abs(np.linalg.norm(N)-1)<1e-10 and
            np.allclose(B,historical_B[:,indices],rtol=0,atol=0) and
            np.array_equal(N,historical_N[indices]),
            model+': source-symbol reordering broke v19 rank or kernel')
    return B,N

def verify_v19_input_sha256(root: Path, witness: dict,
                            required_files: tuple[Path, ...]) -> dict:
    """Recheck original bytes; accept v19's Windows and POSIX path spellings.

    v19 wrote `str(path.relative_to(root))`, which uses the host separator.
    v21 formerly looked up only `.as_posix()` names, an incorrect dictionary
    lookup on Windows.  Normalize the KEYS only, never the SHA256 digests.
    Reject aliased keys, path traversal, missing files, and changed content.
    Check all seven original v19 inputs, not just the two native CG exports.
    """
    root = Path(root).resolve()
    original = witness.get('input_sha256')
    require(isinstance(original, dict) and len(original) == 7,
            'v19 must provide exactly seven independently hashed input files')
    normalized: dict[str, str] = {}
    windows_keys = 0
    for raw, expected in original.items():
        require(isinstance(raw, str) and bool(raw) and
                isinstance(expected, str) and
                re.fullmatch(r'[0-9a-f]{64}', expected) is not None,
                'v19 contains malformed file name or SHA256 digest')
        # This is a source-metadata normalization, NOT a redefinition of a
        # physical tensor/component convention or a rewrite of witness JSON.
        name = raw.replace('\\', '/')
        parts = name.split('/')
        require(all(part not in ('', '.', '..') for part in parts) and
                not name.startswith('/') and ':' not in name and
                not any(ord(c) < 32 or ord(c) == 127 for c in name),
                'v19 contains unsafe or noncanonical project-relative source key: '+repr(raw))
        require(name not in normalized,
                'v19 contains both Windows/POSIX aliases of one source: '+name)
        normalized[name] = expected
        windows_keys += ('\\' in raw)
    for relative, expected in normalized.items():
        path = root.joinpath(*relative.split('/'))
        require(path.is_file(), 'v19 provenance input missing: '+relative)
        actual = digest(path)
        require(actual == expected,
                'v19 provenance SHA256 actually differs for '+relative+
                ' (saved '+expected[:16]+', current '+actual[:16]+')')
    for f in required_files:
        path = Path(f).resolve()
        require(path.is_relative_to(root),
                'native CG evidence must be inside the current repository: '+str(path))
        rel = path.relative_to(root).as_posix()
        require(rel in normalized,
                'v19 has no source hash for native CG: '+rel)
    return {'source_files_sha256_rechecked': len(normalized),
            'legacy_windows_separator_keys_accepted': windows_keys,
            'provenance_hashes_unchanged': True}


def gauge_covariance_recheck(basis:MixedDoubletMassBasis)->dict:
    gauge=float(verify_su2_covariance(basis,seed=102,trials=12))
    mass=quadratic_mass_tensor_transport(basis)
    scale=max(1.,float(np.max(abs(basis.eigenvalues_sq))))
    require(gauge<3e-12 and mass['mass_tensor_diagonalization_residual_gev2']<2e-10*scale,
            'quadratic tensor/SU2 covariance no longer valid')
    return {'quadratic_tensor_transport_max_absolute_gev2':float(mass['mass_tensor_diagonalization_residual_gev2']),
            'su2_covariance_max_absolute':gauge,'quadratic_masses_gev':basis.info()['masses_gev'],
            'explicit_epsilon_conjugation_tested':True}

def audit_model(model:str,basis:MixedDoubletMassBasis,topology:dict,bc:dict,
                v19:dict,*,seed:int=109,trials:int=8)->dict:
    require(model in BC and basis.model==model,'only B_-1/C_-1 supported')
    source=[topology['models'][model]['records'][0]]+bc['models'][model]['records']
    require(len(source)==8,'one topology + seven optional native CGs required')
    converted=[tensor_record(rec) for rec in source]
    symbols=[item[0] for item in converted]
    require(symbols[0]=='lambdaT3' and len(set(symbols))==8,'native Matchete CG naming duplicate')
    B,N=validate_special_matrix(v19,model,symbols[1:])
    old_source_order,native_column_indices=special_source_permutation(v19,model,symbols[1:])
    q=gauge_covariance_recheck(basis)
    all_rot=[];total=0;worst=0.;entries=[]
    for i,(sym,fields,t) in enumerate(converted):
        rotated=rotate_tensor(t,fields,basis.U)
        err=action_residual(t,fields,rotated,basis,seed=seed+41*i,trials=trials)
        worst=max(worst,err);total+=len(rotated)
        all_rot.append(rotated)
        entries.append({'source_matchete_coupling':sym,'source_fields':list(fields),
                        'native_CG_frobenius_norm':float(np.linalg.norm(t)),
                        'rotated_mass_basis_monomials':len(rotated),
                        'mass_basis_fields':[list(x) for x in sorted(rotated)],
                        'plus_Hc_action_relative_residual':err,
                        'numerical_RGBeta_component_tensor_assumed':False})
    # The baseline-six map has one optional source direction invisible to it.
    # Check that this seventh direction is NOT silently discarded by the
    # invertible mass-eigenfield transformation.
    rng=np.random.default_rng(seed+901)
    max_kernel_action=0.;max_kernel_holdout_error=0.
    for i in range(trials):
        old,new=sample_vectors(basis,rng)
        coeff=complex(*rng.normal(size=2))
        original=sum(N[j]*contract(converted[j+1][2],converted[j+1][1],old)
                     for j in range(7))
        moved=sum(N[j]*sum(contract(t,f,new) for f,t in all_rot[j+1].items())
                  for j in range(7))
        a=2*np.real(coeff*original);b=2*np.real(coeff*moved)
        max_kernel_action=max(max_kernel_action,abs(a))
        e=abs(a-b)/max(1.,abs(a),abs(b))
        max_kernel_holdout_error=max(max_kernel_holdout_error,float(e))
        require(e<2e-11,'optional seventh source direction not covariantly transported')
    require(max_kernel_action>1e-6,
            'optional seventh source CG direction appears numerically null: cannot exclude it')
    return {'model':model,'source_native_complex_CG_count':8,
            'source_topology_CGs':1,'source_optional_extra_CGs':7,
            'baseline_RGBeta_special_complex_slots':6,
            'baseline_six_from_seven_rank':int(np.linalg.matrix_rank(B,tol=1e-10)),
            'v19_historical_seven_source_symbol_order':old_source_order,
            'v18_native_seven_source_symbol_order':symbols[1:],
            'v19_column_indices_in_native_v18_order':native_column_indices,
            'source_basis_reordered_by_exact_coupling_symbol':True,
            'optional_seventh_source_direction_action_nonzero':True,
            'optional_seventh_source_direction_worst_holdout_residual':max_kernel_holdout_error,
            'mass_rotation_analytic_and_tensor_checks':q,
            'rotated_mass_basis_monomials_before_dedup':total,
            'max_source_native_action_relative_residual':worst,
            'native_CG_sectors':entries,
            'physical_component_coupling_map_certified':False,
            'physical_scalar_matching_performed':False}

def audit(root:Path,*,topology_path:Path,bc_path:Path,mass_paths:dict[str,Path],
          v19_path:Path,v20_path:Path,trials:int=8)->dict:
    root=Path(root).resolve()
    require(2<=trials<=30,'invalid trial count')
    top=read_json(topology_path);bc=read_json(bc_path)
    native_document(top,root,'TOPOLOGY');native_document(bc,root,'BCEXTRA')
    v19=read_json(v19_path);v20=read_json(v20_path)
    require(v19.get('format')==V19_FORMAT and v19.get('status')==V19_STATUS and
            v19.get('baseline_EFT1_complex_slots')==28 and
            v19.get('certification',{}).get('physical_Wilson_evolution_authorized') is False,
            'v19 source conditioned complex-coupling map invalid')
    v19_provenance=verify_v19_input_sha256(root,v19,(topology_path,bc_path))
    require(v20.get('format')==V20_FORMAT and v20.get('status')==V20_STATUS and
            v20.get('total_EFT1_scalar_quartic_registration_slots')==178 and
            v20.get('native_RGBeta_real_quartics')==150 and
            v20.get('native_RGBeta_complex_quartics')==28 and
            v20.get('certification',{}).get('physical_Wilson_evolution_authorized') is False and
            v20.get('certification',{}).get('BC_minus1_mixed_scalar_mass_basis_certified') is False and
            v20.get('input_sha256',{}).get('complex_bridge_v19')==digest(v19_path),
            'v20 source inventory not verified or v19 source hash differs')
    models={};worst=0.
    for i,tag in enumerate(BC):
        require(tag in mass_paths,'missing existing mass-basis report '+tag)
        basis,_=read_mass_basis(root,tag,mass_paths[tag])
        m=audit_model(tag,basis,top,bc,v19,seed=390+i*500,trials=trials)
        models[tag]=m;worst=max(worst,m['max_source_native_action_relative_residual'])
    return {'format':FORMAT,'status':STATUS,'models_checked':2,
            'native_Matchete_complex_CG_sources_rotated':16,
            'native_topology_CGs_rotated':2,
            'native_BC_optional_extra_CGs_rotated':14,
            'baseline_RGBeta_special_complex_quartic_slots_conditionally_mapped':12,
            'max_complex_quartic_action_transport_relative_residual':worst,
            'v19_input_provenance':v19_provenance,
            'source_sha256':{'v18_topology':digest(topology_path),'v18_bc_extra':digest(bc_path),
                             'v19_complex_bridge':digest(v19_path),'v20_quartic_census':digest(v20_path),
                             **{tag+'_mass_basis':digest(mass_paths[tag]) for tag in BC}},
            'models':models,
            'certification':{
                'native_Matchete_complex_quartic_CG_source_consumed':True,
                'experimental_BC_scalar_doublet_mass_basis_reused':True,
                'exact_epsilon_and_complex_quartic_field_transport_checked':True,
                'matchete_redefined_mass_eigenstate_lagrangian_certified':False,
                'complete_physical_Matchete_RGBeta_coupling_map_certified':False,
                'signed_GammaCCWeyl_certified':False,
                'one_loop_Wilson_ADM_physical_certified':False,
                'physical_Wilson_evolution_authorized':False,
                'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT'}}
