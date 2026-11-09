"""Workstream A v25: cross-check v21/v22 against EXISTING native Matchete B/C eigenfield witnesses.

Read-only. Matchete's ordinary heavy-basis validation is evidence of a
well-formed eigenfield action, not a fermion-first scalar EFT match.
"""
from __future__ import annotations
from hashlib import sha256
from pathlib import Path
from typing import Any
import json
import math
import numpy as np

BC=('B_-1','C_-1')
FORMAT='t3-parallel-a-native-bc-mass-eigen-reconciliation-v25'
STATUS='BCEigenfieldMatcheteNativeBasisAndV21CGTransportConcordant_NoPhysicalMatching'
MASS_FORMAT='t3-mixed-doublet-mass-basis-transport-v1'
PARAM_FORMAT='t3-experimental-mass-eigen-matchete-parameters-v1'
NATIVE_FORMAT='t3-matchete-native-heavy-mass-eigenbasis-v1'
NATIVE_STATUS='MatcheteMassEigenHeavyBasisPassedPhysicalMatchingNotPerformed'
EPS=np.array([[0,1],[-1,0]],dtype=float)

class GateError(ValueError): pass

def require(pred:bool, reason:str)->None:
    if not bool(pred): raise GateError('v25: '+reason)

def load(path:Path)->dict:
    require(path.is_file(),'missing input '+str(path))
    o=json.loads(path.read_text(encoding='utf-8-sig'),parse_constant=lambda x: (_ for _ in ()).throw(GateError('nonfinite JSON constant '+x)))
    require(isinstance(o,dict),'JSON root is not object '+str(path))
    return o

def digest(path:Path)->str: return sha256(path.read_bytes()).hexdigest()

def cp(obj:Any, name:str)->complex:
    require(isinstance(obj,dict) and set(obj)=={'real','imag'},name+': expected complex pair')
    try: z=complex(float(obj['real']),float(obj['imag']))
    except (ValueError,TypeError,OverflowError) as ex: raise GateError(name+': invalid complex pair') from ex
    require(np.isfinite([z.real,z.imag]).all(),name+': nonfinite pair')
    return z

def cm(obj:Any, size:int, name:str)->np.ndarray:
    require(isinstance(obj,dict) and set(obj)=={'real','imag'},name+': expected real/imag matrix')
    try: re=np.array(obj['real'],dtype=float);im=np.array(obj['imag'],dtype=float)
    except (TypeError,ValueError,OverflowError) as ex:raise GateError(name+': nonnumeric entries') from ex
    require(re.shape==im.shape==(size,size) and np.isfinite(re).all() and np.isfinite(im).all(),name+': invalid matrix shape/data')
    return re+1j*im

def relpath(raw:str)->str:
    require(isinstance(raw,str),'expected relative path')
    path=raw.replace('\\','/')
    require(not path.startswith('/') and ':' not in path and all(q not in ('','.','..') for q in path.split('/')), 'unsafe provenance path')
    return path

def recheck_v21_v22(root:Path, v21_path:Path, v22_path:Path)->tuple[dict,dict]:
    # Retain the v22 source-hash guard rather than reimplementing a weaker one.
    from .bc_mass_bridge_gate_v22 import audit as v22_audit, FORMAT as F22, STATUS as S22
    from .bc_mass_eigen_native_bridge_v21 import FORMAT as F21, STATUS as S21
    v21=load(v21_path);v22=load(v22_path)
    require(v21.get('format')==F21 and v21.get('status')==S21,'v21 status not certified')
    require(v22.get('format')==F22 and v22.get('status')==S22,'v22 status not certified')
    base=root/'output/parallel_t3/wilson'
    fresh=v22_audit(root,base/'full_scalar_census_conditional_v20.json',v21_path)
    require(fresh==v22,'saved v22 no longer equals source-guarded recomputation')
    return v21,v22

def reconcile_model(root:Path,model:str,v21:dict, *,base_dir:Path|None=None)->dict:
    require(model in BC,'only B_-1/C_-1')
    source=(root/'output/general_t3_mass_basis' if base_dir is None else base_dir)
    native=root/'output/general_t3_mass_eigen_matchete'
    stem=model[0]+'_m1_final'
    mass_path=source/(stem+'_mass_basis.json')
    param_path=native/(stem+'_parameters.json')
    mat_path=native/(stem+'_matchete.json')
    mass,param,mat=(load(q) for q in (mass_path,param_path,mat_path))
    require(mass.get('format')==MASS_FORMAT and mass.get('model')==model and mass.get('endpoint')=='final',model+': invalid original mass report')
    require(mass.get('physical_mass_matching_performed') is False and
            mass.get('MatcheteRGBetaPhysicalCouplingMapVerified') is False,model+': unsafe mass claim')
    require(param.get('format')==PARAM_FORMAT and param.get('status')=='NumericMassEigenstateBoundaryPreparedMatcheteNotRun' and
            param.get('model')==model and param.get('endpoint')=='final',model+': invalid native parameters')
    require(param.get('source_sha256')==digest(mass_path) and
            relpath(param.get('source_path'))=='output/general_t3_mass_basis/'+stem+'_mass_basis.json',model+': parameter provenance changed')
    require(param.get('dimensions_alpha')==[2,2,1 if model=='B_-1' else 3,-1] and
            param.get('quartic_source_count')==17 and param.get('yukawa_source_count')==2,model+': field/coupling inventory changed')
    for k in ('matchete_mass_basis_verified','physical_matching_performed',
              'matchete_rgbeta_component_map_certified','coefficient_map_certified'):
        require(param.get(k) is False,model+': parameter physical flag '+k)
    require(mat.get('format')==NATIVE_FORMAT and mat.get('status')==NATIVE_STATUS and
            mat.get('model')==model,model+': native Matchete mass-basis certificate missing')
    for k in ('mass_basis_ordinary_check_passed','heavy_mass_basis_detailed_check_passed',
              'all_old_scalar_fields_removed','all_17_quartics_present',
              'both_yukawa_interactions_present','matchete_rotation_uses_exact_symbolic_spurions'):
        require(mat.get(k) is True,model+': missing Matchete native check '+k)
    for k in ('physical_matching_performed','matching_to_two_scalar_thresholds_performed',
              'component_basis_map_certified','matchete_rgbeta_quartic_normalization_certified'):
        require(mat.get(k) is False,model+': unsafe native matching claim '+k)
    require(mat.get('diagonal_heavy_fields')==['NewScalarX','NewScalarY'] and
            mat.get('source_baseline_interactions')==11 and mat.get('source_extra_quartics')==7 and
            mat.get('source_new_interaction_count')==18,model+': native eigenfield inventory invalid')
    interactions=mat.get('source_interactions')
    require(isinstance(interactions,list) and len(interactions)==18 and
            len({x.get('name') for x in interactions})==18 and
            sum(x.get('class')=='AdditionalQuartic' for x in interactions)==7,model+': source interactions incomplete')
    A=cm(param['original_mass_squared_A_gev2'],2,'native A')
    U=cm(param['unitary_U'],2,'native U')
    source_basis=mass.get('basis',{})
    As=cm(source_basis['mass_squared_flavour_matrix_gev2'],2,'saved A')
    Us=cm(source_basis['unitary_flavour_rotation_U'],2,'saved U')
    e=np.asarray(param.get('heavy_mass_squared_gev2',[]),dtype=float)
    es=np.asarray(source_basis.get('masses_squared_gev2',[]),dtype=float)
    require(e.shape==(2,) and es.shape==(2,) and np.isfinite(e).all() and
            np.all(e>0) and e[0]<e[1],model+': eigenvalues malformed')
    scale=max(1.,np.max(abs(A)),np.max(abs(e)))
    require(np.linalg.norm(A-As)/scale<2e-11 and np.linalg.norm(U-Us)<2e-11 and
            np.linalg.norm(e-es)/scale<2e-11,model+': native parameter matrix not source mass basis')
    unitarity=float(np.linalg.norm(U.conj().T@U-np.eye(2)))
    diagonal=float(np.linalg.norm(U.conj().T@A@U-np.diag(e))/scale)
    require(unitarity<3e-11 and diagonal<3e-11,model+': native mass diagonalization failed')
    spurions=mat.get('rotation_spurion_values_from_numerical_U')
    require(isinstance(spurions,dict) and len(spurions)==4,model+': native U spurions missing')
    n_spur=np.array([[cp(spurions['T3MEU11'],'U11'),cp(spurions['T3MEU12'],'U12')],
                      [cp(spurions['T3MEU21'],'U21'),cp(spurions['T3MEU22'],'U22')]])
    require(np.linalg.norm(n_spur-U)<3e-11,model+': native Matchete U spurions disagree')
    for k in ('epsilon_bridge_BB','epsilon_bridge_UU'):
        d=mat.get(k,{})
        require(d.get('native_intertwiner_inputform')=='{{0, 1}, {-1, 0}}' and
                d.get('phase_relative_to_matchete_component_epsilon_verified') is True and
                d.get('physical_RGBeta_component_map_certified') is False,model+': native epsilon convention not verified '+k)
    for k,rnk in (('native_H2_portal_epsilon_tensors',2),('native_cross_scalar_epsilon_tensors',2),('native_S2_self_epsilon_tensor',1)):
        d=mat.get(k,{})
        require(d.get('ExactTensorProjectionVerified') is True and
                d.get('OldInvariantCount')==d.get('NewInvariantCount')==rnk and
                d.get('TransportedRank')==d.get('CombinedBasisRank')==rnk and
                d.get('PhysicalCouplingMapCertified') is False,model+': native tensor rank failure '+k)
    extra=mat.get('native_additional_quartic_tensors',{})
    certs=extra.get('Certificates')
    rows=v21.get('models',{}).get(model,{}).get('native_CG_sectors')
    require(isinstance(rows,list) and len(rows)==8 and rows[0].get('source_matchete_coupling')=='lambdaT3',model+': v21 CG sectors incomplete')
    source_names={x['source_matchete_coupling'] for x in rows[1:]}
    require(isinstance(certs,dict) and len(certs)==7 and
            extra.get('SourceOperatorCount')==extra.get('RecoupledOperatorCount')==7 and
            extra.get('AllProjectedTensorsInvariant') is True and
            extra.get('PhysicalCouplingMapCertified') is False,model+': native seven CG certification incomplete')
    native_names={c.get('SourceCoupling') for c in certs.values()}
    require(native_names==source_names and len(native_names)==7,model+': native seven CG symbol set differs from v21')
    for name,info in certs.items():
        require(info.get('ExactTensorProjectionVerified') is True and
                info.get('SourceTensorIndex') in (1,2) and
                info.get('SourceCoupling') in source_names and
                info.get('NativeInvariantSpanRank')==info.get('CombinedBasisRank')==2,model+': native extra CG projection failed '+name)
    for name,expected in (('native_yukawa1_tensor','1'),('native_yukawa2_epsilon_tensor','-1' if model=='B_-1' else '1')):
        d=mat.get(name,{})
        key='OldToNewCGFactor' if name=='native_yukawa1_tensor' else 'OldTimesEpsilonCGFactor'
        require(d.get(key)==expected and d.get('ExactTensorProjectionVerified') is True and
                d.get('FermionSU2Dimension')==(1 if model=='B_-1' else 3) and
                d.get('PhysicalCouplingMapCertified') is False,model+': native Yukawa epsilon failure '+name)
    require(param.get('epsilon')==EPS.astype(int).tolist(),model+': source epsilon changed')
    plan=param.get('threshold_plan',{})
    masses=np.sqrt(e)
    require(plan.get('decreasing_scale_order')==['NewScalarY','NewScalarX'] and
            plan.get('sequential_matching_performed') is False and
            plan.get('existing_eft1_endpoint_below_both_mass_thresholds') is True,model+': threshold-order evidence changed')
    scale_eft=float(param.get('source_eft1_scale_gev'))
    require(np.isfinite(scale_eft) and 0<scale_eft<masses[0],model+': endpoint not below both thresholds')
    return {'model':model,'mass_squared_gev2':e.tolist(),'mass_gev':masses.tolist(),
            'source_eft1_endpoint_gev':scale_eft,'unitarity_residual':unitarity,
            'mass_diagonalization_relative_residual':diagonal,
            'native_matchete_heavy_mass_basis_verified':True,
            'native_matchete_epsilon_and_seven_extra_invariants_verified':True,
            'source_rotated_16_complex_CGs_in_v21':8,
            'native_cg_source_symbol_reconciliation':sorted(source_names),
            'native_yukawa2_epsilon_factor':-1 if model=='B_-1' else 1,
            'two_scalar_thresholds_above_saved_eft1_endpoint':True,
            'source_sha256':{'mass_basis':digest(mass_path),'native_parameters':digest(param_path),
                             'native_matchete_eigenbasis':digest(mat_path)},
            'physical_fermion_first_scalar_matching_performed':False}

def audit(root:Path,v21_path:Path,v22_path:Path)->dict:
    root=Path(root).resolve();v21_path=Path(v21_path);v22_path=Path(v22_path)
    v21,v22=recheck_v21_v22(root,v21_path,v22_path)
    models={m:reconcile_model(root,m,v21) for m in BC}
    return {'format':FORMAT,'status':STATUS,'models':models,'models_checked':2,
            'native_17_quartic_eigenfield_lagrangians_checked':2,
            'mass_eigen_source_native_CG_concordance':16,
            'source_sha256':{'v21':digest(v21_path),'v22':digest(v22_path)},
            'certification':{'native_Matchete_eigenfield_HeavyMassBasis_passed_for_BC':True,
              'native_Matchete_complex_CG_transport_concordant_with_v21':True,
              'physical_Matchete_RGBeta_coupling_dictionary':False,
              'physical_fermion_first_scalar_threshold_matching':False,
              'signed_GammaCC_to_Weyl':False,'one_loop_dim5_ADM':False,
              'physical_Wilson_evolution_authorized':False,
              'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT'}}
