"""WS4: independent, read-only, full-mass direct-UV p=0 amplitude reference.

No call to core.evaluate, core.takagi, or core.scalar_real_basis is made in
independent_reference().  Heavy propagation is evaluated using the Hermitian
matrix function M^dagger f(M M^dagger), rather than a Takagi diagonalisation.
The scalar Hessian is reconstructed directly from V2 in Cartesian fields.

The reference deliberately DOES NOT certify the provisional input's native
Matchete/RGBeta-to-Weyl coefficient normalisation or sign.  It tests the
physics of the exact loop given those inputs.  Optional comparison with the
production one-loop engine is kept as a separate check.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .native_amplitude_contract import (ContractError, model_dictionary,
                                        source_audit, audit_frozen_inputs)

LOOP = 1.0 / (16.0 * np.pi**2)
FORMAT = 't3-ws4-independent-uv-spectral-reference-v1'
PROVISIONAL = Path('output/general_observable/lo_direct_one_loop/provisional_inputs')
EVIDENCE = {
    'spinor': Path('Reports/output/workstream_a/yukawa_spinor_phase/a0_spinor_phase.json'),
    'yukawa_normalization': Path('Reports/output/workstream_a/yukawa_normalization_bridge/yukawa_normalization_bridge.json'),
    'a0_action': Path('Reports/output/workstream_a/yukawa_a0_action_bridge/a0_action_bridge.json'),
}


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ContractError(reason)


def decode(value: Any, label: str, shape: tuple[int, ...] | None = None) -> np.ndarray:
    if isinstance(value, Mapping) and set(value) == {'real', 'imag'}:
        a = np.asarray(value['real'], float) + 1j * np.asarray(value['imag'], float)
    else:
        a = np.asarray(value, complex)
    if shape is not None:
        require(a.shape == shape, f'{label} shape {a.shape}, expected {shape}')
    require(bool(np.isfinite(a).all()), f'{label}: nonfinite element')
    return np.asarray(a, dtype=np.complex128)


def encode(z: Any) -> dict[str, Any]:
    a = np.asarray(z, complex)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def cartesian_hessian(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """V2=phi^dag A phi + Re(phi^T B phi), phi=(x+i y)/sqrt(2).

    This is a direct real second derivative of the potential, not a call to
    the original scalar mass-eigenstate routine.
    """
    n = a.shape[0]
    require(a.shape == b.shape == (n, n), 'Invalid complex scalar quadratic dimensions')
    unit = max(np.linalg.norm(a), np.linalg.norm(b), 1.)
    require(np.linalg.norm(a-a.conj().T) <= 2e-12*unit, 'A not Hermitian')
    require(np.linalg.norm(b-b.T) <= 2e-12*unit, 'B not symmetric')
    h = np.block([[a.real+b.real, -a.imag-b.imag],
                  [a.imag-b.imag, a.real-b.real]])
    require(np.linalg.norm(h-h.T) <= 2e-12*unit, 'Real Hessian not symmetric')
    return (h+h.T)/2


def potential(a: np.ndarray, b: np.ndarray, phi: np.ndarray) -> float:
    return float(np.real(np.vdot(phi, a@phi) + phi.T@b@phi))


def ratio_function(x: np.ndarray) -> np.ndarray:
    """g(x) = x*log(x)/(x-1), analytic through x=1.

    The exact MSbar p=0 B0 at positive masses is
        B0(a,b;mu) = 1-log(a/mu^2)-g(b/a).
    Only differences in scalar masses are used below, after independently
    checking the UV-pole coefficient vanishes. Thus the leading constant
    and log(mu) are cancelled before floating-point summation.
    """
    x = np.asarray(x, float)
    require(bool(np.all(x>0)) and bool(np.isfinite(x).all()),
            'Nonpositive/nonfinite loop squared-mass ratio')
    t = x-1
    out = np.empty_like(x)
    small = np.abs(t)<1e-4
    if np.any(small):
        v=t[small]
        out[small]=1+v/2-v*v/6+v**3/12-v**4/20+v**5/30-v**6/42
    if np.any(~small):
        z=x[~small]
        out[~small]=z*np.log(z)/(z-1)
    return out


def b0_exact(a: float, b: float, mu: float) -> float:
    """Finite MSbar B0, independently encoded via the exact log primitive."""
    require(min(a,b,mu)>0, 'B0 positive mass condition')
    return float(1-np.log(a/mu**2)-ratio_function(np.array([b/a]))[0])


def _vertex_set(data: Mapping[str, Any], nf: int, ns: int):
    vertices = data['vertices']
    require(isinstance(vertices, list) and vertices, 'Missing UV vertices')
    names=[]; parsed=[]
    for row in vertices:
        name=str(row['name'])
        require(name not in names, 'Duplicate UV source vertex: '+name)
        names.append(name)
        y=decode(row.get('field', np.zeros((3,nf,ns))), name+'.field', (3,nf,ns))
        z=decode(row.get('conjugate', np.zeros((3,nf,ns))), name+'.conjugate', (3,nf,ns))
        parsed.append((y,z))
    return names, parsed


def independent_reference(data: Mapping[str, Any], *, pole_limit: float = 1e-9,
                          check_rephase: bool = False) -> dict[str, Any]:
    """Exact full-propagator UV one-loop, with no heavy-mass expansion.

    A finite radiative mass requires the complete supplied vertex set to
    have a vanishing pole.  If not, no arbitrary finite subtraction is made.
    """
    M=decode(data['fermion_mass_gev'], 'M')
    n=M.shape[0]
    require(M.shape==(n,n), 'Fermion mass not square')
    require(np.linalg.norm(M-M.T)<2e-12*max(np.linalg.norm(M),1.),
            'Fermion mass is not complex symmetric')
    A=decode(data['scalar_A_gev2'], 'A')
    ns=A.shape[0]
    B=decode(data['scalar_B_gev2'],'B',(ns,ns))
    H=cartesian_hessian(A,B)
    scalar_m2,R=np.linalg.eigh(H)
    require(np.min(scalar_m2)>0, 'Scalar spectrum nonpositive: no p=0 finite check')
    C=(R[:ns,:]+1j*R[ns:,:])/np.sqrt(2)
    # Hermitian fermion square; not Takagi.  Its eigenvectors compute an exact
    # mass-dependent propagator M^dagger f(M M^dagger).
    F2=M@M.conj().T
    fm2,Q=np.linalg.eigh((F2+F2.conj().T)*0.5)
    require(np.min(fm2)>0 and np.min(fm2)>1e-15*np.max(fm2),
            'Heavy fermion spectrum singular / numerically unresolved')
    names, yz=_vertex_set(data,n,ns)
    G=[np.einsum('ais,sr->air',y,C,optimize=True)+
       np.einsum('ais,sr->air',z,C.conj(),optimize=True) for y,z in yz]
    # The pole is independent of all scalar eigenmasses; test it before
    # subtracting B0's scalar-independent part.
    pole=np.zeros((3,3),complex)
    natural=0.
    pair_poles={}
    for i,name in enumerate(names):
        for j in range(i,len(names)):
            p=np.zeros((3,3),complex)
            for r in range(2*ns):
                t=G[i][:,:,r]@M.conj().T@G[j][:,:,r].T
                if i!=j:t=t+t.T
                p+=t
                natural+=(np.linalg.norm(G[i][:,:,r])*np.linalg.norm(G[j][:,:,r])
                          *np.linalg.norm(M,ord=2)*(2 if i!=j else 1))
            pole+=p
            pair_poles[f'{name} x {names[j]}']=float(np.linalg.norm(p)*LOOP)
    pole_relative=float(np.linalg.norm(pole)/max(natural,1e-100))
    require(pole_relative < pole_limit,
            f'Nonzero MSbar UV pole relative to natural amplitude: {pole_relative:.3e}')

    reference_m2=float(np.sqrt(np.min(scalar_m2)*np.max(scalar_m2)))
    # Matrix functions diagonal in the Hermitian eigenbasis.  Unlike the
    # production Takagi sum, scalar-independent terms are gone before summing.
    spectral_matrices=[]
    for s2 in scalar_m2:
        gdiff=-(ratio_function(np.array([float(s2)/v for v in fm2]))-
                ratio_function(np.array([reference_m2/v for v in fm2])))
        f=(Q*gdiff[None,:])@Q.conj().T
        spectral_matrices.append(M.conj().T@f)
    pair_matrices={}; individual=[]; total=np.zeros((3,3),complex)
    abs_sum=0.
    for i,name in enumerate(names):
        for j in range(i,len(names)):
            mat=np.zeros((3,3),complex)
            for r,S in enumerate(spectral_matrices):
                term=(G[i][:,:,r]@S@G[j][:,:,r].T)*LOOP
                if i!=j:term=term+term.T
                mat+=term
                abs_sum+=np.linalg.norm(term)
            key=f'{name} x {names[j]}'
            pair_matrices[key]=mat
            total+=mat
    full_norm=float(np.linalg.norm(total))
    asym=float(np.linalg.norm(total-total.T)/max(full_norm,1e-100))
    require(asym<5e-9, 'Independent UV mass is not complex symmetric')
    masses=np.sort(np.linalg.svd(total,compute_uv=False)*1e9)
    scale_mu=float(data.get('renormalization_scale_gev',1e7))
    result={
        'status':'independent_exact_spectral_reference_from_conditional_input',
        'method':'M_dagger f(M M_dagger), real Cartesian scalar Hessian, pole-subtracted B0 differences',
        'subtraction':'MSbar p=0; subtraction only after explicit zero-pole check',
        'input_scale_gev':scale_mu,
        'fermion_mass_singular_values_gev':np.sort(np.sqrt(fm2)).tolist(),
        'scalar_masses_gev':np.sqrt(scalar_m2).tolist(),
        'complex_symmetric_matrix_gev':encode(total),
        'takagi_masses_ev_from_independent_singular_values':masses.tolist(),
        'source_vertex_pair_matrices_gev':{k:encode(v) for k,v in pair_matrices.items()},
        'source_vertex_pair_pole_norms_gev':pair_poles,
        'diagnostics':{
            'uv_pole_relative_to_natural_amplitude':pole_relative,
            'mass_symmetry_relative':asym,
            'sum_of_absolute_component_norms_gev':float(abs_sum),
            'cancellation_factor_absolute_over_total':float(abs_sum/max(full_norm,1e-100)),
            'input_has_scalar_B':bool(np.linalg.norm(B)>0),
            'scalar_cartesian_hessian_symmetry':float(np.linalg.norm(H-H.T)/max(np.linalg.norm(H),1.)),
            'matrix_full_norm_gev':full_norm,
        },
        'physical_source_normalisation_certified':False,
    }
    if check_rephase:
        phases_f=np.exp(1j*np.arange(n)*0.371)
        phases_s=np.exp(1j*np.arange(ns)*0.239)
        transformed=rephase_input(data,phases_f,phases_s)
        other=independent_reference(transformed,pole_limit=pole_limit,check_rephase=False)
        new=decode(other['complex_symmetric_matrix_gev'],'rephased mass',(3,3))
        diff=float(np.linalg.norm(new-total)/max(full_norm,1e-100))
        result['diagnostics']['internal_field_rephasing_relative']=diff
        # Avoid declaring numerical rephasing a failure if input is so
        # cancellation-conditioned that round-off is not under control.
        result['diagnostics']['rephasing_natural_error']=float(
            np.linalg.norm(new-total)/max(abs_sum,1e-100))
    return result


def rephase_input(data: Mapping[str, Any], fermion_phases: np.ndarray,
                  scalar_phases: np.ndarray) -> dict[str, Any]:
    """Passive field replacement psi=Df psi', phi=Ds phi' (all complex)."""
    m=decode(data['fermion_mass_gev'],'M')
    a=decode(data['scalar_A_gev2'],'A')
    b=decode(data['scalar_B_gev2'],'B')
    nf=len(m);ns=len(a)
    require(fermion_phases.shape==(nf,) and scalar_phases.shape==(ns,),
            'Field rephasing dimension mismatch')
    df=np.diag(fermion_phases);ds=np.diag(scalar_phases)
    changed=dict(data)
    changed['fermion_mass_gev']=encode(df.T@m@df)
    changed['scalar_A_gev2']=encode(ds.conj().T@a@ds)
    changed['scalar_B_gev2']=encode(ds.T@b@ds)
    out=[]
    for row in data['vertices']:
        item={'name':row['name']}
        if 'field' in row:
            y=decode(row['field'],row['name']+'.field',(3,nf,ns))
            item['field']=encode(y*fermion_phases[None,:,None]*scalar_phases[None,None,:])
        if 'conjugate' in row:
            z=decode(row['conjugate'],row['name']+'.conjugate',(3,nf,ns))
            item['conjugate']=encode(z*fermion_phases[None,:,None]*scalar_phases.conj()[None,None,:])
        out.append(item)
    changed['vertices']=out
    return changed


def native_portal_symbol_trace(root: Path) -> dict[str, Any]:
    """Source observation: native scalar potential generates indexed Inv names.

    It is NOT legitimate to set lambdaH1==lambdaH1Inv1 solely because their
    names look related; no signed/magnitude conversion is derived here.
    """
    from .native_amplitude_contract import _wl_without_comments
    rels={
        'scalar':'Lagrangian/interactions/ScalarPotential.wl',
        'topology':'Lagrangian/interactions/T3Topology.wl',
        'fields':'Lagrangian/model/T3Fields.wl',
    }
    source={}
    for key,rel in rels.items():
        path=root/rel
        require(path.is_file(),f'Missing {rel}')
        source[key]=_wl_without_comments(path.read_text(encoding='utf-8-sig'))
    scalar=source['scalar'];fields=source['fields'];topology=source['topology']
    required={
        'potential_builds_invariant_specific_portal_family':
            'BuildHiggsPortalInvariantCandidates' in scalar and
            'baseName="lambdaH"<>ToString[which]' in scalar,
        'portal_coupling_symbol_dynamically_has_Inv_index':
            'couplingName=Symbol[couplingBase<>"Inv"<>ToString[index]]' in scalar,
        'topology_calls_general_scalar_potential_family':
            'BuildGeneralScalarPotentialCandidates[model]' in topology,
        'unindexed_lambdaH1_defined':
            'DefineCoupling[lambdaH1,SelfConjugate->True]' in fields,
        'unindexed_lambdaH2_defined':
            'DefineCoupling[lambdaH2,SelfConjugate->True]' in fields,
    }
    require(all(required.values()),'Current native portal symbol source forms changed: '+str(required))
    return {
        'current_source_forms':required,
        'source_level_bare_to_indexed_portal_map_established':False,
        'observed': 'GeneralScalarPotential uses lambdaH1InvN/lambdaH2InvN; fitted input builder reads lambdaH1/lambdaH2',
        'not_inferred':'No signed/magnitude equivalence between bare and invariant-indexed parameters is assumed',
        'next_native_measurement':'Export the scalar quadratic operator with one isolated native lambdaHiInvN, the physical generator basis and H0=v/sqrt(2), then project to A; compare to the benchmark lambdaHi definition',
    }


def optional_native_evidence(root: Path) -> dict[str, Any]:
    results={}
    for name,rel in EVIDENCE.items():
        path=root/rel
        if path.exists():
            raw=path.read_bytes();rec=json.loads(raw.decode('utf-8-sig'))
            results[name]={'available':True,'path':str(rel),
                           'sha256':sha256(raw).hexdigest(),
                           'declared_status':rec.get('status'),
                           'native_physical_cross_package_dictionary_certified':False}
        else:
            results[name]={'available':False,'path':str(rel)}
    return results


def native_operator_export_evidence(root: Path, path: Path | None) -> dict[str, Any]:
    """Verify provenance of the OPTIONAL live Matchete restricted LInt export.

    This certifies only that Matchete assembled the operators in its own
    symbolic action.  It does not set numerical 4-spinor or scalar CG phases.
    """
    if path is None:
        return {'available':False,'physical_component_normalisation_certified':False,
                'reason':'Optional live Matchete operator inventory not supplied'}
    path=path.resolve()
    require(path.is_file(),'Missing requested native Matchete export '+str(path))
    raw=path.read_bytes()
    data=json.loads(raw.decode('utf-8-sig'))
    require(data.get('format')=='t3-ws4-native-matchete-restricted-LInt-operator-export-v1'
            and data.get('status')=='native_restricted_operator_inventory_extracted_not_component_physical',
            'Live Matchete operator export schema/status changed')
    models=data.get('models',[])
    require([m.get('model') for m in models]==['A_+0','B_-1'],
            'Live Matchete operator witness did not inspect the two specified models')
    for m in models:
        f=m.get('native_restricted_LInt_function_names',{})
        require(m.get('native_full_lagrangian_validation') is True and
                all(f.get(x) is True for x in ('topology_y1_in_LInt',
                    'topology_y2_in_LInt','T3_mixing_in_LInt')),
                'Live restricted native topology operator inventory incomplete')
    paths={'LagrangianBuilder':'Lagrangian/model/LagrangianBuilder.wl',
           'ScalarPotential':'Lagrangian/interactions/ScalarPotential.wl',
           'T3Topology':'Lagrangian/interactions/T3Topology.wl'}
    for name,rel in paths.items():
        actual=int(sha256((root/rel).read_bytes()).hexdigest(),16)
        require(int(data.get('source_files',{}).get(name,-1))==actual,
                'Native operator export is stale against local '+rel)
    return {'available':True,'path':str(path),'sha256':sha256(raw).hexdigest(),
            'native_Matchete_restricted_LInt_operators_verified':True,
            'model_native_symbol_profiles':{m['model']:m['native_restricted_LInt_function_names']
                                            for m in models},
            'physical_component_normalisation_certified':False,
            'unresolved_GammaCC_and_RGBeta_phase':True}


def run_all(root: Path, *, rephase: bool = False, compare_production: bool = True,
            input_dir: Path | None = None,
            native_export: Path | None = None) -> dict[str, Any]:
    """Read-only source-guarded all-16 batch; caller owns report file writes."""
    root=root.resolve()
    audit=source_audit(root)
    require(audit['status']=='source_forms_verified',
            'Source forms changed: '+ '; '.join(audit['errors']))
    symbols=native_portal_symbol_trace(root)
    native_evidence=native_operator_export_evidence(root,native_export)
    restricted_input_audit=(audit_frozen_inputs(root) if input_dir is None else
                            {'status':'custom_input_dir_not_frozen_input_sha_audited'})
    directory=(root/PROVISIONAL) if input_dir is None else input_dir.resolve()
    results=[]; output={}
    for row in model_dictionary():
        label=row['model'];path=directory/(label+'.json')
        require(path.is_file(),'Missing provisional input '+str(path))
        payload=path.read_bytes();d=json.loads(payload.decode('utf-8-sig'))
        require(d.get('model')==label, 'Signed model label differs: '+label)
        result=independent_reference(d,check_rephase=rephase)
        result['model']=label
        result['input_file']=str(path.relative_to(root)) if path.is_relative_to(root) else str(path)
        result['input_sha256']=sha256(payload).hexdigest()
        rowout={'model':label, 'independent_masses_ev':result['takagi_masses_ev_from_independent_singular_values'],
                'uv_pole_relative':result['diagnostics']['uv_pole_relative_to_natural_amplitude'],
                'cancellation_factor':result['diagnostics']['cancellation_factor_absolute_over_total'],
                'physical_source_normalisation_certified':False,
                'status':'independent_conditional_UV_reference'}
        if compare_production:
            from Numerical.direct_one_loop.core import evaluate
            old=evaluate(d,retain_individual_terms=False)
            oldM=decode(old['total']['complex_symmetric_matrix_gev'],'production Mnu',(3,3))
            newM=decode(result['complex_symmetric_matrix_gev'],'independent Mnu',(3,3))
            absdelta=float(np.linalg.norm(oldM-newM))
            cond=float(result['diagnostics']['cancellation_factor_absolute_over_total'])
            pairerrs={}
            for pair,z in result['source_vertex_pair_matrices_gev'].items():
                k=decode(z,'independent pair',(3,3))
                other=decode(old['vertex_pair_matrices_gev'][pair],'production pair',(3,3))
                pairerrs[pair]=float(np.linalg.norm(k-other)/max(np.linalg.norm(k),np.linalg.norm(other),1e-100))
            result['production_comparison']={
                'complex_mass_absolute_delta_gev':absdelta,
                'complex_mass_relative_delta':float(absdelta/max(np.linalg.norm(newM),np.linalg.norm(oldM),1e-100)),
                'relative_delta_over_uncancelled_natural_amplitude':float(absdelta/max(result['diagnostics']['sum_of_absolute_component_norms_gev'],1e-100)),
                'max_vertex_pair_relative_delta':max(pairerrs.values(),default=0.),
                'per_vertex_pair_relative_delta':pairerrs,
                'production_neutrino_masses_ev':old['total']['neutrino_masses_ev'],
                'possible_large_roundoff_from_cancellations':bool(cond>1e9),
                'source_physical_certification':False,
            }
            rowout['relative_production_difference']=result['production_comparison']['complex_mass_relative_delta']
        output[label]=result;results.append(rowout)
    return {'format':FORMAT,
            'status':'all_16_independent_conditionally_normalised_uv_references_computed',
            'source_guards':audit, 'native_portal_symbol_trace':symbols,
            'optional_native_evidence':optional_native_evidence(root),
            'live_native_operator_export':native_evidence,
            'restricted_input_structure_audit':restricted_input_audit,
            'model_count':len(results), 'models':results,
            'reports_by_model':output,
            'physical_source_certified_models':0,
            'scale_status':'provisional inputs at mu0, not run to muF',
            'EFT1_tree_interface':'not compared at common scale; formal Gaussian bridge remains separate',
            'hard_one_loop_fermion_threshold_matching':'not implemented',
            'scalar_to_SMEFT_matching':'not implemented',
            'extra_symmetry_allowed_yukawas':'six permitted; set zero by frozen restricted benchmark, not symmetry',
            'unresolved_native_factors':[
                'signed MF and 4-spinor GammaCC map to the same Weyl basis as both Yukawas',
                'RGBeta epsilon/triplet component axes and normalisation against Matchete',
                'bare lambdaH1/lambdaH2 to native lambdaHiInvN scalar-potential identification',
                'source lambdaT3 contraction plus Higgs VEV and potential-to-Hessian sign',
            ]}
