"""Six T3 self-bilinears: exact normal-propagator resummation into existing C5.

Uses the *saved* continuous_v1 scalar/UV-EFT1 state; does NOT repeat RGEs.
The self-conjugate neutral scalar's full quadratic propagator is decomposed
into two real eigenmultiplets. The original one-loop topology, linear in
one NORMAL propagator of that species, is evaluated at both eigenmasses
with exact weights 1/2, retaining the evolved heavy Yukawas and low-energy
SM+Weinberg continuation.

It is not the COMPLETE full-general-scalar C5: anomalous scalar propagators
combined with additional independent HH scalar vertices, and the associated
native Matchete/RGBeta CG normalisations, are NOT evaluated. These must be
separately matched, not inferred from the existing topology kernel.

    python -m scripts.RunT3SelfBilinearNeutrino --all --run-id continuous_v1
    python -m scripts.RunT3SelfBilinearNeutrino --all --run-id continuous_v1 --resume
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from hashlib import sha256
from pathlib import Path
from typing import Any
import numpy as np

from Numerical.experimental.T3SelfBilinearMassBasis import (
    SelfBilinearMassBasis, native_tensors_for_model, evaluate_tensor,
)
from scripts import RunT3FullScalarNeutrino as neutrino
from scripts import RunT3FullScalarContinuous as continuous

ROOT=continuous.ROOT
SELF=('A_+0','A_-2','D_+0','D_-2','E_+0','E_-2')
FORMAT='t3-self-bilinear-normal-propagator-one-loop-v1'
STATUS='CandidateSelfBilinearNormalPropagatorResummedC5AnomalousNotMatched'
MAT_TENSORS=ROOT/'output/general_t3_stage22/matchete_tensors.json'


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError('T3 self-bilinear matching: '+message)


def split_species(model: str) -> str:
    require(model in SELF,'unsupported self-bilinear T3 model')
    return 'S1' if model.endswith('+0') else 'S2'


def candidate_basis(model: str, endpoint: dict) -> SelfBilinearMassBasis:
    """Construct the actual 1/3-dimensional neutral real mass basis."""
    species=split_species(model)
    rep=endpoint['eft1_scalar_endpoint_state']['representation']
    scalar=endpoint['eft1_scalar_endpoint_state']
    dim=rep['d_s1' if species=='S1' else 'd_s2']
    mass=scalar['mS1Sq' if species=='S1' else 'mS2Sq']
    mu=neutrino.complex_pair(endpoint['quadratic_mu2_at_scalar_endpoint_gev2'],'mu2')
    basis=SelfBilinearMassBasis.build(species,dim,mass,mu)
    require(not basis.tachyonic,'nonpositive split scalar mass squared')
    require(endpoint['quadratic_kind']==('self_s1' if species=='S1' else 'self_s2'),
            'quadratic registry does not match the selected self-scalar')
    recorded=np.asarray(endpoint['scalar_endpoint_spectrum']['mass_squared_eigenvalues_gev2'],dtype=float)
    require(np.isfinite(recorded).all(),'nonfinite saved mass spectrum')
    # Both real eigenmultiplets and the unsplit second scalar must appear,
    # with multiplicities matching the full Hessian (2 real per complex component).
    other_dim=rep['d_s2' if species=='S1' else 'd_s1']
    other_sq=float(scalar['mS2Sq' if species=='S1' else 'mS1Sq'])
    expected=np.sort(np.r_[np.repeat(basis.eigen_masses_squared[0],dim),
                              np.repeat(basis.eigen_masses_squared[1],dim),
                              np.repeat(other_sq,2*other_dim)])
    require(recorded.shape==expected.shape and np.allclose(recorded,expected,rtol=2e-10,atol=2e-6),
            'saved 16-model Hessian disagrees with proposed real scalar mass basis')
    return basis


def normal_propagator_c5(kernel: Any, uv: Any, eft: Any,
                         basis: SelfBilinearMassBasis) -> tuple[np.ndarray,np.ndarray,list[np.ndarray]]:
    """Exact normal-propagator channel of the existing one-scalar-line kernel.

    C5_normal = [C5(m_+²) + C5(m_-²)]/2, where the original kernel
    contains one normal propagator for the chosen scalar representation.
    No naive m² -> mean(m²) replacement is performed.
    """
    attr='mS1Sq' if basis.species=='S1' else 'mS2Sq'
    base=neutrino.validated_c5(kernel.evaluate(uv,eft))
    pieces=[]
    for mass_sq in basis.eigen_masses_squared:
        state=replace(eft,**{attr:float(mass_sq)}).validated()
        pieces.append(neutrino.validated_c5(kernel.evaluate(uv,state)))
    resummed=neutrino.validated_c5(0.5*(pieces[0]+pieces[1]))
    return base,resummed,pieces


def _rel(a: np.ndarray,b:np.ndarray) -> float | None:
    den=float(np.linalg.norm(b))
    return float(np.linalg.norm(a-b)/den) if den>0 else None


def _hash(path: Path) -> str:
    require(path.is_file(),'missing required source '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def native_vertex_transport(model: str, basis: SelfBilinearMassBasis) -> dict:
    """Rotate all 5/6 saved extra Matchete quartics as candidate tensors.

    The native-Matchete-to-spherical component identification is still only
    a declared hypothesis. No vertex is inserted into C5 by this operation.
    """
    vertices=native_tensors_for_model(MAT_TENSORS,model,basis)
    records=[]
    rng=np.random.default_rng(9121)
    for item in vertices:
        species=item['species']
        tensor=item['tensor']
        rotated=item['rotated']
        # This confirms field-polynomial transport, not native gauge-CG parity.
        worst=0.0
        for _ in range(3):
            real=rng.normal(size=basis.dim)
            imaginary=rng.normal(size=basis.dim)
            original_self=basis.fields(real,imaginary)
            independent={}
            for species_name,dim in zip(species,tensor.shape):
                core=species_name.removesuffix('dag')
                if core!=basis.species and core not in independent:
                    independent[core]=rng.normal(size=dim)+1j*rng.normal(size=dim)
            old=[]
            for name in species:
                if name.removesuffix('dag')==basis.species:
                    z=original_self
                else:
                    z=independent[name.removesuffix('dag')]
                old.append(z.conjugate() if name.endswith('dag') else z)
            new=dict(independent)
            new[basis.species+'_R']=real
            new[basis.species+'_I']=imaginary
            for core,z in independent.items():
                new[core+'dag']=z.conjugate()
            before=evaluate_tensor(tensor,old)
            after=sum(evaluate_tensor(val,[new[name] for name in names])
                      for names,val in rotated.items())
            err=abs(before-after)/max(1.,abs(before),abs(after))
            worst=max(worst,float(err))
        require(worst<3e-12,model+': rotated quartic tensor fails field-polynomial equality')
        records.append({'matchete_symbol':item['symbol'],
           'source_species':list(species),
           'complex_coupling_plus_hermitian_conjugate':'source Matchete potential; conjugate implicit',
           'component_axes_dimensions':list(tensor.shape),
           'rotation_identity_max_relative_residual':worst,
           'rotated_terms':[{'real_multiplet_fields':list(names),
                             'coeff_tensor_real':np.asarray(val).real.tolist(),
                             'coeff_tensor_imag':np.asarray(val).imag.tolist()}
                            for names,val in sorted(rotated.items())]})
    return {'format':'t3-self-bilinear-rotated-native-extra-vertices-v1',
            'model':model,'source_matchete_tensor_sha256':_hash(MAT_TENSORS),
            'native_matchete_tensor_is_spherical_basis_certified':False,
            'component_map_status':'CandidateMatcheteAxisEqualsDeclaredSphericalAxisHypothesis',
            'extra_quartic_vertices_rotated':len(records),
            'all_vertices_included_in_c5':False,
            'vertices':records}


def calculate(model: str, *, source_root: Path, output_root: Path,
              control_output_id: str='neutrino_v1', resume: bool=False) -> dict:
    require(model in SELF,'only the six A/D/E holomorphic self-bilinears are supported')
    tag=continuous.model_tag(model)
    folder=source_root/tag
    record=folder/'continuous_trajectory.json'
    source=neutrino.verify_continuous(model,record)
    basis=candidate_basis(model,source)
    uv=neutrino.decode_uv(source['uv_fermion_threshold_state'])
    eft=neutrino.decode_eft(source['eft1_scalar_endpoint_state'])
    from scripts.RunT3GeneralScalar48Campaign import matching
    kernel,kernel_path,kind=matching(model)
    upstream_hash=_hash(record)
    kernel_hash=_hash(kernel_path)
    script_hash=_hash(Path(__file__).resolve())
    native_tensor_hash=_hash(MAT_TENSORS)
    existing_control=source_root/control_output_id/(tag+'_neutrino.json')
    control_hash=_hash(existing_control) if existing_control.is_file() else None
    target=output_root/(tag+'_self_bilinear_c5.json')
    if target.exists():
        data=neutrino.load(target)
        require(resume and data.get('format')==FORMAT and data.get('status')==STATUS and
                data.get('source_continuous_sha256')==upstream_hash and
                data.get('source_matching_kernel_sha256')==kernel_hash and
                data.get('source_runner_sha256')==script_hash and
                data.get('source_native_tensor_sha256')==native_tensor_hash and
                data.get('source_control_sha256')==control_hash,
                'existing output not resumable with these exact inputs; select a new --output-id')
        return data
    c5_control,c5_normal,parts=normal_propagator_c5(kernel,uv,eft,basis)
    if control_hash is not None:
        control=neutrino.load(existing_control)
        require(control.get('model')==model and
                control.get('status')=='RestrictedKernelControlSelfBilinearNotMatched' and
                control.get('source_continuous_sha256')==upstream_hash and
                control.get('matching_source_sha256')==kernel_hash,
                'existing no-self-bilinear C5 control has incompatible provenance')
        saved=neutrino.complex_matrix(control['c5_at_1tev_gev_inverse'],'saved control C5')
        require(np.allclose(saved,c5_control,rtol=4e-9,atol=1e-19),
                'recomputed no-mixing C5 fails saved numerical control')
    rotated_vertices=native_vertex_transport(model,basis)
    vertex_file=output_root/(tag+'_rotated_extra_vertices.json')
    neutrino.atomic_json(vertex_file,rotated_vertices)
    low=neutrino.low_energy(model,c5_normal,eft)
    control_low=neutrino.low_energy(model,c5_control,eft)
    mlow=neutrino.complex_matrix(low['mnu_at_100gev_ev'],'low mass')
    m0=neutrino.complex_matrix(control_low['mnu_at_100gev_ev'],'control low mass')
    payload={
        'format':FORMAT,'status':STATUS,'model':model,
        'quadratic_kind':source['quadratic_kind'],
        'source_continuous_report':record.relative_to(ROOT).as_posix(),
        'source_continuous_sha256':upstream_hash,
        'source_matching_kernel':kernel_path.relative_to(ROOT).as_posix(),
        'source_matching_kernel_sha256':kernel_hash,
        'source_runner_sha256':script_hash,
        'source_native_tensor_sha256':native_tensor_hash,
        'rotated_extra_vertex_report':vertex_file.relative_to(ROOT).as_posix(),
        'rotated_extra_vertices':rotated_vertices['extra_quartic_vertices_rotated'],
        'source_control_sha256':control_hash,
        'matching_kernel_kind':kind,
        'mass_basis_species':basis.species,
        'su2_component_dimension':basis.dim,
        'self_mass_coefficient_gev2':{
            'real':float(basis.mu_squared.real),'imag':float(basis.mu_squared.imag)},
        'normal_scalar_mass_squared_eigenvalues_gev2':list(basis.eigen_masses_squared),
        'normal_scalar_mass_eigenvalues_gev':np.sqrt(basis.eigen_masses_squared).tolist(),
        'normal_propagator_eigenstate_weights':[0.5,0.5],
        'original_topology_c5_at_1tev_without_self_mass':neutrino.matrix_payload(c5_control),
        'normal_propagator_c5_by_scalar_mass_eigenstate_at_1tev':[
             neutrino.matrix_payload(z) for z in parts],
        'normal_propagator_resummed_c5_at_1tev':neutrino.matrix_payload(c5_normal),
        'relative_normal_propagator_c5_correction':_rel(c5_normal,c5_control),
        'relative_normal_propagator_low_energy_mnu_correction':_rel(mlow,m0),
        'baseline_takagi_masses_ev_descending':control_low['takagi_masses_ev_descending'],
        'anomalous_propagator_nonzero':bool(basis.mu_squared!=0),
        'anomalous_prop_matching_included':False,
        'additional_independent_quartic_direct_matching_included':False,
        'normal_channel_scalar_bilinear_included':True,
        'physical_matchete_rgbeta_map_certified':False,
        'sequential_wilson_matching_performed':False,
        'full_general_scalar_c5_certified':False,
        'scope':('Exact quadratic mass resummation of the NORMAL scalar propagator '
                 'inside the existing topology one-loop matching expression. '
                 'ANOMALOUS propagator + independent HH scalar vertices and '
                 'their Matchete/RGBeta CG normalization remain unmatched.'),
        **low,
    }
    neutrino.atomic_json(target,payload)
    return payload



def combined_16_summary(*, source_root: Path, output_root: Path,
                        control_output_id: str) -> dict:
    """Merge six normal-resummed cases with ten saved candidate calculations.

    Never overwrites neutrino_v1; each row retains its scope/limitations.
    """
    require(set(SELF).issubset(set(continuous.MODELS)) and len(continuous.MODELS)==16,
            'unexpected global ordinary T3 model inventory')
    rows=[]
    for model in continuous.MODELS:
        tag=continuous.model_tag(model)
        continuous_path=source_root/tag/'continuous_trajectory.json'
        if model in SELF:
            file=output_root/(tag+'_self_bilinear_c5.json')
            allowed={STATUS}
            category='self-bilinear-normal-propagator-only'
        else:
            file=source_root/control_output_id/(tag+'_neutrino.json')
            allowed={'CandidateMassEigenBCEvolvedYukawasSMWeinberg',
                     'CandidateNoAdditionalQuadraticMatchedC5'}
            category='existing-candidate-matching'
        raw=neutrino.load(file)
        # The two upstream report families deliberately use different names
        # for their (False) physical-certification flags.  Do not accept an
        # absent flag: explicitly check the schema for the selected family.
        if model in SELF:
            certification_guard = (
                raw.get('format') == FORMAT and
                raw.get('full_general_scalar_c5_certified') is False and
                raw.get('normal_channel_scalar_bilinear_included') is True and
                raw.get('anomalous_prop_matching_included') is False and
                raw.get('additional_independent_quartic_direct_matching_included') is False
            )
        else:
            certification_guard = (
                raw.get('format') == neutrino.FORMAT and
                raw.get('full_general_scalar_physical_c5_certified') is False and
                raw.get('full_sequential_threshold_matching_performed') is False and
                raw.get('additional_quadratic_c5_omitted') is False
            )
        require(raw.get('model')==model and raw.get('status') in allowed and
                raw.get('source_continuous_sha256')==_hash(continuous_path) and
                certification_guard,
                model+': incompatible source result in 16-model table')
        masses=raw.get('takagi_masses_ev_descending')
        require(isinstance(masses,list) and len(masses)==3 and
                all(np.isfinite(x) and float(x)>=0 for x in masses),
                model+': missing finite low-energy Takagi masses')
        rows.append({'model':model,'status':raw['status'],
                     'matching_category':category,
                     'source_report':file.relative_to(ROOT).as_posix(),
                     'source_report_sha256':_hash(file),
                     'largest_neutrino_mass_ev':masses[0],
                     'middle_neutrino_mass_ev':masses[1],
                     'smallest_neutrino_mass_ev':masses[2],
                     **raw['oscillation_observables_normal_ordering']})
    summary={'format':'t3-continuous-16-model-partially-resummed-neutrino-v1',
             'models_processed':len(rows),
             'models_self_bilinear_normal_channel_resummed':len(SELF),
             'models_existing_candidate_c5':len(rows)-len(SELF),
             'physical_full_general_scalar_matching_complete':False,
             'anomalous_self_scalar_and_extra_quartic_matching_complete':False,
             'sequential_wilson_matching_complete':False,
             'rows':rows}
    neutrino.atomic_json(output_root/'combined_16_summary.json',summary)
    columns=['model','status','matching_category','largest_neutrino_mass_ev',
             'middle_neutrino_mass_ev','smallest_neutrino_mass_ev',
             'delta_m21_sq_ev2','delta_m3l_sq_ev2',
             'sin2_theta12','sin2_theta13','sin2_theta23']
    with (output_root/'combined_16_summary.csv.tmp').open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)
    (output_root/'combined_16_summary.csv.tmp').replace(output_root/'combined_16_summary.csv')
    return summary

def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--all',action='store_true')
    group.add_argument('--model',choices=SELF)
    parser.add_argument('--run-id',required=True,type=continuous.run_id)
    parser.add_argument('--source-root',type=Path,default=continuous.OUTPUT)
    parser.add_argument('--control-output-id',default='neutrino_v1')
    parser.add_argument('--output-id',default='self_bilinear_v1')
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--combine-only',action='store_true',
        help='Rebuild the 16-model summary from already completed model files; do not rerun matching')
    args=parser.parse_args()
    continuous.run_id(args.output_id)
    folder=args.source_root/args.run_id
    require(folder.is_dir(),'missing continuous trajectory directory '+str(folder))
    out=folder/args.output_id
    require(out!=folder/args.control_output_id,'output-id must differ from existing reduced control')
    if args.combine_only:
        require(args.all,'--combine-only requires --all (all six self-bilinear outputs)')
        combined=combined_16_summary(source_root=folder,output_root=out,
                                     control_output_id=args.control_output_id)
        print('Combined candidate results:',combined['models_processed'],
              'models; full general scalar matching NOT certified',flush=True)
        print('Combined:',out/'combined_16_summary.csv')
        return 0
    models=SELF if args.all else (args.model,)
    rows=[]
    for model in models:
        try:
            result=calculate(model,source_root=folder,output_root=out,
                             control_output_id=args.control_output_id,resume=args.resume)
            masses=result['takagi_masses_ev_descending']
            row={'model':model,'status':result['status'],
                 'largest_neutrino_mass_ev':masses[0],
                 'middle_neutrino_mass_ev':masses[1],
                 'smallest_neutrino_mass_ev':masses[2],
                 'relative_normal_c5_delta':result['relative_normal_propagator_c5_correction'],
                 'anomalous_matching_included':False,
                 **result['oscillation_observables_normal_ordering']}
            print(f'{model}: self-bilinear NORMAL channel included; masses={masses}; '
                  f'anomalous+extra quartic matching NOT included',flush=True)
        except (ValueError,KeyError,TypeError,RuntimeError,OSError,OverflowError,
                FloatingPointError,np.linalg.LinAlgError) as exc:
            row={'model':model,'status':'Blocked','error':str(exc)[:600]}
            print(f'{model}: BLOCKED {exc}',flush=True)
        rows.append(row)
        summary={'format':FORMAT,'models_requested':list(models),
                 'models_processed':len(rows),
                 'models_successful':sum(z['status']==STATUS for z in rows),
                 'models_blocked':sum(z['status']=='Blocked' for z in rows),
                 'normal_propagator_self_bilinear_resummation_performed':True,
                 'all_independent_general_scalar_matching_performed':False,
                 'full_sequential_wilson_matching_performed':False,
                 'rows':rows}
        neutrino.atomic_json(out/'summary.json',summary)
        columns=['model','status','largest_neutrino_mass_ev','middle_neutrino_mass_ev',
                 'smallest_neutrino_mass_ev','relative_normal_c5_delta',
                 'delta_m21_sq_ev2','delta_m3l_sq_ev2',
                 'sin2_theta12','sin2_theta13','sin2_theta23','error']
        with (out/'summary.csv.tmp').open('w',newline='',encoding='utf-8') as handle:
            writer=csv.DictWriter(handle,fieldnames=columns,extrasaction='ignore')
            writer.writeheader();writer.writerows(rows)
        (out/'summary.csv.tmp').replace(out/'summary.csv')
    if args.all and not any(z['status']=='Blocked' for z in rows):
        combined=combined_16_summary(source_root=folder,output_root=out,
                                     control_output_id=args.control_output_id)
        print('Combined candidate results:',combined['models_processed'],
              'models; full general scalar matching NOT certified',flush=True)
        print('Combined:',out/'combined_16_summary.csv')
    print('Summary:',out/'summary.json')
    return 2 if any(z['status']=='Blocked' for z in rows) else 0

if __name__=='__main__':
    raise SystemExit(main())
