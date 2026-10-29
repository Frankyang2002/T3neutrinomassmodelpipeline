"""Selected additional-quartic LLSS -> Weinberg source on full-scalar B/C(-1) EFT1 trajectories.

This is a *single source insertion* beyond the existing direct-only sequential
candidate. It reads the previously verified Stage17/19 Matchete/Weinberg
certificates and runs the already-implemented Stage18 source integrator on the
current saved *full-general-scalar* EFT1 background. The heavy-fermion PL LLSS
seed stays frozen at mu_F; the additional quartics run between mu_F and mu_S.
The independently saved direct hierarchical C5 is retained without re-matching.

This is NOT full EFT1 Wilson mixing, complete threshold matching, nor a
certified complete general-scalar neutrino mass calculation.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import sympy as sp

from scripts import RunT3DirectWilsonSequential as direct
from scripts import RunT3FullScalarContinuous as continuous
from scripts import RunT3FullScalarNeutrino as neutrino

ROOT = continuous.ROOT
MODELS = ('B_-1', 'C_-1')
STAGE_NAMES = {'B_-1': 'B_m1', 'C_-1': 'C_m1'}
FORMAT = 't3-bc-extra-quartic-selected-wilson-sequential-candidate-v1'
STATUS = 'SelectedBCExtraQuarticSourceTransportedOnFullScalarEFT1'
DEFAULT_OUTPUT_ID = 'bc_extra_wilson_sequential_v1'
DIRECT_OUTPUT_ID = direct.DEFAULT_OUTPUT_ID
LOOP_FACTOR = 16.0 * math.pi**2
ACTIVE_QUARTICS = ('lambdaHHdagS2S2', 'lambdaHHdagS1barS1bar')


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError('B/C selected extra Wilson: '+message)


def digest(path: Path) -> str:
    path = Path(path)
    require(path.is_file(), 'missing source '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def sources(model: str) -> dict[str, Path]:
    require(model in MODELS, 'only B_-1 and C_-1 supported')
    stage = STAGE_NAMES[model]
    match, rge = direct.model_sources(model)
    tag = continuous.model_tag(model)
    return {
        'stage17': ROOT/'output/general_t3_stage17'/f'{stage}_extra_eft1_scalar_pair_beta.json',
        'stage19': ROOT/'output/general_t3_stage19'/f'{stage}_selected_extra_weinberg_transport.json',
        'seed': ROOT/f'output/full/hypercharge/T3_{model[0]}_alpha_m1/data/eft1_after_F_wilson_seed.json',
        'matching': match,
        'wilson_rge': rge,
        'eft1_beta': continuous.BETA_DIR/f'{tag}_eft1_quadratic_live_rgbeta.json',
    }


def _sym(expression: str) -> sp.Expr:
    return sp.sympify(expression, locals={
        'conjugate': sp.conjugate, 'sqrt': sp.sqrt,
        'I': sp.I, 'MF': sp.Symbol('MF'),
    })


def certify_saved_selected_source(model: str, stage17: Mapping[str, Any],
                                  stage19: Mapping[str, Any], seed: Mapping[str, Any],
                                  wilson: Mapping[str, Any],
                                  matched: Mapping[str, Any]) -> tuple[dict, dict]:
    """Reject a source whose four components or saved full-tensor certificate changed.

    Stage19 previously checked all 30 independent SU(2) components. We verify
    the current four Stage17 expressions *exactly*, not only a status field.
    This is a reuse of the saved certificate, NOT a new full tensor audit.
    """
    from Numerical.experimental.GeneralT3WeinbergProjection import (
        _analytic_mixed_template, verify_production_c5_normalization,
    )
    stage_name = STAGE_NAMES[model]
    require(stage17.get('model')==stage_name and
            stage17.get('status')=='ExtraEFT1ScalarPairMixingDerived' and
            stage17.get('source_convention')==
            'lambda_Matchete=-lambda_RGBeta, same named scalar field basis' and
            stage17.get('full_eft1_wilson_running_implemented') is False,
            model+': invalid saved Stage17 scalar-pair source')
    require(seed.get('Level')==1 and seed.get('IntegratedFields')==['F'] and
            seed.get('ActiveHeavyFields')==['S1','S2'] and
            wilson.get('status')=='Success' and matched.get('status')=='Success',
            model+': fermion-threshold matching artifacts changed')
    require(stage19.get('status')=='SelectedExtraWeinbergSourceTransported' and
            stage19.get('model')==stage_name and
            stage19.get('full_coupled_eft1_wilson_evolution_implemented') is False and
            stage19.get('finite_scalar_matching_implemented') is False and
            stage19.get('complete_two_loop_result') is False,
            model+': unexpected Stage19 source status or completeness')
    certificate = stage19.get('projection_certificate',{})
    require(certificate.get('status')=='FullExtraWeinbergTensorProjectionVerified' and
            certificate.get('model')==stage_name and
            certificate.get('independent_components_verified')==30 and
            certificate.get('mixed_components_matched_to_stage17')==4 and
            certificate.get('all_component_residuals_zero') is True,
            model+': missing full SU(2) projection certificate')
    expected = _sym(certificate['beta_kappa_16pi2_expression'])
    template = _analytic_mixed_template()
    rows = stage17.get('components',[])
    require(isinstance(rows,list) and len(rows)==4 and
            stage17.get('nonzero_weinberg_real_components')==4,
            model+': Stage17 component inventory changed')
    observed = {}
    for row in rows:
        pair = tuple(row['higgs_real_indices'])
        require(pair not in observed and isinstance(row['beta_16pi2_expression'],str),
                model+': repeated or invalid Stage17 component')
        observed[pair] = _sym(row['beta_16pi2_expression'])
    require(set(observed)==set(template), model+': Stage17 Higgs-real tensor support changed')
    for pair, value in template.items():
        require(sp.simplify(observed[pair]-sp.sympify(value)*expected)==0,
                model+': Stage17 source differs from certified full tensor at '+str(pair))
    # Check against the current authoritative hierarchical direct C5 normalization,
    # not merely the earlier Stage19 report text.
    now = verify_production_c5_normalization(wilson,matched)
    historic = stage19.get('normalization_certificate',{})
    require(historic.get('status')=='ProductionDirectC5NormalizationVerified' and
            historic.get('direct_running_log_coefficient_equals_beta_kappa') is True and
            sp.simplify(_sym(historic['beta_kappa_16pi2'])-
                        _sym(now['beta_kappa_16pi2']))==0,
            model+': saved and current direct normalization do not agree')
    return certificate, now


@dataclass(frozen=True)
class EFT1TrajectoryView:
    """Read-only Stage18-compatible view of the native general-scalar trajectory."""
    mu_gev: np.ndarray
    states: tuple[Any,...]
    def state_at_index(self, index: int) -> Any:
        return self.states[index]


def running_scalar_background(model: str, record: Mapping[str,Any],
                              sidecar: Mapping[str,Any], *, samples: int) -> tuple[EFT1TrajectoryView,dict]:
    """Reintegrate *only EFT1* and verify parity with its saved endpoint."""
    from Numerical.experimental.GeneralT3ScalarRGState import (
        GeneralT3ScalarRGState,pack_general_scalar_state,
        unpack_general_scalar_state,
    )
    from scripts import IntermediateEFTRunning as s34
    require(model in MODELS and samples>=17 and samples<=199 and samples%2==1,
            'odd samples in [17,199] required')
    require(sidecar.get('status')=='CandidateQuadraticAndFullQuarticEFT1Integrated' and
            sidecar.get('format')==continuous.qeft.FORMAT and
            sidecar.get('model')==model and record['quadratic_kind']=='mixed',
            model+': saved EFT1 sidecar/mixed mass classification changed')
    registry=continuous.load_registry(continuous.REGISTRY)
    base=s34.decode_intermediate_state(sidecar['initial_baseline'])
    extra={name:neutrino.complex_pair(raw,name)
           for name,raw in sidecar['initial_new_quartics'].items()}
    initial=GeneralT3ScalarRGState(base,extra).validated(registry)
    y0,layout=pack_general_scalar_state(initial,registry)
    initial_mu2=neutrino.complex_pair(sidecar['mu2_initial_gev2'],'initial mixed mass^2')
    y0=np.concatenate((y0,[initial_mu2.real,initial_mu2.imag]))
    beta_path=sources(model)['eft1_beta']
    betas=continuous.beta_inventory(continuous.read(beta_path), model=model,
               stage='EFT1', representation=base.representation, layout=layout,
               registry=registry, kind='mixed')
    result=continuous.integrate_segment(y0,mu_start=s34.START,mu_end=s34.STOP,
               layout=layout,registry=registry,stage='EFT1',kind='mixed',
               betas=betas,points=samples)
    states=tuple(unpack_general_scalar_state(result.y[:layout.size,i], layout,
               registry,uv=False,mu_gev=float(math.exp(result.t[i])))
               for i in range(samples))
    require(all(np.isfinite([complex(getattr(st.baseline,key)).real,
                             complex(getattr(st.baseline,key)).imag]).all()
                for st in states for key in ACTIVE_QUARTICS),
            model+': missing/nonfinite B/C additional portal couplings')
    view=EFT1TrajectoryView(np.exp(np.asarray(result.t,dtype=float)),
                            tuple(st.baseline for st in states))
    final=states[-1]
    previous=neutrino.decode_eft(record['eft1_scalar_endpoint_state'])
    errors={}
    for key in ('mS1Sq','mS2Sq','lambdaT3',*ACTIVE_QUARTICS):
        now=complex(getattr(final.baseline,key))
        prior=complex(getattr(previous,key))
        error=abs(now-prior)/max(abs(prior),1e-12)
        require(error<2e-6,model+': endpoint parity failed for '+key)
        errors[key]=float(error)
    require(set(final.extra)==set(sidecar['final_new_quartics']),
            model+': new-quartic inventory changed')
    for key, now in final.extra.items():
        expected=neutrino.complex_pair(sidecar['final_new_quartics'][key],key)
        require(abs(now-expected)<=2e-6*max(abs(expected),1e-10),
                model+': endpoint new-quartic mismatch '+key)
    final_mu2=complex(result.y[layout.size,-1],result.y[layout.size+1,-1])
    expected_mu2=neutrino.complex_pair(sidecar['mu2_final_gev2'],'final mixed mass^2')
    require(abs(final_mu2-expected_mu2)<=2e-6*max(abs(expected_mu2),1e-10),
            model+': quadratic endpoint parity failed')
    errors['mixed_mu2']=abs(final_mu2-expected_mu2)/max(abs(expected_mu2),1e-12)
    lams=np.asarray([st.baseline.lambdaT3 for st in states],dtype=complex)
    from scipy.integrate import simpson
    topology_integral=complex(simpson(lams,x=result.t))
    return view,{
        'eft1_nfev':int(result.nfev),'samples':samples,
        'source_portals_at_F':{key:continuous.cjson(getattr(view.states[0],key))
                               for key in ACTIVE_QUARTICS},
        'source_portals_at_S':{key:continuous.cjson(getattr(view.states[-1],key))
                               for key in ACTIVE_QUARTICS},
        'lambdaT3_integral':topology_integral,
        'endpoint_parity':errors,
        'beta_path':beta_path,
    }


def selected_source_integral(model: str, stage17: Mapping[str,Any],
                             certificate: Mapping[str,Any], normalization: Mapping[str,Any],
                             view: EFT1TrajectoryView,
                             masses: np.ndarray, y1:np.ndarray,y2:np.ndarray) -> tuple[np.ndarray,dict,dict]:
    """Reuse Stage18 component integration and Stage19 SU2 projection unchanged."""
    from Numerical.experimental.GeneralT3ExtraWilsonTransport import (
        compile_stage17_sources,integrate_source_over_eft1,
    )
    from Numerical.experimental.GeneralT3WeinbergProjection import (
        project_stage18_integral,
    )
    components=compile_stage17_sources(STAGE_NAMES[model],stage17)
    stage18=integrate_source_over_eft1(components,view,masses,y1,y2)
    stage18['model']=STAGE_NAMES[model]
    correction, projection=project_stage18_integral(stage18,STAGE_NAMES[model],
                                                     certificate,normalization)
    return neutrino.validated_c5(correction),stage18,projection


def frozen_selected_source(model: str, stage17: Mapping[str,Any],
                           certificate: Mapping[str,Any], normalization: Mapping[str,Any],
                           view: EFT1TrajectoryView,
                           masses:np.ndarray,y1:np.ndarray,y2:np.ndarray) -> np.ndarray:
    """Independent constant-quartic control using the same verified source tensor."""
    from Numerical.experimental.GeneralT3ExtraWilsonTransport import (
        compile_stage17_sources,frozen_seed_flavor_contractions,
    )
    from Numerical.experimental.GeneralT3WeinbergProjection import (
        _analytic_mixed_template,
    )
    components=compile_stage17_sources(STAGE_NAMES[model],stage17)
    f1,f2=frozen_seed_flavor_contractions(masses,y1,y2)
    first=view.states[0]
    matrices={c.higgs_pair:c.evaluate(first.lambdaHHdagS2S2,
              first.lambdaHHdagS1barS1bar,f1,f2) for c in components}
    template=_analytic_mixed_template()
    denom=sum(abs(v)**2 for v in template.values())
    beta=sum(np.conj(v)*matrices[k] for k,v in template.items())/denom
    return neutrino.validated_c5(beta*math.log(view.mu_gev[-1]/view.mu_gev[0])/LOOP_FACTOR)


def compute(model: str, *,folder:Path,out:Path, direct_output_id: str=DIRECT_OUTPUT_ID,
            samples:int=65,resume:bool=False,reverify_tensors:bool=False) -> dict:
    require(model in MODELS,'unsupported model')
    tag=continuous.model_tag(model)
    continuous_file=folder/tag/'continuous_trajectory.json'
    record=neutrino.verify_continuous(model,continuous_file)
    sidecar_file=folder/tag/record['saved_candidate_eft1_endpoint']
    sidecar=neutrino.load(sidecar_file)
    direct_file=folder/direct_output_id/f'{tag}_direct_wilson.json'
    direct_result=neutrino.load(direct_file)
    require(direct_result.get('format')==direct.FORMAT and
            direct_result.get('status')==direct.STATUS and
            direct_result.get('model')==model and
            direct_result.get('full_sequential_wilson_matching_performed') is False and
            direct_result.get('full_general_scalar_matching_certified') is False,
            model+': previous direct-only matching result incompatible')
    require(direct_result.get('saved_points')==samples,
            model+': previous direct source used different --samples; use matching samples')
    paths=sources(model)
    paths.update({'continuous':continuous_file,'sidecar':sidecar_file,
                  'direct':direct_file,
                  'runner':Path(__file__).resolve(),
                  'continuous_runner':ROOT/'scripts/RunT3FullScalarContinuous.py',
                  'neutrino_runner':ROOT/'scripts/RunT3FullScalarNeutrino.py',
                  'direct_runner':ROOT/'scripts/RunT3DirectWilsonSequential.py'})
    hashes={name:digest(path) for name,path in paths.items()}
    target=out/f'{tag}_extra_wilson.json'
    if target.is_file():
        old=neutrino.load(target)
        require(resume and old.get('format')==FORMAT and old.get('status')==STATUS and
                old.get('input_sha256')==hashes and old.get('saved_points')==samples and
                old.get('reverified_tensor_projection')==reverify_tensors,
                model+': output exists or inputs changed; use --resume or new --output-id')
        return old
    stage17=neutrino.load(paths['stage17'])
    stage19=neutrino.load(paths['stage19'])
    seed=neutrino.load(paths['seed'])
    wilson=neutrino.load(paths['wilson_rge'])
    matched=neutrino.load(paths['matching'])
    certificate,normalization=certify_saved_selected_source(model,stage17,stage19,
                                                             seed,wilson,matched)
    if reverify_tensors:
        from Numerical.experimental.GeneralT3WeinbergProjection import (
            verify_full_seed_weinberg_projection,
        )
        check=verify_full_seed_weinberg_projection(STAGE_NAMES[model],seed,stage17)
        require(sp.simplify(_sym(check['beta_kappa_16pi2_expression'])-
                 _sym(certificate['beta_kappa_16pi2_expression']))==0 and
                 check.get('independent_components_verified')==30,
                model+': live tensor recertification differs from saved certificate')
    uv=neutrino.decode_uv(record['uv_fermion_threshold_state'])
    eft=neutrino.decode_eft(record['eft1_scalar_endpoint_state'])
    mf,y1,y2,heavy_residual=direct.heavy_inputs(uv)
    view,background=running_scalar_background(model,record,sidecar,samples=samples)
    saved_lambda=neutrino.complex_pair(direct_result['lambdaT3_log_integral'],
                                      'saved direct lambdaT3 integral')
    require(abs(saved_lambda-background['lambdaT3_integral'])<=
            5e-7*max(abs(saved_lambda),1e-11),
            model+': EFT1 source integration differs from saved direct source')
    correction,stage18,projection=selected_source_integral(model,stage17,
        certificate,normalization,view,mf,y1,y2)
    control=frozen_selected_source(model,stage17,certificate,normalization,
                                  view,mf,y1,y2)
    baseline=neutrino.complex_matrix(direct_result['c5_at_scalar_threshold_gev_inverse'],
                                    'saved direct C5 at scalar threshold')
    baseline=neutrino.validated_c5(baseline)
    updated=neutrino.validated_c5(baseline+correction)
    low=neutrino.low_energy(model,updated,eft)
    # Run the baseline through the *same* code and same SM endpoint for a clean
    # difference; independent earlier direct masses also provide a parity check.
    baseline_low=neutrino.low_energy(model,baseline,eft)
    earlier=np.asarray(direct_result['takagi_masses_ev_descending'],dtype=float)
    new_base=np.asarray(baseline_low['takagi_masses_ev_descending'],dtype=float)
    require(np.allclose(earlier,new_base,rtol=4e-8,atol=1e-10),
            model+': saved direct-only low-energy mass parity failed')
    low_delta=neutrino.validated_c5(
        neutrino.complex_matrix(low['c5_at_100gev_gev_inverse'],'updated C5 at low')-
        neutrino.complex_matrix(baseline_low['c5_at_100gev_gev_inverse'],'direct C5 at low'))
    rel=float(np.linalg.norm(correction)/max(np.linalg.norm(baseline),1e-300))
    output={
        'format':FORMAT,'status':STATUS,'model':model,'saved_points':samples,
        'input_sha256':hashes,'source_projection_certificate':certificate,
        'normalization_certificate':normalization,
        'source_model':STAGE_NAMES[model],
        'direct_source_result_path':str(direct_file.relative_to(folder)),
        'source_tensor_integration_status':stage18['status'],
        'source_integral_projection':projection,
        'source_convention':stage17['source_convention'],
        'source_quartics_at_F':background['source_portals_at_F'],
        'source_quartics_at_S':background['source_portals_at_S'],
        'direct_C5_at_scalar_threshold_gev_inverse':neutrino.matrix_payload(baseline),
        'selected_extra_delta_C5_at_scalar_threshold_gev_inverse':neutrino.matrix_payload(correction),
        'selected_extra_frozen_quartic_delta_C5_gev_inverse':neutrino.matrix_payload(control),
        'direct_plus_selected_extra_C5_at_scalar_threshold_gev_inverse':neutrino.matrix_payload(updated),
        'selected_extra_delta_C5_at_100gev_gev_inverse':neutrino.matrix_payload(low_delta),
        'relative_extra_C5_vs_direct_at_scalar_threshold':rel,
        'relative_extra_C5_vs_frozen_extra':float(np.linalg.norm(correction-control)/
                                                   max(np.linalg.norm(control),1e-300)),
        'direct_only_masses_ev_descending':earlier.tolist(),
        'updated_masses_ev_descending':low['takagi_masses_ev_descending'],
        'eft1_source_integration_nfev':background['eft1_nfev'],
        'eft1_source_samples':samples,
        'eft1_endpoint_parity':background['endpoint_parity'],
        'heavy_mass_basis_residual':heavy_residual,
        'reverified_tensor_projection':bool(reverify_tensors),
        'stage18_component_integrals':stage18['higgs_real_components'],
        'frozen_fermion_threshold_LLSS_seed':True,
        'additional_quartic_scalar_pair_source_included':True,
        'homogeneous_LLSS_wilson_running_included':False,
        'other_general_quartic_mixing_included':False,
        'finite_one_loop_scalar_threshold_matching_included':False,
        'full_coupled_eft1_wilson_evolution_implemented':False,
        'complete_general_scalar_sequential_matching':False,
        'physical_general_scalar_matching_certified':False,
        'uv_running_repeated':False,
        'eft1_scalar_running_repeated_for_source':True,
        'notes':('Selected Stage17/19 additional quartic scalar-pair source integrated '
                 'on continuous full-scalar EFT1 trajectory, then added to the '
                 'saved direct-only hierarchical C5. This is a single source '
                 'insertion with the heavy LLSS seed frozen at the F threshold; '
                 'NOT a complete EFT1 Wilson beta, finite scalar matching, '
                 'nor certified general-scalar C5 prediction.'),
        **low,
    }
    neutrino.atomic_json(target,output)
    return output


def summary(out:Path, selected:tuple[str,...]) -> dict:
    rows=[]
    for model in selected:
        result=neutrino.load(out/f'{continuous.model_tag(model)}_extra_wilson.json')
        require(result.get('format')==FORMAT and result.get('status')==STATUS and
                result.get('model')==model and
                result.get('complete_general_scalar_sequential_matching') is False,
                model+': summary incompatible result')
        row={
            'model':model,'status':result['status'],
            'relative_extra_C5_vs_direct_at_scalar_threshold':result['relative_extra_C5_vs_direct_at_scalar_threshold'],
            'relative_extra_C5_vs_frozen_extra':result['relative_extra_C5_vs_frozen_extra'],
            'baseline_largest_mnu_ev':result['direct_only_masses_ev_descending'][0],
            'updated_largest_mnu_ev':result['updated_masses_ev_descending'][0],
            'updated_middle_mnu_ev':result['updated_masses_ev_descending'][1],
            'updated_smallest_mnu_ev':result['updated_masses_ev_descending'][2],
            'projection_relative_residual':result['source_integral_projection']['projection_relative_residual'],
        }
        rows.append(row)
    import csv
    payload={'format':FORMAT,'status':'SelectedBCExtraSourceSummary',
             'model_count':len(rows),'complete_general_scalar_sequential_matching':False,
             'selection_scope':'B/C alpha=-1 selected extra scalar-pair source only',
             'rows':rows}
    neutrino.atomic_json(out/'summary.json',payload)
    columns=list(rows[0]) if rows else []
    temp=out/'summary.csv.tmp'
    with temp.open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns)
        writer.writeheader();writer.writerows(rows)
    temp.replace(out/'summary.csv')
    return payload


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    sel=parser.add_mutually_exclusive_group(required=True)
    sel.add_argument('--all',action='store_true')
    sel.add_argument('--model',choices=MODELS)
    parser.add_argument('--run-id',required=True,type=continuous.run_id)
    parser.add_argument('--source-root',type=Path,default=continuous.OUTPUT)
    parser.add_argument('--direct-output-id',default=DIRECT_OUTPUT_ID)
    parser.add_argument('--output-id',default=DEFAULT_OUTPUT_ID)
    parser.add_argument('--samples',type=int,default=65)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--combine-only',action='store_true')
    parser.add_argument('--reverify-tensors',action='store_true',
                        help='repeat expensive full 30-component Stage19 tensor check')
    args=parser.parse_args()
    continuous.run_id(args.output_id)
    continuous.run_id(args.direct_output_id)
    require(args.output_id!=args.direct_output_id and
            args.output_id not in ('neutrino_v1','self_anomalous_v1',
             'self_anomalous_massdep_v1','self_bilinear_v1'),
            'invalid output directory: would overwrite an earlier study')
    require(args.samples>=17 and args.samples<=199 and args.samples%2==1,
            '--samples must be odd integer 17..199')
    require(not args.combine_only or args.all, '--combine-only requires --all')
    folder=args.source_root/args.run_id
    require(folder.is_dir(),'saved continuous run not found: '+str(folder))
    out=folder/args.output_id
    chosen=MODELS if args.all else (args.model,)
    failed=[]
    if not args.combine_only:
        for model in chosen:
            try:
                data=compute(model,folder=folder,out=out,
                        direct_output_id=args.direct_output_id,
                        samples=args.samples,resume=args.resume,
                        reverify_tensors=args.reverify_tensors)
                print(f'{model}: selected extra EFT1 Wilson source OK; '
                      f'|delta C5|/|direct C5|={data["relative_extra_C5_vs_direct_at_scalar_threshold"]:.9g}; '
                      f'vs frozen extra={data["relative_extra_C5_vs_frozen_extra"]:.9g}; '
                      f'masses={data["updated_masses_ev_descending"]}',flush=True)
            except (OSError,ValueError,KeyError,TypeError,RuntimeError,
                    ArithmeticError,np.linalg.LinAlgError) as error:
                failed.append(model)
                print(f'{model}: BLOCKED: {type(error).__name__}: {error}',flush=True)
    if failed:
        print('Blocked models:',', '.join(failed),flush=True)
        return 2
    if args.all:
        completed=summary(out,chosen)
        print('Selected B/C table:',str(out/'summary.csv'),
              'models=',completed['model_count'],flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
