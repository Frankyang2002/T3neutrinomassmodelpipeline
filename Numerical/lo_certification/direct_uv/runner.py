"""All-16 source-gated exact direct UV benchmark; write only owned output tree."""
from __future__ import annotations
import csv
import json
from pathlib import Path
from typing import Any

from Numerical.direct_one_loop.core import evaluate, _complex
from Numerical.direct_one_loop.comparison import compare_same_point
from .certification import (all_models, CLASSES, validate_components,
                            source_evidence, load, digest)
from .channels import partition_exact


def _save(path:Path,value:Any):
    path.parent.mkdir(parents=True,exist_ok=True)
    candidate=path.with_suffix(path.suffix+'.tmp')
    candidate.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    candidate.replace(path)


def _perhaps(path:Path):
    return load(path) if path.is_file() else None


def verify_source_charge_support(structural:dict[str,Any], data:dict[str,Any],
                                 kind:str, alpha:int, charge_entry:dict|None,
                                 cg_entry:dict|None) -> dict[str,Any]:
    """Check source T3 *mixed* charge support, not every Yukawa-only charge.

    Source `charge_resolved_loop_weights` records charges in the simultaneous
    Y1 x Y2 x T3MixCG contraction. Individual Yukawa vertices can have additional
    SU(2)-allowed components that do not participate in this topology. Those
    components remain present in the exact full-UV calculation.
    """
    vertex_rows={v['name']:v for v in structural['vertex_dictionary']}
    y1={float(q) for q in vertex_rows['y1_HC']['charge_channels']}
    y2={float(q) for q in vertex_rows['y2_HC']['charge_channels']}
    observed={float(q) for q in structural['observed_scalar_charge_channels']}
    common=y1 & y2

    expected=(None if charge_entry is None else
              {float(q) for q in charge_entry['internal_charge_channels']})
    cg_charges=(None if cg_entry is None else
                {float(ch['charge']) for ch in cg_entry['channels'] if ch.get('nonzero')})
    # Independent source audits must agree on *topology* channel support.
    if expected is not None and cg_charges is not None and cg_charges != expected:
        raise ValueError('Native CG and source T3 topology charge audits disagree: '
                         f'allowed={sorted(expected)}, weights={sorted(cg_charges)}')
    target=expected if expected is not None else cg_charges
    if target is not None and not target.issubset(common):
        raise ValueError('Missing source T3 topology charge in one or both Yukawas: '
                         f'expected={sorted(target)}, Y1={sorted(y1)}, Y2={sorted(y2)}')
    if not (y1|y2).issubset(observed):
        raise ValueError('Observed vertex charge union is internally inconsistent')

    # Charge-preserving S1/S2 quadratic mixing can be read separately without
    # falsely calling a Yukawa-only charge an additional topology contribution.
    from .certification import charge_bases
    import numpy as np
    d1=CLASSES[kind][0]
    scalar_charges=charge_bases(kind,alpha)['scalar']
    A=_complex(data['scalar_A_gev2'],'scalar A')
    cross=np.abs(A[:d1,d1:])
    cutoff=1e-12*max(float(np.max(cross,initial=0.)),1e-100)
    mixing={float(scalar_charges[i]) for i,j in np.argwhere(cross>cutoff)}
    # Never infer the source lambdaT3 coupling dictionary from the Hessian.
    info={
        'source_yukawa1_charge_support':sorted(y1),
        'source_yukawa2_charge_support':sorted(y2),
        'source_yukawa_union_charge_support':sorted(y1|y2),
        'all_component_vertex_charge_support':sorted(observed),
        'both_yukawa_charge_support':sorted(common),
        'source_expected_internal_charge_channels':None if expected is None else sorted(expected),
        'native_CG_candidate_charge_support':None if cg_charges is None else sorted(cg_charges),
        'source_expected_channel_coverage_checked':expected is not None and expected.issubset(common),
        'native_CG_candidate_channel_coverage_checked':cg_charges is not None and cg_charges.issubset(common),
        'yukawa_only_charge_channels_not_in_T3_contraction':(
            None if target is None else sorted((y1|y2)-target)),
        'candidate_scalar_A_cross_mixing_charge_support':sorted(mixing),
        'source_charge_contract':'T3 weights refer to Y1 x Y2 x scalar-mixing intersection, not Y1/Y2 union',
        'all_original_vertices_retained_in_exact_loop':True,
    }
    if cg_entry is not None:
        info.update({
            'native_source_CG_candidate_weights':{
                str(ch['charge']):str(ch['complex_weight']) for ch in cg_entry['channels']},
            'native_source_CG_candidate_total_coefficient':cg_entry.get('total_group_coefficient'),
            'native_CG_algebraic_evidence_only_not_physical_normalization':True,
            'verified_physical_loop_multiplicity':cg_entry.get('verified_physical_loop_multiplicity'),
        })
    return info


def run(root:Path, output:Path, input_dir:Path|None=None,
        eft_dir:Path|None=None, weak_fraction:float=.01) -> dict:
    root=root.resolve(); output=output.resolve()
    if not output.is_relative_to((root/'output/general_observable/lo_certification_direct_uv').resolve()):
        raise ValueError('Output must remain inside owned output/general_observable/lo_certification_direct_uv')
    if (output/'all_16_certification_status.json').exists() or any((output/'models').glob('*.json')):
        raise FileExistsError('Existing certification output found; use a NEW --output-dir within the owned output tree')
    if input_dir is None:
        input_dir=root/'output/general_observable/lo_direct_one_loop/provisional_inputs'
    source=root/'output/general_observable/lo_direct_one_loop'
    audits={
        'prior':_perhaps(source/'convention_batch1/integrated_convention_audit.json'),
        'portal':_perhaps(source/'source_audits/native_portal_identity_decomposition.json'),
        'charges':_perhaps(source/'source_audits/charge_channels.json'),
        'cg_weights':_perhaps(source/'source_audits/charge_resolved_loop_weights.json'),
    }
    rows=[]
    for short,key in all_models():
        result={'model':short,'status':'not_evaluated','certified':False,
                'benchmark':f'Numerical/benchmarks/fitted_reference/{key}.json',
                'input':str(input_dir/(short+'.json')),
                'blocking':[]}
        try:
            bench_path=root/result['benchmark']
            inp_path=input_dir/(short+'.json')
            if not bench_path.is_file() or not inp_path.is_file():
                raise FileNotFoundError('Missing frozen source benchmark or pre-existing component input')
            bench=load(bench_path); inp=load(inp_path)
            if inp.get('model')!=short: raise ValueError('Input model id disagrees with catalogue')
            kind,alpha=short[0],int(short[2:])
            structural=validate_components(inp,kind,alpha)
            evidence=source_evidence(root,short,key,inp,bench,audits['prior'],audits['portal'],
                                   {'benchmark':bench_path})
            if not all(evidence['checks'][x] for x in
                       ('source_benchmark_matches','source_representation_matches','input_at_source_base_scale')):
                raise ValueError('Source benchmark provenance or UV input scale mismatch')
            # Source charge audit: expected charges describe the nonzero
            # Y1 x Y2 x T3MixCG topology, NOT the union of all vertices.
            rep=bench['representation']
            def matching_entry(audit):
                return next((entry for entry in (audit or {}).get('models',[])
                             if all(entry['model'][k]==rep[r] for k,r in
                                    [('dS1','d_s1'),('dS2','d_s2'),('dF','d_f'),
                                     ('alpha','alpha')])),None)
            structural.update(verify_source_charge_support(
                structural,inp,kind,alpha,
                matching_entry(audits['charges']),
                matching_entry(audits['cg_weights'])))
            direct=evaluate(inp,retain_individual_terms=True)
            if not direct['diagnostics']['uv_finite']:
                raise ValueError('The exact loop has an uncancelled UV pole')
            if not direct['diagnostics']['neutrino_mass_matrix_symmetric']:
                raise ValueError('The full UV neutrino mass is not symmetric')
            channels=partition_exact(inp,kind,alpha,direct,weak_fraction=weak_fraction)
            physically=evidence['physically_certified']
            # Even if an external input sets a certification boolean, no source
            # stage is promoted without independent cross-package witnesses.
            assert not physically
            comparison={'status':'not_requested'}
            if eft_dir is not None:
                p=eft_dir/(short+'.json')
                if p.is_file():
                    comparison=compare_same_point(direct,load(p))
                else:comparison={'status':'same_input_eft_reference_not_supplied'}
            result.update(status='exact_uv_provisional',certified=False,
                          masses_ev=direct['total']['neutrino_masses_ev'],
                          cancellation_ratio=direct['diagnostics']['cancellation_ratio_total_over_sum'],
                          observed_charge_channels=structural['observed_scalar_charge_channels'],
                          blocking=evidence['blocking_physical_conventions'])
            _save(output/'models'/f'{short}.json',{
                'format':'t3-direct-uv-lo-certification-v1',
                'status':'provisional_not_physical_certification',
                'model':short,'representation':bench['representation'],
                'source_files':{
                    'benchmark':str(bench_path.relative_to(root)),
                    'benchmark_sha256':digest(bench_path),
                    'component_input':str(inp_path),
                    'component_input_sha256':digest(inp_path)},
                'EFT_two_thresholds':{
                    'UV_to_EFT1':'integrate ALL heavy BSM fermions together at mu_F',
                    'EFT1_to_SMEFT':'integrate BOTH heavy BSM scalars together at mu_S',
                    'mu_F_gev':bench['scales']['mu_fermion_threshold_gev'],
                    'mu_S_gev':bench['scales']['mu_scalar_threshold_gev'],
                    'UV_reference_mu_gev':bench['base_state']['mu_gev'],
                    'same_scale_as_EFT_threshold':False,
                    'SM_fermions_retained':True},
                'structural':structural,'source_certification':evidence,
                'exact_direct_one_loop':direct,'charge_and_suppression':channels,
                'eft_comparison':comparison})
        except (ValueError,FileNotFoundError,KeyError,IndexError,TypeError,
                RuntimeError,ArithmeticError) as exc:
            result.update(status='blocked',blocking=[type(exc).__name__+': '+str(exc)])
        rows.append(result)
        print(f'{short:5s} {result["status"]} '
              + ('m3={:.6g} eV'.format(result['masses_ev'][-1]) if 'masses_ev' in result else
                 '; '.join(result['blocking'])[:100]),flush=True)
    summary={
        'format':'t3-direct-uv-lo-all16-v1',
        'counts':{
            'catalogue':len(rows),
            'exact_uv_provisional':sum(r['status']=='exact_uv_provisional' for r in rows),
            'physical_source_certified':sum(r['certified'] for r in rows),
            'blocked':sum(r['status']=='blocked' for r in rows)},
        'benchmark_provenance':'frozen fitted_reference at base_state.mu_gev; no invented additional couplings',
        'baseline_unmodified':True,
        'models':rows,
    }
    _save(output/'all_16_certification_status.json',summary)
    output.mkdir(parents=True,exist_ok=True)
    with (output/'all_16_summary.csv').open('w',encoding='utf-8',newline='') as fp:
        wr=csv.writer(fp)
        wr.writerow(['model','status','source_certified','m1_ev','m2_ev','m3_ev','full_over_sum','blocking'])
        for row in rows:
            wr.writerow([row['model'],row['status'],row['certified'],
                         *row.get('masses_ev',[None]*3),row.get('cancellation_ratio'),
                         '; '.join(row['blocking'])])
    return summary
