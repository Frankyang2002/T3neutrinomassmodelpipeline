"""Read-only A(-4) paired *conditional* EFT1-HHSS sensitivity study.

Project path: scripts/RunT3Worker1A4ConditionalIVP.py

Uses the production-selected experimental simultaneous A(-4) runner without editing it.
The Batch3B alternative is NOT a physically certified scalar tensor. Only the
EFT1 Wilson RHS HHSS is substituted; scalar finite matching remains unchanged.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE_OUTPUT = ROOT / 'output/studies/general_scalar_extended/experimental_simultaneous_eft1_seven_term'
DEFAULT_RESPONSE = 'output/audits/batch3b_quartic_fixed_20261111_234059'
SCHEMA = 't3-worker1-A4-conditional-IVP-sensitivity-v1'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as fp:
        for b in iter(lambda: fp.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def complex_matrix(x):
    a = np.asarray(x['real'], dtype=float) + 1j*np.asarray(x['imag'], dtype=float)
    if a.shape != (3,3) or not np.isfinite(a).all():
        raise ValueError('Invalid 3x3 Weinberg coefficient')
    return a


def compare_matrices(a, b):
    """a=alternative, b=baseline. Use complex overlap to remove global phase/scale."""
    a,b = np.asarray(a,complex),np.asarray(b,complex)
    denom = float(np.vdot(b,b).real)
    if denom == 0 or np.linalg.norm(a) == 0:
        return {'status':'undefined_zero_matrix'}
    coefficient = np.vdot(b,a) / denom
    real_coefficient = float(np.vdot(b,a).real / denom)
    nb = float(np.linalg.norm(b))
    na = float(np.linalg.norm(a))
    return {'status':'compared',
            'alternative_over_baseline_norm':na/nb,
            'relative_difference_to_baseline':float(np.linalg.norm(a-b)/nb),
            'best_complex_coefficient':{'real':float(coefficient.real),'imag':float(coefficient.imag)},
            'residual_after_complex_rescaling_relative_to_alternative':float(np.linalg.norm(a-coefficient*b)/na),
            'best_real_coefficient':real_coefficient,
            'residual_after_real_rescaling_relative_to_alternative':float(np.linalg.norm(a-real_coefficient*b)/na),
            'symmetry_baseline':float(np.linalg.norm(b-b.T)/nb),
            'symmetry_alternative':float(np.linalg.norm(a-a.T)/na)}


def mass_angles(case):
    masses = np.array(case.get('masses_ev_ascending', []),dtype=float)
    p = case.get('prediction',{})
    out = {'masses_ev_ascending':masses.tolist(),
           'mass_ratios_m2_over_m3_m1_over_m3':(
               [float(masses[1]/masses[2]),float(masses[0]/masses[2])]
               if masses.shape==(3,) and masses[2]>0 else None),
           'sin2_theta12':p.get('sin2_theta12'),
           'sin2_theta13':p.get('sin2_theta13'),
           'sin2_theta23':p.get('sin2_theta23'),
           'delta_m21_sq_ev2':p.get('delta_m21_sq_ev2'),
           'delta_m3l_sq_ev2':p.get('delta_m3l_sq_ev2')}
    return out


def matrix_stage_fractions(cases):
    """Difference-norm ablation sensitivities, not an additive decomposition."""
    reference = complex_matrix(cases['1111']['c5_at_low_scale'])
    norm = float(np.linalg.norm(reference))
    if norm == 0: return {'status':'undefined_zero_reference'}
    return {stage:float(np.linalg.norm(reference-complex_matrix(cases[key]['c5_at_low_scale']))/norm)
            for stage,key in [('UV_off','0111'),('EFT1_renormalisable_off','1011'),
                              ('EFT1_Wilson_off','1101'),('SMEFT_off','1110'),
                              ('both_EFT1_components_off','1001')]}


def report_pair(baseline, alternative):
    if baseline.get('model')!='A_-4' or alternative.get('model')!='A_-4':
        raise ValueError('Only A_-4 allowed')
    if baseline.get('status')!='Computed' or alternative.get('status')!='Computed':
        raise ValueError('Both runs must be Computed; inspect error.json')
    b,a = baseline['cases'], alternative['cases']
    mandatory={'1111','0111','1001','1110','1011','1101'}
    if not mandatory <= set(a) or not mandatory <= set(b):
        raise ValueError('Missing stage-ablation cases')
    # W=0 must be a strict negative control: the changed HHSS enters only
    # the EFT1 Wilson RHS, not finite scalar matching or renormalisable RGEs.
    negative_controls={}
    for key in sorted(set(a)&set(b)):
        if key[2]=='0':
            source=complex_matrix(b[key]['c5_at_low_scale'])
            conditional=complex_matrix(a[key]['c5_at_low_scale'])
            deviation=float(np.linalg.norm(source-conditional)/max(np.linalg.norm(source),1e-100))
            negative_controls[key]=deviation
            if deviation>1e-9:
                raise ValueError(f'W=0 negative control {key} unexpectedly changed: {deviation}')
    comparison={}
    for key in sorted(set(b)&set(a)):
        x,y=b[key],a[key]
        comparison[key]={
            'low_energy_c5':compare_matrices(complex_matrix(y['c5_at_low_scale']),complex_matrix(x['c5_at_low_scale'])),
            'scalar_threshold_c5':compare_matrices(complex_matrix(y['c5_at_scalar_threshold']),complex_matrix(x['c5_at_scalar_threshold'])),
            'baseline':mass_angles(x),'alternative':mass_angles(y),
            'observable_differences':{k:float(mass_angles(y)[k]-mass_angles(x)[k])
                for k in ('sin2_theta12','sin2_theta13','sin2_theta23','delta_m21_sq_ev2','delta_m3l_sq_ev2')
                if mass_angles(y)[k] is not None and mass_angles(x)[k] is not None}
        }
    return {'schema':SCHEMA,'status':'Computed_Conditional_NotPhysicalCertification',
            'comparison':comparison,
            'W_off_negative_control_relative':negative_controls,
            'stage_fraction_proxy_baseline':matrix_stage_fractions(b),
            'stage_fraction_proxy_alternative':matrix_stage_fractions(a),
            'matching_scope':'Only EFT1 Wilson RHS HHSS changes; scalar finite bubble and renormalisable beta are baseline',
            'physical_alternative_certified':False,
            'absolute_mass_and_flavour_physical_robustness_certified':False}


def execute_variant(args):
    if args.variant not in ('source','batch3b'):
        raise ValueError('Unknown variant')
    from scripts import RunT3SimultaneousEFT1A4 as a4
    if args.variant=='batch3b':
        from Numerical.running.T3CoupledQuarticContext import CandidateQuarticContext
        from Numerical.running.T3ScalarQuarticLinearResponse import symmetry_residual
        response_root=Path(args.response_root).resolve()
        ctx=CandidateQuarticContext.load(ROOT,response_root,'A_-4')
        original_integrator=a4.integrate_simultaneously
        def conditional_integrator(eft_start,beta_eft,source_c,q_basis,mu_s,**kwargs):
            original_tensor=q_basis.tensor
            checks={'points':0,'max_source_candidate_relative':0.,'min_candidate_norm':float('inf')}
            def conditional_tensor(state):
                source=original_tensor(state)
                full=ctx.full_reference(state)
                if full.shape!=(12,12,12,12) or symmetry_residual(full)>1e-9:
                    raise ValueError('Cached conditional quartic is malformed')
                candidate=np.asarray(full[:4,:4,4:,4:],float)
                if candidate.shape!=source.shape or not np.isfinite(candidate).all():
                    raise ValueError('Invalid candidate HHSS')
                checks['points']+=1
                checks['min_candidate_norm']=min(checks['min_candidate_norm'],float(np.linalg.norm(candidate)))
                checks['max_source_candidate_relative']=max(checks['max_source_candidate_relative'],
                    float(np.linalg.norm(candidate-source)/max(np.linalg.norm(source),1e-100)))
                return candidate
            q_basis.tensor=conditional_tensor
            try:
                result=original_integrator(eft_start,beta_eft,source_c,q_basis,mu_s,**kwargs)
            finally:
                q_basis.tensor=original_tensor
            result['conditional_candidate_hhss_evaluations']=checks
            return result
        a4.integrate_simultaneously=conditional_integrator
    argv=['--source-run',args.source_run,'--run-id',args.run_id,'--models','A_-4']
    if args.check_only: argv.append('--check-only')
    status=a4.main(argv)
    if status != 0: return status
    if args.check_only: return 0
    destination=BASE_OUTPUT/args.run_id/'A_m4/cases.json'
    if not destination.is_file():raise FileNotFoundError(destination)
    saved=json.loads(destination.read_text(encoding='utf-8-sig'))
    if saved.get('status')!='Computed':raise RuntimeError('Blocked physics run: inspect error.json')
    # No alteration of saved solver outputs. Write mode provenance separately.
    meta={'schema':SCHEMA,'variant':args.variant,'source_run':args.source_run,
          'run_id':args.run_id,'source_path':str(destination),
          'source_sha256':digest(destination),
          'conditional_tensor_only':args.variant=='batch3b',
          'production_implementation_modified':False,
          'physical_tensor_certified':False}
    if args.variant=='batch3b':
        meta['batch3b_A4_npz_sha256']=digest(Path(args.response_root)/'A_-4/quartic_response.npz')
    path=BASE_OUTPUT/args.run_id/'WORKER1_MODE.json'
    if path.exists():raise FileExistsError(path)
    path.write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    return 0


def _resolve_case(run_id):
    return BASE_OUTPUT/run_id/'A_m4/cases.json'


def report_runs(source_id, candidate_id, out_dir):
    if source_id==candidate_id:
        raise ValueError('Baseline and conditional run IDs must differ')
    bs,cs = _resolve_case(source_id),_resolve_case(candidate_id)
    for p in (bs,cs):
        if not p.is_file():raise FileNotFoundError(p)
    def check_mode(run_id,mode,p):
        meta=BASE_OUTPUT/run_id/'WORKER1_MODE.json'
        if not meta.is_file():raise FileNotFoundError(f'Missing Worker1 provenance: {meta}')
        m=json.loads(meta.read_text(encoding='utf-8'))
        if m['variant']!=mode or digest(p)!=m['source_sha256']:
            raise ValueError('Run mode/hash mismatch')
        return m
    source_meta=check_mode(source_id,'source',bs)
    candidate_meta=check_mode(candidate_id,'batch3b',cs)
    if source_meta['source_run']!=candidate_meta['source_run']:
        raise ValueError('Not identical fitted UV source runs')
    report=report_pair(json.loads(bs.read_text(encoding='utf-8-sig')),
                       json.loads(cs.read_text(encoding='utf-8-sig')))
    report['source_provenance']=source_meta
    report['candidate_provenance']=candidate_meta
    report['source_case_sha256']=digest(bs)
    report['candidate_case_sha256']=digest(cs)
    out=Path(out_dir)
    if out.exists():raise FileExistsError(f'Refusing overwrite: {out}')
    out.mkdir(parents=True)
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    rows=['# A(-4) conditional scalar-quartic IVP sensitivity','','**Not a certified physical tensor alternative.**',
          'Only the EFT1 Wilson HHSS RHS was varied. All UV inputs, renormalisable running, scalar finite matching, and SMEFT conventions remain the same.','',
          '| Quantity | Baseline | Conditional Batch3B |','|---|---:|---:|']
    r=report['comparison']['1111']; b=r['baseline']; a=r['alternative'];c=r['low_energy_c5']
    for key in ('sin2_theta12','sin2_theta13','sin2_theta23','delta_m21_sq_ev2','delta_m3l_sq_ev2'):
        rows.append(f'| {key} | {b.get(key)} | {a.get(key)} |')
    rows += ['',f"Low-energy C5 normalized difference: {c.get('relative_difference_to_baseline')}",
             f"C5 shape residual after optimal complex rescaling: {c.get('residual_after_complex_rescaling_relative_to_alternative')}",
             '', '## Stage freeze proxies (difference norms, not additive fractions)',
             '| Stage disabled | Baseline | Conditional Batch3B |','|---|---:|---:|']
    for key in report['stage_fraction_proxy_baseline']:
        rows.append(f"| {key} | {report['stage_fraction_proxy_baseline'][key]} | {report['stage_fraction_proxy_alternative'][key]} |")
    rows += ['','## Interpretation','This measures sensitivity to an uncertified tensor substitution in the retained simultaneous EFT1 Wilson equations. It does not fix the native/RGBeta coordinate convention or certify the scalar matching.']
    (out/'REPORT.md').write_text('\n'.join(rows)+'\n',encoding='utf-8')
    print('Wrote',out/'REPORT.md',out/'report.json')
    return 0


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--variant',choices=['source','batch3b'],help='Run one A4 physics variant')
    p.add_argument('--source-run',default='thesis16_merged_v1')
    p.add_argument('--run-id',help='Fresh run name; must not exist')
    p.add_argument('--response-root',default=DEFAULT_RESPONSE)
    p.add_argument('--check-only',action='store_true')
    p.add_argument('--compare',action='store_true',help='Compare completed Worker1 source/batch3b runs')
    p.add_argument('--source-id',default='worker1_A4_source_v1')
    p.add_argument('--candidate-id',default='worker1_A4_Batch3B_conditional_v1')
    p.add_argument('--out-dir',default='output/audits/worker1_A4_conditional_IVP_v1')
    p.add_argument('--both',action='store_true',help='Run both variants in fresh separate Python processes, then compare')
    a=p.parse_args(argv)
    if a.compare: return report_runs(a.source_id,a.candidate_id,a.out_dir)
    if a.both:
        for variant,run_id in (('source',a.source_id),('batch3b',a.candidate_id)):
            cmd=[sys.executable,'-X','utf8','-m','scripts.RunT3Worker1A4ConditionalIVP',
                 '--variant',variant,'--run-id',run_id,'--source-run',a.source_run,
                 '--response-root',a.response_root]
            if a.check_only:cmd.append('--check-only')
            subprocess.run(cmd,cwd=ROOT,check=True)
        if not a.check_only:return report_runs(a.source_id,a.candidate_id,a.out_dir)
        return 0
    if not a.variant or not a.run_id:p.error('--variant and --run-id required unless --both or --compare')
    return execute_variant(a)

if __name__=='__main__':
    raise SystemExit(main())
