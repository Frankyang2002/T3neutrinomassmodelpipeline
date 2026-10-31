#!/usr/bin/env python3
"""v23 consolidated B/C(-1) leading one-loop benchmark and certification gate.

Consumes locally passed v18/v19/v21/v22 audits and the archived source files.
Independently evaluates the Euclidean three-mass radial loop integral, verifies
exact direct-vs-sequential C5 matrices, and optionally evolves a frozen-input
REFERENCE using the existing SM+Weinberg RGE solver from MX to 100 GeV.

CRITICAL: The SM-only running background from 1e7 GeV through heavy thresholds
is a deliberately nonphysical proxy; the source UV-EFT1-EFT2 running/parameter
mapping and a separate Feynman-rule derivation of the absolute UV sign are not
certified. Never promote this output to a physical T3 prediction.
"""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
AUD=Path('output/audits')
BASE=Path('output/general_t3_mass_eigen_matchete')
MODELS={'B_m1':('B_-1',1,'T3_dS1_2_dS2_2_dF_1_alpha_m1.json'),
        'C_m1':('C_-1',3,'T3_dS1_2_dS2_2_dF_3_alpha_m1.json')}
V18=AUD/'t3_bc_three_fermion_triangle_v18/summary.json'
V19=AUD/'t3_bc_hard_soft_subtraction_v19/summary.json'
V21=AUD/'t3_bc_topology_closure_v21/summary.json'
V22=AUD/'t3_bc_direct_vs_sequential_full_mass_v22/summary.json'

class CompletionError(ValueError):pass

def need(ok,msg):
    if not ok:raise CompletionError(msg)

def read_json(path):
    need(path.is_file(),'Missing local evidence: '+str(path))
    data=json.loads(path.read_text(encoding='utf-8-sig'))
    need(isinstance(data,dict),'Expected JSON object: '+str(path))
    return data

def digest(path):
    need(path.is_file(),'Required source missing: '+str(path))
    return sha256(path.read_bytes()).hexdigest()

def pair_matrix(values,name):
    import numpy as np
    need(isinstance(values,list) and len(values)==3,name+': expected 3 rows')
    need(all(isinstance(row,list) and len(row)==3 for row in values),name+': expected 3 columns')
    try:
        need(all(isinstance(z,list) and len(z)==2 for row in values for z in row),
             name+': every complex cell requires [real,imag]')
        a=np.asarray([[complex(*z) for z in row] for row in values],dtype=complex)
    except (TypeError,ValueError) as ex:
        raise CompletionError(name+': invalid complex pair') from ex
    need(np.isfinite(a).all(),name+': nonfinite matrix')
    return a

def complex_json(a):
    return [[[float(z.real),float(z.imag)] for z in row] for row in a]

def norm_relative(a,b):
    import numpy as np
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-250))

def masses_from_c5(C5,vev):
    import numpy as np
    mass=-.5*vev**2*C5*1e9
    singular=np.sort(np.linalg.svd(mass,compute_uv=False))
    return {'C5_GeV_inverse':complex_json(C5),
            'mnu_eV':complex_json(mass),
            'singular_masses_eV_ascending':[float(x) for x in singular],
            'delta_m21_squared_eV2':float(singular[1]**2-singular[0]**2),
            'delta_m31_squared_eV2':float(singular[2]**2-singular[0]**2)}

def independent_triangle(A,B,C):
    """Standalone absolutely convergent Euclidean radial integral.

    Phi=-A * int_0^infty dt*t/((t+A)(t+B)(t+C)).
    Minus sign convention inherited from v17, not independently fixed from
    Minkowski Majorana Feynman rules; this test fixes loop-mass dependence.
    """
    import mpmath as mp
    a,b,c=(mp.mpf(str(x)) for x in (A,B,C))
    need(all(x>0 for x in (a,b,c)),'Nonpositive mass squared')
    with mp.workdps(43):
        val=-a*mp.quad(lambda t:t/((t+a)*(t+b)*(t+c)),[0,mp.inf])
        return val

def master_closed(A,B,C):
    """Independent closed form in terms of stable logarithmic divided differences."""
    import mpmath as mp
    a,b,c=(mp.mpf(str(x)) for x in (A,B,C))
    with mp.workdps(43):
        if abs(b-c)<mp.mpf('1e-24')*max(b,c):
            if abs(a-b)<mp.mpf('1e-24')*a:return -mp.mpf('0.5')
            return -a*(a*mp.log(a/b)-(a-b))/(a-b)**2
        def f(x):
            if abs(a-x)<mp.mpf('1e-24')*a:return mp.mpf(1)
            return x*mp.log(a/x)/(a-x)
        return -a*(f(b)-f(c))/(b-c)

def verify_integrals(mf,mx,my):
    import mpmath as mp
    rows=[];max_error=0.
    for M in mf:
        A=float(M)**2;X=mx*mx;Y=my*my
        for label,B,C in (('XX',X,X),('XY',X,Y),('YY',Y,Y)):
            a=independent_triangle(A,B,C)
            b=master_closed(A,B,C)
            rel=float(abs(a-b)/max(mp.mpf('1'),abs(b)))
            need(rel<1e-25,f'Independent triangle master mismatch at M={M}: {label}: {rel}')
            max_error=max(max_error,rel)
            rows.append({'MF_GeV':float(M),'scalar_pair':label,
                         'signed_phi_from_direct_integral':float(a),
                         'signed_phi_from_closed_form':float(b),
                         'relative_residual':rel})
    return rows,max_error

def decode_flavour(obj,label):
    import numpy as np
    need(isinstance(obj,dict) and 'real' in obj and 'imag' in obj,label+': missing matrix components')
    x=np.asarray(obj['real'],dtype=float)+1j*np.asarray(obj['imag'],dtype=float)
    need(x.shape==(3,3) and np.isfinite(x).all(),label+': invalid matrix')
    return x

def reference_sm_at_mx(ref,mx):
    """Proxy SM-only running from the 1e7 GeV reference to MX, K=0.

    This is NOT the actual matched B/C SM background. UV BSM couplings and
    threshold corrections are missing; do not label physical running.
    """
    import numpy as np
    from Numerical.running.SMWeinbergEvolution import SMWeinbergInitialConditions,evolve_sm_weinberg
    sm=ref['base_state']['sm']
    initial=SMWeinbergInitialConditions(
       gY=float(sm['gY']),g2=float(sm['g2']),g3=float(sm['g3']),
       lambdaH=float(sm['lambdaH']),
       ye=decode_flavour(sm['ye'],'ye'),
       yu=decode_flavour(sm['yu'],'yu'),
       yd=decode_flavour(sm['yd'],'yd'),
       K=np.zeros((3,3),complex)).validated()
    source_scale=float(ref['base_state']['mu_gev'])
    need(source_scale==1e7 and 0<mx<source_scale,'Invalid frozen SM reference scale')
    state=evolve_sm_weinberg(initial,source_scale,mx,rtol=1e-9,atol=1e-12)
    need(state.solver_success and np.linalg.norm(state.K)<1e-20,
         'SM-only reference running failed or generated C5 from zero')
    return state

def low_scale_proxy(state,mx,C5,low):
    import numpy as np
    from Numerical.running.SMWeinbergEvolution import SMWeinbergInitialConditions,evolve_sm_weinberg
    initial=SMWeinbergInitialConditions(
       gY=state.gY,g2=state.g2,g3=state.g3,lambdaH=state.lambdaH,
       ye=state.ye,yu=state.yu,yd=state.yd,K=C5).validated()
    need(0<low<mx,'Require 0 < low scale < MX')
    evolved=evolve_sm_weinberg(initial,mx,low,rtol=1e-9,atol=1e-12)
    need(evolved.solver_success and np.isfinite(evolved.K).all(),
         'Low-scale SM-only reference running failed')
    need(norm_relative(evolved.K,evolved.K.T)<1e-12,
         'Running broke symmetric Majorana C5')
    return evolved

def validate_audit_chain(root,v18,v19,v21,v22):
    need(v18.get('format')=='t3-bc-frozen-reference-three-fermion-triangle-v18' and
         v18.get('running_performed') is False and v18.get('physical_C5') is None,
         'v18 reference scope changed')
    need(v19.get('format')=='t3-bc-one-loop-hard-soft-msbar-restricted-v19' and
         v19.get('status')=='RestrictedTriangleHardSoftMSbarSubtractionEstablished' and
         v19.get('full_T3_renormalized_matching_certified') is False,
         'v19 restricted subtraction changed')
    need(v21.get('format')=='t3-bc-one-loop-renormalizable-topology-closure-v21' and
         v21.get('one_loop_diagram_topology_complete_for_archived_UV_under_scope') is True and
         v21.get('full_sequential_physical_matching_certified') is False,
         'v21 archived topology closure changed')
    need(v22.get('format')=='t3-bc-independent-all-heavy-vs-sequential-finite-amplitude-v22' and
         v22.get('all_direct_vs_sequential_full_mass_agree') is True and
         v22.get('full_running_performed') is False and
         v22.get('independent_non_Matchete_UV_Feynman_rule_calculation_performed') is False and
         v22.get('physical_C5') is None,
         'v22 finite matching status inconsistent')
    # A report from a prior local commit must not certify modified matching input.
    source_contract=v21.get('original_UV_source',{}).get('source_files_sha256',{})
    need(isinstance(source_contract,dict) and len(source_contract)==5,
         'Missing archived UV source-hash contract')
    for source_path,expected in source_contract.items():
        rel=Path(source_path)
        need(not rel.is_absolute() and '..' not in rel.parts and
             isinstance(expected,str) and digest(root/rel)==expected,
             'Archived UV source changed since v21: '+source_path)
    need(digest(root/Path('output/scalar_relevance/operator_classification.json'))
         ==v21.get('UV_operator_manifest_sha256'),
         'UV operator manifest changed since v21')
    for stem,(model,df,file) in MODELS.items():
        a18=v18['models'][stem];a21=v21['models'][stem];a22=v22['models'][stem]
        need(a18.get('model')==model and a21.get('model')==model and a22.get('model')==model,
             model+': mismatched audit labels')
        need(a21.get('one_loop_two_trivalent_graphs_excluded_for_archived_UV') is True and
             a21.get('one_loop_one_quartic_channels_exhaustive_for_archived_UV') is True,
             model+': topology branch not complete')
        need(a22.get('direct_vs_sequential_full_mass_agreement') is True and
             a22.get('total_independent_checks')==12 and
             a22.get('direct_vs_sequential_max_relative_residual',1.)<2e-8 and
             a22.get('direct_full_3x3_C5_vs_v18_relative_residual',1.)<2e-8 and
             a22.get('full_physical_C5_certified') is False,
             model+': direct vs sequential finite identity incomplete')
        filenames={
          'UV_model':(BASE/f'{stem}_final_matchete.json',a21.get('UV_archive_sha256')),
          'direct_Matchete':(BASE/'one_loop_probe'/f'{stem}_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt',a22.get('one_loop_direct_coefficient_sha256')),
          'direct_report':(BASE/'one_loop_probe'/f'{stem}_final_one_loop_probe.json',a22.get('all_heavy_direct_json_sha256')),
          'mass_eigen':(BASE/f'{stem}_final_parameters.json',a18.get('mass_basis_sha256')),
          'reference':(Path('Numerical/benchmarks/fitted_reference')/file,a18.get('fitted_reference_sha256')),
        }
        for label,(path,expected) in filenames.items():
            need(isinstance(expected,str) and len(expected)==64 and digest(root/path)==expected,
                 model+': source changed since audits: '+label)
        require_src=v21['source_matched_F_Y_X_sha256_from_v20_verified_unchanged'][stem]
        require_pathmap={
          'F_raw':BASE/'eft1_fermion_first_oneloop'/f'{stem}_final_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt',
          'F_full':BASE/'eft1_fermion_first_oneloop'/f'{stem}_final_eft1_fermion_first_oneloop_eft1_one_loop_correction.inputform.txt',
          'Y_raw':BASE/'eft1_scalar_Y_only_canonical_probe'/f'{stem}_final_y_only_canonical_match_raw_scalar_llhh.inputform.txt',
          'Y_full':BASE/'eft1_scalar_Y_only_canonical_probe'/f'{stem}_final_y_only_canonical_match_scalar_one_loop_correction.inputform.txt',
          'X_raw':BASE/'eft1_scalar_xy_fresh_kernel_v3'/f'{stem}_final_x_from_saved_y_match_x_followup_x_raw_llhh.inputform.txt',
          'X_full':BASE/'eft1_scalar_xy_fresh_kernel_v3'/f'{stem}_final_x_from_saved_y_match_x_followup_x_loop_difference.inputform.txt',
        }
        need(set(require_src)==set(require_pathmap),model+': changed F/Y/X source contract')
        for name,path in require_pathmap.items():
            need(digest(root/path)==require_src[name],model+': F/Y/X matched source changed: '+name)
    return True

def run(root,out,low,skip_running=False):
    import numpy as np
    v18=read_json(root/V18);v19=read_json(root/V19)
    v21=read_json(root/V21);v22=read_json(root/V22)
    validate_audit_chain(root,v18,v19,v21,v22)
    report={
      'format':'t3-bc-consolidated-leading-one-loop-reference-v23',
      'status':'OneLoopTopologyAndFiniteMatchingReproduced_LowScaleReferenceOnly',
      'source_audits_sha256':{str(p):digest(root/p) for p in (V18,V19,V21,V22)},
      'source_checks_passed':True,
      'archived_UV_one_loop_topology_complete_within_scope':True,
      'hard_soft_MSbar_triangle_subtraction_verified':True,
      'direct_Matchete_all_heavy_equals_sequential_full_mass':True,
      'loop_mass_integral_independently_integrated':True,
      'independent_non_Matchete_UV_Feynman_rule_absolute_sign_certified':False,
      'full_UV_EFT1_EFT2_coupling_running_and_finite_scheme_certified':False,
      'SM_only_reference_running_performed':not skip_running,
      'SM_only_background_from_1e7_GeV_is_not_T3_running':True,
      'y1_y2_MF_and_quartics_frozen_from_1e7_GeV_reference':True,
      'physical_C5_certified':False,
      'physical_C5':None,'physical_neutrino_masses':None,
      'models':{},
      'remaining_physics_blockers':[
        'Independent UV Minkowski Feynman-rule derivation of vertex/spinor/global sign and normalization (both Matchete routes use same source conventions)',
        'Full UV-to-EFT1-to-EFT2 renormalizable coupling evolution and threshold map; current SM-only background is a proxy',
        'Complete finite field, parameter, and EOM matching scheme for physical sequential C5 outside the restricted triangle sector',
        'Nonzero independently allowed mixed scalar quadratic if extending beyond archived zero-extra-quadratic source'
      ]}
    for stem,(model,df,filename) in MODELS.items():
        old=v18['models'][stem]; direct=v22['models'][stem]
        ref=read_json(root/Path('Numerical/benchmarks/fitted_reference')/filename)
        need(tuple(ref['representation'][k] for k in ('d_s1','d_s2','d_f','alpha'))==(2,2,df,-1),
             model+': reference representation mismatch')
        need(float(ref['base_state']['mu_gev'])==1e7,model+': reference scale changed')
        masses=[float(x) for x in old['MF_GeV']]
        mx=float(old['MX_GeV']);my=float(old['MY_GeV']);vev=float(old['reference_vev_GeV'])
        need(len(masses)==3 and all(m>my>mx>low>0 for m in masses),model+': invalid threshold hierarchy')
        near=pair_matrix(direct['direct_full_3x3_C5_GeV_inverse_NOT_PHYSICAL'],model+' direct C5')
        previous=pair_matrix(old['strict_and_full_triangle']['full_mass']['reference_C5_3x3_GeV_inverse_NOT_PHYSICAL'],
                             model+' v18 reference C5')
        need(norm_relative(near,near.T)<1e-13,model+': C5 symmetric projection changed')
        c5_diff=norm_relative(near,previous)
        need(c5_diff<2e-8,model+': v22/v18 C5 discrepancy')
        before=masses_from_c5(near,vev)
        oldvals=old['strict_and_full_triangle']['full_mass']['reference_Takagi_singular_masses_eV_ascending_NOT_PHYSICAL']
        need(max(abs(a-b) for a,b in zip(before['singular_masses_eV_ascending'],oldvals))<1e-10,
             model+': old reference masses do not match C5')
        integrals,error=verify_integrals(masses,mx,my)
        row={'model':model,'MF_GeV':masses,'MX_GeV':mx,'MY_GeV':my,
             'independent_loop_integral_relative_residual_max':error,
             'independent_triangle_checks':integrals,
             'C5_direct_vs_v18_relative_residual':c5_diff,
             'direct_vs_sequential_max_relative_residual':direct['direct_vs_sequential_max_relative_residual'],
             'full_mass_matching_scale_reference_GeV':mx,
             'full_mass_fixed_order_at_MX_NOT_PHYSICAL':before,
             'UV_fitted_reference_source_sha256':old['fitted_reference_sha256'],
             'reference_input_Yukawas_at_1e7_GeV_NOT_RUN':True,
             'complete_matched_scalar_EFT_RGE':False,
             'physical_neutrino_masses':None}
        if not skip_running:
            sm_at_mx=reference_sm_at_mx(ref,mx)
            outsm=low_scale_proxy(sm_at_mx,mx,near,low)
            last=masses_from_c5(outsm.K,vev)
            fixed=np.asarray(before['singular_masses_eV_ascending'],float)
            flow=np.asarray(last['singular_masses_eV_ascending'],float)
            row['SM_only_reference_running']={
                 'status':'SM-only illustrative running, NOT T3 matched EFT',
                 'SM_initial_scale_GeV':float(ref['base_state']['mu_gev']),
                 'reference_C5_insertion_scale_GeV':mx,
                 'low_scale_GeV':float(low),
                 'SM_background_UV_to_MX_nfev':sm_at_mx.nfev,
                 'SM_Weinberg_MX_to_low_nfev':outsm.nfev,
                 'low_scale_results_NOT_PHYSICAL':last,
                 'mass_singular_values_relative_change':[float((flow[i]-fixed[i])/max(abs(fixed[i]),1e-250)) for i in range(3)],
                 'one_loop_SM_Weinberg_RGE_applied':True,
                 'correct_UV_and_EFT1_running_applied':False}
            print(f'{model}: direct/sequential PASS; independent loop max residual={error:.2g}; '
                  f'100 GeV SM-only reference masses [eV] '
                  +', '.join(f'{v:.6g}' for v in flow))
        else:
            row['SM_only_reference_running']=None
            print(f'{model}: direct/sequential PASS; integral residual={error:.2g}; SM-only proxy skipped')
        report['models'][stem]=row
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    # Human-readable final snapshot, with explicit certification status.
    lines=['# Workstream C v23 — consolidated B/C benchmark', '',
      '**Status: exploratory leading one-loop reference; NOT physical C5.**', '',
      '| Model | MX [GeV] | MY [GeV] | Fixed-order masses [eV] | SM-only 100 GeV reference masses [eV] |',
      '|---|---:|---:|---|---|']
    csvrows=['model,scale_GeV,kind,m1_eV,m2_eV,m3_eV']
    for stem,item in report['models'].items():
        a=item['full_mass_fixed_order_at_MX_NOT_PHYSICAL']['singular_masses_eV_ascending']
        b=(item['SM_only_reference_running']['low_scale_results_NOT_PHYSICAL']['singular_masses_eV_ascending']
           if item['SM_only_reference_running'] is not None else None)
        fa=', '.join(f'{x:.6g}' for x in a)
        fb=', '.join(f'{x:.6g}' for x in b) if b is not None else '(skipped)'
        lines.append(f"| {item['model']} | {item['MX_GeV']:.4f} | {item['MY_GeV']:.4f} | {fa} | {fb} |")
        csvrows.append(f"{item['model']},{item['MX_GeV']:.10g},fixed_order,"+','.join(f'{x:.16g}' for x in a))
        if b is not None:
            csvrows.append(f"{item['model']},{float(low):.10g},SM_only_proxy,"+','.join(f'{x:.16g}' for x in b))
    lines.extend(['', '## What passed', '',
       'Archived two-Yukawa one-loop topology; restricted MSbar hard/soft triangle subtraction; all-heavy Matchete versus sequential F/Y/X exact finite coefficient; independent Euclidean integral; conditional SM-only reference evolution.',
       '', '## What is still NOT established', ''])
    lines.extend('- '+x for x in report['remaining_physics_blockers'])
    lines.extend(['', '**Do not interpret these frozen-input numbers as measured or independently certified physical neutrino masses.**',''])
    (out.parent/'summary.md').write_text('\n'.join(lines),encoding='utf-8')
    (out.parent/'masses.csv').write_text('\n'.join(csvrows)+'\n',encoding='utf-8')
    print('Saved:',out)
    print('Saved:',out.parent/'summary.md')
    print('Saved:',out.parent/'masses.csv')
    print('PASS: final source-anchored leading one-loop reference; physical C5 NOT CERTIFIED')
    print('No further triangle/matching-reconstruction batches needed for rough numerical reference.')
    print('Remaining physical blockers documented in JSON; running is SM-only illustrative.')

def self_test():
    import numpy as np
    import mpmath as mp
    for a,b,c in ((100.,9.,4.),(100.,9.,9.),(4.,4.,4.),(1e10,1.4e6,1.1e6)):
        x=independent_triangle(a,b,c)
        y=master_closed(a,b,c)
        need(abs(x-y)<mp.mpf('1e-22')*max(1,abs(y)),'Scalar triangle integral wrong')
    toy=np.array([[1.+2j,3.-4j,0],[3.-4j,5.+6j,0],[0,0,7.]],complex)*1e-15
    stats=masses_from_c5(toy,246.22)
    need(len(stats['singular_masses_eV_ascending'])==3 and
         stats['singular_masses_eV_ascending'][0]>0,'Majorana/Takagi conversion failed')
    need(norm_relative(toy,toy.T)==0.,'Majorana coefficient must be symmetric')
    try:pair_matrix([[[1,0]]], 'tampered')
    except CompletionError:pass
    else:raise CompletionError('Malformed complex C5 was accepted')
    fake18={'format':'t3-bc-frozen-reference-three-fermion-triangle-v18','running_performed':False,'physical_C5':None}
    fake19={'format':'t3-bc-one-loop-hard-soft-msbar-restricted-v19','status':'RestrictedTriangleHardSoftMSbarSubtractionEstablished','full_T3_renormalized_matching_certified':False}
    fake21={'format':'t3-bc-one-loop-renormalizable-topology-closure-v21','one_loop_diagram_topology_complete_for_archived_UV_under_scope':True,'full_sequential_physical_matching_certified':False}
    fake22={'format':'t3-bc-independent-all-heavy-vs-sequential-finite-amplitude-v22','all_direct_vs_sequential_full_mass_agree':False,'full_running_performed':False,'independent_non_Matchete_UV_Feynman_rule_calculation_performed':False,'physical_C5':None}
    try:validate_audit_chain(Path('/does_not_exist'),fake18,fake19,fake21,fake22)
    except CompletionError:pass
    else:raise CompletionError('Failure in v22 source gate was accepted')
    print('PASS: independent convergent loop integration and stable closed master')
    print('PASS: source-scope gate, complex symmetric Majorana and mass conversion')
    print('PASS: physical-C5 certification cannot be asserted from these inputs')

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--output',type=Path,default=AUD/'t3_bc_consolidated_leading_reference_v23/summary.json')
    ap.add_argument('--low-scale',type=float,default=100.)
    ap.add_argument('--skip-running',action='store_true',help='Skip only the optional SM-only illustrative running')
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--require-physical',action='store_true')
    a=ap.parse_args()
    try:
        need(not a.require_physical,'Independent UV Feynman-rule normalization, complete EFT running and finite scheme are not yet certified')
        if a.self_test:self_test()
        else:
            root=a.root.resolve()
            sys.path.insert(0,str(root))
            dest=a.output if a.output.is_absolute() else root/a.output
            run(root,dest,a.low_scale,skip_running=a.skip_running)
        return 0
    except (CompletionError,OSError,ValueError,KeyError,IndexError,TypeError,
            ImportError,RuntimeError,AttributeError,ZeroDivisionError) as exc:
        print('FAIL (no physical C5):',exc,file=sys.stderr)
        return 2

if __name__=='__main__':raise SystemExit(main())
