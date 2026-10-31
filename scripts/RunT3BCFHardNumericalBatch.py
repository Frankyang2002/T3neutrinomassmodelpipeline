"""B/C source-anchored F-hard 68+68+68 formal finite sector at the running F threshold.

Numeric extraction only. This is NOT independent renormalised UV/EFT1 matching.
No formal coefficient is inserted into the EFT1 Wilson boundary or added to the
scalar bubble; this avoids double-counting an unknown scheme subtraction.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
from hashlib import sha256
from pathlib import Path
import traceback
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
MODELS={'B_m1':'B_-1','C_m1':'C_-1'}
FIVE=ROOT/'output/studies/t3_full_stack_overnight_v1/bc_full_mass'
AUDIT=ROOT/'output/studies/t3_bc_fermion_threshold_batch_v1/summary.json'
OUT=ROOT/'output/studies/t3_bc_f_hard_numerical_v1'
FORMAT='t3-bc-native-F-hard-three-majorana-formal-numerical-v1'
STATUS='FHardConstantAndLogNumericallyEvaluatedRenormalizationPending'

def require(cond,msg):
    if not cond:raise ValueError('F-hard running batch: '+msg)

def read(path):
    require(path.is_file(),'missing source '+str(path))
    o=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(o,dict),'invalid JSON '+str(path))
    return o

def digest(path):
    require(path.is_file(),'missing hash source '+str(path))
    return sha256(path.read_bytes()).hexdigest()

def external_script(path,alias):
    require(path.is_file(),'missing parser '+str(path))
    spec=importlib.util.spec_from_file_location(alias,path)
    require(spec is not None and spec.loader is not None,'cannot load parser')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def matrix_data(a):
    x=np.asarray(a,complex)
    require(x.shape==(3,3) and np.isfinite(x).all(),'invalid complex matrix')
    return {'real':x.real.tolist(),'imag':x.imag.tolist()}

def validate_audit(audit):
    require(audit.get('status')=='ConditionalUVToEFT1SourceBoundaryAudited' and
            audit.get('physical_C5_certified') is False,'fermion source audit incompatible')
    require(audit.get('EFT2_constructed') is False,'unexpected EFT structure')
    require(audit.get('formal_source_sha256'),'missing native formal-source digest')

def verify_polynomial(source_row,constructed):
    import sympy as sp
    for sector in ('pole','log','constant'):
        saved=source_row['sectors'][sector]['conventional_C5_2x2_without_hbar_MF']
        got=constructed[sector]['conventional_C5_2x2_without_hbar_MF']
        require(len(saved)==len(got)==2,'wrong source polynomial order')
        for i in range(2):
            for j in range(2):
                require(sp.expand(sp.sympify(saved[i][j])-sp.sympify(got[i][j]))==0,
                        f'native F-hard symbolic coefficient changed {sector}[{i},{j}]')

def one_model(stem, five_dir, audit, fformal, parser, v12, *,reference_dir=None):
    from scripts import RunT3BCSameTrajectoryWilson as same
    from Numerical.experimental.T3BCRunningOneScaleTriangle import uv_boundary
    from Numerical.experimental import T3BCFHardNumericalBoundary as numerics
    model=MODELS[stem]
    five_path=five_dir/f'{stem}_uv_to_eft1_coupled_five_mass_v1.json'
    five=read(five_path)
    same.check_five(five,model)
    boundary=uv_boundary(model,five)
    source=fformal['models'][stem]
    require(source['model']==model and source['pole_equals_MF_log_coefficient'] is True,
            'saved source pole/log consistency missing')
    require(source['F_to_EFT1_independent_counterterms_calculated'] is False and
            source['F_hard_finite_physical_matching_certified'] is False,
            'certification changed: re-audit')
    raw=ROOT/source['source_path']
    require(digest(raw)==source['source_sha256'],'raw native F-hard LLHH source changed')
    sectors=parser.split_sectors(raw.read_text(encoding='utf-8-sig'))
    require(parser.exact_pole_log_check(sectors)>0,'empty pole/log source')
    ledgers={};poly={}
    for key in ('pole','log','constant'):
        ledgers[key],counts=parser.make_ledger(sectors[key],model)
        require(counts=={'PR':34,'PL':34},'wrong native chirality inventory')
        poly[key]=parser.construct_polynomials(ledgers[key],v12)
    verify_polynomial(source,poly)
    result,scan=numerics.scale_scan(ledgers,masses=boundary['mf'],
         y1=boundary['y1'],y2=boundary['y2'],U=boundary['scalar_basis'].U,
         quartics=boundary['quartics'],mu_f=boundary['mu_gev'])
    reference=None
    if reference_dir is not None:
        refpath=reference_dir/f'{stem}_same_trajectory_wilson.json'
        if refpath.exists():
            ref=read(refpath)
            require(ref.get('model')==model and ref.get('physical_C5_certified') is False,
                    'same-trajectory reference incompatible')
            refmatrix=ref.get('fermion_grouped_hard_C5_candidate')
            if refmatrix:
                from Numerical.experimental.T3BCRunningOneScaleTriangle import matrix
                c=matrix(refmatrix,'reference F grouped C5')
                reference={'reference_path':str(refpath.relative_to(ROOT)),
                    'reference_sha256':digest(refpath),
                    'formal_over_prior_grouped_hard_norm':float(np.linalg.norm(result['formal_at_mu_f'])/max(np.linalg.norm(c),1e-35)),
                    'both_quantities_are_candidate_in_different_subtractions':True,
                    'not_added_or_subtracted':True}
    def sector(name):
        arr=result['sector_matrix'][name]
        return {'matrix_GeV_inverse':matrix_data(arr),
                'frobenius_norm_GeV_inverse':float(np.linalg.norm(arr)),
                'by_heavy_state_frobenius_norm_GeV_inverse':
                  [float(np.linalg.norm(x)) for x in result['per_heavy_sector'][name]]}
    return {'format':FORMAT,'status':STATUS,'model':model,
        'UV_to_EFT1_to_SM_only_one_intermediate_EFT':True,
        'running_five_mass_source_sha256':digest(five_path),
        'native_F_hard_raw_source_sha256':digest(raw),
        'Majorana_eigenmasses_GeV':boundary['mf'].tolist(),
        'matching_scale_GeV':boundary['mu_gev'],
        'scalar_eigenmasses_at_F_GeV':np.sqrt(boundary['scalar_basis'].eigenvalues_sq).tolist(),
        'source_inventory':{'pole':68,'log':68,'constant':68},
        'formal_native_source_sectors':{k:sector(k) for k in ('pole','log','constant')},
        'log_muF_squared_over_Mr_squared':result['matching_scale_log_per_heavy'].tolist(),
        'formal_constant_plus_log_at_grouped_F_GeV_inverse':matrix_data(result['formal_at_mu_f']),
        'formal_d_dln_muF_at_grouped_F_GeV_inverse':matrix_data(result['derivative_dln_mu_formal']),
        'formal_matching_scale_scan':[{**{k:v for k,v in row.items() if k!='matrix'},
                                      'matrix_GeV_inverse':matrix_data(row['matrix'])} for row in scan],
        'numeric_checks':{'pole_log_relative':result['pole_log_matrix_relative_residual'],
                          'maximum_native_PL_PR_conjugacy_relative':result['max_PR_PL_conjugacy_relative_residual'],
                          'dlog_scale_derivative_relative':result['finite_difference_log_derivative_relative_residual']},
        'optional_previous_grouped_boundary_comparison':reference,
        'corrected_EFT1_Wilson_boundary_computed':False,
        'UV_IR_pole_assignment_computed':False,
        'independent_UV_EFT1_counterterms_computed':False,
        'physical_F_hard_finite_matching':False,
        'physical_C5':None,
        'notes':'This is the source polynomial evaluated at the running UV F threshold; constant + MF logarithm is formal, NOT a renormalized difference of UV and EFT1 amplitudes. No numeric term is added to the baseline Wilson boundary.'}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--five-mass-dir',type=Path,default=FIVE)
    ap.add_argument('--audit',type=Path,default=AUDIT)
    ap.add_argument('--output-dir',type=Path,default=OUT)
    ap.add_argument('--no-previous-reference',action='store_true')
    args=ap.parse_args()
    audit=read(args.audit);validate_audit(audit)
    fpath=args.audit.parent/'f_hard_formal_source.json'
    require(digest(fpath)==audit['formal_source_sha256'],'F-hard source report changed')
    formal=read(fpath)
    require(formal.get('full_physical_matching_certified') is False and formal.get('physical_C5') is None,
            'native formal source flags changed')
    parser=external_script(ROOT/'scripts/BuildT3BCFHardBoundaryInterface.experimental.py','bc_f_hard_native_parser')
    v12=parser.import_v12(ROOT)
    rows={}
    for stem,model in MODELS.items():
        data=one_model(stem,args.five_mass_dir,audit,formal,parser,v12,
            reference_dir=None if args.no_previous_reference else
                  ROOT/'output/studies/t3_bc_same_trajectory_wilson_v1')
        rows[model]=data
        args.output_dir.mkdir(parents=True,exist_ok=True)
        file=args.output_dir/f'{stem}_formal_F_hard_numeric.json'
        file.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        norm=data['formal_native_source_sectors']['constant']['frobenius_norm_GeV_inverse']
        print(model,'formal F-hard finite constant ||C||=',f'{norm:.9e}',
              'pole/log=',f"{data['numeric_checks']['pole_log_relative']:.3g}",
              'counterterms=MISSING',flush=True)
    summary={'format':FORMAT,'status':'NativeFHardThreeMassNumericSourceComputedCountertermsPending',
             'models':{k:{'formal_constant_norm_GeV_inverse':v['formal_native_source_sectors']['constant']['frobenius_norm_GeV_inverse'],
                          'formal_at_grouped_F_norm_GeV_inverse':float(np.linalg.norm(
                              np.asarray(v['formal_constant_plus_log_at_grouped_F_GeV_inverse']['real'])+
                              1j*np.asarray(v['formal_constant_plus_log_at_grouped_F_GeV_inverse']['imag']))),
                          'native_checks':v['numeric_checks']} for k,v in rows.items()},
             'one_intermediate_EFT_only':True,
             'UV_EFT1_counterterms_calculated':False,'full_physical_matching_certified':False,
             'physical_C5':None,'production_files_modified':False}
    dest=args.output_dir/'summary.json'
    dest.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('SUMMARY',summary['status'],'physical C5=NOT CERTIFIED',dest)
    return 0

if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,KeyError,TypeError,OSError,RuntimeError,ImportError) as exc:
        print('F-hard batch BLOCKED:',exc,flush=True)
        traceback.print_exc()
        raise SystemExit(2)
