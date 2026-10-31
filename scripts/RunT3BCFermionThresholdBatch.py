"""B/C F-first UV->EFT1 threshold handoff, preserving uncertified source sectors.

Reuses existing native Matchete tree/one-loop reports and established source
F-hard projector. Does not invent UV/IR counterterms or promote formal finite
terms to renormalized Wilson coefficients.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'output/studies/t3_bc_fermion_threshold_batch_v1'
MASS = ROOT/'output/general_t3_mass_eigen_matchete'
STEMS = {'B_m1':'B_-1','C_m1':'C_-1'}


def require(ok, message):
    if not ok:
        raise ValueError('B/C F-threshold: '+message)


def read(path):
    require(path.is_file(), 'missing '+str(path))
    data=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data,dict),'not a JSON object '+str(path))
    return data


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_script(path,name):
    require(path.is_file(),'missing source '+str(path))
    spec=importlib.util.spec_from_file_location(name,str(path))
    require(spec is not None and spec.loader is not None,'cannot import '+str(path))
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_one(stem,tree,loop,source_dir):
    model=STEMS[stem]
    require(tree.get('model')==model and loop.get('model')==model,'model mismatch')
    require(tree.get('status')=='FermionFirstTreeLLSSExtractedNotSequentialEFT1','tree matching not passed')
    require(loop.get('status')=='FermionFirstOneLoopRawSectorsExtractedNotRenormalized','one-loop extraction not passed')
    for data in (tree,loop):
        require(data.get('fermion_field_absent') is True and data.get('x_scalar_present') is True and data.get('y_scalar_present') is True,
                'incorrect EFT1 active fields')
        require(data.get('physical_matching_performed') is False and data.get('matchete_rgbeta_quartic_normalization_certified') is False,
                'upstream physical status changed: review required')
    require(tree.get('tree_matching_completed') is True and loop.get('one_loop_matching_completed') is True,
            'native matching did not complete')
    require(loop.get('raw_hard_llhh_scheme_renormalized') is False and loop.get('raw_hard_llhh_pole_rge_certified') is False,
            'native one-loop certificate changed')
    require(int(tree.get('llss_term_count',0))>0 and int(loop.get('llss_term_count',0))>0,
            'missing LLSS operators')
    return {'model':model,'status':'NativeTreeAndRawOneLoopVerifiedCountertermsPending',
            'tree_sha256':digest(source_dir[0]),'one_loop_sha256':digest(source_dir[1]),
            'tree_llss_terms':int(tree['llss_term_count']),
            'one_loop_raw_llhh_terms':int(loop['raw_hard_llhh_candidate_term_count']),
            'EFT1_active_fields':['SM','S1','S2'],
            'F_integrated_out':True,
            'EFT1_Wilson_boundary_renormalized':False,
            'physical_C5_certified':False}


def check_formal(report):
    require(report.get('format')=='t3-bc-F-hard-formal-boundary-plus-Y-conditional-ledger-v1',
            'formal source report format mismatch')
    require(report.get('full_physical_matching_certified') is False and report.get('physical_C5') is None,
            'unexpected physical certificate')
    require(report.get('F_hard_constant_and_Y_conditional_are_SEPARATE_boundaries') is True,
            'F/Y double-counting protection missing')
    result={}
    for stem,model in STEMS.items():
        row=report['models'][stem]
        require(row['model']==model and row['raw_sector_counts']=={'pole':68,'log':68,'constant':68},
                'native F-hard term inventory changed')
        require(row['pole_equals_MF_log_coefficient'] is True and
                row['F_to_EFT1_independent_counterterms_calculated'] is False and
                row['F_hard_finite_physical_matching_certified'] is False,
                'F-hard pole/log or counterterm status changed')
        for sector in ('pole','log','constant'):
            entries=row['sectors'][sector]
            matrix=entries['conventional_C5_2x2_without_hbar_MF']
            require(len(matrix)==2 and all(len(r)==2 for r in matrix),
                    'incomplete F-hard polynomial matrix')
        result[model]={'source_sha256':row['source_sha256'],
                       'pole_log_constant_terms':row['raw_sector_counts'],
                       'pole_equals_mass_log':True,
                       'formal_hard_boundary_polynomials':row['sectors'],
                       'counterterms_calculated':False,
                       'renormalized_finite_boundary':False}
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=OUT)
    parser.add_argument('--skip-formal-extraction',action='store_true',
                        help='Use previously generated F-hard report in selected output directory')
    ns=parser.parse_args()
    sys.path.insert(0,str(ROOT))
    from scripts import RunGeneralT3FermionFirstEFT1Tree as t
    from scripts import RunGeneralT3FermionFirstEFT1OneLoop as l
    rows={}
    for stem,model in STEMS.items():
        tree_stem=t._stem(model,'final')
        loop_stem=l._stem(model,'final')
        tp=t.OUTDIR/(tree_stem+'.json')
        lp=l.OUTDIR/(loop_stem+'.json')
        tree,loop=read(tp),read(lp)
        t.validate_tree_report(tree,model,t.OUTDIR,tree_stem)
        l.validate_one_loop_report(loop,model,l.OUTDIR,loop_stem)
        rows[model]=validate_one(stem,tree,loop,(tp,lp))
    output=ns.output_dir.resolve()
    output.mkdir(parents=True,exist_ok=True)
    formal_path=output/'f_hard_formal_source.json'
    if not ns.skip_formal_extraction:
        script=ROOT/'scripts/BuildT3BCFHardBoundaryInterface.experimental.py'
        require(script.is_file(),'missing verified F-hard native source parser')
        cmd=[sys.executable,str(script),'--root',str(ROOT),'--output',str(formal_path)]
        proc=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,errors='replace')
        (output/'f_hard_extraction.log').write_text(proc.stdout+'\nSTDERR:\n'+proc.stderr,encoding='utf-8')
        require(proc.returncode==0,'native F-hard extraction blocked; see '+str(output/'f_hard_extraction.log'))
    formal=check_formal(read(formal_path))
    for model in rows:
        rows[model]['F_hard_source']=formal[model]
        rows[model]['status']='NativeTreeAndFHardSectorsExtractedRenormalizedBoundaryPending'
    report={'format':'t3-bc-fermion-threshold-handoff-batch-v1',
            'status':'ConditionalUVToEFT1SourceBoundaryAudited',
            'stage':'UV -> EFT1 (integrate F; scalars remain)',
            'models':rows,
            'formal_source_sha256':digest(formal_path),
            'missing_for_physical_matching':['UV/IR pole separation',
              'independent UV and EFT1 counterterm subtraction',
              'cross-package Yukawa/quartic normalization certificate',
              'finite field and parameter/EOM shifts',
              'consistent running and scalar threshold matching'],
            'EFT2_constructed':False,'physical_C5_certified':False,
            'production_files_modified':False}
    target=output/'summary.json'
    target.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for model,row in rows.items():
        print(model,'tree LLSS',row['tree_llss_terms'],'raw LLHH',row['one_loop_raw_llhh_terms'],
              'F hard 68+68+68 PASS; finite boundary NOT CERTIFIED')
    print(report['status'])
    print('Output:',target)
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,KeyError,TypeError,OSError,ImportError) as exc:
        print(type(exc).__name__+': '+str(exc),file=sys.stderr)
        raise SystemExit(2)
