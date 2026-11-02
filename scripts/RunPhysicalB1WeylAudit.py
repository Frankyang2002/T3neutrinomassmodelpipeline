#!/usr/bin/env python3
"""Read-only B1 all-16 physical-normalisation audit, explicitly fail closed.

Run in repository root. Reads existing native CG audits, benchmarks and
provisional Weyl inputs; does not invoke the old component builder. Does NOT
promote an agreement with the old provisional builder to source certification.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import traceback
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.physical_normalisation.b1_matchete_weyl import (
    B1Error,model_specs,complex_array,construct_components,four_spinor_witness,
    spurion_charge,relative,require,encode)

BASE=Path('output/general_observable/physical_normalisation/b1_02')
COMPONENT=Path('Reports/output/full/multiplicity/component_support_audit.json')
GENERATOR=Path('Reports/output/full/multiplicity/matchete_generator_analysis.json')
COVARIANCE=Path('Reports/output/full/multiplicity/cg_covariance_audit.json')
FROZEN=Path('output/general_observable/lo_direct_one_loop/provisional_inputs')
BENCHMARK=Path('Numerical/benchmarks/fitted_reference')
# GitHub main at commit 0b210f18... fetched read-only on 2026-11-03.
# Git blob IDs protect against silently assuming different local source rules.
SOURCE_GIT_BLOBS={
 'Lagrangian/interactions/T3Topology.wl':'a45a099a8b766cbad0114a0559cb0bcb80895649',
 'Lagrangian/model/T3Fields.wl':'65467ede6ef06a7743e3ca2878debf777dc94507',
 'Lagrangian/model/LagrangianBuilder.wl':'98d549cda503b751380dc10c4a56212e5312ccd8',
 'Lagrangian/model/SU2Invariants.wl':'dbe0cd490204dec363fcf4278f079b58c8ddcefc',
 'scripts/BuildT3DirectOneLoopInputs.py':'976a5c6303c637b7d1e45942a8dfa06f7e780ac7',
 'Numerical/direct_one_loop/core.py':'abbf261fbf2ca2561482d3e53b84452d1654610e',
 'RGE/running/rgbeta/T3RGBetaModel.wl':'7af7a4ddf99a7e38210543a21b3b1cda203eaf02',
}
HEAD='0b210f18f6a2cfb10fff03059e2ec4a73d2195a4'


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def _git_blob_sha1(data: bytes) -> str:
    """Compute the Git blob object ID for the *supplied bytes*."""
    header = b'blob ' + str(len(data)).encode('ascii') + b'\0'
    return hashlib.sha1(header + data).hexdigest()


def gitblob(p, *, normalize_eol: bool = False):
    """Compute a source blob ID; optionally permit CRLF -> LF only.

    Git checks out LF-tracked files as CRLF under common Windows settings.
    Normalization here is deliberately restricted to line endings: no text
    parsing, whitespace removal, Git filters, or content modifications.
    """
    data = p.read_bytes()
    if normalize_eol:
        # A bare CR is not normal Windows CRLF conversion; fail closed.
        require(b'\r' not in data.replace(b'\r\n', b''),
                'UNSUPPORTED SOURCE DRIFT: lone CR byte in '+str(p))
        data = data.replace(b'\r\n', b'\n')
    return _git_blob_sha1(data)


def read(p):
    require(p.is_file(),'Missing local input: '+str(p))
    return json.loads(p.read_text(encoding='utf-8-sig'))
def git_revision(root):
    try:
        process=subprocess.run(['git','-C',str(root),'rev-parse','HEAD'],capture_output=True,text=True,timeout=6)
        return process.stdout.strip() if process.returncode==0 else None
    except (OSError,subprocess.TimeoutExpired):return None


def guard_source(root):
    checks={}
    for rel,expected in SOURCE_GIT_BLOBS.items():
        p=root/rel
        require(p.is_file(),'Missing local source '+rel)
        raw_blob=gitblob(p)
        canonical_blob=gitblob(p, normalize_eol=True)
        match=canonical_blob==expected
        checks[rel]={
            'git_blob_sha1':canonical_blob,
            'working_tree_raw_git_blob_sha1':raw_blob,
            'sha256':digest(p),
            'github_main_blob_match':match,
            'crlf_normalization_applied':raw_blob!=canonical_blob,
            'source_bytes_match_after_line_ending_normalization':match,
        }
        require(match,'UNSUPPORTED SOURCE DRIFT in '+rel+
                '; changes beyond CRLF/LF conversion rejected. Upload actual local file for updated audit')
    return checks


def rotations(source):
    require(not source.get('errors'),'Invalid generator audit')
    us={1:np.eye(1,dtype=complex)}
    for rec in source['representations']:
        require(rec.get('status')=='generator_algebra_verified','Generator audit not passed')
        d=int(rec['dimension'])
        vals=rec['U_columns_are_weight_eigenvectors']
        u=np.array([[complex(sp.N(sp.sympify(s,locals={'I':sp.I}),18)) for s in row] for row in vals],complex)
        require(u.shape==(d,d) and np.linalg.norm(u.conj().T@u-np.eye(d))<2e-10,'Bad generator rotation')
        us[d]=u
    require(set(us)=={1,2,3},'Missing SU2 rotations')
    return us


def name_for_benchmark(short, dims):
    d1,d2,df,a=dims
    return f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{"m" if a<0 else "p"}{abs(a)}'


def current_native(native_dir):
    if native_dir is None:
        return {'provided':False,'observed':False,'reason':'No fresh native Matchete output supplied'}
    out={}
    for tag in ('A_p0','B_m1','C_m1'):
        path=native_dir/(tag+'.json')
        rec=read(path)
        require(rec.get('format')=='t3-b1-native-restricted-action-v1' and
                rec.get('model')==tag and
                rec.get('status')=='NativeActionTermsIsolated_PhysicalDictionaryPending',
                'Unexpected native output schema/status: '+str(path))
        require(rec.get('physical_source_certified') is False and
                rec.get('source_gammaCC_concrete_matrix_measured') is False,
                'False native physical-certificate promotion: '+str(path))
        checks=rec.get('tests',{})
        require(checks and all(val is True for val in checks.values()),'Native checks not passed')
        for y in ('y1','y2'):
            require(rec['yukawas'][y]['original_term_count']>=1 and
                    rec['yukawas'][y]['hc_term_count']>=1,'Missing genuine Yukawa/HC decomposition')
        out[tag]={'native_terms_observed':True,'sha256':digest(path),
                  'mass_terms':rec.get('mass_term_count'),
                  'yukawa_original_hc_counts':{y:(rec['yukawas'][y]['original_term_count'],
                                                  rec['yukawas'][y]['hc_term_count']) for y in ('y1','y2')}}
    return {'provided':True,'observed':True,'representatives':out,
            'absolute_GammaCC_phase_measured':False}


def compare_one(short,dims,cases,bench,us,inp):
    b1=construct_components(cases[dims],bench,us,short)
    gotM=complex_array(inp['fermion_mass_gev'])
    gotY=np.zeros_like(b1['Y']);gotZ=np.zeros_like(b1['Z'])
    require(len(inp['vertices'])==2,'Restricted model should have two vertices')
    entries={x['name']:x for x in inp['vertices']}
    require(set(entries)=={'y1_HC','y2_HC'},'Vertex labels/orientations changed')
    gotY=complex_array(entries['y2_HC']['field']);gotZ=complex_array(entries['y1_HC']['conjugate'])
    require(gotY.shape==b1['Y'].shape and gotZ.shape==b1['Z'].shape and gotM.shape==b1['M'].shape,'Unexpected frozen shape')
    mass_res=relative(gotM,b1['M'])
    # Production candidate uses + conj(y1,y2); the mathematically derived
    # L_Y has - conj(y1,y2) in this explicitly chosen four-spinor convention.
    # Simultaneous minus on Y,Z is a common external-neutrino rephasing.
    rephaseY=relative(gotY,-b1['Y']);rephaseZ=relative(gotZ,-b1['Z'])
    fixedY=relative(gotY,b1['Y']);fixedZ=relative(gotZ,b1['Z'])
    require(max(mass_res,rephaseY,rephaseZ)<2.e-9,short+': frozen Weyl inputs differ beyond common external field rephase')
    return {'model':short,'fermion_type':b1['fermion_type'],
            'n_heavy_weyl':b1['M'].shape[0],
            'shape_M':list(b1['M'].shape),'shape_YZ':list(b1['Y'].shape),
            'residuals':{'mass':mass_res,'Y_after_global_minus':rephaseY,
                         'Z_after_global_minus':rephaseZ,'Y_fixed_basis':fixedY,
                         'Z_fixed_basis':fixedZ},
            'global_external_neutrino_rephase_is_harmless':True,
            'triplet_metric':encode(b1['metric']) if dims[2]==3 else None,
            'physically_verified':False,
            'pending_B2_basis':True}


def audit(root,native_dir):
    report={'format':'t3-b1-matchete-to-weyl-v1',
      'github_reference_commit':HEAD,'local_git_revision':git_revision(root),
      'execution_utc':datetime.now(timezone.utc).isoformat(),
      'physical_source_certified':False,'physical_B1_scope_complete':False,
      'native_action_observed':False,'relative_physical_factors_verified':False,
      'needs_B2_SU2_dictionary':True,'models':[],'checks':{}}
    report['checks']['source_blobs']=guard_source(root)
    report['checks']['spinor_identity']=four_spinor_witness()
    report['checks']['spurion_dirac']=spurion_charge(False)
    report['checks']['spurion_majorana']=spurion_charge(True)
    report['weyl_rules_mathematically_verified']=True
    native=current_native(native_dir)
    report['native']=native
    report['native_action_observed']=native['observed']
    component=read(root/COMPONENT);gen=read(root/GENERATOR);cov=read(root/COVARIANCE)
    require(not component.get('errors') and len(component.get('cases',[]))==16,'Bad all-16 CG source audit')
    require(cov.get('status')=='cg_generator_covariance_verified','CG Ward audit not accepted')
    us=rotations(gen)
    cases={tuple(int(rec['model'][k]) for k in ('dS1','dS2','dF','alpha')):rec for rec in component['cases']}
    valid={tuple(int(rec['model'][k]) for k in ('dS1','dS2','dF','alpha')):rec.get('all_covariant') for rec in cov['models']}
    require(len(cases)==16 and len(valid)==16 and all(valid.values()),'Native CG covariance unsupported')
    for short,dims in model_specs():
        benchpath=root/BENCHMARK/(name_for_benchmark(short,dims)+'.json')
        frozenpath=root/FROZEN/(short+'.json')
        bench=read(benchpath);inp=read(frozenpath)
        require(inp.get('model')==short,'Frozen model name mismatch')
        require(inp.get('provenance',{}).get('benchmark_sha256')==digest(benchpath),
                'Frozen benchmark provenance mismatch: '+short)
        require(abs(float(inp['renormalization_scale_gev'])-float(bench['base_state']['mu_gev']))<1e-8,
                'Scale mismatch: '+short)
        item=compare_one(short,dims,cases,bench,us,inp)
        item['benchmark_sha256']=digest(benchpath);item['frozen_input_sha256']=digest(frozenpath)
        report['models'].append(item)
    report['status']='All16ConditionalCandidateWeylTensorsReconstructed_PhysicalDictionaryPending'
    report['count_conditional_comparison_passed']=len(report['models'])
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--native-dir',type=Path,default=None)
    p.add_argument('--output-dir',type=Path,default=ROOT/BASE)
    p.add_argument('--self-test',action='store_true')
    args=p.parse_args()
    if args.self_test:
        four_spinor_witness();spurion_charge(False);spurion_charge(True)
        print('B1 chosen-basis Grassmann and spurion identities: PASS; native proof: NOT RUN')
        return 0
    try:
        result=audit(args.project_root,args.native_dir)
        exitcode=0
    except Exception as exc:
        result={'status':'Blocked','error':str(exc),'error_type':type(exc).__name__,
                'github_reference_commit':HEAD,'local_git_revision':git_revision(args.project_root),
                'model_results_not_certified':True,'physical_source_certified':False,
                'native_action_observed':False,'weyl_rules_mathematically_verified':False,
                'relative_physical_factors_verified':False,'needs_B2_SU2_dictionary':True,
                'physical_B1_scope_complete':False}
        exitcode=2
    args.output_dir.mkdir(parents=True,exist_ok=True)
    output=args.output_dir/'all16_component_check.json'
    try:
        with output.open('x',encoding='utf-8') as fd:
            json.dump(result,fd,indent=2,sort_keys=True)
            fd.write('\n')
    except FileExistsError:
        print('REFUSING TO OVERWRITE: '+str(output),file=sys.stderr)
        return 3
    print(result['status'],'models:',len(result.get('models',[])),'physical_source_certified: FALSE')
    if exitcode:print(result.get('error',''),file=sys.stderr)
    return exitcode

if __name__=='__main__':raise SystemExit(main())
