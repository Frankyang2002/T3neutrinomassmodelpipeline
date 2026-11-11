"""Read-only Worker 1 scalar HHSS evidence and controlled saved-IVP comparison.

Never substitutes fitted HHSS signs/factors into the physical solver. Outputs to a
new directory only. Requires numpy and the existing Worker B audit JSON files.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

MODELS=('A_-4','D_-2')
AUDITS={
 'discrepancy':'output/audits/workerb_hhss_discrepancy_v1/report.json',
 'directions':'output/audits/workerb_hhss_direction_v1/report.json',
 'basis':'output/audits/workerb_basis_search_v1/report.json',
 'reconcile':'output/audits/workerb_batch_reconcile_v1/report.json',
 'source':'output/audits/workerb_source_hhss_v1/report.json',
}
CASES={
 'A_-4':'output/studies/general_scalar_extended/experimental_simultaneous_eft1_seven_term/simultaneous_A_m4_v2/A_m4/cases.json',
 'D_-2':'output/studies/general_scalar_extended/experimental_simultaneous_eft1_general_nonmixed/simultaneous_nonmixed13_v1/D_m2/cases.json',
}

def load(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def real_fit(source,candidate):
    s=np.asarray(source,dtype=float).ravel(); c=np.asarray(candidate,dtype=float).ravel()
    if s.shape!=c.shape or not np.all(np.isfinite(s)) or not np.all(np.isfinite(c)):
        raise ValueError('Incompatible or nonfinite tensor')
    n=float(np.linalg.norm(s));
    if n==0: return {'alpha':None,'residual':None}
    a=float(np.vdot(s,c).real/n**2)
    return {'alpha':a,'residual':float(np.linalg.norm(c-a*s)/n),
            'raw_relative':float(np.linalg.norm(c-s)/n)}

def complex_matrix(value):
    if not isinstance(value,dict) or 'real' not in value or 'imag' not in value:
        raise ValueError('Expected project complex matrix encoding')
    m=np.asarray(value['real'],float)+1j*np.asarray(value['imag'],float)
    if m.shape!=(3,3) or not np.all(np.isfinite(m)):
        raise ValueError('Invalid/nonfinite 3x3 C5')
    return m

def observables(c):
    # singular values of a complex symmetric Majorana matrix are Takagi masses.
    masses=np.linalg.svd(c,compute_uv=False)[::-1]
    norm=float(np.linalg.norm(c))
    return {'c5_norm':norm,'takagi_singular_values':masses.tolist(),
            'ratios_m1_m3_m2_m3':[float(x/masses[-1]) if masses[-1] else None for x in masses[:2]],
            'normalized_matrix':((c/norm) if norm else c)}

def compare_cases(base,alternative):
    b=load(base); a=load(alternative)
    if b.get('model')!=a.get('model'):raise ValueError('Different model cases')
    if not isinstance(b.get('cases'),dict) or not isinstance(a.get('cases'),dict):raise ValueError('Missing cases')
    needed={'1111','0111','1001','1110'}
    if not needed.issubset(b['cases']) or not needed.issubset(a['cases']):
        raise ValueError('Stage-freeze cases missing')
    output={}
    for key in sorted(set(b['cases'])&set(a['cases'])):
        bc=b['cases'][key];ac=a['cases'][key]
        cb=complex_matrix(bc['c5_at_low_scale']);ca=complex_matrix(ac['c5_at_low_scale'])
        nb=observables(cb);na=observables(ca)
        scale=real_fit(np.stack((cb.real,cb.imag)),np.stack((ca.real,ca.imag)))
        # Physical shape is invariant under a complex overall phase, unlike a real fit.
        overlap=abs(np.vdot(cb,ca))/(np.linalg.norm(cb)*np.linalg.norm(ca)) if np.linalg.norm(cb)*np.linalg.norm(ca) else None
        pb=bc.get('prediction',{});pa=ac.get('prediction',{})
        fields=['sin2_theta12','sin2_theta13','sin2_theta23','delta_m21_sq_ev2','delta_m3l_sq_ev2']
        changes={f:float(pa[f])-float(pb[f]) for f in fields if f in pb and f in pa}
        output[key]={'low_C5_real_scale_fit':scale,'complex_shape_overlap':float(overlap) if overlap is not None else None,
                     'base_takagi_mass_ratios':nb['ratios_m1_m3_m2_m3'],
                     'alternative_takagi_mass_ratios':na['ratios_m1_m3_m2_m3'],
                     'base_low_C5_norm':nb['c5_norm'],'alternative_low_C5_norm':na['c5_norm'],
                     'prediction_deltas':changes}
    # Fractions use the 1001 case (UV+SMEFT) as no-EFT1 reference, and 1111 as full.
    for label,d in [('baseline',b),('alternative',a)]:
        cfull=complex_matrix(d['cases']['1111']['c5_at_low_scale'])
        cfreeze=complex_matrix(d['cases']['1001']['c5_at_low_scale'])
        output[label+'_eft1_relative_contribution']=float(np.linalg.norm(cfull-cfreeze)/max(np.linalg.norm(cfull),1e-300))
    return output

def report(project, baseline=None,alternative=None, certified=False):
    raw={};missing=[]
    for key,path in AUDITS.items():
        p=project/path
        if p.is_file():raw[key]=load(p)
        else:missing.append(path)
    results={}
    for model in MODELS:
        d=raw.get('discrepancy',{}).get('results',{}).get(model,{})
        dr=raw.get('directions',{}).get('results',{}).get(model,{})
        bs=raw.get('basis',{}).get('results',{}).get(model,{})
        re=raw.get('reconcile',{}).get('results',{}).get(model,{})
        row={'model':model,'status':'saved_audit_available' if d and dr else 'needs_local_audits',
             'overall':d.get('overall'), 'species_blocks':d.get('heavy_species_blocks'),
             'direction_results':[{'coupling':x.get('coupling'),'axis':x.get('axis'),
                'scale':x.get('candidate_over_source_scale'),'residual':x.get('shape_residual_relative_to_source')}
                for x in dr.get('directions',[]) if x.get('source_norm',0)>0],
             'basis_best':bs.get('best_lambdaT3'),
             'basis_span_residual':bs.get('all_direction_span_residual_at_best_rotation'),
             'source_only_directions':re.get('source_only_directions'),
             'physical_normalization_certified':False}
        npz=project/'output/audits/workerb_batch_reconcile_v1'/f'{model}_reconciliation.npz'
        if npz.is_file():
            with np.load(npz,allow_pickle=False) as z:
                row['npz_keys']=list(z.files)
                if 'candidate_hhss' in z and 'source_hhss' in z:
                    c=z['candidate_hhss'];s=z['source_hhss']
                    row['independent_npz_fit']=real_fit(s,c)
                    split=2 if model=='A_-4' else 6
                    row['independent_npz_blocks']={f'{i}{j}':real_fit(s[:,:,si,sj],c[:,:,si,sj])
                        for i,si in [('S1',slice(0,split)),('S2',slice(split,8))]
                        for j,sj in [('S1',slice(0,split)),('S2',slice(split,8))]}
            row['reconciliation_npz_sha256']=digest(npz)
        results[model]=row
    if (baseline is None)!=(alternative is None):raise ValueError('Need both saved cases files')
    if baseline and not certified:raise ValueError('Comparison requires --physical-alternative-certified; do not compare unphysical fitted tensors')
    comparison=compare_cases(baseline,alternative) if baseline else None
    return {'schema':'worker1_scalar_robustness_v1','scope':'read-only; neither solver nor historical files changed',
            'audit_sha256':{key:digest(project/path) for key,path in AUDITS.items() if (project/path).is_file()},
            'missing_local_files':missing,'results':results,
            'physical_alternative_available':bool(comparison),
            'controlled_ivp_comparison':comparison,
            'conclusion':'Tensor structures are not globally proportional. Physical impact on flavour remains unmeasured without a certified alternative and paired coupled-IVP runs.'}

def markdown(r):
    lines=['# Worker 1 scalar quartic robustness — local diagnostic', '',
       '**Status:** diagnostic, not physical certification.', '',
       'EFT order: UV → all fermions out → EFT1 with both scalars → both scalars out → SMEFT.', '',
       '| Model | Raw discrepancy | Best real coefficient | Residual after fit |',
       '|---|---:|---:|---:|']
    for m,v in r['results'].items():
        q=v.get('overall') or {}
        lines.append(f"| {m} | {q.get('relative_to_source','n/a')} | {q.get('least_squares_candidate_over_source','n/a')} | {q.get('residual_after_best_real_scale_relative_to_source','n/a')} |")
        lines.extend(['',f'## {m} blocks','','| Block | Best real coefficient | Rescaled residual |','|---|---:|---:|'])
        for k,b in (v.get('species_blocks') or {}).items():
            lines.append(f"| {k} | {b.get('least_squares_candidate_over_source')} | {b.get('residual_after_best_real_scale_relative_to_source')} |")
        lines.extend(['','### Nonzero coupling directions','','| Coupling | Axis | Coefficient | Shape residual |','|---|---|---:|---:|'])
        for x in v.get('direction_results',[]):
            lines.append(f"| {x['coupling']} | {x['axis']} | {x['scale']} | {x['residual']} |")
    lines.extend(['','## Physical assessment','','A global sign/factor is ruled out as a full explanation. The lambdaT3 direction contains non-proportional structures. The native/RGBeta map and full scalar threshold are not certified; no fitted tensor was injected into the solver.'])
    if r['controlled_ivp_comparison'] is None:
        lines.append('**No certified physical alternative and paired identical-UV IVP runs were supplied; low-energy C5 shape, mixing, and stage fractions cannot be attributed to the tensor discrepancy.**')
    else: lines.append('A certified-alternative paired-results diagnostic is in report.json; verify UV input hashes and matching provenance independently before physical interpretation.')
    if r['missing_local_files']:lines.extend(['','## Missing files','']+['- `'+p+'`' for p in r['missing_local_files']])
    return '\n'.join(lines)+'\n'

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',type=Path,default=Path('.'))
    p.add_argument('--out-dir',type=Path,default=Path('output/audits/worker1_scalar_robustness_v1'))
    p.add_argument('--baseline-cases',type=Path)
    p.add_argument('--alternative-cases',type=Path)
    p.add_argument('--physical-alternative-certified',action='store_true')
    a=p.parse_args(argv)
    r=report(a.project,a.baseline_cases,a.alternative_cases,a.physical_alternative_certified)
    if a.out_dir.exists():raise FileExistsError(f'Refusing to overwrite {a.out_dir}; choose a fresh --out-dir')
    a.out_dir.mkdir(parents=True)
    (a.out_dir/'report.json').write_text(json.dumps(r,indent=2,default=str)+'\n',encoding='utf-8')
    (a.out_dir/'REPORT.md').write_text(markdown(r),encoding='utf-8')
    print('Wrote',a.out_dir/'report.json',a.out_dir/'REPORT.md')
    if r['missing_local_files']:print('Missing existing audit inputs:',*r['missing_local_files'],sep='\n  ')
if __name__=='__main__':main()
