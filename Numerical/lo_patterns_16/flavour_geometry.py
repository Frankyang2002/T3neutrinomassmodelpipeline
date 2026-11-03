"""Chat A Batch 4: invariant heavy-eigenspace and neutrino-mode geometry.

Read-only consumer of saved Batch 1/2/3, with NO UV solver calls. All
interpretations are conditional on the provisional source Weyl dictionary.
Distinct heavy singular-mass eigenspaces are invariant; individual vectors
inside a degenerate eigenspace are not and are deliberately never reported.
"""
from __future__ import annotations

import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary
from Numerical.lo_patterns_16.result_analysis import _source_hash_guard

BASE=Path('output/general_observable/lo_pattern_study16')
STUDY=BASE/'provisional_01'
BATCH02=BASE/'analysis_batch02_01'
BATCH03=BASE/'analysis_batch03_01'
TARGET=BASE/'analysis_batch04_01'
ZERO=1e-100


class FlavourError(ValueError):
    """A source contract, analytical identity or write-safety gate failed."""


def require(test:bool,message:str)->None:
    if not test:raise FlavourError(message)


def load(path:Path)->dict:
    require(path.is_file(),f'Missing source: {path}')
    return json.loads(path.read_text(encoding='utf-8-sig'))


def digest(path:Path)->str:
    return sha256(path.read_bytes()).hexdigest()


def matrix(data:dict,label:str)->np.ndarray:
    """Decode original complex-array representation, reject nonfinite matrices."""
    require(isinstance(data,dict) and 'real' in data and 'imag' in data,label+': invalid encoded complex matrix')
    a=np.asarray(data['real'],dtype=float);b=np.asarray(data['imag'],dtype=float)
    require(a.shape==(3,3) and b.shape==(3,3),label+': expected 3x3 complex matrix')
    require(np.all(np.isfinite(a)) and np.all(np.isfinite(b)),label+': nonfinite entries')
    return a+1j*b


def encode(z:complex)->dict:
    return {'real':float(z.real),'imag':float(z.imag)}


def norm(a:np.ndarray)->float:
    return float(np.linalg.norm(a))


def sv(a:np.ndarray)->np.ndarray:
    """Ascending neutrino masses (eV), never eigenvalues of complex mnu."""
    return np.sort(np.linalg.svd(a,compute_uv=False))


def near(a:np.ndarray,b:np.ndarray,tol:float=1e-8)->bool:
    return norm(a-b)<=tol*max(norm(a),norm(b),ZERO)


def geometry(total:np.ndarray,parts:list[dict],model:str)->dict:
    """Exact singular-vector mode decomposition of a sum of heavy-group matrices.

    For SVD M=U diag(s) V^dag, a_gj=u_j^dag M_g v_j, so
      s_j=sum_g a_gj  (real positive to numerical accuracy).
    Writing G=sum_g ||M_g||F gives the *identity*
      s_j/G = [sum_g |a_gj|/G] [s_j/sum_g |a_gj|].
    First factor diagnoses poor mode projection (geometry); second diagnoses
    destructive complex interference. It is not a fit or approximation.
    """
    require(total.shape==(3,3) and np.all(np.isfinite(total)),model+': invalid mass matrix')
    require(len(parts)>0,model+': empty heavy-mass eigenspace list')
    full=sum((p['matrix_ev'] for p in parts),np.zeros((3,3),complex))
    strength=sum(norm(p['matrix_ev']) for p in parts)
    require(strength>0,model+': all heavy groups are zero')
    require(near(full,total,2e-7),model+': heavy-eigenspace sum does not reproduce full mass matrix')
    require(norm(total-total.T)<1e-7*max(norm(total),ZERO),model+': non-symmetric neutrino matrix')
    for p in parts:
        require(norm(p['matrix_ev']-p['matrix_ev'].T)<1e-7*max(norm(p['matrix_ev']),ZERO),
                model+': a heavy-group matrix is not symmetric')
        require(p['mass_gev']>0 and p['multiplicity']>0,model+': invalid heavy group')
    U,S,Vh=np.linalg.svd(total)
    V=Vh.conj().T
    order=np.argsort(S)
    modes=[]; projections=[]
    for i,idx in enumerate(order):
        u=U[:,idx];v=V[:,idx];mass=float(S[idx])
        amps=[complex(np.vdot(u,p['matrix_ev']@v)) for p in parts]
        amp_sum=sum(amps,0j)
        raw=sum(abs(a) for a in amps)
        require(abs(amp_sum-mass)<3e-8*max(strength,1e-25),
                model+f': failed projected-mode identity for m{i+1}')
        alignment=raw/strength
        coherence=mass/raw if raw>1e-18*strength else None
        require(alignment<=1+3e-8,model+': projection exceeds Frobenius bound')
        if coherence is not None:require(coherence<=1+3e-8,model+': interference ratio exceeds unity')
        identity=mass/strength
        require(abs(alignment*(coherence or 0)-identity)<3e-8,
                model+': geometry-times-interference identity failed')
        if coherence is None:
            mechanism='unresolved_zero_mode'
        elif alignment<coherence/2:
            mechanism='mode_geometry_dominant'
        elif coherence<alignment/2:
            mechanism='heavy_eigenspace_cancellation_dominant'
        else:
            mechanism='mixed_or_neither_dominant'
        dominant=int(np.argmax([abs(z) for z in amps]))
        modes.append({'mode':i+1,'mass_ev':mass,'projected_strength_ev':raw,
                      'group_frobenius_strength_ev':strength,
                      'geometric_projection_fraction':alignment,
                      'coherence_fraction':coherence,
                      'light_mass_fraction_of_group_strength':identity,
                      'product_identity_absolute_error':abs(alignment*(coherence or 0)-identity),
                      'dominant_heavy_group_index':dominant+1,
                      'dominant_group_projection_fraction':abs(amps[dominant])/max(raw,ZERO),
                      'mechanism_diagnostic':mechanism})
        projections.extend({'model':model,'mode':i+1,'heavy_group_index':g+1,
                            'mass_gev':p['mass_gev'],'amplitude_ev':encode(a),
                            'amplitude_real_ev':float(a.real),'amplitude_imag_ev':float(a.imag),
                            'amplitude_absolute_ev':abs(a),
                            'amplitude_fraction_of_projected_strength':abs(a)/max(raw,ZERO)}
                           for g,(p,a) in enumerate(zip(parts,amps)))
    group_rows=[];leave_out=[];cross=[]
    for g,p in enumerate(parts):
        H=p['matrix_ev'];value=norm(H)
        singular=sv(H)
        rank=int(np.sum(singular>1e-8*max(singular[-1],1e-30)))
        group_rows.append({'model':model,'heavy_group_index':g+1,
                           'fermion_mass_gev':p['mass_gev'],
                           'eigenspace_multiplicity':p['multiplicity'],
                           'matrix_frobenius_ev':value,'frobenius_fraction_of_group_strength':value/strength,
                           'relative_numerical_rank':rank,'singular_masses_ev':singular.tolist()})
        remaining=total-H
        leave=sv(remaining)
        leave_out.append({'model':model,'omitted_heavy_group_index':g+1,
                          'fermion_mass_gev':p['mass_gev'],
                          'm1_without_ev':float(leave[0]),'m2_without_ev':float(leave[1]),
                          'm3_without_ev':float(leave[2]),
                          'm1_change_ev':float(leave[0]-modes[0]['mass_ev']),
                          'm2_change_ev':float(leave[1]-modes[1]['mass_ev']),
                          'm3_change_ev':float(leave[2]-modes[2]['mass_ev']),
                          'matrix_fraction_removed':value/max(norm(total),ZERO)})
        for h in range(g+1,len(parts)):
            J=parts[h]['matrix_ev']
            dot=np.vdot(H,J)
            cosine=float(dot.real/max(norm(H)*norm(J),ZERO))
            cross.append({'model':model,'group_i':g+1,'group_j':h+1,
                          'complex_overlap_ev2':encode(dot),
                          'frobenius_overlap_cosine':max(-1.,min(1.,cosine))})
    coh=norm(total)/strength
    singular=sv(total)
    det_from_svd=float(np.prod(singular))
    determinant=float(abs(np.linalg.det(total)))
    require(abs(det_from_svd-determinant)<1e-7*max(det_from_svd,determinant,1e-90),
            model+': determinant/singular-mass identity failed')
    return {'modes':modes,'projections':projections,'heavy_groups':group_rows,
            'leave_one_group_out':leave_out,'group_pair_overlaps':cross,
            'fermion_eigenspace_coherence':coh,
            'fermi_group_strength_ev':strength,
            'matrix_frobenius_ev':norm(total),
            'masses_ev':singular.tolist(),
            'determinant_absolute_ev3':determinant,
            'determinant_scaled_by_frobenius_cubed':determinant/max(norm(total)**3,ZERO),
            'matrix_reconstruction_relative_error':norm(full-total)/max(norm(full),norm(total),ZERO)}


def _write_json(path:Path,obj:dict)->None:
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def _write_csv(path:Path,rows:list[dict],fields:list[str])->None:
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)


def mechanism_explanation(row:dict)->str:
    m=row['lightest_mode'];tag=m['mechanism_diagnostic']
    if tag=='mode_geometry_dominant':
        return ('Its smallest singular direction has weak overlap with the heavy-group '
                'matrices; geometric suppression dominates the diagnostic split.')
    if tag=='heavy_eigenspace_cancellation_dominant':
        return ('Large projections from distinct heavy-mass eigenspaces partially cancel '
                'in the lightest singular direction.')
    if tag=='unresolved_zero_mode':
        return ('The lightest singular direction is below the mode-projection numerical threshold; '
                'the mechanism cannot be uniquely diagnosed from these data.')
    return ('Both projection geometry and inter-eigenspace interference contribute; '
            'neither is separated by more than a factor of two in this diagnostic.')


def report(payload:dict)->str:
    rows=payload['models'];get=lambda model:next(r for r in rows if r['model']==model)
    s=['# Chat A — Batch 4: why the 16 neutrino-mass hierarchies differ','',
       '**Source:** saved Batch-1/2/3 conditional UV results; no new Feynman diagrams, UV runs or EFT matching. '
       '**Physical source normalization remains 0/16 certified.**','',
       '## Analytical separation: flavour geometry vs heavy-state cancellation','',
       'The exact source-conditional heavy-state decomposition in Batch 3 is',
       r'\[M_\nu=\sum_g H_g,\qquad H_g=\frac{1}{16\pi^2}\sum_q\sum_{k\in g}'
       r'K_q(m_k^2)\,\mathrm{sym}\big[Z_{s_q}(M_F^\dagger P_k)Y_{t_q}^{T}\big].\]',
       'Here `g` labels a **distinct heavy singular-mass eigenspace**. Degenerate individual '
       'eigenvectors cannot be assigned a unique contribution. The complete projector is invariant.',
       'For the complex SVD `Mν=U diag(m1,m2,m3) V†`, take singular vectors `u_i`,`v_i` and set',
       r'\[a_{gi}=u_i^\dagger H_g v_i,\qquad m_i=\sum_g a_{gi}.\]',
       'Define `G = Σ_g ||H_g||F`, `A_i = Σ_g |a_gi| / G` (geometric projection fraction) '
       'and `C_i = m_i / Σ_g |a_gi|` (coherence fraction, 0 = cancellation, 1 = aligned). '
       'Then, exactly (up to roundoff),',
       r'\[\boxed{\frac{m_i}{G}=A_i C_i}.\]',
       'This disambiguates a small eigenmass caused by **weak flavour projection** '
       '(`A_i` small) from one caused by **phase cancellation between heavy eigenspaces** '
       '(`C_i` small). No notion of adding group Takagi masses is used.','',
       '## All-16 measured mechanism diagnostics','',
       '| Model | fermion | heavy mass groups | m1/m3 | A1 (geometry) | C1 (coherence) | lightest-mode diagnosis |',
       '|---|---|---:|---:|---:|---:|---|']
    for r in rows:
        m=r['lightest_mode']
        coh=('unresolved' if m['coherence_fraction'] is None else f'{m["coherence_fraction"]:.5g}')
        s.append(f'| {r["model"]} | {r["fermion_type"]} | {r["number_heavy_mass_groups"]} | '
                 f'{r["m1_over_m3"]:.5g} | {m["geometric_projection_fraction"]:.5g} | {coh} | '
                 f'{m["mechanism_diagnostic"]} |')
    for model in ('C_-1','C_+1','C_-3'):
        r=get(model);m=r['lightest_mode']
        s+=['',f'### {model}: {r["fermion_type"]}, m1/m3={r["m1_over_m3"]:.6g}',
            f'- Frozen ||y1||F × ||y2||F = {r["frozen_yukawa_norm_product"]:.7g}; '
            f'|lambdaT3| = {r["frozen_lambdaT3_abs"]:.5g}.',
            f'- The full lightest mass is {r["m1_ev"]:.7g} eV; the corresponding '
            f'group-projection sum of absolute magnitudes is {m["projected_strength_ev"]:.7g} eV.',
            f'- A1={m["geometric_projection_fraction"]:.6g}; '
            f'C1={m["coherence_fraction"] if m["coherence_fraction"] is not None else "unresolved"}; '
            +mechanism_explanation(r),
            '- Exact mode amplitudes and complex phases are listed in `batch04_mode_projections.csv`; '
            'deleting a group is a **diagnostic counterfactual**, not a new UV calculation.']
    s+=['','## Additional cross-family observations','',
        '- **C_-1 versus C_+1 is not a pure Yukawa-controlled comparison:** '
        'the frozen fitted Yukawas differ **and the fermion construction differs** '
        '(neutral Majorana C_-1 versus Dirac C_+1). The group decomposition identifies '
        'numerical mechanisms, but does not isolate their causal origin with all other inputs fixed.',
        '- **A/D matched cases:** compare their mode fractions and group interference, not just '
        'their close `m3`; near equality of selected frozen matrices is not an all-parameter SU(2) identity.',
        '- **E:** Batch 3 found two perfectly aligned scalar-charge channels. The heavy-eigenspace '
        'coherence here is a *different* quantity; aligned charge channels do not preclude '
        'cancellation or weak projections among heavy generations.',
        '- The determinant relation `|det Mν|=m1 m2 m3` provides a complementary '
        'flavour-invariant hierarchy measure; a small determinant need not mean a small loop amplitude.','',
        '## Source and interpretation limits','',
        '- Calculations are on the complex **source-conditional** matrix in eV, not physical-source-certified amplitudes.',
        '- The mode diagnostic threshold compares the two exact suppression factors: '
        '“dominant” means one is at least twice smaller; no model-selection significance is implied.',
        '- A near-degenerate *light-neutrino* singular pair makes individual singular-vector '
        'projections orientation-sensitive; the summed singular masses and group matrices remain unambiguous.',
        '- The direct B1/B2 spinor, symmetry, portal and CG physical normalization is still not certified.',
        '', '## Per-model measured interpretation','']
    for r in rows:
        m=r['lightest_mode']
        c='undefined' if m['coherence_fraction'] is None else f'{m["coherence_fraction"]:.6g}'
        s.append(f'- **{r["model"]}**: `m1/m3={r["m1_over_m3"]:.6g}`, '
                 f'`A1={m["geometric_projection_fraction"]:.6g}`, `C1={c}`. '
                 +mechanism_explanation(r))
    s+=['','The exact `m1/G=A1*C1` identity is the analytic explanation of *how* a '
        'particular frozen texture suppresses the lightest mode. It does not, by itself, '
        'prove which lagrangian coupling made that texture arise.','']
    return '\n'.join(s)+'\n'


def analyse(root:Path, *, study_dir:Path|None=None, batch02_dir:Path|None=None,
            batch03_dir:Path|None=None, output_dir:Path|None=None)->dict:
    root=root.resolve()
    owner=(root/BASE).resolve()
    paths={'study':(study_dir or root/STUDY).resolve(),
           'batch02':(batch02_dir or root/BATCH02).resolve(),
           'batch03':(batch03_dir or root/BATCH03).resolve(),
           'output':(output_dir or root/TARGET).resolve()}
    for key,p in paths.items():
        require(p.is_relative_to(owner) and p!=owner,key+': path outside Chat A output tree')
    require(len(set(paths.values()))==4,'Source and target paths must be distinct')
    target=paths['output']
    require(not target.exists() or (target.is_dir() and not any(target.iterdir())),
            'Output directory already populated: '+str(target))
    summary=load(paths['study']/'all16_summary.json')
    batch2=load(paths['batch02']/'batch02_analysis.json')
    batch3=load(paths['batch03']/'batch03_analytic_explanation.json')
    manifest=load(paths['batch03']/'SOURCE_MANIFEST.json')
    require(summary.get('format')=='t3-conditional-pattern-study-all16-v1',
            'Batch-1 schema changed')
    require(batch2.get('format')=='t3-chatA-readonly-batch02-numerical-patterns-v1',
            'Batch-2 schema changed')
    require(batch3.get('format')=='t3-chatA-batch03-analytic-scalar-kernel-v1',
            'Batch-3 schema changed')
    for key,record in (('Batch-1',summary),('Batch-2',batch2),('Batch-3',batch3)):
        require(record.get('physical_normalisation_certified') is False,
                key+': unexpected physical certification status')
    require(summary['counts']['conditional_evaluated']==16 and summary['counts']['blocked']==0 and
            summary['counts']['physical_source_certified']==0,
            'Batch-1 incomplete or certified unexpectedly')
    require(batch2['counts']['models']==16 and batch2['counts']['source_scan_rows']==288 and
            batch2['counts']['physical_source_certified']==0,
            'Batch-2 incomplete or certified unexpectedly')
    require(batch3['model_count']==16,'Batch-3 incomplete')
    require(batch3.get('safe_for_physical_source_certification') is False,
            'Batch-3 claims physical certification')
    sha1=digest(paths['study']/'all16_summary.json')
    sha2=digest(paths['batch02']/'batch02_analysis.json')
    sha3=digest(paths['batch03']/'batch03_analytic_explanation.json')
    require(batch2.get('input_run_sha256')==sha1 and batch3.get('source_summary_sha256')==sha1 and
            manifest.get('source_summary_sha256')==sha1 and
            manifest.get('batch02_analysis_sha256')==sha2,
            'Batch-1/2/3 provenance SHA256 mismatch')
    require(manifest.get('physical_source_certified') is False,
            'Batch-3 source manifest unexpectedly claims physical certification')
    catalogue=model_dictionary()
    names=[m['model'] for m in catalogue]
    require(len(names)==16 and len(set(names))==16,'Live source model dictionary changed')
    for key,rows in (('Batch-1',summary['models']),('Batch-2',batch2['models']),
                     ('Batch-3',batch3['models'])):
        require([r['model'] for r in rows]==names,key+': all-16 model ordering/coverage changed')
    originals={r['model']:r for r in summary['models']}
    second={r['model']:r for r in batch2['models']}
    third={r['model']:r for r in batch3['models']}
    outputs=[];projrows=[];grouprows=[];leaverows=[];pairrows=[]
    filemanifest={'source_run_sha256':sha1,'source_batch02_sha256':sha2,
                  'source_batch03_sha256':sha3,'original_source_guard_executed':True,
                  'source_revision':summary['provenance']['revision'],
                  'physical_normalisation_certified':False,'models':{}}
    for label in names:
        saved_file=paths['study']/'models'/f'{label}.json'
        last_file=paths['batch03']/'models'/f'{label}.json'
        saved=load(saved_file);last=load(last_file)
        row=third[label];base=originals[label];prev=second[label]
        require(saved.get('model')==label and saved.get('status')=='conditional_uv_evaluated' and
                saved.get('physical_normalisation_certified') is False and
                row.get('physical_normalisation_certified') is False,
                label+': conditional model provenance changed')
        orig_sha=saved['source']
        require(digest(saved_file)==manifest['models'][label]['model_detail_sha256'] and
                orig_sha['benchmark_sha256']==base['benchmark_sha256'] and
                orig_sha['input_sha256']==base['component_input_sha256'] and
                orig_sha['benchmark_sha256']==manifest['models'][label]['benchmark_sha256'] and
                orig_sha['input_sha256']==manifest['models'][label]['weyl_input_sha256'] and
                orig_sha['benchmark_sha256']==prev['benchmark_sha256'] and
                orig_sha['input_sha256']==prev['component_input_sha256'],
                label+': model/benchmark/Weyl SHA256 provenance conflict')
        _source_hash_guard(root,orig_sha,summary,label)
        require(last.get('model')==label and last.get('physical_normalisation_certified') is False and
                last.get('source_SHA256')==orig_sha,
                label+': Batch-3 detail source mismatch')
        require(last.get('summary')==row,label+': Batch-3 per-model summary drift')
        total=matrix(saved['mnu_ev'],label+' original mnu eV')
        last_m=matrix(last['observed_source_conditional_mnu_gev'],label+' Batch3 observed')*1e9
        exact_m=matrix(last['exact_closed_kernel_mnu_gev'],label+' Batch3 exact')*1e9
        require(near(total,last_m,1e-8) and near(total,exact_m,1e-7),
                label+': Batch-3 matrix reconstruction disagrees with source')
        groups=[]
        for idx,g in enumerate(row['heavy_mass_groups']):
            val=float(g['fermion_mass_gev']);mult=int(g['eigenspace_multiplicity'])
            require(math.isfinite(val) and val>0 and mult>=1,label+': invalid heavy group metadata')
            mm=matrix(g['matrix_ev'],label+f' heavy group {idx}')
            require(abs(norm(mm)-float(g['matrix_norm_ev']))<1e-7*max(norm(mm),1e-25),
                    label+': heavy-group source norm drift')
            groups.append({'mass_gev':val,'multiplicity':mult,'matrix_ev':mm})
        require(1<=len(groups)<=24,label+': suspicious heavy-eigenspace group count')
        if groups:
            masses=np.array([p['mass_gev'] for p in groups],dtype=float)
            require(np.all(np.diff(masses)>0),label+': heavy groups not strictly ordered')
        # An independent spectral check of the saved heavy group multiplicities.
        inpath=Path(orig_sha['input'])
        if not inpath.is_absolute():inpath=root/inpath
        inp=load(inpath)
        ff=inp['fermion_mass_gev']
        f=np.asarray(ff['real'],dtype=float)+1j*np.asarray(ff['imag'],dtype=float)
        require(f.ndim==2 and f.shape[0]==f.shape[1] and np.all(np.isfinite(f)),
                label+': invalid source heavy Weyl mass matrix')
        from_groups=np.sort(np.concatenate([np.full(p['multiplicity'],p['mass_gev']) for p in groups]))
        from_mass=np.sort(np.linalg.svd(f,compute_uv=False))
        require(from_groups.shape==from_mass.shape and np.allclose(from_groups,from_mass,rtol=2e-7,atol=0),
                label+': heavy eigenspace masses/multiplicity do not match source M_F')
        info=geometry(total,groups,label)
        ms=np.asarray(saved['takagi_masses_ev_sorted'],dtype=float)
        require(ms.shape==(3,) and np.allclose(info['masses_ev'],ms,rtol=3e-7,atol=1e-10*ms[-1]),
                label+': light-neutrino masses changed')
        require(abs(info['fermion_eigenspace_coherence']-row['heavy_mass_group_coherence'])<2e-7,
                label+': Batch-3 heavy-state coherence drift')
        light=info['modes'][0]
        model_row={'model':label,'family':label[0],'alpha':int(label[2:]),
                   'fermion_type':row['fermion_type'],'m1_ev':float(ms[0]),
                   'm2_ev':float(ms[1]),'m3_ev':float(ms[2]),
                   'm1_over_m3':float(ms[0]/ms[2]),'m2_over_m3':float(ms[1]/ms[2]),
                   'number_heavy_mass_groups':len(groups),
                   'heavy_group_coherence':info['fermion_eigenspace_coherence'],
                   'lightest_mode':light,'modes':info['modes'],
                   'determinant_scaled_by_frobenius_cubed':info['determinant_scaled_by_frobenius_cubed'],
                   'matrix_reconstruction_relative_error':info['matrix_reconstruction_relative_error'],
                   'frozen_yukawa_norm_product':float(row['yukawa_frobenius_product']),
                   'frozen_lambdaT3_abs':float(row['lambdaT3_abs']),
                   'source_conditional_only':True,'physical_normalisation_certified':False}
        outputs.append(model_row)
        projrows.extend(info['projections']);grouprows.extend(info['heavy_groups']);
        leaverows.extend(info['leave_one_group_out']);pairrows.extend(info['group_pair_overlaps'])
        filemanifest['models'][label]={'batch1_model_sha256':digest(saved_file),
                                      'batch03_model_sha256':digest(last_file),
                                      'benchmark_sha256':orig_sha['benchmark_sha256'],
                                      'candidate_weyl_input_sha256':orig_sha['input_sha256']}
    require(len(outputs)==16,'Incomplete frozen model coverage')
    payload={'format':'t3-chatA-batch04-flavour-geometry-v1',
             'basis':'conditional_source_Weyl',
             'model_count':16,'models':outputs,
             'source_sha256':filemanifest,
             'physical_normalisation_certified':False,
             'methods':{'mass_identity':'m_i=sum_g u_i^dag H_g v_i',
                        'factorization':'m_i/G=(sum_g |a_gi|/G)*(m_i/sum_g |a_gi|)',
                        'group_definition':'maximal source heavy singular-mass eigenspaces, not individual degenerate eigenvectors',
                        'mode_diagnosis_rule':'one suppression factor must be at least twice smaller than other',
                        'leave_one_out':'matrix-only counterfactual, not re-fitted benchmark'}}
    # Validate all 16 before the first file write.
    text=report(payload)
    handoff={'format':'t3-chatA-batch04-handoff-v1','models':16,
             'physical_source_certified':False,'safe_for_conditional_interpretation':True,
             'safe_for_physical_source_certification':False,
             'source_manifest_sha256':filemanifest,
             'key_identity':payload['methods']['factorization'],
             'primary_comparison':['C_-1','C_+1','C_-3'],
             'unresolved':['B1 Matchete-Weyl physical spinor dictionary',
                           'B2 gauge/portal/CG physical normalization']}
    target.mkdir(parents=True,exist_ok=True)
    _write_json(target/'batch04_flavour_analysis.json',payload)
    _write_json(target/'SOURCE_MANIFEST.json',filemanifest)
    _write_json(target/'HANDOFF.json',handoff)
    (target/'BATCH04_FLAVOUR_REPORT.md').write_text(text,encoding='utf-8')
    _write_csv(target/'batch04_model_mechanisms.csv',outputs,
               ['model','family','alpha','fermion_type','m1_ev','m2_ev','m3_ev',
                'm1_over_m3','m2_over_m3','number_heavy_mass_groups',
                'frozen_yukawa_norm_product','frozen_lambdaT3_abs','heavy_group_coherence','determinant_scaled_by_frobenius_cubed',
                'matrix_reconstruction_relative_error'])
    _write_csv(target/'batch04_mode_projections.csv',projrows,
               ['model','mode','heavy_group_index','mass_gev','amplitude_real_ev',
                'amplitude_imag_ev','amplitude_absolute_ev','amplitude_fraction_of_projected_strength'])
    _write_csv(target/'batch04_heavy_groups.csv',grouprows,
               ['model','heavy_group_index','fermion_mass_gev','eigenspace_multiplicity',
                'matrix_frobenius_ev','frobenius_fraction_of_group_strength',
                'relative_numerical_rank','singular_masses_ev'])
    _write_csv(target/'batch04_leave_one_group_out.csv',leaverows,
               ['model','omitted_heavy_group_index','fermion_mass_gev','m1_without_ev',
                'm2_without_ev','m3_without_ev','m1_change_ev','m2_change_ev',
                'm3_change_ev','matrix_fraction_removed'])
    _write_csv(target/'batch04_heavy_group_overlaps.csv',pairrows,
               ['model','group_i','group_j','complex_overlap_ev2','frobenius_overlap_cosine'])
    return payload
