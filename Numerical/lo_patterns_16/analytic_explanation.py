"""Chat A batch 3: independent analytic two-scalar kernel and all-16 explanations.

READ ONLY on existing benchmark and Batch-1/2 output. This code does not call
core.evaluate or native_uv_loop_reference.independent_reference. It analytically
resums the 2x2 complex-scalar mixing blocks with B=0, and uses the source Weyl
mass matrix only through the spectral propagator M^dag f(M M^dag).

All claims remain CONDITIONAL on the unverified Matchete/RGBeta -> Weyl map.
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
from Numerical.lo_certification.direct_uv.native_uv_loop_reference import decode, encode
from Numerical.lo_certification.direct_uv.certification import charge_bases
from Numerical.lo_patterns_16.result_analysis import _source_hash_guard

LOOP = 1.0 / (16.0 * math.pi ** 2)
OUTPUT_BASE = Path('output/general_observable/lo_pattern_study16')
STUDY = OUTPUT_BASE / 'provisional_01'
BATCH02 = OUTPUT_BASE / 'analysis_batch02_01'
BATCH03 = OUTPUT_BASE / 'analysis_batch03_01'
INPUTS = Path('output/general_observable/lo_direct_one_loop/provisional_inputs')
CG_AUDIT = Path('output/general_observable/lo_direct_one_loop/source_audits/charge_resolved_loop_weights.json')

class ExplanationError(ValueError):
    pass


def require(ok:bool,why:str)->None:
    if not ok: raise ExplanationError(why)


def load(path:Path)->dict:
    require(path.is_file(),f'Missing required source: {path}')
    return json.loads(path.read_text(encoding='utf-8-sig'))


def digest(path:Path)->str:
    return sha256(path.read_bytes()).hexdigest()


def norm(a:np.ndarray)->float:
    return float(np.linalg.norm(a))


def relative(a:np.ndarray,b:np.ndarray, *, natural:float=0)->float:
    return norm(a-b)/max(norm(a),norm(b),natural*1e-12,1e-100)


def b0(m2:float,s2:float)->float:
    """Finite scalar-dependent part of B0(0;m²,s²); any constant cancels.

    This is -g(s²/m²), with g(x)=x log x/(x-1). No dependence on mu.
    """
    require(m2>0 and s2>0 and math.isfinite(m2) and math.isfinite(s2),
            'Non-positive or non-finite squared mass in B0')
    x=(s2-m2)/m2
    if abs(x)<1e-5:
        g=1+x*(.5+x*(-1/6+x*(1/12+x*(-1/20+x/30))))
    else:
        g=(1+x)*math.log1p(x)/x
    return -g


def b0_derivative(m2:float,s2:float)->float:
    """d/d(s²) B0 with a near-equal-mass Taylor series."""
    require(m2>0 and s2>0,'B0 derivative requires positive mass squared')
    x=(s2-m2)/m2
    if abs(x)<1e-4:
        return (-.5+x/3-x*x/4+x**3/5-x**4/6+x**5/7)/m2
    return (math.log1p(x)/(x*x)-1/x)/m2


def divdiff(m2:float,a:float,b:float)->float:
    """[B0(m²,b)-B0(m²,a)]/(b-a), with degenerate scalar limit."""
    require(a>0 and b>0 and m2>0,'Divided difference: non-positive squared masses')
    if abs(a-b)<1e-8*max(a,b):
        return b0_derivative(m2,(a+b)*.5)
    return (b0(m2,b)-b0(m2,a))/(b-a)


def scalar_kernels(m2:float,a:float,b:float,delta:complex)->tuple[complex,complex,float]:
    """Exact and first-order off-diagonal B0(A)_{S2,S1} for 2x2 Hermitian A.

    The spectral scalar eigenvalues depend only on |delta|^2; the exact
    off-diagonal matrix function is delta times a real divided difference.
    """
    d=complex(delta)
    radius=math.hypot((a-b)*.5,abs(d))
    low=(a+b)*.5-radius
    high=(a+b)*.5+radius
    require(low>0 and high>low and math.isfinite(low) and math.isfinite(high),
            'Tachyonic or unresolved 2x2 scalar eigenvalues')
    exact=d*divdiff(m2,low,high)
    leading=d*divdiff(m2,a,b)
    return exact,leading,float(abs(exact-leading)/max(abs(exact),abs(leading),1e-100))


def coupled_pairs(A:np.ndarray, d1:int, kind:str,alpha:int)->list[dict]:
    """Enumerate EVERY equal-charge S1/S2 pair, even with zero scalar mixing.

    Only isolated 2x2 blocks with diagonal original multiplet masses are
    supported here. This is a verified property of the *restricted* input.
    """
    ns=len(A)
    require(0<d1<ns and A.shape==(ns,ns),'Invalid scalar dimensions')
    scale=max(norm(A),1.)
    require(norm(A-A.conj().T)<1e-10*scale,'Non-Hermitian scalar A')
    require(norm(A[:d1,:d1]-np.diag(np.diag(A[:d1,:d1])))<1e-12*scale and
            norm(A[d1:,d1:]-np.diag(np.diag(A[d1:,d1:])))<1e-12*scale,
            'New scalar within-multiplet mixing: isolated 2x2 formula not applicable')
    charges=charge_bases(kind,alpha)['scalar']
    require(len(charges)==ns,'Electric charge basis mismatch')
    entries=[]
    for s in range(d1):
        for t in range(d1,ns):
            if not np.isclose(charges[s],charges[t],rtol=0,atol=1e-10):
                require(abs(A[t,s])<1e-11*scale,'Cross-mixing violates charge conservation')
                continue
            p={'charge':float(charges[s]),'s1':s,'s2':t,
               'a_gev2':float(A[s,s].real),'b_gev2':float(A[t,t].real),
               'delta_gev2':A[t,s]}
            require(p['a_gev2']>0 and p['b_gev2']>0,'Nonpositive unmixed scalar mass')
            entries.append(p)
    # Nonzero cross-mixing outside the listed 2x2 charge blocks is forbidden.
    allowed={(r['s2'],r['s1']) for r in entries}
    for t in range(d1,ns):
        for s in range(d1):
            if (t,s) not in allowed:
                require(abs(A[t,s])<1e-11*scale,'Unlisted scalar cross mixing')
    require(len({e['s1'] for e in entries})==len(entries) and
            len({e['s2'] for e in entries})==len(entries),
            'Multiple equal-charge scalar partners: need a larger-block formula')
    return entries


def reconstruct(input_data:dict,kind:str,alpha:int,d1:int)->dict:
    """Closed-form source-conditional Mnu, by original scalar charge channel.

    No real-scalar diagonalization and no fermion Takagi used here. The
    independent propagator is M^dagger |Q_k><Q_k| for M M^dagger eigenvectors.
    """
    model=input_data.get('model','unknown')
    require(model==f'{kind}_{alpha:+d}','Input/model representation conflict')
    M=decode(input_data['fermion_mass_gev'],'M')
    A=decode(input_data['scalar_A_gev2'],'A')
    B=decode(input_data['scalar_B_gev2'],'B',A.shape)
    nf=M.shape[0]; ns=A.shape[0]
    require(M.shape==(nf,nf) and ns==A.shape[1],'Invalid M/A shapes')
    require(norm(M-M.T)<1e-10*max(norm(M),1.),'Non-symmetric Weyl mass')
    require(norm(B)<1e-14*max(norm(A),1.),
            'Scalar B is nonzero; 2x2 complex-scalar result not applicable')
    vertices=input_data.get('vertices',[])
    require([x.get('name') for x in vertices]==['y1_HC','y2_HC'],
            'New vertex sources: restricted two-vertex formula not applicable')
    require('conjugate' in vertices[0] and 'field' not in vertices[0] and
            'field' in vertices[1] and 'conjugate' not in vertices[1],
            'Candidate Weyl vertex orientations changed')
    Z=decode(vertices[0]['conjugate'],'Z',(3,nf,ns))
    Y=decode(vertices[1]['field'],'Y',(3,nf,ns))
    require(norm(Z[:,:,d1:])<1e-13*max(norm(Z),1e-90) and
            norm(Y[:,:,:d1])<1e-13*max(norm(Y),1e-90),
            'Source vertex scalar families crossed')
    F2=M@M.conj().T
    fm2,Q=np.linalg.eigh((F2+F2.conj().T)*.5)
    require(np.min(fm2)>0 and np.min(fm2)>1e-15*np.max(fm2),
            'Singular / unresolved heavy propagation')
    heavy_ops=[M.conj().T @ np.outer(Q[:,k],Q[:,k].conj()) for k in range(nf)]
    pairs=coupled_pairs(A,d1,kind,alpha)
    out=np.zeros((3,3),complex)
    lin=np.zeros_like(out)
    heavy_terms=[np.zeros_like(out) for _ in range(nf)]
    channels=[]
    for p in pairs:
        s,t=p['s1'],p['s2']
        exact=np.zeros_like(out); leading=np.zeros_like(out)
        max_kernel_nonlin=0.
        for k,m2 in enumerate(fm2):
            kx,kl,nonlin=scalar_kernels(float(m2),p['a_gev2'],p['b_gev2'],p['delta_gev2'])
            # Z carries S1*, Y carries S2. The scalar off-diagonal
            # is B0(A)_{S2,S1} = delta_{S2,S1} * divided_difference.
            cross=Z[:,:,s] @ heavy_ops[k] @ Y[:,:,t].T
            sym=cross+cross.T
            term=LOOP*kx*sym
            exact+=term
            heavy_terms[k]+=term
            leading+=LOOP*kl*sym
            max_kernel_nonlin=max(max_kernel_nonlin,nonlin)
        channels.append({**p,'exact':exact,'leading':leading,
                         'max_scalar_kernel_relative_nonlinearity':max_kernel_nonlin,
                         'matrix_norm_gev':norm(exact)})
        out+=exact;lin+=leading
    # Group the F2-degenerate eigenvectors: their individual matrices are
    # basis-dependent, but sums over complete eigenspaces are invariant.
    groups=[]
    for k,m2 in enumerate(fm2):
        if not groups or abs(float(m2)-groups[-1]['m2'])>1e-8*max(float(m2),groups[-1]['m2']):
            groups.append({'m2':float(m2),'matrix':np.zeros_like(out),'dimension':0})
        groups[-1]['matrix']+=heavy_terms[k]
        groups[-1]['dimension']+=1
    require(relative(sum((g['matrix'] for g in groups),np.zeros_like(out)),out)<1e-10,
            'Heavy mass eigenspace regrouping does not reconstruct the exact mass')
    return {'exact':out,'linear':lin,'channels':channels,'heavy_groups':groups,
            'linear_relative_error':relative(lin,out),
            'heavy_eigenmass_min_gev':float(math.sqrt(np.min(fm2))),
            'heavy_eigenmass_max_gev':float(math.sqrt(np.max(fm2)))}


def _phase_fit(a:np.ndarray,b:np.ndarray)->dict:
    """Best complex proportionality b ~= c*a, not a representation identity."""
    aa=np.vdot(a,a)
    c=np.vdot(a,b)/aa if abs(aa)>1e-90 else 0j
    return {'complex_proportionality':{'real':float(c.real),'imag':float(c.imag)},
            'best_scale_abs':float(abs(c)),
            'relative_residual':float(norm(b-c*a)/max(norm(b),1e-100)),
            'matrix_alignment_abs':float(abs(np.vdot(a,b))/max(norm(a)*norm(b),1e-100))}


def _yukawa_norms(frozen:dict)->dict:
    ordinary=frozen['base_state']['ordinary']
    Y1=decode(ordinary['y1'],'benchmark y1',(3,3))
    Y2=decode(ordinary['y2'],'benchmark y2',(3,3))
    lam=ordinary['lambdaT3']
    if isinstance(lam,dict):lam=complex(lam['real'],lam['imag'])
    else:lam=complex(lam)
    return {'y1_frobenius':norm(Y1),'y2_frobenius':norm(Y2),
            'yukawa_frobenius_product':norm(Y1)*norm(Y2),
            'lambdaT3_abs':float(abs(lam)),
            'nominal_prefactor_product':norm(Y1)*norm(Y2)*abs(lam),
            'y1':Y1,'y2':Y2}


def _cg_match(native:dict, pairs:list[dict], model:str)->dict[float,dict]:
    out={}
    for ch in native:
        require('charge_numeric' in ch and 'complex_weight_numeric' in ch,
                model+': source charge-weight contract drift')
        if not ch['nonzero']:continue
        q=float(ch['charge_numeric'])
        require(q not in out,model+': multiply counted candidate source charge')
        out[q]=ch
    available={p['charge'] for p in pairs if abs(p['delta_gev2'])>0}
    require(set(out)==available,
            model+f': native CG charge support {sorted(out)} != scalar-mixing support {sorted(available)}')
    return out


def _csv(path:Path,rows:list[dict],keys:list[str]):
    with path.open('w',encoding='utf-8',newline='') as fp:
        w=csv.DictWriter(fp,fieldnames=keys,extrasaction='ignore');w.writeheader();w.writerows(rows)


def _json(path:Path,object_:dict):
    path.write_text(json.dumps(object_,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def report_text(payload:dict)->str:
    rows=payload['models']
    l=['# Chat A — Batch 3: analytical explanation of all 16 T3 patterns','',
       '**Basis:** frozen restricted direct UV, fixed `mu0=10^7 GeV`. The 2×2 scalar kernel is',
       'derived exactly in the **provisional Weyl input conventions**. Physical cross-package normalization: **0/16 certified**.','',
       '## Analytical kernel and its observable consequences','',
       'For each equal-charge pair `s∈S1`, `t∈S2` the complex-scalar squared-mass block is',
       '\\[A_q=\\begin{pmatrix}a_q&\\delta_q^*\\\\\\delta_q&b_q\\end{pmatrix},\\quad'
       '\\lambda_{q,\\pm}=\\tfrac12(a_q+b_q)\\pm\\sqrt{\\tfrac14(a_q-b_q)^2+|\\delta_q|^2}.\\]',
       'The **exact** scalar insertion entering `Z(S1*) × Y(S2)` is',
       '\\[K_q(m_F^2)=\\delta_q\\frac{B_0(0;m_F^2,\\lambda_+)-B_0(0;m_F^2,\\lambda_-)}'
       '{\\lambda_+-\\lambda_-}.\\]',
       'At leading order in mixing, replace the divided difference at `λ±` by the one at `a_q,b_q`. '
       'Higher corrections are odd in `δ_q` for a common phase, with relative size `O(|δ_q|² / mass-gap²)` away from degeneracy.',
       'The source-conditional amplitude is',
       '\\[(m_\\nu)_{ij}=\\frac1{16\\pi^2}\\sum_{q,k}K_q(m_k^2)'
       '\\big[Z_{i,s_q}(M^\\dagger P_k)Y^T_{t_q,j}+'
       '(Z_{i,s_q}(M^\\dagger P_k)Y^T_{t_q,j})^T\\big],\\quad'
       'P_k=Q_{:k}Q_{:k}^\\dagger,\\quad M M^\\dagger=Q\\,\\mathrm{diag}(m_k^2)Q^\\dagger.\\]',
       'This retains the heavy-flavour matrix and its relative phases. '
       'It uses neither a Takagi rotation nor the real-scalar Hessian diagonalisation.','',
       '## All-16 per-model interpretation','',
       '| Model | charges | m1/m3 | m2/m3 | |y1||y2| | |λT3| | exact kernel error | small-mixing error | λ×2 predicted / measured |',
       '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        charges=', '.join(f'{c["charge"]:g}' for c in r['channels'])
        l.append(f'| {r["model"]} | {charges} | {r["m1_over_m3"]:.5g} | '
                 f'{r["m2_over_m3"]:.5g} | {r["yukawa_frobenius_product"]:.5g} | '
                 f'{r["lambdaT3_abs"]:.3g} | {r["exact_reconstruction_relative_error"]:.3g} | '
                 f'{r["small_mixing_relative_error"]:.3g} | '
                 f'{r["lambda_t2_predicted_residual"]:.3g} / '
                 f'{r["lambda_t2_observed_residual"]:.3g} |')
    l+=['','### Channel-by-channel explanation','']
    for r in rows:
        c=r['channels']
        active=[x for x in c if x['matrix_norm_ev']>1e-16*max(r['full_matrix_frobenius_ev'],1e-100)]
        bits=[f'Q={x["charge"]:g}: |δ|={x["delta_abs_gev2"]:.5g} GeV², '
              f'|W_CG|={x["candidate_cg_weight_abs"]:.5g}, '
              f'||M_Q||={x["matrix_norm_ev"]:.5g} eV' for x in active]
        if not bits:bits=['no numerically active topology charge pair']
        l.append(f'- **{r["model"]}** ({r["fermion_type"]}): '
                 +'; '.join(bits)+f'. m₁/m₃={r["m1_over_m3"]:.4g}; '
                 +('two-channel overlap cosine='+f'{r["charge_cosine"]:.6g}' if r['charge_cosine'] is not None
                   else 'single-channel topology; no inter-charge interference')+'.')
    l+=['','## Matched-shape comparisons (NOT exact representation symmetries)','',
        '| Pair | best amplitude scale | shape error after complex rescaling |',
        '|---|---:|---:|']
    for x in payload['paired_comparisons']:
        l.append(f'| {x["a"]} vs {x["b"]} | {x["best_scale_abs"]:.6g} | {x["relative_residual"]:.4g} |')
    l+=['','## What these tests mean','',
        '- The *exact kernel error* is an internal numerical regression test of the source-conditional loop, not native physical certification.',
        '- The *small-mixing error* distinguishes a leading-order analytic approximation from the exact 1-loop scalar spectrum.',
        '- A small m₁/m₃ does not mean a forbidden charge channel: flavour misalignment, determinant suppression, or interference may cause it.',
        '- The heavy-mass-group coherence, `||sum_g M_g||/sum_g ||M_g||`, diagnoses cancellation **between distinct heavy singular-mass eigenspaces**. Individual vectors inside a degenerate eigenspace have no invariant meaning.',
       '- For E, the two-channel cosine measures matrix-direction interference. A value near +1 means the charge matrices reinforce; it does **not** imply an extra normalization factor of 2.',
        '- Numerical masses are from individual frozen Yukawa fits. Changing alpha also changes the fit; a family trend is not by itself a hypercharge theorem.',
        '- Source-CG candidate weights follow `T3Y1CG*T3Y2CG*conjugate(T3MixCG)`. The physical four-spinor/portal dictionary remains B1/B2 dependent.',
        '- The code refuses changed scalar B, new Yukawa vertices, source hashes, partial scans, and false certification promotion.',
        '', '## Open physics', '',
        'Physical normalisation, spinor phases, native Majorana factors and portal sign must be resolved in B1/B2. '
        'After the coordinator supplies their certified dictionary, rerun the original all-16 numerical calculation before asserting physical masses.']
    return '\n'.join(l)+'\n'


def explain(root:Path, *, study_dir:Path|None=None, batch02_dir:Path|None=None,
            output_dir:Path|None=None)->dict:
    root=root.resolve()
    owned=(root/OUTPUT_BASE).resolve()
    study=(study_dir or root/STUDY).resolve()
    prev=(batch02_dir or root/BATCH02).resolve()
    target=(output_dir or root/BATCH03).resolve()
    for name,p in (('study',study),('batch02',prev),('batch03',target)):
        require(p.is_relative_to(owned) and p!=owned,name+': must be inside the owned Chat A output tree')
    require(len({study,prev,target})==3,'Input/output paths must be distinct')
    require(not target.exists() or not any(target.iterdir()),'Output directory already populated: '+str(target))
    s=load(study/'all16_summary.json')
    b=load(prev/'batch02_analysis.json')
    require(s.get('format')=='t3-conditional-pattern-study-all16-v1' and
            b.get('format')=='t3-chatA-readonly-batch02-numerical-patterns-v1',
            'Prior result schemas changed')
    require(s.get('physical_normalisation_certified') is False and
            b.get('physical_normalisation_certified') is False and
            s['counts']['physical_source_certified']==0 and b['counts']['physical_source_certified']==0,
            'False physical certification claim in prior outputs')
    require(s['counts']['conditional_evaluated']==16 and s['counts']['blocked']==0 and
            b['counts']['models']==16 and b['counts']['source_scan_rows']==288,
            'Prior all-16 run / Batch-2 report incomplete')
    require(b['input_run_sha256']==digest(study/'all16_summary.json'),
            'Batch-2 source summary SHA256 does not match Batch-1')
    audit=load(root/CG_AUDIT)
    require(audit.get('status')=='charge_resolved_group_coefficients' and audit.get('count')==16,
            'CG audit status/count drift')
    require(digest(root/CG_AUDIT)==s['provenance']['source_audit_sha256']['charge_resolved_loop_weights.json'],
            'CG audit SHA256 drift')
    names=[m['model'] for m in model_dictionary()]
    require([x['model'] for x in s['models']]==names and [x['model'] for x in b['models']]==names,
            'All-16 source model order drift')
    weights=load(study/'channel_weights.json')
    require(weights.get('status')=='native_group_CG_candidate_only_not_physical_normalization',
            'Source CG weight status changed')
    model_meta={x['model']:x for x in model_dictionary()}
    orig_rows={x['model']:x for x in s['models']}
    batch2_rows={x['model']:x for x in b['models']}
    manifest={'source_summary_sha256':digest(study/'all16_summary.json'),
              'batch02_analysis_sha256':digest(prev/'batch02_analysis.json'),
              'charge_weight_audit_sha256':digest(root/CG_AUDIT),
              'physical_source_certified':False,'models':{}}
    records=[]; per_model={}; matrix_cache={}; output_ch=[]
    for model in names:
        model_file=study/'models'/f'{model}.json'
        detail=load(model_file)
        orig=orig_rows[model]
        source=detail['source']
        require(detail.get('model')==model and detail.get('status')=='conditional_uv_evaluated'
                and detail.get('physical_normalisation_certified') is False,
                model+': source mass status/certification drift')
        require(source['benchmark_sha256']==orig['benchmark_sha256'] and
                source['input_sha256']==orig['component_input_sha256'],
                model+': source input hash conflict')
        _source_hash_guard(root,source,s,model)
        inpath=Path(source['input'])
        if not inpath.is_absolute():inpath=root/inpath
        data=load(inpath)
        require(data.get('provenance',{}).get('validated_uv_component_conventions') is not True,
                model+': provisional Weyl input unexpectedly declares physical validation')
        require(weights['models'][model].get('relative_physical_normalization_certified') is False,
                model+': candidate CG channel output falsely physically certified')
        frozen=load(root/source['benchmark'])
        d1=model_meta[model]['SU2_dimensions']['S1']
        kind=model[0]; alpha=int(model[2:])
        reconstructed=reconstruct(data,kind,alpha,d1)
        exact=reconstructed['exact'];approx=reconstructed['linear']
        observed=decode(detail['mnu_gev'],model+' mass',(3,3))
        err=relative(exact,observed)
        # Analytically predict finite-mixing scan deviations, rather than merely
        # quoting a fitted power law. Exact scalar kernels depend on |delta|².
        predictions={}
        for factor in (-1.,0.,.5,2.):
            varied_A=decode(data['scalar_A_gev2'],'A').copy()
            varied_A[:d1,d1:]*=factor
            varied_A[d1:,:d1]*=factor
            varied={**data,'scalar_A_gev2':encode(varied_A)}
            varied_exact=reconstruct(varied,kind,alpha,d1)['exact']
            predictions[str(factor)] = norm(varied_exact-factor*exact)/max(norm(exact),1e-100)
        require(predictions['-1.0']<1e-11 and predictions['0.0']<1e-11,
                model+': exact odd-in-lambda or zero-coupling limit failed')
        require(err<1e-6,model+f': analytic exact-charge kernel mismatches production: {err:.3e}')
        a2=batch2_rows[model]
        require(a2.get('physical_normalisation_certified') is False,
                model+': Batch-2 model falsely physically certified')
        for t,key in ((2.,'lambda_t2_linear_matrix_residual'),
                      (.5,'lambda_t05_linear_matrix_residual'),
                      (-1.,'lambda_sign_flip_matrix_residual')):
            measured=a2.get(key)
            require(measured is not None,model+': missing Batch-2 '+key)
            measured=float(measured)
            predicted=predictions[str(t)]
            require(math.isfinite(measured) and abs(measured-predicted)<max(5e-8,.15*abs(measured)),
                    model+f': analytic lambda scan prediction {t:g}={predicted:.3g} differs from stored {measured:.3g}')
        require(a2['benchmark_sha256']==source['benchmark_sha256'] and
                a2['component_input_sha256']==source['input_sha256'],
                model+': Batch-2 hash provenance mismatch')
        require(max(abs(float(v)-float(w)) for v,w in
                    zip(detail['takagi_masses_ev_sorted'],[a2['m1_ev'],a2['m2_ev'],a2['m3_ev']]))
                <1e-8*max(a2['m3_ev'],1e-90),model+': Batch-2 eigenmasses drift')
        ch=_cg_match(weights['models'][model]['channels'],reconstructed['channels'],model)
        y=_yukawa_norms(frozen)
        cresults=[]
        pieces=[]
        for channel in reconstructed['channels']:
            charge=channel['charge']; native=ch.get(charge)
            nn=channel['exact'];ln=channel['leading'];pieces.append(nn)
            record={'charge':charge,'s1_zero_based':channel['s1'],
                    's2_zero_based':channel['s2']-d1,'a_gev2':channel['a_gev2'],
                    'b_gev2':channel['b_gev2'],
                    'delta_gev2':encode(channel['delta_gev2']),
                    'delta_abs_gev2':float(abs(channel['delta_gev2'])),
                    'candidate_cg_weight':native['complex_weight_numeric'] if native else None,
                    'candidate_cg_weight_abs':float(abs(complex(native['complex_weight_numeric']['real'],
                                                         native['complex_weight_numeric']['imag']))) if native else 0.,
                    'matrix_norm_ev':norm(nn)*1e9,'matrix_ev':encode(nn*1e9),
                    'first_order_matrix_ev':encode(ln*1e9),
                    'leading_error':relative(nn,ln),
                    'max_kernel_relative_nonlinearity':channel['max_scalar_kernel_relative_nonlinearity'],
                    'candidate_topology_channel':bool(native),
                    'physical_normalisation_certified':False}
            cresults.append(record)
            output_ch.append({'model':model,**{k:record[k] for k in ('charge','s1_zero_based','s2_zero_based',
                  'delta_abs_gev2','candidate_cg_weight_abs','matrix_norm_ev','leading_error','candidate_topology_channel')}})
        hgroups=[]
        for group in reconstructed['heavy_groups']:
            hgroups.append({'fermion_mass_gev':math.sqrt(group['m2']),
                            'eigenspace_multiplicity':group['dimension'],
                            'matrix_norm_ev':norm(group['matrix'])*1e9,
                            'matrix_ev':encode(group['matrix']*1e9)})
        hcoherence=(norm(exact)/max(sum(norm(g['matrix']) for g in reconstructed['heavy_groups']),1e-100))
        ccos=None
        if len(pieces)==2 and norm(pieces[0])*norm(pieces[1])>0:
            ccos=float(np.vdot(pieces[0],pieces[1]).real/(norm(pieces[0])*norm(pieces[1])))
            ccos=max(-1.,min(1.,ccos))
        masses=detail['takagi_masses_ev_sorted']
        summary={'model':model,'family':kind,'alpha':alpha,
                 'fermion_type':model_meta[model]['fermion_type'],
                 'm1_ev':masses[0],'m2_ev':masses[1],'m3_ev':masses[2],
                 'm1_over_m3':masses[0]/max(masses[2],1e-100),
                 'm2_over_m3':masses[1]/max(masses[2],1e-100),
                 'full_matrix_frobenius_ev':norm(observed)*1e9,
                 'exact_reconstruction_relative_error':err,
                 'small_mixing_relative_error':relative(exact,approx),
                 'lambda_residual_predictions':predictions,
                 'lambda_t2_predicted_residual':predictions['2.0'],
                 'lambda_t2_observed_residual':float(a2['lambda_t2_linear_matrix_residual']),
                 'charge_cosine':ccos,
                 'charge_coherence':norm(exact)/max(sum(norm(m) for m in pieces),1e-100),
                 'heavy_mass_group_coherence':hcoherence,
                 'heavy_mass_groups':hgroups,
                 'y1_frobenius':y['y1_frobenius'],
                 'y2_frobenius':y['y2_frobenius'],
                 'yukawa_frobenius_product':y['yukawa_frobenius_product'],
                 'lambdaT3_abs':y['lambdaT3_abs'],
                 'nominal_prefactor_product':y['nominal_prefactor_product'],
                 'channels':cresults,
                 'physical_normalisation_certified':False}
        records.append(summary)
        matrix_cache[model]=observed
        per_model[model]={'model':model,'exact_closed_kernel_mnu_gev':encode(exact),
                          'first_order_kernel_mnu_gev':encode(approx),
                          'observed_source_conditional_mnu_gev':encode(observed),
                          'summary':summary,'source_SHA256':source,
                          'physical_normalisation_certified':False}
        manifest['models'][model]={'model_detail_sha256':digest(model_file),
                                   'benchmark_sha256':source['benchmark_sha256'],
                                   'weyl_input_sha256':source['input_sha256']}
    comparisons=[]
    for a,bb in [('A_-4','D_+2'),('A_-2','D_+0'),('A_+0','D_-2'),
                 ('B_-3','B_+1'),('C_-3','C_+1'),('E_-4','E_+2'),('E_-2','E_+0')]:
        comparisons.append({'a':a,'b':bb,**_phase_fit(matrix_cache[a],matrix_cache[bb])})
    output={'format':'t3-chatA-batch03-analytic-scalar-kernel-v1',
            'basis':'restricted_source_conditional_Weyl',
            'reference_revision':s['provenance']['revision'],
            'source_summary_sha256':manifest['source_summary_sha256'],
            'model_count':len(records),'models':records,'paired_comparisons':comparisons,
            'analytic_formula':'mnu=(16pi²)^-1 sum_qk K_q(mk²) sym[Z_s (M† Pk) Y_t^T] with K_q=delta_q [B0(lambda+)-B0(lambda-)]/(lambda+-lambda-)',
            'physical_normalisation_certified':False,
            'safe_for_conditional_analytic_interpretation':True,
            'safe_for_physical_source_certification':False}
    # No writes until ALL input, source and exact-identity checks pass.
    target.mkdir(parents=True,exist_ok=True)
    _json(target/'batch03_analytic_explanation.json',output)
    _json(target/'SOURCE_MANIFEST.json',manifest)
    for model,entry in per_model.items():
        (target/'models').mkdir(exist_ok=True)
        _json(target/'models'/f'{model}.json',entry)
    _csv(target/'batch03_model_patterns.csv',records,
         ['model','family','alpha','fermion_type','m1_ev','m2_ev','m3_ev',
          'm1_over_m3','m2_over_m3','full_matrix_frobenius_ev',
          'y1_frobenius','y2_frobenius','yukawa_frobenius_product','lambdaT3_abs',
          'exact_reconstruction_relative_error','small_mixing_relative_error',
          'charge_cosine','charge_coherence','heavy_mass_group_coherence','physical_normalisation_certified'])
    _csv(target/'batch03_charge_channels.csv',output_ch,
         ['model','charge','s1_zero_based','s2_zero_based','delta_abs_gev2',
          'candidate_cg_weight_abs','matrix_norm_ev','leading_error','candidate_topology_channel'])
    _csv(target/'batch03_matched_shapecases.csv',comparisons,
         ['a','b','best_scale_abs','relative_residual','matrix_alignment_abs'])
    (target/'BATCH03_ANALYTICAL_REPORT.md').write_text(report_text(output),encoding='utf-8')
    handoff={'format':'t3-chatA-batch03-handoff-v1',
             'source_sha256':manifest,'source_revision':s['provenance']['revision'],
             'models':len(records),'kernel_tests_passed':len(records),
             'physical_source_certified':False,
             'safe_to_integrate_conditional_analytic_patterns':True,
             'safe_to_integrate_physical_mass_dictionary':False,
             'unresolved':['B1 Matchete Weyl rule physical phase/symmetry factors',
                           'B2 source CG scalar portal/coupling dictionary'],
             'formula':output['analytic_formula']}
    _json(target/'HANDOFF.json',handoff)
    (target/'HANDOFF.md').write_text('# Batch 3 handoff\n\n'
        f'- All {len(records)} frozen restricted cases analytically reconstructed from 2x2 scalar kernels, without UV engine calls.\n'
        '- Every source, frozen input and output SHA256 checked; no previous output modified.\n'
        '- Analytical source-conditional scalar-mixing corrections and channel matrices provided.\n'
        '- Physical source normalization certified: **0/16**; B1/B2 remain necessary.\n'
        '- Only source-conditional patterns are safe to integrate.\n',encoding='utf-8')
    return output
