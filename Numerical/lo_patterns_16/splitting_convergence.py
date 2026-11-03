"""Chat A Batch 8: exact finite-mixing scalar-splitting derivative and
transverse heavy-kernel factorization breaking; read only on all prior batches.

This is NOT an EFT calculation or a physical Matchete/RGBeta certification.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.lo_patterns_16 import kernel_factorization as k6
from Numerical.lo_patterns_16 import factorization_theorem as t7
from Numerical.lo_patterns_16 import texture_synthesis as t5
from Numerical.lo_patterns_16.analytic_explanation import (
    b0_derivative, coupled_pairs, divdiff, scalar_kernels,
)
from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary

BASE = k6.BASE
BATCH06 = k6.TARGET
BATCH07 = BASE / 'analysis_batch07_01'
TARGET = BASE / 'analysis_batch08_01'
REFERENCE = k6.REFERENCE_MODEL
TS = (0., 1./16., 1./8., 1./4., 1./2., 1.)
SMALL_TS = (1./16., 1./8., 1./4., 1./2.)
FLOOR = 1e-100


class SplittingConvergenceError(ValueError):
    pass


def require(ok: bool, message: str) -> None:
    if not ok:
        raise SplittingConvergenceError(message)


def read(path: Path) -> dict:
    require(path.is_file(), 'Missing source file: '+str(path))
    value=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(value,dict), 'Invalid source JSON object: '+str(path))
    return value


def norm(a:np.ndarray)->float:
    return float(np.linalg.norm(a))


def flatten(a:np.ndarray)->np.ndarray:
    require(a.shape==(3,3,3) and np.all(np.isfinite(a)), 'Invalid 3×3×3 heavy coefficient tensor')
    return a.reshape(27)


def exact_kernel_derivatives(m2:float,a:float,b:float,delta:complex)->tuple[complex,complex,complex]:
    """K, dK/da, dK/db for K=delta*f[lambda-,lambda+] with fixed delta.

    Unlike derivative_D in Batch 7, these derivatives include *all* odd
    scalar-mixing powers. Valid when the 2×2 scalar eigenmasses are positive
    and resolved. Outputs transform with the phase of delta.
    """
    require(all(math.isfinite(float(x)) and float(x)>0 for x in (m2,a,b)),
            'Positive finite squared masses required')
    delta=complex(delta)
    require(math.isfinite(delta.real) and math.isfinite(delta.imag), 'Invalid scalar mixing')
    sep=math.hypot(a-b,2*abs(delta))
    require(sep>1e-10*max(a,b), 'Unresolved scalar eigenvalue separation')
    lo=(a+b-sep)/2
    hi=(a+b+sep)/2
    require(lo>0 and hi>0, 'Tachyonic scalar eigenvalue')
    D=divdiff(m2,lo,hi)
    # Derivatives of the divided difference wrt its upper and lower nodes.
    d_hi=(b0_derivative(m2,hi)-D)/sep
    d_lo=(D-b0_derivative(m2,lo))/sep
    r=(a-b)/sep
    da=delta*(d_hi*(1+r)/2 + d_lo*(1-r)/2)
    db=delta*(d_hi*(1-r)/2 + d_lo*(1+r)/2)
    return delta*D,da,db


def split_input(inp:dict,d1:int,t:float)->dict:
    """Exact Batch-6 charge-split deformation multiplied by continuous t.

    t=1 reproduces the Batch-6 diagnostic exactly. No other matrix entry is
    changed; the result is *not* a gauge-invariant UV model deformation.
    """
    require(math.isfinite(t) and abs(t)<=1., 'Scalar splitting parameter outside [-1,1]')
    A=t5.decode(inp['scalar_A_gev2'],'A')
    full=t5.decode(k6.deform_scalars(inp,d1,'charge_split_diagnostic')['scalar_A_gev2'],'A1')
    z=A+t*(full-A)
    require(norm(z-z.conj().T)<1e-11*max(norm(z),1.), 'Split matrix not Hermitian')
    require(np.all(np.linalg.eigvalsh(z)>0), 'Split deformation gives tachyonic scalar')
    return {**inp, 'scalar_A_gev2': t5.encode(z)}


def analytic_group_derivative(inp:dict,tensors:dict,model:str,d1:int,
                              *,leading_mixing:bool=False)->dict:
    """Exact dB_g/dt at t=0 using projectors M^dag P_g.

    For leading_mixing=True, use Batch 7 derivative_D (neglect O(delta^3));
    otherwise differentiate the exact 2x2 mixed-scalar loop function.
    """
    A=t5.decode(inp['scalar_A_gev2'],'A')
    end=t5.decode(k6.deform_scalars(inp,d1,'charge_split_diagnostic')['scalar_A_gev2'],'end')
    M=t5.decode(inp['fermion_mass_gev'],'M')
    pairs=coupled_pairs(A,d1,model[0],int(model[2:]))
    groups=k6._groups(M)
    result=[]
    for g in groups:
        op=M.conj().T@g['P']
        group=np.zeros((3,3),complex)
        for p in pairs:
            s,u=p['s1'],p['s2']
            da=float((end[s,s]-A[s,s]).real)
            db=float((end[u,u]-A[u,u]).real)
            d=p['delta_gev2']
            if leading_mixing:
                _,fx,fy=t7.derivative_D(g['m2'],p['a_gev2'],p['b_gev2'])
                fx*=d;fy*=d
            else:
                _,fx,fy=exact_kernel_derivatives(g['m2'],p['a_gev2'],p['b_gev2'],d)
            T=tensors['C1'][:,:,s]@op@tensors['C2'][:,:,u].T
            group+=(da*fx+db*fy)*T
        result.append(group)
    require(len(result)==3, model+': unexpected heavy mass group count')
    return {'slope':np.stack(result),'heavy_masses_gev':[math.sqrt(g['m2']) for g in groups]}


def fit_with_derivative(ref0:np.ndarray,other0:np.ndarray,
                        ref1:np.ndarray,other1:np.ndarray)->dict:
    """Least-squares complex kappa(t), its exact derivative and nonparallel part.

    k(t)=<R(t),O(t)>/<R(t),R(t)> with <x,y>=sum(conj(x)*y).
    D=O-kR, so D'=O'-kR'-k'R. This treats complex derivatives correctly.
    """
    a=flatten(ref0);b=flatten(other0)
    ap=flatten(ref1);bp=flatten(other1)
    n=np.vdot(a,a)
    require(abs(n)>FLOOR,'Zero reference coefficient tensor')
    k=np.vdot(a,b)/n
    nprime=np.vdot(ap,a)+np.vdot(a,ap)
    kp=(np.vdot(ap,b)+np.vdot(a,bp)-k*nprime)/n
    d=b-k*a
    dp=bp-k*ap-kp*a
    require(all(np.isfinite(x) for x in (k,kp)) and np.all(np.isfinite(d)) and np.all(np.isfinite(dp)),
            'Nonfinite transverse factorization coefficient')
    return {'kappa':complex(k),'kappa_derivative':complex(kp),
            'baseline_transverse':d,'transverse_derivative':dp}


def point_fit(reference:np.ndarray,other:np.ndarray)->dict:
    a=flatten(reference);b=flatten(other)
    k=np.vdot(a,b)/np.vdot(a,a)
    d=b-k*a
    return {'kappa':complex(k),'transverse':d,'relative_nonproportionality':norm(d)/max(norm(b),FLOOR)}


def complex_value(z:complex)->dict[str,float]:
    return k6.value(complex(z))


def _csv(path:Path,rows:list[dict],fields:list[str])->None:
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def _slope(xs:dict[float,float],a:float,b:float,noise_floor:float)->float|None:
    u=xs[a];v=xs[b]
    if min(u,v)<=noise_floor or not (math.isfinite(u) and math.isfinite(v)):
        return None
    return math.log(v/u)/math.log(b/a)


def _report(p:dict)->str:
    lines=[
        '# Chat A Batch 8 — Controlled splitting convergence and transverse factorization',
        '',
        '**Source-conditional, physically uncertified: 0/16.** This uses the saved candidate Weyl tensors, '
        'the original scalar split diagnostic and no EFT matching. The split is a numerical '
        'U(1)em-preserving mass deformation, **not** a validated gauge-invariant UV operator.',
        '',
        '## Exact scalar-kernel differentiation',
        '',
        r'For a charge block A_q(t)=[[a_q+t Delta a_q, delta_q^*],[delta_q,b_q+t Delta b_q]], '
        r'let \(\lambda_\pm=(a+b\pm s)/2\), \(s=\sqrt{(a-b)^2+4|\delta|^2}\), '
        r'\(f(x)=B_0(0;M_g^2,x)\), and \(D=(f(\lambda_+)-f(\lambda_-))/s\).',
        '',
        r'\[K_q=\delta D,\qquad D_+=(f^{\prime}(\lambda_+)-D)/s,\quad '
        r'D_-=(D-f^{\prime}(\lambda_-))/s.\]',
        r'\[\partial_aK_q=\delta\left[D_+{1+(a-b)/s\over2}'
        r'+D_-{1-(a-b)/s\over2}\right],\]'
        r'\[\partial_bK_q=\delta\left[D_+{1-(a-b)/s\over2}'
        r'+D_-{1+(a-b)/s\over2}\right].\]',
        '',
        r'The exact slope is \(B_{g,1}=\sum_q(\Delta a_q\partial_aK_q+'
        r'\Delta b_q\partial_bK_q)T_{qg}\). The diagonal mass shifts and '
        'the complex scalar mixing are held fixed as specified, while t varies.',
        '',
        'Crucially, Batch 7 used the **leading-in-mixing** version of this slope. '
        'That is accurate for small mixing, but its difference from the exact slope produces an '
        r'\(O(t|\delta|^3)\) term. Only the exact finite-mixing derivative gives a clean '
        r'\(O(t^2)\) remainder.',
        '',
        '## Nonproportionality: the quantity relevant to factorization',
        '',
        r'Flatten all heavy-generation coefficient matrices into complex vectors '
        r'\(R(t)\) (reference) and \(B(t)\) (model). Define',
        '',
        r'\[\kappa(t)={\langle R(t),B(t)\rangle\over\langle R(t),R(t)\rangle},'
        r'\qquad D(t)=B(t)-\kappa(t)R(t).\]',
        '',
        r'\[\kappa_1={\langle R_1,B_0\rangle+\langle R_0,B_1\rangle'
        r'-\kappa_0(\langle R_1,R_0\rangle+\langle R_0,R_1\rangle)'
        r'\over\langle R_0,R_0\rangle},\]'
        r'\[D_1=B_1-\kappa_0R_1-\kappa_1R_0,\qquad'
        r' D(t)=D_0+tD_1+O(t^2).\]',
        '',
        'This is a complex least-squares projection, and the derivative includes the '
        'change in the *best-fit multiplier*. It compares the factorization-breaking '
        'signal itself rather than disguising it in the much larger total coefficient norm.',
        '',
        '## All-16 scalar-splitting convergence',
        '',
        'The first two errors below are normalized to the t=0 coefficient norm. '
        'They are evaluated at t=1/4. The last is the relative error of the '
        '**transverse increment** against its measured magnitude; n/a indicates a '
        'numerically vanishing transverse signal. A slope of 2 means quadratic convergence.',
        '',
        '| Model | exact-derivative error at t=1/4 | Batch7-leading-derivative error | '
        'O(t²) convergence slope | transverse increment error | t=1 nonparallel residual |',
        '|---|---:|---:|---:|---:|---:|']
    lookup={r['model']:r for r in p['model_summary']}
    at={(r['model'],r['t']):r for r in p['points']}
    for m in p['models']:
        r=lookup[m];s=at[(m,.25)]
        fmt=lambda x:'n/a' if x is None else f'{x:.5g}'
        lines.append(f'| {m} | {s["exact_slope_remainder_relative_base"]:.4e} | '
                     f'{s["leading_slope_remainder_relative_base"]:.4e} | '
                     f'{fmt(r["exact_remainder_order"])} | '
                     f'{fmt(s["transverse_increment_fractional_error"])} | '
                     f'{at[(m,1.)]["nonproportional_residual"]:.4e} |')
    lines += ['',
              '## Summary and interpretation',
              '',
              f'- All {p["model_count"]} baseline and t=1 endpoint tensors match the saved Batch-6 '
              'coefficients within the declared independent reconstruction tolerance.',
              f'- Largest exact finite-mixing slope / central-difference discrepancy: '
              f'{p["max_finite_difference_relative_error"]:.3e}.',
              f'- Maximum exact-derivative t=1/4 full-tensor remainder / baseline norm: '
              f'{p["max_quarter_remainder_relative_base"]:.3e}.',
              f'- Largest saved Batch-6 t=1 nonproportional residual: '
              f'{p["max_endpoint_nonproportionality"]:.3e}.',
              '',
              'A tiny residual of the complete coefficient tensor is not the same as a tiny '
              'relative error in its **transverse** component. The report includes both. '
              'Near-zero transverse channels are reported as noise-limited rather than assigned '
              'spurious power-law exponents.',
              '',
              'These results establish numerical convergence for the *specific* controlled '
              'splitting path. They are not a symbolic identity for an arbitrary T3 scalar '
              'potential. The physical 4-spinor/CG phase dictionary remains unresolved (0/16).',
              '']
    return '\n'.join(lines)


def analyse(root:Path,*,output_dir:Path|None=None)->dict:
    root=root.resolve();base=(root/BASE).resolve();out=(output_dir or root/TARGET).resolve()
    require(out.is_relative_to(base) and out!=base, 'Output must remain inside Chat A analysis directory')
    require(not out.exists() or (out.is_dir() and not any(out.iterdir())),
            'Output directory already populated: '+str(out))

    p6=root/BATCH06;p7=root/BATCH07;pstudy=root/k6.STUDY;p5=root/k6.BATCH05
    b6=read(p6/'batch06_factorization.json');m6=read(p6/'SOURCE_MANIFEST.json')
    h6=read(p6/'HANDOFF.json');b7=read(p7/'batch07_theorem.json')
    m7=read(p7/'SOURCE_MANIFEST.json');h7=read(p7/'HANDOFF.json')
    study=read(pstudy/'all16_summary.json')
    require(b6.get('format')=='t3-chatA-batch06-heavy-kernel-factorization-v1'
            and b7.get('format')=='t3-chatA-batch07-factorization-theorem-v1',
            'Batch-6/7 schema drift')
    require(m6.get('format')=='t3-chatA-batch06-manifest-v1'
            and m7.get('format')=='t3-chatA-batch07-manifest-v1', 'Batch-6/7 manifest schema drift')
    require(h6.get('format')=='t3-chatA-batch06-handoff-v1'
            and h7.get('format')=='t3-chatA-batch07-handoff-v1', 'Batch-6/7 handoff schema drift')
    require(all(x.get('physical_normalisation_certified') is False for x in (b6,b7,m6,m7,h6,h7,study)),
            'Physical certification status changed: refusing provisional analysis')
    require(b6.get('safe_for_physical_source_certification') is False and
            b7.get('safe_for_physical_source_certification') is False and
            b7.get('physical_source_certified_count')==0 and
            study.get('counts',{}).get('physical_source_certified')==0,
            'Uncertified calculation was mislabelled physically certified')
    require(h6.get('source_manifest_sha256')==k6.digest(p6/'SOURCE_MANIFEST.json')
            and h7.get('manifest_sha256')==k6.digest(p7/'SOURCE_MANIFEST.json'),
            'Handoff manifest hash changed')
    w=b6['source_sha256']
    require(m6['source_sha256']==w and
            k6.digest(pstudy/'all16_summary.json')==w['study'] and
            k6.digest(p5/'batch05_synthesis.json')==w['batch05'] and
            k6.digest(p5/'SOURCE_MANIFEST.json')==w['batch05_manifest'],
            'Batch-1/5/6 source hash chain drift')
    require(b7['source_sha256']['batch06_factorization']==k6.digest(p6/'batch06_factorization.json')
            and b7['source_sha256']['batch06_manifest']==k6.digest(p6/'SOURCE_MANIFEST.json')
            and m7['upstream_sha256']==b7['source_sha256'],
            'Batch-7 source hash chain drift')
    models=[r['model'] for r in model_dictionary()]
    require(len(models)==16 and set(b6['models'])==set(b7['models'])==set(models),
            '16-model source catalogue mismatch')
    require(len(b6['coefficient_comparisons'])==80 and len(b7['correction_tests'])==16,
            'Incomplete earlier calculations')
    saved=b6['kernel_coefficients']
    require(set(saved)==set(k6.SCENARIOS),'Batch-6 scenarios altered')
    for scenario in ('baseline','charge_split_diagnostic'):
        require(set(saved[scenario])==set(models),'Batch-6 endpoint missing model')
    raw={};srcinfo={};savedcheck=[]
    for model in models:
        record=read(pstudy/'models'/f'{model}.json')
        require(record.get('physical_normalisation_certified') is False and
                record.get('status')=='conditional_uv_evaluated',
                model+': provisional benchmark status changed')
        source=record['source']
        inp_path=t5._resolve_local(root,source['input'])
        benchmark=t5._resolve_local(root,source['benchmark'])
        require(k6.digest(inp_path)==source['input_sha256']==m6['models'][model]['input_sha256'] and
                k6.digest(benchmark)==source['benchmark_sha256']==m6['models'][model]['benchmark_sha256'] and
                k6.digest(pstudy/'models'/f'{model}.json')==m6['models'][model]['study_model_sha256'],
                'Source SHA drift '+model)
        inp=read(inp_path);frozen=read(benchmark)
        tensor=t5.extract_source_tensors(inp,frozen,model)
        d1=next(row for row in model_dictionary() if row['model']==model)['SU2_dimensions']['S1']
        lam=frozen['base_state']['ordinary']['lambdaT3']
        lam=abs(complex(lam['real'],lam['imag'])) if isinstance(lam,dict) else abs(complex(lam))
        require(lam>0 and math.isfinite(lam),'Missing scalar mixing '+model)
        controlled=t5.scalar_mixing_scale(inp,d1,0.10/lam)
        base_coeff=k6.coefficient_groups(controlled,tensor,model,d1)
        base=np.stack([g['B'] for g in base_coeff['groups']])
        expected=np.stack([t5.decode(g['B'],model+' baseline') for g in
                           saved['baseline'][model]['heavy_groups']])
        base_err=norm(base-expected)/max(norm(base),FLOOR)
        require(base_err<5e-8,model+': baseline coefficient reconstruction changed')
        terminal=k6.coefficient_groups(split_input(controlled,d1,1.),tensor,model,d1)
        end=np.stack([g['B'] for g in terminal['groups']])
        end_saved=np.stack([t5.decode(g['B'],model+' saved t=1') for g in
                            saved['charge_split_diagnostic'][model]['heavy_groups']])
        endpoint_err=norm(end-end_saved)/max(norm(end),FLOOR)
        require(endpoint_err<5e-8,model+': t=1 no longer reproduces Batch 6')
        derivative=analytic_group_derivative(controlled,tensor,model,d1)['slope']
        leading=analytic_group_derivative(controlled,tensor,model,d1,leading_mixing=True)['slope']
        h=1e-4
        plus=k6.coefficient_groups(split_input(controlled,d1,h),tensor,model,d1)
        minus=k6.coefficient_groups(split_input(controlled,d1,-h),tensor,model,d1)
        plus=np.stack([g['B'] for g in plus['groups']])
        minus=np.stack([g['B'] for g in minus['groups']])
        finite_diff=(plus-minus)/(2*h)
        # Compare relative to the derivative signal; allow numerical differentiation noise.
        fd_error=norm(finite_diff-derivative)/max(norm(derivative),1e-13*norm(base),FLOOR)
        require(fd_error<5e-4,model+f': exact split derivative inconsistent with central difference: {fd_error:g}')
        raw[model]={'base':base,'slope':derivative,'leading':leading,'end':end,
                    'source':controlled,'tensors':tensor,'d1':d1}
        srcinfo[model]={'input_sha256':source['input_sha256'],
                        'benchmark_sha256':source['benchmark_sha256'],
                        'study_model_sha256':m6['models'][model]['study_model_sha256']}
        savedcheck.append({'model':model,'baseline_saved_reconstruction_relative_error':base_err,
                           'endpoint_saved_reconstruction_relative_error':endpoint_err,
                           'central_difference_slope_relative_error':fd_error,
                           'exact_vs_leading_slope_relative_difference':norm(derivative-leading)/max(norm(derivative),FLOOR)})
    fit_slopes={m:fit_with_derivative(raw[REFERENCE]['base'],raw[m]['base'],
                                      raw[REFERENCE]['slope'],raw[m]['slope']) for m in models}
    points=[];summary=[]
    for t in TS:
        evaluated={}
        for model in models:
            obj=raw[model]
            if t==0.:arr=obj['base']
            elif t==1.:arr=obj['end']
            else:
                z=k6.coefficient_groups(split_input(obj['source'],obj['d1'],t),
                                        obj['tensors'],model,obj['d1'])
                arr=np.stack([g['B'] for g in z['groups']])
            evaluated[model]=arr
        ref=evaluated[REFERENCE]
        for model in models:
            obj=raw[model];a=obj['base'];b=evaluated[model]
            scale=max(norm(a),FLOOR)
            f=point_fit(ref,b)
            fit=fit_slopes[model]
            t0=fit['baseline_transverse']
            t1=fit['transverse_derivative']
            signal=norm(f['transverse']-t0)
            defect_remainder=norm(f['transverse']-(t0+t*t1))
            relative_signal=signal/max(norm(flatten(a)),FLOOR)
            # Ignore meaningless relative errors on transverse channels below 1e-9
            # of the full tensor (A family has nearly exact factorization).
            detectable=relative_signal>1e-9
            points.append({'model':model,'t':t,
                'exact_slope_remainder_relative_base':norm(b-(a+t*obj['slope']))/scale,
                'leading_slope_remainder_relative_base':norm(b-(a+t*obj['leading']))/scale,
                'change_relative_base':norm(b-a)/scale,
                'nonproportional_residual':f['relative_nonproportionality'],
                'transverse_increment_relative_base':relative_signal,
                'transverse_predicted_increment_relative_base':t*norm(t1)/scale,
                'transverse_remainder_relative_base':defect_remainder/scale,
                'transverse_increment_fractional_error':defect_remainder/max(signal,FLOOR) if detectable else None,
                'transverse_detectable':detectable,
                'kappa_real':f['kappa'].real,'kappa_imag':f['kappa'].imag})
    by={(r['model'],r['t']):r for r in points}
    for model in models:
        r={t:by[(model,t)]['exact_slope_remainder_relative_base'] for t in TS}
        remord=_slope(r,.125,.25,1e-13)
        trans={t:by[(model,t)]['transverse_increment_relative_base'] for t in TS}
        transord=_slope(trans,.125,.25,1e-10)
        modelpoints=[by[(model,t)] for t in TS]
        summary.append({'model':model,
                        'exact_remainder_order':remord,
                        'transverse_signal_order':transord,
                        'baseline_nonproportionality':by[(model,0.)]['nonproportional_residual'],
                        'endpoint_nonproportionality':by[(model,1.)]['nonproportional_residual'],
                        'max_full_tensor_remainder':max(x['exact_slope_remainder_relative_base'] for x in modelpoints),
                        'kappa0':complex_value(fit_slopes[model]['kappa']),
                        'kappa1':complex_value(fit_slopes[model]['kappa_derivative']),
                        'baseline_transverse_relative':norm(fit_slopes[model]['baseline_transverse'])/max(norm(raw[model]['base']),FLOOR),
                        'transverse_derivative_relative':norm(fit_slopes[model]['transverse_derivative'])/max(norm(raw[model]['base']),FLOOR)})
    for s in summary:
        order=s['exact_remainder_order']
        if order is not None:
            require(1.5<order<2.5,s['model']+f': exact splitting remainder not O(t²): slope={order:g}')
        if s['transverse_signal_order'] is not None:
            require(.6<s['transverse_signal_order']<1.4,
                    s['model']+': transverse splitting signal not first-order in t')
    payload={'format':'t3-chatA-batch08-splitting-convergence-v1',
             'models':models,'model_count':16,'sampled_t':list(TS),
             'source_conditional_only':True,'physical_normalisation_certified':False,
             'safe_for_physical_source_certification':False,'physical_source_certified_count':0,
             'reference_model':REFERENCE,
             'source_sha256':{'batch06':k6.digest(p6/'batch06_factorization.json'),
                              'batch06_manifest':k6.digest(p6/'SOURCE_MANIFEST.json'),
                              'batch07':k6.digest(p7/'batch07_theorem.json'),
                              'batch07_manifest':k6.digest(p7/'SOURCE_MANIFEST.json'),
                              'study':w['study']},
             'model_sources':srcinfo,'points':points,'model_summary':summary,
             'derivative_checks':savedcheck,
             'max_finite_difference_relative_error':max(x['central_difference_slope_relative_error'] for x in savedcheck),
             'max_quarter_remainder_relative_base':max(by[(m,.25)]['exact_slope_remainder_relative_base'] for m in models),
             'max_endpoint_nonproportionality':max(by[(m,1.)]['nonproportional_residual'] for m in models)}
    output_report=_report(payload)
    out.mkdir(parents=True,exist_ok=True)
    (out/'BATCH08_SPLITTING_CONVERGENCE.md').write_text(output_report,encoding='utf-8')
    (out/'batch08_convergence.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    _csv(out/'batch08_scaling_points.csv',points,
         ['model','t','change_relative_base','exact_slope_remainder_relative_base',
          'leading_slope_remainder_relative_base','nonproportional_residual',
          'transverse_increment_relative_base','transverse_predicted_increment_relative_base',
          'transverse_remainder_relative_base','transverse_increment_fractional_error',
          'transverse_detectable','kappa_real','kappa_imag'])
    _csv(out/'batch08_model_summary.csv',summary,
         ['model','exact_remainder_order','transverse_signal_order',
          'baseline_nonproportionality','endpoint_nonproportionality',
          'max_full_tensor_remainder','kappa0','kappa1',
          'baseline_transverse_relative','transverse_derivative_relative'])
    _csv(out/'batch08_derivative_validation.csv',savedcheck,
         ['model','baseline_saved_reconstruction_relative_error',
          'endpoint_saved_reconstruction_relative_error',
          'central_difference_slope_relative_error','exact_vs_leading_slope_relative_difference'])
    # Optional independent, non-interactive figures; the exact CSV/JSON diagnostics
    # are always emitted, even without matplotlib installed.
    pictures=[]
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        pass
    else:
        chosen=['A_-2','B_-1','C_+1','D_+2','E_-2']
        fig,ax=plt.subplots(figsize=(7.4,4.7))
        for m in chosen:
            xs=[x['t'] for x in points if x['model']==m and x['t']>0]
            ys=[x['exact_slope_remainder_relative_base'] for x in points if x['model']==m and x['t']>0]
            ax.loglog(xs,ys,marker='o',label=m)
        ax.set(xlabel='Split fraction t',ylabel='Exact tangent remainder / ||B(0)||',
               title='Scalar-splitting quadratic remainder (conditional)')
        ax.grid(True,which='both',alpha=.2);ax.legend(fontsize='small')
        fig.tight_layout();fig.savefig(out/'batch08_quadratic_convergence.png',dpi=160)
        plt.close(fig);pictures.append('batch08_quadratic_convergence.png')

        fig,ax=plt.subplots(figsize=(7.4,4.7))
        for m in ['B_-1','C_+1','D_+2','E_-2']:
            data=[x for x in points if x['model']==m and x['t']>0]
            ax.plot([x['t'] for x in data],
                    [x['transverse_increment_relative_base'] for x in data],
                    marker='o',label=m)
        ax.set(xlabel='Split fraction t',ylabel='Transverse increment / ||B(0)||',
               title='Onset of factorization breaking (conditional)')
        ax.grid(True,alpha=.2);ax.legend(fontsize='small')
        fig.tight_layout();fig.savefig(out/'batch08_transverse_response.png',dpi=160)
        plt.close(fig);pictures.append('batch08_transverse_response.png')
    manifest={'format':'t3-chatA-batch08-manifest-v1',
              'physical_normalisation_certified':False,
              'source_sha256':payload['source_sha256'],
              'generated_files':['BATCH08_SPLITTING_CONVERGENCE.md','batch08_convergence.json',
                                 'batch08_scaling_points.csv','batch08_model_summary.csv',
                                 'batch08_derivative_validation.csv']+pictures}
    (out/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    handoff={'format':'t3-chatA-batch08-handoff-v1',
             'all16_evaluated':16,'physical_normalisation_certified':False,
             'physically_certified_count':0,'report':'BATCH08_SPLITTING_CONVERGENCE.md',
             'source_manifest_sha256':k6.digest(out/'SOURCE_MANIFEST.json')}
    (out/'HANDOFF.json').write_text(json.dumps(handoff,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return payload
