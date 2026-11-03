"""Chat A Batch 7: source-conditional factorization theorem and first corrections.

Read-only on Batches 1/5/6 and source inputs. No production UV solver calls,
no EFT and no physical Matchete/RGBeta normalization certificate.
"""
from __future__ import annotations

import ast
import csv
import json
import math
from pathlib import Path
from typing import Any
import numpy as np

from Numerical.lo_patterns_16 import kernel_factorization as k6
from Numerical.lo_patterns_16 import texture_synthesis as t5
from Numerical.lo_patterns_16.analytic_explanation import (
    b0, b0_derivative, divdiff, coupled_pairs, scalar_kernels
)
from Numerical.lo_certification.direct_uv.certification import charge_bases
from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary

BASE=k6.BASE
STUDY=k6.STUDY
BATCH05=k6.BATCH05
BATCH06=k6.TARGET
TARGET=BASE/'analysis_batch07_01'
AUDIT=Path('output/general_observable/lo_direct_one_loop/source_audits/charge_resolved_loop_weights.json')
REF=k6.REFERENCE_MODEL
SCENARIOS=k6.SCENARIOS
TINY=1e-100


class TheoremError(ValueError):
    pass


def require(ok:bool,msg:str)->None:
    if not ok:raise TheoremError(msg)


def load(p:Path)->dict:
    require(p.is_file(),'Missing required source: '+str(p))
    x=json.loads(p.read_text(encoding='utf-8-sig'))
    require(isinstance(x,dict),'Invalid JSON object '+str(p))
    return x


def digest(p:Path)->str:
    return k6.digest(p)


def fro(x:np.ndarray)->float:
    return float(np.linalg.norm(x))


def flatten(groups:list[dict],key:str='B')->np.ndarray:
    return np.concatenate([g[key].ravel() for g in groups])


def derivative_D(m2:float,a:float,b:float)->tuple[float,float,float]:
    """D=f[a,b], dD/da and dD/db, where f=B0(0;m²,s²).

    This branch requires nondegenerate a,b; near degeneracy is supported by
    the exact B0-divided difference elsewhere, but not by this derivative form.
    """
    require(m2>0 and a>0 and b>0 and abs(a-b)>1e-5*max(a,b),
            'Nondegenerate positive scalar masses required for derivative formula')
    D=divdiff(m2,a,b)
    fa=b0_derivative(m2,a);fb=b0_derivative(m2,b)
    da=(fa-D)/(a-b)
    db=(D-fb)/(a-b)
    require(all(math.isfinite(v) for v in (D,da,db)), 'Nonfinite scalar divided difference')
    return D,da,db


def cubic_E(m2:float,a:float,b:float)->float:
    """Coefficient of delta*|delta|² in off-diagonal f(A): f[a,a,b,b]."""
    D,_,_=derivative_D(m2,a,b)
    return (b0_derivative(m2,a)+b0_derivative(m2,b)-2*D)/(a-b)**2


def kernel_series(m2:float,a:float,b:float,delta:complex)->tuple[complex,complex]:
    D,_,_=derivative_D(m2,a,b)
    lin=delta*D
    cubic=delta*(D+abs(delta)**2*cubic_E(m2,a,b))
    return lin,cubic


def scalar_split_correction(inp:dict,tensors:dict,model:str,d1:int)->dict:
    """First-order charge-splitting correction at baseline scalar masses.

    Delta B = sum_q delta_q*(da d_a D + db d_b D) * (C1 MdagP C2^T)
    for each degenerate heavy mass eigenspace. No UV loop re-evaluation.
    """
    A=t5.decode(inp['scalar_A_gev2'],'A');M=t5.decode(inp['fermion_mass_gev'],'M')
    Aq=t5.decode(k6.deform_scalars(inp,d1,'charge_split_diagnostic')['scalar_A_gev2'],'A_split')
    pairs=coupled_pairs(A,d1,model[0],int(model[2:]))
    original=[]; correction=[]
    for g in k6._groups(M):
        op=M.conj().T@g['P'];B=np.zeros((3,3),complex);db=np.zeros((3,3),complex)
        for p in pairs:
            s,u=p['s1'],p['s2'];a,b=p['a_gev2'],p['b_gev2'];delta=p['delta_gev2']
            d,dda,ddb=derivative_D(g['m2'],a,b)
            da=float((Aq[s,s]-A[s,s]).real);ds=float((Aq[u,u]-A[u,u]).real)
            factor=tensors['C1'][:,:,s] @ op @ tensors['C2'][:,:,u].T
            x,_,_=scalar_kernels(g['m2'],a,b,delta)
            B+=x*factor
            db+=delta*(da*dda+ds*ddb)*factor
        original.append(B);correction.append(db)
    return {'base':np.stack(original),'correction':np.stack(correction),
            'linearized':np.stack(original)+np.stack(correction)}


def large_mix_series(inp:dict,tensors:dict,model:str,d1:int,multiplier:float=20.)->dict:
    """Exact source CG tensors with leading and cubic scalar-mixing kernels."""
    A=t5.decode(inp['scalar_A_gev2'],'A');M=t5.decode(inp['fermion_mass_gev'],'M')
    pairs=coupled_pairs(A,d1,model[0],int(model[2:]))
    linear=[];cubic=[]
    for g in k6._groups(M):
        op=M.conj().T@g['P'];bl=np.zeros((3,3),complex);bc=np.zeros((3,3),complex)
        for p in pairs:
            s,u=p['s1'],p['s2'];a,b=p['a_gev2'],p['b_gev2'];delta=p['delta_gev2']*multiplier
            lin,third=kernel_series(g['m2'],a,b,delta)
            factor=tensors['C1'][:,:,s] @ op @ tensors['C2'][:,:,u].T
            bl+=lin*factor;bc+=third*factor
        linear.append(bl);cubic.append(bc)
    return {'linear':np.stack(linear),'cubic':np.stack(cubic)}


def compare_coefficient_arrays(reference:np.ndarray,other:np.ndarray)->dict:
    require(reference.shape==other.shape and reference.size and fro(reference)>TINY,
            'Coefficient tensor shape or reference invalid')
    k=np.vdot(reference,other)/np.vdot(reference,reference)
    err=fro(other-k*reference)/max(fro(other),TINY)
    return {'kappa':k6.value(complex(k)), 'residual':err}


def _simple_expr(expr:str)->complex:
    """Evaluate a small, explicit source-CG arithmetic grammar, never eval()."""
    def walk(t):
        if isinstance(t,ast.Constant) and type(t.value) in (int,float,complex):return complex(t.value)
        if isinstance(t,ast.Name) and t.id=='I':return 1j
        if isinstance(t,ast.UnaryOp) and isinstance(t.op,(ast.UAdd,ast.USub)):
            z=walk(t.operand);return z if isinstance(t.op,ast.UAdd) else -z
        if isinstance(t,ast.BinOp) and isinstance(t.op,(ast.Add,ast.Sub,ast.Mult,ast.Div,ast.Pow)):
            a,b=walk(t.left),walk(t.right)
            if isinstance(t.op,ast.Add):return a+b
            if isinstance(t.op,ast.Sub):return a-b
            if isinstance(t.op,ast.Mult):return a*b
            if isinstance(t.op,ast.Div):return a/b
            require(abs(b.imag)<1e-14 and abs(b.real)<=8,'Unsupported CG exponent')
            return a**b
        if isinstance(t,ast.Call) and isinstance(t.func,ast.Name) and t.func.id=='sqrt' and len(t.args)==1 and not t.keywords:
            return walk(t.args[0])**.5
        raise TheoremError('Unsupported source-CG arithmetic '+expr)
    z=walk(ast.parse(expr,mode='eval').body)
    require(math.isfinite(z.real) and math.isfinite(z.imag),'Nonfinite CG expression')
    return z


def audit_weight_lookup(payload:dict,catalogue:dict[str,dict])->dict:
    require(payload.get('status')=='charge_resolved_group_coefficients' and payload.get('count')==16,
            'Source CG weight audit schema/count drift')
    require(payload.get('contraction')=='Y1 * Y2 * conjugate(T3MixCG) using same physical internal weight indices',
            'Source CG weight contraction changed')
    rows={}
    for item in payload.get('models',[]):
        m=f"{item['class']}_{int(item['model']['alpha']):+d}"
        require(m in catalogue and m not in rows,'Unknown/duplicate source CG model '+m)
        dim=catalogue[m]['SU2_dimensions']
        require(tuple(item['model'][k] for k in ('dS1','dS2','dF'))==
                tuple(dim[k] for k in ('S1','S2','F')),'Source CG representation drift '+m)
        cs=charge_bases(m[0],int(m[2:]))
        for ch in item['channels']:
            i={name:int(ch['indices_weight_basis_one_based'][name])-1 for name in ('S1','S2','F')}
            require(all(0<=i[name]<dim[name] for name in ('S1','S2','F')),
                    'Source CG component index invalid '+m)
            charge=float(_simple_expr(ch['charge']).real)
            require(all(abs(q-charge)<1e-9 for q in (
                cs['scalar'][i['S1']],cs['scalar'][dim['S1']+i['S2']],
                (dim['F']-1)/2-i['F']+(int(m[2:])+1)/2)),
                'Source CG charge mismatch '+m)
            product=_simple_expr(ch['y1'])*_simple_expr(ch['y2'])*_simple_expr(ch['mixing']).conjugate()
            audited=_simple_expr(ch['complex_weight'])
            require(abs(product-audited)<1e-10*max(abs(product),abs(audited),1e-100),
                    'Source CG product/orientation mismatch '+m)
        rows[m]=item
    require(set(rows)==set(catalogue),'Incomplete source CG audit 16-model coverage')
    return rows


def expected_model_weights(audit:dict,models:list[str])->dict:
    catalogue={r['model']:r for r in model_dictionary()}
    require(set(models)==set(catalogue) and len(models)==16,'Expected original 16-model catalogue')
    checked=audit_weight_lookup(audit,catalogue)
    by={}
    for model in models:
        x=checked[model]
        active=[z for z in x['channels'] if z['nonzero']]
        require(len(active)==(2 if model[0]=='E' else 1),model+': source charge multiplicity changed')
        weights=[_simple_expr(z['complex_weight']) for z in active]
        summed=sum(weights)
        require(abs(summed-_simple_expr(x['total_group_coefficient']))<1e-12,
                model+': source CG total inconsistent')
        by[model]={'family':model[0], 'alpha':int(model[2:]),
                   'charges':[float(_simple_expr(z['charge']).real) for z in active],
                   'weight_formulas':[z['complex_weight'] for z in active],
                   'weight_total_formula':x['total_group_coefficient'],
                   'weight_total':k6.value(summed),'weight_magnitude':float(abs(summed)),
                   'active_charge_count':len(active)}
    ref=by[REF]['weight_magnitude']
    require(ref>0,'Reference group weight vanished')
    for x in by.values():x['magnitude_relative_to_ref']=x['weight_magnitude']/ref
    return by


def _flat_saved(src:dict,sc:str,model:str)->np.ndarray:
    rows=src['kernel_coefficients'][sc][model]['heavy_groups']
    require(len(rows)==3,'Missing heavy groups '+model+' '+sc)
    return np.stack([t5.decode(g['B'],model+' B',(3,3)) for g in rows])


def _csv(p:Path,rows:list[dict],cols:list[str])->None:
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)


def _report(p:dict)->str:
    rows=p['cg_and_factorization']; corrections=p['correction_tests']
    maximum=lambda sc:max(r['coefficient_residual'] for r in p['scenario_summary'] if r['scenario']==sc)
    txt=[
        '# Chat A Batch 7 — Conditional factorization theorem and first corrections','',
        '**Status:** derived and numerically checked within the frozen candidate Weyl dictionary. '
        'Physical Matchete/RGBeta normalization and phases are **not certified (0/16)**. '
        'No EFT matching or production-code changes.','',
        '## 1. Source-convention input and exact matrix identity','',
        r'The candidate source vertices are \(Z=Y_1^*C_1\) and \(Y=Y_2^*C_2\), '
        r'and the restricted scalar quadratic has \(B_{SS}=0\) with isolated 2×2 equal-charge blocks. '
        r'For a full heavy-degenerate eigenspace projector \(P_g\),', '',
        r'\[\boxed{B^{(R)}_{g,rs}=\sum_q K_q(M_g^2)\,C^{(1)}_{r,I,s_q}'
        r'(M_F^\dagger P_g)_{IJ}C^{(2)}_{s,J,t_q}}\]', '',
        r'\[K_q=\delta_q\,{f_g(\lambda_{q+})-f_g(\lambda_{q-})\over'
        r'\lambda_{q+}-\lambda_{q-}},\quad f_g(x)=B_0(0;M_g^2,x).\]', '',
        r'\[\boxed{M_\nu^{(R)}={1\over16\pi^2}\sum_g\operatorname{Sym}'
        r'\left(Y_1^*B_g^{(R)}Y_2^{*T}\right)},\quad'
        r'\operatorname{Sym}(X)=X+X^T.\]', '',
        'The coefficients B_g are 3×3 **source-generation** matrices, not neutrino-mass matrices. '
        'Their independent heavy-group projections avoid arbitrary phase choices inside degenerate fermion eigenspaces.','',
        '## 2. Fixed-spectrum factorization theorem','',
        r'**Theorem.** If for all heavy groups \(B_g^{(R)}=\kappa_R B_g^{(0)}\) '
        r'with one common complex \(\kappa_R\), then for *any* common 3×3 Yukawa pair', '',
        r'\[M_\nu^{(R)}(Y_1,Y_2)=\kappa_R M_\nu^{(0)}(Y_1,Y_2),\]','',
        'as follows immediately from linearity of the symmetrized bilinear contraction. '
        'The statement is exact for the specified coefficient tensors; a small measured '
        'residual only establishes approximate factorization at the checked spectra.','',
        r'For errors \(E_g=B_g^{(R)}-\kappa_RB_g^{(0)}\), the Frobenius bound is', '',
        r'\[\|\Delta M_\nu\|_F\le {2\over16\pi^2}\|Y_1\|_2\|Y_2\|_2'
        r'\sum_g\|E_g\|_F.\]', '',
        'This is a deterministic bound for arbitrary common Yukawa matrices, '
        'not a fitted-texture observation. The theorem does not assert a physical equality '
        'between differently normalized theories.','',
        '## 3. Source CG-weight table (16 models)','',
        'The listed weights use the **source audit** contraction '
        '`Y1 * Y2 * conjugate(T3MixCG)`, not a certified four-spinor vertex rule. '
        'The source-audit signed phase must **not** automatically be identified '
        'with the fitted Weyl coefficient phase κ.','',
        '| Model | Q of active S1/S2 pairs | Source weight(s) | Total weight | '
        '|weight|/|A_-2| | measured |κ| | baseline tensor ε |',
        '|---|---|---|---|---:|---:|---:|']
    for r in rows:
        txt.append(f'| {r["model"]} | '+', '.join(f'{q:g}' for q in r['charges'])+' | '+
                   ', '.join(r['weight_formulas'])+f' | {r["weight_total_formula"]} | '
                   f'{r["magnitude_relative_to_ref"]:.9f} | '
                   f'{r["measured_kappa_magnitude"]:.9f} | {r["baseline_residual"]:.3e} |')
    txt.extend(['','**CG magnitude versus measured tensor scale:** largest absolute ratio mismatch '
                f'{p["max_cg_kappa_magnitude_disagreement"]:.3e}. '
                'The finite mismatch reflects conditional scalar-kernel effects/numerical precision '
                'at this benchmark; neither source audit nor fitted κ certifies a physical phase.','',
                '## 4. First factorization-breaking terms','',
                r'Write \(a_q=a+\epsilon_{a,q}\), \(b_q=b+\epsilon_{b,q}\) and '
                r'\(\delta_q\) for the cross-scalar mixing. '
                r'For nondegenerate unmixed scalar squared masses define', '',
                r'\[D_g(a,b)={f_g(a)-f_g(b)\over a-b},\quad '
                r'D_{g,a}={f_g\!\prime(a)-D_g\over a-b},\quad '
                r'D_{g,b}={D_g-f_g\!\prime(b)\over a-b}.\]','',
                r'The exact off-diagonal two-scalar loop kernel has the expansion', '',
                r'\[\boxed{K_{qg}=\delta_q\{D_g+'
                r'\epsilon_{a,q}D_{g,a}+\epsilon_{b,q}D_{g,b}'
                r'+|\delta_q|^2E_g+\cdots\}},\]', '',
                r'\[E_g={f_g\!\prime(a)+f_g\!\prime(b)-2D_g\over(a-b)^2}'
                r'=f_g[a,a,b,b].\]', '',
                'The terms retained in this displayed expression are first order in the '
                'charge-dependent diagonal shifts and cubic in mixing; mixed terms such as '
                r'\(\delta|\delta|^2\epsilon\) and second-order splittings are omitted.','',
                'For each heavy group and charge, let '
                r'\(T_{qg,rs}=C^{(1)}_{r,I,s_q}(M_F^\dagger P_g)_{IJ}C^{(2)}_{s,J,t_q}\). '
                'Then the leading splitting correction is', '',
                r'\[\boxed{\Delta B_{g,rs}^{\rm split} '
                r'=\sum_q\delta_q(\epsilon_{a,q}D_{g,a}'
                r'+\epsilon_{b,q}D_{g,b})T_{qg,rs}},\]', '',
                r'and the first nonlinear-mixing term is '
                r'\(\Delta B_g^{\rm cubic}=\sum_q\delta_q|\delta_q|^2E_gT_{qg}\). '
                'These corrections need not be proportional to the reference tensor '
                'if they weight different charge sectors or heavy groups differently.','',
                '## 5. Numerical evidence and perturbative checks','',
                '| Scenario | Largest B_g factorization residual (all 16) |',
                '|---|---:|'])
    for sc in SCENARIOS:
        xs=[x for x in p['scenario_summary'] if x['scenario']==sc]
        worst=max(xs,key=lambda x:x['coefficient_residual'])
        txt.append(f'| {sc} (worst {worst["model"]}) | {worst["coefficient_residual"]:.6e} |')
    txt.extend(['','The following compares the analytic **truncated** first-correction '
                'formulae with exact saved Batch-6 heavy-group coefficient tensors. '
                'All errors are relative Frobenius errors; a perturbative approximation '
                'is not assumed exact under the deliberate stress tests.','',
                '| Model | Charge-split first-order error / exact | '
                '|ΔB_split| / |B_baseline| | Large-mix leading error | '
                'Large-mix cubic error |',
                '|---|---:|---:|---:|---:|'])
    for r in corrections:
        txt.append(f'| {r["model"]} | {r["charge_split_linearized_relative_error"]:.3e} | '
                   f'{r["split_correction_fraction_of_baseline"]:.3e} | '
                   f'{r["large_mix_linear_relative_error"]:.3e} | '
                   f'{r["large_mix_cubic_relative_error"]:.3e} |')
    txt.extend(['','## 6. Domain of validity and limitations','',
                '- The 16 input states use a restricted scalar quadratic: no scalar-B terms, '
                'and the same unmixed S1 and S2 masses for components of each multiplet.',
                '- With identical charge-sector kernels and compatible heavy-generation '
                'contractions, representation dependence can collapse to an overall '
                'candidate CG coefficient. Scalar degeneracy alone does **not** establish '
                'the necessary tensor identity; Batch 6 checked it numerically.',
                '- The heavy-generation tilt changes three eigenmasses while keeping the '
                'spectral projectors; it is a mass-matrix experiment, not a UV renormalization.',
                '- Large scalar mixing resums higher odd powers in δ; the cubic expansion '
                'is conditional on small enough |δ| and nondegenerate unmixed masses.',
                '- Charge splitting deliberately violates the multiplet-degenerate '
                'restricted mass assumption; it is **not** asserted to arise from an '
                'SU(2)-invariant T3 potential.',
                '- Frozen Yukawa differences explain the normalized spectra, but the '
                'CG coefficients in this report are not standalone physical amplitudes.',
                '- The source audit uses a separate CG contraction convention; its signs '
                'cannot be equated to the candidate Weyl coefficient signs without B1/B2.',
                '- No source physical normalization certification is claimed: **0/16**.', '',
                '## 7. Scientific synthesis','',
                'At matched Yukawas and within the original restricted scalar-sector '
                'assumptions, the saved heavy-generation coefficient tensors of all 16 '
                'models are nearly proportional, with baseline errors at most '
                f'{p["max_baseline_residual"]:.3e}. The candidate CG magnitude pattern '
                'is (A,B,D):1, C:1/√3, E:√2 relative to A_-2. '
                'Consequently the *normalized* neutrino mass spectrum is controlled by '
                'the common Yukawa flavour texture at the checked spectrum. '
                'The leading analytically identified violations are charge-dependent '
                'scalar mass shifts and nonlinear scalar mixing. This is a '
                '**source-conditional numerical-plus-algebraic statement**, '
                'not a universal UV theorem or EFT result.',''])
    return '\n'.join(txt)


def analyse(root:Path,*,output_dir:Path|None=None)->dict:
    root=root.resolve();owner=(root/BASE).resolve();out=(output_dir or root/TARGET).resolve()
    require(out.is_relative_to(owner) and out != owner,
            'Output must remain inside Chat A analysis directory')
    require(not out.exists() or (out.is_dir() and not any(out.iterdir())),
            'Output directory already populated: '+str(out))

    p6=root/BATCH06;p5=root/BATCH05;pstudy=root/STUDY
    src=load(p6/'batch06_factorization.json');manifest=load(p6/'SOURCE_MANIFEST.json')
    handoff=load(p6/'HANDOFF.json');study=load(pstudy/'all16_summary.json')
    audit=load(root/AUDIT)
    require(src.get('format')=='t3-chatA-batch06-heavy-kernel-factorization-v1' and
            manifest.get('format')=='t3-chatA-batch06-manifest-v1' and
            handoff.get('format')=='t3-chatA-batch06-handoff-v1', 'Batch 6 schema mismatch')
    require(src.get('physical_normalisation_certified') is False and
            manifest.get('physical_normalisation_certified') is False and
            handoff.get('physical_normalisation_certified') is False and
            src.get('safe_for_physical_source_certification') is False,
            'False/unclear physical certification status')
    require(study.get('physical_normalisation_certified') is False and
            study.get('counts',{}).get('physical_source_certified')==0,
            'Study physical certification status drift')
    require(handoff['source_manifest_sha256']==digest(p6/'SOURCE_MANIFEST.json'),
            'Batch 6 manifest hash chain drift')
    want=src['source_sha256']
    require(manifest.get('source_sha256')==want and
            digest(pstudy/'all16_summary.json')==want['study'] and
            digest(p5/'batch05_synthesis.json')==want['batch05'] and
            digest(p5/'SOURCE_MANIFEST.json')==want['batch05_manifest'],
            'Batch 1/5/6 upstream SHA chain drift')
    models=[r['model'] for r in model_dictionary()]
    require(len(models)==16 and sorted(src['models'])==sorted(models) and
            len(src['coefficient_comparisons'])==16*len(SCENARIOS),
            'Incomplete all-16 scenario result')
    weights=expected_model_weights(audit,models)
    cmp={(r['scenario'],r['model']):r for r in src['coefficient_comparisons']}
    require(set(cmp)=={(sc,m) for sc in SCENARIOS for m in models},
            'Missing/duplicated model–scenario coefficient rows')
    guard={}
    for m in models:
        study_model=load(pstudy/'models'/f'{m}.json')
        require(study_model.get('status')=='conditional_uv_evaluated' and
                study_model.get('physical_normalisation_certified') is False,
                'Source conditional result status drift '+m)
        f=study_model['source'];inp=t5._resolve_local(root,f['input'])
        frozen=t5._resolve_local(root,f['benchmark'])
        require(digest(inp)==f['input_sha256']==manifest['models'][m]['input_sha256'] and
                digest(frozen)==f['benchmark_sha256']==manifest['models'][m]['benchmark_sha256'] and
                digest(pstudy/'models'/f'{m}.json')==manifest['models'][m]['study_model_sha256'],
                'Source SHA drift '+m)
        guard[m]=(load(inp),load(frozen))

    rows=[];test_rows=[]
    for m in models:
        w=weights[m];baseline=cmp[('baseline',m)]
        kappa=baseline['absolute_scale']
        require(math.isfinite(kappa) and kappa>0 and baseline['coefficient_residual']>=0,
                'Malformed coefficient scale for '+m)
        rows.append({'model':m,**w,'measured_kappa_magnitude':kappa,
                     'baseline_residual':baseline['coefficient_residual'],
                     'weight_ratio_minus_measured_scale':w['magnitude_relative_to_ref']-kappa})
        inp,frozen=guard[m]
        ts=t5.extract_source_tensors(inp,frozen,m)
        dimension=next(r for r in model_dictionary() if r['model']==m)['SU2_dimensions']['S1']
        lam=frozen['base_state']['ordinary']['lambdaT3']
        lam=abs(complex(lam['real'],lam['imag'])) if isinstance(lam,dict) else abs(complex(lam))
        require(lam>0,'Zero scalar mixing strength '+m)
        controlled=t5.scalar_mixing_scale(inp,dimension,0.10/lam)
        original=k6.coefficient_groups(controlled,ts,m,dimension)
        B0=np.stack([g['B'] for g in original['groups']]);stored=_flat_saved(src,'baseline',m)
        require(fro(B0-stored)/max(fro(stored),TINY)<5e-8,
                'Baseline coefficient tensor drift '+m)

        split=scalar_split_correction(controlled,ts,m,dimension)
        s_exact=_flat_saved(src,'charge_split_diagnostic',m)
        m_series=large_mix_series(controlled,ts,m,dimension)
        l_exact=_flat_saved(src,'large_scalar_mixing',m)
        spliterr=fro(split['linearized']-s_exact)/max(fro(s_exact),TINY)
        linerr=fro(m_series['linear']-l_exact)/max(fro(l_exact),TINY)
        cubicerr=fro(m_series['cubic']-l_exact)/max(fro(l_exact),TINY)
        require(all(np.isfinite(x) and x>=0 for x in (spliterr,linerr,cubicerr)),
                'Nonfinite correction approximation '+m)
        test_rows.append({'model':m,'charge_split_linearized_relative_error':spliterr,
                          'split_correction_fraction_of_baseline':fro(split['correction'])/max(fro(B0),TINY),
                          'large_mix_linear_relative_error':linerr,
                          'large_mix_cubic_relative_error':cubicerr,
                          'cubic_improves_over_linear':cubicerr<=linerr*(1+1e-4)+1e-12})
        require(cubicerr<=linerr*(1+1e-3)+1e-10, 'Cubic mixing correction did not improve '+m)
        require(spliterr<0.10,'First-order charge-splitting approximation outside working range '+m)

    maxdiff=max(abs(x['weight_ratio_minus_measured_scale']) for x in rows)
    require(maxdiff<2e-5,'CG weight magnitudes no longer track coefficient scales')
    scenario_rows=[{'scenario':sc,'model':m,
                    'coefficient_residual':cmp[(sc,m)]['coefficient_residual'],
                    'max_group_fractional_spread':cmp[(sc,m)]['max_group_fractional_spread']}
                   for sc in SCENARIOS for m in models]
    payload={'format':'t3-chatA-batch07-factorization-theorem-v1',
             'source_conditional_only':True,'physical_normalisation_certified':False,
             'safe_for_physical_source_certification':False,
             'models':models,'model_count':16,'reference_model':REF,
             'source_sha256':{'batch06_factorization':digest(p6/'batch06_factorization.json'),
                              'batch06_manifest':digest(p6/'SOURCE_MANIFEST.json'),
                              'study':want['study'],'batch05':want['batch05'],
                              'cg_weight_audit':digest(root/AUDIT)},
             'cg_and_factorization':rows,'scenario_summary':scenario_rows,
             'correction_tests':test_rows,
             'max_cg_kappa_magnitude_disagreement':maxdiff,
             'max_baseline_residual':max(r['baseline_residual'] for r in rows),
             'theorem_scope':'Restricted scalar B=0, isolated equal-charge scalar blocks, '
                            'candidate Weyl vertices and same Yukawa pair at the tested spectra',
             'physical_source_certified_count':0}
    report=_report(payload)
    out.mkdir(parents=True,exist_ok=True)
    (out/'BATCH07_FACTORISATION_THEOREM.md').write_text(report,encoding='utf-8')
    (out/'batch07_theorem.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    _csv(out/'batch07_cg_weight_table.csv',rows,
         ['model','family','alpha','charges','weight_formulas','weight_total_formula',
          'magnitude_relative_to_ref','measured_kappa_magnitude','baseline_residual',
          'weight_ratio_minus_measured_scale'])
    _csv(out/'batch07_first_corrections.csv',test_rows,
         ['model','charge_split_linearized_relative_error','split_correction_fraction_of_baseline',
          'large_mix_linear_relative_error','large_mix_cubic_relative_error',
          'cubic_improves_over_linear'])
    _csv(out/'batch07_scenario_residuals.csv',scenario_rows,
         ['scenario','model','coefficient_residual','max_group_fractional_spread'])
    manifest_out={'format':'t3-chatA-batch07-manifest-v1',
                  'physical_normalisation_certified':False,
                  'upstream_sha256':payload['source_sha256'],
                  'generated_files':['BATCH07_FACTORISATION_THEOREM.md','batch07_theorem.json',
                                     'batch07_cg_weight_table.csv','batch07_first_corrections.csv',
                                     'batch07_scenario_residuals.csv']}
    (out/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest_out,indent=2)+'\n',encoding='utf-8')
    handoff={'format':'t3-chatA-batch07-handoff-v1','all16_derived':16,
             'physical_normalisation_certified':False,'physically_certified_count':0,
             'max_baseline_residual':payload['max_baseline_residual'],
             'max_cg_scale_disagreement':maxdiff,
             'report':'BATCH07_FACTORISATION_THEOREM.md',
             'manifest_sha256':digest(out/'SOURCE_MANIFEST.json')}
    (out/'HANDOFF.json').write_text(json.dumps(handoff,indent=2)+'\n',encoding='utf-8')
    return payload
