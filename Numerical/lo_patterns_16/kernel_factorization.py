"""Chat A Batch 6: source-conditional CG x heavy-loop coefficient factorization.

NO EFT. NO UV source changes. NO physical Matchete/RGBeta certification.
Computes Yukawa-INDEPENDENT rank-2 flavour coefficient tensors in each
heavy-mass eigenspace and charge channel from the saved candidate Weyl input.
All mutations are in-memory counterfactuals; all source files are read-only.
"""
from __future__ import annotations

import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.lo_patterns_16 import texture_synthesis as t5
from Numerical.lo_patterns_16.analytic_explanation import (
    LOOP, coupled_pairs, scalar_kernels, reconstruct,
)

BASE = Path('output/general_observable/lo_pattern_study16')
STUDY = BASE / 'provisional_01'
BATCH03 = BASE / 'analysis_batch03_01'
BATCH04 = BASE / 'analysis_batch04_01'
BATCH05 = BASE / 'analysis_batch05_01'
TARGET = BASE / 'analysis_batch06_01'
REFERENCE_MODEL = 'A_-2'
TEXTURES = ('A_-2', 'C_-1', 'C_+1', 'E_-2')
SCENARIOS = ('baseline', 'scalar_common_shift', 'heavy_generation_tilt',
             'large_scalar_mixing', 'charge_split_diagnostic')
EPS = 1e-100


class FactorizationError(ValueError):
    pass


def require(ok: bool, why: str) -> None:
    if not ok:
        raise FactorizationError(why)


def read(p: Path) -> dict:
    require(p.is_file(), 'Missing source file: ' + str(p))
    x = json.loads(p.read_text(encoding='utf-8-sig'))
    require(isinstance(x, dict), 'Expected JSON object: ' + str(p))
    return x


def digest(p: Path) -> str:
    return sha256(p.read_bytes()).hexdigest()


def norm(x: np.ndarray) -> float:
    return float(np.linalg.norm(x))


def rel(x: np.ndarray, y: np.ndarray) -> float:
    return norm(x-y) / max(norm(x), norm(y), EPS)


def value(z: complex) -> dict:
    return {'real': float(complex(z).real), 'imag': float(complex(z).imag)}


def _groups(M: np.ndarray) -> list[dict]:
    require(M.ndim == 2 and M.shape[0] == M.shape[1] and np.all(np.isfinite(M)),
            'Invalid fermion mass matrix')
    require(rel(M,M.T)<1e-10, 'Fermion mass matrix must be symmetric')
    h = M @ M.conj().T
    w,u = np.linalg.eigh((h+h.conj().T)/2.)
    require(np.min(w)>0 and np.min(w)>1e-15*np.max(w), 'Singular heavy mass spectrum')
    result = []
    for k,m2 in enumerate(w):
        if not result or abs(float(m2)-result[-1]['m2'])>1e-8*max(float(m2),result[-1]['m2']):
            result.append({'m2': float(m2), 'P': np.zeros_like(M), 'multiplicity': 0})
        result[-1]['P'] += np.outer(u[:,k], u[:,k].conj())
        result[-1]['multiplicity'] += 1
    for g in result:
        require(rel(g['P'], g['P']@g['P'])<2e-8, 'Heavy group projector not idempotent')
    return result


def deform_heavy_generations(inp: dict, factors: tuple[float, ...]) -> dict:
    """Scale Takagi mass groups while retaining Weyl symmetry and projectors.

    For symmetric M, P_g M = M P_g^T; real scales keep M' symmetric.
    This is not a renormalization-group evolution or Yukawa refit.
    """
    M=t5.decode(inp['fermion_mass_gev'],'M')
    groups=_groups(M)
    require(len(factors)==len(groups) and all(math.isfinite(v) and v>0 for v in factors),
            'Heavy deformation requires one positive factor per eigenmass group')
    new=sum((f*g['P']@M for f,g in zip(factors,groups)),np.zeros_like(M))
    require(rel(new,new.T)<1e-9, 'Spectral heavy deformation broke Weyl mass symmetry')
    return {**inp,'fermion_mass_gev':t5.encode((new+new.T)/2)}


def deform_scalars(inp: dict, d1: int, scenario: str) -> dict:
    A=t5.decode(inp['scalar_A_gev2'],'scalar A')
    require(A.ndim==2 and A.shape[0]==A.shape[1] and 0<d1<len(A),
            'Invalid scalar matrix dimensions')
    new=A.copy()
    if scenario=='scalar_common_shift':
        for i in range(len(A)):
            new[i,i] = A[i,i] * (1.28 if i<d1 else .78)
    elif scenario=='large_scalar_mixing':
        new[:d1,d1:] *= 20.
        new[d1:,:d1] *= 20.
    elif scenario=='charge_split_diagnostic':
        # Formally U(1)_em-preserving mass splittings, not a claim that the
        # given restricted SU(2)-invariant Lagrangian realizes these masses.
        for i in range(len(A)):
            j=i if i<d1 else i-d1
            new[i,i] = A[i,i] * (1. + (.13 if i<d1 else .11) * (j-1.))
    elif scenario=='baseline' or scenario=='heavy_generation_tilt':
        pass
    else:
        raise FactorizationError('Unknown scenario '+scenario)
    require(rel(new,new.conj().T)<1e-12, 'Scalar deformation broke Hermiticity')
    require(all(new[i,i].real>0 for i in range(len(new))), 'Scalar diagonal tachyon')
    return {**inp,'scalar_A_gev2':t5.encode(new)}


def deform(inp: dict, d1: int, scenario: str) -> dict:
    z=deform_scalars(inp,d1,scenario)
    if scenario=='heavy_generation_tilt':
        z=deform_heavy_generations(z,(.71,1.19,1.64))
    return z


def coefficient_groups(inp: dict, tensors: dict, model: str, d1: int) -> dict:
    """Yukawa-independent B_g[r,s] and B_(q,g)[r,s].

    Contract *original* source candidate C1,C2 with M^dagger P_g and exact
    off-diagonal B0(A)_{S2,S1}. No fitted flavour matrix enters this function.
    """
    require(inp.get('model')==model, 'Model label mismatch')
    M=t5.decode(inp['fermion_mass_gev'],'M')
    A=t5.decode(inp['scalar_A_gev2'],'A')
    B=t5.decode(inp['scalar_B_gev2'],'B',A.shape)
    require(norm(B)<1e-14*max(norm(A),1.), 'Restricted scalar B must be zero')
    require(tensors['C1'].shape==(3,len(M),len(A)) and
            tensors['C2'].shape==(3,len(M),len(A)), 'Extracted source tensor dimensions changed')
    pairs=coupled_pairs(A,d1,model[0],int(model[2:]))
    groups=_groups(M)
    out=[]; channels={}
    for g in groups:
        m2=g['m2']; op=M.conj().T@g['P']
        matrix=np.zeros((3,3),complex)
        leading=np.zeros((3,3),complex)
        for pair in pairs:
            q=pair['charge'];s=pair['s1'];u=pair['s2']
            x,y,nonlin=scalar_kernels(m2,pair['a_gev2'],pair['b_gev2'],pair['delta_gev2'])
            # r,j : original generation indices for Y1 and Y2, NOT light flavour.
            coeff=tensors['C1'][:,:,s] @ op @ tensors['C2'][:,:,u].T
            term=x*coeff; matrix+=term;leading+=y*coeff
            channels.setdefault(str(q),[]).append({'m2':float(m2),'B':t5.encode(term),
                                                   'coefficient_norm':norm(term),
                                                   'kernel_exact':value(x),
                                                   'kernel_leading':value(y),
                                                   'relative_nonlinearity':nonlin})
        out.append({'m2':m2, 'multiplicity':g['multiplicity'],
                    'B':matrix, 'B_leading':leading})
    require(len(out)==3, model+': expected 3 distinct heavy mass groups; inspect source')
    return {'groups':out,'charges':channels,'pair_count':len(pairs),
            'mass_gev':[math.sqrt(g['m2']) for g in out],
            'group_multiplicities':[g['multiplicity'] for g in out]}


def assemble(coeff: dict, y1: np.ndarray, y2: np.ndarray, *, leading: bool=False) -> dict:
    require(y1.shape==y2.shape==(3,3) and np.all(np.isfinite(y1)) and np.all(np.isfinite(y2)),
            'Yukawa source must be finite 3x3')
    individual=[]
    for g in coeff['groups']:
        z=np.conj(y1)@(g['B_leading'] if leading else g['B'])@np.conj(y2).T
        individual.append(LOOP*(z+z.T))
    return {'matrix':sum(individual,np.zeros((3,3),complex)), 'groups':individual}


def _fit(a: np.ndarray, b: np.ndarray) -> dict:
    require(a.shape==b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)),
            'Coefficient comparison shape/finite mismatch')
    require(norm(a)>1e-80, 'Reference coefficient tensor is zero')
    c=complex(np.vdot(a,b)/np.vdot(a,a))
    return {'complex_scale':value(c),'absolute_scale':abs(c),
            'matrix_residual':norm(b-c*a)/max(norm(b),EPS)}


def compare_coeff(reference: dict, other: dict) -> dict:
    r=reference['groups'];o=other['groups']
    require(len(r)==len(o)==3, 'Nonmatching heavy-group count')
    for x,y in zip(r,o):
        require(abs(x['m2']-y['m2'])<1e-7*max(x['m2'],y['m2']),
                'Different heavy masses across compared models')
    flat_a=np.concatenate([g['B'].ravel() for g in r])
    flat_b=np.concatenate([g['B'].ravel() for g in o])
    full=_fit(flat_a,flat_b)
    per=[_fit(x['B'],y['B']) for x,y in zip(r,o)]
    c=complex(full['complex_scale']['real'],full['complex_scale']['imag'])
    spread=max(abs(complex(k['complex_scale']['real'],k['complex_scale']['imag'])-c)/max(abs(c),EPS)
               for k in per)
    return {**full,'per_group':per,'max_group_scale_fractional_spread':spread,
            'max_group_shape_residual':max(p['matrix_residual'] for p in per)}


def _csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open('w', newline='', encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def _report(payload: dict) -> str:
    lines=[
        '# Chat A Batch 6 — CG-weighted heavy-kernel factorization', '',
        '**Conditional only:** candidate Matchete/RGBeta→Weyl source dictionary; '
        'no physical normalization certificate and no EFT calculation. '
        'Do not interpret counterfactual masses as fitted or gauge-invariant models.', '',
        '## Source-derived algebra', '',
        r'The extracted candidate Weyl vertices obey '
        r'\(Z_{iIs}=\sum_r Y_{1,ir}^{*}C^{(1)}_{rIs}\), '
        r'\(Y_{jIt}=\sum_sY_{2,js}^{*}C^{(2)}_{sIt}\). '
        'The two-scalar kernel and heavy-mass spectral projector give', '',
        r'\[B^{(R)}_{g,rs}=\sum_{q}K_q(M_g^2)\, '
        r'C^{(1)}_{r,I,s_q}(M^\dagger P_g)_{IJ}C^{(2)}_{s,J,t_q},\]', '',
        r'\[K_q(M^2)=\delta_q\frac{B_0(0;M^2,\lambda_{q+})-'
        r'B_0(0;M^2,\lambda_{q-})}{\lambda_{q+}-\lambda_{q-}},\]', '',
        r'\[M_\nu={1\over16\pi^2}\sum_g\left('
        r'Y_1^* B^{(R)}_gY_2^{*T} + (Y_1^*B^{(R)}_gY_2^{*T})^T\right).\]', '',
        'Here each B is a **3×3 source-generation coefficient matrix** '
        '(not a neutrino mass matrix). B_g is evaluated before choosing any Yukawa matrices. '
        'The heavy group projector P_g includes the *whole degenerate mass eigenspace*.', '',
        '## Exact factorization criterion', '',
        r'If \(B^{(R)}_g=\kappa_R B^{(0)}_g\) with the **same** complex '
        r'\(\kappa_R\) for every heavy group g, then '\
        r'\(M_\nu^{(R)}=\kappa_R M_\nu^{(0)}\) for **arbitrary** Y1,Y2 '
        'at those scalar and heavy mass parameters. No Yukawa sampling is needed to prove '
        'this finite-dimensional statement. A tiny single-spectrum residual does not prove '
        'factorization for all mass parameters.', '',
        'For small scalar mixing, each charge contribution has '
        r'\(K_q=\delta_qD(M_g^2,a_q,b_q)+O(\delta_q|\delta_q|^2)\). '
        'Common unmixed scalar masses and charge-independent mixing magnitudes can make '
        'the common loop function factor out of the CG contractions; this is a '
        '**conditional sufficient scenario**, not an assertion about every possible T3 scalar potential.', '',
        '## Scenarios and interpretation', '',
        '| Scenario | Change | Interpretation |', '|---|---|---|',
        '| baseline | Matched |lambdaT3|=0.10, original scalar and heavy spectra | Source-conditional benchmark comparison |',
        '| scalar_common_shift | S1 diagonal +28%, S2 diagonal −22% | Common SU(2)-multiplet masses shifted |',
        '| heavy_generation_tilt | Three heavy Takagi mass groups ×(0.71,1.19,1.64) | Mass eigenvalues changed with invariant projectors |',
        '| large_scalar_mixing | S1–S2 scalar mixing ×20 | Stress nonlinear exact kernel; not a fitted point |',
        '| charge_split_diagnostic | Multiplet component diagonals deliberately split | Outside original restricted SU(2)-degenerate mass ansatz; **not a derived UV potential** |', '',
        'Every scenario uses the exact two-complex-scalar kernel and leaves vertices/CG tensors unchanged. '
        'Counterfactual scenarios do not alter original outputs.', '',
        '## All-16 Yukawa-independent coefficient comparisons', '',
        'Reference is `A_-2`; each coefficient fit uses **all nine complex entries in all three heavy groups**, '
        'not only the final neutrino matrix. `eps` is the relative Frobenius residual, '
        'and `spread` is the maximum relative disagreement of group-specific multipliers with the common multiplier.', '',
        '| Model | baseline eps | baseline spread | common-shift eps | heavy-tilt eps | large-mix eps | charge-split eps |',
        '|---|---:|---:|---:|---:|---:|---:|']
    lookup={(r['model'],r['scenario']):r for r in payload['coefficient_comparisons']}
    symmetry={(r['model_a'],r['model_b'],r['scenario']):r for r in payload['pairwise_symmetry_tests']}
    for model in payload['models']:
        item=lambda sc:lookup[(model,sc)]
        lines.append(f'| {model} | {item("baseline")["coefficient_residual"]:.3e} | '
                     f'{item("baseline")["max_group_fractional_spread"]:.3e} | '
                     f'{item("scalar_common_shift")["coefficient_residual"]:.3e} | '
                     f'{item("heavy_generation_tilt")["coefficient_residual"]:.3e} | '
                     f'{item("large_scalar_mixing")["coefficient_residual"]:.3e} | '
                     f'{item("charge_split_diagnostic")["coefficient_residual"]:.3e} |')
    lines += ['', '## Heavy-group multipliers at the frozen common-mixing spectrum', '',
              'A single κ across all three heavy groups is the test; κ itself is a '
              '**candidate basis-dependent complex coefficient ratio**, not a physical CG normalization.', '',
              '| Model | Re(κ) | Im(κ) | Re(κ1) | Re(κ2) | Re(κ3) |',
              '|---|---:|---:|---:|---:|---:|']
    for model in payload['models']:
        row=lookup[(model,'baseline')]
        vals=[p['complex_scale']['real'] for p in row['per_group_scales']]
        lines.append(f'| {model} | {row["complex_scale"]["real"]:.9g} | '
                     f'{row["complex_scale"]["imag"]:.9g} | '
                     + ' | '.join(f'{v:.9g}' for v in vals) + ' |')
    lines += ['', '## Mirror/symmetry partner tests at B_g level', '',
              'Unlike the Batch-5 final-matrix fit, this compares all heavy-group '
              'coefficient entries and tests the result under spectral deformations.', '',
              '| Pair | baseline eps | heavy tilt eps | scalar shift eps | charge split eps |',
              '|---|---:|---:|---:|---:|']
    for a,b in payload['matched_pairs']:
        lines.append(f'| {a}/{b} | {symmetry[(a,b,"baseline")]["coefficient_residual"]:.3e} | '
                     f'{symmetry[(a,b,"heavy_generation_tilt")]["coefficient_residual"]:.3e} | '
                     f'{symmetry[(a,b,"scalar_common_shift")]["coefficient_residual"]:.3e} | '
                     f'{symmetry[(a,b,"charge_split_diagnostic")]["coefficient_residual"]:.3e} |')
    lines += ['', '## Exact versus small-mixing group tensors', '',
              'The file `batch06_linear_approximation.csv` records the relative '
              'Frobenius error of the first-order scalar-mixing approximation '
              'for each model and stress scenario. The exact kernel is used for '
              'all other results.', '',
              '## Same-texture full neutrino matrix checks', '',
              'For each coefficient set, two reference Yukawa pairs are inserted into the '\
              'algebraic expression for Mnu. The real frozen source was also reconstructed '
              'and checked against Batch 1/3/5. Matrix alignment across representations is '
              'a consequence of B_g proportionality, not an independent UV measurement.', '',
              '## E charge-channel equality test', '',
              'Baseline E has two constructively aligned charge sectors (Batches 2/3). '
              'The deliberate component diagonal split changes their scalar kernel functions; '
              'the reported per-charge coefficient norms show if that equality survives. '
              'These are **coefficient** norms before choosing a Yukawa texture.', '',
              '| E model | scenario | charge 1 norm | charge 2 norm | norm ratio |',
              '|---|---|---:|---:|---:|']
    for row in payload['E_channel_diagnostics']:
        lines.append(f'| {row["model"]} | {row["scenario"]} | '
                     f'{row["charge_1_norm"]:.6e} | {row["charge_2_norm"]:.6e} | '
                     f'{row["norm_ratio"]:.8g} |')
    lines += ['', '## Interpretation guardrails', '',
              '- A small epsilon for one frozen spectrum proves only a finite-dimensional '
              'identity of the saved candidate Weyl tensors at that spectrum.',
              '- If epsilon is small for **all** mass deformations, that is stronger evidence '
              'for structural factorization, not a proof for an arbitrary potential.',
              '- If a deformation raises epsilon, it provides an explicit counterexample '
              'to extrapolating the baseline identity beyond its assumptions.',
              '- Different models retain different fermion representations and '
              'Dirac/Majorana constructions. No shared physical coupling convention is implied.',
              '- The charge-split diagnostic is deliberately outside the restricted '
              'source multiplet-degenerate scalar mass ansatz; its gauge-invariant UV '
              'completion has **not** been established.',
              '- Numerical precision at ~1e−14 is not a physical amplitude certificate.',
              '- The 0/16 B1/B2 normalization certificate count is unchanged.', '']
    return '\n'.join(lines)


def analyse(root: Path, *, output_dir: Path | None = None,
            batch05_dir: Path | None = None, make_plots: bool = True) -> dict:
    """Read-only verified all-16 sources; emit a NEW non-overwriting report directory."""
    from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary
    from Numerical.lo_patterns_16.result_analysis import _source_hash_guard

    root=root.resolve()
    owner=(root/BASE).resolve()
    b5_dir=(batch05_dir or root/BATCH05).resolve()
    target=(output_dir or root/TARGET).resolve()
    require(all(x.is_relative_to(owner) and x != owner for x in (b5_dir,target)) and
            b5_dir!=target,'Sources/output must be distinct and contained within Chat A results')
    require(not target.exists() or (target.is_dir() and not any(target.iterdir())),
            'Output directory already populated: '+str(target))

    study=read(root/STUDY/'all16_summary.json')
    b5=read(b5_dir/'batch05_synthesis.json')
    manifest=read(b5_dir/'SOURCE_MANIFEST.json')
    require(study.get('format')=='t3-conditional-pattern-study-all16-v1' and
            b5.get('format')=='t3-chatA-batch05-analytical-yukawa-texture-v1' and
            manifest.get('format')=='t3-chatA-batch05-source-manifest-v1',
            'Source schema/format drift')
    require(study.get('physical_normalisation_certified') is False and
            b5.get('physical_normalisation_certified') is False and
            manifest.get('physical_normalisation_certified') is False and
            b5.get('safe_for_physical_source_certification') is False,
            'Unexpected physical certification claim')
    require(b5.get('model_count')==16 and study['counts']['conditional_evaluated']==16 and
            study['counts']['blocked']==0 and study['counts']['physical_source_certified']==0,
            'Incomplete source run')
    require(manifest['study_sha256']==digest(root/STUDY/'all16_summary.json') and
            manifest['batch03_sha256']==digest(root/BATCH03/'batch03_analytic_explanation.json') and
            manifest['batch04_sha256']==digest(root/BATCH04/'batch04_flavour_analysis.json'),
            'Batch5 study/Batch3/Batch4 hash chain drift')
    catalogue=model_dictionary();models=[r['model'] for r in catalogue]
    require([r['model'] for r in b5['models']]==models and
            [r['model'] for r in study['models']]==models and
            sorted(manifest['models'])==sorted(models), 'Source model catalogue mismatch')
    meta={r['model']:r for r in catalogue}
    sources={};textures={};source_hashes={}
    for model in models:
        base=read(root/STUDY/'models'/f'{model}.json')
        src=base['source'];p_in=t5._resolve_local(root,src['input'])
        p_b=t5._resolve_local(root,src['benchmark'])
        require(base.get('model')==model and base.get('status')=='conditional_uv_evaluated' and
                base.get('physical_normalisation_certified') is False,
                model+': source model status drift')
        require(digest(p_in)==src['input_sha256']==manifest['models'][model]['candidate_weyl_input_sha256'] and
                digest(p_b)==src['benchmark_sha256']==manifest['models'][model]['benchmark_sha256'] and
                digest(root/STUDY/'models'/f'{model}.json')==manifest['models'][model]['study_model_sha256'],
                model+': source file SHA mismatch')
        _source_hash_guard(root,src,study,model)
        inp=read(p_in);frozen=read(p_b)
        tensors=t5.extract_source_tensors(inp,frozen,model)
        d1=meta[model]['SU2_dimensions']['S1']
        lam=frozen['base_state']['ordinary']['lambdaT3']
        lam=abs(complex(lam['real'],lam['imag'])) if isinstance(lam,dict) else abs(complex(lam))
        require(lam>0,'Frozen lambdaT3 cannot be zero')
        controlled=t5.scalar_mixing_scale(inp,d1,.10/lam)
        sources[model]={'input':controlled,'d1':d1,'tensors':tensors,'frozen':frozen}
        textures[model]=(tensors['y1'],tensors['y2'])
        source_hashes[model]={'input_sha256':src['input_sha256'],
                              'benchmark_sha256':src['benchmark_sha256'],
                              'study_model_sha256':manifest['models'][model]['study_model_sha256']}
    comparison_rows=[];details=[];charge_rows=[];masscheck_rows=[];model_summary=[]
    raw_coefficients={};linear_rows=[];pairwise=[]
    matched=[('A_-4','D_+2'),('A_-2','D_+0'),('A_+0','D_-2'),
             ('B_-3','B_+1'),('C_-3','C_+1'),('E_-4','E_+2'),('E_-2','E_+0')]
    for scenario in SCENARIOS:
        scenario_data={}
        for model in models:
            info=sources[model]
            inp=deform(info['input'],info['d1'],scenario)
            coeff=coefficient_groups(inp,info['tensors'],model,info['d1'])
            scenario_data[model]=coeff
            raw_coefficients.setdefault(scenario,{})[model]={
                'heavy_groups':[{'mass_gev':math.sqrt(g['m2']),
                                 'multiplicity':g['multiplicity'],
                                 'B':t5.encode(g['B']), 'B_leading':t5.encode(g['B_leading'])}
                                for g in coeff['groups']],
                'charge_channel_kernels':coeff['charges']}
            group_exact=np.concatenate([g['B'].ravel() for g in coeff['groups']])
            group_linear=np.concatenate([g['B_leading'].ravel() for g in coeff['groups']])
            linear_rows.append({'scenario':scenario,'model':model,
                                'B_group_linear_relative_error':rel(group_linear,group_exact)})
            # Directly test the independent B_g algebra against the already
            # checked 2x2 scalar + fermion matrix-function reconstruction.
            if scenario=='baseline':
                y1,y2=textures[model]
            else:
                y1,y2=textures['C_+1']
            direct=reconstruct(t5.substitute_textures(inp,info['tensors'],y1,y2),
                               model[0],int(model[2:]),info['d1'])
            bytensor=assemble(coeff,y1,y2)
            err=rel(direct['exact'],bytensor['matrix'])
            require(err<5e-8,model+' '+scenario+f': B_g tensor/Mnu reconstruction failed ({err:g})')
            for a,b in zip(direct['heavy_groups'],bytensor['groups']):
                require(rel(a['matrix'],b)<5e-8,
                        model+' '+scenario+': heavy-eigenspace individual matrix mismatch')
            # Independently recheck original mass matrix with fitted Yukawas
            # and original scalar lambda (only the special case E differs).
            if scenario=='baseline':
                if model[0]!='E':
                    expected=t5.decode(read(root/STUDY/'models'/f'{model}.json')['mnu_gev'],'saved Mnu')
                    require(rel(bytensor['matrix'],expected)<2e-7,
                            model+': original Batch1 matrix changed')
                b5data=next(z for z in b5['models'] if z['model']==model)
                for ref in TEXTURES:
                    yy1,yy2=textures[ref]
                    mass=np.sort(np.linalg.svd(assemble(coeff,yy1,yy2)['matrix']*1e9,compute_uv=False))
                    saved=b5data['common_mixing_controlled'][ref]
                    require(np.max(np.abs(mass-np.array([saved['m1_ev'],saved['m2_ev'],saved['m3_ev']])))<
                            2e-5*max(mass[-1],1e-20),
                            model+' '+ref+': Batch5 common texture spectra do not match')
            masscheck_rows.append({'scenario':scenario,'model':model,
                                   'exact_tensor_reconstruction_error':err,
                                   'heavy_groups':len(coeff['groups']),
                                   'scalar_pairs':coeff['pair_count']})
            if model[0]=='E' and scenario in ('baseline','charge_split_diagnostic','large_scalar_mixing'):
                entries=[]
                for q,ks in coeff['charges'].items():
                    entries.append((q,math.sqrt(sum(z['coefficient_norm']**2 for z in ks))))
                entries=[x for x in entries if x[1]>1e-14*max(z[1] for z in entries)]
                require(len(entries)==2,model+': expected two E source coefficient charge channels')
                charge_rows.append({'model':model,'scenario':scenario,
                                    'charge_1':entries[0][0],'charge_2':entries[1][0],
                                    'charge_1_norm':entries[0][1],
                                    'charge_2_norm':entries[1][1],
                                    'norm_ratio':entries[0][1]/max(entries[1][1],EPS)})
        # Check full coefficient array before any optional flavour insertion.
        reference=scenario_data[REFERENCE_MODEL]
        for left,right in matched:
            pairfit=compare_coeff(scenario_data[left],scenario_data[right])
            pairwise.append({'scenario':scenario,'model_a':left,'model_b':right,
                             'coefficient_residual':pairfit['matrix_residual'],
                             'max_heavy_multiplier_spread':pairfit['max_group_scale_fractional_spread'],
                             'complex_scale':pairfit['complex_scale']})
        for model in models:
            other=scenario_data[model]
            fit=compare_coeff(reference,other)
            comparison_rows.append({'scenario':scenario,'model':model,
                                    'reference_model':REFERENCE_MODEL,
                                    'coefficient_residual':fit['matrix_residual'],
                                    'max_group_fractional_spread':fit['max_group_scale_fractional_spread'],
                                    'max_single_group_shape_residual':fit['max_group_shape_residual'],
                                    'complex_scale':fit['complex_scale'],
                                    'absolute_scale':fit['absolute_scale'],
                                    'per_group_scales':fit['per_group'],
                                    'fermion_mass_groups_gev':other['mass_gev'],
                                    'multiplicities':other['group_multiplicities']})
            for refname in ('C_+1','C_-1'):
                yy1,yy2=textures[refname]
                vref=assemble(reference,yy1,yy2)['matrix']
                vmod=assemble(other,yy1,yy2)['matrix']
                mfit=_fit(vref,vmod)
                details.append({'scenario':scenario,'model':model,'texture_reference':refname,
                                'neutrino_matrix_proportional_residual':mfit['matrix_residual'],
                                'matrix_multiplier':mfit['complex_scale'],
                                'ratio_m1_over_m3':float(np.linalg.svd(vmod,compute_uv=False)[-1]/
                                                         max(np.linalg.svd(vmod,compute_uv=False)[0],EPS))})
    require(len(comparison_rows)==16*len(SCENARIOS) and
            len(masscheck_rows)==16*len(SCENARIOS), 'All-16 scenario coverage incomplete')
    by={(p['model'],p['scenario']):p for p in comparison_rows}
    for model in models:
        baseline=by[(model,'baseline')]
        model_summary.append({'model':model,
                              'baseline_coefficient_residual':baseline['coefficient_residual'],
                              'baseline_heavy_scale_spread':baseline['max_group_fractional_spread'],
                              'largest_stress_coefficient_residual':max(by[(model,s)]['coefficient_residual']
                                                                     for s in SCENARIOS if s!='baseline'),
                              'baseline_algebraic_factorization_at_1e-9':baseline['coefficient_residual']<1e-9 and
                                baseline['max_group_fractional_spread']<1e-9})
    payload={'format':'t3-chatA-batch06-heavy-kernel-factorization-v1',
             'source_conditional_only':True,'physical_normalisation_certified':False,
             'safe_for_physical_source_certification':False,
             'source_revision':study['provenance']['revision'],
             'models':models,'model_count':16,
             'reference_model':REFERENCE_MODEL,
             'scenarios':list(SCENARIOS),
             'source_sha256':{'study':digest(root/STUDY/'all16_summary.json'),
                              'batch05':digest(b5_dir/'batch05_synthesis.json'),
                              'batch05_manifest':digest(b5_dir/'SOURCE_MANIFEST.json')},
             'coefficient_comparisons':comparison_rows,
             'flavour_matrix_controls':details,
             'E_channel_diagnostics':charge_rows,
             'tensor_reconstruction_checks':masscheck_rows,
             'linear_approximation':linear_rows,
             'matched_pairs':matched,
             'pairwise_symmetry_tests':pairwise,
             'kernel_coefficients':raw_coefficients,
             'model_summary':model_summary}
    markdown=_report(payload)
    # No write happens until ALL input / algebraic checks have succeeded.
    target.mkdir(parents=True,exist_ok=True)
    (target/'BATCH06_ANALYTICAL_FACTORISATION.md').write_text(markdown,encoding='utf-8')
    (target/'batch06_factorization.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    _csv(target/'batch06_coefficient_tests.csv',
         [{k:v for k,v in row.items() if k not in ('complex_scale','per_group_scales',
                                                  'fermion_mass_groups_gev','multiplicities')}
          for row in comparison_rows],
         ['scenario','model','reference_model','coefficient_residual',
          'max_group_fractional_spread','max_single_group_shape_residual','absolute_scale'])
    _csv(target/'batch06_symmetry_pairs.csv',
         [{k:v for k,v in x.items() if k!='complex_scale'} for x in pairwise],
         ['scenario','model_a','model_b','coefficient_residual','max_heavy_multiplier_spread'])
    _csv(target/'batch06_linear_approximation.csv',linear_rows,
         ['scenario','model','B_group_linear_relative_error'])
    _csv(target/'batch06_flavour_matrix_tests.csv',
         [{k:v for k,v in row.items() if k!='matrix_multiplier'} for row in details],
         ['scenario','model','texture_reference','neutrino_matrix_proportional_residual','ratio_m1_over_m3'])
    _csv(target/'batch06_E_charge_channels.csv',charge_rows,
         ['model','scenario','charge_1','charge_2','charge_1_norm','charge_2_norm','norm_ratio'])
    _csv(target/'batch06_tensor_reconstruction.csv',masscheck_rows,
         ['scenario','model','exact_tensor_reconstruction_error','heavy_groups','scalar_pairs'])
    provenance={'format':'t3-chatA-batch06-manifest-v1',
                'physical_normalisation_certified':False,
                'source_sha256':payload['source_sha256'],'models':source_hashes,
                'tests':{'source_inputs_checked':16,'coefficient_scenarios':80,
                         'flavour_matrix_controls':len(details),
                         'physical_source_certified':0}}
    (target/'SOURCE_MANIFEST.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    handoff={'format':'t3-chatA-batch06-handoff-v1','all16_evaluated':16,
             'scenarios':len(SCENARIOS),'physical_normalisation_certified':False,
             'source_manifest_sha256':digest(target/'SOURCE_MANIFEST.json'),
             'primary_report':'BATCH06_ANALYTICAL_FACTORISATION.md',
             'open_certification':'B1/B2 source physical Weyl normalization/phase'}
    (target/'HANDOFF.json').write_text(json.dumps(handoff,indent=2)+'\n',encoding='utf-8')
    if make_plots:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            pass
        else:
            x=np.arange(len(models));fig,ax=plt.subplots(figsize=(13,5))
            for sc in SCENARIOS:
                ax.plot(x,[max(by[(m,sc)]['coefficient_residual'],1e-17) for m in models],
                        '.-',label=sc)
            ax.set_yscale('log');ax.set_xticks(x,models,rotation=55)
            ax.set_ylabel('Yukawa-independent 3-group coefficient residual')
            ax.set_title('Kernel factorisation across original and altered mass spectra (conditional)')
            ax.legend(ncol=2,fontsize=8);fig.tight_layout()
            fig.savefig(target/'batch06_factorization_stress.png',dpi=170);plt.close(fig)
    return payload
