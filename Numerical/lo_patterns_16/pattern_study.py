"""Isolated, fail-closed all-16 restricted T3 *conditional* UV pattern study.

Imports existing production numerical kernels without changing them. This module
never certifies a native Matchete/RGBeta->Weyl dictionary. No EFT imports.
"""
from __future__ import annotations

import ast
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.direct_one_loop.core import _complex, encoded, evaluate
from Numerical.lo_certification.direct_uv.native_uv_loop_reference import (
    independent_reference, decode,
)
from Numerical.lo_certification.direct_uv.native_amplitude_contract import (
    model_dictionary, source_audit,
)
from Numerical.lo_certification.direct_uv.certification import (
    charge_bases, validate_components, all_models,
)
from Numerical.lo_certification.direct_uv.channels import partition_exact

SOURCE_AUDITS = (
    'component_support_audit.json', 'charge_resolved_loop_weights.json',
    'cg_covariance_audit.json', 'matchete_generator_analysis.json',
    'charge_channels.json',
)
INPUT_REL = Path('output/general_observable/lo_direct_one_loop/provisional_inputs')
AUDIT_REL = Path('output/general_observable/lo_direct_one_loop/source_audits')
OUTPUT_REL = Path('output/general_observable/lo_pattern_study16/provisional_01')
BASELINE_REV = '0b210f18f6a2cfb10fff03059e2ec4a73d2195a4'
T_VALUES = (-1., 0., .5, 1., 2.)
HEAVY_VALUES = (.5, 1., 2.)


class StudyError(ValueError):
    pass


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def jload(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8-sig'))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def source_tree_state(root: Path) -> dict:
    """Record actual local checkout; do not equate it with the remote baseline."""
    import subprocess
    try:
        cmd = subprocess.run(['git','-C',str(root),'rev-parse','HEAD'],
                             capture_output=True, text=True, check=True, timeout=8)
        head = cmd.stdout.strip()
        state = subprocess.run(['git','-C',str(root),'status','--porcelain'],
                               capture_output=True, text=True, check=True, timeout=8)
        dirty = bool(state.stdout.strip())
    except (OSError, subprocess.SubprocessError):
        head, dirty = None, None
    return {'local_head': head, 'local_worktree_dirty': dirty,
            'reviewed_remote_main_head_at_implementation': BASELINE_REV,
            'local_head_matches_reviewed_remote': head == BASELINE_REV if head else None}


def audit_weight_lookup(payload: dict, catalogue: dict[str, dict]) -> dict:
    if payload.get('status') != 'charge_resolved_group_coefficients' or payload.get('count') != 16:
        raise StudyError('Source CG weight audit status/count changed')
    if payload.get('contraction') != 'Y1 * Y2 * conjugate(T3MixCG) using same physical internal weight indices':
        raise StudyError('Source CG weight contraction or orientation changed')
    out = {}
    for item in payload['models']:
        label = f"{item['class']}_{int(item['model']['alpha']):+d}"
        if label not in catalogue or label in out:
            raise StudyError('Unexpected / duplicate CG audit model '+label)
        dim = catalogue[label]['SU2_dimensions']
        if (item['model']['dS1'], item['model']['dS2'], item['model']['dF']) != (dim['S1'],dim['S2'],dim['F']):
            raise StudyError('CG audit representation mismatch '+label)
        for ch in item['channels']:
            for name, d in (('S1',dim['S1']),('S2',dim['S2']),('F',dim['F'])):
                n = int(ch['indices_weight_basis_one_based'][name])
                if not 1 <= n <= d:
                    raise StudyError('Bad CG component index '+label)
            for name in ('y1','y2','mixing','complex_weight','charge'):
                if name not in ch:
                    raise StudyError('Missing source CG weight data '+label)
            # Source audit's exact algebraic convention, not an arbitrary weight
            # assigned to a charge: Y1*Y2*conjugate(T3MixCG).
            product = (_simple_expr(ch['y1'])*_simple_expr(ch['y2'])*
                       _simple_expr(ch['mixing']).conjugate())
            audited = _simple_expr(ch['complex_weight'])
            if abs(product-audited)>1e-10*max(abs(product),abs(audited),1e-100):
                raise StudyError('Source CG product/complex phase mismatch '+label)
            base = charge_bases(label[0],int(item['model']['alpha']))
            charge = float(_simple_expr(ch['charge']).real)
            i1=int(ch['indices_weight_basis_one_based']['S1'])-1
            i2=int(ch['indices_weight_basis_one_based']['S2'])-1
            fi=int(ch['indices_weight_basis_one_based']['F'])-1
            q1=base['scalar'][i1];q2=base['scalar'][dim['S1']+i2]
            qf=base['source_fermion_weights'][fi]
            if not (np.isclose(q1,charge) and np.isclose(q2,charge) and np.isclose(qf,charge)):
                raise StudyError('CG indexed component charge mismatch '+label)
        if len(item['channels']) < 1:
            raise StudyError('No charge weights in source audit '+label)
        out[label] = item
    if set(out) != set(catalogue):
        raise StudyError('CG audit does not cover exact 16-model catalogue')
    return out


def _simple_expr(text: str) -> complex:
    """Restricted expression evaluator for *audit numbers*, never arbitrary eval()."""
    names = {'I':1j, 'pi':math.pi}
    fun = {'sqrt':lambda x:x**.5}
    def walk(node):
        if isinstance(node,ast.Constant) and type(node.value) in (int,float,complex): return complex(node.value)
        if isinstance(node,ast.Name) and node.id in names: return names[node.id]
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,(ast.UAdd,ast.USub)):
            v=walk(node.operand);return v if isinstance(node.op,ast.UAdd) else -v
        if isinstance(node,ast.BinOp) and isinstance(node.op,(ast.Add,ast.Sub,ast.Mult,ast.Div,ast.Pow)):
            a,b=walk(node.left),walk(node.right)
            if isinstance(node.op,ast.Add):return a+b
            if isinstance(node.op,ast.Sub):return a-b
            if isinstance(node.op,ast.Mult):return a*b
            if isinstance(node.op,ast.Div):return a/b
            if abs(b.imag)>1e-14 or abs(b.real)>8:raise StudyError('Unsafe exponent')
            return a**b
        if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id in fun and len(node.args)==1 and not node.keywords:
            return fun[node.func.id](walk(node.args[0]))
        raise StudyError('Unsupported source CG expression: '+text)
    ans = walk(ast.parse(text, mode='eval').body)
    if not (math.isfinite(ans.real) and math.isfinite(ans.imag)):
        raise StudyError('Nonfinite source CG coefficient')
    return ans


def _candidate_group_channels(item:dict) -> list[dict]:
    return [{**ch, 'charge_numeric':float(_simple_expr(ch['charge']).real),
             'complex_weight_numeric':{'real':_simple_expr(ch['complex_weight']).real,
                                       'imag':_simple_expr(ch['complex_weight']).imag},
             'status':'native_CG_candidate_not_physical_loop_normalisation'}
            for ch in item['channels']]


def _norm(a) -> float:
    return float(np.linalg.norm(a))


def _mass(ans:dict) -> np.ndarray:
    return decode(ans['total']['complex_symmetric_matrix_gev'], 'mass', (3,3))


def _copy_encoded_array(raw:Any, factor:complex=1) -> dict:
    return encoded(_complex(raw,'tensor')*factor)


def transformed_input(data:dict, mode:str, t:float, d_s1:int) -> dict:
    """Numerical intervention on frozen Weyl tensors, NOT a new source-derived model."""
    if not math.isfinite(t):raise StudyError('Nonfinite scan parameter')
    o = dict(data)
    o['vertices'] = [dict(v) for v in data['vertices']]
    o['provenance'] = dict(data.get('provenance') or {})
    o['provenance']['scan_intervention']={'mode':mode,'parameter':t,
       'conditional_on_provisional_mapping':True}
    if mode in ('y1','y2'):
        name = 'y1_HC' if mode == 'y1' else 'y2_HC'
        matches = [v for v in o['vertices'] if v['name']==name]
        if len(matches)!=1:raise StudyError('Cannot find unique scan Yukawa vertex '+name)
        for field in ('field','conjugate'):
            if field in matches[0]:matches[0][field] = _copy_encoded_array(matches[0][field], t)
    elif mode == 'lambdaT3':
        A = _complex(data['scalar_A_gev2'],'A').copy()
        # Original restricted input builder isolates lambdaT3 in cross-S1/S2 A.
        A[:d_s1,d_s1:]*=t
        A[d_s1:,:d_s1]*=t
        o['scalar_A_gev2'] = encoded(A)
    elif mode == 'heavy_mass':
        if t <= 0:raise StudyError('Nonpositive heavy-mass scan factor')
        o['fermion_mass_gev'] = _copy_encoded_array(data['fermion_mass_gev'], t)
    else:
        raise StudyError('Unknown scan mode '+mode)
    return o


def _compare(ans:dict, ref:dict) -> dict:
    m=_mass(ans)
    q=decode(ref['complex_symmetric_matrix_gev'],'reference mass',(3,3))
    diff=_norm(m-q)
    scale=max(_norm(q),_norm(m),1e-100)
    val=diff/scale
    return {'absolute_matrix_error_gev':diff,'relative_matrix_error':val,
            'absolute_over_reference_sum_component_norms':diff/max(
                float(ref['diagnostics']['sum_of_absolute_component_norms_gev']),1e-100),
            'reference_cancellation_conditioning':ref['diagnostics']['cancellation_factor_absolute_over_total']}


def _scans(inp:dict, baseline:dict, dim:int, model:str) -> list[dict]:
    rows=[]
    m0=_mass(baseline)
    norm0=_norm(m0)
    for mode, ts in (('y1',T_VALUES),('y2',T_VALUES),('lambdaT3',T_VALUES),
                     ('heavy_mass',HEAVY_VALUES)):
        for t in ts:
            row={'model':model,'mode':mode,'t':t,'status':'blocked',
                 'm1_ev':None,'m2_ev':None,'m3_ev':None,
                 'matrix_scaling_residual':None,'reason':None,
                 'prediction':'exact linear for y1/y2 if no same-vertex terms' if mode in ('y1','y2')
                     else 'leading O(lambda) for small mixing; higher odd powers' if mode=='lambdaT3'
                     else 'no universal power law'}
            try:
                uv=transformed_input(inp,mode,t,dim)
                ans=evaluate(uv,retain_individual_terms=False)
                if not ans['diagnostics']['uv_finite']:
                    raise StudyError('UV pole failure after scan variation')
                ms=ans['total']['neutrino_masses_ev']
                if len(ms)!=3 or not np.isfinite(ms).all():raise StudyError('Invalid scan eigenmasses')
                m=_mass(ans)
                row.update(status='conditional_uv_scanned',m1_ev=ms[0],m2_ev=ms[1],m3_ev=ms[2],
                           mass_norm_gev=_norm(m), cancellation_factor=ans['diagnostics']['cancellation_factor_sum_over_total'])
                if mode in ('y1','y2','lambdaT3'):
                    row['matrix_scaling_residual']=_norm(m-t*m0)/max(norm0,1e-100)
            except (ValueError,KeyError,RuntimeError,TypeError,ArithmeticError,StudyError) as exc:
                row['reason']=type(exc).__name__+': '+str(exc)
            rows.append(row)
    return rows


def _charge_rows(inp:dict, dim:int, item:dict, channels:dict, kind:str,alpha:int) -> dict:
    """Return *source-component* support plus exact charge-partition matrices.

    Avoid treating Takagi eigenstate indices as definite electric-charge labels.
    """
    A=_complex(inp['scalar_A_gev2'],'A')
    q=charge_bases(kind,alpha)['scalar']
    y1=next(v for v in inp['vertices'] if v['name']=='y1_HC')
    y2=next(v for v in inp['vertices'] if v['name']=='y2_HC')
    z=_complex(y1['conjugate'],'y1 conjugate')
    y=_complex(y2['field'],'y2 field')
    audit=_candidate_group_channels(item)
    charged_support={float(c['charge_numeric']) for c in audit if c['nonzero']}
    records=[]
    # Every equal-charge pair that could mix across scalar multiplets is enumerated,
    # including null group channels. Complex phases of Weyl vertices are retained.
    for s in range(dim):
        for t in range(dim,A.shape[0]):
            if not np.isclose(q[s],q[t]):continue
            charge=float(q[s]);mix=A[t,s]
            zslice=z[:,:,s];yslice=y[:,:,t]
            record={'charge':charge,'S1_zero_based':s,'S2_zero_based':t-dim,
                    'candidate_mix_A_S2S1_gev2':{'real':mix.real,'imag':mix.imag},
                    'y1_tensor_frobenius':_norm(zslice),'y2_tensor_frobenius':_norm(yslice),
                    'simultaneous_component_support':bool(abs(mix)>1e-16 and _norm(zslice)>1e-16 and _norm(yslice)>1e-16),
                    'in_native_group_audit':charge in charged_support,
                    'group_coefficients':[c['complex_weight'] for c in audit if c['charge_numeric']==charge]}
            records.append(record)
    predicted={str(c['charge_numeric']) for c in audit if c['nonzero']}
    supported={str(r['charge']) for r in records if r['simultaneous_component_support']}
    if predicted!=supported:
        raise StudyError('Candidate CG topology charge support differs from source Weyl tensor/mixing intersection: '
                         +str({'cg':sorted(predicted),'tensor':sorted(supported)}))
    partition=channels['charge_channel_matrices']
    for c in audit:
        lab=f"{c['charge_numeric']:g}|{c['charge_numeric']:g}"
        if lab not in partition:raise StudyError('Charge partition lacks CG topology charge '+lab)
    return {'source_scalar_equal_charge_cross_pairs_including_zeros':records,
            'native_group_audit_channels':audit,
            'group_weight_total':item['total_group_coefficient'],
            'charge_partition_direct_exact':channels,
            'prediction_checked_against_input_CG_support':True,
            'warning':'Source charge labels are scalar-component basis labels. The exact charge partition includes all vertex interference and is not a product of CG weights times one common loop function.'}


def _reconstruct_individual(ans:dict) -> float:
    terms=ans.get('individual_terms')
    if terms is None:raise StudyError('Individual terms were not retained')
    m=sum((_complex(t['matrix_gev'],'individual',(3,3)) for t in terms),np.zeros((3,3),complex))
    return _norm(m-_mass(ans))/max(_norm(_mass(ans)),
            1e-12*float(ans['diagnostics']['sum_individual_frobenius_gev']),1e-90)


def _validate_dictionary(path:Path | None, root:Path) -> dict:
    """A user-supplied document is evidence, not self-validating authority.

    There is no agreed B1/B2 certificate schema at this revision. Reject any
    purported promotion, rather than accepting a JSON boolean as proof.
    """
    if path is None:return {'status':'not_supplied','physical_source_certified':False}
    if not path.is_file():raise StudyError('Requested physical dictionary does not exist: '+str(path))
    obj=jload(path)
    if obj.get('physical_source_certified') is True or obj.get('physical_normalisation_certified') is True:
        raise StudyError('Physical certification cannot be promoted by unverified JSON boolean; coordinator witness verification is not implemented')
    raise StudyError('B1/B2 coordinator certificate schema and independent validation not yet available: file is NOT a certified input')


def _plots(out:Path, rows:list[dict], scans:list[dict]):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return ['matplotlib unavailable: plots not produced']
    active=[x for x in rows if x['status']=='conditional_uv_evaluated']
    if not active:return ['No benchmark outputs: plots not produced']
    labels=[x['model'] for x in active]
    values=[max(x['m3_ev'],1e-70) for x in active]
    fig,ax=plt.subplots(figsize=(11,4))
    ax.bar(range(len(labels)),values)
    ax.set_yscale('log');ax.set_xticks(range(len(labels)),labels,rotation=55,ha='right')
    ax.set_ylabel('Largest Takagi mass [eV], provisional')
    ax.set_title('Frozen restricted T3 benchmarks — conditional normalisation')
    fig.tight_layout();fig.savefig(out/'frozen_m3_conditional.png',dpi=150);plt.close(fig)
    dat=[r for r in scans if r['mode']=='lambdaT3' and r['status']=='conditional_uv_scanned' and r['t']>=0]
    if dat:
        fig,ax=plt.subplots(figsize=(8,5))
        for label in labels:
            sub=sorted((r for r in dat if r['model']==label),key=lambda r:r['t'])
            if sub:
                ax.plot([r['t'] for r in sub],[r['m3_ev'] for r in sub],marker='.',linewidth=.8,label=label)
        ax.set_xlabel('lambdaT3 multiplication t');ax.set_ylabel('Largest mass [eV]')
        ax.set_title('Controlled scalar-mixing intervention (conditional)')
        ax.legend(fontsize=7,ncol=4);fig.tight_layout()
        fig.savefig(out/'lambda_scan_conditional.png',dpi=150);plt.close(fig)
    return []


def run_study(root:Path, output:Path|None=None, input_dir:Path|None=None,
              physical_dictionary:Path|None=None, *, do_scans:bool=True,
              do_plots:bool=True) -> dict:
    root=root.resolve()
    output=(output or root/OUTPUT_REL).resolve()
    owned=(root/'output/general_observable/lo_pattern_study16').resolve()
    if not output.is_relative_to(owned) or output==owned:
        raise StudyError('Output must be a NEW directory beneath '+str(owned))
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Output already exists and is nonempty; select a fresh --output-dir')
    input_dir=(input_dir or root/INPUT_REL).resolve()
    cert=_validate_dictionary(physical_dictionary,root)
    models=model_dictionary();lookup={m['model']:m for m in models}
    if len(lookup)!=16:raise StudyError('Source catalogue not 16 distinct model names')
    if [r[0] for r in all_models()]!=[m['model'] for m in models]:
        raise StudyError('Source model inventories disagree')
    guard=source_audit(root)
    # All results fail closed on source drift. A report still contains 16 blocked rows.
    audit_path=root/AUDIT_REL/'charge_resolved_loop_weights.json'
    audit = jload(audit_path) if audit_path.is_file() else None
    audit_error=None
    try:
        if guard['status']!='source_forms_verified':raise StudyError('Source forms drift: '+'; '.join(guard['errors']))
        if not audit:raise StudyError('Missing charge-resolved CG audit')
        weight_map=audit_weight_lookup(audit,lookup)
    except (StudyError,KeyError,TypeError,ValueError) as err:
        audit_error=str(err);weight_map={}
    audit_hashes={p:(digest(root/AUDIT_REL/p) if (root/AUDIT_REL/p).is_file() else None)
                  for p in SOURCE_AUDITS}
    # Do not claim the remote revision equals the user's actual Windows checkout.
    provenance={'revision':source_tree_state(root), 'source_audit':guard,
                'source_audit_sha256':audit_hashes,
                'source_audits_verified_for_this_input':audit_error is None,
                'input_dir':str(input_dir),'physical_dictionary':cert}
    rows=[];scan_rows=[];weight_rows={}
    for short,key in all_models():
        bfile=root/'Numerical/benchmarks/fitted_reference'/f'{key}.json'
        pfile=input_dir/f'{short}.json'
        row={'model':short,'family':short[0],'alpha':int(short[2:]),
             'dimensions':lookup[short]['SU2_dimensions'],
             'dS1':lookup[short]['SU2_dimensions']['S1'],
             'dS2':lookup[short]['SU2_dimensions']['S2'],
             'dF':lookup[short]['SU2_dimensions']['F'],
             'candidate_topology_charge_channel_count':None,
             'fermion_type':lookup[short]['fermion_type'],
             'status':'blocked','reason':None,'m1_ev':None,'m2_ev':None,'m3_ev':None,
             'benchmark_sha256':digest(bfile) if bfile.is_file() else None,
             'component_input_sha256':digest(pfile) if pfile.is_file() else None,
             'physical_normalisation_certified':False}
        try:
            if audit_error:raise StudyError(audit_error)
            if not bfile.is_file():raise StudyError('Missing fitted frozen benchmark '+str(bfile))
            if not pfile.is_file():raise StudyError('Missing provisional component input '+str(pfile))
            b=jload(bfile);d=jload(pfile)
            if d.get('model')!=short:raise StudyError('Input model label mismatch')
            if d.get('provenance',{}).get('benchmark_sha256')!=digest(bfile):
                raise StudyError('Input not generated from these benchmark bytes')
            pmeta=d['provenance']
            for provenance_key, expected_file in (
                ('cg_support', 'component_support_audit.json'),
                ('generator_analysis', 'matchete_generator_analysis.json')):
                name_key=provenance_key+'_path'
                sha_key=provenance_key+'_sha256'
                original=str(pmeta.get(name_key,'')).replace('\\','/')
                if original != (AUDIT_REL/expected_file).as_posix():
                    raise StudyError('Provisional input audit source path drift: '+name_key)
                if pmeta.get(sha_key)!=audit_hashes[expected_file]:
                    raise StudyError('Provisional input source audit SHA256 drift: '+sha_key)
            rep=b['representation']
            cls=lookup[short]['SU2_dimensions']
            if (rep['d_s1'],rep['d_s2'],rep['d_f'],rep['alpha'])!=(cls['S1'],cls['S2'],cls['F'],row['alpha']):
                raise StudyError('Frozen benchmark representation mismatch')
            if not math.isclose(float(d['renormalization_scale_gev']),float(b['base_state']['mu_gev']),rel_tol=1e-12):
                raise StudyError('Source benchmark base scale mismatch')
            if not math.isclose(float(d['renormalization_scale_gev']),1e7,rel_tol=1e-12):
                raise StudyError('Unexpected frozen benchmark reference scale')
            if d['provenance'].get('validated_uv_component_conventions') is True:
                raise StudyError('Provisional input unexpectedly claims certified conventions')
            structural=validate_components(d,short[0],row['alpha'])
            ans=evaluate(d,retain_individual_terms=True)
            ref=independent_reference(d)
            if not ans['diagnostics']['uv_finite'] or not ans['diagnostics']['neutrino_mass_matrix_symmetric']:
                raise StudyError('Production mass fails UV finiteness/symmetry')
            cmp=_compare(ans,ref)
            if cmp['relative_matrix_error']>1e-6 and cmp['absolute_over_reference_sum_component_norms']>1e-9:
                raise StudyError('Independent spectral reference disagrees with production')
            sumerr=_reconstruct_individual(ans)
            if sumerr>1e-6:raise StudyError('Individual term reconstruction failed')
            parts=partition_exact(d,short[0],row['alpha'],ans)
            ch=_charge_rows(d,cls['S1'],weight_map[short],parts,short[0],row['alpha'])
            weight_rows[short]={'model':short,'dimensions':cls,'source_weight_audit_sha256':audit_hashes['charge_resolved_loop_weights.json'],
                                'channels':ch['native_group_audit_channels'],
                                'total_group_weight':ch['group_weight_total'],
                                'relative_physical_normalization_certified':False}
            masses=ans['total']['neutrino_masses_ev']
            if len(masses)!=3 or any(not math.isfinite(x) or x<0 for x in masses):
                raise StudyError('Invalid mass eigenvalues')
            row.update(status='conditional_uv_evaluated',m1_ev=masses[0],m2_ev=masses[1],m3_ev=masses[2],
                       candidate_topology_charge_channel_count=sum(bool(c['nonzero']) for c in
                           ch['native_group_audit_channels']),
                       source_reference_relative_error=cmp['relative_matrix_error'],
                       source_reference_absolute_over_sum=cmp['absolute_over_reference_sum_component_norms'],
                       pole_relative_production=ans['diagnostics']['relative_uv_pole_coefficient'],
                       pole_relative_reference=ref['diagnostics']['uv_pole_relative_to_natural_amplitude'],
                       cancellation_factor=ans['diagnostics']['cancellation_factor_sum_over_total'],
                       min_scalar_mass_gev=min(ans['scalar_eigenmasses_gev']),
                       max_scalar_mass_gev=max(ans['scalar_eigenmasses_gev']),
                       min_fermion_mass_gev=min(ans['fermion_eigenmasses_gev']),
                       max_fermion_mass_gev=max(ans['fermion_eigenmasses_gev']))
            detail={'format':'t3-lo-pattern-study-v1','model':short,'status':row['status'],
                    'physical_normalisation_certified':False,'source':{'benchmark':str(bfile.relative_to(root)),
                      'benchmark_sha256':row['benchmark_sha256'],'input':str(pfile),
                      'input_sha256':row['component_input_sha256'],'audit_hashes':audit_hashes},
                    'mu_gev':d['renormalization_scale_gev'],
                    'mnu_gev':ans['total']['complex_symmetric_matrix_gev'],
                    'mnu_ev':encoded(_mass(ans)*1e9),
                    'takagi_masses_ev_sorted':masses,
                    'scalar_mass_gev':ans['scalar_eigenmasses_gev'],
                    'fermion_mass_gev':ans['fermion_eigenmasses_gev'],
                    'production_diagnostics':ans['diagnostics'],
                    'independent_spectral_diagnostics':ref['diagnostics'],
                    'comparison':cmp,'individual_terms_reconstruction_relative':sumerr,
                    'production_individual_eigenstate_term_count':len(ans['individual_terms']),
                    'source_vertex_pair_matrices_gev':ans['vertex_pair_matrices_gev'],
                    'component_selection':structural,'source_charge_and_group':ch,
                    'scope':'restricted, direct UV p=0, frozen mu0, no EFT',
                    'unresolved':['four-spinor to Weyl physical mapping',
                        'cross-package CG/portal/T3 phase and normalisation']}
            write_json(output/'models'/f'{short}.json',detail)
            if do_scans:scan_rows.extend(_scans(d,ans,cls['S1'],short))
        except (StudyError,ValueError,KeyError,IndexError,FileNotFoundError,TypeError,
                RuntimeError,ArithmeticError,OverflowError) as exc:
            row['reason']=type(exc).__name__+': '+str(exc)
            write_json(output/'models'/f'{short}.json',{'model':short,'status':'blocked',
                        'physical_normalisation_certified':False,'reason':row['reason'],
                        'benchmark_sha256':row['benchmark_sha256'],
                        'component_input_sha256':row['component_input_sha256']})
        rows.append(row)
        print(f"{short:<5} {row['status']:<27} {row['reason'] or ('m3='+str(row['m3_ev'])+' eV')}",flush=True)
    fields=['model','family','alpha','dS1','dS2','dF','fermion_type',
            'candidate_topology_charge_channel_count','status','reason','m1_ev','m2_ev','m3_ev',
            'source_reference_relative_error','cancellation_factor','benchmark_sha256',
            'component_input_sha256','physical_normalisation_certified']
    with (output/'all16_summary.csv').open('w',newline='',encoding='utf-8') as fp:
        wr=csv.DictWriter(fp,fieldnames=fields,extrasaction='ignore');wr.writeheader();wr.writerows(rows)
    scan_fields=['model','mode','t','status','reason','prediction','m1_ev','m2_ev','m3_ev',
                 'matrix_scaling_residual','mass_norm_gev','cancellation_factor']
    with (output/'scans.csv').open('w',newline='',encoding='utf-8') as fp:
        wr=csv.DictWriter(fp,fieldnames=scan_fields,extrasaction='ignore');wr.writeheader();wr.writerows(scan_rows)
    write_json(output/'channel_weights.json',{'format':'t3-conditional-CG-weights-v1',
                'audit_hashes':audit_hashes,'models':weight_rows,
                'status':'native_group_CG_candidate_only_not_physical_normalization'})
    # A separate deliberately *group-only* common-kernel comparison. This is
    # never presented as a calculation of physical common-input neutrino mass.
    baseline_group=abs(_simple_expr('-2*sqrt(3)/3'))
    common=[]
    for m in models:
        short=m['model']
        if short in weight_rows:
            amplitude=sum((_simple_expr(c['complex_weight']) for c in
                           weight_rows[short]['channels'] if c['nonzero']),0j)
            common.append({'model':short,'candidate_weight_real':float(amplitude.real),
                           'candidate_weight_imag':float(amplitude.imag),
                           'candidate_weight_abs':float(abs(amplitude)),
                           'ratio_to_A_group_weight_abs':float(abs(amplitude)/baseline_group),
                           'channel_count':sum(bool(c['nonzero']) for c in weight_rows[short]['channels']),
                           'comparison':'common identical loop kernel and Yukawa/lambda factors; source CG ONLY',
                           'physically_certified':False})
    with (output/'common_kernel_group_only.csv').open('w',newline='',encoding='utf-8') as fp:
        columns=['model','candidate_weight_real','candidate_weight_imag','candidate_weight_abs',
                 'ratio_to_A_group_weight_abs','channel_count','comparison','physically_certified']
        writer=csv.DictWriter(fp,fieldnames=columns)
        writer.writeheader();writer.writerows(common)
    write_json(output/'common_kernel_group_only.json',{
        'assumptions':['identical heavy masses, identical scalar squared masses and identical nonzero off-diagonal scalar mass insertion in matched charge channels',
                       'identical Yukawa and mixing coefficients after choosing source component conventions',
                       'all internal loop kernels identical across channels and model families',
                       'this comparison does NOT establish physically consistent coupling normalisation across representations'],
        'status':'formal_common_kernel_group_weights_only',
        'rows':common,'physical_normalisation_certified':False})
    counts={'catalogue':16,'conditional_evaluated':sum(r['status']=='conditional_uv_evaluated' for r in rows),
            'blocked':sum(r['status']=='blocked' for r in rows),
            'physical_source_certified':0,'scan_passes':sum(r['status']=='conditional_uv_scanned' for r in scan_rows),
            'scan_blocked':sum(r['status']=='blocked' for r in scan_rows)}
    report={'format':'t3-conditional-pattern-study-all16-v1','counts':counts,
            'provenance':provenance,'models':rows,
            'physical_normalisation_certified':False,
            'notes':['No native B1/B2 physical mapping established',
                     '16/16 numerical evaluation does not imply 16/16 source certification',
                     'Historical strict 7 provisional / 9 blocked is not recomputed here',
                     'Yukawa scans act on source-derived candidate Weyl vertex arrays; lambda scan on off-diagonal scalar A',
                     'Zero/small masses may reflect charge selection, small coupling, or destructive interference; do not conflate']}
    write_json(output/'all16_summary.json',report)
    plot_errors=_plots(output,rows,scan_rows) if do_plots else []
    handoff=f"""# Chat A handoff — restricted all-16 direct UV

- Reviewed GitHub main: `{BASELINE_REV}`; actual local HEAD: `{provenance['revision']['local_head']}`; dirty worktree: `{provenance['revision']['local_worktree_dirty']}`.
- This run: `{counts['conditional_evaluated']}/16` conditional UV evaluated, `{counts['blocked']}/16` blocked, `0/16` physically certified.
- Input provider: `{input_dir}`; frozen input and audit SHA256 are recorded in per-model JSON.
- Exact conditional computation: existing `core.evaluate` in canonical left-Weyl/complex-scalar basis; independent reference uses `M^dagger f(M M^dagger)` and scalar Cartesian Hessian.
- Source-CG rule (as audited): `T3Y1CG * T3Y2CG * conjugate(T3MixCG)` in common source-weight indices; NOT certified cross-package normalisation.
- Yukawa scans: conditional linearity for topology-only cross terms, tested numerically; lambda scan varies only `A[S1,S2]` off-diagonal; mass scan scales heavy `M`.
- Critical unresolved: native four-spinor chirality/phase/normalisation, RGBeta/Matchete scalar portals/T3 sign and component-axis mapping, physical fermion multiplicities.
- This package **is not safe for asserting physical neutrino masses**. It is safe for numerical regression and source-conditional structure. No EFT/running/threshold matching.
- Plot warnings: {plot_errors}.
- No previously existing output directory overwritten.
"""
    (output/'HANDOFF.md').write_text(handoff,encoding='utf-8')
    write_json(output/'HANDOFF.json',{'format':'t3-chat-A-handoff-v1',
              'reviewed_remote_sha':BASELINE_REV,'local_revision':provenance['revision'],
              'checks':counts,'source_audit_status':guard['status'],
              'unresolved':['B1 Weyl spinor dictionary','B2 SU2/RGBeta/Matchete dictionary'],
              'safe_to_integrate_for_physical_mass':False,
              'safe_to_integrate_for_provisional_numerics':counts['blocked']==0,
              'physical_source_certified':False})
    return report
