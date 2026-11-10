"""T3 Batch 4: conditional numerical finite *native* Matchete F-only LLHH terms.

This is an ALTERNATIVE one-loop MF boundary, replacing (not adding to) the
Batch-1 analytic triangle. All PL Laurent kernels from the recorded 16-model
native one-loop MatchedEFT are evaluated, including extra quartics and B/C
X/Y mixing. Conditional operation: formally remove the 1/epsilon pole and
retain 1+log(muF**2/Mk**2). THIS IS NOT an independently verified UV-hard
minus EFT1-soft subtraction or an independently established GammaCC/Weyl map.

This module never interprets the result as physical MSbar matching. The exact
source inventory, each native quartic symbol, numerical heavy mass, and
benchmark assumptions are recorded for reproducibility.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Mapping

import numpy as np

MF=100_000.0
LOOP=1/(16*math.pi**2)
TRIOS=Path('output/general_observable/lo_certification_spinor/pole_log_e515d82d118284ba_bridge01')
V3=Path('output/parallel_t3/wilson/additional_quartic_dictionary_conditional_v3.json')
B_RUN=Path('output/parallel_t3/background/B_all16_qualitative_v1')
MIXED={'B_-1','C_-1'}
EXPECTED_PL={'A_-4':1,'A_-2':2,'A_+0':2,'B_-3':1,'B_-1':34,'B_+1':1,
             'C_-3':1,'C_-1':34,'C_+1':1,'D_-2':2,'D_+0':2,'D_+2':1,
             'E_-4':1,'E_-2':2,'E_+0':2,'E_+2':1}
SOURCE_FORMAT='t3-f-first-source-resolved-pole-log-bridge-v1'
MAPPING_FORMAT='t3-parallel-a-extra-quartic-polynomial-dictionary-conditional-v3'
NATIVE_SCHEMA='t3-all16-native-pole-dropped-finite-candidate-batch4-v1'
_SOURCE=re.compile(r'^([+-])\s+\((.*?)hbar\*(.*)\)/(\(.*\)|Coupling\[MF.*\])$', re.DOTALL)
_COUPLING=re.compile(r'(?P<bar>Bar\[)?Coupling\[(?P<name>[a-zA-Z][a-zA-Z0-9]*),')
_YUK={'y1','y2','y1Conj','y2Conj'}
_MIX=re.compile(r'^T3MEU([12])([12])$')
_EXTRA=re.compile(r'^lambdaT3ExtraQ\d+I\d+$')


def digest(path:Path)->str:
    return sha256(Path(path).read_bytes()).hexdigest()


def read_json(path:Path)->dict:
    obj=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(obj,dict):raise ValueError('Expected object: '+str(path))
    return obj


def relative(x:np.ndarray,y:np.ndarray)->float:
    return float(np.linalg.norm(np.asarray(x)-np.asarray(y))/max(np.linalg.norm(y),1e-100))


def _balanced_close(text:str,opening:int)->int:
    """Find closing ] for one nested Matchete head (no eval/exec)."""
    depth=0
    for i in range(opening,len(text)):
        if text[i]=='[':depth+=1
        elif text[i]==']':
            depth-=1
            if depth==0:return i
    raise ValueError('Unbalanced native Matchete coupling head')


def parse_scalar_number(s:str)->float:
    """Only exact numeric syntax present in the tracked Laurent kernels."""
    v=s.strip().rstrip('*')
    if not v:return 1.0
    pieces=v.split('*')
    out=1.0
    for part in pieces:
        if re.fullmatch(r'\d+(?:/\d+)?',part):
            a,*b=part.split('/')
            y=float(a)/(float(b[0]) if b else 1.0)
        elif re.fullmatch(r'Sqrt\[\d+(?:/\d+)?\]',part):
            y=math.sqrt(parse_scalar_number(part[5:-1]))
        else:raise ValueError('Unsupported native scalar prefactor: '+repr(part))
        out*=y
    if not math.isfinite(out) or out<=0:raise ValueError('Invalid source prefactor')
    return out


@dataclass(frozen=True)
class ParsedKernel:
    overall:float
    yukawas:tuple[tuple[str,bool],...]
    scalars:tuple[tuple[str,bool,int],...]
    source:str


def parse_native_kernel(s:str)->ParsedKernel:
    """Narrow lossless evaluator for the *actual* 88 PL finite kernels.

    Source H/L/epsilon/spinor structures are checked as a unique invariant,
    NOT silently parsed into a different operator. Only scalar prefactors,
    coupling monomials, and flavour labels enter the numerical evaluation.
    """
    if not isinstance(s,str):raise ValueError('Kernel must be InputForm text')
    m=_SOURCE.fullmatch(s.strip())
    if not m:raise ValueError('Unexpected native finite kernel syntax')
    sign, pref, body, den=m.groups()
    if not body.endswith(']') or 'DiracProduct[GammaCC, Proj[-1]]' not in body:
        raise ValueError('Not the native PL GammaCC Weinberg structure')
    if body.count('Field[H, Scalar,')!=2 or body.count('Field[l, Fermion,')!=2:
        raise ValueError('Not exactly two Higgs and two lepton fields')
    if body.count('CG[Bar[eps[SU2L]]')!=2 or 'CG[eps[SU2L]' in body:
        raise ValueError('Noncanonical native PL epsilon invariant')
    if body.count('hbar') or 'Log[' in body or 'ϵ' in body or 'ε' in body:
        raise ValueError('Kernel was not separated into constant/pole/log')
    # Exactly one heavy MF denominator with optional native numerical factor.
    if den.startswith('(') and den.endswith(')'):den=den[1:-1]
    md=re.fullmatch(r'(?:(Sqrt\[3\]|3)\*)?Coupling\[MF, (?:\{\}|\{Index\[d\$\$1, NFlavor\]\}), 0\]',den)
    if not md:raise ValueError('Unexpected heavy mass factor or index in native kernel')
    normal=parse_scalar_number(pref)/parse_scalar_number(md.group(1) or '')
    if sign=='-':normal=-normal
    result=[]
    spans=[]
    for x in _COUPLING.finditer(body):
        opening=body.index('[',x.start()+len(x.group('bar') or ''))
        last=_balanced_close(body,opening)
        barred=bool(x.group('bar'))
        end=last+1
        if barred:
            if end>=len(body) or body[end]!=']':raise ValueError('Unbalanced Bar[Coupling]')
            end+=1
        power_match=re.match(r'\^(\d+)',body[end:])
        power=int(power_match.group(1)) if power_match else 1
        if not 1<=power<=8:raise ValueError('Unsupported coupling multiplicity')
        end+=len(power_match.group(0)) if power_match else 0
        name=x.group('name')
        if name in _YUK:
            actual=body[opening+1:last]
            required=re.fullmatch(
                r'\w+, \{Index\[d\$\$([12]), Flavor\], Index\[d\$\$1, NFlavor\]\}, 0',
                actual)
            if required is None:
                raise ValueError('Unexpected Yukawa source heavy/external flavour indices')
        if name not in _YUK and name!='lambdaT3' and not _EXTRA.fullmatch(name) and not _MIX.fullmatch(name):
            raise ValueError('Unmapped native coupling '+name)
        result.append((name,barred,power))
        spans.append((x.start(),end))
    if not spans or any(b>c for (_,b),(c,_) in zip(spans,spans[1:])):
        raise ValueError('Overlapping or missing native couplings')
    # Enforce a fixed source polynomial grammar: only the other recognized
    # factor heads remain. Unsupported extra numerators must not be dropped.
    rest=body
    for start,end in reversed(spans):rest=rest[:start]+rest[end:]
    atoms=[part for part in rest.split('*') if part.strip()]
    if len(atoms)!=5:
        raise ValueError('Unexpected native numerator factor inventory: '+str(len(atoms)))
    if sum(x.startswith('CG[Bar[eps[SU2L]]') for x in atoms)!=2 or \
       sum(x.startswith('Field[H, Scalar,') for x in atoms)!=2 or \
       sum(x.startswith('NCM[Transp[Field[l, Fermion,') for x in atoms)!=1:
        raise ValueError('Unknown native external tensor/CG factor')
    for atom in atoms:
        if atom.startswith('CG[') and not re.fullmatch(
           r'CG\[Bar\[eps\[SU2L\]\], \{Bar\[Index\[d\$\$\d+, SU2L\[fund\]\]\], Bar\[Index\[d\$\$\d+, SU2L\[fund\]\]\]\}\]',atom):
            raise ValueError('Unexpected native epsilon CG')
        if atom.startswith('Field[') and not re.fullmatch(
           r'Field\[H, Scalar, \{Index\[d\$\$\d+, SU2L\[fund\]\]\}, \{\}\]',atom):
            raise ValueError('Unexpected native Higgs field factor')
        if atom.startswith('NCM['):
            if not atom.endswith(']') or atom.count('Field[l, Fermion,')!=2 or \
               atom.count('DiracProduct[GammaCC, Proj[-1]]')!=1 or \
               not atom.startswith('NCM[Transp[Field[l, Fermion,'):
                raise ValueError('Unexpected native PL fermion/spinor structure')
    yuk=[];scalar=[]
    for name,barred,power in result:
        if name in _YUK:
            yuk.extend([(name,barred)]*power)
        else:scalar.append((name,barred,power))
    if len(yuk)!=2 or len([x for x in scalar if x[0].startswith('lambda')])!=1:
        raise ValueError('Source must have precisely two Yukawas and one quartic')
    return ParsedKernel(normal,tuple(yuk),tuple(scalar),s)


def matchete_quartics(entry:Mapping, native_state:Mapping[str,np.ndarray])->tuple[dict,dict]:
    """Read full V3 conditional linear map r=A m and solve it, no zipped aliases."""
    mats=list(entry.get('matchete_symbols',[]));rgbs=list(entry.get('rgbeta_symbols',[]))
    if len(mats)!=len(set(mats)) or len(rgbs)!=len(set(rgbs)):
        raise ValueError('Duplicated quartic dictionary symbols')
    if not mats:return {},{'mapped_extra_count':0,'matrix_condition':None}
    raw=entry['map_rgbeta_from_matchete']
    A=np.asarray(raw['real'],float)+1j*np.asarray(raw['imag'],float)
    if A.shape!=(len(rgbs),len(mats)) or len(rgbs)!=len(mats):
        raise ValueError('V3 conditional quartic dictionary not square')
    if any(x not in native_state for x in rgbs):
        raise ValueError('Running RGBeta lacks extra native quartic coordinate: '+str(sorted(set(rgbs)-set(native_state))))
    rhs=np.asarray([complex(np.asarray(native_state[x]).item()) for x in rgbs])
    if not np.isfinite(rhs).all():raise ValueError('Nonfinite native quartic input')
    cond=float(np.linalg.cond(A))
    if not math.isfinite(cond) or cond>1e10:raise ValueError('Singular conditional quartic map')
    m=np.linalg.solve(A,rhs)
    if relative(A@m,rhs)>1e-9:raise ValueError('Conditional Matchete/RGBeta quartic reconstruction failure')
    return dict(zip(mats,m)),{'mapped_extra_count':len(mats),'matrix_condition':cond,
        'basis':'V3 source-conditioned polynomial mapping; physical coupling normalization not independently fixed'}


def evaluate_kernels(kernels:list[ParsedKernel], masses:np.ndarray,
                     yukawas:Mapping[str,np.ndarray],quartics:Mapping[str,complex],
                     mixing_U:np.ndarray|None, *,mu_f:float=MF,
                     native_spinor_factor:float=1.0) -> tuple[np.ndarray,np.ndarray,list[dict]]:
    """Return [all-kernel,flavour], [kernel,heavy-index,flavour], audit rows.

    Multiplying by native_spinor_factor explicitly selects an UNVERIFIED
    GammaCC-to-canonical-Weyl overall factor; default 1 reproduces the
    existing conditional convention, not an independent amplitude derivation.
    """
    mf=np.asarray(masses,float)
    if mf.shape!=(3,) or not np.isfinite(mf).all() or np.any(mf<=0):
        raise ValueError('Need three positive same-scale heavy eigenmasses')
    if mu_f<=0 or not np.isfinite(native_spinor_factor):
        raise ValueError('Invalid matching scale or declared spinor factor')
    U=None if mixing_U is None else np.asarray(mixing_U,complex)
    if U is not None and (U.shape!=(2,2) or not np.allclose(U.conj().T@U,np.eye(2),atol=1e-9)):
        raise ValueError('Invalid source native X/Y U matrix')
    per=np.zeros((len(kernels),3,3,3),complex)
    rows=[]
    for n,kernel in enumerate(kernels):
        scalar=complex(kernel.overall)
        channel=[]
        for name,barred,power in kernel.scalars:
            if name.startswith('lambda'):
                if name not in quartics:raise ValueError('Unmapped named quartic '+name)
                val=complex(quartics[name])
                channel.append(name)
            else:
                match=_MIX.fullmatch(name)
                if match is None or U is None:raise ValueError('Native mixing spurion not available: '+name)
                val=complex(U[int(match[1])-1,int(match[2])-1])
            if barred:val=val.conjugate()
            scalar*=val**power
        y=[]
        for name,barred in kernel.yukawas:
            if name not in yukawas:raise ValueError('Missing heavy-basis Yukawa '+name)
            value=np.asarray(yukawas[name],complex)
            if value.shape!=(3,3) or not np.isfinite(value).all():
                raise ValueError('Invalid heavy-basis Yukawa '+name)
            y.append(value.conj() if barred else value)
        for mass_index,heavy_mass in enumerate(mf):
            pref=native_spinor_factor*LOOP*scalar/float(heavy_mass)
            pref*=1+math.log(float(mu_f)**2/float(heavy_mass)**2)
            w=pref*np.outer(y[0][:,mass_index],y[1][:,mass_index])
            per[n,mass_index]=(w+w.T)/2
        rows.append({'source_kernel_sha256':sha256(kernel.source.encode()).hexdigest(),
                     'yukawa_pair':[[name,bool(bar)] for name,bar in kernel.yukawas],
                     'native_quartic':channel[0],
                     'native_scalar_prefactor_re':float(scalar.real),
                     'native_scalar_prefactor_im':float(scalar.imag),
                     'by_fermion_norm':[float(np.linalg.norm(x)) for x in per[n]]})
    total=per.sum(axis=(0,1))
    if not np.isfinite(total).all() or not np.allclose(total,total.T,rtol=1e-10,atol=1e-25):
        raise ValueError('Nonfinite or nonsymmetric native finite C5')
    return total,per,rows


def native_all16_boundary(project:Path,model:str,*,native_spinor_factor:float=1.0)->tuple[np.ndarray,np.ndarray,dict]:
    """Evaluate one model with actual SHA-validated local B, V3 and Laurent sources."""
    from Numerical.parallel_t3.background.api import Background,load_formal_boundary
    from Numerical.parallel_t3.background.contracts import tag
    from Numerical.parallel_t3.wilson.fermion_hard_triangle_batch1 import _uv_at_mf,complex_matrix
    from Numerical.matching.GeneralT3FlavourC5 import prepare_heavy_flavour

    root=Path(project).resolve()
    if model not in EXPECTED_PL:raise ValueError('Model outside all-16 set')
    slug=tag(model)
    trios_path=root/TRIOS/slug/'trios.json'
    source=read_json(trios_path)
    if (source.get('format')!=SOURCE_FORMAT or source.get('model')!=model or
            source.get('PL_groups')!=EXPECTED_PL[model] or
            source.get('PR_groups')!=EXPECTED_PL[model] or
            source.get('exact_kernel_groups')!=2*EXPECTED_PL[model] or
            source.get('raw_pole_residue_equals_finite_constant_equals_log_coefficient') is not True or
            source.get('EFT1_soft_subtraction_performed') is not False or
            source.get('counterterm_subtraction_performed') is not False or
            source.get('physical_MSbar_hard_coefficient') is not None or
            source.get('master_Weyl_phase') is not None):
        raise ValueError('Unexpected native one-loop source or physical-completion claim')
    all_records=source.get('records',[])
    pl=[x for x in all_records if x.get('chirality')=='PL']
    pr=[x for x in all_records if x.get('chirality')=='PR']
    if len(pl)!=EXPECTED_PL[model] or len(pr)!=len(pl):
        raise ValueError('All source PL/PR Laurent groups required')
    seen=set();kernels=[]
    for row in all_records:
        k=row.get('kernel_inputform','')
        if (row.get('kernel_sha256')!=sha256(k.encode()).hexdigest() or
                row.get('coefficient_factor')!='1 + 1/ϵ + Log[μbar2/MF^2]' or
                row.get('pole_residue_equals_constant_kernel') is not True or
                row.get('log_scale_coefficient_equals_pole_residue') is not True):
            raise ValueError('Changed native Laurent triplet source')
        ids=row.get('source_indices',{})
        if set(ids)!={'constant','pole','scale_log'} or any(not isinstance(i,int) or i in seen for i in ids.values()):
            raise ValueError('Incomplete native Laurent term indices')
        seen.update(ids.values())
        if row['chirality']=='PL':kernels.append(parse_native_kernel(k))
    dictpath=root/V3
    v3=read_json(dictpath)
    if v3.get('format')!=MAPPING_FORMAT or v3.get('model_count')!=16 or v3.get('n_failed_models')!=0:
        raise ValueError('V3 conditional Matchete quartic dictionary changed')
    entry=v3['models'][model]
    base=root/B_RUN/slug
    bg=Background.load(base)
    formal,fmeta=load_formal_boundary(base)
    if (fmeta.get('model')!=model or bg.manifest.get('model')!=model or
            fmeta.get('active_BSM_fields_below_MF')!=['S1','S2'] or
            fmeta.get('all_three_heavy_fermions_decoupled_at_single_MF') is not True):
        raise ValueError('B manifest/model/common heavy fermion threshold mismatch')
    if bg.manifest['input_scales_gev']['F']!=MF:
        raise ValueError('B UV and EFT1 trajectories use different MF')
    native_state={k:np.asarray(bg.arrays['eft1__'+k][0]) for k in bg.manifest['beta_symbol_order_eft1']}
    extra,mapaudit=matchete_quartics(entry,native_state)
    quartics={'lambdaT3':complex(np.asarray(native_state['lambdaT3']).item()),**extra}
    used={x[0] for k in kernels for x in k.scalars if x[0].startswith('lambda')}
    if not used<=set(quartics):raise ValueError('Source requires missing quartics: '+repr(used-set(quartics)))
    need={name for kernel in kernels for name,_ in kernel.yukawas}|{'y1','y2'}
    uv=_uv_at_mf(bg,model)
    masses,ys,heavy_audit=prepare_heavy_flavour(uv,need)
    saved_masses=np.asarray(formal['heavy_mass_eigenvalues_gev'],float)
    if not np.allclose(masses,saved_masses,rtol=2e-8,atol=5e-5):
        raise ValueError('Native one-loop source evolved heavy mass eigenbasis does not match B F-tree')
    mixed=None
    if model in MIXED:
        meta=fmeta.get('mixed_native_diagnostics',{})
        if 'mixing_U_muS' not in meta:raise ValueError('Missing original B/C MF mixed source U')
        mixed=complex_matrix(meta['mixing_U_muS'])
    c5,per,records=evaluate_kernels(kernels,masses,ys,quartics,mixed,
                                    native_spinor_factor=native_spinor_factor)
    channels={}
    for i,row in enumerate(records):
        name=row['native_quartic']
        channels[name]=channels.get(name,np.zeros((3,3),complex))+per[i].sum(axis=0)
    if relative(sum(channels.values(),np.zeros((3,3),complex)),c5)>1e-11:
        raise ValueError('Native quartic C5 channels do not close')
    data={
        'schema':NATIVE_SCHEMA,'model':model,
        'native_pole_drop_is_physical_MSbar_subtraction':False,
        'EFT1_soft_subtraction_performed':False,'counterterm_subtraction_performed':False,
        'native_gamma_cc_to_master_Weyl_signed_factor':float(native_spinor_factor),
        'native_gamma_cc_factor_independently_certified':False,
        'benchmark_hbar':LOOP,'matching_scale_gev':MF,
        'operator_inputform_factor':'1 + 1/epsilon + Log[mubar2/Mk2]',
        'numerical_FORMAL_pole_dropped_factor':'1 + Log[MF2/Mk2]',
        'one_common_MF_all_fermions_out':True,'both_scalars_retained_at_MF':True,
        'native_PL_kernel_count':len(kernels),'native_PR_excluded_count':len(pr),
        'native_quartic_channels':{name:{'real':a.real.tolist(), 'imag':a.imag.tolist(),
                         'norm':float(np.linalg.norm(a))} for name,a in channels.items()},
        'kernel_records':records,
        'source_sha256':{'native_trios':digest(trios_path),'conditional_v3':digest(dictpath),
                         'B_background_npz':digest(base/'background.npz'),
                         'B_fermion_boundary_manifest':digest(base/'fermion_boundary_manifest.json')},
        'heavy_eigenmasses_gev':np.asarray(masses,float).tolist(),
        'heavy_flavour_basis_diagnostics':str(heavy_audit),
        'extra_quartic_map':mapaudit,
        'RGBeta_native_to_Matchete_quartic_coefficients_PHYSICAL_mapping_certified':False,
        'actual_UV_hard_subtraction_certified':False,
        'full_one_loop_matching_complete':False,
        'physics_status':'conditional_native_finite_source_scenario_only'}
    return c5,per,data
