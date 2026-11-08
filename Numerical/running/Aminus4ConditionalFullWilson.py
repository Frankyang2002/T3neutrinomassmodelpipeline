"""A(-4) source-backed *conditional* 45-Weyl/12-real-scalar Wilson transport.

UV --integrate out ALL F--> EFT1(SM+S1+S2) --integrate out BOTH S--> SMEFT.
This program stops at the END OF EFT1. It neither performs scalar matching nor
identifies the signed Matchete GammaCC -> Weyl conversion. The rho=+1 source
boundary is a declared convention, not a physically certified C(M_F).

No production RGE or historical inputs are modified. Uses the already tested
Fonseca six-channel *candidate* and the corrected RGBeta 17-beta EFT1 samples.
The source Wilson beta's full gauge-dependent extension is not independently
loop-certified. Only global rho -> -rho is inferred by real-linearity; an
arbitrary complex normalization phase is NOT generally factorizable because
beta(C) contains C.conj().
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import Any, Mapping
import json

import numpy as np

BASIS = Path('output/a_minus4_next_batch/fresh_index_production/running/A_m4_eft1_tensor_basis.json')
TENSORS = Path('output/a_minus4_next_batch/fresh_index_production/running/A_m4_eft1_tensors.npz')
BETA = Path('output/a_minus4_next_batch/fresh_index_production/A_m4_fixed_17beta_compatible.json')
PRODUCTION_SUMMARY = Path('output/a_minus4_next_batch/fresh_index_production/summary.json')
BOUNDARY = Path('output/general_observable/a_minus4_boundary_batch1/declared_convention_boundary.json')
MASTER = Path('output/general_observable/a_minus4_boundary_batch1/master_basis_candidate.json')
FORMAT = 'a-minus4-declared-conditional-full-wilson-eft1-v1'
SOURCE_HASH = '147b25f172ebcc3644e3ae1be8a7686fe66e1b7c1673ca17dd10e444ce0fff0d'
MODES = ('fonseca_eq22_six_term',)


class ConditionalBlocked(ValueError):
    """Changed or uncertified numerical input; do not silently coerce."""


def require(test: bool, message: str) -> None:
    if not test:
        raise ConditionalBlocked(message)


def file_sha(path: Path) -> str:
    require(path.is_file(), 'Missing prerequisite file: '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    require(path.is_file(), 'Missing prerequisite JSON: '+str(path))
    x = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(x,dict), 'Expected JSON object: '+str(path))
    return x


def _orbit_size(i: int, j: int, a: int, b: int) -> int:
    return (1 if i==j else 2)*(1 if a==b else 2)


def _put_orbit(out: np.ndarray, i: int, j: int, a: int, b: int, z: complex):
    """Indices are 0-based. Fail if source coordinates would double-count."""
    for p,q in ((i,j),(j,i)):
        for u,v in ((a,b),(b,a)):
            old=out[p,q,u,v]
            require(old==0 or abs(old-z) < 1e-26,
                    'Inconsistent declared permutation orbit')
            out[p,q,u,v]=z


def declared_tensor(report: Mapping[str,Any], fi: Mapping[str,int], *, rho: int=1) -> np.ndarray:
    """36 local-lepton, real-scalar independent coordinates -> full EFT1 C.

    The lepton local ordering is L[1,1],L[1,2],L[2,1],L[2,2],L[3,1],L[3,2].
    It is NOT the global 45-Weyl ordering; the first 36 global Weyl slots
    belong to quarks in the numerical EFT1 tensor basis.
    """
    require(type(rho) is int and rho in (-1,1), 'Choose rho=+1 or -1')
    require(report.get('format')=='t3-a-minus4-declared-weyl-uv-coupling-dictionary-v1'
            and report.get('model')=='A_-4'
            and report.get('independent_coordinate_count')==36,
            'Wrong declared Matchete/UV boundary source')
    require(report.get('native_matchete_to_declared_Weyl_factor') is None and
            report.get('physical_C_ijab_at_MF') is None and
            report.get('physical_certification') is False and
            report.get('declared_gamma_cc_choice',{}).get('native_GammaCC_to_C0_measured') is None,
            'Signed native normalization status changed; requires independent review')
    require(report.get('master_full_shape')==[45,45,12,12],
            'Declared boundary shape does not match full EFT1')
    pos=report.get('declared_rho_plus_one_independent_coordinates')
    neg=report.get('rho_minus_one_fixed_native_y1_alternative_coordinates')
    require(isinstance(pos,list) and isinstance(neg,list) and len(pos)==len(neg)==36,
            'Missing full signed 36-coordinate declared inventory')
    key=lambda row:(row['i'],row['j'],row['a'],row['b'])
    pmap={key(row):row for row in pos};nmap={key(row):row for row in neg}
    require(len(pmap)==len(nmap)==36 and set(pmap)==set(nmap),
            'Duplicated/mismatched declared independent coordinate keys')
    local=[int(fi[f'L[{p},{w}]']) for p in (1,2,3) for w in (1,2)]
    require(len(set(local))==6 and all(0<=k<45 for k in local),
            'Numerical 45-Weyl basis has invalid lepton placement')
    out=np.zeros((45,45,12,12),np.complex128)
    for k,plus in pmap.items():
        i,j,a,b=(int(v) for v in k)
        require(1<=i<=j<=6 and 5<=a<=6 and 7<=b<=12,
                'Declared support exceeds LL[S1bar][S2]')
        minus=nmap[k]
        orbit=_orbit_size(i,j,a,b)
        require(plus.get('distinct_orbit_size')==minus.get('distinct_orbit_size')==orbit,
                'Wrong full-sum orbit normalization')
        z=complex(float(plus['value_re_gev_inv']),float(plus['value_im_gev_inv']))
        zm=complex(float(minus['value_re_gev_inv']),float(minus['value_im_gev_inv']))
        require(np.isfinite(z) and np.isfinite(zm) and abs(z+zm)<=1e-12*max(abs(z),1e-30),
                'rho branch coordinates are not exact sign alternatives')
        _put_orbit(out,local[i-1],local[j-1],a-1,b-1,rho*z)
    require(np.count_nonzero(out)>0 and np.isfinite(out).all(),
            'Declared Wilson tensor empty/nonfinite')
    return out


def compare_master_source(declared: Mapping[str,Any], master: Mapping[str,Any],
                          c: np.ndarray, fi: Mapping[str,int]) -> float:
    """Source-independent master candidate must reproduce declared rho=+1."""
    require(master.get('model')=='A_-4' and master.get('physical_certification') is False and
            master.get('physical_C_ijab_at_MF') is None,
            'Master candidate unexpectedly physically certified or changed')
    rows=master.get('conditional_master_weight_basis_independent_entries')
    require(isinstance(rows,list) and len(rows)==36,
            'Native-master triplet candidate lacks 36 independent entries')
    source=np.zeros_like(c)
    local=[int(fi[f'L[{p},{w}]']) for p in (1,2,3) for w in (1,2)]
    seen=set()
    for row in rows:
        i,j,a,b=(int(row[k]) for k in ('i','j','a','b'))
        k=(i,j,a,b)
        require(k not in seen and 1<=i<=j<=6 and 5<=a<=6 and 7<=b<=12,
                'Master candidate contains unexpected coordinate')
        seen.add(k)
        require(int(row['distinct_orbit_size'])==_orbit_size(i,j,a,b),
                'Master candidate orbit-size mismatch')
        z=complex(float(row['coefficient_real_if_rho_1']),
                  float(row['coefficient_imag_if_rho_1']))
        _put_orbit(source,local[i-1],local[j-1],a-1,b-1,z)
    scale=max(np.max(np.abs(source)),1e-300)
    residual=float(np.max(np.abs(source-c))/scale)
    require(residual<3e-12,'Declared signed boundary differs from source-anchored master CG')
    return residual


def _encode(z: np.ndarray) -> dict:
    return {'real':np.asarray(z).real.tolist(),'imag':np.asarray(z).imag.tolist()}


def _ward_max(c:np.ndarray, samples:list, names:tuple[str,...], *, run_index:int=0) -> dict:
    from RGE.running.intermediate.GeneralEFT1SevenTermFlow import gauge_ward_residual
    result={}
    for name in names:
        g,pairs=samples[run_index].gauge[names.index(name)]
        result[name]=float(max(gauge_ward_residual(c,t,s) for t,s in pairs))
    return result


def load_contract(root:Path, *, basis_path:Path=BASIS, tensors_path:Path=TENSORS,
                  beta_path:Path=BETA, summary_path:Path=PRODUCTION_SUMMARY,
                  boundary_path:Path=BOUNDARY, master_path:Path=MASTER) -> dict:
    """Read-only verification: source-digest, 17 betas, samples, basis and C0.

    Imports the existing repository's inspected EFT1 adapter only here; pure
    declaration algebra above remains directly unit-testable without Wolfram.
    """
    from Numerical.running import Aminus4EFT1Couplings as eft
    from RGE.running.intermediate.GeneralEFT1SevenTermFlow import Couplings, TensorBasis

    root=Path(root).resolve()
    def loc(p:Path)->Path: return p if p.is_absolute() else root/p
    paths={k:loc(v) for k,v in {'basis':basis_path,'tensors':tensors_path,'beta':beta_path,
           'production_summary':summary_path,'boundary':boundary_path,'master':master_path}.items()}
    require(file_sha(root/'RGE/running/rgbeta/T3RGBetaModel.wl')==SOURCE_HASH,
            'Corrected production RGBeta source hash differs; do not reuse stale trajectory')
    source=read_json(paths['production_summary'])
    require(source.get('model')=='A_-4' and source.get('native_17_parameter_beta_integrity') is True
            and source.get('ode_renormalisable_couplings_7_scales_computed') is True
            and source.get('physical_full_Wilson_ODE_enabled') is False,
            'Regenerated A(-4) production audit is missing/uncertified')
    beta=read_json(paths['beta'])
    eft.validate_beta_payload(beta)
    require(source.get('new_beta_sha256')==file_sha(paths['beta']),
            'Corrected beta payload is not the previously verified production beta')
    bg=read_json(paths['basis'])
    require(bg.get('model')=='A_-4' and
            bg.get('matchete_scalar_and_Wilson_basis_map_certified') is False and
            bg.get('physical_full_Wilson_running_ready') is False,
            'Incorrect or prematurely certified EFT1 tensor basis')
    require(bg.get('fermion_basis')==list(eft.F_LABELS) and
            bg.get('real_scalar_basis')==list(eft.S_LABELS),
            'Incorrect fermion/scalar basis ordering')
    require(bg.get('rge_payload_sha256')==file_sha(paths['beta']),
            'Tensor samples were not generated with the corrected 17-beta payload')
    prov=eft.source_contract(root,beta)['source_sha256']
    require(bg.get('provenance',{}).get('source_sha256')==prov,
            'EFT1 tensor samples do not match current executable source')
    declared=read_json(paths['boundary'])
    master=read_json(paths['master'])
    require(abs(float(declared.get('mu_F_gev',0))-float(bg.get('mu_gev',-1)))<1e-8,
            'Different fermion matching threshold scales')
    sha_values=set(declared.get('source_files_sha256',{}).values())
    require(file_sha(paths['master']) in sha_values,
            'Declared operator dictionary is not tied to the source master candidate')
    c0=declared_tensor(declared,eft.FI,rho=1)
    cross=compare_master_source(declared,master,c0,eft.FI)
    scales=np.asarray(bg.get('sampled_scales_gev',[]),float)
    rows=bg.get('running_couplings')
    require(scales.shape==(7,) and np.all(np.isfinite(scales)) and
            np.all(np.diff(scales)<0) and float(scales[0])>float(scales[-1])>0 and
            isinstance(rows,list) and len(rows)==7,
            'Corrected production ODE must contain seven descending samples')
    require(set(rows[0])==set(beta['report_betas']) and all(set(x)==set(rows[0]) for x in rows),
            'Saved running couplings do not contain all 17 EFT1 parameters')
    with np.load(paths['tensors'],allow_pickle=False) as data:
        require(set(('mu_gev','Y','lambda4'))<=set(data.files),
                'Saved NPZ lacks complete tensor samples')
        require(np.allclose(data['mu_gev'],scales,rtol=1e-12,atol=1e-8),
                'NPZ/JSON scale grid mismatch')
        y=np.asarray(data['Y']);lam=np.asarray(data['lambda4'])
        require(y.shape==(7,45,45,12) and lam.shape==(7,12,12,12,12),
                'Saved EFT1 samples have wrong complete tensor dimensions')
        base_g=eft.gauge_generators()
        for name in ('SU3','SU2','U1'):
            for rep in ('fermion','scalar'):
                key=f'generator_{name}_{rep}'
                require(key in data.files and
                        np.array_equal(data[key],np.stack(base_g[name][rep])),
                        'Saved gauge representation changed: '+key)
    for k in (0,3,6):
        point=eft.tensor_point(float(scales[k]),rows[k])
        require(np.max(abs(point.Y-y[k]))<1e-12 and
                np.max(abs(point.lambda4-lam[k]))<1e-11,
                'Saved 45-Weyl tensors disagree with reconstructed corrected couplings')
    names=('SU3','SU2','U1')
    input_samples=[]
    for k in range(7):
        g=tuple((float(rows[k][{'SU3':'g3','SU2':'g2','U1':'gY'}[name]]),
                 tuple(zip(base_g[name]['fermion'],base_g[name]['scalar'])))
                 for name in names)
        input_samples.append(Couplings(y[k],lam[k],g))
    basis=TensorBasis(45,12,
        tuple((eft.FI[f'L[{p},1]'],eft.FI[f'L[{p},2]']) for p in (1,2,3)),
        (0,1,2,3))
    return {'basis':basis,'c0':c0,'scales':scales[::-1].copy(),
            'samples':input_samples[::-1], 'boundary_sha256':file_sha(paths['boundary']),
            'master_sha256':file_sha(paths['master']), 'fixed_beta_sha256':file_sha(paths['beta']),
            'production_summary_sha256':file_sha(paths['production_summary']),
            'tensor_npz_sha256':file_sha(paths['tensors']),
            'tensor_basis_sha256':file_sha(paths['basis']),
            'crosswalk_relative':cross, 'declared_operator':declared,
            'source_files':{k:str(v.relative_to(root)) if v.is_relative_to(root) else str(v)
                            for k,v in paths.items()}, 'gauge_names':names}


def run(root:Path,output_root:Path, *, points:int=3,
        rtol:float=5e-5,atol:float=1e-20,method:str='RK23') -> dict:
    """Solve one rho=+1 ODE; infer rho=-1 by exact real-linearity."""
    from RGE.running.intermediate.GeneralEFT1SevenTermFlow import (
        evolve, c5_projection, gauge_ward_residual, symmetry_residual,
        total_beta, SymmetricCoordinates)
    from RGE.running.intermediate.GeneralEFT1SixTermReferenceCandidate import six_term_beta

    root=Path(root).resolve()
    out=Path(output_root)
    if not out.is_absolute():out=root/out
    require(not out.exists(),'Refusing to overwrite existing output: '+str(out))
    require(2<=points<=11 and 1e-8<=rtol<=1e-2 and 1e-24<=atol<=1e-15,
            'Invalid conditional numerical tolerance/point count')
    require(method in ('RK23','RK45','DOP853'), 'Unsupported Wilson numerical method')
    ctx=load_contract(root)
    basis=ctx['basis']; c0=ctx['c0']; samples=ctx['samples']; x=ctx['scales']
    mu_f=float(x[-1]); mu_s=float(x[0])
    initial_max=max(_ward_max(c0,samples,ctx['gauge_names'],run_index=-1).values())
    require(initial_max<5e-11,'Declared rho=+1 source C fails gauge Ward invariance')
    b0=total_beta(six_term_beta(c0,samples[-1]))
    bm=total_beta(six_term_beta(-c0,samples[-1]))
    anti_res=float(np.max(abs(b0+bm))/max(np.max(abs(b0)),1e-300))
    require(anti_res<1e-11,'Six-term flow violates C -> -C real-linearity')
    init_beta_ward=max(_ward_max(b0,samples,ctx['gauge_names'],run_index=-1).values())
    require(init_beta_ward<1e-7,'Six-term candidate beta fails initial gauge Ward test')
    result=evolve(basis,c0,mu_f,mu_s,x,samples,
                  points=points,rtol=rtol,atol=atol,method=method,
                  beta_mode='fonseca_eq22_six_term')
    end=result['final']
    last_ward=_ward_max(end,samples,ctx['gauge_names'],run_index=0)
    require(max(last_ward.values())<5e-7,
            'Conditional Wilson evolution fails endpoint gauge Ward invariance')
    require(max(result['symmetry_residuals'])<1e-8,
            'Conditional Wilson beta violates independent pair symmetries')
    coords=SymmetricCoordinates(basis)
    require(coords.size==80730 and result['trajectory'].shape==(points,80730),
            'Full Wilson symmetric coordinate inventory not retained')
    ll_f=c5_projection(c0,basis)
    ll_s=c5_projection(end,basis)
    require(max(abs(ll_f).ravel())<1e-25,'Source F threshold unexpectedly contains LLHH')
    from RGE.running.intermediate.GeneralEFT1SevenTermFlow import gauge_ward_residual
    report={
        'format':FORMAT,'model':'A_-4',
        'status':'CONDITIONAL_DECLARED_C0_FULL_EFT1_WILSON_RUNNING_NOT_PHYSICAL',
        'equation':'16pi² dC/dln(mu) = six noncrossed terms of Fonseca Eq.22 with candidate gauge extension',
        'beta_mode':'fonseca_eq22_six_term',
        'full_beta_loop_certified':False,
        'authority':'declared_rho_plus_one_conditional_only',
        'native_matchete_to_declared_weyl_factor':None,
        'rho_plus_and_minus_compatible_with_abstract_native_cc':True,
        'rho_minus_endpoint_equals_negative_rho_plus_endpoint_by_real_linearity':True,
        'complex_multiplier_cannot_be_factorised_from_antilinear_RGE':True,
        'full_complex_Wilson_shape':[45,45,12,12],
        'independent_symmetric_complex_coordinates':int(coords.size),
        'fermion_field_content':'All 45 left-Weyl SM species, no heavy F',
        'scalar_field_content':'H,S1,S2; all 12 real components until mu_S',
        'tree_F_only_boundary_nonzero_full_tensor_entries':int(np.count_nonzero(c0)),
        'C_at_MF_norm_gev_inv':float(np.linalg.norm(c0)),
        'C_at_MS_norm_gev_inv':float(np.linalg.norm(end)),
        'C_relative_norm_change':float((np.linalg.norm(end)-np.linalg.norm(c0))/np.linalg.norm(c0)),
        'full_tensor_relative_change':float(np.linalg.norm(end-c0)/np.linalg.norm(c0)),
        'direct_LLHH_projection_at_MF_not_SMEFT_C5_gev_inv':_encode(ll_f),
        'induced_LLHH_projection_at_MS_not_SMEFT_C5_gev_inv':_encode(ll_s),
        'rho_minus_induced_LLHH_projection_at_MS':_encode(-ll_s),
        'max_abs_induced_LLHH_not_SMEFT_C5_gev_inv':float(np.max(abs(ll_s))),
        'full_F_only_boundary_phase_branch_was_not_native_measured':True,
        'hard_fermion_one_loop_boundary_included':False,
        'scalar_thresholds_integrated_simultaneously':False,
        'physical_scalar_mass_eigenbasis_and_CG_map_certified':False,
        'matched_SMEFT_C5_available':False,
        'neutrino_masses_available':False,
        'EFT_order':'UV -> integrate out ALL fermions -> EFT1(SM+S1+S2) -> integrate out ALL scalars simultaneously -> SMEFT',
        'input_provenance':{k:v for k,v in ctx.items() if k.endswith('_sha256') or k=='source_files'},
        'source_declared_boundary_vs_master_CG_relative_residual':ctx['crosswalk_relative'],
        'initial_branch_real_linearity_residual':anti_res,
        'gauge_ward_initial_max_relative':float(initial_max),
        'gauge_ward_initial_beta_max_relative':float(init_beta_ward),
        'gauge_ward_endpoint_max_relative':last_ward,
        'pair_symmetry_ODE_max_relative':float(max(result['symmetry_residuals'])),
        'fixed_order_full_RGE_relative_difference':result['fixed_order_difference_relative'],
        'ODE_solver':result['solver'],
        'mu_F_gev':mu_f,'mu_S_gev':mu_s,
        'trajectory_scales_gev':result['scales_gev'].tolist(),
        'generated_zero_boundary_independent_coordinates':int(len(result['generated_indices'])),
        'npz_file':'full_conditional_C_trajectory.npz',
    }
    out.mkdir(parents=True)
    np.savez_compressed(out/'full_conditional_C_trajectory.npz',
        indices_1based=result['coordinates'],
        scales_gev=result['scales_gev'],
        independent_C_trajectory_rho_plus_gev_inv=result['trajectory'],
        C_full_rho_plus_MF_gev_inv=c0,C_full_rho_plus_MS_gev_inv=end,
        C_full_rho_minus_MS_gev_inv=-end)
    (out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return report
