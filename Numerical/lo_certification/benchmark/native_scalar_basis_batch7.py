"""Batch 7: direct native-source scalar quadratic audit and UV-weight-basis handoff.

Conditional diagnostics only. The native quartic V=-L is evaluated directly at
H0=v/sqrt(2), without using the symmetrized HHSS tensor to form the Hessian.
The scalar SU(2) generator weight rotation is taken from the original audit.
Original Batch4/5/6 outputs and underlying source files are read-only.
"""
from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any

import numpy as np

from Numerical.lo_certification.benchmark import common_input_batch4 as b4
from Numerical.lo_certification.benchmark import common_point_batch5 as b5
from Numerical.lo_certification.benchmark import native_portal_batch6 as b6

ROOT = Path(__file__).resolve().parents[3]
OWNED = Path('output/general_observable/lo_certification_benchmark')
FORMAT = 't3-a-m4-native-source-scalar-weight-basis-audit-batch7-v1'
STAGES = ('UV', 'F', 'S')
BRANCHES = ('adjoint_alias', 'ordinary_alias')
EXPECTED_HHSS = {2: 'lambdaH1Inv1', 3: 'lambdaT3',
                 4: 'lambdaT3', 5: 'lambdaH2Inv2'}
FIELD_DIMS = {'H': 2, 'NewScalar1': 1, 'NewScalar2': 3}


class Batch7Error(ValueError):
    pass


def require(ok: Any, why: str) -> None:
    if not bool(ok):
        raise Batch7Error(why)


def read(path: Path) -> dict:
    require(path.is_file(), f'Missing input: {path}')
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), 'Expected JSON object: '+str(path))
    return obj


def read_real_hessian(path: Path, label: str) -> np.ndarray:
    """Read Batch 6's bare 8x8 JSON matrix (not its object-format reports).

    Batch 6 writes `h.tolist()` and `delta.tolist()` directly via json.dumps.
    Do not send those arrays through `read()`, which deliberately only accepts
    JSON objects for source documents and summary reports.
    """
    require(path.is_file(), f'Missing {label}: {path}')
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, list) and len(obj) == 8 and
            all(isinstance(row, list) and len(row) == 8 for row in obj),
            f'{label}: expected bare 8x8 JSON array from Batch6: {path}')
    require(all(type(x) in (int, float) and math.isfinite(float(x))
                for row in obj for x in row),
            f'{label}: nonfinite or non-real Hessian entry: {path}')
    matrix = np.asarray(obj, dtype=float)
    require(matrix.shape == (8, 8) and
            np.linalg.norm(matrix - matrix.T) <=
            1e-10 * max(np.linalg.norm(matrix), 1.0),
            f'{label}: nonsymmetric Hessian: {path}')
    return matrix


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def decode(obj: dict, shape: tuple[int, ...], label: str) -> np.ndarray:
    require(isinstance(obj, dict) and {'real', 'imag'} <= set(obj),
            f'{label}: explicit real/imag JSON required')
    value = np.asarray(obj['real'], float)+1j*np.asarray(obj['imag'], float)
    require(value.shape == shape and np.isfinite(value).all(),
            f'{label}: wrong/nonfinite shape')
    return value


def encode(z: np.ndarray) -> dict:
    a = np.asarray(z, complex)
    require(np.isfinite(a).all(), 'Cannot encode nonfinite matrix')
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def neutral_vev_hessian(seed: dict, couplings: dict[str, Any], vev: float) -> tuple[np.ndarray, dict]:
    """Differentiate the original native Matchete L4 terms after H0=v/sqrt2.

    Direct source contraction: select SU2 components of each Field[], evaluate
    each registered native CG at repeated Matchete indices, insert two neutral
    Higgs VEVs and expand only the two remaining scalar legs. This bypasses
    the quartic-to-lambda_ABab factorial conversion used in Batch6.

    Real heavy axes: (S1_R,S1_I,S2_1R,S2_1I,S2_2R,S2_2I,S2_3R,S2_3I).
    """
    require(math.isfinite(vev) and vev > 0, 'Invalid neutral Higgs VEV')
    require(seed.get('IntegratedFields') == ['F'] and
            seed.get('ActiveHeavyFields') == ['S1','S2'],
            'Native input is not the all-fermion-first EFT1 scalar source')
    terms = seed.get('ScalarQuarticTerms')
    require(isinstance(terms,list) and
            len(terms) == seed.get('ScalarQuarticTermCount'),
            'Unverified source scalar inventory')
    hhss = [t for t in terms if t.get('HiggsFieldCount') == 2 and
            t.get('Scalar1FieldCount',0)+t.get('Scalar2FieldCount',0) == 2]
    require({t['Index'] for t in hhss} == set(EXPECTED_HHSS) and len(hhss)==4,
            'Native HHSS term inventory changed; no guessed extra channels')
    require(set(couplings) == set(EXPECTED_HHSS.values()),
            'Native HHSS coupling values incomplete/extra')
    from RGE.running.intermediate import ScalarOnlyTensorAdapters as source
    registry = source.load_cg_registry(seed)
    htot = np.zeros((8,8), complex)
    per_term: dict[str,dict] = {}
    for term in hhss:
        idx = int(term['Index'])
        raw = term['TermInputForm']
        legs = source._find_scalar_legs(raw)
        calls = source._find_cg_calls(raw)
        name, conjugated = source._find_coupling(raw)
        require(name == EXPECTED_HHSS[idx],
                'Native term coupling or term ordering changed: '+str(idx))
        require(len(legs) == 4 and sum(leg.name == 'H' for leg in legs) == 2,
                'Wrong source H/H/S/S field assignment')
        require(all(leg.name in FIELD_DIMS for leg in legs),
                'Unknown native scalar field')
        prefactor = -complex(source.sp.N(source._numeric_prefactor(raw), 32))
        z = b6.value(couplings[name], name)
        if conjugated:
            z = z.conjugate()
        h = np.zeros((8,8),complex)
        ranges = [range(1, FIELD_DIMS[leg.name]+1) for leg in legs]
        for components in product(*ranges):
            assignment: dict[tuple[str,str],int] = {}
            consistent = True
            heavy = []
            amplitude = prefactor*z
            for leg,comp in zip(legs,components):
                if leg.dummy is not None:
                    key = (leg.dummy,leg.rep)
                    if key in assignment and assignment[key] != comp:
                        consistent = False
                        break
                    assignment[key] = comp
                if leg.name == 'H':
                    if comp != 2:
                        consistent = False
                        break
                    amplitude *= vev/math.sqrt(2)
                else:
                    vector = np.zeros(8,complex)
                    base = 0 if leg.name == 'NewScalar1' else 2
                    pos = base+2*(comp-1)
                    vector[pos] = 1/math.sqrt(2)
                    vector[pos+1] = (-1j if leg.conjugated else 1j)/math.sqrt(2)
                    heavy.append(vector)
            if not consistent:
                continue
            require(len(heavy) == 2, 'Missing heavy scalar legs')
            for call in calls:
                amplitude *= complex(source.sp.N(source._cg_value(call,registry,assignment),32))
            h += amplitude*(np.outer(heavy[0],heavy[1]) +
                            np.outer(heavy[1],heavy[0]))
        require(np.isfinite(h).all(), 'Native term Hessian nonfinite')
        per_term[str(idx)] = {'native_coupling':name,
                              'source_term_sha256':sha256(raw.encode('utf8')).hexdigest(),
                              'hessian_frobenius_gev2':float(np.linalg.norm(h)),
                              'matrix_imag_max_gev2':float(np.max(abs(h.imag)))}
        htot += h
    require(np.max(abs(htot.imag)) <= 3e-10*max(np.max(abs(htot)),1e-30),
            'Source V=-L is not real after native HC pair summation')
    require(np.linalg.norm(htot-htot.T) <= 3e-11*max(np.linalg.norm(htot),1e-30),
            'Directly differentiated native Hessian not symmetric')
    return htot.real, per_term


def weight_rotation(triplet_columns: np.ndarray) -> np.ndarray:
    """Native scalar triplet field = U_triplet * scalar weight eigenfield.

    Same passive replacement used by BuildT3DirectOneLoopInputs._tensor:
    plain fields use U columns; conjugate fields use U* columns.
    """
    t = np.asarray(triplet_columns,complex)
    require(t.shape == (3,3) and np.isfinite(t).all(), 'Invalid native triplet weight rotation')
    require(np.linalg.norm(t.conj().T@t-np.eye(3)) < 1e-11,
            'Native triplet generator eigenvectors are not unitary')
    u = np.eye(4,dtype=complex)
    u[1:,1:] = t
    return u


def rotate_native_complex_matrices(a: np.ndarray, b: np.ndarray,
                                   u: np.ndarray) -> tuple[np.ndarray,np.ndarray]:
    a,b,u = (np.asarray(x,complex) for x in (a,b,u))
    require(a.shape == b.shape == u.shape == (4,4), 'Scalar A/B/rotation shape changed')
    require(np.linalg.norm(a-a.conj().T) < 1e-11*max(np.linalg.norm(a),1.), 'Native A not Hermitian')
    require(np.linalg.norm(b-b.T) < 1e-11*max(np.linalg.norm(b),1.), 'Native B not symmetric')
    require(np.linalg.norm(u.conj().T@u-np.eye(4)) < 1e-11, 'Nonunitary source rotation')
    aa = u.conj().T@a@u
    bb = u.T@b@u
    require(np.linalg.norm(aa-aa.conj().T) < 1e-11*max(np.linalg.norm(aa),1.),
            'Weight-basis A not Hermitian')
    require(np.linalg.norm(bb-bb.T) < 1e-11*max(np.linalg.norm(bb),1.),
            'Weight-basis B not symmetric')
    require(np.linalg.norm(u@aa@u.conj().T-a) < 2e-11*max(np.linalg.norm(a),1.),
            'Reverse native mass rotation failed')
    require(np.linalg.norm(u.conj()@bb@u.conj().T-b) < 2e-11*max(np.linalg.norm(b),1.),
            'Reverse native holomorphic mass rotation failed')
    return aa,bb


def matrix_relative(a: np.ndarray, b: np.ndarray, floor: float=1e-100) -> float:
    a,b = np.asarray(a,complex),np.asarray(b,complex)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(),
            'Relative matrix input mismatch')
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),floor))


def preflight(root: Path, batch4_id: str, batch5_id: str, batch6_id: str) -> dict:
    root=root.resolve()
    for run_id in (batch4_id,batch5_id,batch6_id):
        require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,50}',run_id) is not None,
                'Unsafe source run ID')
    b6.preflight(root,batch4_id,batch5_id)
    f4=root/OWNED/('batch4_common_input_'+batch4_id)
    f5=root/OWNED/('batch5_uv_eft_'+batch5_id)
    f6=root/OWNED/('batch6_portal_'+batch6_id)
    o6=read(f6/'summary.json')
    require(o6.get('format') == b6.FORMAT and
            o6.get('status') == 'NATIVE_HHSS_HESSIAN_CANDIDATE_ALIASES_UNCERTIFIED' and
            o6.get('physical_C5_certified') is False,
            'Completed uncertified Batch6 source missing')
    require(o6.get('batch4_summary_SHA256') == digest(f4/'summary.json') and
            o6.get('batch5_summary_SHA256') == digest(f5/'summary.json'),
            'Batch6 differs from source Batch4/5')
    require(o6.get('source_quartic_seed_SHA256') == digest(root/b4.QUARTIC_SEED),
            'Native source quartic seed has changed since Batch6')
    require(set(o6.get('scenarios',{})) == set(STAGES) and
            all(set(o6['scenarios'][s])==set(BRANCHES) for s in STAGES),
            'Batch6 did not finish all six branches')
    original=read(root/b5.INPUT_UV)
    _case, us, _paths, _cov = b5.check_audits(root,original)
    u = weight_rotation(np.asarray(us[3],complex))
    return {'model':'A_-4','status':'PREFLIGHT_OK_NO_OUTPUT_WRITTEN',
            'stages':list(STAGES),'branches':list(BRANCHES),
            'native_triplet_rotation_unitarity_residual':float(np.linalg.norm(u.conj().T@u-np.eye(4))),
            'source_fingerprint':digest(root/b4.QUARTIC_SEED),
            'source_predecessors':{'batch4':str(f4),'batch5':str(f5),'batch6':str(f6)},
            'scalar_basis_rotation_required_for_direct_UV_evaluator':True,
            'physical_certification':False,'output':None}


def calculate(root: Path, batch4_id: str, batch5_id: str, batch6_id: str,
              run_id: str) -> dict:
    root=root.resolve()
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,50}',run_id) is not None,
            'Unsafe new run ID')
    target=root/OWNED/('batch7_native_basis_'+run_id)
    require(not target.exists(),'Refusing to overwrite existing Batch7 output: '+str(target))
    preflight(root,batch4_id,batch5_id,batch6_id)
    from Numerical.lo_certification.direct_uv.native_uv_loop_reference import (
        independent_reference, decode as uv_decode)
    from Numerical.direct_one_loop.core import evaluate

    f4=root/OWNED/('batch4_common_input_'+batch4_id)
    f5=root/OWNED/('batch5_uv_eft_'+batch5_id)
    f6=root/OWNED/('batch6_portal_'+batch6_id)
    o4,o5,o6 = (read(x/'summary.json') for x in (f4,f5,f6))
    frozen=read(root/b4.BENCH)
    running=read(f4/'fitted_running.json')
    seed=read(root/b4.QUARTIC_SEED)
    original=read(root/b5.INPUT_UV)
    _case, us, _paths, _cov=b5.check_audits(root,original)
    u=weight_rotation(np.asarray(us[3],complex))
    vev=float(frozen['numerical']['vev_gev'])
    eft=decode(o4['conditional_mnu_eV_per_N'],(3,3),'Batch4 EFT mass')
    require(np.linalg.norm(eft)>0, 'Cannot compare zero EFT mass')

    output = {}
    records = {}
    for stage in STAGES:
        state=b6.scalar_stage(frozen,running,stage)
        m1=b6.value(state['mS1Sq'],'mS1Sq')
        m2=b6.value(state['mS2Sq'],'mS2Sq')
        require(m1.imag==0 and m2.imag==0 and min(m1.real,m2.real)>0,
                'Invalid real unbroken scalar mass values')
        bare=np.diag([m1.real]*2+[m2.real]*6)
        oldinput=read(f5/(stage+'_direct_input.json'))
        old_A=uv_decode(oldinput['scalar_A_gev2'],'Batch5 A',(4,4))
        old_B=uv_decode(oldinput['scalar_B_gev2'],'Batch5 B',(4,4))
        by_branch={}
        for branch in BRANCHES:
            couplings=b6.candidate_couplings(state,branch)
            shift,terms=neutral_vev_hessian(seed,couplings,vev)
            h_native=bare+shift
            h6=read_real_hessian(f6/(stage+'_'+branch)/'Hessian_real_interleaved_GeV2.json',
                                 'Batch6 full scalar Hessian')
            d6=read_real_hessian(f6/(stage+'_'+branch)/'delta_Hessian_from_HHSS_GeV2.json',
                                 'Batch6 HHSS Hessian shift')
            require(h6.shape==d6.shape==(8,8),'Batch6 Hessian dimension changed')
            source_vs_b6=matrix_relative(shift,d6)
            # An independently found mismatch is an audit RESULT, not a value to overwrite.
            # Retain native source directly as authority for this branch while flagging
            # any disagreement with the previous quartic/factorial implementation.
            hessian_agrees_with_batch6 = (source_vs_b6 < 4e-8 and
                                            matrix_relative(h_native,h6) < 3e-11)
            require(np.linalg.eigvalsh(h_native).min()>0,
                    'Native source gives nonpositive scalar mass eigenvalues')
            a_native,b_native=b6.hessian_to_complex(h_native)
            a_weight,b_weight=rotate_native_complex_matrices(a_native,b_native,u)
            original6=read(f6/(stage+'_'+branch)/'input.json')
            a_native_vs_batch6 = matrix_relative(a_native,uv_decode(original6['scalar_A_gev2'],'Batch6 native A',(4,4)))
            b_native_vs_batch6 = matrix_relative(b_native,uv_decode(original6['scalar_B_gev2'],'Batch6 native B',(4,4)))
            hessian_agrees_with_batch6 = (hessian_agrees_with_batch6 and
                a_native_vs_batch6 < 4e-11 and b_native_vs_batch6 < 4e-11)
            aligned=deepcopy(original6)
            aligned['scalar_A_gev2']=encode(a_weight)
            aligned['scalar_B_gev2']=encode(b_weight)
            aligned['provenance'].update({
                'batch7_weight_basis_rotation_applied':True,
                'batch7_original_native_Lagrangian_differentiated':True,
                'batch7_basis_rotation_source':'matchete_generator_analysis.json U_columns_are_weight_eigenvectors',
                'batch7_reference_quartic_alias_hypothesis':branch,
                'batch7_spinor_normalization_certified':False,
                'batch7_complete_UV_EFT_matching_certified':False,
                'batch7_native_HHSS_seed_SHA256':digest(root/b4.QUARTIC_SEED),
            })
            uv=independent_reference(aligned,check_rephase=False)
            prod=evaluate(aligned,retain_individual_terms=False)
            nu=uv_decode(uv['complex_symmetric_matrix_gev'],'Batch7 UV',(3,3))
            nu_prod=uv_decode(prod['total']['complex_symmetric_matrix_gev'],'production UV',(3,3))
            independent_delta=matrix_relative(nu,nu_prod)
            require(independent_delta < 4e-7, 'Independent direct-UV evaluation disagrees with production')
            require(prod['diagnostics']['uv_finite'] is True and
                    uv['diagnostics']['uv_pole_relative_to_natural_amplitude']<1e-9,
                    'UV pole cancellation failed for rotated scalar input')
            m6=uv_decode(read(f6/(stage+'_'+branch)/'independent.json')['complex_symmetric_matrix_gev'],
                         'Batch6 UV mass',(3,3))
            m5=uv_decode(read(f5/(stage+'_direct_independent.json'))['complex_symmetric_matrix_gev'],
                         'Batch5 UV mass',(3,3))
            residual=b5.comparison(nu,eft)
            prior=b5.comparison(m6,eft)
            projected=matrix_relative(a_weight,old_A)
            record={**residual,
                'stage':stage,'branch':branch,'higgs_vev_gev':vev,
                'source_direct_hessian_vs_Batch6_quartic_relative':source_vs_b6,
                'source_native_scalar_A_vs_Batch6_relative':a_native_vs_batch6,
                'source_native_scalar_B_vs_Batch6_relative':b_native_vs_batch6,
                'source_direct_Hessian_agrees_with_Batch6':hessian_agrees_with_batch6,
                'source_native_vs_legacy_Batch5_weight_A_relative':projected,
                'source_weight_B_frobenius_gev2':float(np.linalg.norm(b_weight)),
                'source_weight_vs_legacy_Batch5_B_relative':matrix_relative(b_weight,old_B),
                'scalar_source_to_weight_rotation_U':encode(u),
                'Batch6_mass_matrix_relative_difference':matrix_relative(nu,m6),
                'Batch5_mass_matrix_relative_difference':matrix_relative(nu,m5),
                'Batch6_matrix_shape_residual':prior['residual_after_one_complex_multiplier'],
                'UV_masses_eV_descending':np.linalg.svd(nu*1e9,compute_uv=False).tolist(),
                'Batch5_masses_eV_descending':np.linalg.svd(m5*1e9,compute_uv=False).tolist(),
                'Batch6_masses_eV_descending':np.linalg.svd(m6*1e9,compute_uv=False).tolist(),
                'native_source_term_Hessians':terms,
                'independent_vs_production_UV_relative':independent_delta,
                'uv_pole_relative':uv['diagnostics']['uv_pole_relative_to_natural_amplitude'],
                'physical_C5_certified':False,
            }
            by_branch[branch]=record
            records[(stage,branch)]={'aligned_UV_input.json':aligned,
                                   'independent_direct_UV.json':uv,
                                   'production_direct_UV.json':prod,
                                   'native_Hessian_interleaved_GeV2.json':h_native.tolist(),
                                   'native_source_HHSS_Hessian_shift_GeV2.json':shift.tolist(),
                                   'weight_basis_scalar_A_gev2.json':encode(a_weight),
                                   'weight_basis_scalar_B_gev2.json':encode(b_weight)}
        output[stage]=by_branch

    all_hessian_agree=all(output[s][b]['source_direct_Hessian_agrees_with_Batch6']
                          for s in STAGES for b in BRANCHES)
    summary={
        'format':FORMAT,'model':'A_-4',
        'status':('NATIVE_SOURCE_HESSIAN_WEIGHT_BASIS_RECONCILED_CONDITIONAL'
                  if all_hessian_agree else
                  'NATIVE_SOURCE_HESSIAN_DISAGREEMENT_EXPLICIT_NOT_CERTIFIED'),
        'all_six_native_source_Hessians_agree_with_Batch6':all_hessian_agree,
        'source':{'batch4_summary_SHA256':digest(f4/'summary.json'),
                  'batch5_summary_SHA256':digest(f5/'summary.json'),
                  'batch6_summary_SHA256':digest(f6/'summary.json'),
                  'native_HHSS_seed_SHA256':digest(root/b4.QUARTIC_SEED)},
        'native_source_evaluated_without_lambda_tensor_factorial':True,
        'native_triplet_weight_rotation_from_existing_SU2_generator_audit':True,
        'H0_vev':'H0=v/sqrt(2); source V4=-L4',
        'matching_sequence':'UV -> integrate out ALL heavy F at muF -> EFT1(S1,S2) -> integrate out BOTH scalars together at muS -> SMEFT',
        'physics_note':'UV is a broken-phase loop; EFT is an unbroken fixed-order LLSS-induced conditional result. No full one-loop equivalence asserted.',
        'branch_aliases_not_certified':dict(b6.BRANCHES),
        'scale_status':'UV scalar stage is frozen muUV, F is muF, S is EFT1 muS hybrid. Fermion parameters are muF for all.',
        'scenarios':output,
        'gammaCC_to_master_Weyl_signed_factor':None,
        'source_rgbeta_quartic_aliases_certified':False,
        'fermion_hard_one_loop_matching_complete':False,
        'physical_C5_certified':False,
        'complete_neutrino_mass_certified':False,
    }
    target.parent.mkdir(parents=True,exist_ok=True)
    staging=Path(tempfile.mkdtemp(prefix='.batch7_native_stage_',dir=target.parent))
    try:
        (staging/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf8')
        report=['# Batch 7 — A_-4 native Lagrangian Hessian / SU(2) scalar weight basis',
                '', '**Conditional diagnostic, not a physical UV-to-SMEFT matching certificate.**','',
                'Original Matchete HHSS terms are contracted directly at H0=v/sqrt(2), then differentiated analytically in 8 real scalar coordinates. The direct-source Hessian is checked against the symmetrized Batch6 source tensor; any disagreement is recorded explicitly, not silently corrected. The triplet unitary from the existing generator audit rotates scalar matrices to the basis used for the direct-UV Yukawa vertices.',
                '', '| Scalar scale | Alias branch | UV masses (eV; descending) | Complex-scale shape residual | Batch5 weight-A relative | Batch6 UV-matrix relative |',
                '|---|---|---|---:|---:|---:|']
        for stage in STAGES:
            for branch in BRANCHES:
                row=output[stage][branch]
                report.append('| '+stage+' | '+branch+' | '+', '.join(f'{v:.9g}' for v in row['UV_masses_eV_descending'])+
                              f" | {row['residual_after_one_complex_multiplier']:.8g} | {row['source_native_vs_legacy_Batch5_weight_A_relative']:.8g} | {row['Batch6_mass_matrix_relative_difference']:.8g} |")
        report += ['', 'Absolute source-to-Weyl GammaCC normalization is unknown.',
                   'Matchete-to-RGBeta quartic aliases are candidate-only, and the full hard-fermion matching is absent.',
                   'The two scalar multiplets remain active together in EFT1 and are integrated out together at one muS.', '']
        (staging/'report.md').write_text('\n'.join(report),encoding='utf8')
        for (stage,branch),outputs in records.items():
            folder=staging/(stage+'_'+branch)
            folder.mkdir()
            for name,doc in outputs.items():
                (folder/name).write_text(json.dumps(doc,indent=2,allow_nan=False)+'\n',encoding='utf8')
        require(not target.exists(),'Output became occupied while Batch7 was computing')
        staging.rename(target)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return {'model':'A_-4','status':summary['status'],'output':str(target),
            'all_source_Hessians_reproduce_Batch6':all_hessian_agree,
            'scenarios':{s:{b:{'residual':round(output[s][b]['residual_after_one_complex_multiplier'],12),
                                  'UV_masses_ev':output[s][b]['UV_masses_eV_descending'],
                                  'scalar_A_difference_from_Batch5':output[s][b]['source_native_vs_legacy_Batch5_weight_A_relative']}
                           for b in BRANCHES} for s in STAGES},
            'physical_C5_certified':False}
