"""Batch 9 v11: compare F-first and scalar-only LLHH poles in ONE Warsaw normalization.

Read-only. Reuse the archived B/C F-first one-loop LLHH source (68 pole,
68 log, 68 finite terms), and the completed Batch 8 scalar one-loop C5.
Evaluate both in the same *frozen* muF couplings and 3-diagonal-MF lift.
Do NOT infer a missing inherited pole as minus the scalar pole. The archived
native F expression, with the project's neutral SU(2)/Warsaw projector,
must independently produce it. A mismatched/missing archive blocks. The already committed spinor-certification
F-one-loop raw source may be used instead, but only after its native JSON,
inspection SHA256, and source inventory have been checked.

Batch 8 stores the symmetrized ORDERED Matchete A5 coefficient, whereas the
archived F-first projector already returns the half-normalized symmetric
C5=A5+A5.T. v11 verifies the source definition in RunMatching.wl and
converts the saved scalar numeric contribution by a factor of TWO before
comparing pole residues or adding finite matrices. This does not modify old
Batch 8 artifacts and never assumes the pole cancellation succeeds.

A passing matrix equality is conditional fixed-order cancellation of these
specific raw sectors, NOT complete UV/EFT renormalization, EFT1 RGE, or
physical C5. No Matchete or fermion matching is invoked.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import sys
from typing import Any

import numpy as np
import sympy as sp

from Numerical.running.T3MixedBatch6Numerics import (
    encode_complex, evaluate_hard, metrics, numerical_inputs,
)
from scripts import RunT3MixedFixedOrderBatch6 as batch6
from scripts import RunT3MixedScalarHardBatch6 as native

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT/'output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop'
# Existing repository native F-first one-loop export, independently archived.
# Do not confuse it with an EFT1 tree export or scalar-hard one-loop input.
CERTIFIED_NATIVE_F_ARCHIVE = ROOT/('output/general_observable/lo_certification_spinor/'
                                   'native_hard_e515d82d118284ba_native01')
# Exact committed native evidence fingerprints (verified against inspection.json).
KNOWN_NATIVE_F_SHA256 = {
    'B_-1':'dceaa8c40af8eb2885715a2fc9719188c22220b020ab75c3487e641203337b0c',
    'C_-1':'c043dac85e99f99ec88ebc474137c6ff95eca923621c25ea6868999851f24e32',
}
OUTPUT = ROOT/'output/general_sequential_pilot/mixed_inherited_pole_batch9'
STATUS = 'CommonCouplingFInheritedAndScalarRawPolesCancelled_Conditional'
FORMAT = 't3-mixed-f-inherited-scalar-common-pole-v1'
MODELS = ('B_-1','C_-1')
QSET = {'lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1'}


class PoleGateError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PoleGateError(message)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    require(path.is_file(), 'Missing source JSON: '+str(path))
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj,dict), 'Not a JSON object: '+str(path))
    return obj


def verify_scalar_extraction_convention(root: Path) -> dict:
    """Source-gate the factor converting Batch 8's ordered A5 into C5.

    `ExtractWeinbergCoefficient` returns its ordered `Coefficient` and reports
    separately that C5=A5+A5.T for the 1/2-normalized Weinberg operator. The
    archived Batch 8 numerical code symmetrized its ordered A5 instead of
    applying this full flavour conversion. Refuse to apply factor 2 if the
    actual on-disk Matchete wrapper does not have the checked convention.
    """
    path=root/'Lagrangian'/'RunMatching.wl'
    require(path.is_file(), 'Cannot establish scalar A5/C5 convention without '+str(path))
    content=path.read_text(encoding='utf-8-sig')
    anchor='ExtractWeinbergCoefficient[eft_] := Module['
    require(content.count(anchor)==1,
            'Matchete ordered Weinberg extractor missing or ambiguous')
    definition=content[content.index(anchor):]
    doc=content[:content.index(anchor)]
    require("Matchete's ordered coefficient A5" in doc and
            'C5^(ij) = A5^(ij) + A5^(ji)' in doc and
            '"Coefficient" -> coefficient' in definition and
            '"CanonicalCoefficient" -> canonicalCoefficient' in definition and
            bool(re.search(r'holomorphicTerms\s*=\s*Select\[terms,\s*!FreeQ\[#,\s*Proj\[-1\]\]\s*&\s*\];', definition)) and
            'coefficient = CompactWeinbergCoefficient[holomorphicTerms]' in definition,
            'Unsupported ExtractWeinbergCoefficient A5 versus C5 convention; no factor assumed')
    return {
        'source_path':str(path.resolve()),
        'source_sha256':digest(path),
        'ordered_matchete_extractor_returns_A5':True,
        'half_normalized_warsaw_definition':'C5 = A5 + Transpose[A5]',
        'batch8_old_numeric_symmetrization':'Asym = (A5 + Transpose[A5])/2',
        'factor_from_saved_batch8_Asym_to_C5':2,
        'historical_batch8_results_not_rewritten':True,
    }


def decode_matrix(obj: Any, label: str) -> np.ndarray:
    require(isinstance(obj,dict) and 'real' in obj and 'imag' in obj,
            label+': missing complex matrix')
    a=np.asarray(obj['real'],dtype=float)+1j*np.asarray(obj['imag'],dtype=float)
    require(a.shape==(3,3) and np.isfinite(a).all(),label+': bad matrix shape/value')
    require(np.linalg.norm(a-a.T)<2e-10*max(np.linalg.norm(a),1e-40),
            label+': Majorana C5 not symmetric')
    return a


def relerr(a: np.ndarray, b: np.ndarray, floor: float=1e-40) -> float:
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),floor))


def archived_source_path(archive_root: Path, model: str) -> Path:
    require(model in MODELS,'Unsupported model')
    stem=native.tag(model)
    return archive_root/f'{stem}_final_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt'



def recovered_native_f_archive_paths(root: Path, model: str) -> dict[str, Path]:
    """The repository's pre-existing native F-first one-loop evidence (read-only)."""
    require(model in MODELS, 'Unsupported model')
    stem = native.tag(model)
    folder = root / stem
    prefix = f'{stem}_native_one_loop'
    return {
        'source': folder / (prefix + '_raw_hard_llhh_sector.inputform.txt'),
        'report': folder / (prefix + '.json'),
        'inspection': folder / 'inspection.json',
    }


def inspect_recovered_native_archive(root: Path, model: str) -> dict:
    """Require real F-one-loop Matchete provenance and raw-file SHA256.

    The archived source contains symbolic U11..U22 spurions. Its old
    numerical U values are provenance only; current Batch8 U is substituted
    by evaluate_polynomials, not copied from the archived report.
    """
    paths = recovered_native_f_archive_paths(root, model)
    for label, path in paths.items():
        require(path.is_file(), model+': native F-one-loop recovery missing '+label+
                ': '+str(path))
    inspection = read_json(paths['inspection'])
    report = read_json(paths['report'])
    stem = native.tag(model)
    require(inspection.get('status') == 'MATCHeTE_RAW_MIXED_LOOP_LLHH_EXTRACTED_NOT_MSbar'
            and inspection.get('model') == model
            and inspection.get('driver_type') == 'mixed_BC'
            and inspection.get('native_exit_code') == 0
            and inspection.get('one_loop_F_native_run_returned') is True
            and inspection.get('hard_sector_renormalized') is False
            and inspection.get('explicit_hard_LLHH_candidate_terms') == 204
            and inspection.get('physical_wilson_coefficient') is None
            and inspection.get('raw_hard_llhh_file') == paths['source'].name,
            model+': unsupported native F-one-loop inspection certificate')
    require(inspection.get('raw_hard_llhh_sha256') == digest(paths['source'])
            and digest(paths['source']) == KNOWN_NATIVE_F_SHA256[model]
            and inspection.get('native_report_sha256') == digest(paths['report']),
            model+': native F-one-loop archive source/report SHA256 mismatch')
    require(report.get('format') == 't3-mass-eigen-fermion-first-eft1-oneloop-v1'
            and report.get('status') == 'FermionFirstOneLoopRawSectorsExtractedNotRenormalized'
            and report.get('model') == model
            and report.get('matching_approach') == 'FermionFirstTreeAndOneLoopScalarsRemainLight'
            and report.get('matching_EFT_order') == 5
            and report.get('matching_loop_order') == 1
            and report.get('source_quartic_count') == 17
            and report.get('source_yukawa_count') == 2
            and report.get('integrated_heavy_fields') == ['NewFermion']
            and report.get('scalar_light_fields') == ['NewScalarX','NewScalarY']
            and report.get('llss_term_count') == (20 if model == 'B_-1' else 24)
            and report.get('raw_hard_llhh_candidate_term_count') == 204
            and report.get('raw_hard_llhh_file') == paths['source'].name
            and report.get('tree_matching_completed') is True
            and report.get('one_loop_matching_completed') is True
            and report.get('loop_difference_fermion_absent') is True
            and report.get('raw_hard_llhh_scheme_renormalized') is False
            and report.get('raw_hard_llhh_pole_rge_certified') is False
            and report.get('physical_C5_computed') is False,
            model+': unsupported native F-first loop report')
    # The full native report is not a physical renormalization certificate.
    return {'archive': paths['source'],
            'provenance': 'ExistingNativeSpinorCertificationFOneLoop',
            'inspection_sha256': digest(paths['inspection']),
            'native_report_sha256': digest(paths['report']),
            'source_sha256': digest(paths['source']),
            'current_U_substituted_for_symbolic_native_spurions': True,
            'native_archive_U_is_not_current_Batch8_U': None,
            'source_one_loop_pole_rge_certified': False}


def resolve_f_first_archive(archive_root: Path, model: str) -> dict:
    """Original F-first source first; verified existing repository fallback.

    An explicit non-default --archive-root is never silently overridden.
    """
    primary = archived_source_path(archive_root, model)
    if primary.is_file():
        return {'archive':primary,
                'provenance':'OriginalFermionFirstOneLoopArchive',
                'source_sha256':digest(primary),
                'source_one_loop_pole_rge_certified':False}
    if archive_root.resolve() != ARCHIVE.resolve():
        require(False, 'Archived F-FIRST ONE-LOOP LLHH source missing at explicitly selected '
                '--archive-root: '+str(primary))
    # This is an independent successful native one-loop F export committed
    # under lo_certification_spinor, not the historical general_t3 path.
    candidate = recovered_native_f_archive_paths(CERTIFIED_NATIVE_F_ARCHIVE, model)
    require(candidate['source'].is_file(),
            'No F-first ONE-LOOP LLHH source at '+str(primary)+' or '+
            str(candidate['source'])+'; TREE EOMEFT is insufficient')
    return inspect_recovered_native_archive(CERTIFIED_NATIVE_F_ARCHIVE, model)

def load_native_f_polynomials(root: Path, source: Path, model: str) -> dict:
    """Use the verified legacy parser, not a new CG/epsilon convention."""
    require(source.is_file(), 'Missing archived F-FIRST ONE-LOOP source: '+str(source)
            +'; tree EOMEFT alone cannot supply this inheritance')
    import importlib.util
    script=root/'scripts/BuildT3BCFHardBoundaryInterface.experimental.py'
    require(script.is_file(), 'Missing repository F-boundary parser: '+str(script))
    spec=importlib.util.spec_from_file_location('t3bc_f_first_archived_batch9',script)
    require(spec is not None and spec.loader is not None,
            'Cannot load repository F-boundary parser')
    v13=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(v13)
    v13.check_project_imports(root)
    v12=v13.import_v12(root)
    raw=source.read_text(encoding='utf-8-sig')
    sectors=v13.split_sectors(raw)  # fail closed: exactly 68/68/68
    v13.exact_pole_log_check(sectors)  # exact opaque-factor check
    polys={}
    for key in ('pole','log','constant'):
        ledger, counts=v13.make_ledger(sectors[key],model)
        require(counts=={'PL':34,'PR':34},model+': F boundary chirality inventory')
        pl=v12.polynomial_coefficients(ledger,'PL')
        pr=v12.polynomial_coefficients(ledger,'PR')
        require(len(pl)==3 and all(sp.expand(
            v12.polynomial_conjugate(x)-y)==0 for x,y in zip(pr,pl)),
            model+': F inherited PL/PR conjugation failed')
        polys[key]=tuple(pl)
    require(all(sp.expand(p-l)==0 for p,l in zip(polys['pole'],polys['log'])),
            model+': F pole/log coefficients differ after SU(2) reduction')
    return {'poly':polys,'archive_path':str(source.resolve()),'archive_sha256':digest(source),
            'terms':{'pole':68,'log':68,'constant':68},
            'F_finite_constant_equals_pole':all(sp.expand(c-p)==0
               for c,p in zip(polys['constant'],polys['pole']))}


def evaluate_polynomials(polys: dict, inp: dict) -> dict[str,np.ndarray]:
    """Freeze all native mixing spurions/quartics at Batch 8 muF inputs."""
    require(set(inp['couplings'])==QSET,
            'All three native quartics required; conditional RGBeta map must be explicit')
    substitutions={}
    for i in range(2):
        for j in range(2):
            name=f'T3MEU{i+1}{j+1}'
            z=complex(inp['U'][i,j])
            substitutions[sp.Symbol(name)]=z
            substitutions[sp.Symbol('bar_'+name)]=z.conjugate()
    for name,z in inp['couplings'].items():
        z=complex(z)
        substitutions[sp.Symbol(name)]=z
        substitutions[sp.Symbol('bar_'+name)]=z.conjugate()
    result={}
    for sector,vector in polys.items():
        require(len(vector)==3,sector+': expected three source Yukawa combinations')
        values=[]
        for p in vector:
            require(not p.free_symbols.difference(substitutions),
                    'Unbound native F inherited spurions: '+str(p.free_symbols-set(substitutions)))
            z=complex(sp.N(p.subs(substitutions),18))
            require(math.isfinite(z.real) and math.isfinite(z.imag),
                    'Nonfinite projected F source polynomial')
            values.append(z)
        result[sector]=np.asarray(values,dtype=complex)
    require(relerr(result['pole'],result['log'])<1e-12,
            'Numerical F pole and log polynomials disagree')
    return result


def coefficient_matrix(vec: np.ndarray) -> np.ndarray:
    """The established v12/v13 convention: A+A.T, *not* A+A.conj().T."""
    require(vec.shape==(3,),'Expected three source Yukawa-pair coefficients')
    return np.asarray([[2*vec[0],vec[1]],[vec[1],2*vec[2]]],complex)


def f_inherited_matrices(inputs: dict, vectors: dict[str,np.ndarray],
                         mu: float) -> tuple[np.ndarray,np.ndarray,list[dict]]:
    """Raw F-first one-loop pole and finite plus MF-log at a common scale.

    Every F eigenmass enters an independently evaluated 1/M_r term. This is
    still an inferred lift of the original ONE unindexed MF source; not a
    newly native matched three-MF one-loop calculation.
    """
    require(mu>0 and math.isfinite(mu),'Invalid common comparison scale')
    y1,y2=inputs['y1'],inputs['y2']
    require(y1.shape==(3,3) and y2.shape==(3,3),'Expected 3x3 Yukawa matrices')
    pf=np.zeros((3,3),complex)
    finite=np.zeros((3,3),complex)
    records=[]
    for r,m in enumerate(inputs['mf']):
        mass=float(m)
        require(mass>0 and math.isfinite(mass),'Bad MF')
        L=np.column_stack((y1[:,r].conjugate(),y2[:,r].conjugate()))
        hbar_over_m=1/(16*math.pi**2*mass)
        p=hbar_over_m*L@coefficient_matrix(vectors['pole'])@L.T
        logarithm=2*math.log(mu/mass)
        f=hbar_over_m*L@coefficient_matrix(
            vectors['constant']+logarithm*vectors['log'])@L.T
        pf+=p; finite+=f
        records.append({'fermion_eigenstate':r+1,'MF_gev':mass,
                        'pole_residue_c5_gev_inverse':encode_complex(p),
                        'finite_plus_MF_log_c5_gev_inverse':encode_complex(f),
                        'log_mu2_over_MF2':logarithm})
    require(relerr(pf,pf.T)<1e-12 and relerr(finite,finite.T)<1e-12,
            'F inherited Majorana C5 must be symmetric')
    return pf,finite,records


def assess_poles(scalar_pole: np.ndarray, inherited_pole: np.ndarray,
                 atol_relative: float=2e-8) -> dict:
    """Never set the inherited residue equal to minus the scalar residue."""
    relative=float(np.linalg.norm(scalar_pole+inherited_pole) /
                   max(np.linalg.norm(scalar_pole),np.linalg.norm(inherited_pole),1e-40))
    return {'common_pole_residual_relative':relative,
            'common_pole_cancellation_verified_in_raw_projected_sectors':relative<atol_relative,
            'absolute_frobenius_residual_gev_inverse':float(np.linalg.norm(scalar_pole+inherited_pole)),
            'hard_pole_norm_gev_inverse':float(np.linalg.norm(scalar_pole)),
            'inherited_F_pole_norm_gev_inverse':float(np.linalg.norm(inherited_pole))}


def verify_inputs(model: str, hard_run: str, fixed_run: str,
                  archive_root: Path) -> dict:
    """Read-only provenance: gate both native and numerical archive results."""
    source,fingerprints=batch6.verify_hard(model,hard_run)
    d=native.OUTPUT/hard_run/native.tag(model)
    hard_report=read_json(d/'hard_scalar_report.json')
    recon=read_json(d/'hard_scalar_reconstruction_audit.json')
    gate=read_json(native.OUTPUT/hard_run/'component_gate.json')
    require(hard_report.get('component_equivalent_native_reconstruction_used') is True
            and hard_report.get('full_EFT1_CheckLagrangian_passed') is True
            and hard_report.get('independent_three_mass_native_matching') is False,
            model+': not a valid gated Batch 8 native source')
    require(recon.get('component_equivalence_proven') is True and
            recon.get('independent_matchete_check_passed') is True,
            model+': native reconstructed EFT1 not independently checked')
    require(gate.get('component_equivalence_verified') is True and
            gate.get('full_pole_cancellation_certified') is False and
            gate.get('physical_C5_certified') is False,
            model+': component certificate unavailable/incompatible')
    fixed=read_json(batch6.OUTPUT/fixed_run/native.tag(model)/'result.json')
    require(fixed.get('model')==model and
            fixed.get('status')=='NativeHardThreeFlavorPoleSeparated_NotInheritedPoleCancelled'
            and fixed.get('quartic_map_is_conditional') is True and
            fixed.get('physical_C5_certified') is False and
            fixed.get('checks',{}).get('three_distinct_MF_lift_inferred_from_common_MF_symbol') is True,
            model+': not the requested 3MF conditional Batch 8 numerical result')
    require(fixed.get('hard_scalar_source_sha256')==digest(source) and
            fixed.get('source_fingerprints',{}).get('hard_candidate_sha256')==digest(source),
            model+': Batch 8 native scalar source changed')
    ps=native.paths(model)
    require(sha256(ps['running'].read_bytes()).hexdigest()==
            fixed.get('source_fingerprints',{}).get('source_fingerprints',{}).get('running')
            and digest(ps['parameters'])==
            fixed.get('source_fingerprints',{}).get('source_fingerprints',{}).get('parameters'),
            model+': current running parameters differ from the completed fixed-order run')
    archive_source=resolve_f_first_archive(archive_root,model)
    archive=archive_source['archive']
    return {'hard_source':source,'hard_stage_fingerprints':fingerprints,
            'fixed_result':fixed,'running':read_json(ps['running']),
            'parameters':read_json(ps['parameters']),
            'archive':archive,'archive_evidence':archive_source,
            'native_gate_sha256':digest(native.OUTPUT/hard_run/'component_gate.json'),
            'fixed_result_sha256':digest(batch6.OUTPUT/fixed_run/native.tag(model)/'result.json')}


def calculate_one(model: str, evidence: dict, root: Path,
                  *, vev: float=246.22) -> dict:
    src=evidence['hard_source']; prior=evidence['fixed_result']
    native_f=load_native_f_polynomials(root,evidence['archive'],model)
    source_convention=verify_scalar_extraction_convention(root)
    inp=numerical_inputs(model,evidence['running'],evidence['parameters'],
                         conditional_quartic_map=True)
    mu=float(prior['muS_gev'])
    require(mu>0 and math.isfinite(mu),'Invalid frozen scalar matching scale')
    # Re-evaluate the native Batch 8 hard source rather than trusting JSON
    # numbers and do not modify the saved result file.
    recomputed=evaluate_hard(model,src,evidence['running'],evidence['parameters'],
                             conditional_quartic_map=True,muS=mu,vev_gev=vev)
    hard_pole=decode_matrix(recomputed['hard_pole_residue_c5_gev_inverse'],'hard pole')
    hard_finite=decode_matrix(recomputed['hard_finite_c5_gev_inverse'],'hard finite')
    require(relerr(hard_pole,decode_matrix(prior['hard_pole_residue_c5_gev_inverse'],
            'saved pole'))<2e-11 and
            relerr(hard_finite,decode_matrix(prior['hard_finite_c5_gev_inverse'],
            'saved finite'))<2e-11,
            model+': current native hard does not reproduce the Batch 8 numerical output')
    projected=evaluate_polynomials(native_f['poly'],inp)
    f_pole,f_fin,per_f=f_inherited_matrices(inp,projected,mu)
    # Source convention: ExtractWeinbergCoefficient[hard]["Coefficient"] is
    # ORDERED A5. Batch 8 recorded Asym=(A5+A5.T)/2, but the F one-loop
    # projector already computes half-normalized C5=A5+A5.T. Therefore
    # C5_scalar=2*Asym_scalar, for BOTH pole and finite coefficients.
    scalar_c5_pole=2*hard_pole
    scalar_c5_finite=2*hard_finite
    pole_check=assess_poles(scalar_c5_pole,f_pole)
    pole_check['ordered_batch8_unconverted_relative_residual']=assess_poles(
        hard_pole,f_pole)['common_pole_residual_relative']
    # An independent scale shift checks the logarithmic slope of BOTH
    # archived F-first and current scalar hard expressions.
    mu_shift=mu*math.exp(0.051)
    require(mu_shift<inp['muF'],'Scale-shift check exceeds frozen F threshold')
    h2=evaluate_hard(model,src,evidence['running'],evidence['parameters'],
                     conditional_quartic_map=True,muS=mu_shift,vev_gev=vev)
    f_pole2,f_fin2,_=f_inherited_matrices(inp,projected,mu_shift)
    require(relerr(f_pole,f_pole2)<1e-12,'F inherited pole should not depend on mu')
    finite_sum=scalar_c5_finite+f_fin
    finite_sum2=2*decode_matrix(h2['hard_finite_c5_gev_inverse'],
                               'shifted ordered scalar A5')+f_fin2
    scale_error=relerr(finite_sum,finite_sum2)
    passed=(pole_check['common_pole_cancellation_verified_in_raw_projected_sectors']
            and scale_error<2e-8)
    result={'format':FORMAT,'model':model,
            'status':STATUS if passed else 'ProjectedRawPoleOrScaleMismatch_NotCombined',
            'source_archived_F_first_one_loop':{
                'path':native_f['archive_path'],'sha256':native_f['archive_sha256'],
                'sector_terms':native_f['terms'],
                'SU2_neutral_PL_PR_conjugacy_checked':True,
                'pole_log_coefficients_equal':True,
                'finite_constant_equals_pole':native_f['F_finite_constant_equals_pole']},
            'source_batch8':{'hard_scalar_sha256':digest(src),
                'native_gate_sha256':evidence['native_gate_sha256'],
                'fixed_result_sha256':evidence['fixed_result_sha256']},
            'scalar_ordered_to_symmetric_C5_convention':source_convention,
            'comparison_scale_gev':mu,'source_couplings_frozen_at':'muF',
            'three_independent_MF_inferred_not_native_one_loop_matched':True,
            'native_quartic_to_RGBeta_map_certified':False,
            'hard_scalar_ordered_A5_sym_pole_residue_gev_inverse':encode_complex(hard_pole),
            'hard_scalar_pole_residue_c5_gev_inverse':encode_complex(scalar_c5_pole),
            'archived_F_inherited_pole_residue_c5_gev_inverse':encode_complex(f_pole),
            'sum_pole_residue_c5_gev_inverse':encode_complex(scalar_c5_pole+f_pole),
            'source_pole_check':pole_check,
            'common_scale_shift_gev':mu_shift,
            'common_scale_finite_matrix_relative_change':scale_error,
            'per_fermion':per_f,
            'UV_EFT_matching_counterterms_independently_verified':False,
            'EFT1_RGE_evolution_performed':False,
            'full_UV_EFT1_scheme_matching_certified':False,
            'complete_physical_C5':None,
            'physical_C5_certified':False}
    if passed:
        result.update({
            'conditional_raw_hard_finite_ordered_A5_sym_gev_inverse':encode_complex(hard_finite),
            'conditional_raw_hard_finite_c5_gev_inverse':encode_complex(scalar_c5_finite),
            'conditional_raw_F_inherited_finite_plus_log_c5_gev_inverse':encode_complex(f_fin),
            'conditional_common_scale_sum_c5_gev_inverse':encode_complex(finite_sum),
            'conditional_common_scale_sum_singular_values_ev':
                metrics(finite_sum,vev)['conditional_takagi_masses_ev_ascending'],
            'conditional_common_scale_mu_independence_test_passed':True})
    return result


def write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main(argv: list[str]|None=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hard-run',default='mixed_hard_08')
    p.add_argument('--fixed-run',default='mixed_fixed_08')
    p.add_argument('--run-id',default='mixed_inherited_09')
    p.add_argument('--models',default='B_-1,C_-1')
    p.add_argument('--archive-root',type=Path,default=ARCHIVE)
    p.add_argument('--check-only',action='store_true')
    a=p.parse_args(argv)
    for name in (a.hard_run,a.fixed_run,a.run_id):
        if not re.fullmatch(r'[A-Za-z0-9_-]+',name):
            p.error('Run IDs must be alphanumeric/_/-')
    models=tuple(s.strip() for s in a.models.split(','))
    if not models or len(set(models))!=len(models) or any(s not in MODELS for s in models):
        p.error('Only B_-1,C_-1 once each')
    ready={};errors={}
    for model in models:
        try:
            e=verify_inputs(model,a.hard_run,a.fixed_run,a.archive_root)
            ready[model]=e
            print(model+': READY, Batch 8 and '+e['archive_evidence']['provenance']+
                  '; source='+str(e['archive']),flush=True)
        except Exception as exc:
            errors[model]=str(exc)
            print(model+': BLOCKED: '+str(exc),file=sys.stderr,flush=True)
    if a.check_only:
        return 0 if not errors else 2
    out=OUTPUT/a.run_id
    require(not out.exists(),'Output already exists; select a new run ID: '+str(out))
    out.mkdir(parents=True)
    statuses={}
    for model in models:
        d=out/native.tag(model)
        if model not in ready:
            statuses[model]='BlockedMissingInput'
            write_json(d/'error.json',{'model':model,'status':'Blocked',
                'reason':errors[model],'physical_C5_certified':False})
            continue
        try:
            result=calculate_one(model,ready[model],ROOT)
            write_json(d/'result.json',result)
            statuses[model]=result['status']
            chk=result['source_pole_check']
            print(f'{model}: {result["status"]}; F + scalar raw pole residual='
                  f'{chk["common_pole_residual_relative"]:.3e}; '
                  f'common-scale change={result["common_scale_finite_matrix_relative_change"]:.3e}',
                  flush=True)
            if result['status']==STATUS:
                print('  conditional common-scale singular values [eV]: '+
                      str(result['conditional_common_scale_sum_singular_values_ev']),flush=True)
        except Exception as exc:
            statuses[model]='BlockedAnalysisError'
            write_json(d/'error.json',{'model':model,'status':'Blocked',
                'reason':str(exc),'physical_C5_certified':False})
            print(model+': BLOCKED: '+str(exc),file=sys.stderr,flush=True)
    success=all(statuses.get(model)==STATUS for model in models)
    write_json(out/'summary.json',{'format':FORMAT,'status':
        'ConditionalRawPoleAndCommonScaleChecksPassed' if success else 'IncompleteOrUncancelled',
        'models':statuses,'hard_run':a.hard_run,'fixed_run':a.fixed_run,
        'archived_F_one_loop_source_reused':all(m in ready for m in models),
        'archived_native_F_source_fallback_supported':True,
        'scalar_A5_to_C5_conversion_required':True,
        'scalar_A5_to_C5_factor_checked_for_all_model_results':all(
            (out/native.tag(model)/'result.json').is_file() for model in models),
        'fermion_Match_recomputed':False,
        'scalar_Match_recomputed':False,'ordinary_models_rerun':False,
        'complete_physical_C5':None,'physical_C5_certified':False,
        'renormalized_UV_EFT_counterterms_certified':False,
        'scope':'F-first raw LLHH C5 and Batch8 scalar ordered A5 converted to half-normalized C5; frozen common couplings'})
    print('Evidence report: '+str(out/'summary.json'),flush=True)
    return 0 if success else 2


if __name__=='__main__':
    try:
        raise SystemExit(main())
    except Exception as exc:
        print('BATCH9 ERROR: '+str(exc),file=sys.stderr)
        raise SystemExit(2)
