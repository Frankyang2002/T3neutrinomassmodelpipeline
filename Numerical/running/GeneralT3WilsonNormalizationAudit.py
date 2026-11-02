"""All-16 Matchete LLSS -> master polynomial *normalization audit*.

This audits source projector/chirality, correct lepton/flavour dummy pairing,
permutation combinatorics, and one-generation master polynomial convention.
It explicitly does NOT identify Matchete's four-component GammaCC PL bilinear
with the master two-component Weyl bilinear or certify an absolute sign/factor.

In particular, the old B/C mixed converter mislabeled the conjugate Proj[1]
sector as PL. Historical mixed outputs must fail this audit until rebuilt.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import math
from pathlib import Path
import re
from typing import Any, Mapping

import numpy as np

from Numerical.running.GeneralT3WilsonPolynomialBatch import (
    FORMAT as BRIDGE_FORMAT, STATUS as BRIDGE_STATUS,
    canonical_key, distinct_orbit, build_polynomial_bridge,
    master_weinberg_polynomial_check,
)

FORMAT = 't3-16-matchete-to-master-normalization-audit-v1'
STATUS = 'PolynomialAndPLStructureVerified_SpinorNormalizationOpen'
PL = re.compile(r'Proj\s*\[\s*-1\s*\]')
PR = re.compile(r'Proj\s*\[\s*1\s*\]')
TRANSPOSE_PLAIN = re.compile(r'Transp\s*\[\s*Field\s*\[\s*l\s*,\s*Fermion')
TRANSPOSE_BAR = re.compile(r'Transp\s*\[\s*Bar\s*\[\s*Field\s*\[\s*l\s*,\s*Fermion')


def _require(ok: bool, why: str) -> None:
    if not ok:
        raise ValueError(why)


def digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def identify_chirality(text: str) -> str:
    """Validate physical PL/PR bilinear against the native RunMatching.wl rule."""
    _require(isinstance(text, str) and 'GammaCC' in text,
             'Term lacks Matchete charge-conjugation bilinear')
    pl = bool(PL.search(text))
    pr = bool(PR.search(text))
    _require(pl != pr, 'Term has missing/ambiguous native chirality projector')
    if pl:
        _require(bool(TRANSPOSE_PLAIN.search(text)) and
                 not bool(TRANSPOSE_BAR.search(text)),
                 'Proj[-1] term is not the plain L-transpose Majorana operator')
        return 'PL'
    _require(bool(TRANSPOSE_BAR.search(text)) and
             not bool(TRANSPOSE_PLAIN.search(text)),
             'Proj[1] term is not the conjugate barred-lepton operator')
    return 'PR'


def verify_native_convention(wolfram_source: str) -> dict:
    """Check exact declared Matchete chirality/C5 rules; NOT GammaCC matrix."""
    text = str(wolfram_source)
    _require(re.search(r'"Proj\[-1\]"\]\s*,\s*"PL"', text) is not None,
             'RunMatching.wl no longer maps Proj[-1] to PL')
    _require(re.search(r'"Proj\[1\]"\]\s*,\s*"PR"', text) is not None,
             'RunMatching.wl no longer maps Proj[1] to PR')
    _require('"L_EFT superset (1/2) C5^(ij) O_ij' in text and
             'C5^(ij) = A5^(ij) + A5^(ji)' in text,
             'Expected reported Weinberg C5 convention absent or changed')
    return {
        'matchete_pl_projector':'Proj[-1]',
        'matchete_pr_conjugate_projector':'Proj[1]',
        'reported_weinberg_coefficient':'L_EFT contains (1/2) C5^(ij) O_ij + h.c.',
        'reported_ordered_A5_to_C5':'C5^(ij)=A5^(ij)+A5^(ji)',
        'the_reported_C5_rule_is_not_assumed_for_LLSS':True,
        'GammaCC_to_two_component_spinor_phase_certified':False,
    }


def verify_source_chiralities(model: str, ordered: Mapping[str,Any],
                              source: Mapping[str,Any]) -> dict:
    terms = source.get('operator_terms')
    _require(isinstance(terms, list) and bool(terms),
             f'{model}: missing matched Matchete source term inventory')
    counts = Counter()
    problems = []
    for t in terms:
        raw = t.get('source_term_inputform')
        if not isinstance(raw,str):
            problems.append('source term '+str(t.get('source_index'))+' lacks exact InputForm')
            continue
        actual = identify_chirality(raw)
        declared = t.get('chirality', 'PL' if model not in ('B_-1','C_-1') else None)
        if declared != actual:
            problems.append(f"source term {t.get('source_index')} tagged {declared}, native Matchete is {actual}")
        counts[actual] += 1
        _require('Coupling[MF,' in raw,
                 f'{model}: matched tree term lacks its heavy mass source')
    _require(not problems, f'{model}: native chirality selection mismatch: '+ '; '.join(problems[:6]))
    expect = int(ordered['tree_operator_count_PL'])
    _require(counts['PL'] == expect, f'{model}: PL count does not match tensor')
    if model in ('B_-1','C_-1'):
        _require(counts['PR'] == expect,
                 f'{model}: incomplete native PL/PR Hermitian-conjugate inventory')
    else:
        _require(counts['PR'] == 0, f'{model}: unexpected PR source term in PL boundary')
    return {
        'native_PL_matched_term_count':counts['PL'],
        'native_PR_conjugate_term_count':counts['PR'],
        'selected_tensor_is_holomorphic_PL':True,
        'Hermitian_conjugate_not_summed_into_PL_twice':True,
        'native_chirality_structure_verified':True,
    }


def verify_saved_polynomial(model: str, ordered: Mapping[str,Any],
                            source: Mapping[str,Any], bridge: Mapping[str,Any]) -> dict:
    _require(bridge.get('format')==BRIDGE_FORMAT and bridge.get('status')==BRIDGE_STATUS
             and bridge.get('model')==model, f'{model}: invalid polynomial bridge')
    _require(bridge.get('overall_matchete_spinor_normalization_certified') is False,
             f'{model}: unexpected source spinor normalization claim')
    regenerated = build_polynomial_bridge(model,ordered,source)
    actual = bridge.get('independent_symmetric_coordinates')
    expected = regenerated['independent_symmetric_coordinates']
    _require(isinstance(actual,list) and len(actual)==len(expected),
             f'{model}: symmetric coordinate inventory changed')
    orbit_hist = Counter()
    for row, calculated in zip(actual,expected):
        _require(row['coordinate']==calculated['coordinate'] and
                 row['orbit_size']==calculated['orbit_size'],
                 f'{model}: coordinate/orbit mismatch')
        _require(row['orbit_size']==len(distinct_orbit(tuple(row['coordinate']))),
                 f'{model}: incorrect scalar/fermion identical-field multiplicity')
        x=complex(row['coefficient_real'],row['coefficient_imag'])
        y=complex(calculated['coefficient_real'],calculated['coefficient_imag'])
        _require(all(map(math.isfinite,(x.real,x.imag))) and
                 abs(x-y)<=1e-12*max(abs(y),1e-30),
                 f'{model}: normalized polynomial coefficient differs from regenerated source')
        orbit_hist[str(row['orbit_size'])]+=1
    _require(regenerated['polynomial_identity']['preserved_to_numerical_precision'],
             f'{model}: source/symmetric polynomial reconstruction failed')
    return {
        'independent_symmetric_coordinate_count':len(actual),
        'distinct_permutation_orbit_histogram':dict(sorted(orbit_hist.items())),
        'full_symmetric_tensor_requires_no_extra_combinatorial_factor':True,
        'formal_polynomial_reconstruction_residual':
            regenerated['polynomial_identity']['relative_monomial_difference'],
        'master_two_component_spinor_identification_certified':False,
    }


def audit_model(model: str, ordered: Mapping[str,Any],
                source: Mapping[str,Any], bridge: Mapping[str,Any]) -> dict:
    chirality = verify_source_chiralities(model, ordered,source)
    polynomial=verify_saved_polynomial(model,ordered,source,bridge)
    return {
        'format':FORMAT,'model':model,'status':STATUS,
        'mu_F_gev':bridge['mu_F_gev'],
        'native_matchete_chirality':chirality,
        'master_polynomial_coordinate_check':polynomial,
        'normalization_findings':{
            'symmetric_orbit_division_already_applied':True,
            'do_not_multiply_by_scalar_identical_field_factor_again':True,
            'do_not_add_the_PR_Hermitian_conjugate_to_the_PL_boundary':True,
            'absolute_matchete_GammaCC_to_master_Weyl_bilinear_factor':'UNRESOLVED',
            'absolute_matchete_GammaCC_to_master_Weyl_bilinear_sign':'UNRESOLVED',
            'matchete_to_master_spinor_normalization_certified':False,
            'physical_master_rge_dictionary_certified':False,
            'three_family_gamma_A_B_constructed':False,
            'physical_eft1_wilson_running_performed':False,
        },
    }


def master_evidence(root: Path) -> dict:
    repo = Path(root)
    files={
        'matchete_matching':'Lagrangian/RunMatching.wl',
        'master_weinberg_adapter':'RGE/matching/WeinbergTensorAdapter.py',
        'master_wilson_rge':'RGE/general/WilsonTensorRGE.py',
    }
    results={}
    for label,rel in files.items():
        path=repo/rel
        _require(path.is_file(),f'Missing repository convention file {rel}')
        results[label]={'project_path':rel,'sha256':digest(path)}
    conventions=verify_native_convention((repo/files['matchete_matching']).read_text(encoding='utf-8-sig'))
    master=master_weinberg_polynomial_check()
    return {'files':results,'source_declared_conventions':conventions,
            'internal_master_polynomial_example':master,
            'physical_spinor_normalization_certified':False}
