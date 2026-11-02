"""Independent full SU(2) LLSS-response verification and conditional scalar handoff.

All real T3 sources are read-only.  This is a formal Matchete PL-polynomial
cross-check, never a proof of absolute GammaCC/Weyl normalization.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import Any
import json
import math

import numpy as np

from .core import project_weinberg
from .fermion_boundary_bridge import SourceBoundary, read_json, require, digest
from .large_batch import ScalarBackground, _rotate, parse_matrix, _beta_project

FORMAT = 't3-eft1-full-su2-source-response-v1'
STATUS = 'SourceFactorisedFullSU2_NoQuartic'
SUMMARY_FORMAT = 't3-eft1-full-su2-source-response-batch-v1'
SUMMARY_STATUS = 'FullSU2SourceResponsesOnly_NotPhysicalC5'
EXPECTED_CONTRACTION = ('beta[p,q,i,j,H_A,H_B] = sum_{a<=b}(2-delta_ab)*'
                        'lambda[H_A,H_B,a,b]*C[p,q,i,j,a,b]')
EXPECTED_WEINBERG = '(1/2) C5[p,q] (nu_p H0 - e_p H+) (nu_q H0 - e_q H+)'


def parse_complex(raw: Any, shape: tuple[int, ...], name: str) -> np.ndarray:
    require(isinstance(raw, dict) and set(raw) == {'real', 'imag'},
            f'{name}: expected explicit real/imag parts')
    real = np.asarray(raw['real'], dtype=float)
    imag = np.asarray(raw['imag'], dtype=float)
    require(real.shape == shape and imag.shape == shape,
            f'{name}: invalid array shape (expected {shape})')
    value = real + 1j * imag
    require(np.isfinite(value).all(), f'{name}: nonfinite source coefficient')
    return value


def _close(a: np.ndarray, b: np.ndarray, label: str) -> float:
    error = float(np.max(np.abs(a - b)))
    scale = max(float(np.max(np.abs(a))), float(np.max(np.abs(b))), 1e-30)
    require(error <= 2e-9*scale + 1e-28,
            f'{label}: full SU2 source differs from fermion boundary (max {error:.4g})')
    return error


def validate_su2_response(path: Path, source: SourceBoundary) -> dict:
    """Verify every independently exported PL source, SU2 component and scalar pair.

    Only unordered scalar pairs with nonzero formal coefficient are exported;
    a pair can cancel in the sum even when individual PL sources are nonzero.
    """
    o = read_json(path)
    require(o.get('format') == FORMAT and o.get('status') == STATUS,
            'unexpected full-SU2 source format/status')
    require(o.get('model') == source.model and o.get('source_boundary_sha256') == source.hash,
            'full-SU2 response is not tied to the current fermion release')
    require(float(o.get('mu_F_gev', -1)) == source.mu_f_gev and
            o.get('mu_S_gev') is None and
            o.get('N_GammaCC_to_Weyl') is None,
            'full-SU2 response is muF-only with undefined GammaCC factor')
    require(o.get('source_count') == len(source.source_indices),
            'source count differs from upstream boundary')
    for key in ('absolute_spinor_normalization_certified', 'physical_C5_computed',
                'complete_one_loop_operator_inventory_certified',
                'one_loop_fermion_hard_matching_included',
                'finite_joint_scalar_matching_included'):
        require(o.get(key) is False, f'{key}: false certification required')
    require(o.get('full_su2_contraction_convention') == EXPECTED_CONTRACTION,
            'SU2 response tensor convention changed')
    require(o.get('weinberg_tensor_convention') == EXPECTED_WEINBERG,
            'Weinberg convention changed; do not silently reinterpret')

    rows = o.get('source_resolved_blocks')
    require(isinstance(rows, list) and len(rows) == len(source.source_indices),
            'missing or repeated source SU2 blocks')
    n = len(source.labels)
    maxdiff = 0.0
    pairs_verified = 0
    source_nonzero = 0
    source_cancelled_pairs: dict[tuple[int, int], np.ndarray] = {}
    seen_source_indices: set[int] = set()
    for s, row in enumerate(rows):
        idx = row.get('matchete_source_index')
        require(idx == source.source_indices[s] and idx not in seen_source_indices,
                'Matchete PL source list order/identity changed')
        seen_source_indices.add(idx)
        entries = row.get('scalar_pairs')
        require(isinstance(entries, list), f'PL {idx}: scalar pair list missing')
        expected_pairs = {
            (a,b) for a in range(n) for b in range(a,n)
            if bool(np.any(source.by_source[s,...,a,b] != 0))
        }
        seen: set[tuple[int, int]] = set()
        for entry in entries:
            pair = entry.get('heavy_pair_one_based')
            require(isinstance(pair,list) and len(pair)==2 and
                    all(type(x) is int for x in pair), f'PL {idx}: malformed scalar pair')
            a,b = pair[0]-5,pair[1]-5
            require(0 <= a <= b < n and (a,b) not in seen,
                    f'PL {idx}: repeated/out-of-range/noncanonical scalar pair')
            seen.add((a,b))
            multiplicity = 1 if a==b else 2
            require(entry.get('ordered_pair_multiplicity') == multiplicity,
                    f'PL {idx} pair {(a,b)}: wrong ordered multiplicity')
            actual = parse_complex(entry.get('LLSS_C_pqij_per_NGamma_GeV_inv'),
                                   (3,3,2,2), f'PL {idx} pair {(a,b)}')
            expected = source.by_source[s,...,a,b]
            maxdiff = max(maxdiff, _close(actual, expected,
                         f'PL {idx} pair {(a,b)}'))
            source_cancelled_pairs[(a,b)] = (source_cancelled_pairs.get(
                (a,b), np.zeros((3,3,2,2), complex)) + actual)
            pairs_verified += 1
            if bool(np.any(actual!=0)):
                source_nonzero += 1
        require(seen == expected_pairs,
                f'PL {idx}: full-SU2 response omitted or added scalar pairs; '
                f'missing={sorted(expected_pairs-seen)}, extra={sorted(seen-expected_pairs)}')
    for (a,b), tensor in source_cancelled_pairs.items():
        maxdiff = max(maxdiff, _close(tensor, source.raw[...,a,b],
                         f'combined scalar pair {(a,b)}'))
    return {
        'full_su2_source_file_sha256': digest(path),
        'source_boundary_sha256': source.hash,
        'source_count_verified': len(rows),
        'unordered_source_scalar_pairs_verified': pairs_verified,
        'nonzero_source_scalar_pairs': source_nonzero,
        'distinct_pair_count_before_cancellation': len(source_cancelled_pairs),
        'distinct_pair_count_after_cancellation': sum(
            bool(np.any(source.raw[...,a,b] != 0)) for a,b in source_cancelled_pairs),
        'max_component_abs_error_gev_inverse': maxdiff,
        'full_two_lepton_SU2_component_parity': True,
        'physical_GammaCC_normalization_certified': False,
        'full_one_loop_certified': False,
    }


def validate_su2_manifest(dir_path: Path, release_dir: Path, *,models: tuple[str,...]) -> dict:
    report = read_json(dir_path/'summary.json')
    require(report.get('format') == SUMMARY_FORMAT and
            report.get('status') == SUMMARY_STATUS,
            'full-SU2 batch summary version/status changed')
    require(report.get('model_count') == 16 and
            report.get('fermion_release_summary_sha256') == digest(release_dir/'summary.json'),
            'full-SU2 source manifest disagrees with 16-model fermion source release')
    require(report.get('quartic_input_kind') == 'none' and
            report.get('physical_C5_computed') is False and
            report.get('source_spinor_normalization_certified') is False,
            'source SU2 response has been mislabeled as physical')
    rows = report.get('models')
    require(isinstance(rows,list) and len(rows)==16, 'expected all 16 source SU2 records')
    indexed={}
    for row in rows:
        model=row.get('model')
        require(isinstance(model,str) and model not in indexed,
                'invalid or duplicate model in SU2 manifest')
        require(row.get('status') == STATUS, 'SU2 model is not in expected source-only state')
        indexed[model]=row
    require(set(models).issubset(indexed), 'missing requested model in full SU2 manifest')
    return indexed


def independent_source_beta_check(source: SourceBoundary,
                                  background: ScalarBackground,
                                  *, relative_tol: float = 2e-9) -> dict:
    """Contract all source scalar pairs directly, then compare with rotated beta.

    This check is independent of the per-source quartic-pivot (H0R-only) worker.
    It tests all four Higgs coordinates and both lepton SU2 indices at once.
    The source carries no certified spinor factor; the beta remains per unit N.
    """
    U=background.U
    q=background.quartic
    q_src=np.einsum('ABkl,ak,bl->ABab',q,U,U,optimize=True)
    C=source.raw
    # Distinct unordered a,b: C_ab = C_ba, q_ab = q_ba.
    full=np.zeros((3,3,2,2,4,4),complex)
    for a in range(len(source.labels)):
        for b in range(a,len(source.labels)):
            full += (1 if a==b else 2)*np.einsum(
                'fgij,AB->fgijAB',C[...,a,b],q_src[...,a,b],optimize=True)
    direct, direct_res=project_weinberg(full,strict=False)
    rotated,_=_beta_project(_rotate(C,U),q,strict=False)
    diff=float(np.max(abs(direct-rotated)))
    scale=max(float(np.max(abs(direct))),float(np.max(abs(rotated))),1e-30)
    require(diff <= relative_tol*scale+1e-28,
            f'full-SU2 unordered scalar-pair beta differs from rotated tensor: {diff}')
    return {'full_SU2_unordered_pair_beta_parity':True,
            'beta_C5_16pi2_per_N_max_abs_error_gev_inverse':diff,
            'weinberg_projection_residual':direct_res,
            'quartic_mass_at_single_muS':background.mu_s,
            'conditional_on_unknown_GammaCC_normalization':True}


def audit_signed_evidence(native: Path, signed: Path,
                          source_release_dir: Path) -> dict:
    """Read existing evidence without inferring GammaCC phase from contact algebra.

    Missing reports remain explicit; a blocked native Feynman rule report does
    not invalidate the source-polynomial audit, but does block physical C5.
    """
    result={'absolute_matchete_to_weyl_factor':None,
            'physical_spinor_normalization_certified':False,
            'source_boundary_release_sha256':digest(source_release_dir/'summary.json'),
            'blocking_items':[]}
    for label,path in (('native_feynman',native),('grassmann_contact',signed)):
        if not path.is_file():
            result[label]={'status':'evidence_file_missing','path':str(path)}
            result['blocking_items'].append(f'{label} evidence missing')
            continue
        doc=read_json(path)
        result[label]={'status':doc.get('status'), 'sha256':digest(path),
                       'path':str(path),'format':doc.get('format')}
        if label=='native_feynman':
            require(doc.get('format')=='t3-native-matchete-fr-evidence-v1',
                    'unsupported native Feynman-rule evidence format')
            result[label]['blocked_phase']=doc.get('blocked_phase')
            result[label]['reason']=doc.get('reason')
            if doc.get('physical_master_wilson_normalization_certified') is not True:
                result['blocking_items'].append('native signed UV-to-EFT1 four-point amplitude not certified')
        else:
            require(doc.get('format')=='t3-signed-zero-momentum-grassmann-diagnostic-v1',
                    'unsupported signed Grassmann check format')
            require(doc.get('input_release_summary_sha256')==result['source_boundary_release_sha256'],
                    'signed Grassmann evidence is not for current F-boundary release')
            if doc.get('assumed_GammaCC_phase_certified') is not True:
                result['blocking_items'].append('declared Weyl contact does not fix native GammaCC phase')
    result['blocking_items'].append('independent numerical Matchete-to-master signed normalization certificate absent')
    return result
