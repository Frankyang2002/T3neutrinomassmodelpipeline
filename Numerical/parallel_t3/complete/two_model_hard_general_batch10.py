"""Batch 10: source-backed A(-4)/D(-2) hard + general-HHSS numerical handoff.

The F-hard entry is the *raw* ordered Matchete LLHH regular part, NOT a
UV-minus-EFT1-soft renormalized Wilson boundary.  The D(-2) extra HHSS
channels use the existing exact SU(2)-invariant reference basis; the
RGBeta-to-Matchete coupling coordinates and GammaCC/Weyl factor remain open.

No result of this module may be called a physical full one-loop C5.
The O(hbar) LLSS scalar-pair logarithm and simultaneous finite two-scalar
bubble are included exactly ONCE; the LLSS self-running contributes to LLHH
at O(hbar**2), not at the one-loop order claimed here.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any, Mapping
import json
import math

import numpy as np

from Numerical.parallel_t3.complete.two_model_joint_scalar_batch9 import (
    B1_SUBDIR, B2_SUBDIR, MODELS, load_snapshot, insist,
)

B9 = Path('output/parallel_t3/thesis_full_general_joint_scalar_batch9_v1')
LAURENT = Path('output/parallel_t3/full_general_batch2/thesis_full_general_batch2_laurent_v2')
EXTRA_D = Path('output/general_observable/lo_certification_couplings/'
               'ws2_quartic_completion_v1/D_m2/exact_extra_HHSS_reference_tensors.json')
REGISTRY = Path('RGE/general/generated/general_rgbeta_quartic_registry.json')
FORMAT = 't3-thesis-hard-general-hhss-batch10-v1'
REF_FORMAT = 't3-ws2-extra-hhss-exact-su2-reference-space-v1'
EXPECTED_D = {
    (0, 2, 0, 1, 0, 1): 'lambdaT3GeneralQ020101I1',
    (1, 1, 0, 0, 0, 2): 'lambdaT3GeneralQ110002I1',
}
LOOP = 16.0 * math.pi**2


def digest(path: Path) -> str:
    insist(Path(path).is_file(), 'Missing source file: ' + str(path))
    return sha256(Path(path).read_bytes()).hexdigest()


def matrix_json(a: np.ndarray) -> dict:
    x = np.asarray(a, dtype=np.complex128)
    insist(x.shape == (3, 3) and np.isfinite(x).all(), 'Expected finite 3x3 matrix')
    return {'real': x.real.tolist(), 'imag': x.imag.tolist()}


def matrix(value: Any, name: str = 'coefficient', *, symmetric: bool = True) -> np.ndarray:
    """Decode both actual source-specific complex 3x3 JSON representations.

    Batch 9 serialises one ``{'real': 3x3, 'imag': 3x3}`` matrix object.
    ``numerical_unsubtracted_source`` instead calls ``encode_native``, which
    serialises a 3x3 list of ``{'real': x, 'imag': y}`` entries.  Preserve
    all complex phases and the ordered flavour indices before any optional
    (complex *transpose*, never Hermitian) symmetrisation by the caller.
    """
    if isinstance(value, Mapping):
        insist(set(value) == {'real', 'imag'},
               name + ': expected a real/imag matrix object')
        try:
            real = np.asarray(value['real'], dtype=float)
            imag = np.asarray(value['imag'], dtype=float)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(name + ': invalid scalar values in matrix object') from exc
        insist(real.shape == imag.shape == (3, 3),
               name + ': real/imag matrix components must both be 3x3')
        a = real + 1j * imag
    elif isinstance(value, list):
        # The existing native-JSON encoder operates recursively per element.
        insist(len(value) == 3 and all(
            isinstance(row, list) and len(row) == 3 and all(
                isinstance(z, Mapping) and set(z) == {'real', 'imag'}
                for z in row) for row in value),
            name + ': expected a 3x3 list of real/imag complex entries')
        try:
            real = np.asarray([[z['real'] for z in row] for row in value], dtype=float)
            imag = np.asarray([[z['imag'] for z in row] for row in value], dtype=float)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(name + ': invalid native complex-entry values') from exc
        a = real + 1j * imag
    else:
        raise ValueError(name + ': unknown source complex-matrix JSON encoding')
    insist(a.shape == (3, 3) and np.isfinite(a).all(),
           name + ': coefficient must be a finite complex 3x3 matrix')
    if symmetric:
        insist(np.allclose(a, a.T, rtol=1e-9, atol=1e-23),
               name + ' not complex-symmetric')
    return a


def symmetric(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=complex)
    insist(a.shape == (3, 3) and np.isfinite(a).all(), 'Bad ordered 3x3 tensor')
    return a + a.T


def _coeff(value: Any, label: str) -> complex:
    """Decode explicitly supplied RGBeta scalar without ignoring its phase."""
    if isinstance(value, dict) and set(value) == {'real', 'imag'}:
        z = complex(float(value['real']), float(value['imag']))
    elif isinstance(value, (int, float, complex, np.number)) and not isinstance(value, bool):
        z = complex(value)
    else:
        raise ValueError(label + ': invalid numerical scalar encoding')
    insist(np.isfinite([z.real, z.imag]).all(), label + ' must be finite scalar')
    return z


def _reference_number(s: str) -> float:
    """Parse source-authored exact radical; no arbitrary executable code."""
    import sympy as sp
    from sympy.parsing.sympy_parser import parse_expr
    insist(isinstance(s, str) and s and len(s) <= 160 and
           all(c in '0123456789+-*/() .sqrt' for c in s),
           'Non-numeric exact SU2 tensor coefficient')
    expr = parse_expr(s, local_dict={'sqrt': sp.sqrt},
                      global_dict={'Integer': sp.Integer, 'Rational': sp.Rational,
                                   'Float': sp.Float, 'Pow': sp.Pow, 'Mul': sp.Mul,
                                   'Add': sp.Add}, evaluate=True)
    insist(bool(expr.is_number and expr.is_real and expr.is_finite),
           'Non-real exact SU2 tensor coefficient')
    return float(expr.evalf(17))


def reference_hhss_tensor(rows: list[dict], z: complex, *, kind: str) -> np.ndarray:
    """Explicit fourth derivative q[H_A,H_B,S_a,S_b], original real basis.

    Each source row supplies an ordered (A,B,a,b) component. Insert the
    bosonic A/B and a/b exchanges once, preserving all relative signs.
    The source convention for its complex invariant is V4=-cP-c*P*.
    """
    insist(kind in ('Re', 'Im') and isinstance(rows, list), 'Invalid reference channel')
    scale = z.real if kind == 'Re' else z.imag
    q = np.zeros((4, 4, 8, 8), float)
    seen: dict[tuple[int, int, int, int], float] = {}
    for row in rows:
        key = row.get('H1_H2_S1_S2_0based')
        insist(isinstance(key, list) and len(key) == 4 and
               all(type(k) is int for k in key) and
               all(0 <= key[i] < (4 if i < 2 else 8) for i in range(4)),
               'Invalid reference four-point tensor coordinate')
        a, b, c, d = key
        val = scale * _reference_number(row['value_exact'])
        for orb in {(a,b,c,d),(b,a,c,d),(a,b,d,c),(b,a,d,c)}:
            if orb in seen:
                insist(math.isclose(seen[orb], val, rel_tol=1e-12, abs_tol=1e-12),
                       'Conflicting exact reference components: ' + str(orb))
            seen[orb] = val
            q[orb] = val
    return q


def extra_hhss_from_source(doc: dict, registry: dict,
                            scalar_couplings: Mapping[str, Any]
                            ) -> tuple[dict[str, np.ndarray], dict]:
    """D2 extra channels; never silently drop an allowed HHSS coupling.

    Four other steps of this model's full generality remain separate:
    source CG/RGBeta normalization, finite F matching, complete scalar
    one-loop threshold inventory, and signed GammaCC spinor translation.
    """
    insist(doc.get('format') == REF_FORMAT and doc.get('model') == 'D_-2',
           'Wrong D2 exact reference geometry')
    row = registry.get('models', {}).get('D_-2', {})
    new = row.get('new_rgbeta_symbols', [])
    source_map = {
        tuple(entry['counts']): entry['name'] for entry in new
        if sum(entry['counts'][:2]) == 2 and sum(entry['counts'][2:]) == 2
    }
    insist(source_map == EXPECTED_D, 'D2 full-general HHSS registry changed')
    sectors = doc.get('monomial_sectors', [])
    insist(len(sectors) == len(EXPECTED_D), 'D2 extra-HHSS reference channels incomplete')
    channels: dict[str, np.ndarray] = {}
    details = {}
    for sector in sectors:
        counts = tuple(sector.get('counts_H_BarH_S1_BarS1_S2_BarS2', []))
        insist(counts in EXPECTED_D, 'Unexpected D2 HHSS reference sector')
        name = EXPECTED_D[counts]
        symbols = sector.get('source_registry_symbols', [])
        insist(symbols == [{'rgbeta_symbol': name, 'source': 'new'}],
               'Extra HHSS source-to-RGBeta registry changed')
        references = sector.get('exact_reference_real_component_tensor_basis', [])
        insist(len(references) == 1 and references[0].get('id') == 'reference_invariant_1',
               'D2 extra HHSS requires a unique exact reference invariant')
        insist(name in scalar_couplings and name not in channels,
               'Missing nonzero-capable D2 HHSS coupling ' + name)
        z = _coeff(scalar_couplings[name], name)
        ref = references[0]
        q = reference_hhss_tensor(ref['Re_c_unit_tensor_components'], z, kind='Re')
        q += reference_hhss_tensor(ref['Im_c_unit_tensor_components'], z, kind='Im')
        insist(np.isfinite(q).all() and np.max(abs(q-q.swapaxes(0,1))) < 1e-11 and
               np.max(abs(q-q.swapaxes(2,3))) < 1e-11,
               'D2 reference tensor not real/boson-symmetric')
        channels[name] = q
        details[name] = {'coupling': {'real': z.real, 'imag': z.imag},
                         'source_component_entries': {
                            'Re': len(ref['Re_c_unit_tensor_components']),
                            'Im': len(ref['Im_c_unit_tensor_components'])},
                         'Frobenius_original_real_basis': float(np.linalg.norm(q)),
                         'reference_coordinate_basis_is_canonical_SU2_not_certified_RGBeta': True}
    insist(set(channels) == set(EXPECTED_D.values()), 'D2 extra HHSS loss')
    return channels, details


def matrix_metrics(a: np.ndarray, b: np.ndarray) -> dict:
    """Complex Frobenius diagnostics without fitting away missing counterterms."""
    a, b = np.asarray(a, complex), np.asarray(b, complex)
    insist(a.shape == b.shape == (3, 3) and np.isfinite(a).all() and np.isfinite(b).all(),
           'Incompatible coefficient matrices')
    den = max(float(np.vdot(b,b).real), 1e-90)
    projection = complex(np.vdot(b, a) / den)
    difference = a - b
    return {'A_norm': float(np.linalg.norm(a)),
            'B_norm': float(np.linalg.norm(b)),
            'A_minus_B_norm': float(np.linalg.norm(difference)),
            'relative_to_B': float(np.linalg.norm(difference) /
                                   max(np.linalg.norm(b), 1e-100)),
            'complex_fit_coefficient_diagnostic_only':
                {'real': projection.real, 'imag': projection.imag},
            'fit_is_not_applied': True}


def _soft_response(c: np.ndarray, q: np.ndarray, mass2: np.ndarray,
                   mf: float, ms: float) -> dict:
    from Numerical.joint_scalar_matching.core import (
        project_weinberg, scalar_bubble_matching,
    )
    amplitude = np.einsum('pqijab,ABab->pqijAB', c, q, optimize=True)
    beta, ward_beta = project_weinberg(amplitude, strict=False)
    running = math.log(ms/mf) / LOOP * beta
    loop = scalar_bubble_matching(c, q, mass2, ms, reference_scale=mf)
    bubble, ward_bubble = project_weinberg(loop['amplitude'], strict=False)
    return {'beta': beta, 'leading_log': running, 'bubble': bubble,
            'total': running + bubble,
            'projection_residual_beta': float(ward_beta),
            'projection_residual_bubble': float(ward_bubble)}


def a4_seven_channel_diagnostic(snapshot, boundary) -> dict:
    """Evaluate the *complete seven-channel candidate* master beta at muF.

    This uses a formal identity native->RGBeta heavy-scalar basis. It is NOT
    substituted for Batch9's native-source scalar log, because the identity
    CG mapping and absolute signed GammaCC factor are not established.
    """
    from Numerical.running.Aminus4EFT1Couplings import (
        BETA_NAMES, F_LABELS, S_LABELS, tensor_point,
    )
    from Numerical.running.T3EFT1MasterSourceBridge import embed_formal_source
    from Numerical.parallel_t3.wilson.adm import beta, SYMMETRIC_TERMS
    from RGE.running.intermediate.GeneralEFT1SevenTermFlow import (
        TensorBasis, c5_projection,
    )
    f = snapshot.running.mu_f
    s = snapshot.running.mu_s
    native = snapshot.trajectory.at('eft1', f)
    point = tensor_point(f, {name: native[name] for name in BETA_NAMES})
    full, bridge = embed_formal_source(
        boundary.raw, model='A_-4', fermion_labels=F_LABELS,
        scalar_labels=S_LABELS, heavy_source_from_target=None,
        spinor_signed_factor=None,
    )
    gauge = tuple((float(point.couplings[name]), tuple(zip(
        point.gauge[kind]['fermion'], point.gauge[kind]['scalar'])))
        for name, kind in (('gY','U1'),('g2','SU2'),('g3','SU3')))
    parts, backend = beta(full, point.Y, point.lambda4, gauge,
                          variant='fonseca_eq22_six_term')
    insist(set(parts) == set(SYMMETRIC_TERMS), 'Seven-channel ADM inventory changed')
    fm = {label: i for i, label in enumerate(F_LABELS)}
    leptons = tuple((fm[f'L[{p},1]'], fm[f'L[{p},2]']) for p in (1,2,3))
    basis = TensorBasis(45, 12, leptons, (0,1,2,3))
    delta = {name: c5_projection(parts[name], basis)*math.log(s/f)/LOOP
             for name in SYMMETRIC_TERMS}
    for name, arr in delta.items():
        insist(np.isfinite(arr).all() and
               np.allclose(arr, arr.T, rtol=1e-9, atol=1e-23),
               'Non-symmetric provisional A4 ADM projected term: ' + name)
    return {'status': 'SevenChannelCandidateProjected_FixedOrder_NotPhysical',
            'variant': 'fonseca_eq22_six_term',
            'backend': backend, 'model': 'A_-4',
            'formal_bridge': bridge,
            'numerical_channels': {k: matrix_json(v) for k,v in delta.items()},
            'projected_total': matrix_json(sum(delta.values(), np.zeros((3,3),complex))),
            'projected_total_Frobenius': float(np.linalg.norm(sum(delta.values(),np.zeros((3,3),complex)))),
            'includes_only_O_hbar_evolution_of_tree_source': True,
            'gamma_to_weyl_certified': False,
            'matchete_to_rgbeta_scalar_CG_basis_certified': False,
            'crossed_term_for_selected_variant_is_explicit_zero_candidate': True,
            'complete_full_general_D2_ADM_not_claimed': True}


def run_one(root: Path, model: str, output: Path, *,
            include_a4_seven: bool = True,
            illustrative_smeft: bool = False,
            ) -> dict:
    """Use saved successful Batch9 + frozen source data; write only Batch10."""
    from Numerical.parallel_t3.complete.hard_llhh_source import numerical_unsubtracted_source
    from Numerical.joint_scalar_matching.large_batch import load_background
    from Numerical.running.GeneralT3DynamicFermionWilson import (
        uv_from_json, nonmixed_boundary,
    )
    from Numerical.parallel_t3.background.general_batch1_api import read_json

    root = Path(root).resolve()
    output = Path(output).resolve()
    insist(model in MODELS, 'Unknown model')
    b9_dir = root / B9 / model
    b9_path = b9_dir / 'result.json'
    b9 = read_json(b9_path)
    insist(b9.get('model') == model and
           b9.get('status') == 'SourceBasedLLSSScalarLoopComputed_ConditionalNotPhysicalC5' and
           b9.get('only_scalar_pair_log_and_joint_scalar_bubble_computed') is True and
           b9.get('physical_C5_computed') is False,
           model + ': Batch9 successful conditional scalar calculation required')
    snap = load_snapshot(root, model)
    insist(b9.get('file_sha256') == snap.provenance,
           model + ': Batch9 source provenance differs from current Batch1/2 sources')
    uv = uv_from_json(model, snap.uv_values, mu_f=snap.running.mu_f)
    boundary, f_basis = nonmixed_boundary(model, uv, repo_root=root,
                                          mu_f=snap.running.mu_f)
    insist(boundary.raw.shape == (3,3,2,2,8,8),
           model + ': F-tree tensor has unexpected scalar dimension')
    st = b9['LLSS_response_per_unit_unknown_signed_GammaCC_Weyl_factor']
    baseline_log = matrix(st['C5_EFT1_induced_leading_log_per_N_gev_inverse'], 'B9 log')
    baseline_bubble = matrix(st['C5_joint_scalar_finite_bubble_per_N_gev_inverse'], 'B9 bubble')
    baseline = matrix(st['C5_LLSS_total_one_loop_per_N_gev_inverse'], 'B9 sum')
    insist(np.allclose(baseline, baseline_log + baseline_bubble,
                       rtol=2e-9, atol=1e-23), 'Batch9 double-count/missing scalar contribution')
    bg_name = 'A_m4' if model == 'A_-4' else 'D_m2'
    bg = load_background(b9_dir/(bg_name+'.json'), boundary)
    insist(bg.model == model and math.isclose(bg.mu_f, 1e5) and
           math.isclose(bg.mu_s, 1e3), 'Scalar threshold background does not match F source')
    output.mkdir(parents=True, exist_ok=True)
    extra_channels = {}
    extra_outputs = {}
    extra_reference_hash = None
    d_extra_record = {'status': 'NotApplicable_A4_NoExtraHHSS'}
    if model == 'D_-2':
        source = root/EXTRA_D
        extra_reference_hash = digest(source)
        registry = root/REGISTRY
        ref = read_json(source)
        reg = read_json(registry)
        extra, d_details = extra_hhss_from_source(ref, reg, snap.running.couplings)
        from Numerical.joint_scalar_matching.special_mass_v6 import build_special_basis
        basis = build_special_basis(model, snap.running, boundary)
        R = np.asarray(basis['R_original_from_boundary'], dtype=float)
        insist(R.shape == (8,8) and np.allclose(R.T@R, np.eye(8), rtol=1e-11, atol=1e-11),
               'D2 original->split-real rotation invalid')
        insist(np.allclose(basis['mass2_boundary'], bg.mass2, rtol=1e-11, atol=1e-7),
               'D2 split-real scalar mass matrix has drifted since Batch9')
        for name, q_old in extra.items():
            q_new = np.einsum('ABab,ai,bj->ABij', q_old, R, R, optimize=True)
            out = _soft_response(boundary.raw, q_new, bg.mass2, bg.mu_f, bg.mu_s)
            extra_channels[name] = q_new
            extra_outputs[name] = out
        if extra_channels:
            np.savez_compressed(output/'D_m2_extra_HHSS_reference_channels.npz',
                                **extra_channels)
        d_extra_record = {
            'status': 'TwoExtraHHSSReferenceChannelsComputed_CouplingMapUncertified',
            'explicit_HHSS_names': sorted(extra_channels),
            'additional_allowed_non_HHSS_scalar_quartics_not_in_this_one_loop_kernel':
                sorted(set(snap.extras_nonzero)-set(extra_channels)),
            'reference_tensor_original_basis_sha256': extra_reference_hash,
            'registry_sha256': digest(registry),
            'source_details': d_details,
            'reused_split_real_mass_basis': True,
            'source_to_basis_rotation_status': basis['mass_basis_assumption'],
            'native_RGBeta_extra_quartic_normalization_verified': False,
            'additional_HHSS_projection_residuals': {
                name: {'beta': float(out['projection_residual_beta']),
                       'bubble': float(out['projection_residual_bubble'])}
                for name,out in extra_outputs.items()},
        }
    zero = np.zeros((3,3), complex)
    extra_beta = sum((z['beta'] for z in extra_outputs.values()), zero.copy())
    extra_log = sum((z['leading_log'] for z in extra_outputs.values()), zero.copy())
    extra_bubble = sum((z['bubble'] for z in extra_outputs.values()), zero.copy())
    extra_total = extra_log + extra_bubble
    scalar_total = baseline + extra_total

    raw_ledger = read_json(snap.b2_folder/'one_loop_operator_ledger.json')
    at_f = read_json(snap.b1_folder/'fermion_matching_input.json')
    hard = numerical_unsubtracted_source(
        model, raw_ledger, at_f['native_uv_rgbeta_couplings_at_muF'],
        mu_f=snap.running.mu_f)
    h_regular = symmetric(matrix(hard['baseline_ordered_holomorphic_regular_GeVinv'],
                                 'F-hard regular', symmetric=False))
    h_pole = symmetric(matrix(hard['baseline_ordered_holomorphic_pole_GeVinv'],
                              'F-hard pole', symmetric=False))
    # Explicitly DO NOT cancel a pole with an unknown source counterterm.
    scalar_beta = matrix(st['C5_LLSS_beta_16pi2_per_N_gev_inverse'], 'B9 beta') + extra_beta
    candidate_soft_pole = -0.5 / LOOP * scalar_beta
    pole_ledger = matrix_metrics(h_pole, -candidate_soft_pole)
    # Formal source composition is a *diagnostic assumption*, not an
    # MSbar finite C5.  Include the raw hard regular once, and each
    # first-order EFT1 log / scalar B0 once. Never add LLSS running again.
    formal_total = h_regular + scalar_total
    for name, a in (('hard regular', h_regular), ('hard pole', h_pole),
                    ('extra LLSS', extra_total), ('scalar total', scalar_total),
                    ('formal composition', formal_total)):
        insist(np.isfinite(a).all() and
               np.allclose(a, a.T, rtol=2e-9, atol=1e-23),
               name + ': non-symmetric or nonfinite')

    channels_report = None
    if model == 'A_-4' and include_a4_seven:
        try:
            channels_report = a4_seven_channel_diagnostic(snap, boundary)
            a_7 = matrix(channels_report['numerical_channels']['scalar_pair'],
                         'A4 beta scalar pair')
            channels_report['scalar_pair_vs_B9_native_source_log'] = matrix_metrics(
                a_7, baseline_log)
            channels_report['comparison_not_used_to_redefine_B9_or_fit_CG_factor'] = True
        except (ValueError, RuntimeError, FloatingPointError,
                ImportError, TypeError, KeyError) as exc:
            # The source-candidate seven-channel ADM is unverified. A failure
            # must remain visible but cannot invalidate independently obtained
            # D HHSS, A scalar B0, or native F-hard source contributions.
            channels_report = {
                'status': 'CandidateSevenChannelProjectionBlocked',
                'reason': type(exc).__name__ + ': ' + str(exc),
                'replaced_by_scalar_pair_proxy': False,
                'F_hard_and_scalar_sources_unaffected': True,
            }

    record = {
        'format': FORMAT, 'model': model,
        'status': 'HardAndGeneralScalarSourceComponentsComputed_ConditionalNotPhysicalC5',
        'order': [['F'], ['S1', 'S2']],
        'scales_GeV': {'UV':1e7, 'F':1e5, 'S':1e3, 'low':100},
        'F_hard': hard,
        'F_hard_holomorphic_Weinberg_formal_symmetrized_regular': matrix_json(h_regular),
        'F_hard_holomorphic_Weinberg_formal_symmetrized_pole': matrix_json(h_pole),
        'F_hard_pole_vs_scalar_beta_pole': pole_ledger,
        'pole_compensation_not_assumed': True,
        'D2_extra_HHSS': d_extra_record,
        'extra_HHSS_channels': {
            name: {'beta_per_N': matrix_json(out['beta']),
                   'leading_log_per_N': matrix_json(out['leading_log']),
                   'joint_scalar_bubble_per_N': matrix_json(out['bubble']),
                   'scalar_total_per_N': matrix_json(out['total']),
                   'weinberg_projection_residual_beta':out['projection_residual_beta'],
                   'weinberg_projection_residual_bubble':out['projection_residual_bubble']}
            for name,out in extra_outputs.items()},
        'C5_B9_baseline_scalar_per_N_GeVinv': matrix_json(baseline),
        'C5_D2_additional_HHSS_scalar_per_N_GeVinv': matrix_json(extra_total),
        'C5_new_total_scalar_per_N_GeVinv': matrix_json(scalar_total),
        'C5_FORMAL_hard_plus_scalar_unsubtracted_source_sum_GeVinv': matrix_json(formal_total),
        'formal_total_Frobenius':float(np.linalg.norm(formal_total)),
        'F_hard_raw_pole_remains':True,
        'physical_full_MSbar_C5_computed':False,
        'source_spinor_factor_fixed':False,
        'matchete_RGBeta_all_coupling_normalizations_verified':False,
        'D2_nonHHSS_general_quartics_retained_in_running_not_in_this_scalar_pair_vertex':True,
        'A4_full_seven_channel_candidate_at_muF':channels_report,
        'D2_full_seven_channel_tensor_beta_evaluated':False,
        'full_general_18_real_fit_performed':False,
        'source_basis_meta':f_basis,
        'provenance_sha256':{
            'B9_result':digest(b9_path),
            'B9_background':digest(b9_dir/(bg_name+'.json')),
            'F_operator_ledger':digest(snap.b2_folder/'one_loop_operator_ledger.json'),
            'F_boundary':digest(snap.b1_folder/'fermion_matching_input.json'),
            **({'D2_extra_HHSS_reference':extra_reference_hash} if model=='D_-2' else {})
        },
        'unresolved_physics':[
            'Source-exact UV-minus-EFT1-soft finite fermion hard matching and counterterms',
            'Absolute signed GammaCC-to-Weyl normalization, native/RGBeta invariant map',
            'D2 full general scalar one-loop threshold operators beyond LLSS bubble',
            'Complete dimensional-five operator closure, including non-tree Wilson sources',
            'Genuine recomputed-trajectory 18-real-UV-Yukawa fit',
        ],
    }
    if illustrative_smeft:
        # This output is NOT a physical prediction. It is a useful numerical
        # sensitivity exercise under the explicit *unverified* N_Gamma=+1
        # identification of the raw hard and scalar source expressions.
        from Numerical.parallel_t3.complete.full_general_completion import (
            existing_evidence, finish_smeft,
        )
        config = json.loads((root/'configs/thesis_numerics/t3_thesis_A4_Dm2_pilot.json')
                            .read_text(encoding='utf-8-sig'))
        target = json.loads((root/config['fitting']['target'])
                            .read_text(encoding='utf-8-sig'))
        try:
            ev = existing_evidence(model, root/B1_SUBDIR/model, root/B2_SUBDIR/model,
                                   root/LAURENT/model)
            result = finish_smeft(ev, formal_total, vev=float(config['vev_GeV']),
                                  target=target)
            record['illustrative_SMEFT_and_neutrino_observables'] = result
            record['illustrative_low_energy_status'] = 'ComputedNotPhysical'
        except (ValueError, RuntimeError, FloatingPointError) as exc:
            # Do not throw away successful hard/soft decomposition when the
            # OPTIONAL phenomenology of this unverified source point fails.
            record['illustrative_low_energy_status'] = 'BlockedNotPhysical'
            record['illustrative_low_energy_reason'] = type(exc).__name__ + ': ' + str(exc)
        record['illustrative_assumption'] = (
            'UNVERIFIED: holomorphic raw GammaCC hard regular and Batch9/10 '
            'scalar source share a signed factor N_GammaCC=+1; raw one-loop '
            'pole/subtraction and other operators have NOT been completed. '
            'These masses/chi2 are numerical examples, NOT physical predictions.')
    return record
