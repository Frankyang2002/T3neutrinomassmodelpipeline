"""Fail-closed all-16 native Laurent / EFT1-soft subtraction *readiness* audit.

This module does not claim Matchete's F-only MatchedEFT is an un-subtracted
UV amplitude. No native EFT1-soft graphs, counterterms, or physical C5 are
inferred from the algebraic equality of native pole and log coefficients.

The only evaluated EFT1 loop ingredient is the universal zero-momentum
scalar two-point master B0(0;m_a^2,m_b^2), not a full LLHH diagram.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any

from .boundary import BoundaryError, MODELS, require, tag
from .pole_log_bridge import FORMAT as BRIDGE_FORMAT, audit_all, expressions_for
from .raw_llhh import sha_bytes

FORMAT = 't3-f-first-soft-counterterm-readiness-v1'
NATIVE_FORMAT = 't3-f-only-raw-one-loop-native-evidence-v1'
EXPECTED_TOTALS = {'kernels': 176, 'PL': 88, 'PR': 88,
                   'pole': 176, 'finite_constant': 176, 'scale_log': 176}
SCALAR_BUBBLE_FORMULA = (
    'B0(0;a,b;mu2)_MSbar_finite = '
    '1 - [a*Log[a/mu2] - b*Log[b/mu2]]/(a-b); '
    'equal-mass limit = -Log[a/mu2]'
)
SOURCE_SEMANTICS = (
    'Matchete F-only one-loop matched-EFT output; not independently verified '
    'as the un-subtracted UV 1PI amplitude. Do not subtract EFT1 twice.'
)
POLE_CONVENTION = (
    'CONDITIONAL ONLY: if R_native were independently verified as the bare UV '
    'pole residue, an amplitude convention A_UV - A_EFT1_soft + A_CT '
    'would require R_native - R_EFT1_soft + R_CT = 0.'
)


def _read_json(path: Path) -> dict:
    require(path.is_file(), f'Missing JSON: {path}')
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), f'Expected JSON object: {path}')
    return obj


def _inside(root: Path, p: Path) -> Path:
    candidate = p.resolve()
    require(candidate.is_relative_to(root.resolve()), f'Path outside project: {candidate}')
    return candidate


def bubble_b0_finite(a_gev2: float, b_gev2: float, mu2_gev2: float) -> float:
    """Real MSbar finite B0(0;a,b) for nonnegative scalar mass squared.

    Convention: B0 = 1/epsilon-bar + finite, with finite =
    -Integral_0^1 dx ln[(x*a+(1-x)*b)/mu2]. This does NOT specify a
    loop i/(16pi^2) factor, EFT1 contraction or operator sign.
    Stable near a=b via log1p; a=b=0 is infrared singular.
    """
    values = (a_gev2, b_gev2, mu2_gev2)
    require(all(isinstance(v, (float, int)) and not isinstance(v, bool)
                and math.isfinite(v) for v in values), 'Nonfinite/invalid B0 input')
    require(a_gev2 >= 0 and b_gev2 >= 0 and mu2_gev2 > 0,
            'B0 needs nonnegative scalar mass-squared and positive scale-squared')
    a, b = sorted((float(a_gev2), float(b_gev2)), reverse=True)
    require(a > 0, 'B0(0;0,0) has an infrared divergence')
    if b == 0:
        return 1.0 - math.log(a / mu2_gev2)
    r = (a - b) / a
    if r == 0.0:
        return -math.log(a / mu2_gev2)
    if r < 1.0e-6:
        # (1-r)*(-log(1-r))/r = 1 - r/2 - r^2/6 - r^3/12 - ...
        q = 1.0 - r/2.0 - r*r/6.0 - r*r*r/12.0 - r**4/20.0
    else:
        q = (1.0-r) * (-math.log1p(-r)) / r
    return 1.0 - math.log(a / mu2_gev2) - q


def scalar_master_table(mass_sq: list[float], mu_f_gev: float,
                        model: str, input_sha256: str) -> dict:
    require(isinstance(mass_sq, list) and len(mass_sq) == 2,
            'Two scalar mass eigenvalues required')
    require(isinstance(mu_f_gev, (float, int)) and not isinstance(mu_f_gev, bool)
            and math.isfinite(mu_f_gev) and mu_f_gev > 0, 'Invalid mu_F')
    m = [float(x) for x in mass_sq]
    require(all(math.isfinite(x) and x > 0 for x in m),
            'Nonpositive scalar mass eigenvalue squared')
    values = []
    for i, j in ((0, 0), (0, 1), (1, 1)):
        finite = bubble_b0_finite(m[i], m[j], mu_f_gev**2)
        values.append({'scalar_pair': ['X','Y'][i] + ['X','Y'][j],
                       'mass_a_squared_gev2': m[i], 'mass_b_squared_gev2': m[j],
                       'B0_MSbar_finite_at_mu_F': finite,
                       'B0_1_over_epsilon_bar_pole': 1,
                       'd_B0_finite_d_ln_mu_bar_squared': 1})
    return {
        'model': model, 'native_mass_basis_input_sha256': input_sha256,
        'mass_eigen_squared_gev2': m,
        'mu_F_gev': mu_f_gev, 'scalar_pairs': values,
        'evaluated_object': 'Universal scalar B0 master only; not EFT1 operator coefficient',
        'B0_formula': SCALAR_BUBBLE_FORMULA,
        'dimensional_regularization_pole_convention': '+1/epsilon-bar',
        'loop_normalization_and_operator_symmetry_factors_applied': False,
        'EFT1_soft_diagram_calculated': False,
    }


def _validate_mass_basis(root: Path, plan: dict, model: str, mu_f: float) -> dict:
    name = plan.get('mass_basis_input')
    digest = plan.get('mass_basis_input_sha256')
    require(isinstance(name, str) and bool(name) and isinstance(digest, str)
            and re.fullmatch('[a-f0-9]{64}', digest),
            f'{model}: missing matching-verified mass-basis input')
    path = _inside(root, root / Path(name.replace('\\', '/')))
    require(path.is_file() and sha_bytes(path.read_bytes()) == digest,
            f'{model}: mass-basis source hash mismatch')
    data = _read_json(path)
    require(data.get('model') == model and
            data.get('format') == 't3-experimental-mass-eigen-matchete-parameters-v1' and
            data.get('status') == 'NumericMassEigenstateBoundaryPreparedMatcheteNotRun',
            f'{model}: unexpected native mass-basis schema')
    return scalar_master_table(data.get('heavy_mass_squared_gev2'), mu_f, model, digest)


def validate_bridge_snapshot(repo_root: Path, native_dir: Path, inventory_dir: Path,
                             bridge_dir: Path) -> tuple[dict, dict[str, dict], dict]:
    """Recompute upstream from native hashed sources and compare every bridge export."""
    root = repo_root.resolve()
    owned = (root / 'output/general_observable/lo_certification_spinor').resolve()
    for d in (native_dir, inventory_dir, bridge_dir):
        require(d.resolve().is_relative_to(owned) and d.resolve().is_dir(),
                'Inputs must be existing owned spinor output directories')
    require(len({native_dir.resolve(), inventory_dir.resolve(), bridge_dir.resolve()}) == 3,
            'Native, inventory, and bridge directories must be distinct')
    calculated, models = audit_all(root, native_dir, inventory_dir)
    manifest = _read_json(bridge_dir / 'summary.json')
    require(manifest == calculated, 'Committed bridge summary differs from reconstructed native evidence')
    require(calculated.get('totals') == EXPECTED_TOTALS,
            'Unexpected kernel counts in this exact 16-model release')
    for model in MODELS:
        d = bridge_dir / tag(model)
        recorded = _read_json(d / 'trios.json')
        require(recorded == models[model], f'{model}: bridge trios mismatch')
        for chirality in ('PL', 'PR'):
            e = expressions_for(models[model], chirality)
            for kind, expected in e.items():
                path = d / f'{chirality}_{kind}.inputform.txt'
                require(path.is_file() and path.read_text('utf8') == expected + '\n',
                        f'{model}: bridge expression mismatch {chirality}/{kind}')
    native = _read_json(native_dir / 'summary.json')
    require(native.get('format') == NATIVE_FORMAT and
            native.get('n_models_with_raw_native_one_loop_evidence') == 16 and
            native.get('hard_one_loop_MSbar_certified') is False and
            native.get('native_GammaCC_to_master_Weyl_phase') is None,
            'Native status unexpectedly changed')
    return calculated, models, native


def build_readiness(root: Path, bridge_summary: dict, models: dict[str, dict],
                    native_summary: dict, native_sha: str, bridge_sha: str) -> tuple[dict, dict]:
    """Create all-16 subtraction constraints, without inventing soft/CT data."""
    root = root.resolve()
    plans = native_summary.get('model_plans')
    require(isinstance(plans, list) and len(plans) == 16 and
            {x.get('model') for x in plans} == set(MODELS),
            'Native plan lacks the 16 model inputs')
    by_model = {p['model']: p for p in plans}
    mu_f = bridge_summary.get('mu_F_gev')
    require(type(mu_f) in (int, float) and math.isfinite(mu_f) and mu_f > 0,
            'Invalid fermion threshold')
    scalar_masters: dict[str, Any] = {}
    ledgers: dict[str, Any] = {}
    counts = Counter()
    for model in MODELS:
        source = models[model]
        records = []
        for t in source['records']:
            chirality = t['chirality']
            require(chirality in ('PL', 'PR'), 'Invalid native chirality')
            kernel = t['kernel_inputform']
            require(sha256(kernel.encode('utf8')).hexdigest() == t['kernel_sha256'],
                    'Modified native pole kernel')
            records.append({
                'chirality': chirality,
                'native_kernel_sha256': t['kernel_sha256'],
                'native_pole_residue_inputform': kernel,
                'native_constant_inputform': kernel,
                'native_mass_log_inputform': t['mass_log_inputform'],
                'source_indices': t['source_indices'],
                'raw_pole_residue_ratio_in_native_kernel_basis': '1',
                'raw_scale_log_derivative_d_ln_mu_bar_squared_ratio': '1',
                'external_EFT1_soft_pole_residue_inputform': None,
                'external_EFT1_soft_finite_inputform': None,
                'native_counterterm_pole_inputform': None,
                'native_counterterm_finite_inputform': None,
                'EFT1_soft_operator_and_CG_contraction_verified': False,
                'native_matching_vs_bare_UV_amplitude_semantics_verified': False,
                'counterterm_MSbar_scheme_verified': False,
                'physical_hard_one_loop_wilson': None,
            })
            counts['kernels'] += 1
            counts[chirality] += 1
        ledgers[model] = {
            'format': FORMAT, 'model': model,
            'status': 'NATIVE_POLE_RESIDUES_KNOWN_SOFT_AND_CT_UNCOMPUTED',
            'source_inventory_terms_sha256': source['source_inventory_terms_sha256'],
            'native_kernel_count': len(records), 'records': records,
            'physical_MSbar_hard_matching': None,
            'physical_C5': None,
        }
        if model in ('B_-1', 'C_-1'):
            scalar_masters[model] = _validate_mass_basis(root, by_model[model], model, mu_f)
        else:
            require(by_model[model].get('mass_basis_input') is None,
                    f'{model}: unexpected scalar input ignored')
            scalar_masters[model] = {
                'model': model, 'status': 'SYMBOLIC_SCALAR_MASTER_ONLY',
                'scalar_masses_squared_gev2': None,
                'scalar_masses_from_native_run_available': False,
                'B0_formula': SCALAR_BUBBLE_FORMULA,
                'd_B0_finite_d_ln_mu_bar_squared': 1,
                'EFT1_soft_diagram_calculated': False,
            }
    require(dict(counts) == {'kernels': 176, 'PL': 88, 'PR': 88},
            'Unexpected 16-model pole population')
    payload = {
        'format': FORMAT, 'status': 'ALL16_NATIVE_POLE_AND_SCALAR_BUBBLE_MASTERS_READY_SOFT_CT_BLOCKED',
        'mu_F_gev': mu_f, 'models': 16, 'totals': dict(counts),
        'source_release_sha256': bridge_summary['source_release_summary_sha256'],
        'source_native_summary_sha256': native_sha,
        'source_pole_log_bridge_summary_sha256': bridge_sha,
        'native_loop_source_interpretation': SOURCE_SEMANTICS,
        'conditional_pole_balance_equation': POLE_CONVENTION,
        'scalar_bubble_basis': 'B0(0,m_a^2,m_b^2) with a,b active at mu_F',
        'scalar_bubble_MSbar_formula': SCALAR_BUBBLE_FORMULA,
        'mixed_models_with_hashed_mass_eigenvalues': ['B_-1','C_-1'],
        'other_models_scalar_mass_inputs_not_in_native_batch': 14,
        'only_all_BSM_fermions_integrated_at_mu_F': True,
        'both_BSM_scalars_remain_light_in_EFT1': True,
        'SM_fermions_remain_light': True,
        'integrate_both_BSM_scalars_together_only_at_later_mu_S': True,
        'native_one_loop_source_already_soft_subtracted': None,
        'native_CT_witness_available': False,
        'EFT1_soft_loop_matching_amplitude_computed': False,
        'physical_MSbar_hard_matching_certified': False,
        'absolute_matchete_gamma_cc_to_master_weyl_factor': None,
        'physical_C5': None,
        'physical_SMEFT_scalar_threshold_performed': False,
        'needs': [
            'Identify Matchete Match[UV,LoopOrder->1] hard/soft and CT semantics; do NOT double-subtract',
            'Compute signed, SU2-contracted EFT1 soft LLHH loops with both scalars active',
            'Compute native counterterm and pole cancellation in one scheme',
            'Check wavefunction, coupling and mixed-mass-basis contributions as required',
            'Fix independent Matchete GammaCC to master two-component Weyl dictionary',
        ],
    }
    return payload, {'ledgers': ledgers, 'scalar_masters': scalar_masters}
