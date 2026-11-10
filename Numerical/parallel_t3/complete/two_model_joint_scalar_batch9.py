"""A(-4)/D(-2) integrated, source-derived *conditional* LLSS scalar response.

This fills a real numerical handoff missing from the thesis-specific full-general
Batch1/Batch2 workflow.  It combines the saved, continuously run native UV/EFT1
couplings with source-exact F-tree LLSS and the canonical joint-scalar B0 kernel.
It is deliberately NOT a finite physical one-loop C5: the F-hard boundary,
complete EFT1 ADM, full scalar-soft inventory and signed GammaCC/Weyl map are
not supplied by this operation.

No GitHub writes.  All model-specific numerical files come from the local repo.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any
import json
import math

import numpy as np

MODELS = ('A_-4', 'D_-2')
B1_SUBDIR = Path('output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1')
B2_SUBDIR = Path('output/parallel_t3/full_general_batch2/thesis_full_general_batch2_raw_v1')
SYMBOLIC_SUBDIR = Path('output/full16_symbolic')
# These archived tree seeds are byte-identical to the current thesis Batch2
# F-tree source for A_-4 and D_-2, as inspected in repository main.
# Verify again on the user's local checkout before evaluating their tensors.
ORIGINAL_F_TREE_SEED = {
    'A_-4': Path('output/general_observable/fermion_first_all16_v1/A_m4/data/eft1_after_F_wilson_seed.json'),
    'D_-2': Path('output/general_observable/fermion_first_self_fix_v2/D_m2/data/eft1_after_F_wilson_seed.json'),
}
FORMAT = 't3-thesis-source-derived-llss-joint-scalar-batch9-v1'
UNCOMPUTED = (
    'Finite renormalised full one-loop F-hard Wilson matching (UV minus EFT1 soft)',
    'Complete seven-channel EFT1 Wilson evolution (only scalar-pair leading log here)',
    'Complete scalar-threshold one-loop hard/soft and inherited operator inventory',
    'Absolute signed Matchete GammaCC-to-Weyl normalisation',
    'Full-general physical C5(100 GeV), Takagi masses, PMNS, NuFIT fit',
)


class IntegrationBlocked(ValueError):
    """Missing or incompatible source data; no synthetic result is substituted."""


def insist(ok: bool, explanation: str) -> None:
    if not ok:
        raise IntegrationBlocked(explanation)


def sha(path: Path) -> str:
    insist(path.is_file(), f'Missing source: {path}')
    return sha256(path.read_bytes()).hexdigest()


def real_scalar(value: Any, name: str, *, positive: bool = False) -> float:
    """Convert explicitly finite real RGBeta scalar, never discard an imaginary part."""
    a = np.asarray(value, dtype=np.complex128)
    insist(a.shape == () and np.isfinite(a).all(), f'{name}: expected finite scalar')
    z = complex(a)
    insist(abs(z.imag) <= 1e-12 * max(1.0, abs(z.real)),
           f'{name}: non-negligible imaginary component {z.imag}')
    x = float(z.real)
    insist(math.isfinite(x) and (not positive or x > 0),
           f'{name}: non-finite or non-positive scalar')
    return x


def quartic_source_path(root: Path, model: str) -> Path:
    """Original (unsplit) native HHSS CG source, NOT the D split-field seed.

    D's split-field Batch2 quartic export omits unsplit NewScalar2 channels;
    source_quartic_v4 instead converts the original-field Matchete quartics
    into the split named matching basis with the existing real rotation.
    """
    insist(model in MODELS, f'Unsupported model {model}')
    family, alpha = model.split('_')
    tag = ('m' + alpha[1:]) if alpha.startswith('-') else ('p' + alpha[1:])
    return root / SYMBOLIC_SUBDIR / f'T3_{family}_alpha_{tag}' / 'data' / 'eft1_after_F_scalar_quartic_seed.json'


@dataclass(frozen=True)
class SourceSnapshot:
    model: str
    trajectory: Any
    uv_values: dict[str, Any]
    running: Any
    b1_folder: Path
    b2_folder: Path
    boundary_seed: Path
    quartic_seed: Path
    provenance: dict[str, str]
    extras_nonzero: dict[str, dict[str, float]]


def load_snapshot(root: Path, model: str) -> SourceSnapshot:
    """Read exact saved scale endpoints with checked Batch1 input file hashes."""
    from Numerical.parallel_t3.background.general_batch1_api import load_batch1
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import RunningSnapshot

    root = Path(root).resolve()
    insist(model in MODELS, f'Unsupported model {model}')
    b1 = root / B1_SUBDIR / model
    b2 = root / B2_SUBDIR / model
    trajectory = load_batch1(b1)
    insist(trajectory.model == model, f'{model}: Batch1 model mismatch')
    scales = trajectory.manifest['scales_gev']
    insist(scales == {'UV': 1e7, 'F': 1e5, 'S': 1e3, 'low': 100},
           f'{model}: unexpected RG scale hierarchy')
    at_f = trajectory.fermion_matching_input()
    at_s = trajectory.scalar_matching_input()
    insist(at_f.get('all_heavy_fermions_to_integrate_out_together') is True and
           at_f.get('both_heavy_scalar_multiplets_remain_light') is True and
           at_s.get('both_heavy_scalar_multiplets_to_integrate_out_together') is True,
           f'{model}: threshold decoupling order is inconsistent')
    uv_values = at_f['native_uv_rgbeta_couplings_at_muF']
    from Numerical.parallel_t3.background.general_batch1_api import decode_native
    for key in ('MF', 'y1', 'y2'):
        insist(key in uv_values, f'{model}: missing F-endpoint {key}')
        a = np.asarray(decode_native(uv_values[key]), dtype=np.complex128)
        insist(a.shape == (3, 3) and np.isfinite(a).all(),
               f'{model}: bad 3x3 F-endpoint {key}')
    eft = trajectory.at_scalar_threshold()
    native_at_s = at_s.get('native_eft1_rgbeta_couplings_at_muS', {})
    insist(set(native_at_s) == set(eft),
           f'{model}: saved scalar handoff and RG endpoint inventories differ')
    for name in eft:
        observed = np.asarray(decode_native(native_at_s[name]), dtype=np.complex128)
        insist(np.allclose(observed, eft[name], rtol=1e-13, atol=1e-11),
               f'{model}: scalar endpoint {name} differs from source trajectory')

    # RunningSnapshot supplies values in the native RGBeta coordinate system.
    # Complex phases are retained. The Matchete coordinate *identification* is
    # the existing v5 candidate map and must not be reported as established.
    numeric: dict[str, Any] = {}
    for name, value in eft.items():
        a = np.asarray(value, dtype=np.complex128)
        if a.shape != ():
            continue   # SM Yukawa matrices are not scalar couplings
        z = complex(a)
        insist(np.isfinite([z.real, z.imag]).all(),
               f'{model}: nonfinite scalar RG coordinate {name}')
        numeric[name] = ({'real': float(z.real), 'imag': float(z.imag)}
                         if abs(z.imag) > 1e-12 * max(1.0, abs(z.real))
                         else float(z.real))
    m1 = real_scalar(eft['mS1Sq'], 'mS1Sq', positive=True)
    m2 = real_scalar(eft['mS2Sq'], 'mS2Sq', positive=True)
    if model == 'D_-2':
        insist('muT3GeneralS2Sq' in numeric and 'y2Conj' in uv_values,
               'D(-2): missing scalar self-bilinear or extra Yukawa')
        extras = sorted(k for k in eft if k.startswith('lambdaT3GeneralQ'))
        insist(len(extras) == 5, 'D(-2): expected five independently evolved extra quartics')
    else:
        insist(not any(k.startswith('lambdaT3GeneralQ') for k in eft),
               'A(-4): unexpected general quartics')
    extras_nonzero = {}
    for name in sorted(k for k in eft if k.startswith('lambdaT3GeneralQ')):
        z = complex(np.asarray(eft[name]))
        if z != 0:
            extras_nonzero[name] = {'real': float(z.real), 'imag': float(z.imag)}

    running = RunningSnapshot(
        model=model, source=b1/'scalar_matching_input.json',
        sha256=sha(b1/'scalar_matching_input.json'),
        mu_f=float(scales['F']), mu_s=float(scales['S']),
        scalar_mass1_sq=m1, scalar_mass2_sq=m2,
        couplings=numeric,
        evidence={'source': 'Batch1 general UV/EFT1 RGBeta trajectory at actual muS',
                  'native_matchete_scalar_coupling_aliases_certified': False},
    )
    native_status = b2/'native_matching/native_status.json'
    native_seed = b2/'native_matching/data/eft1_after_F_wilson_seed.json'
    for path in (native_status, native_seed):
        insist(path.is_file(), f'{model}: missing saved full-general F match {path}')
    status = json.loads(native_status.read_text(encoding='utf-8-sig'))
    doc = json.loads(native_seed.read_text(encoding='utf-8-sig'))
    insist(status.get('model') == model and status.get('integrated_fields') == ['F']
           and status.get('heavy_flavour_count') == 3,
           f'{model}: full-general raw F-stage status invalid')
    insist(doc.get('IntegratedFields') == ['F'] and
           doc.get('ActiveHeavyFields') == ['S1', 'S2'] and
           doc.get('TreeWilsonTermCount') == (2 if model == 'A_-4' else 8),
           f'{model}: full-general F-tree operator inventory changed')
    original_tree_seed = root / ORIGINAL_F_TREE_SEED[model]
    insist(sha(native_seed) == sha(original_tree_seed),
           f'{model}: archived F-tree native CG source differs from current full-general thesis Batch2')
    original_seed = quartic_source_path(root, model)
    provenance = {p.name + ':' + str(i): sha(p) for i, p in enumerate((
        b1/'manifest.json', b1/'fermion_matching_input.json',
        b1/'scalar_matching_input.json', native_status, native_seed,
        original_tree_seed, original_seed
    ))}
    return SourceSnapshot(model, trajectory, uv_values, running,
                          b1, b2, native_seed, original_seed,
                          provenance, extras_nonzero)


def execute_model(root: Path, model: str, output: Path) -> dict[str, Any]:
    """Calculate existing source-exact F-tree -> scalar pair RG + joint B0.

    This is a fixed-order O(hbar) *subset*. No physical C5 or neutrino mass
    is returned because the F-hard subtracted coefficient and other one-loop
    channels are not yet implemented in this combined two-model evaluator.
    """
    from Numerical.running.GeneralT3DynamicFermionWilson import (
        uv_from_json, nonmixed_boundary,
    )
    from Numerical.joint_scalar_matching.source_quartic_v4 import build_source_channels
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import build_ordinary_background
    from Numerical.joint_scalar_matching.special_mass_v6 import build_special_background
    from Numerical.joint_scalar_matching.large_batch import (
        evaluate_conditional, load_background, parse_matrix,
    )

    root = Path(root).resolve()
    output = Path(output).resolve()
    snap = load_snapshot(root, model)
    uv = uv_from_json(model, snap.uv_values, mu_f=snap.running.mu_f)
    boundary, f_info = nonmixed_boundary(model, uv, repo_root=root,
                                         mu_f=snap.running.mu_f)
    insist(boundary.model == model and boundary.raw.shape[:4] == (3,3,2,2),
           f'{model}: invalid source-exact F-tree tensor')
    insist(boundary.raw.shape[-2:] == (8,8),
           f'{model}: eight heavy real scalar components required')
    channels = build_source_channels(snap.quartic_seed, model)
    insist(channels.model == model,
           f'{model}: Matchete source HHSS quartics model mismatch')
    output.mkdir(parents=True, exist_ok=True)
    if model == 'A_-4':
        background = build_ordinary_background(channels, boundary,
                                               snap.running, output)
    else:
        background = build_special_background(channels, boundary,
                                              snap.running, output)
    bg = load_background(Path(background['background_json']), boundary)
    response = evaluate_conditional(boundary, bg)
    k = parse_matrix(response['C5_LLSS_total_one_loop_per_N_gev_inverse'])
    insist(np.isfinite(k).all() and
           np.allclose(k, k.T, rtol=1e-9, atol=1e-23),
           f'{model}: non-symmetric conditional scalar-pair response')
    return {
        'format': FORMAT,
        'status': 'SourceBasedLLSSScalarLoopComputed_ConditionalNotPhysicalC5',
        'model': model,
        'order': [['F'], ['S1', 'S2']],
        'scales_GeV': {'UV': 1e7, 'F': 1e5, 'S': 1e3, 'low': 100},
        'source_F_tree_wilson_terms': 2 if model == 'A_-4' else 8,
        'source_F_basis': f_info,
        'native_general_quartics_nonzero_at_S_not_discarded_from_RGE': snap.extras_nonzero,
        'scalar_background': background,
        'LLSS_response_per_unit_unknown_signed_GammaCC_Weyl_factor': response,
        'LLSS_response_Frobenius_GeV_inverse': float(np.linalg.norm(k)),
        'only_scalar_pair_log_and_joint_scalar_bubble_computed': True,
        'full_EFT1_seven_term_Wilson_RG_computed': False,
        'F_hard_one_loop_matching_included': False,
        'full_one_loop_scalar_matching_included': False,
        'physical_C5_computed': False,
        'neutrino_masses_computed': False,
        'full_general_18_real_fit_computed': False,
        'source_native_Matchete_to_RGBeta_mapping_certified': False,
        'known_missing_physics': list(UNCOMPUTED),
        'file_sha256': snap.provenance,
        'original_source_quartic_file': str(snap.quartic_seed.relative_to(root)),
        'F_wilson_seed_file': str(snap.boundary_seed.relative_to(root)),
    }
