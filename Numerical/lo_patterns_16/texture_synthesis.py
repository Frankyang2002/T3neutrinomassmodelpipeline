"""Chat A Batch 5 — source-conditional Yukawa-tensor controls and flavour phasors.

Read-only consumer of Batches 1, 3, 4 and their frozen source inputs.
Does not call the production UV solver, alter a benchmark or certify physical
normalisation. All generated files go to a new, protected Batch-5 directory.
"""
from __future__ import annotations

import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

BASE = Path('output/general_observable/lo_pattern_study16')
STUDY = BASE / 'provisional_01'
BATCH03 = BASE / 'analysis_batch03_01'
BATCH04 = BASE / 'analysis_batch04_01'
TARGET = BASE / 'analysis_batch05_01'
REFERENCES = ('A_-2', 'C_-1', 'C_+1', 'E_-2')
ZERO = 1e-100


class SynthesisError(ValueError):
    """Invalid model provenance, tensor identity or analytical result."""


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise SynthesisError(msg)


def read(path: Path) -> dict:
    require(path.is_file(), 'Required source file missing: ' + str(path))
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), 'Source must be a JSON object: ' + str(path))
    return obj


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def decode(value: Any, label: str, shape: tuple[int, ...] | None = None) -> np.ndarray:
    require(isinstance(value, dict) and set(('real', 'imag')) <= set(value), label + ': invalid complex encoding')
    re = np.asarray(value['real'], dtype=float)
    im = np.asarray(value['imag'], dtype=float)
    require(re.shape == im.shape and (shape is None or re.shape == shape), label + ': bad array shape')
    require(np.all(np.isfinite(re)) and np.all(np.isfinite(im)), label + ': nonfinite entries')
    return re + 1j * im


def encode(value: np.ndarray | complex) -> dict:
    z = np.asarray(value, dtype=complex)
    return {'real': z.real.tolist(), 'imag': z.imag.tolist()}


def f_norm(mat: np.ndarray) -> float:
    return float(np.linalg.norm(mat))


def relative(a: np.ndarray, b: np.ndarray) -> float:
    return f_norm(a - b) / max(f_norm(a), f_norm(b), ZERO)


def singular_masses(matrix_ev: np.ndarray) -> list[float]:
    require(matrix_ev.shape == (3, 3) and np.all(np.isfinite(matrix_ev)), 'Invalid neutrino matrix')
    return np.sort(np.linalg.svd(matrix_ev, compute_uv=False)).tolist()


def mode_summary(exact: np.ndarray, groups: list[dict], model: str) -> dict:
    """Re-use tested Batch 4 singular-mode decomposition, not sum of group masses."""
    from Numerical.lo_patterns_16.flavour_geometry import geometry
    data = geometry(exact * 1e9,
                    [{'matrix_ev': g['matrix'] * 1e9,
                      'mass_gev': math.sqrt(g['m2']),
                      'multiplicity': g['dimension']} for g in groups], model)
    light = data['modes'][0]
    return {'m1_ev': data['masses_ev'][0], 'm2_ev': data['masses_ev'][1],
            'm3_ev': data['masses_ev'][2],
            'm1_over_m3': data['masses_ev'][0] / max(data['masses_ev'][-1], ZERO),
            'm2_over_m3': data['masses_ev'][1] / max(data['masses_ev'][-1], ZERO),
            'm1_geometry_A': light['geometric_projection_fraction'],
            'm1_coherence_C': light['coherence_fraction'],
            'm1_mechanism': light['mechanism_diagnostic'],
            'group_strength_ev': data['fermi_group_strength_ev'],
            'phasors': data['projections'][:len(groups)],
            'max_group_matrix_norm_ev': max(g['matrix_frobenius_ev'] for g in data['heavy_groups'])}


def _fitted_yukawas(frozen: dict, model: str) -> tuple[np.ndarray, np.ndarray]:
    require('base_state' in frozen and 'ordinary' in frozen['base_state'],
            model + ': missing frozen benchmark ordinary inputs')
    p = frozen['base_state']['ordinary']
    return decode(p['y1'], model + ' y1', (3, 3)), decode(p['y2'], model + ' y2', (3, 3))


def _extract_channel_tensor(y: np.ndarray, vertex: np.ndarray, label: str) -> tuple[np.ndarray, dict]:
    """Find CG/Weyl tensor C in vertex[i,...]=sum_r conj(y[i,r]) C[r,...].

    A valid full-rank frozen texture is required: under rank deficiency the
    channel tensor cannot be identified from the frozen vertices alone.
    """
    require(y.shape == (3, 3) and vertex.ndim == 3 and vertex.shape[0] == 3,
            label + ': incompatible frozen Yukawa and vertex dimensions')
    require(np.all(np.isfinite(y)) and np.all(np.isfinite(vertex)), label + ': nonfinite input')
    A = np.conj(y)
    sv = np.linalg.svd(A, compute_uv=False)
    require(sv[-1] > 1e-10 * sv[0] and sv[0] > 0,
            label + ': frozen Yukawa rank/condition prevents an identifiable CG tensor')
    C = np.linalg.solve(A, vertex.reshape(3, -1)).reshape((3,) + vertex.shape[1:])
    reconstructed = np.einsum('ir,rfs->ifs', A, C)
    error = relative(reconstructed, vertex)
    require(error < 5e-9, label + f': vertex does not factorize as fitted Yukawa × CG ({error:g})')
    return C, {'relative_reconstruction_error': error,
               'yukawa_condition_number': float(sv[0] / sv[-1]),
               'tensor_norm': f_norm(C)}


def extract_source_tensors(inp: dict, frozen: dict, model: str) -> dict:
    """Extract tensor maps, respecting the checked conjugate/field orientation."""
    require(inp.get('model') == model, model + ': input model label mismatch')
    vertices = inp.get('vertices')
    require(isinstance(vertices, list) and len(vertices) == 2, model + ': unexpected vertex count')
    require(vertices[0].get('name') == 'y1_HC' and vertices[1].get('name') == 'y2_HC'
            and 'conjugate' in vertices[0] and 'field' not in vertices[0]
            and 'field' in vertices[1] and 'conjugate' not in vertices[1],
            model + ': unexpected Weyl vertex orientation; cannot infer tensors')
    y1, y2 = _fitted_yukawas(frozen, model)
    z = decode(vertices[0]['conjugate'], model + ' Z')
    y = decode(vertices[1]['field'], model + ' Y', z.shape)
    c1, diag1 = _extract_channel_tensor(y1, z, model + ' y1->Z')
    c2, diag2 = _extract_channel_tensor(y2, y, model + ' y2->Y')
    return {'C1': c1, 'C2': c2, 'y1': y1, 'y2': y2,
            'tensor_reconstruction': {'y1_HC': diag1, 'y2_HC': diag2}}


def substitute_textures(inp: dict, tensors: dict, y1: np.ndarray, y2: np.ndarray) -> dict:
    """Controlled Yukawa substitution, all other candidate Weyl inputs fixed."""
    require(y1.shape == y2.shape == (3, 3), 'Controlled Yukawas must be 3x3')
    require(np.all(np.isfinite(y1)) and np.all(np.isfinite(y2)), 'Nonfinite counterfactual Yukawa')
    z = np.einsum('ir,rfs->ifs', np.conj(y1), tensors['C1'])
    y = np.einsum('ir,rfs->ifs', np.conj(y2), tensors['C2'])
    v = [{'name': 'y1_HC', 'conjugate': encode(z)},
         {'name': 'y2_HC', 'field': encode(y)}]
    return {**inp, 'vertices': v}


def scalar_mixing_scale(inp: dict, d1: int, factor: float) -> dict:
    """Scale only the S1/S2 off-diagonal A blocks, preserving Hermiticity."""
    require(math.isfinite(factor) and factor > 0, 'Scalar mixing scale must be finite and positive')
    A = decode(inp['scalar_A_gev2'], 'scalar A')
    require(A.ndim == 2 and A.shape[0] == A.shape[1] and 0 < d1 < len(A),
            'Scalar A dimensions do not match model')
    B = A.copy()
    B[:d1, d1:] *= factor
    B[d1:, :d1] *= factor
    require(relative(B, B.conj().T) < 1e-10, 'Scaled scalar mass matrix not Hermitian')
    return {**inp, 'scalar_A_gev2': encode(B)}


def kernel_eval(inp: dict, model: str, d1: int) -> dict:
    from Numerical.lo_patterns_16.analytic_explanation import reconstruct
    return reconstruct(inp, model[0], int(model[2:]), d1)


def phasor_explanation(summary: dict, expected: dict, model: str) -> dict:
    """Explain observed lightest mode using **complex** heavy-group amplitudes.

    No phases or directions are defined for heavy vectors within degenerate
    mass eigenspaces; only summed eigenspace matrices are analyzed.
    """
    parts = summary['phasors']
    complex_a = [complex(x['amplitude_ev']['real'], x['amplitude_ev']['imag']) for x in parts]
    observed = complex(sum(complex_a))
    m1 = summary['m1_ev']
    require(abs(observed - m1) < 2e-7 * max(summary['group_strength_ev'], m1, ZERO),
            model + ': complex mode phasors do not sum to m1')
    amplitudes = sum(abs(z) for z in complex_a)
    require(abs(summary['m1_coherence_C'] - m1 / amplitudes) < 3e-7
            if amplitudes > 1e-17 * summary['group_strength_ev'] else summary['m1_coherence_C'] is None,
            model + ': coherence does not reproduce phasor sum')
    pairs = []
    for j in range(len(complex_a)):
        for k in range(j + 1, len(complex_a)):
            z, w = complex_a[j], complex_a[k]
            cos = float((z * w.conjugate()).real / max(abs(z) * abs(w), ZERO))
            pairs.append({'g': j + 1, 'h': k + 1,
                          'phase_cosine': max(-1., min(1., cos)),
                          'destructive_cross_term_ev2': 2 * (z * w.conjugate()).real,
                          'pair_residual_fraction': abs(z + w) / max(abs(z) + abs(w), ZERO)})
    pairs.sort(key=lambda p: p['destructive_cross_term_ev2'])
    return {'model': model, 'mode': 1, 'm1_ev': m1,
            'phasor_sum_ev': encode(observed),
            'sum_absolute_projected_ev': amplitudes,
            'cancelled_fraction': max(0., 1. - m1 / max(amplitudes, ZERO)),
            'phasors': [{'heavy_group': i + 1,
                         'heavy_mass_gev': float(parts[i]['mass_gev']),
                         'amplitude_ev': encode(z), 'magnitude_ev': abs(z),
                         'relative_phase_deg': float(math.degrees(math.atan2(z.imag, z.real))),
                         'fraction_of_abs_sum': abs(z) / max(amplitudes, ZERO)}
                        for i, z in enumerate(complex_a)],
            'strongest_destructive_pair': pairs[0] if pairs else None,
            'all_pairs': pairs,
            'original_Batch4_coherence': expected['lightest_mode']['coherence_fraction']}


def _resolve_local(root: Path, recorded: str) -> Path:
    """Use the local checkout, not an obsolete absolute path from a past run."""
    p = Path(recorded)
    if p.is_absolute():
        require(p.is_relative_to(root), 'Source path points outside this checkout: ' + str(p))
        return p
    q = (root / p).resolve()
    require(q.is_relative_to(root), 'Relative source path escapes checkout')
    return q


def _write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def _write_csv(path: Path, records: list[dict], fields: list[str]) -> None:
    with path.open('w', newline='', encoding='utf-8') as fp:
        w = csv.DictWriter(fp, fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(records)


def report(payload: dict) -> str:
    allrows = payload['models']
    lookup = {r['model']: r for r in allrows}
    blocks = [
        '# Chat A — Batch 5: analytical origin of the 16 neutrino-mass patterns', '',
        '**Scope:** frozen, source-conditional T3 Weyl amplitudes. All B1/B2 physical normalization '
        'is unresolved. No EFT matching, no refit, no production UV-solver calls.', '',
        '## Exact component and flavour decomposition', '',
        'The saved candidate Weyl vertices are factorized into **frozen Yukawas** and '
        'model-specific component tensors:', '',
        r'\[Z_{iIs}=\sum_r y_{1,ir}^{*} C^{(1)}_{rIs},\qquad '
        r'Y_{iIt}=\sum_r y_{2,ir}^{*} C^{(2)}_{rIt}.\]', '',
        'The complete source-conditional mass is rebuilt as the sum of charge and distinct '
        'heavy singular-mass eigenspace terms:', '',
        r'\[M_\nu=\frac{1}{16\pi^2}\sum_{q,g}\sum_{k\in g}'
        r'K_q(m_k^2)\operatorname{sym}'
        r'\bigl[Z_{s_q}(M_F^\dagger P_k)Y_{t_q}^{T}\bigr],\]', '',
        r'\[K_q(m^2)=\delta_q\frac{B_0(m^2,\lambda_{q+})-B_0(m^2,\lambda_{q-})}'
        r'{\lambda_{q+}-\lambda_{q-}}.\]', '',
        'The scalar formula is exact for the restricted isolated 2×2 S1/S2 charge blocks '
        'already verified by Batch 3. It is not a claim about missing B1/B2 phases or symmetry factors.', '',
        '## Frozen physical input versus controlled *counterfactual* inputs', '',
        'For each model we leave its SU(2), hypercharge, heavy Weyl matrix, scalar masses '
        'and candidate component tensors fixed, but substitute exactly the same '
        '3×3 complex Yukawa matrices from four reference frozen benchmarks: '
        '`A_-2`, `C_-1`, `C_+1`, `E_-2`. Each replacement is performed at the '
        '**vertex level**, with every conjugation preserved. None is a refitted benchmark. '
        'A separate common-|lambda| comparison rescales E scalar mixing from 0.05 to 0.10; '
        'this is only a controlled source-conditional comparison.', '',
        '## Observed mass pattern and complex-mode cancellation', '',
        '| Model | Fermion | m1/m3 | A1 | C1 | cancelled projected strength | strongest opposing pair |',
        '|---|---|---:|---:|---:|---:|---|',
    ]
    for r in allrows:
        o = r['observed']
        p = r['phasor_explanation']
        pair = p['strongest_destructive_pair']
        pairtext = ('n/a' if pair is None else
                    f'{pair["g"]}–{pair["h"]} (cos={pair["phase_cosine"]:.3g})')
        blocks.append(f'| {r["model"]} | {r["fermion_type"]} | {o["m1_over_m3"]:.5g} | '
                      f'{o["m1_geometry_A"]:.5g} | {o["m1_coherence_C"]:.5g} | '
                      f'{p["cancelled_fraction"]:.3%} | {pairtext} |')
    blocks.extend(['', '## Controlled reference textures, common scalar mixing parameter', '',
                   'For this comparison, each reference provides **both** `y1` and `y2`. '
                   'The ratio m1/m3 is invariant under an overall nonzero matrix scaling '
                   'but not under a change of texture or Majorana/Dirac construction. '
                   'All entries use |lambdaT3|=0.10, rescaling E by two.', '',
                   '| Model | own frozen m1/m3 | A_-2 texture | C_-1 texture | C_+1 texture | E_-2 texture |',
                   '|---|---:|---:|---:|---:|---:|'])
    for r in allrows:
        v=r['common_mixing_controlled']
        blocks.append('| ' + r['model'] + ' | ' + f'{r["observed"]["m1_over_m3"]:.5g}'
                      + ' | ' + ' | '.join(f'{v[x]["m1_over_m3"]:.5g}' for x in REFERENCES) + ' |')
    blocks.extend(['', '## C triplet family: controlled texture swaps and interpolation', '',
                   'C_-3 and C_+1 use charged Dirac heavy fermions; C_-1 has a neutral '
                   'Majorana heavy fermion. The following Yukawa changes do **not** alter '
                   'this fundamental difference.', '',
                   '| Frozen model | substitution | m1/m3 | A1 | C1 | m1 (eV) |',
                   '|---|---|---:|---:|---:|---:|'])
    for m in ('C_-3','C_-1','C_+1'):
        r=lookup[m]
        for label,val in r['C_family_swaps'].items():
            blocks.append(f'| {m} | {label} | {val["m1_over_m3"]:.5g} | '
                          f'{val["m1_geometry_A"]:.5g} | {val["m1_coherence_C"]:.5g} | '
                          f'{val["m1_ev"]:.6g} |')
    blocks.extend(['', '### Measured C-family heavy-eigenspace phasors', '',
                   'Each row is the **complex** contribution projected on the full lightest-neutrino SVD mode. '
                   'The sum is a real positive `m1` within numerical accuracy; the individual '
                   'terms need not be positive.', '',
                   '| Frozen model | heavy group | mass [GeV] | Re(a_g1) [eV] | Im(a_g1) [eV] | angle [deg] |',
                   '|---|---:|---:|---:|---:|---:|'])
    for model in ('C_-3','C_-1','C_+1'):
        for ph in lookup[model]['phasor_explanation']['phasors']:
            amp=ph['amplitude_ev']
            blocks.append(f'| {model} | {ph["heavy_group"]} | {ph["heavy_mass_gev"]:.7g} | '
                          f'{amp["real"]:.8g} | {amp["imag"]:.8g} | '
                          f'{ph["relative_phase_deg"]:.5g} |')
    blocks.extend(['', '### C_-1 versus C_+1: algebraic suppression factors', '',
                   r'\[m_1=G A_1 C_1,\quad G=\sum_g\|H_g\|_F,\quad '
                   r'A_1=\frac{\sum_g|a_{g1}|}{G},\quad '
                   r'C_1=\frac{|\sum_g a_{g1}|}{\sum_g|a_{g1}|}.\]', ''])
    q=lookup['C_-1']['observed'];p=lookup['C_+1']['observed']
    blocks.extend([
       f'The frozen ratio m1(C_-1)/m1(C_+1) = {q["m1_ev"]/max(p["m1_ev"],ZERO):.8g} '
       f'factorizes as [{q["group_strength_ev"]/max(p["group_strength_ev"],ZERO):.6g}] '
       f'(overall group strength) × '
       f'[{q["m1_geometry_A"]/max(p["m1_geometry_A"],ZERO):.6g}] '
       f'(mode projection) × '
       f'[{q["m1_coherence_C"]/max(p["m1_coherence_C"],ZERO):.6g}] '
       '(heavy-eigenspace coherence).',
       'This is an *exact numerical factorization*, not an independent theory '
       'prediction or a causal attribution to Majorana versus Dirac structure.', ''])
    for model,other in (('C_+1','C_-1'),('C_-1','C_+1')):
        path=[v for v in payload['C_texture_interpolations'] if v['model']==model]
        best=min(path,key=lambda z:z['m1_over_m3'])
        blocks.append(f'- {model} → {other} fitted texture: frozen m1/m3 = '
                      f'{path[0]["m1_over_m3"]:.7g}; endpoint = '
                      f'{path[-1]["m1_over_m3"]:.7g}; smallest **sampled** ratio '
                      f'{best["m1_over_m3"]:.7g} at t={best["t"]:.1f}.')
    blocks.extend(['', 'The CSV `batch05_C_texture_interpolation.csv` evaluates 11 '
                   'intermediate Yukawa pairs in each direction. Its extrema and the '
                   'mode projection/coherence diagnostics reveal whether the C_+1 '
                   'cancellation persists after the C_-1 Yukawa texture is substituted. '
                   'A small mass at an intermediate point is a **counterfactual** algebraic '
                   'cancellation and is not a fitted or naturalness statement.', '',
                   '## Pairwise A/D, B, C, E patterns', '',
                   'To avoid mistaking a fit for a group-theory identity, both the original '
                   'frozen mass-matrix proximity and the controlled-common-texture proximity '
                   'are reported in `batch05_matched_pairs.csv`. A near-perfect alignment '
                   'of the former alone does not establish representation equivalence. '
                   'E has two aligned scalar-charge channels in the saved inputs; this is '
                   'separate from the three heavy-eigenspace interference phases.', '',
                   '### Controlled matrix proportionality tests', '',
                   'A small residual below means both *complex* 3x3 mass matrices '
                   'are related by a single fitted complex scalar for the indicated '
                   'shared source Yukawa texture. It tests a conditional benchmark relation, '
                   'not an operator-level SU(2) identity.', '',
                   '| Pair | frozen residual | common A_-2 | common C_-1 | common C_+1 | common E_-2 |',
                   '|---|---:|---:|---:|---:|---:|'])
    for p in payload['matched_pairs']:
        l=p['model_a'];r=p['model_b']
        blocks.append(f'| {l} / {r} | {p["frozen_matrix_proportional_residual"]:.3g} | '
                      + ' | '.join(f'{p["common_reference_matrix_shapes"][t]["proportional_matrix_residual"]:.3g}'
                                   for t in REFERENCES) + ' |')
    blocks.extend(['', '## Model-specific explanations', ''])
    for r in allrows:
        o=r['observed'];p=r['phasor_explanation'];z=r['tensor_reconstruction']
        blocks.extend([f'### {r["model"]}: {r["fermion_type"]}', '',
                       f'- Source masses: ({o["m1_ev"]:.7g}, {o["m2_ev"]:.7g}, '
                       f'{o["m3_ev"]:.7g}) eV. m1/m3={o["m1_over_m3"]:.5g}; '
                       f'candidate charge channels={r["active_charge_channels"]}.',
                       f'- Lightest mode: A1={o["m1_geometry_A"]:.5g}, '
                       f'C1={o["m1_coherence_C"]:.5g}; '
                       f'{p["cancelled_fraction"]:.3%} of absolute projected strength cancels.',
                       f'- Original source vertices reconstructed with relative errors '
                       f'{z["y1_HC"]["relative_reconstruction_error"]:.2e} and '
                       f'{z["y2_HC"]["relative_reconstruction_error"]:.2e}.',
                       f'- At common scalar mixing, C_+1 reference texture gives '
                       f'm1/m3={r["common_mixing_controlled"]["C_+1"]["m1_over_m3"]:.5g}; '
                       f'C_-1 reference texture gives '
                       f'm1/m3={r["common_mixing_controlled"]["C_-1"]["m1_over_m3"]:.5g}.',
                       '- The heavy eigenspace amplitudes (magnitudes and phases) are '
                       'recorded in `batch05_mode_phasors.csv`; the reference-texture '
                       'counterfactuals isolate flavour input changes only.'])
        blocks.append('')
    blocks.extend(['## Interpretation constraints', '',
                   '- A common Yukawa **matrix** is not a universal physical effective coupling '
                   'across different SU(2) representations. The source-conditional group '
                   'coefficients remain embedded in the extracted component tensors.',
                   '- Identical heavy singular-mass eigenvalues do not imply identical '
                   'fermion propagator numerators or Majorana/Dirac contractions.',
                   '- Source Yukawas used to extract the tensors must be identifiable '
                   '(full-rank); ambiguous tensors trigger an error rather than being inferred.',
                   '- Phasor angles refer to projections onto singular vectors of the *full* '
                   'Mν. They are not physical CP phases and should not be compared as '
                   'generation eigenvector phases within degenerate heavy subspaces.',
                   '- None of the controlled substitutions is used to certify the source '
                   'Lagrangian, and no EFT calculation is performed.', ''])
    return '\n'.join(blocks)


def analyse(root: Path, *, study_dir: Path | None = None,
            batch03_dir: Path | None = None, batch04_dir: Path | None = None,
            output_dir: Path | None = None, make_plots: bool = True) -> dict:
    """Validate full source lineage and produce Batch 5 without overwriting files."""
    from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary
    from Numerical.lo_patterns_16.result_analysis import _source_hash_guard

    root=root.resolve()
    owner=(root / BASE).resolve()
    paths={'study':(study_dir or root/STUDY).resolve(),
           'batch03':(batch03_dir or root/BATCH03).resolve(),
           'batch04':(batch04_dir or root/BATCH04).resolve(),
           'output':(output_dir or root/TARGET).resolve()}
    for key,p in paths.items():
        require(p.is_relative_to(owner) and p!=owner,
                key+': not a contained Chat A output path')
    require(len(set(paths.values()))==4, 'Source and target directories overlap')
    target=paths['output']
    require(not target.exists() or (target.is_dir() and not any(target.iterdir())),
            'Output directory already populated: '+str(target))

    summary=read(paths['study'] / 'all16_summary.json')
    b3=read(paths['batch03'] / 'batch03_analytic_explanation.json')
    m3=read(paths['batch03'] / 'SOURCE_MANIFEST.json')
    b4=read(paths['batch04'] / 'batch04_flavour_analysis.json')
    m4=read(paths['batch04'] / 'SOURCE_MANIFEST.json')
    require(summary.get('format')=='t3-conditional-pattern-study-all16-v1', 'Batch1 format drift')
    require(b3.get('format')=='t3-chatA-batch03-analytic-scalar-kernel-v1', 'Batch3 format drift')
    require(b4.get('format')=='t3-chatA-batch04-flavour-geometry-v1', 'Batch4 format drift')
    require(summary.get('physical_normalisation_certified') is False and
            b3.get('physical_normalisation_certified') is False and
            b4.get('physical_normalisation_certified') is False and
            m3.get('physical_source_certified') is False and
            m4.get('physical_normalisation_certified') is False,
            'Unexpected physical certification flags')
    require(summary['counts']['conditional_evaluated']==16 and
            summary['counts']['blocked']==0 and summary['counts']['physical_source_certified']==0 and
            b3['model_count']==16 and b4['model_count']==16,
            'Incomplete conditional all-16 source reports')
    sha1=digest(paths['study']/'all16_summary.json')
    sha3=digest(paths['batch03']/'batch03_analytic_explanation.json')
    sha4=digest(paths['batch04']/'batch04_flavour_analysis.json')
    require(b3.get('source_summary_sha256')==sha1 and
            b4['source_sha256']['source_run_sha256']==sha1 and
            b4['source_sha256']['source_batch03_sha256']==sha3 and
            m3.get('source_summary_sha256')==sha1 and
            m4.get('source_batch03_sha256')==sha3,
            'Source summary/B3/B4 lineage SHA mismatch')
    catalogue=model_dictionary()
    names=[item['model'] for item in catalogue]
    require(len(names)==16 and len(set(names))==16 and
            [r['model'] for r in summary['models']]==names and
            [r['model'] for r in b3['models']]==names and
            [r['model'] for r in b4['models']]==names,
            'Missing/reordered model catalogue')
    bymeta={x['model']:x for x in catalogue}
    by3={x['model']:x for x in b3['models']}
    by4={x['model']:x for x in b4['models']}
    sources={};textures={};records=[];phasor_rows=[];control_rows=[];interpolation_rows=[]
    controlled_matrix_cache={}
    source_manifest={'format':'t3-chatA-batch05-source-manifest-v1',
                     'study_sha256':sha1,'batch03_sha256':sha3,'batch04_sha256':sha4,
                     'source_revision':summary['provenance']['revision'],
                     'physical_normalisation_certified':False,'models':{}}
    for model in names:
        base=read(paths['study']/'models'/f'{model}.json')
        original=by3[model]
        expected=by4[model]
        src=base['source']
        require(base.get('model')==model and base.get('status')=='conditional_uv_evaluated' and
                base.get('physical_normalisation_certified') is False and
                original.get('physical_normalisation_certified') is False,
                model+': model status changed')
        require(digest(paths['study']/'models'/f'{model}.json')==
                m3['models'][model]['model_detail_sha256']==
                m4['models'][model]['batch1_model_sha256'] and
                src['benchmark_sha256']==m3['models'][model]['benchmark_sha256'] and
                src['input_sha256']==m4['models'][model]['candidate_weyl_input_sha256'],
                model+': model source hash drift')
        _source_hash_guard(root, src, summary, model)
        input_path=_resolve_local(root,src['input'])
        bench_path=_resolve_local(root,src['benchmark'])
        require(digest(input_path)==src['input_sha256'] and
                digest(bench_path)==src['benchmark_sha256'],
                model+': source file SHA mismatch')
        inp=read(input_path)
        frozen=read(bench_path)
        d1=bymeta[model]['SU2_dimensions']['S1']
        extraction=extract_source_tensors(inp,frozen,model)
        benchmark_y1,benchmark_y2=extraction['y1'],extraction['y2']
        check=substitute_textures(inp,extraction,benchmark_y1,benchmark_y2)
        require(relative(decode(check['vertices'][0]['conjugate'],'Z'),
                         decode(inp['vertices'][0]['conjugate'],'Z'))<5e-9 and
                relative(decode(check['vertices'][1]['field'],'Y'),
                         decode(inp['vertices'][1]['field'],'Y'))<5e-9,
                model+': fitted tensor reconstruction failed')
        calculated=kernel_eval(check,model,d1)
        expected_matrix=decode(base['mnu_gev'],model+' baseline mnu',(3,3))
        reported_matrix=decode(read(paths['batch03']/'models'/f'{model}.json')
                               ['exact_closed_kernel_mnu_gev'], model+' Batch3 mnu',(3,3))
        require(relative(calculated['exact'],expected_matrix)<2e-7 and
                relative(calculated['exact'],reported_matrix)<2e-7,
                model+': extracted tensor baseline disagrees with saved mass matrix')
        # Regrouping invariant: the eigenspace-wise matrices must match Batch 3.
        claimed=original['heavy_mass_groups']
        got=calculated['heavy_groups']
        require(len(got)==len(claimed), model+': heavy group count drift')
        for x,y in zip(got,claimed):
            require(abs(math.sqrt(x['m2'])-float(y['fermion_mass_gev']))<
                    3e-7*max(float(y['fermion_mass_gev']),ZERO) and
                    x['dimension']==y['eigenspace_multiplicity'] and
                    relative(x['matrix']*1e9,decode(y['matrix_ev'],'saved heavy group',(3,3)))<2e-7,
                    model+': heavy group source/analytical contract drift')
        observed=mode_summary(calculated['exact'],calculated['heavy_groups'],model)
        require(abs(observed['m1_over_m3']-float(expected['m1_over_m3']))<1e-6 and
                abs(observed['m1_coherence_C']-float(expected['lightest_mode']['coherence_fraction']))<2e-6,
                model+': Batch4 lightest-mode coherence drift')
        phasors=phasor_explanation(observed,expected,model)
        require(abs(phasors['original_Batch4_coherence']-observed['m1_coherence_C'])<2e-6,
                model+': phasor coherence source drift')
        sources[model]={'input':inp,'d1':d1,'frozen':frozen,'tensors':extraction,
                        'baseline_mnu':calculated['exact']}
        textures[model]=(benchmark_y1,benchmark_y2)
        entry={'model':model,'family':model[0],'alpha':int(model[2:]),
               'fermion_type':bymeta[model]['fermion_type'],
               'active_charge_channels':sum(1 for x in calculated['channels'] if x['matrix_norm_gev']>1e-17*f_norm(calculated['exact'])),
               'observed':{k:v for k,v in observed.items() if k!='phasors'},
               'phasor_explanation':phasors,
               'tensor_reconstruction':extraction['tensor_reconstruction'],
               'common_mixing_controlled':{},'C_family_swaps':{},
               'conditional_only':True,'physical_normalisation_certified':False}
        records.append(entry)
        phasor_rows.extend({'model':model,**{k:v for k,v in ph.items() if k!='amplitude_ev'},
                            'real_ev':ph['amplitude_ev']['real'],'imag_ev':ph['amplitude_ev']['imag'],
                            'observed_m1_ev':phasors['m1_ev'],
                            'observed_coherence':observed['m1_coherence_C']}
                           for ph in phasors['phasors'])
        source_manifest['models'][model]={'benchmark_sha256':src['benchmark_sha256'],
                                          'candidate_weyl_input_sha256':src['input_sha256'],
                                          'study_model_sha256':digest(paths['study']/'models'/f'{model}.json')}
    # Deliberately wait until every frozen source is verified before using
    # other models' fitted textures as counterfactual inputs.
    controlled_mixing = .10
    for rec in records:
        model=rec['model'];base=sources[model]
        controlled_matrix_cache[model]={}
        fit_y1,fit_y2=textures[model]
        for name in REFERENCES:
            y1,y2=textures[name]
            modified=substitute_textures(base['input'],base['tensors'],y1,y2)
            native_lam=base['frozen']['base_state']['ordinary']['lambdaT3']
            native_lam=abs(complex(native_lam['real'],native_lam['imag'])) if isinstance(native_lam,dict) else abs(complex(native_lam))
            require(native_lam>0,model+': zero reference lambdaT3')
            altered=scalar_mixing_scale(modified,base['d1'],controlled_mixing/native_lam)
            m=kernel_eval(altered,model,base['d1'])
            datum=mode_summary(m['exact'],m['heavy_groups'],model)
            controlled_matrix_cache[model][name]=m['exact']
            rec['common_mixing_controlled'][name]={k:v for k,v in datum.items() if k!='phasors'}
            control_rows.append({'model':model,'texture_reference':name,
                                 'lambdaT3_abs':controlled_mixing,
                                 **{k:v for k,v in datum.items() if k not in ('phasors','m1_mechanism')}})
        if model[0]=='C':
            recipes={'frozen':(fit_y1,fit_y2),
                     'C_-1_both':textures['C_-1'],
                     'C_+1_both':textures['C_+1'],
                     'C_-1_y1_only':(textures['C_-1'][0],fit_y2),
                     'C_-1_y2_only':(fit_y1,textures['C_-1'][1]),
                     'C_+1_y1_only':(textures['C_+1'][0],fit_y2),
                     'C_+1_y2_only':(fit_y1,textures['C_+1'][1])}
            for label,(y1,y2) in recipes.items():
                result=kernel_eval(substitute_textures(base['input'],base['tensors'],y1,y2),
                                   model,base['d1'])
                meas=mode_summary(result['exact'],result['heavy_groups'],model)
                rec['C_family_swaps'][label]={k:v for k,v in meas.items() if k!='phasors'}
    # Continuous interpolation in both directions, preserving each model's
    # fermion nature and original scalar/Higgs/CG inputs.
    for model,other in (('C_+1','C_-1'),('C_-1','C_+1')):
        source=sources[model]
        y1a,y2a=textures[model];y1b,y2b=textures[other]
        for step in range(11):
            t=step/10.
            yy1=(1-t)*y1a+t*y1b
            yy2=(1-t)*y2a+t*y2b
            result=kernel_eval(substitute_textures(source['input'],source['tensors'],yy1,yy2),
                               model,source['d1'])
            datum=mode_summary(result['exact'],result['heavy_groups'],model)
            interpolation_rows.append({'model':model,'target_texture':other,'t':t,
                **{k:v for k,v in datum.items() if k not in ('phasors','m1_mechanism')}})
    pairs=[]
    matched=[('A_-4','D_+2'),('A_-2','D_+0'),('A_+0','D_-2'),
             ('B_-3','B_+1'),('C_-3','C_+1'),('E_-4','E_+2'),('E_-2','E_+0')]
    indexed={x['model']:x for x in records}
    for left,right in matched:
        a=indexed[left];b=indexed[right]
        x=sources[left]['baseline_mnu'];y=sources[right]['baseline_mnu']
        pa=(np.vdot(x,y)/max(np.vdot(x,x).real,ZERO))
        original_rel=relative(pa*x,y)
        # Apply common reference texture to BOTH models; this isolates the
        # conditional source representation/fermion/scalar effects.
        common_ratio={name: b['common_mixing_controlled'][name]['m1_over_m3'] /
                      max(a['common_mixing_controlled'][name]['m1_over_m3'],ZERO)
                      for name in REFERENCES}
        common_shapes={}
        for name in REFERENCES:
            p=controlled_matrix_cache[left][name]
            q=controlled_matrix_cache[right][name]
            z=np.vdot(p,q)/max(np.vdot(p,p).real,ZERO)
            common_shapes[name]={'best_complex_scale':encode(z),
                                 'proportional_matrix_residual':relative(z*p,q),
                                 'frob_norm_ratio_b_over_a':f_norm(q)/max(f_norm(p),ZERO)}
        pairs.append({'model_a':left,'model_b':right,
                      'frozen_complex_scale_best_fit':encode(pa),
                      'frozen_matrix_proportional_residual':original_rel,
                      'frozen_ratio_m1_over_m3':b['observed']['m1_over_m3']/max(a['observed']['m1_over_m3'],ZERO),
                      'common_reference_m1_over_m3_ratios':common_ratio,
                      'common_reference_matrix_shapes':common_shapes})
    require(len(records)==16 and len(control_rows)==16*len(REFERENCES) and
            len(interpolation_rows)==22 and len(phasor_rows)>=48,
            'Incomplete all-16 control matrix coverage')
    payload={'format':'t3-chatA-batch05-analytical-yukawa-texture-v1',
             'source_conditional_only':True,'physical_normalisation_certified':False,
             'safe_for_physical_source_certification':False,
             'model_count':16,'reference_texture_models':list(REFERENCES),
             'common_scalar_mixing_abs':controlled_mixing,
             'methods':{'tensor_factorization':'Z[i,I,s]=sum_r conjugate(y1[i,r])*C1[r,I,s]; Y similarly',
                        'mass_kernel':'Batch3 exact finite 2x2 scalar divided difference and invariant heavy eigenspace projector',
                        'counterfactual':'replace only candidate Yukawa vertices; preserve original fermion/scalar/CG inputs',
                        'phasor':'Batch4 u1-dagger H_g v1 complex projected amplitudes',
                        'trajectory':'linear complex Yukawa interpolation; not a fit or RGE'},
             'models':records,'matched_pairs':pairs,
             'C_texture_interpolations':interpolation_rows,
             'source_manifest_sha256':{'study':sha1,'batch03':sha3,'batch04':sha4}}
    text=report(payload)  # All computations completed before creating output.
    target.mkdir(parents=True,exist_ok=True)
    _write_json(target/'batch05_synthesis.json',payload)
    _write_json(target/'SOURCE_MANIFEST.json',source_manifest)
    handoff={'format':'t3-chatA-batch05-handoff-v1',
             'model_count':16,'controlled_reference_scenarios':len(control_rows),
             'C_family_single_texture_swaps':21,'interpolation_points':len(interpolation_rows),
             'physical_source_certified':False,'safe_for_physical_source_certification':False,
             'source_manifest_sha256':source_manifest,
             'primary_result':'complex heavy-state phasors and source-conditioned fixed-texture controls',
             'open_certification':'B1/B2 Matchete/RGBeta-to-Weyl phases and normalization'}
    _write_json(target/'HANDOFF.json',handoff)
    (target/'BATCH05_ANALYTICAL_SYNTHESIS.md').write_text(text,encoding='utf-8')
    _write_csv(target/'batch05_mode_phasors.csv',phasor_rows,
               ['model','heavy_group','heavy_mass_gev','real_ev','imag_ev','magnitude_ev',
                'relative_phase_deg','fraction_of_abs_sum','observed_m1_ev','observed_coherence'])
    _write_csv(target/'batch05_common_texture_controls.csv',control_rows,
               ['model','texture_reference','lambdaT3_abs','m1_ev','m2_ev','m3_ev',
                'm1_over_m3','m2_over_m3','m1_geometry_A','m1_coherence_C','group_strength_ev'])
    _write_csv(target/'batch05_C_texture_interpolation.csv',interpolation_rows,
               ['model','target_texture','t','m1_ev','m2_ev','m3_ev',
                'm1_over_m3','m2_over_m3','m1_geometry_A','m1_coherence_C','group_strength_ev'])
    _write_csv(target/'batch05_matched_pairs.csv',
               [{'model_a':p['model_a'],'model_b':p['model_b'],
                 'frozen_matrix_proportional_residual':p['frozen_matrix_proportional_residual'],
                 'frozen_ratio_m1_over_m3':p['frozen_ratio_m1_over_m3'],
                 **{f'common_{k}_m1ratio_b_over_a':v for k,v in p['common_reference_m1_over_m3_ratios'].items()},
                 **{f'common_{k}_matrix_residual':p['common_reference_matrix_shapes'][k]['proportional_matrix_residual']
                    for k in REFERENCES}}
                for p in pairs],
               ['model_a','model_b','frozen_matrix_proportional_residual',
                'frozen_ratio_m1_over_m3']+
                [f'common_{k}_m1ratio_b_over_a' for k in REFERENCES]+
                [f'common_{k}_matrix_residual' for k in REFERENCES])
    if make_plots:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            pass  # Plotting is optional; JSON/CSV/Markdown are complete.
        else:
            names=[r['model'] for r in records]
            x=np.arange(len(names))
            fig,ax=plt.subplots(figsize=(12,5))
            ax.plot(x,[r['observed']['m1_over_m3'] for r in records], 'o-',label='Frozen Yukawas')
            for label in REFERENCES:
                ax.plot(x,[r['common_mixing_controlled'][label]['m1_over_m3'] for r in records],
                        '.--',label=label+' reference')
            ax.set_xticks(x, names, rotation=60)
            ax.set_yscale('log');ax.set_ylabel(r'$m_1/m_3$');ax.legend(ncol=3,fontsize=8)
            ax.set_title('Frozen versus controlled Yukawa textures (conditional)')
            fig.tight_layout();fig.savefig(target/'batch05_controlled_hierarchy.png',dpi=170);plt.close(fig)
            fig,ax=plt.subplots(figsize=(8,5))
            for model in ('C_-1','C_+1'):
                series=[v for v in interpolation_rows if v['model']==model]
                ax.plot([v['t'] for v in series],[v['m1_over_m3'] for v in series],'.-',label=model)
            ax.set_xlabel('Fraction of opposite C texture');ax.set_ylabel(r'$m_1/m_3$')
            ax.set_yscale('log');ax.legend();ax.set_title('C Yukawa interpolation (conditional)')
            fig.tight_layout();fig.savefig(target/'batch05_C_texture_trajectory.png',dpi=170);plt.close(fig)
            fig,ax=plt.subplots(figsize=(8,5))
            for model in ('C_-3','C_-1','C_+1'):
                pp=indexed[model]['phasor_explanation']['phasors']
                ax.plot([p['amplitude_ev']['real'] for p in pp],
                        [p['amplitude_ev']['imag'] for p in pp],'.-',label=model)
            ax.set_xlabel('Re($a_{g1}$) [eV]');ax.set_ylabel('Im($a_{g1}$) [eV]')
            ax.axhline(0,linewidth=.8);ax.axvline(0,linewidth=.8)
            ax.legend();ax.set_title('C family lightest-mode heavy phasors (conditional)')
            fig.tight_layout();fig.savefig(target/'batch05_C_complex_phasors.png',dpi=170);plt.close(fig)
    return payload
