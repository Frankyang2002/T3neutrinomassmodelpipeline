#!/usr/bin/env python3
"""Read-only Workstream A UV closure stage for 16 ordinary T3 points.

Independent Sym^n(SU2) singlet counting for *all* scalar monomials of degree
2/3/4, exact comparison with GeneralScalarVerifiedQuartics.wl, numerical
quadratic and generic interaction-tensor transport for all 8 exceptional
mass bases, and optionally native Matchete free/interaction term Z2 checks for
all 16.  Prior eigenbasis outputs are consumed as evidence, never generated.

The certificate is LIMITED: this does not implement a Matchete discrete-Z2
registry, verify every possible Lorentz fermion operator, equate symbolic
Matchete spurions with UV masses, or certify RGBeta matching normalization.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement, product
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

WL = ROOT / 'Lagrangian' / 'RunT3UVZ2TermAudit.wl'
CATALOG = ROOT / 'Lagrangian' / 'interactions' / 'GeneralScalarVerifiedQuartics.wl'
SELFDIR = ROOT / 'Reports' / 'output' / 'workstream_a' / 'six_self_eigenbasis'
BCMF_DIR = ROOT / 'Reports' / 'output' / 'workstream_a'

DIMS = {'A': (1, 3, 2), 'B': (2, 2, 1), 'C': (2, 2, 3),
        'D': (3, 1, 2), 'E': (3, 3, 2)}
ALPHAS = {'A': (-4, -2, 0), 'B': (-3, -1, 1), 'C': (-3, -1, 1),
          'D': (-2, 0, 2), 'E': (-4, -2, 0, 2)}
MODELS = tuple(f'{family}_{"m" if a < 0 else "p"}{abs(a)}'
               for family, aa in ALPHAS.items() for a in aa)
SELF = tuple(f'{f}_{tag}' for f in ('A', 'D', 'E') for tag in ('m2', 'p0'))
BC = ('B_m1', 'C_m1')
SPECIES = ('H', 'Hdag', 'S1', 'S1dag', 'S2', 'S2dag')
BASELINE = (
    (0, 0, 2, 2, 0, 0), (0, 0, 0, 0, 2, 2),
    (1, 1, 1, 1, 0, 0), (1, 1, 0, 0, 1, 1),
    (0, 0, 1, 1, 1, 1), (2, 0, 1, 0, 0, 1),
)
FLAGS = ('T3_EXPERIMENTAL_FLAVOR_MASSES',
         'T3_EXPERIMENTAL_GENERAL_FERMIONS',
         'T3_EXPERIMENTAL_GENERAL_QUARTICS')
Z2_SUCCESS = 'PassedActualBSMTermParityAndMatchete'
SELF_SUCCESS = 'MatcheteSelfEigenHeavyBasisPassedMatchingNotPerformed'
BCMF_SUCCESS = 'PassedFinalEigenbasisIndexedMFAndMatchete'


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise ValueError(msg)


def name_alpha(model: str) -> tuple[str, int]:
    f, tag = model.split('_')
    return f, (-1 if tag.startswith('m') else 1) * int(tag[1:])


def conjugate(counts: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(counts[i] for i in (1, 0, 3, 2, 5, 4))


def canonical(counts: tuple[int, ...]) -> tuple[int, ...]:
    return min(counts, conjugate(counts))


BASELINE_KEYS = frozenset(canonical(c) for c in BASELINE)


@lru_cache(None)
def sym_weights(dim: int, degree: int) -> Counter[int]:
    """SU2 weights in Sym^degree(V_dim), including Bose statistics."""
    weights = range(1 - dim, dim, 2)
    return Counter(sum(choice) for choice in combinations_with_replacement(weights, degree))


@lru_cache(None)
def singlet_multiplicity(dims: tuple[int, ...],
                         counts: tuple[int, ...]) -> int:
    weights = Counter({0: 1})
    for dim, degree in zip(dims, counts):
        if not degree:
            continue
        next_weights: Counter[int] = Counter()
        for w, a in weights.items():
            for v, b in sym_weights(dim, degree).items():
                next_weights[w + v] += a * b
        weights = next_weights
    # For SU2, multiplicity of J=0 equals multiplicity of weight 0 minus
    # multiplicity of weight 2 (weights measured as twice T3).
    result = weights[0] - weights[2]
    require(result >= 0, 'Negative SU2 singlet multiplicity')
    return result


def scalar_operators(model: str) -> dict[int, dict[tuple[int, ...], int]]:
    family, alpha = name_alpha(model)
    ds1, ds2, _ = DIMS[family]
    dims = (2, 2, ds1, ds1, ds2, ds2)
    twice_y = (1, -1, alpha, -alpha, alpha + 2, -alpha - 2)
    results = {}
    for degree in (2, 3, 4):
        out = {}
        for choice in combinations_with_replacement(range(6), degree):
            counts = tuple(choice.count(i) for i in range(6))
            if not sum(counts[2:]) or sum(counts[2:]) % 2:
                continue  # no BSM fields or excluded by exact Z2
            if sum(n * y for n, y in zip(counts, twice_y)):
                continue
            key = canonical(counts)
            if key in out:
                continue
            number = singlet_multiplicity(dims, key)
            if number:
                out[key] = number
        results[degree] = out
    return results


def parse_wolfram_manifest(path: Path) -> dict[str, dict[tuple[int, ...], int]]:
    """Narrow, fail-closed parser for the source's single-line WL table rows."""
    content = path.read_text(encoding='utf-8-sig')
    start = content.index('T3ExtendedVerifiedQuarticTable = <|')
    end = content.index('|>;', start)
    body = content[start:end]
    rows = {}
    for line in body.splitlines():
        if '->' not in line:
            continue
        match = re.fullmatch(r'\s*"(\d+_\d+_-?\d+)"\s*->\s*(\{.*\}),?\s*', line)
        require(match is not None, 'Unexpected Wolfram manifest row: ' + line[:120])
        key, data = match.groups()
        data = ast.literal_eval(data.replace('{', '[').replace('}', ']'))
        require(key not in rows, 'Repeated Wolfram structure key ' + key)
        records = {}
        for monomial, multiplicity in data:
            mono = tuple(monomial)
            require(len(mono) == 6 and all(isinstance(n, int) and n >= 0 for n in mono),
                    'Malformed manifest monomial')
            require(mono not in records and sum(mono) == 4 and multiplicity > 0,
                    'Duplicate or invalid manifest entry')
            records[mono] = multiplicity
        rows[key] = records
    require(len(rows) == 13, f'Expected 13 distinct representation/hypercharge structures, got {len(rows)}')
    return rows


def model_manifest_key(model: str) -> str:
    family, alpha = name_alpha(model)
    return f'{DIMS[family][0]}_{DIMS[family][1]}_{alpha}'


def audit_operator_basis() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    manifest = parse_wolfram_manifest(CATALOG)
    rows = []
    for model in MODELS:
        spectrum = scalar_operators(model)
        extra = {k: v for k, v in spectrum[4].items() if k not in BASELINE_KEYS}
        original = manifest[model_manifest_key(model)]
        q_diff = {
            'missing_from_matchete_manifest': [list(k) for k in sorted(extra.keys() - original.keys())],
            'unexpected_in_matchete_manifest': [list(k) for k in sorted(original.keys() - extra.keys())],
            'multiplicity_mismatch': [list(k) for k in extra.keys() & original.keys()
                                      if extra[k] != original[k]],
        }
        # Baseline normalizations are NOT checked by this structural comparison.
        quadratics = spectrum[2]
        extra_bilinears = sum(quadratics.values()) - 2  # two normal masses
        ok = all(not v for v in q_diff.values()) and not spectrum[3] and extra_bilinears in (0, 1)
        rows.append({'model': model, 'passed': ok,
                     'quartic_slots': sum(extra.values()),
                     'independent_extra_bilinears': extra_bilinears,
                     'allowed_scalar_cubic_singlets': sum(spectrum[3].values()),
                     'extra_quartic_structures': [
                         {'species_counts': list(k), 'multiplicity': v}
                         for k, v in sorted(extra.items())],
                     'manifest_comparison': q_diff})
    totals = {'models': len(rows), 'quartic_slots': sum(x['quartic_slots'] for x in rows),
              'bilinear_slots': sum(x['independent_extra_bilinears'] for x in rows),
              'scalar_cubic_slots': sum(x['allowed_scalar_cubic_singlets'] for x in rows),
              'all_models_match_manifest': all(x['passed'] for x in rows)}
    return rows, totals


def complex_rand(rng: np.random.Generator, shape: tuple[int, ...]) -> np.ndarray:
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def contract(tensor: np.ndarray, vectors: list[np.ndarray]) -> complex:
    out = np.asarray(tensor, dtype=complex)
    for v in vectors:
        out = np.tensordot(out, v, axes=(0, 0))
    return complex(out)


def generic_self_transport(model: str, mu: complex, *, trials: int = 3) -> dict[str, Any]:
    """Every relevant scalar source structure + representative Yukawa tensors.

    Complex neutral field is native Cartesian: S=e^(i arg(mu)/2)(R+iI)/sqrt2.
    This identity transports arbitrary invariant tensors without claiming
    their native Clebsch-Gordan components have been matched to RGBeta.
    """
    f, a = name_alpha(model)
    d1, d2, dferm = DIMS[f]
    neutral = 'S2' if a == -2 else 'S1'
    dim = d2 if neutral == 'S2' else d1
    require(dim in (1, 3), 'Self neutral field must be a real SU2 irrep')
    phase = np.exp(0.5j * np.angle(mu))
    scale = 1 / np.sqrt(2)
    src = scalar_operators(model)[4]
    source_specs: list[tuple[str, ...]] = []
    # Include baseline scalar quartics (independent from extra structures).
    for counts in src:
        names = tuple(s for s, n in zip(SPECIES, counts) for _ in range(n))
        if neutral in names or (neutral + 'dag') in names:
            source_specs.extend([names] * src[counts])
    # Scalar in every T3 Yukawa (fermion tensors are arbitrary spectators).
    source_specs.extend([
        ('Ldag', 'F', neutral), ('Ldag', 'F', neutral + 'dag')
    ])
    rng = np.random.default_rng(33021 + len(model) + dim)
    max_error = 0.0
    checked = 0
    for names in source_specs:
        dims = {'H': 2, 'Hdag': 2, 'S1': d1, 'S1dag': d1,
                'S2': d2, 'S2dag': d2, 'Ldag': 2, 'F': dferm}
        size = tuple(dims[x] for x in names)
        for _ in range(trials):
            r = rng.normal(size=dim)
            im = rng.normal(size=dim)
            s = phase * (r + 1j * im) * scale
            other = complex_rand(rng, (d1 if neutral == 'S2' else d2,))
            vectors: dict[str, np.ndarray] = {'H': complex_rand(rng, (2,)),
                'Ldag': complex_rand(rng, (2,)), 'F': complex_rand(rng, (dferm,))}
            vectors['Hdag'] = vectors['H'].conj()
            for label in ('S1', 'S2'):
                vec = s if label == neutral else other
                vectors[label] = vec
                vectors[label + 'dag'] = vec.conj()
            vectors[neutral + '_R'] = r
            vectors[neutral + '_I'] = im
            tensor = complex_rand(rng, size)
            before = contract(tensor, [vectors[k] for k in names])
            opts = []
            for name in names:
                if name in (neutral, neutral + 'dag'):
                    barred = name.endswith('dag')
                    coeff = phase.conjugate() if barred else phase
                    opts.append(((neutral+'_R', coeff*scale),
                                 (neutral+'_I', coeff*(-1j if barred else 1j)*scale)))
                else:
                    opts.append(((name, 1.0),))
            after = 0j
            for choice in product(*opts):
                coeff = np.prod([x[1] for x in choice])
                after += coeff * contract(tensor, [vectors[x[0]] for x in choice])
            error = abs(before - after) / max(1., abs(before), abs(after))
            max_error = max(max_error, float(error))
            checked += 1
    require(max_error < 1e-11, f'{model}: self tensor transport error {max_error}')
    return {'tested_source_slots_including_baseline_and_yukawa': len(source_specs),
            'random_complex_tensor_trials': checked,
            'maximum_relative_residual': max_error,
            'native_cartesian_basis_not_cross_package_certified': True}


def read_json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), f'JSON object required: {path}')
    return obj


def physical_self_benchmark(model: str, report: dict[str, Any]) -> dict[str, Any]:
    from Lagrangian.RunSixSelfBilinearEigenAudit import (
        evaluate_report, native_cartesian_benchmark, parse_native_factor, native_dim)
    token = model.replace('_m2', '_-2').replace('_p0', '_+0')
    original = evaluate_report(report, token)
    require(original['passed'], 'Prior self eigenbasis report does not certify passing')
    factor = parse_native_factor(report['native_bilinear_metric_factor_inputform'])
    lam = (120.0**2) * np.exp(0.37j)
    mu = -lam * factor  # potential = - Lagrangian for constant scalars
    transported = generic_self_transport(model, mu)
    masses = [float(900.0**2 + 2*abs(mu)), float(900.0**2 - 2*abs(mu))]
    numerical = native_cartesian_benchmark(native_dim(token), factor)
    return {'passed': True, 'native_factor': str(factor),
            'eigenmass_squared_gev2': masses,
            'physical_relation_verified_for_numerical_benchmark': True,
            'potential_diagonalization': numerical,
            'interaction_tensor_transport': transported,
            'symbolic_matchete_mass_spurions_identified': False}


def generic_bc_benchmark(model: str) -> dict[str, Any]:
    from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
        MixedDoubletMassBasis, verify_tensor_transport, quadratic_mass_tensor_transport)
    token = model.replace('_m1', '_-1')
    basis = MixedDoubletMassBasis.build(token, 930.0**2, 1110.0**2,
                                       (82.0**2) * np.exp(0.41j))
    quartics = scalar_operators(model)[4]
    tests = []
    for counts, n in sorted(quartics.items()):
        fields = tuple(s for s, k in zip(SPECIES, counts) for _ in range(k))
        for _ in range(n):
            tests.append({'fields': list(fields),
                          'max_relative_residual': verify_tensor_transport(fields, basis.U,
                                                                           seed=125+len(tests), trials=3)})
    mass = quadratic_mass_tensor_transport(basis)
    rng = np.random.default_rng(1003)
    errors = []
    for _ in range(20):
        x, y = complex_rand(rng, (2,)), complex_rand(rng, (2,))
        s1, s2 = basis.original_fields(x, y)
        recovered = basis.eigen_fields(s1, s2)
        require(np.allclose(recovered, [x, y], rtol=1e-11, atol=1e-9),
                'B/C inverse transformation failed')
        a, b = basis.old_potential(s1, s2), basis.new_potential(x, y)
        errors.append(abs(a-b) / max(1, abs(a), abs(b)))
    maximum = float(max(errors))
    require(maximum < 1e-11, 'B/C physical mass potential discrepancy')
    return {'passed': True, 'mass_tensor_transport': mass,
            'random_potential_max_relative_residual': maximum,
            'interaction_tensor_slots_tested': len(tests),
            'interaction_tensor_max_relative_residual': max(x['max_relative_residual'] for x in tests),
            'source_matchete_tensor_components_used': False}


def bc_saved_mass_relation(model: str) -> dict[str, Any]:
    """Check the *actual* B/C saved U, source A and eigenmasses, if available.

    The existing Matchete B/C driver reads exactly this prepared file.
    This verifies its numerical mass map without changing the source or
    claiming that the symbolic Matchete spurions are constrained.
    """
    p = ROOT / 'output' / 'general_t3_mass_eigen_matchete' / (model + '_final_parameters.json')
    require(p.is_file(), 'Missing original B/C prepared parameters: ' + str(p))
    obj = read_json(p)
    key = model.replace('_m1', '_-1')
    require(obj.get('format') == 't3-experimental-mass-eigen-matchete-parameters-v1' and
            obj.get('model') == key and
            obj.get('status') == 'NumericMassEigenstateBoundaryPreparedMatcheteNotRun',
            'Wrong B/C saved parameter schema')
    def complex_matrix(name: str) -> np.ndarray:
        raw = obj[name]
        m = np.asarray(raw['real'], dtype=float) + 1j * np.asarray(raw['imag'], dtype=float)
        require(m.shape == (2, 2) and np.isfinite(m).all(), 'Malformed matrix '+name)
        return m
    U = complex_matrix('unitary_U')
    A = complex_matrix('original_mass_squared_A_gev2')
    eig = np.asarray(obj['heavy_mass_squared_gev2'], dtype=float)
    require(eig.shape == (2,) and np.isfinite(eig).all() and (eig > 0).all(),
            'Invalid B/C eigenmasses')
    scale = max(1., float(np.max(abs(A))))
    unitary = float(np.max(abs(U.conj().T @ U - np.eye(2))))
    diagonal = float(np.max(abs(U.conj().T @ A @ U - np.diag(eig))))
    hermitian = float(np.max(abs(A-A.conj().T)))
    require(unitary < 2e-10 and diagonal < 2e-9*scale and hermitian < 2e-9*scale,
            'Saved B/C physical eigenmass map not diagonal')
    return {'passed': True, 'parameter_source': str(p),
            'source_sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
            'unitarity_max_error': unitary,
            'diagonal_mass_max_error_gev2': diagonal,
            'hermiticity_max_error_gev2': hermitian,
            'eigenmass_squared_gev2': eig.tolist(),
            'symbolic_matchete_spurions_identified': False}


def prior_eigen_evidence(model: str) -> dict[str, Any]:
    if model in SELF:
        path = SELFDIR / (model.split('_')[0] + '_' + model.split('_')[1] + '.json')
        if not path.is_file():
            return {'status': 'MissingPriorEigenReport', 'passed': False, 'path': str(path)}
        try:
            r = read_json(path)
            return {'status': r.get('status'), 'passed': r.get('status') == SELF_SUCCESS,
                    'path': str(path), 'native_numerical': physical_self_benchmark(model, r)}
        except (ValueError, ImportError, KeyError, TypeError, RuntimeError) as e:
            return {'status': 'PriorEigenReportVerificationFailed', 'passed': False,
                    'path': str(path), 'diagnostic': str(e)}
    if model in BC:
        path = BCMF_DIR / (model + '_mf_eigenbasis_audit.json')
        if not path.is_file():
            return {'status': 'MissingPriorBCMFReport', 'passed': False, 'path': str(path)}
        try:
            r = read_json(path)
            require(r.get('status') == BCMF_SUCCESS, 'B/C MF and heavy eigenbasis not both passed')
            return {'status': r['status'], 'passed': True,
                    'path': str(path), 'native_numerical': generic_bc_benchmark(model),
                    'saved_physical_mass_map': bc_saved_mass_relation(model)}
        except (ValueError, ImportError, KeyError, TypeError, RuntimeError) as e:
            return {'status': 'PriorBCReportVerificationFailed', 'passed': False,
                    'path': str(path), 'diagnostic': str(e)}
    return {'status': 'NotExceptional', 'passed': True}


def run_z2(model: str, out_dir: Path, command: str, timeout: int) -> dict[str, Any]:
    output = out_dir / 'native_z2' / (model + '.json')
    log = out_dir / 'native_z2' / (model + '.log')
    output.parent.mkdir(parents=True, exist_ok=True)
    family, tag = model.split('_')
    cmd = [command, '-file', str(WL), str(output), family, tag]
    env = dict(os.environ)
    env.update({name: '1' for name in FLAGS})
    env['T3_EXPERIMENTAL_GENERAL_QUADRATICS'] = '0'
    print(f'Z2 Matchete {model}', flush=True)
    # Never accept a stale JSON from a previous subprocess run.
    output.unlink(missing_ok=True)
    try:
        proc = subprocess.run(cmd, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, timeout=timeout, check=False)
        logs = proc.stdout
        log.write_text(logs.decode('utf-8', 'replace'), encoding='utf-8')
        data = read_json(output)
        valid = (proc.returncode == 0 and data.get('format') == 't3-workstream-a-term-parity-v1'
                 and data.get('model') == model and data.get('status') == Z2_SUCCESS)
        return {'status': data.get('status', 'Missing'), 'passed': bool(valid),
                'exit_code': proc.returncode, 'data': data, 'log': str(log), 'path': str(output)}
    except (OSError, subprocess.TimeoutExpired, ValueError, KeyError) as exc:
        log.write_text(str(exc), encoding='utf-8')
        return {'status': 'NativeZ2RunFailed', 'passed': False,
                'diagnostic': str(exc), 'log': str(log), 'path': str(output)}


def collect_z2(model: str, out_dir: Path) -> dict[str, Any]:
    p = out_dir / 'native_z2' / (model + '.json')
    if not p.is_file():
        return {'status': 'NotRun', 'passed': False, 'path': str(p)}
    try:
        obj = read_json(p)
        ok = obj.get('format') == 't3-workstream-a-term-parity-v1' and obj.get('model') == model
        return {'status': obj.get('status'), 'passed': bool(ok and obj.get('status') == Z2_SUCCESS),
                'data': obj, 'path': str(p)}
    except (ValueError, OSError, KeyError) as exc:
        return {'status': 'InvalidNativeZ2Report', 'passed': False,
                'diagnostic': str(exc), 'path': str(p)}


def write_results(out_dir: Path, basis_rows: list[dict[str, Any]],
                  totals: dict[str, Any], z2_rows: dict[str, dict[str, Any]]) -> dict[str, Any]:
    rows = []
    for basis in basis_rows:
        model = basis['model']
        eigen = prior_eigen_evidence(model)
        z2 = dict(z2_rows[model])
        if z2['passed']:
            data = z2['data']
            expected_y = 1 if model in SELF else 0
            count_ok = (data.get('extra_quartics_included') == basis['quartic_slots'] and
                        data.get('extra_yukawas_included') == expected_y and
                        data.get('extra_quadratics_included') == 0 and
                        data.get('matchete_full_check_passed') is True and
                        data.get('indexed_mf_in_original_free_lagrangian') is True and
                        data.get('discrete_z2_explicitly_registered_in_matchete') is False)
            if not count_ok:
                z2['passed'] = False
                z2['status'] = 'NativeParityInventoryOrScopeMismatch'
        rows.append({'model': model, 'independent_scalar_basis': basis,
                     'z2_parity_of_actual_unrotated_bsm_terms': z2,
                     'previous_heavy_eigenbasis': eigen,
                     'complete_workstream_a_certified': False})
    result = {
        'format': 't3-workstream-a-uv-closure-stage-v1',
        'status': 'StructuralAndNumericalStagePassed' if
            totals['all_models_match_manifest'] and
            all(x['z2_parity_of_actual_unrotated_bsm_terms']['passed'] for x in rows) and
            all(x['previous_heavy_eigenbasis']['passed'] for x in rows) else
            'PartialOrFailedEvidence',
        'source': 'IndependentSU2SymmetricPowerCounter+ExistingMatcheteZ2AndEigenReports',
        'catalog_file': str(CATALOG),
        'catalog_sha256': hashlib.sha256(CATALOG.read_bytes()).hexdigest(),
        'z2_driver_sha256': hashlib.sha256(WL.read_bytes()).hexdigest(),
        'scalar_basis_totals': totals,
        'native_z2_pass_count': sum(x['z2_parity_of_actual_unrotated_bsm_terms']['passed'] for x in rows),
        'prior_eigenbasis_pass_count': sum(x['previous_heavy_eigenbasis']['passed'] for x in rows if x['model'] in SELF+BC),
        'remaining_requirements': [
            'Symbolic original-to-eigenmass and phase parameter relation in native Matchete',
            'Exact independent full interaction CG normalization and Matchete/RGBeta map',
            'Explicitly registered Z2 selection rule, not just parity of assembled terms',
            'Independent exhaustive Lorentz/fermion operator completeness and flavour tensors',
            'Physical EFT matching / threshold normalization (outside Workstream A)'
        ],
        'workstream_a_full_uv_completeness_certified': False,
        'rows': rows,
    }
    (out_dir / 'uv_closure_stage.json').write_text(json.dumps(result, indent=2,
        allow_nan=False, default=lambda o: float(o) if isinstance(o, np.floating) else str(o))+'\n',encoding='utf-8')
    lines = ['# Workstream A — UV closure stage', '',
             f"Outcome: **{result['status']}**. This is not a full UV certificate.", '',
             '| Model | Extra quartics | Extra bilinear | Scalar structure | Actual Z2 terms | Eigenbasis |',
             '|---|---:|---:|---|---|---|']
    for row in rows:
        s = row['independent_scalar_basis']; z = row['z2_parity_of_actual_unrotated_bsm_terms']; e = row['previous_heavy_eigenbasis']
        lines.append(f"| {row['model']} | {s['quartic_slots']} | {s['independent_extra_bilinears']} | "
                     f"{'PASS' if s['passed'] else 'FAIL'} | {'PASS' if z['passed'] else z['status']} | "
                     f"{'PASS' if e['passed'] else e['status']} |")
    lines += ['', '## Scope and boundaries',
              '', 'Independent scalar counting uses Bose-symmetrized SU(2) Sym^n weight multiplicities, U(1)Y, and the assumed exact Z2.',
              'Actual-terms parity uses Matchete-built BSM interactions and free terms with extra bilinears disabled, so B/C can be built in their original basis.',
              'Eight native mass eigenbases are certified only by the earlier user-generated Matchete reports; this stage adds numerical potential and arbitrary-tensor transport checks.',
              'The six self-scalar Matchete eigenmass spurions remain independent symbolic masses; the numerical relation is conditional on the source potential sign and native factor.',
              '', '## Not certified', '']
    lines.extend('- '+x for x in result['remaining_requirements'])
    lines.append('')
    (out_dir / 'uv_closure_stage.md').write_text('\n'.join(lines),encoding='utf-8')
    return result


def self_test() -> None:
    require(len(MODELS) == 16, 'Wrong class enumeration')
    row, totals = audit_operator_basis()
    require(totals == {'models': 16, 'quartic_slots': 50, 'bilinear_slots': 8,
                       'scalar_cubic_slots': 0, 'all_models_match_manifest': True},
            'Independent scalar basis disagrees with existing 16-point manifest: '+str(totals))
    # Verify all 8 mass-basis calculations without needing the saved Matchete reports.
    for name in SELF:
        trial = generic_self_transport(name, (80**2) * np.exp(.29j), trials=2)
        require(trial['maximum_relative_residual'] < 1e-11, 'Self operator transport failed')
    try:
        import Numerical.experimental.GeneralT3MixedScalarMassBasis  # noqa: F401
    except ImportError:
        print('NOTE: B/C repository source unavailable in this isolated test runtime')
    else:
        for name in BC:
            b = generic_bc_benchmark(name)
            require(b['passed'], 'Mixed scalar tensor transport failed')
    # Structural proof: an SU2 doublet's antisymmetric self tensor vanishes
    # under exact Bose symmetrization, while 2⊗2 contains one singlet.
    require(singlet_multiplicity((2,), (2,)) == 0, 'Bose projection missing')
    require(singlet_multiplicity((2, 2), (1, 1)) == 1, 'SU2 singlet count incorrect')
    require('Needs["Matchete`"]' in WL.read_text(encoding='utf-8'), 'Native parity driver missing')
    print('PASS: 16-model independent scalar operator basis (50 quartic, 8 bilinear, 0 cubic)')
    print('PASS: six self-scalar arbitrary interaction-tensor transport benchmarks')
    print('NOTE: B/C numerical tests run when the repository Numerical package is available')
    print('NOTE: neither Wolfram kernel nor user-generated eigenbasis reports inspected')


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    modes = p.add_mutually_exclusive_group(required=True)
    modes.add_argument('--self-test', action='store_true')
    modes.add_argument('--static-only', action='store_true',
                       help='Independent enumeration + aggregate existing reports; no Wolfram runs')
    modes.add_argument('--all', action='store_true',
                       help='Execute native parity audits for all 16 models, aggregate evidence')
    p.add_argument('--output-dir', type=Path,
                   default=ROOT / 'Reports' / 'output' / 'workstream_a' / 'uv_closure_stage')
    p.add_argument('--wolframscript', default='wolframscript')
    p.add_argument('--timeout-seconds', type=int, default=2400)
    p.add_argument('--overwrite', action='store_true', help='Allow overwriting own stage reports')
    args = p.parse_args()
    if args.self_test:
        self_test()
        return 0
    require(args.timeout_seconds > 0, 'Timeout must be positive')
    basis, totals = audit_operator_basis()
    require(totals['all_models_match_manifest'], 'Scalar basis differs from current source')
    out_dir = args.output_dir.resolve()
    protected = [out_dir/'uv_closure_stage.json', out_dir/'uv_closure_stage.md']
    if args.all:
        protected += [out_dir/'native_z2'/(model + ext) for model in MODELS for ext in ('.json', '.log')]
    conflicts = [str(p) for p in protected if p.exists()]
    require(not conflicts or args.overwrite, 'Results exist; back up first then pass --overwrite: '+str(conflicts[:4]))
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.all:
        for model in MODELS:
            result = run_z2(model, out_dir, args.wolframscript, args.timeout_seconds)
            print(' ', result['status'], 'PASS' if result['passed'] else 'FAIL', flush=True)
            if result['status'] == 'NativeZ2RunFailed' and ('No such file' in result.get('diagnostic', '') or
                                                          'not found' in result.get('diagnostic', '')):
                print('STOP: Wolfram executable unavailable; remaining native checks skipped')
                break
    z2 = {model: collect_z2(model, out_dir) for model in MODELS}
    result = write_results(out_dir, basis, totals, z2)
    print(f"RESULT: {result['status']} | scalar basis 16/16 | "
          f"native Z2 {result['native_z2_pass_count']}/16 | "
          f"prior heavy eigenbasis {result['prior_eigenbasis_pass_count']}/8")
    print('JSON:', out_dir/'uv_closure_stage.json')
    print('REPORT:', out_dir/'uv_closure_stage.md')
    return 0 if result['status'] == 'StructuralAndNumericalStagePassed' else 6


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError, ImportError) as exc:
        print('ERROR:', exc, file=sys.stderr)
        sys.exit(2)
