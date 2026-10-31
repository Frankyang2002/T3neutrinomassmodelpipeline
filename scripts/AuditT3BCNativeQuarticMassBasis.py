"""Audit actual Matchete extra-quartic tensors in B/C(-1) scalar mass basis.

Read-only audit of the saved Stage22 *native* Matchete SU(2) tensors and the
already verified candidate mixed-doublet mass-basis report. Rotates each
independent tensor with the full doublet epsilon map, checks numerical
polynomials at complex field samples, and retains both independent contractions
in the 110101 block. No coupling equality with RGBeta, no native Matchete
field registry update, no one-loop Match[], and no sequential EFT matching.

    python -u -X utf8 -m scripts.AuditT3BCNativeQuarticMassBasis --all --run-id bc_native_vertices_v1
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import re

import numpy as np
import sympy as sp

from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
    EPS, MixedDoubletMassBasis, field_map,
)

ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT / 'output/general_t3_stage22/matchete_tensors.json'
MASS_BASIS = ROOT / 'output/general_t3_mass_basis'
OUTPUT = ROOT / 'output/audits/t3_bc_native_quartic_mass_basis'
MODELS = ('B_-1', 'C_-1')
EXPECTED = {
    'lambdaT3ExtraQ110101I1': ('H', 'Hdag', 'S1dag', 'S2dag'),
    'lambdaT3ExtraQ110101I2': ('H', 'Hdag', 'S1dag', 'S2dag'),
    'lambdaT3ExtraQ000112I1': ('S1dag', 'S2', 'S2dag', 'S2dag'),
    'lambdaT3ExtraQ000202I1': ('S1dag', 'S1dag', 'S2dag', 'S2dag'),
    'lambdaT3ExtraQ001201I1': ('S1', 'S1dag', 'S1dag', 'S2dag'),
    'lambdaT3ExtraQ020020I1': ('Hdag', 'Hdag', 'S2', 'S2'),
    'lambdaT3ExtraQ020200I1': ('Hdag', 'Hdag', 'S1dag', 'S1dag'),
}


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('T3 B/C native quartic mass-basis audit: ' + message)


def read_json(path: Path) -> dict:
    require(path.is_file(), 'missing source: ' + str(path))
    value = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(value, dict), 'expected JSON object: ' + str(path))
    return value


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def complex_matrix(obj: dict, shape: tuple[int, ...]) -> np.ndarray:
    require(isinstance(obj, dict) and set(obj) == {'real', 'imag'}, 'complex array schema')
    value = np.asarray(obj['real'], dtype=float) + 1j*np.asarray(obj['imag'], dtype=float)
    require(value.shape == shape and np.isfinite(value).all(), 'bad complex matrix')
    return value


def array_payload(a: np.ndarray) -> dict:
    value = np.asarray(a, dtype=complex)
    require(np.isfinite(value).all(), 'nonfinite component tensor')
    return {'shape': list(value.shape), 'real': value.real.ravel().tolist(),
            'imag': value.imag.ravel().tolist()}


def native_number(value: str) -> complex:
    """Parse the restricted numeric InputForm emitted by saved native exporter."""
    require(isinstance(value, str) and
            re.fullmatch(r'[0-9+\-*/^()\[\]Sqrt\s]+', value) is not None,
            'unexpected native tensor component: ' + repr(value))
    expr = sp.sympify(value.replace('Sqrt[', 'sqrt(').replace(']', ')').replace('^', '**'),
                      locals={'sqrt': sp.sqrt})
    require(not expr.free_symbols and expr.is_finite is True,
            'nonnumeric native tensor component')
    z = complex(sp.N(expr, 35))
    require(np.isfinite([z.real, z.imag]).all(), 'nonfinite native component')
    return z


def native_tensors(source: dict, model: str) -> dict[str, tuple[tuple[str, ...], np.ndarray]]:
    require(source.get('format') == 'stage22-exact-tensors-v1' and
            source.get('source') == 'MatcheteInvariantTensors' and
            source.get('basis') == 'Matchete', 'not an exported Matchete tensor source')
    records = source.get('models', {}).get(model, {}).get('tensors')
    require(isinstance(records, list) and len(records) == 7,
            model + ': expected exactly seven independent quartics')
    result = {}
    for record in records:
        name = record.get('symbol')
        require(name in EXPECTED and name not in result,
                model + ': missing/duplicate/unknown native invariant ' + str(name))
        species = tuple(record.get('species', ()))
        require(species == EXPECTED[name] and record.get('dimensions') == [2]*4 and
                record.get('potential_prefactor') == '1',
                model + '/' + name + ': native invariant contract changed')
        counts = [sum(f == s for f in species)
                  for s in ('H', 'Hdag', 'S1', 'S1dag', 'S2', 'S2dag')]
        require(record.get('counts') == counts, model + '/' + name + ': species mismatch')
        components = record.get('components')
        require(isinstance(components, list) and len(components) == 16,
                model + '/' + name + ': missing 2^4 native entries')
        tensor = np.asarray([native_number(z) for z in components], dtype=complex).reshape((2,)*4)
        require(np.linalg.norm(tensor) > 0, model + '/' + name + ': zero native tensor')
        result[name] = species, tensor
    require(set(result) == set(EXPECTED), model + ': incomplete extra quartic inventory')
    return result


def load_basis(path: Path, model: str) -> MixedDoubletMassBasis:
    report = read_json(path)
    require(report.get('format') == 't3-mixed-doublet-mass-basis-transport-v1' and
            report.get('status') == 'ExperimentalMassBasisTensorTransportNoPhysicalMatching' and
            report.get('model') == model and report.get('endpoint') == 'final' and
            report.get('physical_mass_matching_performed') is False and
            report.get('MatcheteRGBetaPhysicalCouplingMapVerified') is False,
            model + ': unexpected source mass-basis certificate')
    data = report.get('basis', {})
    A = complex_matrix(data.get('mass_squared_flavour_matrix_gev2'), (2, 2))
    U = complex_matrix(data.get('unitary_flavour_rotation_U'), (2, 2))
    require(np.allclose(A, A.conj().T, atol=1e-8, rtol=1e-12), 'mass matrix not Hermitian')
    basis = MixedDoubletMassBasis.build(model, float(A[0, 0].real),
                                       float(A[1, 1].real), complex(A[0, 1]))
    require(not basis.tachyonic and not basis.degenerate and
            np.allclose(basis.A, A, atol=1e-8, rtol=1e-12) and
            np.allclose(U.conj().T @ A @ U, np.diag(basis.eigenvalues_sq),
                        rtol=1e-11, atol=1e-6) and
            np.allclose(U, basis.U, rtol=1e-11, atol=1e-11),
            model + ': saved scalar rotation differs from mass matrix')
    saved = np.asarray(data.get('masses_squared_gev2', []), dtype=float)
    require(saved.shape == (2,) and np.allclose(saved, basis.eigenvalues_sq,
                                                rtol=1e-11, atol=1e-6),
            model + ': saved eigenvalues mismatch')
    return basis


def rotate_tensor(tensor: np.ndarray, species: tuple[str, ...], U: np.ndarray
                 ) -> dict[tuple[str, ...], np.ndarray]:
    """Compute T'_{jklm}=T_{abcd} M0_{aj} M1_{bk} M2_{cl} M3_{dm}."""
    require(tensor.shape == (2, 2, 2, 2) and len(species) == 4,
            'expected a fourth-order doublet tensor')
    rotated = {}
    for choice in product(*(field_map(field, U) for field in species)):
        labels = tuple(name for name, _ in choice)
        matrices = [matrix for _, matrix in choice]
        term = np.einsum('abcd,ai,bj,ck,dl->ijkl', tensor, *matrices, optimize=True)
        rotated[labels] = rotated.get(labels, np.zeros((2,)*4, dtype=complex)) + term
    return rotated


def verify_polynomial(species: tuple[str, ...], tensor: np.ndarray,
                      rotated: dict[tuple[str, ...], np.ndarray],
                      basis: MixedDoubletMassBasis, *, trials: int = 12,
                      seed: int = 2081) -> float:
    """Compare native old-basis polynomial with explicit rotated components."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(trials):
        rand = lambda: rng.normal(size=2) + 1j*rng.normal(size=2)
        h, x, y = rand(), rand(), rand()
        s1, s2 = basis.original_fields(x, y)
        old = {'H':h, 'Hdag':h.conjugate(), 'S1':s1, 'S1dag':s1.conjugate(),
               'S2':s2, 'S2dag':s2.conjugate()}
        new = {'H':h, 'Hdag':h.conjugate(), 'X':x, 'Xdag':x.conjugate(),
               'Y':y, 'Ydag':y.conjugate()}
        before = np.einsum('abcd,a,b,c,d->', tensor, *(old[s] for s in species))
        after = sum(np.einsum('abcd,a,b,c,d->', t, *(new[s] for s in labels))
                    for labels, t in rotated.items())
        residual = abs(before-after)/max(1.0, abs(before), abs(after))
        worst = max(worst, float(residual))
        # The native potential has lambda*T + h.c.; test a genuinely complex coupling.
        lam = 0.23 + 0.37j
        real_residual = (abs(2*np.real(lam*before)-2*np.real(lam*after)) /
                         max(1., abs(2*np.real(lam*before)), abs(2*np.real(lam*after))))
        worst = max(worst, float(real_residual))
    require(worst < 3e-12, 'native quartic polynomial transport failed')
    return worst


def shared_independent_block(a: dict, b: dict) -> dict:
    """Verify 110101 I1/I2 remain independent at an identical new field signature."""
    candidates = set(a) & set(b)
    require(bool(candidates), 'no shared field labels for independent contractions')
    key = max(sorted(candidates),
              key=lambda k: min(np.linalg.norm(a[k]), np.linalg.norm(b[k])))
    vectors = np.stack([a[key].ravel(), b[key].ravel()])
    singular = np.linalg.svd(vectors, compute_uv=False)
    require(singular[0] > 1e-12 and singular[-1] / singular[0] > 1e-9,
            '110101 two independent native SU2 contractions collapsed')
    return {'field_signature':list(key), 'rank':2,
            'small_to_large_singular_value_ratio':float(singular[-1]/singular[0])}


def run_model(model: str, *, native: Path = NATIVE, mass_basis: Path = MASS_BASIS,
              output: Path = OUTPUT, trials: int = 12) -> dict:
    require(model in MODELS, 'only B_-1 and C_-1')
    require(1 <= trials <= 100, 'invalid trial count')
    mass_path = mass_basis / (model.replace('-', 'm') + '_final_mass_basis.json')
    basis = load_basis(mass_path, model)
    sources = native_tensors(read_json(native), model)
    reports = {}
    transformed = {}
    worst = 0.0
    for i, (name, (species, source)) in enumerate(sources.items()):
        rotated = rotate_tensor(source, species, basis.U)
        residual = verify_polynomial(species, source, rotated, basis, trials=trials,
                                     seed=20900+i)
        worst = max(worst, residual)
        transformed[name] = rotated
        reports[name] = {
            'native_source_species':list(species),
            'native_source_component_norm':float(np.linalg.norm(source)),
            'rotated_terms':[
                {'species':list(labels), 'component_tensor':array_payload(t)}
                for labels, t in sorted(rotated.items())],
            'polynomial_and_hc_max_relative_residual':residual,
            'native_potential_prefactor': '1',
            'rgbeta_coupling_normalization_certified':False,
        }
    independent = shared_independent_block(transformed['lambdaT3ExtraQ110101I1'],
                                            transformed['lambdaT3ExtraQ110101I2'])
    report = {
        'format':'t3-bc-native-matchete-quartic-mass-basis-audit-v1',
        'status':'NativeMatcheteExtraQuarticTensorTransportChecked',
        'model':model, 'source_stage22_matchete_tensors':native.relative_to(ROOT).as_posix()
                 if native.is_relative_to(ROOT) else str(native),
        'source_mass_basis_report':mass_path.relative_to(ROOT).as_posix()
                 if mass_path.is_relative_to(ROOT) else str(mass_path),
        'sha256':{'stage22_matchete_tensors':digest(native),
                  'mixed_mass_basis_report':digest(mass_path)},
        'native_extra_quartic_count':7,
        'max_native_polynomial_relative_error':worst,
        'two_contraction_110101':independent,
        'quartic_tensors':reports,
        'all_source_vertex_contractions_rotated':False,
        'all_native_extra_quartics_rotated':True,
        'matchete_rgbeta_coupling_map_certified':False,
        'matchete_redefined_mass_basis_lagrangian_validated':False,
        'mass_eigenstate_matchete_match_performed':False,
        'sequential_eft_wilson_matching_performed':False,
        'scope':('Actual saved native Matchete EXTRA quartic tensors only; '
                 'the baseline scalar quartics, Yukawa CGs, and RGBeta-native '
                 'coupling normalization have not been certified by this audit.'),
    }
    output.mkdir(parents=True, exist_ok=True)
    target = output / (model.replace('-', 'm') + '_native_quartic_mass_basis.json')
    require(not target.exists(), 'output exists: ' + str(target))
    target.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    selector = parser.add_mutually_exclusive_group(required=True)
    selector.add_argument('--all', action='store_true')
    selector.add_argument('--model', choices=MODELS)
    parser.add_argument('--run-id', required=True)
    parser.add_argument('--trials', type=int, default=12)
    ns = parser.parse_args()
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', ns.run_id) is not None,
            'invalid run ID')
    destination = OUTPUT / ns.run_id
    require(not destination.exists(), 'run ID already used: ' + str(destination))
    chosen = MODELS if ns.all else (ns.model,)
    failures = []
    rows = []
    for model in chosen:
        try:
            row = run_model(model, output=destination, trials=ns.trials)
            rows.append({'model':model, 'status':row['status'],
                         'max_native_polynomial_relative_error':row['max_native_polynomial_relative_error'],
                         'independent_contraction_rank':row['two_contraction_110101']['rank']})
            print(f'{model}: native tensors 7/7; 110101 rank=2; '
                  f'polynomial residual={row["max_native_polynomial_relative_error"]:.3e}',
                  flush=True)
        except (ValueError, TypeError, OSError, RuntimeError, KeyError,
                np.linalg.LinAlgError) as exc:
            failures.append({'model':model, 'error':str(exc)})
            print(f'{model}: BLOCKED {exc}',flush=True)
    destination.mkdir(parents=True, exist_ok=True)
    summary = {'format':'t3-bc-native-quartic-mass-basis-batch-v1',
               'status':'SourceNativeExtraQuarticTransportChecked' if not failures else 'Incomplete',
               'rows':rows,'errors':failures,
               'physical_matching_certified':False}
    (destination/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',
                                            encoding='utf-8')
    print('Summary:', destination/'summary.json',flush=True)
    return 2 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
