"""Exact source-level Y/LLHH pole, logarithm and neutral-flavour ledger.

This module does NOT calculate an MS-bar counterterm, renormalize Matchete's
one-loop output, or certify a physical Weinberg coefficient.  The negative of
an observed pole is only the REQUIRED NET cancellation residue, not a counterterm
computed from the Lagrangian.  Do not use the formal finite part as physical C5.

Inputs are the archived native B(-1)/C(-1) Y-only Matchete expressions.
The algebraic conventions and parsers are reused from the inspected project
modules, rather than defining another set of gamma/epsilon conventions.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import zipfile

from Numerical.experimental.T3BCMatcheteFiniteTerms import (
    split_inputform, additive_terms, MASS_LOGS,
)
from Numerical.experimental.T3BCPoleLogConsistency import (
    compare_residue_and_log, terms,
)
from Numerical.experimental.T3BCNativeYLLHHFlavorKernel import (
    _record_from_term, reduce_finite_y_llhh, check_exact_chirality_conjugacy,
)

BASE = 'output/general_t3_mass_eigen_matchete/eft1_scalar_Y_only_canonical_probe'
MODELS = {'B_-1': 'B_m1', 'C_-1': 'C_m1'}


class LedgerError(ValueError):
    pass


def _hash(data: bytes) -> str:
    return sha256(data).hexdigest()


def _archive_text(archive: zipfile.ZipFile, manifest: dict, path: str) -> str:
    meta = manifest.get('files', {}).get(path)
    if not isinstance(meta, dict):
        raise LedgerError('Unmanifested native input: ' + path)
    blob = archive.read(path)
    if _hash(blob) != meta.get('sha256') or len(blob) != meta.get('size_bytes'):
        raise LedgerError('Native source hash/length mismatch: ' + path)
    return blob.decode('utf-8').strip()


def _canonical_key(d: dict) -> tuple:
    return (
        d['chirality'], d['quartic'], d['quartic_conjugated'],
        d['yukawa_i'], d['yukawa_j'], d['yukawas_conjugated'],
        d['mixing'], d['sqrt3_exponent'],
    )


def _project(source: str, model: str, count: int) -> dict:
    """Extend the existing strict 32-term neutral projector to 50 pole terms.

    Each term is still validated by the project's exact Matchete parser.
    No new field, CG, Majorana, or coupling conventions are introduced.
    """
    raw = terms(source)
    if len(raw) != count:
        raise LedgerError(f'{model}: expected {count} terms, got {len(raw)}')
    reduced: dict[tuple, Fraction] = defaultdict(Fraction)
    for term in raw:
        rec, coeff = _record_from_term(term, model)
        reduced[_canonical_key(rec)] += coeff
    entries = []
    for key, val in sorted(reduced.items()):
        if not val:
            continue
        ch, q, qbar, yi, yj, ybar, mix, sqrt3 = key
        entries.append({
            'chirality': ch,
            'quartic': q,
            'quartic_conjugated': qbar,
            'yukawa_i': yi,
            'yukawa_j': yj,
            'yukawas_conjugated': ybar,
            'mixing': [
                {'name': name, 'conjugated': bar, 'exponent': exp}
                for name, bar, exp in mix
            ],
            'sqrt3_exponent': sqrt3,
            'coefficient_numerator': val.numerator,
            'coefficient_denominator': val.denominator,
        })
    if not entries:
        raise LedgerError(f'{model}: component residue vanished unexpectedly')
    projected = {
        'model': model,
        'source_additive_terms': len(raw),
        'one_common_factor': 'hbar/Coupling[MF,{},0]',
        'neutral_SU2_slice': 'H_component_2__lepton_component_1__epsilon12_plus_1',
        'monomials': entries,
        'chirality_monomial_counts': {
            ch: sum(x['chirality'] == ch for x in entries) for ch in ('PL', 'PR')
        },
        'projection_type': 'ordered_neutral_component_not_physical_C5',
        'Majorana_Wilson_operator_normalization_certified': False,
        'MSbar_counterterms_certified': False,
        'physical_matching_certified': False,
    }
    if not check_exact_chirality_conjugacy(projected):
        raise LedgerError(f'{model}: pole chirality conjugacy failed')
    return projected


def _tree_llhh_term_count(source: str) -> int:
    """Identify an explicit tree-level Weinberg-like term in Matchete output."""
    count = 0
    for part in additive_terms(source):
        if ('GammaCC' in part and part.count('Field[l, Fermion') == 2 and
                part.count('Field[H, Scalar') == 2 and
                'Field[NewScalarX' not in part and
                'Field[NewScalarY' not in part):
            count += 1
    return count


def build(archive_path: Path, output_dir: Path) -> dict:
    if not archive_path.is_file():
        raise LedgerError('Missing native archive: ' + str(archive_path))
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        'format': 't3-bc-native-y-llhh-renormalization-ledger-v1',
        'interpretation': 'Native hard matching pole constraints only; NOT CT evaluation',
        'models': {},
        'physical_matching_certified': False,
        'physical_c5': None,
    }
    with zipfile.ZipFile(archive_path) as z:
        manifest = json.loads(z.read('MANIFEST.json'))
        if manifest.get('format') != 't3-bc-partial-xy-native-evidence-v3':
            raise LedgerError('Wrong native archive format')
        for model, prefix in MODELS.items():
            stem = f'{BASE}/{prefix}_final_y_only_canonical_match'
            raw_path = stem + '_raw_scalar_llhh.inputform.txt'
            tree_path = stem + '_sm_tree.inputform.txt'
            status_path = stem + '.json'
            raw = _archive_text(z, manifest, raw_path)
            tree = _archive_text(z, manifest, tree_path)
            status = json.loads(_archive_text(z, manifest, status_path))
            if (status.get('model') != model or
                status.get('tree_matching_completed') is not True or
                status.get('one_loop_matching_completed') is not True or
                status.get('hard_pole_counterterms_combined') is not False):
                raise LedgerError(model + ': stage/source status inconsistent')
            split = split_inputform(raw, stage='Y', sector='LLHH')
            if (len(split.terms), len(split.pole_terms),
                len(split.log_terms), len(split.constant_terms)) != (132, 50, 50, 32):
                raise LedgerError(model + ': changed Y/LLHH Laurent inventory')
            exact_check = compare_residue_and_log(split.pole_residue,
                                                   split.pole_free, 'Y')
            if not exact_check['pole_equals_log_coefficient_in_opaque_factor_algebra']:
                raise LedgerError(model + ': pole/log mismatch')
            pole = _project(split.pole_residue, model, 50)
            # Setting the native logarithm to 1 extracts its coefficient,
            # not a finite-part or a renormalization prescription.
            log_source = ''.join(t.replace(MASS_LOGS['Y'], '1')
                                 for t in split.log_terms)
            log_kernel = _project(log_source, model, 50)
            if pole['monomials'] != log_kernel['monomials']:
                raise LedgerError(model + ': pole/log neutral-component mismatch')
            finite = reduce_finite_y_llhh(split.scale_equal_mass, model)
            if not check_exact_chirality_conjugacy(finite):
                raise LedgerError(model + ': finite chirality conjugacy failed')
            tree_llhh_count = _tree_llhh_term_count(tree)
            if tree_llhh_count:
                raise LedgerError(model + ': unexpectedly found tree LLHH terms')
            # A formal *requirement*. There is no independently calculated
            # EFT1/EFT2 counterterm source in this archive.
            residue = split.pole_residue
            required_pole = '-((' + residue + ')/ϵ)'
            scale_slope = '2*(' + residue + ')'  # explicit log derivative only
            paths = {
                'raw_pole_residue': f'{prefix}_y_llhh_raw_pole_residue.inputform.txt',
                'required_net_pole_cancellation': f'{prefix}_y_llhh_required_net_pole_cancellation.NOT_A_COUNTERTERM.inputform.txt',
                'explicit_scale_derivative': f'{prefix}_y_llhh_explicit_d_dlnmubar_at_fixed_parameters.inputform.txt',
                'formal_finite_at_muY': f'{prefix}_y_llhh_formal_finite_at_muY.NOT_PHYSICAL_C5.inputform.txt',
                'pole_flavour_kernel': f'{prefix}_y_llhh_pole_neutral_flavour.json',
                'constant_flavour_kernel': f'{prefix}_y_llhh_constant_neutral_flavour.json',
            }
            output = {
                'raw_pole_residue': residue,
                'required_net_pole_cancellation': required_pole,
                'explicit_scale_derivative': scale_slope,
                'formal_finite_at_muY': split.scale_equal_mass,
            }
            for name, expr in output.items():
                (output_dir / paths[name]).write_text(expr + '\n', encoding='utf-8')
            for name, obj in [('pole_flavour_kernel', pole),
                              ('constant_flavour_kernel', finite)]:
                (output_dir / paths[name]).write_text(
                    json.dumps(obj, indent=2, ensure_ascii=False) + '\n',
                    encoding='utf-8',
                )
            summary['models'][model] = {
                'native_raw_sha256': manifest['files'][raw_path]['sha256'],
                'tree_y_llhh_explicit_operator_count': tree_llhh_count,
                'raw_terms': 132,
                'pole_terms': 50,
                'log_terms': 50,
                'constant_terms': 32,
                'opaque_pole_monomials': exact_check['nonzero_unique_opaque_pole_monomials'],
                'neutral_pole_monomials': len(pole['monomials']),
                'neutral_finite_monomials': len(finite['monomials']),
                'pole_log_opaque_and_neutral_identities_verified': True,
                'explicit_d_dlnmubar_at_fixed_parameters_equals_2_residue': True,
                'net_counterterm_residue_required_but_not_computed': True,
                'MSbar_renormalized_pole_cancellation_verified': False,
                'EFT1_EFT2_counterterms_supplied': False,
                'physical_C5': None,
                'outstanding_independent_checks': [
                    'EFT1_EFT2_operator_counterterm_difference',
                    'UV_vs_IR_pole_assignment',
                    'field_EOM_and_finite_scheme_consistency',
                    'Majorana_GammaCC_to_C5_normalization',
                    'EFT1_running_and_sequential_X_inheritance',
                ],
                'files': paths,
            }
    out = output_dir / 'summary.json'
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n',
                   encoding='utf-8')
    return summary


def self_test() -> None:
    # No fake amplitudes or counterterms: only local numeric/shape tests.
    ex = {
        'chirality': 'PR', 'quartic': 'lambdaT3',
        'quartic_conjugated': False,
        'yukawa_i': 'y1', 'yukawa_j': 'y2',
        'yukawas_conjugated': False,
        'mixing': (('T3MEU11', False, 1),),
        'sqrt3_exponent': 0,
    }
    assert _canonical_key(ex) == (
        'PR', 'lambdaT3', False, 'y1', 'y2', False,
        (('T3MEU11', False, 1),), 0,
    )
    assert Fraction(-2, 4) == Fraction(-1, 2)
    assert _tree_llhh_term_count('1 + 2') == 0
    import io
    fake = b'Y/LLHH'
    manifest = {'files': {'evidence.txt': {
        'size_bytes': len(fake), 'sha256': _hash(fake),
    }}}
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, 'w') as z:
        z.writestr('evidence.txt', fake)
    with zipfile.ZipFile(stream) as z:
        assert _archive_text(z, manifest, 'evidence.txt') == 'Y/LLHH'
        tampered = {'files': {'evidence.txt': {'size_bytes': len(fake),
                                              'sha256': '0' * 64}}}
        try:
            _archive_text(z, tampered, 'evidence.txt')
        except LedgerError:
            pass
        else:
            raise AssertionError('Tampered manifest accepted')
    print('PASS: exact keys, archive manifest hash guard, no CT certified')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            return 0
        if not args.archive or not args.output_dir:
            parser.error('--archive and --output-dir are required')
        report = build(args.archive, args.output_dir)
        for m, item in report['models'].items():
            print(f"{m}: 50 pole/log terms -> "
                  f"{item['neutral_pole_monomials']} neutral monomials; "
                  'required CT pole only; physical C5=NOT CERTIFIED')
        print('PASS: exact hard Y/LLHH pole ledger; report=' +
              str(args.output_dir / 'summary.json'))
        return 0
    except (LedgerError, ValueError, KeyError, RuntimeError,
            zipfile.BadZipFile) as exc:
        print('FAIL (no physical matching certificate): ' + str(exc))
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
