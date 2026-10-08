"""Exact electric-charge spectra for T3 models; no unverified CG basis mapping.

Input: component_support_audit.json from AuditT3ComponentMultiplicity.
Output: per-model SU(2) weight/charge spectra and a basis-mapping warning.

The Matchete d=3 CG array index is *not* asserted to equal ordered T3.
The spectra do NOT assign physical neutrino-mass multiplicities.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


DEFAULT_INPUT = Path('Reports/output/full/multiplicity/component_support_audit.json')
DEFAULT_OUTPUT = Path('Reports/output/full/multiplicity/charge_spectrum_audit.json')


def _fmt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f'{x.numerator}/{x.denominator}'


def charge_spectrum(dimension: int, hypercharge: Fraction) -> list[dict]:
    """List abstract isospin weights m=j,...,-j and Q=m+Y."""
    if not isinstance(dimension, int) or dimension < 1:
        raise ValueError('SU(2) dimension must be a positive integer')
    j = Fraction(dimension - 1, 2)
    return [{'m': _fmt(j - k), 'Q': _fmt(j - k + hypercharge)}
            for k in range(dimension)]


def audit_case(case: dict) -> dict:
    m = case['model']
    d1, d2, df, alpha = (int(m[k]) for k in ('dS1', 'dS2', 'dF', 'alpha'))
    y = {'S1': Fraction(alpha, 2), 'S2': Fraction(alpha + 2, 2),
         'F': Fraction(alpha + 1, 2), 'L': Fraction(-1, 2), 'H': Fraction(1, 2)}
    dimensions = {'S1': d1, 'S2': d2, 'F': df, 'L': 2, 'H': 2}
    fields = {
        name: {'dimension': dimensions[name], 'hypercharge': _fmt(y[name]),
               'weights_and_charges': charge_spectrum(dimensions[name], y[name]),
               'conjugate_charges': [_fmt(-Fraction(x['Q']))
                                     for x in charge_spectrum(dimensions[name], y[name])]}
        for name in ('L', 'H', 'F', 'S1', 'S2')
    }
    return {'model': m, 'fields': fields,
            'neutral_components': {name: [v['m'] for v in item['weights_and_charges'] if v['Q'] == '0']
                                   for name, item in fields.items()},
            'input_tensor_nonzero_entries': {name: int(v['nonzero_component_count'])
                                            for name, v in case['vertices'].items()},
            'physical_charge_basis_contraction': 'not_computed',
            'verified_loop_multiplicity': None,
            'compensation_factor': None}


def audit(payload: dict) -> dict:
    cases = payload.get('cases')
    if not isinstance(cases, list) or not cases:
        raise ValueError('Input must contain nonempty cases array')
    models = [audit_case(case) for case in cases]
    keys = [(p['model']['dS1'], p['model']['dS2'], p['model']['dF'], p['model']['alpha'])
            for p in models]
    if len(set(keys)) != len(keys):
        raise ValueError('Duplicated T3 model in input')
    return {'status': 'charge_spectra_verified_component_basis_unresolved',
            'count': len(models), 'models': models,
            'convention': 'Q = T3 + Y; SU(2) abstract weights j,j-1,...,-j; conjugates have opposite charges',
            'warning': ('These are spectra, not CG-array component labels. Matchete triplet indices may be '
                        'in a real Cartesian basis. Without the exact generator and basis rotation, '
                        'the neutrino loop cannot be charge-selected from CG array positions; '
                        'no multiplicity factor is assigned.')}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    p.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = p.parse_args(argv)
    if not args.input.is_file():
        p.error(f'Input not found: {args.input} (use --input to specify location)')
    result = audit(json.loads(args.input.read_text(encoding='utf-8-sig')))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(f'Charge spectra: {result["count"]} models; multiplicities assigned: 0')
    print('Triplet CG basis mapping: unresolved; no scalar coupling changed')
    print(f'JSON: {args.output.resolve()}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
