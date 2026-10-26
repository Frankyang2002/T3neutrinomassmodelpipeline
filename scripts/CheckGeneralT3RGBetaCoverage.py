"""Audit STRUCTURAL (not normalised) Matchete/RGBeta quartic coverage.

The two C5-relevant B/C(-1) structures are already present in RGBeta by
field content. THIS DOES NOT MEAN their coupling normalisations have been
matched. For generic model points the omitted extra structures do not have
RGBeta vertices or UV beta functions yet.

Usage:
    python scripts/CheckGeneralT3RGBetaCoverage.py --output output/general_rgbeta_coverage.json
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'docs/general_scalar_trial_operator_manifest.json'
RGBETA = ROOT / 'RGE/running/rgbeta/T3RGBetaModel.wl'

# This is a field-content match, NOT a basis, phase, factor, CG, or beta mapping.
# One of the two mixed H Hdag S1dag S2dag contractions is still missing.
B_MINUS_ONE_STRUCTURAL_MAP: dict[tuple[int, ...], list[str]] = {
    (0, 0, 0, 1, 1, 2): ['lambdaS1barS2S2bar2'],
    (0, 0, 0, 2, 0, 2): ['lambdaS1bar2S2bar2'],
    (0, 0, 1, 2, 0, 1): ['lambdaS1S1bar2S2bar'],
    (0, 2, 0, 0, 2, 0): ['lambdaHHdagS2S2'],
    (0, 2, 0, 2, 0, 0): ['lambdaHHdagS1barS1bar'],
    (1, 1, 0, 1, 0, 1): ['lambdaHHdagS1barS2barCross'],
}


def audit(manifest: dict, source: str) -> dict:
    models: dict[str, dict] = {}
    totals = {'additional_matchete_quartics': 0,
              'rgbeta_structurally_covered_quartics': 0,
              'rgbeta_missing_quartics': 0}
    for model, body in sorted(manifest['models'].items()):
        is_bc_m1 = model in ('B_-1', 'C_-1')
        rows = []
        for structure in body['extra_quartic_structures']:
            counts = tuple(structure['counts'])
            total = int(structure['invariants'])
            candidates = (B_MINUS_ONE_STRUCTURAL_MAP.get(counts, [])
                          if is_bc_m1 else [])
            found = [name for name in candidates if re.search(
                r'\b' + re.escape(name) + r'\b', source)]
            # There may be *more* SU2 invariants than mapped RGBeta tensors.
            covered = min(total, len(found))
            totals['additional_matchete_quartics'] += total
            totals['rgbeta_structurally_covered_quartics'] += covered
            totals['rgbeta_missing_quartics'] += total - covered
            rows.append({
                'structure': structure['field_structure'],
                'counts': counts,
                'matchete_independent_contractions': total,
                'rgbeta_symbol_candidates': found,
                'rgbeta_structural_contractions': covered,
                'rgbeta_unrepresented_contractions': total - covered,
                'normalisation_verified': False,
            })
        models[model] = {
            'matchete_extra_quartics': int(body['extra_quartic_complex_couplings']),
            'structural_coverage': rows,
        }
    return {
        'scope': 'All 16 ordinary supported models: only extra scalar quartics',
        'status': 'StructuralInventoryOnly',
        'totals': totals,
        'models': models,
        'caveats': [
            'Exact Matchete CG vs RGBeta invariant-tensor normalisation is not verified.',
            'The extra RGBeta quartics are not implemented for all other model classes.',
            'Only the B/C(-1) RGBeta special quartics are recognised by the map.',
            'The complex scalar mass mixing operators are not included in RGBeta or matching.',
        ]
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('output/general_rgbeta_coverage.json'))
    args = parser.parse_args()
    if not RGBETA.is_file():
        parser.error(f'RGBeta source missing: {RGBETA}; execute inside full repository')
    source = RGBETA.read_text(encoding='utf-8')
    manifest = json.loads(MANIFEST.read_text(encoding='utf-8'))
    data = audit(manifest, source)
    destination = args.output if args.output.is_absolute() else ROOT / args.output
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    print('Extra Matchete quartics:', data['totals']['additional_matchete_quartics'])
    print('RGBeta structural analogues:', data['totals']['rgbeta_structurally_covered_quartics'])
    print('Missing RGBeta structures:', data['totals']['rgbeta_missing_quartics'])
    print('Normalisations verified: 0 (requires explicit invariant-tensor comparison)')
    print('Report:', destination)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
