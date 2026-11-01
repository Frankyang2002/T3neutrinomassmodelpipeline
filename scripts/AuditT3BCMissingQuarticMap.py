"""Audit the two B/C(-1) missing Matchete/RGBeta quartic maps.

Read-only, source-evidence-based. This reports CONDITIONAL signed factors;
it does not claim the native component basis or physical map is certified.
Run: python -X utf8 -m scripts.AuditT3BCMissingQuarticMap
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sympy as sp
from scripts.QuarticPrefactorAudit import INPUTS, run as audit_native

ROOT = Path(__file__).resolve().parents[1]
PAIRS = {
    'lambdaT3ExtraQ020020I1': 'lambdaHHdagS2S2',
    'lambdaT3ExtraQ020200I1': 'lambdaHHdagS1barS1bar',
}
MODELS = ('B_-1', 'C_-1')


def run(folder: Path) -> dict:
    missing = [str(folder / name) for name in INPUTS if not (folder / name).is_file()]
    if missing:
        raise FileNotFoundError('Native tensor evidence unavailable; no factors may be inferred. Missing:\n  '
                                + '\n  '.join(missing))
    evidence = audit_native(base=folder)
    rows = evidence['blocks']
    found = []
    for model in MODELS:
        for mat, rgb in PAIRS.items():
            matches = [r for r in rows if r['model'] == model and
                       mat in r['matchete_couplings'] and rgb in r['rgbeta_couplings']]
            if len(matches) != 1:
                raise ValueError(f'{model}: unique native geometry block missing for {mat}/{rgb}')
            block = matches[0]
            mnames = block['matchete_couplings']
            rnames = block['rgbeta_couplings']
            a = sp.Matrix([[sp.sympify(x) for x in row] for row in block['r_in_m']])
            # Geometry: invariant_RGB = invariant_MAT * A.  L_MAT = -V_RGB.
            # Consequently lambda_MAT = - A * lambda_RGB (for aligned fields).
            mapping = -a
            if mapping.shape != (len(mnames), len(rnames)):
                raise ValueError('Inconsistent coupling-map shape')
            idx = mnames.index(mat)
            for j, other in enumerate(rnames):
                if other != rgb and sp.simplify(mapping[idx, j]) != 0:
                    raise ValueError(f'{model}/{mat}: coupled block; not a single-factor map')
            factor = sp.simplify(mapping[idx, rnames.index(rgb)])
            if factor == 0:
                raise ValueError(f'{model}/{mat}: zero geometric mapping')
            found.append({'model': model, 'matchete': mat, 'rgbeta': rgb,
                          'conditional_matchete_from_rgbeta_factor': str(factor),
                          'geometry_class': block['geometry_class'],
                          'physical_mapping_certified': False})
    return {'status': 'ConditionalNativeTensorMapAudited',
            'physical_mapping_certified': False,
            'source_directory': str(folder), 'entries': found,
            'note': 'Factors depend on component-basis and RGBeta potential-sign identification; '
                    'not permission to substitute into a physical C5 calculation.'}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--native-dir', type=Path, default=ROOT/'output/general_t3_stage22')
    ap.add_argument('--output-json', type=Path,
                    default=ROOT/'output/audits/t3_bc_missing_quartic_map.json')
    args = ap.parse_args()
    result = run(args.native_dir)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n',
                                encoding='utf-8')
    for row in result['entries']:
        print(f"{row['model']}: {row['matchete']} = ({row['conditional_matchete_from_rgbeta_factor']}) "
              f"* {row['rgbeta']} [CONDITIONAL]")
    print('Physical map NOT certified; output:', args.output_json)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
