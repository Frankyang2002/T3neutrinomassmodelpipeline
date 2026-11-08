"""Audit A(-4) global RGBeta quartic projection before using it in EFT1.

SOURCE: RGE/running/rgbeta/RunGeneralT3RGBeta.wl EFT1 1 3 2 m4 A_-4

The general runner uses BetaTerm[Quartic,1], unlike the production
T3RGBetaEFT1OneLoopBetas individual-coupling extractions. A_-4 has no
extra quartics, but two *existing* quartics have four S2 legs and
must be projected together.  A successful status flag is not a certificate:
all 55 unit/pair quadratic scalar beta polynomials are tested against the
independent real-scalar fourth-derivative contraction before any compatible
beta input is emitted.  No physical Matchete coupling map is implied.

This file never changes the native RGBeta sources, earlier beta exports,
existing ODE results, Matchete/Wilson inputs, or threshold matching.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.running import Aminus4EFT1Couplings as eft
from Numerical.running import Aminus4EFT1Audit as audit

QUARTICS = frozenset({
    'lambdaH', 'lambdaS1', 'lambdaS2', 'lambdaS2Adj',
    'lambdaH1', 'lambdaH2', 'lambdaH2Adj', 'lambda12', 'lambdaT3',
})
NONQUARTICS = eft.BETA_NAMES - QUARTICS
COORDINATES = audit.SCALAR_NAMES


class ProjectionError(ValueError):
    """Unverified native export, source change or incompatible beta input."""


def require(condition, reason: str) -> None:
    if not condition:
        raise ProjectionError(reason)


def read_json(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(payload, dict), f'{path}: expected JSON object')
    return payload


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def check_source(root: Path) -> dict:
    paths = (
        'RGE/running/rgbeta/RunGeneralT3RGBeta.wl',
        'RGE/running/rgbeta/T3RGBetaModel.wl',
        'RGE/general/generated/A_-4_extra_quartics.wl',
    )
    src = {}
    for rel in paths:
        p = root / rel
        require(p.is_file(), f'Missing executable source: {rel}')
        src[rel] = p.read_text(encoding='utf-8-sig')
    general, model, extra = (src[p] for p in paths)
    require('BetaTerm[Quartic, 1]' in general and
            'GeneralScalarQuarticBetaBasis' in general and
            'RGBetaFullQuarticProjection' in general,
            'Global quartic projection source has changed')
    require('T3RGBetaEFT1OneLoopBetas[]' in general and
            'T3RGBetaBuildEFT1[dS1,dS2,dF,alpha]' in general,
            'Native A(-4) EFT1 build path has changed')
    require('T3RGBetaAddRGClosedQuartics[dS1, dS2, dF, alpha]' in model and
            'T3RGBetaAddBasicQuartics[dS1, dS2]' in model,
            'Native A(-4) scalar potential build path has changed')
    require('T3GeneralQuarticExpectedCount = 0;' in extra and
            'T3GeneralQuarticReusedCount = 0;' in extra,
            'A(-4) additional quartic inventory changed; do not assume zero')
    return {p: sha(root / p) for p in paths}


def validate_native(original: dict, projected: dict) -> dict:
    """Validate source stage and full inventory independently of status names."""
    eft.validate_beta_payload(original)
    require(projected.get('status') == 'Success', 'Global RGBeta export did not succeed')
    meta = projected.get('metadata')
    require(isinstance(meta, dict), 'Missing global projection metadata')
    baseline = original['metadata']
    for key in ('SharedScalar', 'dS1', 'dS2', 'dF', 'alpha', 'YS1', 'YS2',
                'ReportBetaConvention'):
        require(meta.get(key) == baseline.get(key), f'Metadata changed: {key}')
    require(meta.get('GeneralScalarStage') == 'EFT1' and
            meta.get('GeneralScalarModel') == eft.MODEL and
            meta.get('GeneralScalarQuarticBetaBasis') == 'RGBetaFullQuarticProjection',
            'Not the native A(-4) global EFT1 quartic projection')
    require(meta.get('GeneralScalarAdditionalQuartics') == 0 and
            meta.get('GeneralScalarExistingStructuralCount') == 0 and
            meta.get('GeneralScalarNewBetaKeys') == [],
            'A(-4) does not permit extra quartic symbols in this exporter')
    projected_keys = meta.get('GeneralScalarExistingQuarticBetasProjected')
    require(isinstance(projected_keys, list) and
            set(projected_keys) == QUARTICS and
            len(projected_keys) == len(QUARTICS),
            'Global projection did not cover all nine complex/real quartic symbols')
    beta = projected.get('report_betas')
    require(isinstance(beta, dict) and set(beta) == eft.BETA_NAMES,
            'Global RGBeta beta inventory must exactly equal all 17 active EFT1 parameters')
    require(all(isinstance(v,str) and v.strip() and '$Failed' not in v and
                'BetaTerm[' not in v for v in beta.values()),
            'Projected beta export contains unresolved expressions')
    require(all(original['report_betas'][name].replace(' ', '') ==
                beta[name].replace(' ', '') for name in NONQUARTICS),
            'Global projection unexpectedly changed gauge, SM Yukawa or scalar mass betas')
    require(not (set(beta) & {'MF', 'y1', 'y2', 'y1Conj', 'y2Conj'}),
            'Integrated-out heavy fermion retained in EFT1')
    return {
        'quartic_symbols': sorted(QUARTICS),
        'nonquartic_symbols_unchanged': sorted(NONQUARTICS),
        'changed_quartic_expressions': sorted(n for n in QUARTICS
            if original['report_betas'][n].replace(' ','') != beta[n].replace(' ','')),
    }


def compatibility_payload(projected: dict) -> dict:
    """Add only existing EFT1 stage metadata required by the old solver.

    This is NOT a tensor-basis conversion or physical Matchete certificate.
    """
    out = dict(projected)
    out['metadata'] = dict(projected['metadata'])
    out['metadata']['IntegratedField'] = 'F'
    out['metadata']['ActiveBSMFields'] = ['S1','S2']
    # Do not set 'format': this is the repository's standard 17-entry EFT1
    # state layout, not the experimental quadratic-mass extension.
    eft.validate_beta_payload(out)
    return out


def unit_point(names: tuple[str,...]) -> dict:
    p = {
        'gY': 0., 'g2': 0., 'g3': 0.,
        'yu': np.zeros((3,3), dtype=complex),
        'yd': np.zeros((3,3), dtype=complex),
        'ye': np.zeros((3,3), dtype=complex),
        'mS1Sq': 1.0e6, 'mS2Sq': 1.0e6,
    }
    p.update({k:(0j if k=='lambdaT3' else 0.) for k in eft.SCALAR_COUPLINGS})
    for name in names:
        if name == 'lambdaT3.real':p['lambdaT3'] += 1.
        elif name == 'lambdaT3.imag':p['lambdaT3'] += 1j
        else:p[name] += 1.
    return p


def quadratic_probes() -> list[tuple[str,...]]:
    """Ten squares plus 45 distinct cross terms determine a quadratic map."""
    require(len(COORDINATES) == 10 and len(set(COORDINATES)) == 10,
            'Changed scalar tensor coordinates')
    return [(name,) for name in COORDINATES] + list(combinations(COORDINATES,2))


def verify_scalar_polynomial(payload: dict, *, comparison=None) -> dict:
    if comparison is None:
        comparison = audit.scalar_beta_comparison
    failures = []
    max_closure = 0.
    for probe in quadratic_probes():
        row = comparison(unit_point(probe), payload)
        closure = row['scalar_loop_projection_closure_relative']
        max_closure = max(max_closure, float(closure))
        if (not row['reference_and_native_scalar_loop_match'] or
                closure >= 1.e-9):
            failures.append({
                'nonzero_unit_couplings': list(probe),
                'mismatched_betas': row['mismatched_coupling_betas'],
                'closure_relative': float(closure),
                'components': {name: row['coupling_beta_comparison'][name]
                               for name in row['mismatched_coupling_betas']},
            })
    return {
        'test': 'all 10 real-coupling squares and 45 pair cross terms',
        'scalar_quadratic_probe_count': 55,
        'maximum_reference_projection_closure_relative': max_closure,
        'failures': failures,
        'native_and_reference_scalar_quadratic_polynomials_agree':not failures,
    }


def evaluate(*, root: Path, old_path: Path, projected_path: Path,
             report_path: Path, compatible_path: Path) -> dict:
    report_path = report_path.resolve()
    compatible_path = compatible_path.resolve()
    for target in (report_path,compatible_path):
        require(target not in (old_path.resolve(),projected_path.resolve()),
                'Refusing to overwrite a native RGBeta source input')
    require(report_path != compatible_path, 'Distinct report/compatible paths required')
    source_hashes = check_source(root)
    original = read_json(old_path)
    projected = read_json(projected_path)
    inventory = validate_native(original, projected)
    compatible = compatibility_payload(projected)
    scalar = verify_scalar_polynomial(compatible)
    report = {
        'model': eft.MODEL,
        'format': 'a-m4-global-rgbeta-quartic-projection-audit-v1',
        'status': 'AllScalarPolynomialsAgree' if scalar[
            'native_and_reference_scalar_quadratic_polynomials_agree'] else
            'ProjectedScalarPolynomialsStillDisagree',
        'native_global_projection_source_sha256': source_hashes,
        'old_individual_beta_export_sha256': sha(old_path),
        'new_global_beta_export_sha256': sha(projected_path),
        'inventory': inventory,
        'scalar': scalar,
        'compatible_beta_file_created':bool(not scalar['failures']),
        'matchete_to_rgbeta_scalar_basis_certified':False,
        'full_seven_term_Wilson_running_certified':False,
        'beta_export_interpretation':
            'No fitted factors; global RGBeta quartic projection checked against '
            'reference 12-real-scalar one-loop scalar-only polynomial, not '
            'the native Matchete invariant tensor or full gauge/Yukawa quartic beta.',
    }
    report_path.parent.mkdir(parents=True,exist_ok=True)
    report_path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    if not scalar['failures']:
        compatible_path.parent.mkdir(parents=True,exist_ok=True)
        compatible_path.write_text(json.dumps(compatible,indent=2)+'\n',encoding='utf-8')
    print('Global quartic beta expressions changed:',
          ', '.join(inventory['changed_quartic_expressions']) or '(none)')
    print('Independent scalar polynomial probes:', 55)
    print('Disagreeing probes:', len(scalar['failures']))
    if scalar['failures']:
        for item in scalar['failures'][:8]:
            print('  ', '+'.join(item['nonzero_unit_couplings']),
                  ':', ', '.join(item['mismatched_betas']))
        print('No compatible beta file written: investigate native projection.')
    else:
        print('Scalar-loop polynomial agrees with reference in all ten directions.')
        print('Compatible beta file:', compatible_path)
    print('Report:',report_path)
    print('Matchete-to-RGBeta invariant map certified: NO')
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original-beta',type=Path,required=True)
    parser.add_argument('--projected-beta',type=Path,required=True)
    parser.add_argument('--audit-output',type=Path,required=True)
    parser.add_argument('--compatible-output',type=Path,required=True)
    args=parser.parse_args(argv)
    report=evaluate(root=ROOT,old_path=args.original_beta,
                    projected_path=args.projected_beta,
                    report_path=args.audit_output,
                    compatible_path=args.compatible_output)
    return 0 if report['compatible_beta_file_created'] else 1


if __name__=='__main__':
    raise SystemExit(main())
