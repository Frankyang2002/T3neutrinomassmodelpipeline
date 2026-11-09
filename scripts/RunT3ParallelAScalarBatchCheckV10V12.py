"""Fail-closed, read-only resume verification for the v10-v12 scalar batch.

Used only to check previously generated files, not to create/modify physics
outputs. A saved v11/v12 witness is independently recalculated before reuse.
Native v10 smoke and completed witnesses are revalidated against CURRENT local
sources and installed RGBeta provenance, never trusted from their names alone.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys


class ResumeBlocked(ValueError):
    pass


def compare_saved_to_fresh(saved, fresh, *, where: str = 'root') -> None:
    """Compare entire nested diagnostic records; tolerate float roundoff only."""
    if type(saved) is not type(fresh):
        raise ResumeBlocked(f'{where}: stored/fresh value type differs')
    if isinstance(fresh, dict):
        if set(saved) != set(fresh):
            raise ResumeBlocked(f'{where}: stored/fresh keys differ')
        for key in fresh:
            compare_saved_to_fresh(saved[key], fresh[key], where=f'{where}.{key}')
    elif isinstance(fresh, list):
        if len(saved) != len(fresh):
            raise ResumeBlocked(f'{where}: stored/fresh list lengths differ')
        for i, (s, f) in enumerate(zip(saved, fresh)):
            compare_saved_to_fresh(s, f, where=f'{where}[{i}]')
    elif isinstance(fresh, (int, float)) and not isinstance(fresh, bool):
        if not math.isfinite(float(saved)) or not math.isfinite(float(fresh)):
            raise ResumeBlocked(f'{where}: NaN or infinity is not allowed')
        if not math.isclose(float(saved), float(fresh), rel_tol=1e-10, abs_tol=1e-11):
            raise ResumeBlocked(f'{where}: saved numeric evidence differs from fresh audit')
    elif saved != fresh:
        raise ResumeBlocked(f'{where}: saved evidence differs from fresh audit')


def read(path: Path) -> dict:
    if not path.is_file():
        raise ResumeBlocked(f'missing prior diagnostic: {path}')
    record = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(record, dict):
        raise ResumeBlocked(f'not a JSON object: {path}')
    return record


def verify_python_audit(root: Path, which: str) -> dict:
    if which == 'v11':
        from Numerical.parallel_t3.wilson.scalar_nonportal_rank_v11 import audit
        path = root / 'output/parallel_t3/wilson/nonportal_scalar_invariants_conditional_v11.json'
    elif which == 'v12':
        from Numerical.parallel_t3.wilson.scalar_topology_action_v12 import audit
        path = root / 'output/parallel_t3/wilson/topology_scalar_action_conditional_v12.json'
    else:
        raise ResumeBlocked(f'unsupported mathematical audit: {which}')
    fresh = audit(root=root)
    compare_saved_to_fresh(read(path), fresh, where=which)
    return {'which': which, 'path': str(path), 'status': fresh['status']}


def verify_smoke(native: dict, root: Path) -> dict:
    from Numerical.parallel_t3.wilson.native_quartic_projectors_v10 import (
        FORMAT, SOURCE_FILES, DIMENSIONS, expected_symbols, file_hash,
    )
    import numpy as np
    if native.get('format') != FORMAT or native.get('smoke_only') is not True:
        raise ResumeBlocked('v10 smoke: wrong native format or not smoke')
    if native.get('status') != 'NativeRGBetaQuarticProjectorSmoke3ModelsPassed':
        raise ResumeBlocked('v10 smoke: native contractions did not pass')
    if native.get('registered_model_count') != 3 or native.get('native_Lam_evaluated') is not True or \
       native.get('native_quartic_projectors_evaluated') is not True:
        raise ResumeBlocked('v10 smoke: missing actual native execution evidence')
    if set(native.get('source_sha256', {})) != set(SOURCE_FILES):
        raise ResumeBlocked('v10 smoke: source hash inventory changed')
    for name in SOURCE_FILES:
        p = root / name
        if not p.is_file() or native.get('source_sha256', {}).get(name) != file_hash(p):
            raise ResumeBlocked('v10 smoke: stale source hash ' + name)
    ep = native.get('installed_RGBeta_main_file')
    if not isinstance(ep, str) or not Path(ep).is_file():
        raise ResumeBlocked('v10 smoke: installed RGBeta entrypoint missing')
    if native.get('installed_RGBeta_main_file_sha256') != file_hash(Path(ep)):
        raise ResumeBlocked('v10 smoke: installed RGBeta entrypoint changed')
    models = native.get('models')
    if not isinstance(models, dict) or set(models) != {'A_-4','B_-1','E_-2'}:
        raise ResumeBlocked('v10 smoke: wrong model list')
    for tag, row in models.items():
        expected = expected_symbols(tag)
        if row.get('ok') is not True or row.get('registry_exact_source_inventory') is not True:
            raise ResumeBlocked(f'v10 smoke: {tag} not passed')
        d1, d2, df = DIMENSIONS[tag[0]]
        if (row.get('dimensions') != [d1,d2,df] or
                row.get('alpha') != int(tag.split('_')[1]) or
                row.get('registered_quartics_count') != len(expected)):
            raise ResumeBlocked(f'v10 smoke: {tag} representation/count changed')
        inventory = row.get('inventory', [])
        if len(inventory) != len(expected) or {x.get('symbol') for x in inventory} != set(expected):
            raise ResumeBlocked(f'v10 smoke: {tag} inventory changed')
        for rec in inventory:
            sym = rec['symbol']
            if (rec.get('native_self_conjugate') is not expected[sym] or
                    rec.get('has_installed_RGBeta_projector') is not True or
                    not (1 <= rec.get('native_unique_arrangement_count', 0) <= 24) or
                    ''.join(rec.get('native_coupling_inputform','').split()) !=
                    f'6*Matrix[{sym}][##1]&'):
                raise ResumeBlocked(f'v10 smoke: {tag}/{sym} native registration changed')
        diagonals = row.get('native_real_self_projectors', [])
        if {x.get('symbol') for x in diagonals} != {k for k, v in expected.items() if v}:
            raise ResumeBlocked(f'v10 smoke: {tag} real diagonal coverage missing')
        if any(x.get('native_response', {}).get('ok') is not True or
               abs(x['native_response'].get('value', 999) - 1.) > 2e-8
               for x in diagonals):
            raise ResumeBlocked(f'v10 smoke: {tag} non-unit diagonal')
        blocks = row.get('native_portal_response_blocks', [])
        if len(blocks) != 2:
            raise ResumeBlocked(f'v10 smoke: {tag} portal response incomplete')
        for which, block in enumerate(blocks, start=1):
            dim = d1 if which == 1 else d2
            names = [f'lambdaH{which}'] + ([f'lambdaH{which}Adj'] if dim>1 else [])
            raw = block.get('native_response_matrix', [])
            n = len(names)
            if (block.get('symbols') != names or block.get('which') != which or
                    block.get('dimension') != dim or len(raw) != n or
                    any(len(z) != n for z in raw)):
                raise ResumeBlocked(f'v10 smoke: {tag}/S{which} response dimension changed')
            if any(x.get('ok') is not True for group in raw for x in group):
                raise ResumeBlocked(f'v10 smoke: {tag}/S{which} native response failed')
            response = np.array([[x.get('value') for x in group] for group in raw],float)
            if (not np.all(np.isfinite(response)) or
                    np.max(abs(np.diag(response)-1)) > 2e-8 or
                    np.linalg.svd(response, compute_uv=False)[-1] < 1e-8):
                raise ResumeBlocked(f'v10 smoke: {tag}/S{which} response singular')
        if (row.get('native_real_self_projectors_all_unit') is not True or
                row.get('native_portal_response_blocks_evaluated') is not True or
                row.get('native_numeric_S2_components_measured') is not False):
            raise ResumeBlocked(f'v10 smoke: {tag} coverage/certification changed')
    if native.get('physical_Wilson_evolution_authorized') is not False:
        raise ResumeBlocked('v10 smoke: unsupported physics authorization')
    return {'which': 'v10-smoke', 'status': native['status']}


def verify_native(root: Path, which: str) -> dict:
    base = root / 'output/parallel_t3/wilson'
    if which == 'v10-smoke':
        path = base / 'rgbeta_quartic_projectors_native_v10_smoke.json'
        result = verify_smoke(read(path), root)
    elif which == 'v10-native':
        from Numerical.parallel_t3.wilson.native_quartic_projectors_v10 import validate_native_evidence
        path = base / 'rgbeta_quartic_projectors_native_v10.json'
        result = validate_native_evidence(read(path), root)
        result = {'which': 'v10-native', 'status': 'native full evidence independently validated',
                  'registered_quartic_slots': result['registered_quartic_slots']}
    elif which == 'v10-final':
        from Numerical.parallel_t3.wilson.native_quartic_projectors_v10 import audit_from_paths
        path = base / 'native_quartic_projectors_conditional_v10.json'
        fresh = audit_from_paths(root,base/'rgbeta_quartic_projectors_native_v10.json',
                                 base/'portal_real_action_conditional_v9.json')
        compare_saved_to_fresh(read(path), fresh, where='v10-final')
        result = {'which': 'v10-final', 'status': fresh['status']}
    else:
        raise ResumeBlocked('unsupported native audit: ' + which)
    result['path'] = str(path)
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=Path('.'))
    parser.add_argument('--which', required=True,
                        choices=('v11','v12','v10-smoke','v10-native','v10-final'))
    args = parser.parse_args(argv)
    try:
        root = args.project_root.resolve()
        report = (verify_python_audit(root,args.which) if args.which in {'v11','v12'}
                  else verify_native(root,args.which))
        print(f"{args.which}: saved output independently VERIFIED and reusable: {report['path']}")
        return 0
    except (ResumeBlocked,ValueError,OSError,KeyError,TypeError,IndexError,
            ArithmeticError,ImportError) as exc:
        print(f'{args.which}: saved output BLOCKED: {exc}',file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
