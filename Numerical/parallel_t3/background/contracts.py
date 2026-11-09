"""Fail-closed provenance and I/O contracts for Workstream B.

Only imports stdlib and numpy. The input symbolic results are never regenerated.
The GitHub baseline identifies public source for human comparison, whereas the
per-model local SHA256 fields are the actual reuse gate.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

import numpy as np

FAMILIES = {
    'A': ((1, 3, 2), (-4, -2, 0)),
    'B': ((2, 2, 1), (-3, -1, 1)),
    'C': ((2, 2, 3), (-3, -1, 1)),
    'D': ((3, 1, 2), (-2, 0, 2)),
    'E': ((3, 3, 2), (-4, -2, 0, 2)),
}
MODELS = tuple(f'{k}_{a:+d}' for k, (_, alphas) in FAMILIES.items() for a in alphas)
MIXED = frozenset(('B_-1', 'C_-1'))
SELF = frozenset(('A_-2', 'A_+0', 'D_-2', 'D_+0', 'E_-2', 'E_+0'))
HEAVY = frozenset(('MF', 'y1', 'y2', 'y1Conj', 'y2Conj'))
BASE = frozenset(('gY', 'g2', 'g3', 'yu', 'yd', 'ye', 'lambdaH',
                  'mS1Sq', 'mS2Sq', 'lambdaS1', 'lambdaS2',
                  'lambdaH1', 'lambdaH2', 'lambda12', 'lambdaT3'))
BETA_FORMAT = 'general-t3-quadratic-rgbeta-experimental-v1'
BETA_CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'
PUBLIC_HEAD_CHECKED = 'a944e0030b8ef24b23342d99e0d8f33cef3cffd5'
DEFAULT_BATCH1 = 'output/batch1_general_contract/batch1_native_v1'


class ContractError(ValueError):
    """Input is missing, stale or incompatible; never synthesize missing data."""


def tag(model: str) -> str:
    if model not in MODELS:
        raise ContractError(f'Unknown production T3 model {model!r}')
    return model.replace('-', 'm').replace('+', 'p')


def hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    data = Path(path).read_bytes()
    if data.startswith((b'\xff\xfe', b'\xfe\xff')):
        text = data.decode('utf-16')
    else:
        text = data.decode('utf-8-sig')
    obj = json.loads(text)
    if not isinstance(obj, dict):
        raise ContractError(f'Expected JSON object: {path}')
    return obj


def save_json(path: Path, obj: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.partial')
    temporary.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False)
                         + '\n', encoding='utf-8')
    temporary.replace(path)


def check_beta(model: str, stage: str, doc: Mapping[str, Any]) -> list[str]:
    if model not in MODELS or stage not in ('UV', 'EFT1'):
        raise ContractError(f'Invalid model/stage: {model}/{stage}')
    fam, alpha = model.split('_')
    meta = doc.get('metadata', {})
    needed = dict(zip(('dS1', 'dS2', 'dF'), FAMILIES[fam][0]))
    needed.update(alpha=int(alpha), GeneralScalarModel=model,
                  GeneralScalarStage=stage, SharedScalar=False,
                  ReportBetaConvention=BETA_CONVENTION)
    if doc.get('format') != BETA_FORMAT or doc.get('status') != 'SuccessCandidateMassBeta':
        raise ContractError(f'{model}/{stage}: not a successful full-general candidate export')
    if not isinstance(meta, dict) or any(meta.get(k) != v for k, v in needed.items()):
        raise ContractError(f'{model}/{stage}: representation/stage/beta convention mismatch')
    if meta.get('GeneralScalarBasisStatus') != 'CandidateNativeRGBetaUnmappedToMatchete':
        raise ContractError(f'{model}/{stage}: unexpected RGBeta basis status')
    betas = doc.get('report_betas')
    if not isinstance(betas, dict) or not betas or any(
            not isinstance(k, str) or not isinstance(v, str) or not v.strip()
            for k, v in betas.items()):
        raise ContractError(f'{model}/{stage}: absent or malformed beta expressions')
    keys = set(betas)
    if not BASE <= keys:
        raise ContractError(f'{model}/{stage}: missing required betas {sorted(BASE-keys)}')
    if stage == 'UV' and not {'MF', 'y1', 'y2'} <= keys:
        raise ContractError(f'{model}: incomplete UV heavy-fermion beta sector')
    if stage == 'EFT1' and keys & HEAVY:
        raise ContractError(f'{model}: fermion betas survive in EFT1: {sorted(keys & HEAVY)}')
    if set(doc.get('beta_symbols', [])) != keys:
        raise ContractError(f'{model}/{stage}: beta_symbols disagrees with report_betas')
    if not set(meta.get('GeneralScalarNewBetaKeys', [])) <= keys:
        raise ContractError(f'{model}/{stage}: missing extra general quartic betas')
    return list(betas)


def check_pair(model: str, uv: Mapping[str, Any], eft: Mapping[str, Any]) -> tuple[list[str], list[str]]:
    a, b = check_beta(model, 'UV', uv), check_beta(model, 'EFT1', eft)
    if set(a) - HEAVY != set(b):
        raise ContractError(f'{model}: nonfermion beta inventory changed across F threshold; '
                            f'UV-only={sorted(set(a)-set(b))}, EFT-only={sorted(set(b)-set(a))}')
    return a, b


def verify_batch1_provenance(root: Path, batch_root: Path, model: str,
                              record: Mapping[str, Any], uv_path: Path,
                              eft_path: Path, *, strict_sources: bool = True) -> dict:
    """Check actual local source SHA256 and cached symbolic output SHA256.

    This deliberately requires the source used for Batch 1, not just a Git
    commit. Never accept a mismatched cache because the public HEAD matches.
    """
    if record.get('model') != model:
        raise ContractError(f'{model}: Batch 1 model record identifies another model')
    expected = record.get('local_source_sha256')
    if not isinstance(expected, dict) or not expected:
        raise ContractError(f'{model}: Batch 1 has no source fingerprint')
    observed: dict[str, dict] = {}
    for relative, reference in expected.items():
        p = root / relative
        if not p.is_file():
            raise ContractError(f'{model}: missing source file {relative}')
        actual = hash_file(p)
        is_b_input = (relative.startswith('Lagrangian/') or
                      relative.startswith('RGE/running/rgbeta/') or
                      relative.startswith('RGE/general/generated/'))
        observed[relative] = {'sha256': actual, 'batch1_sha256': reference,
                              'same_as_batch1': actual == reference,
                              'required_for_B_reuse': is_b_input}
        # Workstream A can legitimately edit Wilson-only source files in
        # parallel. Such changes do not invalidate renormalisable betas or
        # the original Matchete F-tree seed; record them but do not block B.
        if strict_sources and is_b_input and actual != reference:
            raise ContractError(f'{model}: local source changed from successful Batch 1: '
                                f'{relative}; cached betas/matching must be regenerated')
    for stage, p in (('uv', uv_path), ('eft1', eft_path)):
        reference = record.get('rgbeta', {}).get(stage, {})
        if reference.get('status') != 'SuccessCandidateMassBeta':
            raise ContractError(f'{model}: Batch 1 {stage} export not recorded successful')
        if hash_file(p) != reference.get('sha256'):
            raise ContractError(f'{model}: Batch 1 cached {stage} JSON has changed')
    if model == 'A_-4':
        pot = root / 'Lagrangian/interactions/ScalarPotential.wl'
        if 'T3HiggsPortalHermitianPhase' not in pot.read_text(encoding='utf-8-sig'):
            raise ContractError('A_-4: anti-Hermitian portal correction absent')
        portal = record.get('f_threshold', {}).get('a4_portal', {})
        if not (portal.get('inv1_survives') and portal.get('inv2_survives')):
            raise ContractError('A_-4: Batch 1 is not based on corrected Matchete portal')
    return {'model': model, 'batch1_root':str(batch_root.relative_to(root))
            if batch_root.is_relative_to(root) else str(batch_root),
            'local_sources': observed,
            'uv_export_sha256': hash_file(uv_path),
            'eft1_export_sha256': hash_file(eft_path),
            'model_record_sha256': hash_file(batch_root/'models'/f'{tag(model)}.json'),
            'required_B_source_fingerprints_equal_batch1': all(
                v['same_as_batch1'] for v in observed.values() if v['required_for_B_reuse']),
            'public_main_last_inspected': PUBLIC_HEAD_CHECKED}


def encode_complex(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return encode_complex(value.tolist())
    if isinstance(value, (np.complexfloating, complex)):
        z = complex(value)
        return {'real': float(z.real), 'imag': float(z.imag)}
    if isinstance(value, np.generic):
        return encode_complex(value.item())
    if isinstance(value, (list, tuple)):
        return [encode_complex(v) for v in value]
    if isinstance(value, dict):
        return {str(k): encode_complex(v) for k, v in value.items()}
    return value


def check_finite(*arrays: np.ndarray) -> None:
    if any(not np.isfinite(np.asarray(a)).all() for a in arrays):
        raise ContractError('Nonfinite numerical state/tensor')


def safe_run_id(value: str) -> str:
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,79}', value):
        raise ContractError('run_id must be alphanumeric, underscore or dash')
    return value
