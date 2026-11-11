"""Fail-closed Stage 22 -> full-quartic EFT1 production boundary.

A response JSON is NOT certified merely because it declares a provenance string.
This module requires an explicit post-Stage22 response-to-certificate binding,
including a SHA256 of the exact response bytes and a per-model verification.
It does not construct or certify physical tensors independently.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Mapping

from Numerical.running.T3RunningQuarticTensor import (
    QuarticEvidenceError, RunningQuarticBasis, expected_scalar_symbols,
)

CERT_STATUS = 'ExactTensorAndExplicitPrefactorAlgebraPassed'
CERT_BLOCK = 'ExactAlgebraPassed'
BINDING_FORMAT = 't3-quartic-response-certificate-binding-v1'


def _load_json(path: Path) -> dict:
    try:
        obj = json.loads(path.read_text(encoding='utf-8-sig'))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise QuarticEvidenceError(f'Unable to load {path}: {exc}') from exc
    if not isinstance(obj, dict):
        raise QuarticEvidenceError(f'Expected JSON object: {path}')
    return obj


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def registry_baseline_symbols(registry: Mapping, model: str,
                              baseline_symbols: set[str]) -> set[str]:
    """Require registry extra couplings and all supplied baseline couplings."""
    return expected_scalar_symbols(registry, model, baseline_symbols)


def audit_gate(registry_path: Path, stage22_path: Path,
               response_path: Path, binding_path: Path,
               model: str, baseline_symbols: set[str]) -> dict:
    """Check all pieces of evidence, without integrating or changing outputs."""
    registry = _load_json(registry_path)
    stage22 = _load_json(stage22_path)
    response = _load_json(response_path)
    binding = _load_json(binding_path)
    expected = registry_baseline_symbols(registry, model, baseline_symbols)
    if stage22.get('normalisation_certified') is not True or \
       stage22.get('certification_status') != CERT_STATUS:
        raise QuarticEvidenceError('Stage 22 has not certified native coupling normalization')
    cert_models = stage22.get('certification', {}).get('models', {})
    row = cert_models.get(model)
    if not isinstance(row, dict) or row.get('status') != CERT_BLOCK:
        raise QuarticEvidenceError(f'Stage 22 model {model} not certified')
    blocks = row.get('blocks', [])
    if not isinstance(blocks, list) or any(b.get('status') != CERT_BLOCK for b in blocks):
        raise QuarticEvidenceError(f'Stage 22 {model}: incomplete certified block list')
    if row.get('slots') != registry['models'][model]['total_additional_quartics']:
        raise QuarticEvidenceError(f'Stage 22 {model}: slot inventory mismatch')
    if binding.get('format') != BINDING_FORMAT or binding.get('model') != model:
        raise QuarticEvidenceError('Incorrect response certificate binding format/model')
    if binding.get('status') != 'IndependentlyVerifiedResponseMapping':
        raise QuarticEvidenceError('Response mapping not independently verified')
    if binding.get('stage22_sha256') != _sha(stage22_path):
        raise QuarticEvidenceError('Stage 22 report changed after response verification')
    if binding.get('response_sha256') != _sha(response_path):
        raise QuarticEvidenceError('Response tensor bytes changed after verification')
    if binding.get('registry_sha256') != _sha(registry_path):
        raise QuarticEvidenceError('Quartic registry changed after verification')
    checks = binding.get('checks', {})
    required_checks = ('field_component_alignment','potential_prefactors',
                       'rgbeta_to_real_cartesian','coupling_direction_completeness')
    if not isinstance(checks, dict) or any(checks.get(k) is not True for k in required_checks):
        raise QuarticEvidenceError('Missing independent real-Cartesian response checks')
    basis = RunningQuarticBasis.from_payload(response, expected)
    if basis.model != model:
        raise QuarticEvidenceError('Tensor response model mismatch')
    return {'model':model, 'status':'ResponseEvidenceBound',
            'response_sha256':_sha(response_path),
            'stage22_sha256':_sha(stage22_path),
            'registry_sha256':_sha(registry_path),
            'coupling_count':len(expected),
            'response_directions':len(basis.directions),
            'basis':basis}


def integrate_certified_eft1(adapter, start_ren, initial_wilson, *,
                             registry_path: Path, stage22_path: Path,
                             response_path: Path, binding_path: Path,
                             baseline_symbols: set[str], **solver_options):
    """A full-quartic production entry point gated on immutable evidence.

    Keeps the seven-term RHS and single DOP853 joint IVP of Batch 4.
    Existing unguarded synthetic solver remains available for tests.
    """
    from Numerical.running.T3ExtendedSimultaneousEFT1 import integrate_extended_eft1
    evidence = audit_gate(Path(registry_path), Path(stage22_path),
                          Path(response_path), Path(binding_path),
                          adapter.model, set(baseline_symbols))
    if set(adapter.required_quartics) != set(evidence['basis'].symbols):
        raise QuarticEvidenceError('Adapter/state quartic list differs from certified basis')
    if solver_options.get('check_hhss') is False:
        raise QuarticEvidenceError('HHSS validation cannot be disabled in production')
    result = integrate_extended_eft1(adapter, start_ren, initial_wilson,
                                     evidence['basis'], **solver_options)
    result['quartic_normalisation_certified'] = True
    result['quartic_response_sha256'] = evidence['response_sha256']
    result['quartic_stage22_sha256'] = evidence['stage22_sha256']
    return result
