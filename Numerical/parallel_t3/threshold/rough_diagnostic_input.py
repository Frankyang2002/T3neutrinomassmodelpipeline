"""C-owned, explicit non-Weinberg projection *diagnostic* adapter.

No input of Workstreams 1/2 is altered.  C's optional redundancy key
``C5_EFT1_ms`` is removed only from a disposable C-owned NPZ if the source
full LLHH tensor fails Weinberg purity AND that redundant C5 differs from
its least-squares projection.  The discrepancy is fully documented, not
presented as a valid SMEFT coefficient.  Genuine benchmark mode NEVER calls
this function.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import numpy as np

from .rough_adapter import PreparedPair
from .pipeline import W_FORMAT, _load_pair
from .matching import master_weinberg_template, extract_lepton_sectors, project_master, symmetric3
from .models import n_real, normalize_model


def make_projection_only_pair(model: str, pair: PreparedPair, output: Path,
                              *, observed_purity_residual: float) -> tuple[PreparedPair, dict]:
    """Prepare C copy solely for ``--allow-projection-residual`` execution.

    The original Wilson tensor remains byte-identical.  No leading-log,
    finite matching or flavour coefficient is generated here.  The optional
    independent C5 redundancy cross-check is skipped *only* when it is
    incompatible with the already impure LLHH tensor.  Its original value,
    hashes and disagreement are retained in metadata and a JSON sidecar.
    """
    model = normalize_model(model)
    if not np.isfinite(observed_purity_residual) or observed_purity_residual <= 1e-8:
        return pair, {'performed': False, 'reason': 'Pure input, or non-Weinberg projection not required'}
    w, arrays, original_manifest_sha, original_npz_sha = _load_pair(pair.wilson, W_FORMAT)
    if normalize_model(w['model']) != model:
        raise ValueError('Projection diagnostic model disagrees with A source')
    if 'C5_EFT1_ms' not in arrays:
        return pair, {'performed': False, 'reason': 'No redundant C5_EFT1_ms supplied'}
    n = n_real(model)
    hh, _ = extract_lepton_sectors(arrays['C_psi2phi2_ms'],
        np.asarray(w['lepton_indices'], int), np.asarray(w['higgs_source_indices'], int),
        np.asarray(w['heavy_source_indices'], int), number_of_real_scalars=n+4)
    c5, purity = project_master(hh, enforce_pure=False)
    if not np.isclose(purity, observed_purity_residual, rtol=1e-7, atol=1e-10):
        raise ValueError('Diagnostic purity residual differs from reloaded source; possible stale data')
    supplied = symmetric3(arrays['C5_EFT1_ms'], 'separate C5_EFT1_ms')
    mismatch = float(np.linalg.norm(supplied - c5) /
                     max(float(np.linalg.norm(supplied)), float(np.linalg.norm(c5)),1e-100))
    if mismatch <= 1e-8:
        return pair, {'performed': False, 'reason': 'Source redundant C5 agrees with LS projector',
                      'relative_disagreement': mismatch}
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    # Always use the original adapted C tensor; do not try to replace LLHH
    # with the projected Weinberg polynomial here.
    saved = {k:v for k,v in arrays.items() if k != 'C5_EFT1_ms'}
    npz = output/'wilson_at_MS.npz'
    np.savez_compressed(npz, **saved)
    metadata = {**w, 'tensor_npz':'wilson_at_MS.npz',
        'tensor_npz_sha256': sha256(npz.read_bytes()).hexdigest(),
        'nonweinberg_projection_diagnostic_only': True,
        'source_redundant_C5_not_used_for_diagnostic_only': True}
    meta_path = output/'wilson.json'
    meta_path.write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    audit = {
        'format': 't3-nonweinberg-projection-input-override-v1',
        'classification': 'NonWeinbergProjectedDiagnostic_ONLY',
        'performed': True, 'model': model,
        'reason': 'Full evolved LLHH has non-Weinberg components; separate C5 is not its LS projection',
        'original_A_manifest_sha256': original_manifest_sha,
        'original_A_npz_sha256': original_npz_sha,
        'temporary_C_npz_sha256': metadata['tensor_npz_sha256'],
        'original_redundant_C5_EFT1_ms': {'real': supplied.real.tolist(), 'imag':supplied.imag.tolist()},
        'LS_projected_C5_gev_inverse': {'real': c5.real.tolist(), 'imag':c5.imag.tolist()},
        'relative_disagreement': mismatch,
        'Weinberg_purity_residual': purity,
        'original_wilson_tensor_changed': False,
        'independent_redundancy_check_deactivated_for_diagnostic_only': True,
        'cannot_be_classified_as_gauge_invariant_SMEFT_match': True,
    }
    (output/'redundant_C5_diagnostic.json').write_text(
        json.dumps(audit,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    transformed = PreparedPair(meta_path, pair.background,
                              {**pair.provenance,
                               'projected_diagnostic_input_only': audit})
    return transformed, audit
