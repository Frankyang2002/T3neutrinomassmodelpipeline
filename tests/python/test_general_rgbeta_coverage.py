"""Structural coverage only; no proof of tensor normalisation or beta correctness."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.CheckGeneralT3RGBetaCoverage import audit, B_MINUS_ONE_STRUCTURAL_MAP  # noqa


def _manifest():
    return json.loads((ROOT / 'docs/general_scalar_trial_operator_manifest.json').read_text())


def test_mapping_is_exactly_six_special_bm1_structures():
    body = _manifest()['models']['B_-1']
    keys = {tuple(entry['counts']) for entry in body['extra_quartic_structures']}
    assert keys == set(B_MINUS_ONE_STRUCTURAL_MAP)
    assert sum(x['invariants'] for x in body['extra_quartic_structures']) == 7


def test_rgbeta_structural_counts_not_normalisation():
    source = '\n'.join(sum(B_MINUS_ONE_STRUCTURAL_MAP.values(), []))
    result = audit(_manifest(), source)
    assert result['totals'] == {
        'additional_matchete_quartics': 50,
        'rgbeta_structurally_covered_quartics': 12,
        'rgbeta_missing_quartics': 38,
    }
    assert not any(row['normalisation_verified']
                   for value in result['models'].values()
                   for row in value['structural_coverage'])
    assert sum(row['rgbeta_unrepresented_contractions']
               for row in result['models']['B_-1']['structural_coverage']) == 1
    assert sum(row['rgbeta_unrepresented_contractions']
               for row in result['models']['C_-1']['structural_coverage']) == 1


def test_missing_rg_source_must_not_claim_present():
    result = audit(_manifest(), '')
    assert result['totals']['rgbeta_structurally_covered_quartics'] == 0
    assert result['totals']['rgbeta_missing_quartics'] == 50
