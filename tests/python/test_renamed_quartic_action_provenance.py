"""Regression checks for filename-only provenance compatibility in saved data."""
from copy import deepcopy

import pytest

from scripts.NativeQuarticActionAudit import (
    HISTORICAL_EXPORT_PROVENANCE,
    RENAMED_EXPORT_PROVENANCE,
    verify_stored_stage23,
)


def records():
    fresh = {
        'source_verified': {
            'matchete_source': RENAMED_EXPORT_PROVENANCE,
            'upstream_rgbeta_source': 'unchanged',
        },
        'source_prefactors_verified': True,
        'signed_physical_coupling_map_certified': False,
        'totals': {'models': 16, 'slots': 50, 'blocks': 46},
        'blocks': [{'model': 'A_+0', 'rank': 1}],
    }
    saved = deepcopy(fresh)
    saved['source_verified']['matchete_source'] = HISTORICAL_EXPORT_PROVENANCE
    return saved, fresh


def test_exact_report_is_valid():
    _, fresh = records()
    verify_stored_stage23(fresh, fresh)


def test_only_historical_export_path_change_is_valid():
    saved, fresh = records()
    verify_stored_stage23(saved, fresh)


def test_changed_physics_block_is_rejected():
    saved, fresh = records()
    saved['blocks'][0]['rank'] = 2
    with pytest.raises(ValueError, match='stale'):
        verify_stored_stage23(saved, fresh)


def test_other_provenance_change_is_rejected():
    saved, fresh = records()
    saved['source_verified']['upstream_rgbeta_source'] = 'different source'
    with pytest.raises(ValueError, match='stale'):
        verify_stored_stage23(saved, fresh)


def test_unknown_source_rename_is_rejected():
    saved, fresh = records()
    saved['source_verified']['matchete_source'] = 'an unverified legacy path'
    with pytest.raises(ValueError, match='stale'):
        verify_stored_stage23(saved, fresh)


def test_missing_source_metadata_is_rejected():
    saved, fresh = records()
    saved.pop('source_verified')
    with pytest.raises(ValueError, match='stale'):
        verify_stored_stage23(saved, fresh)
