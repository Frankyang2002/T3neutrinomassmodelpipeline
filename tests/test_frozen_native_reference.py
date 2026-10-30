from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]

def test_reference_complete():
    d=json.loads((ROOT/'docs/FROZEN_LAMBDA_TARGET_REFERENCE.json').read_text())
    assert len(d['target']['observable_names'])==5
    assert len(d['target']['covariance'])==5
    assert d['future_sigma']['sin2_theta12']>0

def test_native_uses_original_template_and_overlays():
    src=(ROOT/'scripts/BuildT3FrozenLambdaNative.py').read_text()
    assert 'extract_lambda_template' in src
    assert "'future_sigma':ref['future_sigma']" in src
    assert "'target':target" in src
    assert 'coverage.csv' in src
