"""Optional integration smoke test against the user's existing Batch 1 outputs.
Run with T3_B_RUN_NATIVE_TESTS=1; never generates Wolfram symbolic files.
"""
from __future__ import annotations
import json
import os
from pathlib import Path
import pytest
from Numerical.parallel_t3.background.contracts import (
    MODELS, DEFAULT_BATCH1, check_pair, hash_file, read_json, tag)

ROOT = Path(__file__).resolve().parents[1]

@pytest.mark.skipif(os.getenv('T3_B_RUN_NATIVE_TESTS')!='1',
                    reason='opt-in native Batch 1 reuse check')
@pytest.mark.parametrize('model',MODELS)
def test_all_beta_pairs_native(model):
    b=ROOT/DEFAULT_BATCH1/'rgbeta'/tag(model)
    u=read_json(b/'uv.json'); e=read_json(b/'eft1.json')
    a,c=check_pair(model,u,e)
    assert len(a)>len(c)

@pytest.mark.skipif(os.getenv('T3_B_RUN_NATIVE_TESTS')!='1',
                    reason='opt-in native Batch 1 reuse check')
def test_corrected_a4_source_present():
    p=ROOT/'Lagrangian/interactions/ScalarPotential.wl'
    assert 'T3HiggsPortalHermitianPhase' in p.read_text(encoding='utf-8-sig')
