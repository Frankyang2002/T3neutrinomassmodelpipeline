"""Reconstruction tests; these do not execute physics or LaTeX."""
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'BuildCombinedAnalyticalReports.py'


def _module():
    spec = importlib.util.spec_from_file_location('BuildCombinedAnalyticalReports', SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_all_models_are_uniquely_named():
    m = _module()
    keys = [m.model_keys(x) for x in m.MODELS]
    assert len(keys) == 16
    assert len({a for a, b in keys}) == 16
    assert m.model_keys(('E',3,3,2,2)) == ('T3-E-p2','T3_dS1_3_dS2_3_dF_2_alpha_p2')


def test_reconstruct_all_saved_run_records(tmp_path):
    m = _module()
    for index, model in enumerate(m.MODELS):
        short, long = m.model_keys(model)
        folder = tmp_path / (short if index % 2 else long)
        (folder / 'T3_physical' / 'data').mkdir(parents=True)
        (folder / 'T3_physical' / 'data' / 'uv_rgbeta_rge.json').write_text('{}')
        (folder / 't3_model_comparison.json').write_text(json.dumps([{
            'UVRGEFile': 'data/uv_rgbeta_rge.json',
            'UVRGEStatus': 'Success', 'EFTStages': [],
        }]))
    records = m.reconstruct_records(tmp_path)
    assert len(records) == 16
    assert records[0].name == 'T3-A-m4'
    assert records[-1].name == 'T3-E-p2'


def test_missing_model_rejected(tmp_path):
    with pytest.raises(FileNotFoundError, match='Missing models'):
        _module().reconstruct_records(tmp_path)
