"""Regression: negative alpha must not be interpreted as WolframScript options.

Run from the repository root with:
    python -m pytest -q tests/python/test_general_t3_eft1_alpha_argument.py

No Wolfram installation, RGBeta calculation, or Stage-33 data is required.
"""
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import RunGeneralT3EFT1Quartics as runner


CASES = [
    ('A_+0', 1, 3, 2, 0, '0'),
    ('A_-2', 1, 3, 2, -2, 'm2'),
    ('A_-4', 1, 3, 2, -4, 'm4'),
    ('B_+1', 2, 2, 1, 1, '1'),
    ('B_-1', 2, 2, 1, -1, 'm1'),
    ('B_-3', 2, 2, 1, -3, 'm3'),
    ('C_+1', 2, 2, 3, 1, '1'),
    ('C_-1', 2, 2, 3, -1, 'm1'),
    ('C_-3', 2, 2, 3, -3, 'm3'),
    ('D_+0', 3, 1, 2, 0, '0'),
    ('D_+2', 3, 1, 2, 2, '2'),
    ('D_-2', 3, 1, 2, -2, 'm2'),
    ('E_+0', 3, 3, 2, 0, '0'),
    ('E_+2', 3, 3, 2, 2, '2'),
    ('E_-2', 3, 3, 2, -2, 'm2'),
    ('E_-4', 3, 3, 2, -4, 'm4'),
]


@pytest.mark.parametrize('model,d_s1,d_s2,d_f,alpha,token', CASES)
def test_all_sixteen_models_pass_cli_safe_alpha(
    monkeypatch, tmp_path, model, d_s1, d_s2, d_f, alpha, token
):
    script = tmp_path / 'RunGeneralT3RGBeta.wl'
    script.touch()
    monkeypatch.setattr(runner, 'RGBETA_RUNNER', script)
    monkeypatch.setattr(runner, 'ROOT', tmp_path)
    seen = []

    def fake_run(command, **kwargs):
        seen.append((command, kwargs))
        Path(command[-1]).write_text('{}', encoding='utf-8')
        return SimpleNamespace(returncode=0, stdout='ok', stderr='')

    monkeypatch.setattr(runner.subprocess, 'run', fake_run)
    output = tmp_path / f'{model}_eft1.json'
    rep = SimpleNamespace(d_s1=d_s1, d_s2=d_s2, d_f=d_f, alpha=alpha)
    runner.export_fresh_eft1_betas(model, rep, output, wolframscript='wolframscript')

    assert output.is_file()
    assert len(seen) == 1
    args, opts = seen[0]
    assert len(args) == 10  # executable, -file, script, 7 script arguments
    assert args[:4] == ['wolframscript', '-file', str(script), 'EFT1']
    assert args[4:8] == [str(d_s1), str(d_s2), str(d_f), token]
    assert args[8:] == [model, str(output)]
    assert opts['cwd'] == tmp_path
    assert opts['capture_output'] and opts['text']


def test_failed_wolfram_command_reports_invocation(monkeypatch, tmp_path):
    script = tmp_path / 'RunGeneralT3RGBeta.wl'
    script.touch()
    monkeypatch.setattr(runner, 'RGBETA_RUNNER', script)
    monkeypatch.setattr(runner, 'ROOT', tmp_path)
    monkeypatch.setattr(
        runner.subprocess, 'run',
        lambda *args, **kwargs: SimpleNamespace(returncode=2, stdout='Usage:', stderr=''),
    )
    with pytest.raises(ValueError, match=r'A_-2.*exit=2') as exc:
        runner.export_fresh_eft1_betas(
            'A_-2', SimpleNamespace(d_s1=1, d_s2=3, d_f=2, alpha=-2),
            tmp_path / 'failed.json', wolframscript='wolframscript',
        )
    assert "'m2'" in str(exc.value)
    assert 'Usage:' in str(exc.value)
