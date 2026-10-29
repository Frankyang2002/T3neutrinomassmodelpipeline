"""Offline tests for the descriptive-name migration tool."""
import importlib.util
from pathlib import Path
import subprocess

TOOL = Path(__file__).resolve().parents[2] / 'tools/RenameNumberedResearchFiles.py'
SPEC = importlib.util.spec_from_file_location('local_name_migration', TOOL)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_production_modules_replaced():
    assert mod.rename_name('scripts/Stage31CoupledUVRunning.py') == 'scripts/CoupledUVRunning.py'
    assert mod.rename_name('scripts/Stage34IntermediateEFTRunning.py') == 'scripts/IntermediateEFTRunning.py'
    assert mod.rename_name('scripts/RunGeneralT3Stage18.py') == 'scripts/RunGeneralT3EFT1WilsonSourceTransport.py'
    assert mod.rename_name('scripts/RunGeneralT3Stage10.ps1') == 'scripts/RunGeneralT3MatchingRGSourcePilots.ps1'


def test_wolfram_and_validation_names():
    assert mod.rename_name('scripts/ExportGeneralT3Stage31OneLoop.wl') == 'scripts/ExportGeneralT3AllModelsUVOneLoop.wl'
    assert mod.rename_name('scripts/Stage22RGBetaNativeAccess.wl') == 'scripts/RGBetaNativeAccess.wl'
    assert mod.rename_name('scripts/VerifyGeneralT3Stage24BCPostFix.py') == 'scripts/VerifyGeneralT3BCQuarticCovarianceCorrection.py'


def test_docs_and_tests_get_named_by_physics():
    assert mod.rename_name('docs/GENERAL_T3_STAGE17_STATUS.md') == 'docs/GENERAL_T3_EFT1_SCALAR_PAIR_MIXING.md'
    assert mod.rename_name('docs/GENERAL_T3_STAGE22_110101_SU2_COVARIANCE.md') == 'docs/GENERAL_T3_TENSOR_110101_SU2_COVARIANCE.md'
    assert mod.rename_name('tests/python/test_general_t3_stage36_final_weinberg.py') == 'tests/python/test_general_t3_final_weinberg.py'


def test_other_artifacts_not_named_or_edited():
    assert mod.rename_name('output/general_t3_stage31/results.json') is None
    assert mod.rename_name('scripts/RunGeneralT3ScalarCampaign.py') is None
    assert not mod.eligible('output/studies/scalar_stage31.json')
    assert not mod.eligible('RGE/general/generated/scalar.json')
    assert not mod.eligible('docs/archive/GENERAL_T3_STAGE22_TENSOR_CERTIFICATION.md')


def test_references_replace_but_historical_numbers_are_preserved():
    pairs = mod.replacement_pairs({
        'scripts/Stage31CoupledUVRunning.py': 'scripts/CoupledUVRunning.py',
        'docs/GENERAL_T3_STAGE31_16_MODEL_UV_RUNNING.md': 'docs/GENERAL_T3_16_MODEL_UV_RUNNING.md',
    })
    text = (b'from scripts import Stage31CoupledUVRunning\n'
            b'python -m scripts.Stage31CoupledUVRunning\n'
            b'docs/GENERAL_T3_STAGE31_16_MODEL_UV_RUNNING.md\n'
            b'output/general_t3_stage31/run.json\n'
            b'Stage 31: historical experiment\n')
    new = mod.replace_text(text, pairs, 'README.md').decode()
    assert 'from scripts import CoupledUVRunning' in new
    assert 'python -m scripts.CoupledUVRunning' in new
    assert 'docs/GENERAL_T3_16_MODEL_UV_RUNNING.md' in new
    assert 'output/general_t3_stage31/run.json' in new
    assert 'Stage 31: historical experiment' in new


def test_document_top_heading_descriptive():
    source = b'# Stage 17: additional scalar-pair mixing\n\nStage 16 input\n'
    new = mod.replace_text(source, [], 'docs/GENERAL_T3_STAGE17_STATUS.md').decode()
    assert new.startswith('# General T3: additional scalar-pair mixing')
    assert 'Stage 16 input' in new


def test_local_git_plan_and_external_backup(tmp_path, monkeypatch):
    root = tmp_path/'repo'
    root.mkdir()
    subprocess.run(['git','init','-q'],cwd=root,check=True)
    subprocess.run(['git','config','user.name','Local Test'],cwd=root,check=True)
    subprocess.run(['git','config','user.email','test@localhost'],cwd=root,check=True)
    files = {
        'scripts/Stage31CoupledUVRunning.py': b'DEF = 3\n',
        'scripts/RunGeneralT3Stage31.ps1': b'python -m scripts.Stage31CoupledUVRunning\n',
        'tests/python/test_general_t3_stage31_coupled_uv_running.py': b'from scripts import Stage31CoupledUVRunning\n',
        'docs/GENERAL_T3_STAGE31_16_MODEL_UV_RUNNING.md': b'# Stage 31: running models\n',
        'README.md': b'See docs/GENERAL_T3_STAGE31_16_MODEL_UV_RUNNING.md\n',
        'output/general_t3_stage31/data.json': b'{"stage31": "immutable"}\n',
    }
    for path, data in files.items():
        p = root/path
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)
    subprocess.run(['git','add','.'],cwd=root,check=True)
    subprocess.run(['git','commit','-qm','fixture'],cwd=root,check=True)
    plan = mod.make_plan(root)
    assert len(plan['rename']) == 4
    assert len(plan['updates']) == 5
    assert 'output/general_t3_stage31/data.json' not in plan['updates']
    monkeypatch.setattr(mod, 'EXPECTED_HEAD', mod.head(root))
    backup = tmp_path/'external'
    manifest = mod.apply_plan(root, plan, backup)
    assert len(manifest['renames']) == 4
    assert (root/'scripts/CoupledUVRunning.py').read_bytes() == b'DEF = 3\n'
    assert not (root/'scripts/Stage31CoupledUVRunning.py').exists()
    assert (root/'README.md').read_text().strip().endswith('GENERAL_T3_16_MODEL_UV_RUNNING.md')
    assert (root/'output/general_t3_stage31/data.json').read_bytes() == files['output/general_t3_stage31/data.json']
    assert (backup/'originals/README.md').read_bytes() == files['README.md']
    assert (backup/'originals/scripts/Stage31CoupledUVRunning.py').read_bytes() == files['scripts/Stage31CoupledUVRunning.py']


def test_distinct_native_quartic_action_and_action_identity_test_names():
    assert mod.rename_name(
        'tests/python/test_general_t3_stage25_quartic_action.py'
    ) == 'tests/python/test_general_t3_native_quartic_action.py'
    assert mod.rename_name(
        'tests/python/test_general_t3_stage26_quartic_action.py'
    ) == 'tests/python/test_general_t3_quartic_action_identity.py'


def test_real_collision_regression_and_reference_rewriting(tmp_path):
    root = tmp_path / 'repo'
    root.mkdir()
    subprocess.run(['git', 'init', '-q'], cwd=root, check=True)
    subprocess.run(['git', 'config', 'user.name', 'Local Test'], cwd=root, check=True)
    subprocess.run(['git', 'config', 'user.email', 'test@localhost'], cwd=root, check=True)

    s25 = 'tests/python/test_general_t3_stage25_quartic_action.py'
    s26 = 'tests/python/test_general_t3_stage26_quartic_action.py'
    sources = {
        s25: b'def test_native_action(): pass\n',
        s26: b'def test_action_identity(): pass\n',
        'README.md': (f'python -m pytest {s25} {s26}\n').encode(),
    }
    for rel, data in sources.items():
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    subprocess.run(['git', 'add', '.'], cwd=root, check=True)
    subprocess.run(['git', 'commit', '-qm', 'fixture'], cwd=root, check=True)
    plan = mod.make_plan(root)
    assert len(plan['rename']) == 2
    assert len(set(plan['rename'].values())) == 2
    assert plan['rename'][s25] == 'tests/python/test_general_t3_native_quartic_action.py'
    assert plan['rename'][s26] == 'tests/python/test_general_t3_quartic_action_identity.py'
    new_readme = plan['updates']['README.md'][1].decode()
    assert 'test_general_t3_native_quartic_action.py' in new_readme
    assert 'test_general_t3_quartic_action_identity.py' in new_readme
    assert 'test_general_t3_stage25_quartic_action.py' not in new_readme
    assert 'test_general_t3_stage26_quartic_action.py' not in new_readme
