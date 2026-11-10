"""Isolated fixture-driven validation of final all-16 reproducer.

No Matchete installation or Windows benchmark outputs required. The fixtures
match the exact saved Batch-7 keys and 14-input source manifest contract.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.wilson.final_repro_batch8 import (
    B7, B7_SCHEMA, B7_STATUS, MODELS, SOURCE_CODE,
    SCENARIOS, FACTORS, ReproError, c5, collect, digest,
    source_paths, tag, verify_report, write_report,
)


def jwrite(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, sort_keys=True, allow_nan=False), encoding='utf-8')


def zwrite(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **data)


def diag(a):
    return np.diag(np.asarray(a, dtype=complex))


def build_saved_fixture(root: Path):
    root.mkdir(parents=True, exist_ok=True)
    for rel in SOURCE_CODE:
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('# synthetic validation fixture; NOT actual project implementation\n', encoding='utf-8')
    summary = {
        'schema': B7_SCHEMA, 'completed_count': 16, 'blocked_count': 0,
        'mass_scan_run': True, 'vev_gev': 246.22, 'mu_low_gev': 173.0,
        'results': {},
    }
    for idx, model in enumerate(MODELS):
        slug = tag(model)
        fac = 1 + idx*0.025
        hard = diag([2e-12, 4e-12, 7e-12])*fac
        running = hard*0.05
        bubble = hard*0.1
        reference = {
            'analytic_triangle': hard+running+bubble,
            'native_formal': hard*0.90+running+bubble,
            'native_as_UV_minus_soft_hypothesis': hard*1.12+running+bubble,
        }
        allpaths = source_paths(root, model)
        for name, path in allpaths.items():
            if name.endswith('_json'):
                jwrite(path, {'model': model, 'fixture': True})
            elif name == 'B1_hard':
                zwrite(path, {'C5_hard_F_triangle_MF': hard})
            elif name == 'B2_strict':
                zwrite(path, {'C5_EFT1_strict_running_delta': running})
            elif name == 'B5_MS_C5':
                zwrite(path, {'C5_scalar_bubble_MS': bubble, 'C5_SMEFT_MS': reference['native_as_UV_minus_soft_hypothesis']})
            elif name == 'B2_MS_C5':
                zwrite(path, {'C5_SMEFT_MS': reference['analytic_triangle']})
            elif name == 'B4_MS_C5':
                zwrite(path, {'C5_SMEFT_MS': reference['native_formal']})
            else:
                zwrite(path, {'fixture': np.asarray([idx])})
        source_hash = {k: digest(p) for k,p in allpaths.items()}
        folder = root/B7/slug
        folder.mkdir(parents=True, exist_ok=True)
        ct = hard*0.01
        native = ct*1.5
        zwrite(folder/'source_poles_and_counterterms.npz', {
            'C5_EFT1_scalar_total': ct,
            'C5_B0_determined_scalar_CT': ct * (1 + 2e-15),
            'C5_native_raw_pole_total': native,
        })
        scanned = {}
        byfactor = {}
        for factor in FACTORS:
            group = {}
            for scenario in SCENARIOS:
                # Fixed MF couplings, illustrative dependence on subtraction scale.
                upper = reference[scenario] * (1 + 0.02*np.log(factor))
                lower = upper * 0.95
                masses = np.sort(np.linalg.svd(-0.5*246.22**2*lower,
                                                 compute_uv=False))*1e9
                key = f'{scenario}_scale_{factor:g}'
                scanned[key+'_C5_SMEFT_MS'] = upper
                scanned[key+'_C5_at_low'] = lower
                group[scenario] = {
                    'scenario': scenario, 'subtraction_scale_factor': factor,
                    'mass_results_conditional_not_physical_MSbar': True,
                    'fixed_MF_source_couplings': True,
                    'masses_eV': masses.tolist(),
                    'dm21_sq_eV2': float(masses[1]**2-masses[0]**2),
                    'dm31_sq_eV2': float(masses[2]**2-masses[0]**2),
                    'mixing_defined': True,
                    'sin2_theta12': 0.31, 'sin2_theta13': 0.022, 'sin2_theta23': 0.55,
                    'C5_SMEFT_MS_norm_GeV_inverse': float(np.linalg.norm(upper)),
                }
            byfactor[str(factor)] = group
        zwrite(folder/'matching_scale_scenario_matrices.npz', scanned)
        (folder/'conditional_mass_and_mixing_scan.csv').write_text('header\nfixture\n', encoding='utf-8')
        ledger = {
            'schema': B7_SCHEMA, 'model': model, 'status': B7_STATUS,
            'source_SHA256': source_hash,
            'EFT_order': 'UV -> ALL F out -> EFT1 BOTH S -> BOTH S out -> SMEFT',
            'integrated_out_fermions_at_MF': 3, 'fermion_threshold_count': 1,
            'integrated_out_scalars_at_MS': 2, 'scalar_threshold_count': 1,
            'MF_GeV': 100_000., 'MS_GeV': 1000.,
            'mu_low_GeV': 173., 'vev_GeV': 246.22,
            'matching_scale_and_low_energy': {'by_factor': byfactor},
            'data_files_SHA256': {k: digest(folder/k) for k in (
                'source_poles_and_counterterms.npz',
                'matching_scale_scenario_matrices.npz',
                'conditional_mass_and_mixing_scan.csv')},
        }
        file = folder/'source_pole_and_matching_ledger.json'
        jwrite(file, ledger)
        original = byfactor['1.0']['analytic_triangle']['masses_eV']
        summary['results'][model] = {
            'model': model, 'status': B7_STATUS,
            'low_energy_running_computed': True,
            'native_vs_explicit_scalar_CT_relative': float(np.linalg.norm(native-ct)/np.linalg.norm(ct)),
            'EFT1_scalar_CT_vs_B0_relative': float(np.linalg.norm(ct-ct*(1+2e-15))/np.linalg.norm(ct*(1+2e-15))),
            'native_vs_candidate_seven_channel_ADM_relative': 0.4,
            'native_PL_kernel_count': 1 if model != 'B_-1' else 34,
            'EFT1_source_quartic_count': 5,
            'nominal_analytic_masses_eV': original,
            'output_sha256': digest(file),
        }
    jwrite(root/B7/'summary.json', summary)
    return root


@pytest.fixture
def fixture_root(tmp_path):
    return build_saved_fixture(tmp_path/'repo')


def test_all16_collect_and_write_verify_repeatable(fixture_root):
    rows, overview, matrices, details = collect(fixture_root)
    assert len(rows) == 144 and len(overview) == 16
    assert len(matrices) == 16*(3*3*2+5)
    assert len({r['model'] for r in overview}) == 16
    assert all(r['prescription'] == 'CONDITIONAL_NOT_PHYSICAL_MSbar' for r in rows)
    manifest = write_report(fixture_root)
    assert manifest['completed_models'] == 16
    assert manifest['nominal_scenario_rows'] == 48
    assert len(manifest['saved_inputs_and_batch7_outputs']['provenance']) == 16
    assert verify_report(fixture_root)['schema'] == manifest['schema']
    again = write_report(fixture_root)
    assert again['report_file_SHA256'] == manifest['report_file_SHA256']
    destination = fixture_root/'output/parallel_t3/final_repro_all16_batch8_v1'
    assert (destination/'all16_scenarios.csv').read_text().count('\n') == 145
    assert (destination/'all16_nominal.csv').read_text().count('\n') == 49


def test_modified_source_hash_refuses_report(fixture_root):
    p = source_paths(fixture_root, 'B_-1')['B1_hard_json']
    p.write_text('{"tampered": true}', encoding='utf-8')
    with pytest.raises(ReproError, match='source SHA256 mismatch'):
        write_report(fixture_root)
    assert not (fixture_root/'output/parallel_t3/final_repro_all16_batch8_v1').exists()


def test_missing_model_refuses_report(fixture_root):
    summary = fixture_root/B7/'summary.json'
    data = json.loads(summary.read_text())
    del data['results']['C_-1']
    jwrite(summary, data)
    with pytest.raises(ReproError, match='missing, duplicated or unexpected models'):
        collect(fixture_root)


def test_wrong_mass_from_low_c5_refuses_report(fixture_root):
    target = fixture_root/B7/tag('A_-4')/'source_pole_and_matching_ledger.json'
    ledger = json.loads(target.read_text())
    ledger['matching_scale_and_low_energy']['by_factor']['1.0']['analytic_triangle']['masses_eV'][0] *= 1.2
    jwrite(target, ledger)
    summary = fixture_root/B7/'summary.json'
    s = json.loads(summary.read_text())
    s['results']['A_-4']['output_sha256'] = digest(target)
    jwrite(summary, s)
    with pytest.raises(ReproError, match='masses vs low C5'):
        collect(fixture_root)


def test_missing_subtraction_scale_refuses_report(fixture_root):
    target = fixture_root/B7/tag('C_-1')/'source_pole_and_matching_ledger.json'
    ledger = json.loads(target.read_text())
    del ledger['matching_scale_and_low_energy']['by_factor']['2.0']
    jwrite(target, ledger)
    sfile = fixture_root/B7/'summary.json'
    s = json.loads(sfile.read_text())
    s['results']['C_-1']['output_sha256'] = digest(target)
    jwrite(sfile, s)
    with pytest.raises(ReproError, match='incomplete 3-scale output'):
        collect(fixture_root)


def test_report_hash_detects_tampering(fixture_root):
    write_report(fixture_root)
    path = fixture_root/'output/parallel_t3/final_repro_all16_batch8_v1/all16_overview.csv'
    path.write_text('tampered\n', encoding='utf-8')
    with pytest.raises(ReproError, match='Final report content changed'):
        verify_report(fixture_root)


def test_matrix_rejects_antisymmetric_total():
    with pytest.raises(ReproError, match='not flavour symmetric'):
        c5(np.array([[1,2,0],[-2,1,0],[0,0,1]],complex), 'bad')


def test_refuse_to_write_within_existing_source_batch(fixture_root):
    with pytest.raises(ReproError, match='must not overwrite'):
        write_report(fixture_root, B7/'my_report')


def test_one_command_collect_only_and_verify_only(fixture_root):
    import subprocess
    script = ROOT/'scripts/RunT3FinalRepro.py'
    p = subprocess.run([sys.executable, str(script), '--root', str(fixture_root),
                        '--collect-only', '--skip-tests'], capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    assert '16/16 models' in p.stdout
    v = subprocess.run([sys.executable, str(script), '--root', str(fixture_root),
                        '--verify-only'], capture_output=True, text=True)
    assert v.returncode == 0, v.stdout + v.stderr
    assert 'VERIFIED: 16/16 models' in v.stdout


def test_verify_detects_changed_code_without_new_report(fixture_root):
    write_report(fixture_root)
    source = fixture_root/SOURCE_CODE[0]
    source.write_text('# updated batch 7 runner\n', encoding='utf-8')
    with pytest.raises(ReproError, match='Installed physics/report code differs'):
        verify_report(fixture_root)


def test_verify_detects_new_batch7_summary_even_if_source_files_unchanged(fixture_root):
    write_report(fixture_root)
    old = fixture_root/B7/'summary.json'
    summary = json.loads(old.read_text())
    summary['executed_in_this_call'] = ['A_-4']
    jwrite(old, summary)
    with pytest.raises(ReproError, match='Batch-7 source or results changed'):
        verify_report(fixture_root)
