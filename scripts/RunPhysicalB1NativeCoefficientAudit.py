#!/usr/bin/env python3
"""B1 independent structural action coefficient/chiral audit (no production edits).

Reads three NATIVE Matchete action exports and the PREVIOUS completed B1 all-16
conditional audit. Does not regenerate frozen tensors. Writes a new evidence
folder; refuses overwriting. Optional abstract GammaCC witness can be supplied.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Numerical.lo_certification.physical_normalisation.b1_native_action_contract import (  # noqa E402
    ActionContractError, require, verify_record,
)
from Numerical.lo_certification.physical_normalisation.b1_matchete_weyl import (  # noqa E402
    four_spinor_witness, spurion_charge, model_specs,
)
from scripts.RunPhysicalB1WeylAudit import guard_source, git_revision  # noqa E402


MODELS = ('A_p0', 'B_m1', 'C_m1')
AUDIT_STATUS = 'All16ConditionalCandidateWeylTensorsReconstructed_PhysicalDictionaryPending'
STATUS = 'ThreeNativeActionsCoefficientAndChiralityVerified_AbsoluteCCAndB2Pending'
CC_STATUS = 'NativeAbstractCCAlgebraVerified_NumericOrientationPending'


def load(path: Path) -> dict:
    require(path.is_file(), 'Missing evidence: '+str(path))
    value = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(value, dict), 'Expected JSON object: '+str(path))
    return value


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_all16(value: dict) -> dict:
    require(value.get('status') == AUDIT_STATUS and
            value.get('count_conditional_comparison_passed') == 16 and
            value.get('physical_source_certified') is False,
            'Missing successful 16/16 CONDITIONAL comparison; source certification must remain false')
    require(value.get('native_action_observed') is True,
            'Prior all-16 run did not observe three native actions')
    models = value.get('models')
    require(isinstance(models, list) and len(models) == 16, 'Expected 16 model records')
    expected = {m for m, _ in model_specs()}
    require({m.get('model') for m in models} == expected,
            'All-16 set has missing or duplicate models')
    results = {}
    for item in models:
        require(item.get('physically_verified') is False, 'False source certificate in all16 model')
        data = item.get('residuals')
        require(isinstance(data, dict), 'Missing model relative tensor residuals')
        required = ('mass', 'Y_after_global_minus', 'Z_after_global_minus')
        import math
        for key in required:
            val = data.get(key)
            require(type(val) in (int, float) and math.isfinite(val) and 0 <= val < 2e-9,
                    f"{item['model']}: invalid residual {key}={val!r}")
        results[item['model']] = {k: data[k] for k in required}
    return {'conditional_models': 16,
            'max_mass_relative_residual': max(x['mass'] for x in results.values()),
            'max_Y_relative_residual': max(x['Y_after_global_minus'] for x in results.values()),
            'max_Z_relative_residual': max(x['Z_after_global_minus'] for x in results.values()),
            'models': results}


def check_cc(value: dict) -> dict:
    require(value.get('format') == 't3-b1-native-cc-abstract-v2' and
            value.get('status') == CC_STATUS and value.get('model') == 'A_p0',
            'Unexpected native CC witness schema/status')
    tests = value.get('checks')
    require(isinstance(tests, dict) and len(tests) == 6 and all(x is True for x in tests.values()),
            'One or more installed abstract CC identities failed')
    for key in ('chosen_numeric_gamma_matrix_test_performed_here',
                'installed_numeric_gamma_cc_sign_measured',
                'cross_package_spinor_phase_certified', 'physical_source_certified'):
        require(value.get(key) is False, 'Native abstract witness cannot claim '+key)
    manifest = value.get('installed_symbol_manifest', {})
    require(isinstance(manifest, dict) and str(manifest.get('symbol_full_name', '')).endswith('GammaCC'),
            'Installed GammaCC metadata missing')
    return {'abstract_native_cc_checks': 6,
            'installed_numeric_gamma_cc_sign_measured': False,
            'installed_symbol_manifest': manifest}


def build(root: Path, native: Path, prior: Path, cc: Path | None) -> dict:
    inventory = guard_source(root)
    actions = {}
    for model in MODELS:
        src = native/(model+'.json')
        actions[model] = {'sha256': sha(src), **verify_record(load(src), model)}
    prior_record = load(prior)
    conditional = check_all16(prior_record)
    # Bind the three structurally parsed native files to the EXACT files that
    # the prior successful 16-model reconstruction audited.
    previous_native = prior_record.get('native', {}).get('representatives', {})
    require(isinstance(previous_native, dict) and set(previous_native) == set(MODELS),
            'Prior 16/16 audit has no matching three-native provenance')
    for model in MODELS:
        require(previous_native[model].get('sha256') == actions[model]['sha256'],
                model+': native JSON changed since the prior 16/16 audit')
    installed_cc = check_cc(load(cc)) if cc else {'status': 'NOT_EXECUTED',
        'installed_numeric_gamma_cc_sign_measured': False}
    four = four_spinor_witness()
    require(four.get('exact_chosen_basis_checked') is True and
            four.get('absolute_native_GammaCC_phase_measured') is False,
            'Chosen-basis mathematical witness wrongly promoted')
    # This field is intentionally a physical, observable-invariant statement,
    # independent of any unmeasured overall heavy-fermion spinor phase.
    return {
        'format': 't3-b1-native-source-factor-audit-v2',
        'status': STATUS,
        'execution_utc': datetime.now(timezone.utc).isoformat(),
        'local_git_revision': git_revision(root),
        'guarded_source_blobs': inventory,
        'native_action_witnesses': actions,
        'native_actions_structurally_verified_count': 3,
        'mass_coefficient_contract': {'Dirac_A_p0': '-1', 'Majorana_B_m1': '-1/2',
                                      'Majorana_C_m1': '-1/2'},
        'yukawa_original_and_hc_coefficient_contract': 'positive 1 (each), in all three actions',
        'previous_conditional_16': conditional,
        'previous_conditional_audit_sha256': sha(prior),
        'native_cc_abstract_probe': installed_cc,
        'chosen_basis_algebra': four,
        'spurion': {'Dirac': spurion_charge(False),'Majorana': spurion_charge(True)},
        'checks': {'native_mass_yukawa_and_hc_structure': True,
                   'restricted_flags_checked': True,
                   'conditional_all16_residuals_small': True,
                   'exact_chosen_basis_weyl_algebra': True},
        'limitations': [
            'Native symbolic GammaCC not anchored to a concrete 4x4 gamma matrix sign',
            'No physical Matchete-RGBeta SU2 CG component axes or relative phases',
            'Three representatives do not independently measure the remaining 13 native actions',
            'Frozen all16 inputs are conditional candidate, not external Feynman rule certification',
            'No physical one-loop amplitude/source normalisation certificate was issued'],
        'relative_physical_factors_verified': False,
        'physical_source_certified': False,
        'physical_source_certified_models': 0,
        'physical_B1_scope_complete': False,
        'safe_to_integrate_into_production': False,
        'EFT_or_running_performed': False}


def handoff(report: dict) -> str:
    a = report['previous_conditional_16']
    return ('# B1 next-batch: native source coefficient and chirality witness\n\n'
      f'- Local Git revision: `{report["local_git_revision"]}`. Source byte equivalence checked against pinned Git blobs.\n'
      '- Native symbolic actions A_+0, B_-1, C_-1: **3/3 structurally checked** (original and HC for both Yukawas, spinor/projector orientation, mass and symmetry prefactors).\n'
      '- Extracted mass coefficient: Dirac `-1`, Majorana `-1/2` **once**; Yukawa originals and HC `+1` for each.\n'
      f'- Prior conditional 16/16 maximum relative residuals: mass `{a["max_mass_relative_residual"]:.4g}`, '
      f'Y `{a["max_Y_relative_residual"]:.4g}`, Z `{a["max_Z_relative_residual"]:.4g}`.\n'
      '- A native GammaCC probe, if supplied, verifies symbolic Matchete identities only; an overall numerical sign relative to the declared Weyl representation remains **unmeasured**.\n'
      '- Relative CG phases/normalisations from B2 remain **unmeasured** in this B1 audit. No source-certified neutrino mass coefficient.\n'
      '- Common neutrino sign reversal of both Y and Z is a field rephasing; not a physical factor-of-two discrepancy.\n'
      '- Physical certifications: **0/16**, even after all checks pass.\n'
      '- Safe to integrate into production: **NO**. Native mass/Yukawa structural coefficients and chosen-basis identities may be used as source-grounded inputs to the later joint physical-dictionary step.\n'
      '- No production files, frozen inputs, GitHub, EFT, or RGE changed.\n')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root', type=Path, default=ROOT)
    ap.add_argument('--native-dir', type=Path, required=True)
    ap.add_argument('--conditional-audit', type=Path, required=True)
    ap.add_argument('--cc-abstract-json', type=Path)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    root = args.project_root.resolve()
    def resolve(p): return p if p.is_absolute() else root / p
    out = resolve(args.output_dir)
    if out.exists():
        print('BLOCKED: output folder already exists: '+str(out), file=sys.stderr)
        return 3
    try:
        r = build(root, resolve(args.native_dir), resolve(args.conditional_audit),
                  resolve(args.cc_abstract_json) if args.cc_abstract_json else None)
    except (ActionContractError, ValueError, KeyError, TypeError, OSError) as exc:
        print('BLOCKED (no output written): '+str(exc), file=sys.stderr)
        return 2
    out.mkdir(parents=True, exist_ok=False)
    (out/'native_action_coefficients.json').write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',
                                                      encoding='utf-8')
    (out/'HANDOFF.md').write_text(handoff(r),encoding='utf-8')
    print(STATUS+'; native=3/3; prior_conditional=16/16; physical_source_certified=FALSE')
    print('Evidence: '+str(out/'native_action_coefficients.json'))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
