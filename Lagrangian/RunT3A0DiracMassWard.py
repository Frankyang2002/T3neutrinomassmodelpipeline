#!/usr/bin/env python3
"""Workstream A: actual RGBeta A(0) Dirac mass SU(2) covariance gate.

Runs one native RGBeta registration probe on the *unchanged production
model functions*, and performs the SU(2) component Ward identity independently.
If F and Fc are both registered as SU2L[fund], their holomorphic left-Weyl
Dirac mass is invariant only for epsilon, not for a two-fundamental delta:

  (T^a)^T C + C T^a = 0,  a=1,2,3.

A successful AUDIT can therefore report a verified NONINVARIANT coupling;
it must never confuse successful AddFermionMass registration with a passed Ward
identity. Does not modify any production or RGE files.

Usage (repository root):
  python Lagrangian/RunT3A0DiracMassWard.py --self-test
  python Lagrangian/RunT3A0DiracMassWard.py --all

Writes Reports/output/workstream_a/a0_dirac_mass_ward/...
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import sys
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path('Lagrangian/RunT3A0DiracMassRGBetaNative.wl')
SOURCE = Path('RGE/running/rgbeta/T3RGBetaModel.wl')
OUTPUT = Path('Reports/output/workstream_a/a0_dirac_mass_ward')
NATIVE = OUTPUT / 'native_rgbeta/A_p0.json'
REPORT = OUTPUT / 'a0_dirac_mass_ward.json'
NATIVE_STATUS = 'NativeA0FandFcMassTensorCaptured'
STATUS_FAIL = 'NativeA0DiracMassSU2WardViolationConfirmed'
STATUS_PASS = 'NativeA0DiracMassSU2WardPassed'


class AuditError(RuntimeError):
    pass


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise AuditError(reason)


def matrices() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    return (sp.Matrix([[0, 1], [1, 0]]) / 2,
            sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2,
            sp.Matrix([[1, 0], [0, -1]]) / 2)


def residuals(c: sp.Matrix, conjugate_second: bool = False) -> list[sp.Matrix]:
    require(c.shape == (2, 2), 'Expected a 2x2 mass contraction tensor')
    return [sp.simplify(t.T * c + c * (-sp.conjugate(t) if conjugate_second else t))
            for t in matrices()]


def zero(m: sp.Matrix) -> bool:
    return all(sp.simplify(x) == 0 for x in m)


def source_guard(root: Path) -> str:
    p = root / SOURCE
    require(p.is_file(), 'Missing current RGBeta model source: ' + str(p))
    s = p.read_text(encoding='utf-8')
    compact = re.sub(r'\s+', '', re.sub(r'\(\*.*?\*\)', '', s, flags=re.S))
    # Guard the specific model-definition claims with the installed local file.
    for token in (
        'T3RGBetaAddHeavyFermion[dF_Integer,yF_]',
        'T3RGBetaGaugeRep[dF,yF]',
        'T3RGBetaGaugeRep[dF,-yF]',
        'AddFermionMass[MF,{F,Fc}',
        'del[T3RGBetaRep[dF],#1,#2]',
        'Chirality->Right',
    ):
        require(token in compact, 'RGBeta source changed; inspect '+str(p)+' (missing '+token+')')
    return sha256(p.read_bytes()).hexdigest()


def check_native(r: dict[str, Any]) -> str:
    require(r.get('format') == 't3-workstream-a-native-a0-mass-gauge-tensor-v1'
            and r.get('status') == NATIVE_STATUS and r.get('model') == 'A_p0',
            'Native RGBeta probe did not capture the actual A(0) mass tensor')
    for field, expected in (
        ('F_gauge_rep_inputform', '{U1Y[1/2], SU2L[fund]}'),
        ('Fc_gauge_rep_inputform', '{U1Y[-1/2], SU2L[fund]}'),
        ('MF_fields_inputform', '{F, Fc}'),
        ('MF_chirality', 'Right'),
        ('MF_internal_mass_tensor_prefactor', '2'),
    ):
        require(r.get(field) == expected,
                'Native A(0) mass registration changed: '+field+' = '+str(r.get(field)))
    require(r.get('both_fields_in_SU2_fundamental') is True
            and r.get('F_and_Fc_have_opposite_U1_hypercharge') is True
            and r.get('native_invariant_checked_by_rgbeta') is False
            and r.get('component_ward_evaluated_by_native_package') is False
            and r.get('mass_beta_evaluated') is False
            and r.get('production_files_modified') is False,
            'Unexpected native mass-probe metadata')
    inv = r.get('MF_group_invariant_kind')
    require(inv in ('delta', 'epsilon', 'other'), 'Native mass invariant kind missing')
    expected = {
        'delta': 'del[SU2L[fund], dm1, dm2]',
        'epsilon': 'eps[SU2L[fund], dm1, dm2]',
    }
    if inv in expected:
        require(r.get('MF_group_invariant_inputform') == expected[inv],
                'Native invariant expression does not match its stated kind')
    return inv


def analyze(r: dict[str, Any], source_hash: str) -> dict[str, Any]:
    invariant = check_native(r)
    require(invariant != 'other',
            'Native MF invariant not recognized: cannot map to a physical tensor')
    delta, eps = sp.eye(2), sp.Matrix([[0, 1], [-1, 0]])
    c = delta if invariant == 'delta' else eps
    ward = residuals(c)
    passed = all(zero(x) for x in ward)
    # Independent controls: only epsilon is invariant between two SU2 fundamentals.
    require(all(zero(x) for x in residuals(eps)),
            'Internal SU2 epsilon Ward reference is wrong')
    require(not all(zero(x) for x in residuals(delta)),
            'Internal SU2 delta Ward negative control unexpectedly passed')
    require(all(zero(x) for x in residuals(delta, conjugate_second=True)),
            'Fundamental versus antifundamental delta control unexpectedly failed')
    require(passed == (invariant == 'epsilon'),
            'Native tensor classification inconsistent with exact SU2 Ward identity')
    return {
        'format': 't3-workstream-a-a0-native-dirac-mass-ward-v1',
        'status': STATUS_PASS if passed else STATUS_FAIL,
        'model': 'A_p0',
        'registered_F_representation': r['F_gauge_rep_inputform'],
        'registered_Fc_representation': r['Fc_gauge_rep_inputform'],
        'registered_MF_invariant': r['MF_group_invariant_inputform'],
        'native_MF_internal_prefactor': '2',
        'local_rgbeta_source_sha256': source_hash,
        'a0_hypercharge_sum_zero': True,
        'native_add_fermion_mass_registered': True,
        'actual_su2_mass_ward_passed': passed,
        'nonzero_SU2_generator_residuals': [i + 1 for i, x in enumerate(ward) if not zero(x)],
        'su2_residual_matrices': [[[str(sp.simplify(v)) for v in row]
                                  for row in x.tolist()] for x in ward],
        'su2_epsilon_positive_control_passed': True,
        'su2_delta_negative_control_detected': True,
        'fundamental_antifundamental_delta_control_passed': True,
        'source_model_is_physically_mass_invariant': passed,
        'a0_yukawa_spinor_phase_mass_transport_certified': False,
        'physical_matchete_rgbeta_coupling_dictionary_certified': False,
        'production_files_modified': False,
        'next_gate': ('Mass tensor already SU2 covariant; proceed to one common spinor/mass phase map'
                      if passed else
                      'Before phase transport, correct the production SU2 mass invariant with a separate validated change; then rerun Ward and beta calculations'),
    }


def execute(root: Path) -> dict[str, Any]:
    source_hash = source_guard(root)
    out = root / OUTPUT
    (out / 'native_rgbeta').mkdir(parents=True, exist_ok=True)
    native = root / NATIVE
    if native.exists():
        native.unlink()  # prevent stale JSON being accepted
    cmd = ['wolframscript', '-file', str(root / SCRIPT), str(native)]
    try:
        proc = subprocess.run(cmd, cwd=root, capture_output=True, text=True,
                              errors='replace', check=False)
    except OSError as exc:
        raise AuditError('Could not launch wolframscript: ' + str(exc)) from exc
    log = out / 'native_rgbeta/A_p0.log'
    log.write_text(proc.stdout + '\nSTDERR:\n' + proc.stderr, encoding='utf-8')
    if proc.returncode != 0:
        details = ''
        if native.is_file():
            try:
                nr = json.loads(native.read_text(encoding='utf-8'))
                details = ': '+str(nr.get('failed_step', ''))+' '+str(nr.get('reason', ''))
            except (OSError, ValueError):
                pass
        raise AuditError(f'Native RGBeta mass registration failed (exit {proc.returncode})'
                         +details+'; see '+str(log))
    require(native.is_file(), 'Wolfram returned success but native JSON missing: '+str(native))
    result = analyze(json.loads(native.read_text(encoding='utf-8')), source_hash)
    (root / REPORT).write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


def self_test() -> None:
    raw = {
        'format': 't3-workstream-a-native-a0-mass-gauge-tensor-v1',
        'status': NATIVE_STATUS, 'model': 'A_p0',
        'F_gauge_rep_inputform': '{U1Y[1/2], SU2L[fund]}',
        'Fc_gauge_rep_inputform': '{U1Y[-1/2], SU2L[fund]}',
        'both_fields_in_SU2_fundamental': True,
        'F_and_Fc_have_opposite_U1_hypercharge': True,
        'MF_fields_inputform': '{F, Fc}', 'MF_chirality': 'Right',
        'MF_group_invariant_inputform': 'del[SU2L[fund], dm1, dm2]',
        'MF_group_invariant_kind': 'delta',
        'MF_internal_mass_tensor_prefactor': '2',
        'native_invariant_checked_by_rgbeta': False,
        'component_ward_evaluated_by_native_package': False,
        'mass_beta_evaluated': False,
        'production_files_modified': False,
    }
    test = analyze(raw, 'fixture')
    require(test['status'] == STATUS_FAIL and
            test['nonzero_SU2_generator_residuals'] == [1, 3],
            'Expected native delta mass violation not detected')
    good = dict(raw, MF_group_invariant_kind='epsilon',
                MF_group_invariant_inputform='eps[SU2L[fund], dm1, dm2]')
    require(analyze(good, 'fixture')['actual_su2_mass_ward_passed'] is True,
            'Epsilon positive control failed')
    for name, value in (
        ('Fc_gauge_rep_inputform','{U1Y[-1/2], Bar[SU2L[fund]]}'),
        ('both_fields_in_SU2_fundamental', False),
        ('MF_internal_mass_tensor_prefactor', '1'),
    ):
        bad = dict(raw, **{name: value})
        try:
            analyze(bad, 'fixture')
        except AuditError:
            pass
        else:
            raise AuditError('Did not reject corrupt native fixture: '+name)
    print('PASS: exact SU2 mass Ward; delta fails generators 1,3; epsilon passes; conjugate-rep delta passes; negative controls')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--repo-root', type=Path, default=ROOT)
    args = ap.parse_args()
    try:
        if args.self_test:
            self_test()
            return 0
        require(args.all, 'Specify --all or --self-test')
        r = execute(args.repo_root.resolve())
        print(r['status'])
        print('Native RGBeta MF tensor: '+r['registered_MF_invariant'])
        print('Independent SU(2) Ward residuals: '+str(r['nonzero_SU2_generator_residuals']))
        print('Mass phase transport / physical dictionary: NOT CERTIFIED')
        print('Report: '+str(args.repo_root.resolve()/REPORT))
        return 0 if r['actual_su2_mass_ward_passed'] else 4
    except (AuditError, ValueError, KeyError, OSError, TypeError) as exc:
        print('A(0) mass gauge covariance audit BLOCKED: '+str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
