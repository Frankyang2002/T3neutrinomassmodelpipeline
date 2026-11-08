"""Read-only A(-4) F-only Wilson master-basis and gauge-only beta audit.

This reconstructs the actual unscaled native 60-coordinate source tensor,
transports it through the source-anchored complex triplet matrix, independently
reconstructs the 36-coordinate master tensor, and checks the full 45-Weyl
SU(3)xSU(2)xU(1) Ward constraints. It compares a separate gauge-only tensor
contraction against the existing provisional six-term EFT1 implementation.

NO physical Matchete GammaCC/Weyl normalization is supplied or inferred.
NO matching, full RGE integration, or RGBeta correction is performed.

PowerShell (from repository root):
  python scripts/RunAminus4FullGaugeWardAudit.py --output output/a_minus4_next_batch/gauge_ward_audit.json
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import subprocess
from math import log, pi
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.running import Aminus4EFT1Couplings as eft
from RGE.running.intermediate.GeneralEFT1SevenTermFlow import Couplings, total_beta
from RGE.running.intermediate.GeneralEFT1SixTermReferenceCandidate import six_term_beta

NATIVE = Path('output/general_observable/lo_certification_spinor/source_e515d82d118284ba/A_m4/boundary.json')
MASTER = Path('output/general_observable/a_minus4_boundary_batch1/master_basis_candidate.json')
QUARTIC = Path('output/a_minus4_eft1_batch2/A_m4_projected_raw.json')
QUARTIC_AUDIT = Path('output/a_minus4_eft1_batch2/A_m4_projection_audit.json')

class AuditFailure(ValueError):
    pass


def check(condition, message):
    if not condition:
        raise AuditFailure(message)


def load(root, relative):
    path = root / relative
    check(path.is_file(), f'missing source {path}')
    raw = path.read_bytes()
    return json.loads(raw.decode('utf-8-sig')), sha256(raw).hexdigest()


def fill_orbit(c, i, j, a, b, value):
    for ii, jj in {(i, j), (j, i)}:
        for aa, bb in {(a, b), (b, a)}:
            old = c[ii, jj, aa, bb]
            check(abs(old) < 1e-25 or abs(old - value) < 1e-25,
                  f'inconsistent duplicate coefficient at {(ii, jj, aa, bb)}')
            c[ii, jj, aa, bb] = value


def native_tensor(data):
    check(data.get('operator') == 'LLSS', 'not the F-only LLSS source')
    check(data.get('index_order') == ['lepton_su2_i', 'lepton_flavour_p',
                                    'lepton_su2_j', 'lepton_flavour_q',
                                    'real_scalar_a', 'real_scalar_b'],
          'source coordinate order changed')
    rows = data.get('summed_unscaled')
    check(isinstance(rows, list) and len(rows) == 60, 'native source inventory changed')
    c = np.zeros((45, 45, 12, 12), dtype=complex)
    for row in rows:
        i, p, j, q, a, b = row['coordinate']
        check(all(1 <= x <= top for x, top in
                  ((i, 2), (j, 2), (p, 3), (q, 3), (a, 12), (b, 12))),
              'invalid native coordinate')
        check(a in (5, 6) and b in range(7, 13), 'unexpected source scalar support')
        ii = eft.FI[f'L[{p},{i}]']
        jj = eft.FI[f'L[{q},{j}]']
        val = complex(row['coefficient_real'], row['coefficient_imag'])
        orbit_size = len({(x, y, u, v) for x, y in ((ii, jj), (jj, ii))
                          for u, v in ((a-1, b-1), (b-1, a-1))})
        check(row['orbit_size'] == orbit_size, 'native full-sum orbit size changed')
        fill_orbit(c, ii, jj, a-1, b-1, val)
    return c


def master_tensor(data):
    rows = data.get('conditional_master_weight_basis_independent_entries')
    check(isinstance(rows, list) and len(rows) == 36, 'master-coordinate inventory changed')
    c = np.zeros((45, 45, 12, 12), dtype=complex)
    leptons = [eft.FI[f'L[{p},{w}]'] for p in range(1, 4) for w in (1, 2)]
    for row in rows:
        i, j, a, b = (row[k] for k in ('i', 'j', 'a', 'b'))
        check(1 <= i <= j <= 6 and 5 <= a <= 6 and 7 <= b <= 12,
              'unexpected master coordinate support')
        ii, jj = leptons[i-1], leptons[j-1]
        val = complex(row['coefficient_real_if_rho_1'],
                      row['coefficient_imag_if_rho_1'])
        orbit_size = len({(x, y, u, v) for x, y in ((ii, jj), (jj, ii))
                          for u, v in ((a-1, b-1), (b-1, a-1))})
        check(row['distinct_orbit_size'] == orbit_size,
              'master full-sum orbit size changed')
        fill_orbit(c, ii, jj, a-1, b-1, val)
    return c


def realification(matrix):
    u = np.asarray(matrix, complex)
    check(u.ndim == 2 and u.shape[0] == u.shape[1], 'non-square complex scalar map')
    size = len(u)
    ret = np.zeros((2*size, 2*size), float)
    for i in range(size):
        for j in range(size):
            z = u[i, j]
            ret[2*i:2*i+2, 2*j:2*j+2] = [[z.real, -z.imag],
                                          [z.imag, z.real]]
    return ret


def transform_native_to_master(c, r):
    # NATIVE S2 = N^dagger MASTER S2 => R(NATIVE) = R(N)^T R(MASTER).
    check(r.shape == (6, 6), 'triplet realification is not 6x6')
    result = np.zeros_like(c)
    result[:, :, :, 6:12] = np.einsum('ijac,cb->ijab',
                                     c[:, :, :, 6:12], r.T, optimize=True)
    result[:, :, 6:12, :] = result[:, :, :, 6:12].swapaxes(2, 3)
    return result


def ward(c, tf, theta):
    return (np.einsum('ki,kjab->ijab', tf, c, optimize=True)
            + np.einsum('kj,ikab->ijab', tf, c, optimize=True)
            + np.einsum('ca,ijcb->ijab', theta, c, optimize=True)
            + np.einsum('cb,ijac->ijab', theta, c, optimize=True))


def relative_max(a, b):
    return float(np.max(np.abs(a-b)) / max(np.max(np.abs(b)), 1e-300))


def ward_checks(c, gauge):
    scale = max(float(np.max(np.abs(c))), 1e-300)
    return {name: [float(np.max(np.abs(ward(c, tf, theta))) / scale)
                   for tf, theta in zip(gauge[name]['fermion'], gauge[name]['scalar'])]
            for name in ('U1', 'SU2', 'SU3')}


def gauge_only_reference(c, tf_list, theta_list):
    """Independent gauge-only transcription; g=1, Y=lambda=0.

    The four gauge contributions are evaluated directly from the sums in
    arXiv:2512.16890v1 Eq. (4.85), but with the fermion Casimir transposed
    for coefficient-index covariance. This is not a second loop derivation.
    """
    pieces = {k: np.zeros_like(c) for k in
              ('scalar_pair', 'mixed_yukawa_gauge',
               'scalar_anomalous_dimension', 'fermion_anomalous_dimension')}
    for tf, theta in zip(tf_list, theta_list):
        raw = -2 * np.einsum('ac,bd,ijcd->ijab', theta, theta, c, optimize=True)
        pieces['scalar_pair'] += (raw + raw.swapaxes(2, 3)) / 2
        raw = 4 * np.einsum('lj,bc,ilac->ijab', tf, theta, c, optimize=True)
        pieces['mixed_yukawa_gauge'] += (raw + raw.swapaxes(0, 1)
                                          + raw.swapaxes(2, 3)
                                          + raw.transpose(1, 0, 3, 2))
        raw = np.einsum('bd,ijad->ijab', -4 * theta @ theta, c, optimize=True)
        pieces['scalar_anomalous_dimension'] += raw + raw.swapaxes(2, 3)
        raw = np.einsum('jk,ikab->ijab', -3 * (tf @ tf).T, c, optimize=True)
        pieces['fermion_anomalous_dimension'] += raw + raw.swapaxes(0, 1)
    return pieces


def gauge_sector_report(c, ff, ss, leptons, sector):
    tfs = [tf[np.ix_(leptons, leptons)] for tf in ff]
    theta = list(ss)
    independent = gauge_only_reference(c, tfs, theta)
    expected = sum(independent.values(), np.zeros_like(c))
    couplings = Couplings(np.zeros((6, 6, 12), complex),
                          np.zeros((12,)*4, float),
                          ((1., tuple(zip(tfs, theta))),))
    production = six_term_beta(c, couplings)
    actual = total_beta(production)
    if sector != 'SU3':
        check(float(np.max(np.abs(expected))) > 0, 'expected gauge beta vanished')
    if float(np.max(np.abs(expected))) > 0:
        discrepancy = relative_max(actual, expected)
        check(discrepancy < 5e-11, f'six-term gauge disagreement in {sector}')
    else:
        discrepancy = float(np.max(np.abs(actual)))
        check(discrepancy < 1e-20, 'gauge-neutral SU3 Wilson term unexpectedly runs')
    for name, values in independent.items():
        check(relative_max(production[name], values) < 5e-11 if np.max(abs(values))
              else np.max(abs(production[name])) < 1e-20,
              f'gauge channel {name} differs in {sector}')
    check(np.max(abs(production['crossed_yukawa'])) < 1e-20,
          'six-term crossed-Yukawa term unexpectedly nonzero')
    norm_c = np.vdot(c, c)
    eigenvalue = np.vdot(c, expected)/norm_c
    eigen_res = float(np.linalg.norm(expected-eigenvalue*c)/
                      max(np.linalg.norm(expected), 1e-300))
    check(eigen_res < 5e-11, f'A(-4) Wilson not an eigenvector of {sector} gauge beta')
    w = max(float(np.max(np.abs(ward(expected,tf,ts))) /
                 max(np.max(np.abs(expected)),1e-300))
            for tf,ts in zip(tfs,theta))
    check(w < 5e-11, f'gauge beta fails {sector} Ward closure')
    return {'gauge_beta_over_g_squared_eigenvalue': {'real':float(eigenvalue.real),
                                                     'imag':float(eigenvalue.imag)},
            'eigenvector_relative_residual': eigen_res,
            'implementation_vs_independent_contraction_relative': discrepancy,
            'beta_Ward_max_relative': w,
            'channel_coefficients': {k: float((np.vdot(c,v)/norm_c).real)
                                     for k,v in independent.items()}}


def local_git_revision(root: Path):
    try:
        result = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=root,
                                capture_output=True, text=True, timeout=5,
                                check=True)
        return result.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None


def report(root: Path):
    native, native_sha = load(root,NATIVE)
    master, master_sha = load(root,MASTER)
    check(native.get('model') == master.get('model') == 'A_-4', 'different models')
    check(native.get('normalization',{}).get('certified') is False,
          'native normalization state changed; review before auditing')
    check(master.get('physical_C_ijab_at_MF') is None and
          master.get('physical_certification') is False,
          'master physical boundary status changed; review before auditing')
    cn = native_tensor(native)
    cm = master_tensor(master)
    mat = master['native_to_master_triplet_complex_matrix']
    n = np.array(mat['real']) + 1j*np.array(mat['imag'])
    check(n.shape == (3,3), 'invalid triplet intertwiner')
    unitarity = float(np.max(abs(n.conj().T@n-np.eye(3))))
    check(unitarity < 5e-12, 'triplet map not unitary')
    r = realification(n)
    check(relative_max(r, np.array(master['native_to_master_triplet_real_matrix'])) < 5e-12,
          'source complex and real triplet maps inconsistent')
    transported = transform_native_to_master(cn, r)
    source_transport = relative_max(transported, cm)
    check(source_transport < 5e-11, 'source native-to-master tensor transport failed')
    norm_err = abs(np.linalg.norm(cn)-np.linalg.norm(cm))/np.linalg.norm(cm)
    check(norm_err < 5e-11, 'source/master tensor norm not preserved')
    gauge=eft.gauge_generators()
    naive = ward_checks(cn, gauge)
    converted = ward_checks(cm, gauge)
    check(max(v for values in converted.values() for v in values) < 5e-11,
          'converted 45-Weyl A(-4) Wilson fails Ward closure')
    check(max(naive['SU2']) > 1e-2,
          'raw native triplet tensor unexpectedly already master SU2-invariant')
    leptons = [eft.FI[f'L[{p},{w}]'] for p in (1,2,3) for w in (1,2)]
    c6 = cm[np.ix_(leptons, leptons, range(12), range(12))]
    channel = {name:gauge_sector_report(c6,gauge[name]['fermion'],
                                         gauge[name]['scalar'],leptons,name)
               for name in ('U1','SU2','SU3')}
    expected = {'U1':-13.5, 'SU2':-4.5, 'SU3':0.0}
    for name, value in expected.items():
        check(abs(channel[name]['gauge_beta_over_g_squared_eigenvalue']['real']-value)<5e-11,
              f'A(-4) gauge-only eigenvalue changed for {name}')
    # An illustrative fixed-g leading log; NOT a physical full EFT1 trajectory.
    gY, g2 = .36, .65
    logratio = log(1e3/1e5)
    gauge_log = (expected['U1']*gY*gY + expected['SU2']*g2*g2)*logratio/(16*pi*pi)

    # Re-evaluate native RGBeta polynomial against reference *now*, not just
    # trust the saved audit's status. The relevant functions already exist in repo.
    from Numerical.running.Aminus4EFT1Audit import scalar_beta_comparison
    from scripts.VerifyAminus4GlobalQuarticProjection import (
        unit_point, quadratic_probes, compatibility_payload)
    beta_raw,_ = load(root,QUARTIC)
    beta = compatibility_payload(beta_raw)
    failures = []
    for probe in quadratic_probes():
        comp = scalar_beta_comparison(unit_point(probe), beta)
        if not comp['reference_and_native_scalar_loop_match']:
            failures.append({'unit_directions':list(probe),
                             'mismatched_components':comp['mismatched_coupling_betas'],
                             'details':{k:comp['coupling_beta_comparison'][k]
                                        for k in comp['mismatched_coupling_betas']}})
    saved,_ = load(root,QUARTIC_AUDIT)
    check(saved['scalar']['scalar_quadratic_probe_count'] == len(list(quadratic_probes())),
          'saved and executable quartic probe counts differ')
    check(len(saved['scalar']['failures']) == len(failures),
          'current quartic source results disagree with saved source audit')
    quartic_summary = {'independently_rerun_quadratic_probes':len(list(quadratic_probes())),
                       'mismatched_probe_count':len(failures),
                       'first_failure_examples':failures[:4],
                       'physical_rgbeta_quartic_basis_certified':False}
    return {
        'format':'a-minus4-source-backed-gauge-quartic-batch-v1',
        'model':'A_-4',
        'github_source_inspection_baseline':'54857efafc7261f0e4c2bd0802f23e61dc3c165a',
        'local_git_head':local_git_revision(root),
        'inputs_sha256':{str(NATIVE):native_sha,str(MASTER):master_sha},
        'native_coordinate_count':len(native['summed_unscaled']),
        'master_coordinate_count':len(master['conditional_master_weight_basis_independent_entries']),
        'source_tensor_frobenius_norm_gev_inv':float(np.linalg.norm(cm)),
        'master_transport_relative_max_error':source_transport,
        'master_transport_relative_norm_error':float(norm_err),
        'triplet_unitarity_max_error':unitarity,
        'unconverted_native_using_master_generators_Ward':naive,
        'converted_master_Ward':converted,
        'gauge_only_six_term_checks':channel,
        'illustrative_fixed_gauge_leading_log':{'gY':gY,'g2':g2,'muF_gev':1e5,
                                                 'muS_gev':1e3,
                                                 'delta_C_over_C_first_order':gauge_log,
                                                 'full_physical_RGE':False},
        'quartic_global_rgbeta_reference_comparison':quartic_summary,
        'physical_Matchete_Weyl_normalization_verified':False,
        'complete_gauge_one_loop_loop_derivation_verified':False,
        'full_EFT1_coupling_trajectory_certified':False,
        'physical_C5_calculated':False,
        'source_files_modified':False,
    }


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--output',type=Path,
                   default=Path('output/a_minus4_next_batch/gauge_ward_audit.json'))
    args=p.parse_args(argv)
    root=args.root.resolve()
    result=report(root)
    out=args.output if args.output.is_absolute() else root/args.output
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Native 60 -> master 36 transport max relative:',
          f'{result["master_transport_relative_max_error"]:.3e}')
    print('Converted SU2 Ward max:',
          f'{max(result["converted_master_Ward"]["SU2"]):.3e}')
    print('U1 and SU2 gauge-only coefficients:',
          result['gauge_only_six_term_checks']['U1']['gauge_beta_over_g_squared_eigenvalue']['real'],
          result['gauge_only_six_term_checks']['SU2']['gauge_beta_over_g_squared_eigenvalue']['real'])
    print('RGBeta scalar-polynomial mismatched probes:',
          result['quartic_global_rgbeta_reference_comparison']['mismatched_probe_count'])
    print('Physical A(-4) Wilson boundary and complete one-loop EFT1: NOT CERTIFIED')
    print('Audit JSON:',out)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
