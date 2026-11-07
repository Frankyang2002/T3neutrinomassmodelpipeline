"""Nondegenerate NF=3 lift of the *saved* ordinary-T3 hard scalar template.

This is deliberately restricted to the exact source-verified D(+2)/A(-4)
Matchete expression.  It computes the finite hard scalar contribution and
its direct EFT1 LLSS log, not the inherited one-loop sector or complete C5.

UV -> integrate out ALL F -> EFT1(S1,S2) -> integrate out ALL scalars -> SMEFT.
The input template is O(hbar), hbar = 1/(16*pi*pi). Its residual 1/epsilon
is EXPORTED, not silently discarded as a proved renormalized observable.
"""
from __future__ import annotations

from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

MODELS = ('D_+2', 'A_-4')
TAGS = {'D_+2': 'D_p2', 'A_-4': 'A_m4'}
# Whitespace-insensitive fingerprint of the exact six-term Mathematica
# expression supplied by the two successful local Wolfram matching stages.
TEMPLATE_SHA256 = 'b9fbc4c11e366602e912b665cf3cb08ecaf7c280aa2998b93036e7c851ec395f'
LOOP = 1.0 / (16 * math.pi**2)
DEFAULT_VEV_GEV = 246.22


class HardFlavorError(ValueError):
    """The requested flavor lift is not justified by the supplied source."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise HardFlavorError(message)


def read_json(path: Path) -> dict[str, Any]:
    path = Path(path)
    require(path.is_file(), f'Missing {path}')
    try:
        obj = json.loads(path.read_text(encoding='utf-8-sig'))
    except (OSError, ValueError) as exc:
        raise HardFlavorError(f'Unreadable JSON {path}: {exc}') from exc
    require(isinstance(obj, dict), f'Expected JSON object: {path}')
    return obj


def file_sha(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def encode_complex(a: np.ndarray) -> dict:
    a = np.asarray(a, np.complex128)
    require(a.shape == (3, 3) and np.isfinite(a).all(), 'Invalid 3x3 complex C5')
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def decode_complex(raw: Any, name: str) -> np.ndarray:
    require(isinstance(raw, dict) and 'real' in raw and 'imag' in raw,
            f'{name}: expected real/imag complex matrix')
    try:
        value = np.asarray(raw['real'], float) + 1j * np.asarray(raw['imag'], float)
    except (ValueError, TypeError) as exc:
        raise HardFlavorError(f'{name}: malformed complex matrix') from exc
    require(value.shape == (3, 3) and np.isfinite(value).all(),
            f'{name}: expected finite 3x3 matrix')
    return value


def typed_complex_matrix(raw: Any, name: str) -> np.ndarray:
    require(isinstance(raw, list) and len(raw) == 3, f'{name}: expected 3 rows')
    result = np.empty((3, 3), np.complex128)
    for i, row in enumerate(raw):
        require(isinstance(row, list) and len(row) == 3, f'{name}: expected 3 columns')
        for j, value in enumerate(row):
            require(isinstance(value, dict) and {'real', 'imag'} <= set(value),
                    f'{name}[{i},{j}]: missing real/imag')
            try:
                result[i, j] = complex(float(value['real']), float(value['imag']))
            except (TypeError, ValueError) as exc:
                raise HardFlavorError(f'{name}: invalid numeric entry') from exc
    require(np.isfinite(result).all(), f'{name}: nonfinite entries')
    return result


def scalar_complex(raw: Any, name: str) -> complex:
    if isinstance(raw, dict):
        require({'real', 'imag'} <= set(raw), f'{name}: missing complex components')
        value = complex(float(raw['real']), float(raw['imag']))
    else:
        require(not isinstance(raw, bool), f'{name}: invalid Boolean')
        value = complex(raw)
    require(np.isfinite(value.real) and np.isfinite(value.imag), f'{name}: nonfinite')
    return value


def positive_real(raw: Any, name: str) -> float:
    require(not isinstance(raw, bool), f'{name}: invalid Boolean')
    try:
        value = float(raw)
    except (ValueError, TypeError) as exc:
        raise HardFlavorError(f'{name}: expected real number') from exc
    require(math.isfinite(value) and value > 0, f'{name}: must be finite and positive')
    return value


def mass_kernel(m1_sq: float, m2_sq: float, mu: float) -> float:
    r"""1 + (m1^2 log(mu^2/m1^2) - m2^2 log(mu^2/m2^2))/(m1^2-m2^2).

    Analytic equal-mass limit: log(mu^2/m^2).  Evaluated stably when
    m1^2 is very close to m2^2 without cancellation of large numbers.
    """
    a = positive_real(m1_sq, 'MS1 squared')
    b = positive_real(m2_sq, 'MS2 squared')
    scale = positive_real(mu, 'scalar matching scale')
    t = (b - a) / a
    if abs(t) < 1e-5:
        # (1+t)*log(1+t)/t = 1+t/2-t^2/6+t^3/12-t^4/20+...
        shift = -t/2 + t*t/6 - t**3/12 + t**4/20 - t**5/30
        return float(math.log(scale**2 / a) + shift)
    return float(1.0 + math.log(scale**2 / a)
                 - (1.0 + t) * math.log1p(t) / t)


def verify_hard_stage(folder: Path, model: str, *, source_root: Path | None = None) -> dict:
    require(model in MODELS, f'Model not supported by this source: {model}')
    folder = Path(folder)
    report = read_json(folder/'hard_scalar_report.json')
    audit = read_json(folder/'hard_scalar_input_audit.json')
    provenance = read_json(folder/'provenance.json')
    require(report.get('Status') == 'HardScalarCandidateExported' and
            report.get('Model') == model and
            report.get('HardOneLoopMatched') is True and
            report.get('HardHbarSectorExtracted') is True and
            report.get('PerHeavyMassEigenvalueTemplate') is True and
            report.get('FermionMatchingRecomputed') is False,
            f'{model}: not a verified successful native one-eigenvalue hard stage')
    require(report.get('CompleteC5Assembled') is False and
            report.get('SU2NormalizationIndependentlyCertified') is False and
            report.get('ThreeHeavyMassNumericalSumImplemented') is False,
            f'{model}: template verification flags unexpected')
    require(audit.get('Status') == 'Passed' and audit.get('Model') == model and
            audit.get('CheckLagrangianPassed') is True and
            audit.get('DiagonalNFlavorMassIndexEliminatedForTemplate') is True and
            audit.get('PerEigenvalueHeavyMassSumRequired') is True and
            audit.get('TreeWilsonTermCount') == 2 and
            audit.get('OddHeavyScalarDegreeAbsent') is True,
            f'{model}: scalar threshold input audit not passed')
    require(provenance.get('model') == model and provenance.get('status') == 'Ready' and
            provenance.get('three_heavy_mass_numerical_sum_required') is True and
            provenance.get('indexed_fermion_mass_factors') == 2,
            f'{model}: mismatched template provenance')
    path = folder/'hard_c5_candidate.txt'
    require(path.is_file(), f'{model}: missing hard candidate')
    content = path.read_text(encoding='utf-8-sig')
    normalized = ''.join(content.split())
    digest = sha256(normalized.encode('utf-8')).hexdigest()
    require(digest == TEMPLATE_SHA256,
            f'{model}: unexpected Mathematica hard C5 expression {digest}; '
            'this restricted analytic lift MUST NOT guess a new expression')
    require(normalized.count('hbar') == 6 and
            normalized.count('Coupling[MF,{},0]') == 6 and
            normalized.count('ϵ') == 2 and
            'Coupling[lambdaT3,{},0]' in normalized and
            'Log[μbar2/Coupling[MS1,{},0]^2]' in normalized and
            'Log[μbar2/Coupling[MS2,{},0]^2]' in normalized,
            f'{model}: unexpected loop/pole/scale source structure')
    if source_root is not None:
        fps = provenance.get('source_fingerprints', {})
        require(isinstance(fps, dict) and len(fps) >= 5,
                f'{model}: no recorded fermion source fingerprints')
        root = Path(source_root) / TAGS[model]
        for relative, expected in fps.items():
            if relative == 'matching_wl':
                continue   # This is the archived matching code, not the F-only release.
            local = root / relative
            require(local.is_file() and file_sha(local) == expected,
                    f'{model}: saved F-only source has changed: {local}')
    return {'status': 'VerifiedOneEigenvalueTemplate',
            'model': model, 'source_expression_sha256_normalized': digest,
            'source_expression_sha256_raw': file_sha(path),
            'source_hard_report_sha256': file_sha(folder/'hard_scalar_report.json'),
            'source_provenance_sha256': file_sha(folder/'provenance.json')}


def running_parameters(doc: dict, model: str) -> dict:
    require(doc.get('model') == model and
            doc.get('status') == 'NumericalRenormalisableRunningOnly',
            f'{model}: not the expected illustrative full-general RG snapshot')
    scales = doc.get('scales_gev')
    require(isinstance(scales, dict), f'{model}: scales absent')
    mu_f = positive_real(scales.get('F'), 'muF')
    mu_s = positive_real(scales.get('S'), 'muS')
    require(mu_f > mu_s, f'{model}: matching hierarchy muF > muS required')
    uv = doc.get('uv_at_fermion_threshold')
    scalar = doc.get('eft1_at_scalar_threshold')
    scalar_f = doc.get('eft1_at_fermion_threshold')
    require(all(isinstance(x, dict) for x in (uv, scalar, scalar_f)),
            f'{model}: incomplete renormalisable RG snapshots')
    y1 = typed_complex_matrix(uv.get('y1'), 'y1(muF)')
    y2 = typed_complex_matrix(uv.get('y2'), 'y2(muF)')
    matrix = typed_complex_matrix(uv.get('MF'), 'MF(muF)')
    diag = np.diag(matrix)
    norm = max(float(np.linalg.norm(diag)), 1.0)
    require(np.linalg.norm(matrix-np.diag(diag)) < 1e-12*norm,
            f'{model}: MF is not diagonal; a mass-basis Yukawa rotation is required')
    require(np.max(np.abs(diag.imag)) < 1e-12*norm and np.all(diag.real > 0),
            f'{model}: diagonal MF must have positive real eigenvalues')
    mf = diag.real.copy()
    require(float(np.max(mf)/np.min(mf)) > 1.00001,
            f'{model}: test specifically expects nondegenerate MF')
    m1sq = positive_real(scalar.get('mS1Sq'), 'mS1Sq(muS)')
    m2sq = positive_real(scalar.get('mS2Sq'), 'mS2Sq(muS)')
    l = scalar_complex(scalar.get('lambdaT3'), 'lambdaT3(muS)')
    return dict(mu_f=mu_f, mu_s=mu_s, y1=y1, y2=y2, mf=mf,
                m1sq=m1sq, m2sq=m2sq, lam=l,
                m1sq_f=positive_real(scalar_f.get('mS1Sq'), 'mS1Sq(muF)'),
                m2sq_f=positive_real(scalar_f.get('mS2Sq'), 'mS2Sq(muF)'),
                lam_f=scalar_complex(scalar_f.get('lambdaT3'), 'lambdaT3(muF)'))


def source_pole_matrix(y1: np.ndarray, y2: np.ndarray, mf: np.ndarray,
                       lam: complex) -> tuple[np.ndarray, list[np.ndarray]]:
    """Only the exact holomorphic ybar1*ybar2 term in the source template.

    The source term has external (i,j), and symmetric flavour projection uses
    0.5*(i,j + j,i); its off-diagonal overall normalization is provisional.
    """
    rows = []
    for r in range(3):
        unsym = np.outer(np.conj(y1[:, r]), np.conj(y2[:, r]))
        sym = .5*(unsym + unsym.T)
        rows.append((-2.0 / math.sqrt(3))*LOOP*lam*sym/mf[r])
    return sum(rows, np.zeros((3, 3), complex)), rows


def metrics(c5: np.ndarray, vev: float) -> dict:
    require(c5.shape == (3, 3) and np.isfinite(c5).all(), 'C5 nonfinite')
    norm = max(float(np.linalg.norm(c5)), 1e-40)
    symmetry = float(np.linalg.norm(c5-c5.T)/norm)
    require(symmetry <= 1e-11, f'flavour C5 is not symmetric: {symmetry}')
    return {'c5_frobenius_gev_inverse': float(np.linalg.norm(c5)),
            'c5_symmetry_relative_residual': symmetry,
            'conditional_takagi_masses_ev_ascending':
            np.sort(np.linalg.svd(c5, compute_uv=False) * vev**2/2*1e9).tolist()}


def evaluate(*, model: str, stage_folder: Path, running_doc: dict,
             vev_gev: float = DEFAULT_VEV_GEV,
             source_root: Path | None = None) -> dict:
    """Numerical fixed-order LLSS *subset*: native finite scalar hard + direct log.

    Full C5 requires inherited first-threshold one-loop and independently
    established MSbar counterterms. Do not add archived conditional LLSS to
    this result: it contains the same log and finite scalar contribution.
    """
    source = verify_hard_stage(stage_folder, model, source_root=source_root)
    p = running_parameters(running_doc, model)
    vev = positive_real(vev_gev, 'Higgs vev')
    pole, by_flavour = source_pole_matrix(p['y1'], p['y2'], p['mf'], p['lam'])
    kernel = mass_kernel(p['m1sq'], p['m2sq'], p['mu_s'])
    finite_hard = pole * kernel
    direct_log = (-2.0*math.log(p['mu_s']/p['mu_f']))*pole
    subtotal = finite_hard + direct_log

    # The scale test varies ONLY the arbitrary scalar matching mu, holding
    # F-boundary Yukawa/couplings and scalar masses fixed. In a fixed-order
    # truncation the hard + direct log has zero derivative in ln(mu).
    probe_mu = p['mu_s']*math.exp(.03125)
    shifted = (pole*mass_kernel(p['m1sq'], p['m2sq'], probe_mu)
               - 2*pole*math.log(probe_mu/p['mu_f']))
    denominator = max(float(np.linalg.norm(subtotal)),
                      float(np.linalg.norm(finite_hard)),
                      float(np.linalg.norm(direct_log)), 1e-40)
    scale_residual = float(np.linalg.norm(shifted-subtotal)/denominator)
    require(scale_residual < 2e-11,
            f'{model}: hard/direct-log matching scale cancellation failed {scale_residual}')
    closure = float(np.linalg.norm(subtotal-finite_hard-direct_log)/denominator)
    require(closure < 1e-12, f'{model}: finite plus log does not close')
    # Separate scheme benchmark with muF-frozen scalar couplings. It is NOT
    # added to the muS evolved result. Difference represents renormalisable
    # running/resummation and unverified full EFT1 effects.
    frozen_pole, _ = source_pole_matrix(p['y1'],p['y2'],p['mf'],p['lam_f'])
    frozen_total = frozen_pole * (mass_kernel(p['m1sq_f'],p['m2sq_f'],p['mu_s'])
                                  - 2*math.log(p['mu_s']/p['mu_f']))
    return {'format':'t3-ordinary-native-hard-three-nondegenerate-flavour-subset-v1',
            'model':model, 'status':'HardThreeFlavorAndDirectLogComputed_NotFullC5',
            'sequence':'UV -> integrate ALL F -> EFT1(S1,S2) -> integrate ALL S -> SMEFT',
            'source_template':source,
            'matching_scales_gev':{'muF':p['mu_f'],'muS':p['mu_s']},
            'heavy_fermion_mass_eigenvalues_gev':p['mf'].tolist(),
            'scalar_evolved_at_muS':{'mS1_gev':math.sqrt(p['m1sq']),
                   'mS2_gev':math.sqrt(p['m2sq']),
                   'lambdaT3':{'real':p['lam'].real,'imag':p['lam'].imag}},
            'source_formula':{
                'hbar':'1/(16*pi^2)',
                'template_prefactor':'-2/sqrt(3) * hbar * lambdaT3',
                'bar_yukawa_convention':'complex-conjugate y1 and y2; lambdaT3 NOT conjugated',
                'flavour_sum':'sum_r conj(y1_ir) conj(y2_jr) / MF_r',
                'symmetry_projection':'(unsym_ij + unsym_ji)/2; offdiagonal physical normalization not independently certified',
                'pole':'P/epsilon',
                'finite_scalar_hard':'P * (1 + (m1sq*ln(muS^2/m1sq)-m2sq*ln(muS^2/m2sq))/(m1sq-m2sq))',
                'direct_EFT1_log':'-2*P*ln(muS/muF)',
                'log_running_provenance':'D(+2) and A(-4) source general EFT1 direct-running expressions agree in sign and coefficient; full inherited operators not included',
            },
            'raw_unsubtracted_hard_1_over_epsilon_residue_gev_inverse':encode_complex(pole),
            'hard_finite_scalar_gev_inverse':encode_complex(finite_hard),
            'direct_eft1_llss_log_gev_inverse':encode_complex(direct_log),
            'c5_hard_finite_plus_direct_log_subset_gev_inverse':encode_complex(subtotal),
            'frozen_scalar_parameter_subset_gev_inverse':encode_complex(frozen_total),
            'per_heavy_eigenstate':[{'r':r+1,'MF_gev':float(p['mf'][r]),
                'hard_pole_residue_gev_inverse':encode_complex(component),
                'hard_finite_gev_inverse':encode_complex(component*kernel),
                'direct_log_gev_inverse':encode_complex(-2*math.log(p['mu_s']/p['mu_f'])*component)}
                for r,component in enumerate(by_flavour)],
            'checks':{
                'independent_Wolfram_scalar_hard_template':True,
                'nondegenerate_three_eigenstate_sum_performed':True,
                'hard_finite_plus_direct_log_closure_residual':closure,
                'hard_plus_direct_log_matching_scale_variation_residual':scale_residual,
                'bare_hard_1_over_epsilon_pole_present':bool(np.linalg.norm(pole)>0),
                'pole_residue_norm_gev_inverse':float(np.linalg.norm(pole)),
                'rge_log_coefficient_over_hard_pole':'-2 for ln(muS/muF)',
                'actual_inherited_one_loop_added':False,
                'independent_MSbar_counterterm_proved':False,
                'physical_SU2_and_spinor_normalization_certified':False,
                'complete_sequential_C5_certified':False,
            },
            'hard_finite_metrics':metrics(finite_hard,vev),
            'hard_plus_running_subset_metrics':metrics(subtotal,vev),
            'vev_gev':vev,
            'interpretation':('Independent scalar hard one-loop finite term evaluated for three nondegenerate '
                'heavy masses, plus corresponding leading EFT1 direct logarithm. '
                'Pole is separately retained. This is a conditional LLSS sector, NOT a complete '
                'or independently renormalized physical C5; it overlaps earlier archived LLSS outputs.')}


def compare_archived_v6(result: dict, archived_path: Path) -> dict:
    """Independent existing real-component LLSS kernel parity; never add twice."""
    archived = read_json(archived_path)
    require(archived.get('model') == result['model'] and
            archived.get('complete_one_loop_C5_certified') is False,
            'Unexpected archived provisional LLSS result')
    response = archived.get('response')
    require(isinstance(response, dict) and response.get('model') == result['model'],
            'Wrong reference response')
    pairs = {
        'hard_finite':('hard_finite_scalar_gev_inverse',
                        'C5_joint_scalar_finite_bubble_per_N_gev_inverse'),
        'eft1_direct_log':('direct_eft1_llss_log_gev_inverse',
                            'C5_EFT1_induced_leading_log_per_N_gev_inverse'),
        'hard_plus_direct_log':('c5_hard_finite_plus_direct_log_subset_gev_inverse',
                                'C5_LLSS_total_one_loop_per_N_gev_inverse')
    }
    errors = {}
    for key,(local,remote) in pairs.items():
        a = decode_complex(result[local],local)
        b = decode_complex(response[remote],remote)
        denominator = max(float(np.linalg.norm(a)),float(np.linalg.norm(b)),1e-40)
        errors[key] = float(np.linalg.norm(a-b)/denominator)
    return {'status':'SourceBackedLLSSKernelAgreement'
                    if max(errors.values())<1e-9 else 'SourceBackedLLSSKernelMismatch',
            'reference_path':str(archived_path),
            'reference_sha256':file_sha(archived_path),
            'relative_frobenius_discrepancies':errors,
            'warning':'Comparison only; archived LLSS must NOT be added again to this hard+log subtotal',
            'absolute_SU2_spinor_normalization_proved':False}
