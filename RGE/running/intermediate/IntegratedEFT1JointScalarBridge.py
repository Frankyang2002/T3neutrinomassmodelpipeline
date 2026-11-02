"""Read-only all-F-first EFT1 seven-term source + simultaneous both-scalar bubble.

Explicitly *conditional* per unit unknown Matchete GammaCC/Weyl factor.  Does
not claim the UV hard one-loop boundary or full dimension-five operator closure.
This bridge reuses the repository's seven-term master and joint-scalar bubble.
"""
from __future__ import annotations

from hashlib import sha256
from itertools import permutations
from math import exp, log, pi
from pathlib import Path
from typing import Any, Mapping
import json
import numpy as np

from RGE.running.intermediate.SourceBackedEFT1Batch import (
    COUPLING_FORMAT, _validate_couplings, evaluate as seven_term_evaluate,
    required_couplings, read_json,
)
from RGE.running.intermediate.GeneralEFT1C5Mixing import (
    LOOP, TERMS, calculate_master_beta, _read_quartic, load_boundary, encode_matrix,
)
from RGE.running.intermediate.FullSU2SourceResponse import project_weinberg
from RGE.running.intermediate.ProvisionalSourceKernel import (
    SourceKernelError, kernel_for_boundary,
)

FORMAT = 't3-eft1-seven-term-joint-two-scalar-conditional-v1'
CONTRACT_FORMAT = 't3-eft1-seven-term-joint-two-scalar-contract-v1'
WARD_TOLERANCE = 1e-8


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SourceKernelError(message)


def complex_matrix(value: Mapping[str, Any], label: str) -> np.ndarray:
    require(isinstance(value, Mapping) and set(value) == {'real', 'imag'},
            f'{label}: real/imag complex matrix required')
    try:
        z = np.asarray(value['real'], float) + 1j * np.asarray(value['imag'], float)
    except (ValueError, TypeError) as exc:
        raise SourceKernelError(f'{label}: invalid numerical matrix') from exc
    require(z.shape == (3, 3) and np.isfinite(z).all(),
            f'{label}: finite 3x3 matrix required')
    return z


def _symmetric(z: np.ndarray, label: str) -> None:
    require(z.shape == (3, 3) and np.isfinite(z).all(), f'{label}: invalid 3x3 matrix')
    require(np.linalg.norm(z-z.T) < 1e-10*max(np.linalg.norm(z), 1e-30),
            f'{label}: neutrino flavours must be complex symmetric (not Hermitian)')


def _same(a: np.ndarray, b: np.ndarray, label: str, *, rtol: float=2e-10) -> float:
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(),
            f'{label}: dimensions/finite values disagree')
    error = float(np.max(abs(a-b))) if a.size else 0.0
    scale = max(float(np.max(abs(a))) if a.size else 0.,
                float(np.max(abs(b))) if b.size else 0., 1e-30)
    require(error <= rtol * scale + 1e-28,
            f'{label}: source and scalar background differ by {error:g} (scale {scale:g})')
    return error


def _real_channels(couplings: Mapping[str, Any], ns: int) -> dict[str, np.ndarray]:
    """Expand fully symmetric 4-real-scalar quartics, summing no channels here."""
    result = {}
    for channel in couplings['quartic_channels']:
        name = channel['label']
        require(name not in result, 'Duplicate quartic channel')
        result[name] = _read_quartic(channel['components'], ns)
    return result


def _master_source_heavy(context) -> np.ndarray:
    """Reorder C[i,j,a,b] -> C[p,q,lepton_i,lepton_j,heavy_a,heavy_b]."""
    n = context.ns - 4
    out = np.zeros((3,3,2,2,n,n),complex)
    total = sum((c for _,c in context.operators),
                np.zeros((context.nf,context.nf,context.ns,context.ns), complex))
    for p in range(3):
        for q in range(3):
            for i in range(2):
                for j in range(2):
                    out[p,q,i,j] = total[context.leptons[p][i],
                                          context.leptons[q][j],4:,4:]
    return out


def _project(amplitude: np.ndarray, name: str, *, strict: bool=True) -> tuple[np.ndarray,float]:
    r = project_weinberg(amplitude)
    mat = complex_matrix(r['projected_C5_beta_16pi2_per_NGamma'], name)
    _symmetric(mat, name)
    if strict:
        require(r['SU2_Ward_tensor_algebra_pass'] and
                r['weinberg_residual_relative'] <= WARD_TOLERANCE,
                f'{name}: full LLHH tensor fails Weinberg SU2 Ward residual '
                f"{r['weinberg_residual_relative']:.3g}; no physical C5 projection")
    return mat,float(r['weinberg_residual_relative'])


def _contract(C: np.ndarray, q: np.ndarray) -> np.ndarray:
    return np.einsum('fgijab,ABab->fgijAB',C,q,optimize=True)


def _rotate(C: np.ndarray, U: np.ndarray) -> np.ndarray:
    return np.einsum('fgijab,ak,bl->fgijkl', C, U, U, optimize=True)


def _backrotate(q_matching: np.ndarray, U: np.ndarray) -> np.ndarray:
    return np.einsum('ABkl,ak,bl->ABab',q_matching,U,U,optimize=True)


def _bubble(C_matching: np.ndarray, q_matching: np.ndarray, bg, *, mu:float|None=None):
    from Numerical.joint_scalar_matching.core import scalar_bubble_matching
    return scalar_bubble_matching(C_matching, q_matching, bg.mass2,
                                  bg.mu_s if mu is None else mu,
                                  reference_scale=bg.reference_scale)


def _bubble_c5(amplitude: np.ndarray, label: str, *, strict: bool) -> tuple[np.ndarray,float]:
    """Project in the *master* (1/2)C5 convention, not joint module's O5 convention.

    The joint module's Weinberg template is exactly TWICE the master template.
    Equivalently master_C5 = 2 * joint_module_projected_coefficient.
    """
    from Numerical.joint_scalar_matching.core import project_weinberg as project_joint
    master,r = _project(amplitude,label,strict=strict)
    joint, r_joint = project_joint(amplitude, strict=False)
    _same(master, 2*joint, f'{label}: joint-to-master operator factor of two')
    require(abs(r-r_joint)<=2e-11,
            f'{label}: different Weinberg tensor geometry, not just normalization')
    return master,r


def contracts(boundary: Mapping[str, Any], boundary_sha: str) -> dict:
    c = required_couplings(boundary,boundary_sha)
    n = boundary['real_scalar_dimension']-4
    return {'format':CONTRACT_FORMAT, 'model':boundary['model'],
            'status':'SourceValidated_ExternalNumericalInputsMissing',
            'source_boundary_sha256':boundary_sha,
            'mu_F_gev':boundary['mu_F_gev'],
            'required_master_couplings_format':COUPLING_FORMAT,
            'required_master_couplings':c,
            'required_joint_scalar_background':{
                'format':'t3-joint-scalar-conditional-background-v1',
                'model':boundary['model'],
                'source_boundary_sha256':boundary_sha,
                'numerical_array_shapes':{
                    'source_to_matching_real_U':[n,n],
                    'mass_squared_gev2_muS':[n,n],
                    'lambda_HHSS_muS':[4,4,n,n],
                    'lambda_HHSS_channels_muS':['k',4,4,n,n],
                },
                'required_channel_consistency':'quartic_channels labels and fully normalized HHSS components must match the master inputs after source/matching basis rotation',
                'required_scalar_mass_origin':'unbroken_electroweak or explicitly declared EWSB_resummed',
                'one_common_scalar_threshold':True,
            },
            'physical_source_spinor_factor':None,
            'complete_one_loop_UV_to_SMEFT_certified':False,
            'note':'This contract is not a fabricated mass/quartic/gauge input. No numerical C5 yet.'}


def evaluate(boundary_path: Path, couplings_path: Path, background_path: Path,
             *, trajectory_points: int=21) -> dict:
    """Combine 7-term frozen one-loop source and ONE joint scalar bubble, per NΓ.

    Same explicit couplings at muF/muS are required for the tested fixed-order
    O(hbar) scale-cancellation identity. A full running EFT1 background needs an
    independent higher-order transport and a revised contract.
    """
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.large_batch import load_background

    boundary_path = Path(boundary_path)
    boundary = read_json(boundary_path)
    bsha = sha256(boundary_path.read_bytes()).hexdigest()
    kernel_for_boundary(boundary, source_sha256=bsha)
    source = load_source(boundary_path, expected_model=boundary['model'])
    require(source.hash == bsha and source.normalization_certified is False,
            'Unexpected/certified spinor boundary: re-audit GammaCC externally')
    couplings = read_json(couplings_path)
    payload, channel_inventory = _validate_couplings(couplings,boundary,bsha)
    ctx = load_boundary(payload,allow_synthetic=True)
    bg = load_background(background_path,source)
    require(bg.model == source.model and
            abs(bg.mu_s/ctx.mu_s-1) <= 1e-12 and
            abs(bg.mu_f/ctx.mu_f-1) <= 1e-12,
            'Scalar background and master RGE must have exactly the same muF and muS')
    require(bg.mass_origin=='unbroken_electroweak',
            'Ward-certified dimension-five projection requires unbroken-EW masses')
    require(couplings['coverage_assertion'].get('fixed_order_quartic_identical_muF_muS') is True,
            'This fixed-order bridge requires explicit frozen quartic assertion')
    n = len(source.labels)
    require(ctx.ns == n + 4 and bg.quartic.shape==(4,4,n,n),
            'All real components of both heavy scalar multiplets required')
    _same(_master_source_heavy(ctx), source.raw, 'Matchete PL source -> master Weyl components')

    U = bg.U
    qsource = ctx.quartic[:4,:4,4:,4:]
    qtgt = _backrotate(bg.quartic,U)
    qerr = _same(qsource,qtgt,'complete HHSS quartet master-source to muS matching basis')
    qchannels = _real_channels(couplings,ctx.ns)
    require(len(qchannels)==len(bg.names) and set(qchannels)==set(bg.names),
            'Every named quartic channel must be preserved, with no duplications')
    channel_q = {}
    channel_err = {}
    for name,qt in zip(bg.names,bg.channels):
        q = qchannels[name][:4,:4,4:,4:]
        channel_err[name]=_same(q,_backrotate(qt,U),f'quartic channel {name}')
        channel_q[name]=q

    seven = seven_term_evaluate(boundary,bsha,couplings,trajectory_points=trajectory_points)
    require(seven.get('N_GammaCC_to_master_weyl') is None and
            seven.get('physical_C5_computed') is False,
            'Upstream diagnostics falsely assert physical normalization')
    beta_full = complex_matrix(seven['conditional_total_beta_C5_16pi2_per_NGamma'],'7term total beta')
    delta_full = complex_matrix(seven['conditional_total_delta_C5_log_per_NGamma_gev_inv'],
                                '7term total log')
    # Reject a neutrino-only pivot masquerading as an SU2-invariant total.
    full_master_beta = np.zeros((ctx.nf,ctx.nf,ctx.ns,ctx.ns),complex)
    for _,wilson in ctx.operators:
        terms = calculate_master_beta(ctx,wilson)
        for term in TERMS: full_master_beta += terms[term]
    weak_beta = np.zeros((3,3,2,2,4,4),complex)
    for p in range(3):
        for q in range(3):
            for i in range(2):
                for j in range(2):
                    weak_beta[p,q,i,j] = full_master_beta[ctx.leptons[p][i],
                             ctx.leptons[q][j]][np.ix_(ctx.higgs,ctx.higgs)]
    full_proj,ward_full = _project(weak_beta,'all seven terms combined',strict=True)
    _same(full_proj,beta_full,'7term full SU2 tensor versus master projected beta')
    beta_quartic,ward_q = _project(_contract(source.raw,qsource),'quartic beta',strict=True)
    beta_master_q = complex_matrix(seven['quartic_only_SU2_projection']
                ['projected_C5_beta_16pi2_per_NGamma'],'master quartic beta')
    _same(beta_quartic,beta_master_q,'7term master and independent scalar quartic beta')
    _same(LOOP*log(bg.mu_s/bg.mu_f)*beta_full,delta_full,'seven-term fixed-order logarithm')

    cmatch = _rotate(source.raw,U)
    bubble = _bubble(cmatch,bg.quartic,bg)
    bubble_c5,ward_b = _bubble_c5(bubble['amplitude'],'joint scalar bubble',strict=True)
    _same(beta_quartic,_project(_contract(cmatch,bg.quartic),'rotated beta')[0],
          'quartic beta changes under orthogonal scalar rotation')
    # Check the derivative only for the scalar-quartic channel. The other six
    # anomalous-dimension terms have their own light-field RGE scale partition.
    shifted_mu = bg.mu_s*exp(0.03125)
    shifted = _bubble(cmatch,bg.quartic,bg,mu=shifted_mu)
    shifted_c5,_ = _bubble_c5(shifted['amplitude'],'shifted scalar bubble',strict=True)
    residual = shifted_c5-bubble_c5 + LOOP*log(shifted_mu/bg.mu_s)*beta_quartic
    scale = max(np.linalg.norm(bubble_c5),
                abs(LOOP*log(shifted_mu/bg.mu_s))*np.linalg.norm(beta_quartic),1e-30)
    scale_error = float(np.linalg.norm(residual)/scale)
    require(scale_error<=2e-9, f'quartic log versus bubble scale cancellation fails: {scale_error:g}')

    per_sources=[]
    for i,sid in enumerate(source.source_indices):
        C = source.by_source[i]
        qbeta,qres = _project(_contract(C,qsource),f'PL {sid} quartic',strict=False)
        Cb = _rotate(C,U)
        bm = _bubble(Cb,bg.quartic,bg)
        bm5,bres = _bubble_c5(bm['amplitude'],f'PL {sid} bubble',strict=False)
        row = seven['source_resolved_seven_term_results'][i]
        require(row['source_id']==f'Matchete_PL_{sid}', 'Matchete PL source order changed')
        d5=complex_matrix(row['delta_C5_log_gev_inverse'],f'PL {sid} delta')
        per_sources.append({'source_index':sid,'quartic_beta_master_per_N':encode_matrix(qbeta),
           'seven_term_delta_C5_per_N_gev_inv':encode_matrix(d5),
           'joint_scalar_bubble_master_per_N_gev_inv':encode_matrix(bm5),
           'combined_increment_per_N_gev_inv':encode_matrix(d5+bm5),
           'individual_full_SU2_residuals':[qres,bres]})
    per_quartics=[]
    for name,q in channel_q.items():
        qbeta,qres = _project(_contract(source.raw,q),f'{name} beta',strict=False)
        matching_q = bg.channels[list(bg.names).index(name)]
        bm = _bubble(cmatch,matching_q,bg)
        bm5,bres = _bubble_c5(bm['amplitude'],f'{name} bubble',strict=False)
        per_quartics.append({'label':name,'quartic_beta_master_per_N':encode_matrix(qbeta),
             'quartic_delta_C5_log_per_N_gev_inv':encode_matrix(LOOP*log(bg.mu_s/bg.mu_f)*qbeta),
             'joint_scalar_bubble_master_per_N_gev_inv':encode_matrix(bm5),
             'combined_quartic_increment_per_N_gev_inv':encode_matrix(
                 LOOP*log(bg.mu_s/bg.mu_f)*qbeta+bm5),
             'individual_full_SU2_residuals':[qres,bres]})
    for field,total,key in [(per_sources,delta_full,'seven_term_delta_C5_per_N_gev_inv'),
                            (per_sources,bubble_c5,'joint_scalar_bubble_master_per_N_gev_inv'),
                            (per_quartics,beta_quartic,'quartic_beta_master_per_N'),
                            (per_quartics,bubble_c5,'joint_scalar_bubble_master_per_N_gev_inv')]:
        _same(sum((complex_matrix(x[key],key) for x in field),np.zeros((3,3),complex)),
              total,f'{key} source/channel recomposition')
    for mat,name in ((beta_full,'7term beta'),(delta_full,'7term log'),
                     (bubble_c5,'bubble'),(delta_full+bubble_c5,'combined increment')):
        _symmetric(mat,name)
    return {'format':FORMAT,'status':'ConditionalAllSevenTermLogPlusJointScalarBubble_PerUnknownGamma',
        'model':source.model,'source_boundary_sha256':source.hash,
        'coupling_input_sha256':sha256(Path(couplings_path).read_bytes()).hexdigest(),
        'background_json_sha256':sha256(Path(background_path).read_bytes()).hexdigest(),
        'background_npz_sha256':sha256((Path(background_path).parent / read_json(background_path)['tensor_npz']).read_bytes()).hexdigest(),
        'mu_F_gev':bg.mu_f,'mu_S_gev':bg.mu_s,'scheme':'MSbar',
        'N_GammaCC_to_master_Weyl':None,'physical_C5_computed':False,
        'C5_muF_direct_hard_matching':'NOT PROVIDED (not assumed zero)',
        'one_loop_input_order':'C_LLSS^(tree) times one renormalizable loop',
        'scalar_threshold_count':1,'both_heavy_scalar_multiplets_integrated_together':True,
        'C5_7term_beta_16pi2_per_N':encode_matrix(beta_full),
        'C5_quartic_beta_16pi2_per_N':encode_matrix(beta_quartic),
        'C5_EFT1_7term_log_per_N_gev_inv':encode_matrix(delta_full),
        'C5_joint_scalar_bubble_master_per_N_gev_inv':encode_matrix(bubble_c5),
        'C5_LLSS_induced_increment_to_SMEFT_per_N_gev_inv':encode_matrix(delta_full+bubble_c5),
        'C5_quartic_only_increment_per_N_gev_inv':encode_matrix(
              LOOP*log(bg.mu_s/bg.mu_f)*beta_quartic+bubble_c5),
        'per_Matchete_PL_source':per_sources,
        'per_independent_quartic_channel':per_quartics,
        'per_scalar_mass_eigenpair':[
          {'pair':k,'bubble_C5_master_per_N_gev_inv':encode_matrix(
              _bubble_c5(v,k,strict=False)[0])}
          for k,v in bubble['pair_amplitudes'].items()],
        'scalar_eigenstate_masses_gev':np.sqrt(bubble['eigenvalues_sq_gev2']).tolist(),
        'full_SU2_seven_term_residual':ward_full,
        'full_SU2_quartic_residual':ward_q,
        'full_SU2_joint_bubble_residual':ward_b,
        'quartic_muS_log_cancellation_relative_residual':scale_error,
        'master_vs_background_HHSS_max_abs_error':qerr,
        'channel_max_abs_errors':channel_err,
        'seven_term_provenance':seven['conditional_master_beta_by_rge_term'],
        'physical_complete_one_loop_certified':False,
        'remaining_blockers':['signed Matchete GammaCC to canonical Weyl C5 normalization',
            'independently certified physical EFT1 quartic/Yukawa/gauge dictionary',
            'hard one-loop fermion matching and independent direct C5 at muF',
            'other one-loop EFT1 operator insertions and finite threshold counterterms'],
        'note':('Seven-term *source* log and one joint scalar bubble are added once, in the '
            'master (1/2)C5 convention. The joint matcher uses a Weinberg template '
            'twice as large; its projected bubble coefficient is multiplied by two. '
            'All results are increments per an UNKNOWN N_Gamma and omit independent C5(muF).')}
