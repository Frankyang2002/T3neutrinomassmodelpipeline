"""Strict, lossless Stream-1/2 -> existing Workstream-C input adapter.

The producer's numerical tensors, not its display names, establish the bridge.
Native RGBeta beta-coordinate trajectories alone CANNOT determine the full
Matchete-normalised real quartic tensor or the scalar coordinate rotation.
Those objects must be supplied independently with explicit conventions.

Only C-owned temporary files are written; no A/B product is modified.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
import json
import math
import re
from typing import Any

import numpy as np

from .models import n_real, normalize_model, scalar_labels
from .pipeline import B_FORMAT, W_FORMAT, Q_CONVENTION, W_CONVENTION


class AdapterBlocked(ValueError):
    """Explicitly unavailable or inconsistent numerical producer input."""


@dataclass(frozen=True)
class PreparedPair:
    wilson: Path
    background: Path
    provenance: dict


def _sha(path: Path) -> str:
    h = sha256()
    with path.open('rb') as stream:
        for buf in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(buf)
    return h.hexdigest()


def _read_input(path: Path, *, nominal_npz: str) -> tuple[dict, dict, dict]:
    if not path.is_file():
        raise AdapterBlocked(f'Missing producer manifest: {path}')
    try:
        meta = json.loads(path.read_text(encoding='utf-8-sig'))
    except (ValueError, UnicodeError) as exc:
        raise AdapterBlocked(f'Invalid manifest {path}: {exc}') from exc
    if not isinstance(meta, dict):
        raise AdapterBlocked(f'{path}: expected JSON object')
    npzname = meta.get('tensor_npz', nominal_npz)
    if not isinstance(npzname, str) or not re.fullmatch(r'[A-Za-z0-9_.-]+\.npz', npzname):
        raise AdapterBlocked('Input NPZ must be a sibling basename')
    npz = path.with_name(npzname)
    if not npz.is_file():
        raise AdapterBlocked(f'Missing producer tensor NPZ: {npz}')
    digest = _sha(npz)
    for key in ('tensor_npz_sha256', 'npz_sha256', 'background_npz_sha256'):
        if key in meta and meta[key] != digest:
            raise AdapterBlocked(f'{key}: failed original producer checksum')
    outputs = meta.get('outputs', {})
    if isinstance(outputs, dict) and npz.name in outputs and outputs[npz.name] != digest:
        raise AdapterBlocked('Producer output-hash record disagrees with NPZ')
    with np.load(npz, allow_pickle=False) as loaded:
        arrays = {key: loaded[key].copy() for key in loaded.files}
    return meta, arrays, {'json_path':str(path), 'npz_path':str(npz),
                          'json_sha256':_sha(path), 'npz_sha256':digest}


def _numerical_array(value: Any, label: str, *, shape: tuple[int, ...] | None = None,
                     real: bool = False) -> np.ndarray:
    if value is None:
        raise AdapterBlocked(f'Missing required numerical array: {label}')
    a = np.asarray(value)
    if a.dtype.kind not in 'fiuc':
        raise AdapterBlocked(f'{label}: expected numeric real/complex array')
    if shape is not None and a.shape != shape:
        raise AdapterBlocked(f'{label}: expected {shape}; received {a.shape}')
    if not np.isfinite(a).all():
        raise AdapterBlocked(f'{label}: nonfinite data')
    if real and np.max(abs(a.imag)) > 1e-10 * max(1., float(np.max(abs(a.real)))):
        raise AdapterBlocked(f'{label}: imaginary components in real-field tensor')
    return np.asarray(a.real if real else a, dtype=np.float64 if real else np.complex128)


def _identity_fields(meta: dict, model: str, *, label: str) -> dict:
    """Require genuine shared provenance rather than deriving matching IDs from model."""
    if normalize_model(meta.get('model', '')) != model:
        raise AdapterBlocked(f'{label} model mismatch')
    ids = {}
    for key in ('matching_context_id', 'matched_seed_id', 'lagrangian_source_id'):
        value = meta.get(key)
        if not isinstance(value, str) or not value.strip():
            raise AdapterBlocked(f'{label} missing {key}; cannot assert same A/B input state')
        ids[key] = value
    mf = meta.get('mu_f_gev', meta.get('MF_gev'))
    ms = meta.get('mu_s_gev', meta.get('MS_gev'))
    try:
        mf, ms = float(mf), float(ms)
    except (TypeError, ValueError) as exc:
        raise AdapterBlocked(f'{label} missing explicit MF and MS scale') from exc
    if not (math.isfinite(mf) and math.isfinite(ms) and mf > ms > 0):
        raise AdapterBlocked(f'{label}: invalid MF > MS > 0')
    if abs(ms - 1e3) > 1e-9:
        raise AdapterBlocked(f'{label}: rough benchmark requires MS=1000 GeV, got {ms}')
    if meta.get('scheme') != 'MSbar' or meta.get('units') != 'GeV':
        raise AdapterBlocked(f'{label}: requires declared MSbar and GeV conventions')
    if type(meta.get('synthetic_fixture')) is not bool:
        raise AdapterBlocked(f'{label}: synthetic_fixture must be declared explicitly')
    return {**ids, 'mu_f_gev':mf, 'mu_s_gev':ms,
            'synthetic_fixture':meta['synthetic_fixture']}


def _same_input(w: dict, b: dict) -> None:
    for key in ('matching_context_id','matched_seed_id','lagrangian_source_id',
                'mu_f_gev','mu_s_gev','synthetic_fixture'):
        if w[key] != b[key]:
            raise AdapterBlocked(f'Cross-stream mismatch: {key} ({w[key]!r} != {b[key]!r})')
    if not w['synthetic_fixture'] and any('synthetic' in w[k].lower() for k in
                                         ('matching_context_id','matched_seed_id','lagrangian_source_id')):
        raise AdapterBlocked('Synthetic provenance cannot be relabelled as a genuine model')


def _source_wilson(meta: dict, z: dict, model: str, identity: dict) -> tuple[dict,dict]:
    n = n_real(model)
    fmt = meta.get('format')
    if fmt == W_FORMAT:
        if meta.get('canonical_wilson_convention') != W_CONVENTION:
            raise AdapterBlocked('Workstream A does not declare the C no-factorial Wilson basis')
        if 'C_psi2phi2_ms' not in z:
            raise AdapterBlocked('Canonical Wilson endpoint lacks full C_psi2phi2_ms')
        data = {k:np.asarray(v) for k,v in z.items()}
        out = meta.copy()
    elif fmt == 't3-parallel-a-wilson-output-v1':
        if (meta.get('all_F_integrated_out_at_MF') is not True or
            meta.get('both_S1_S2_remain_until_MS') is not True or
            meta.get('no_posthoc_fixed_order_leading_log') is not True or
            meta.get('C5_mutual_mixing_present') is not True):
            raise AdapterBlocked('Legacy A output is not a verified fermion-first full EFT1 endpoint')
        basis = meta.get('basis', {})
        if not isinstance(basis, dict):
            raise AdapterBlocked('Missing A full Weyl/real scalar basis')
        fermion_labels = basis.get('fermion_labels')
        scalars = basis.get('scalar_labels')
        if not isinstance(fermion_labels,list) or not isinstance(scalars,list):
            raise AdapterBlocked('Missing A ordered fermion/scalar label lists')
        if len(set(fermion_labels)) != len(fermion_labels) or len(set(scalars)) != len(scalars):
            raise AdapterBlocked('A basis contains duplicated field labels')
        if len(scalars) != 4+n or scalars[:4] != ['H[1,R]','H[1,I]','H[2,R]','H[2,I]']:
            raise AdapterBlocked('A Higgs basis is not canonical H(real,imag) interleaved')
        try:
            lep = [[fermion_labels.index(f'L[{p},{w}]') for w in (1,2)] for p in (1,2,3)]
        except ValueError as exc:
            raise AdapterBlocked('A full tensor has no complete three-generation lepton doublets') from exc
        c = _numerical_array(z.get('C_full'), 'A/C_full',
                             shape=(len(fermion_labels),)*2+(4+n,)*2)
        if meta.get('C5_projection_definition') != 'C5[p,q]=4*C[L[p,1],L[q,1],H[2,R],H[2,R]]':
            raise AdapterBlocked('A C5 normalization not known to agree with Workstream C')
        out = {**identity, 'format':W_FORMAT,'tensor_npz':'wilson_at_MS.npz',
               'canonical_wilson_convention':W_CONVENTION,
               'wilson_evolved_to_ms':True,'includes_eft1_weinberg':True,
               'legacy_direct_llss_to_c5_log_added':False,
               'scalar_basis_id':basis.get('basis_id'),
               'scalar_source_real_labels':scalars[4:],
               'lepton_indices':lep,
               'higgs_source_indices':list(range(4)),
               'heavy_source_indices':list(range(4,4+n)),
               'wilson_adm_certified':False, 'wilson_dim5_closure_certified':False,
               'signed_GammaCC_to_Weyl_certified':False,
               'source_A_format':fmt}
        data = {'C_psi2phi2_ms':c}
        if 'C5_EFT1_MS' in z:
            data['C5_EFT1_ms']=_numerical_array(z['C5_EFT1_MS'],'A/C5_EFT1_MS',shape=(3,3))
        # A's optional linearised leading log is diagnostic only and is NOT inserted.
    else:
        raise AdapterBlocked(f'Unrecognised A format {fmt!r}. Export a declared C-v1 endpoint')
    # The existing threshold preflight checks full Weyl pair symmetry and C5 purity.
    return out, data


def _end_of_eft1(z: dict, *, name: str, ms: float) -> np.ndarray:
    """Numerically select MS endpoint of a genuine descending RGBeta trajectory."""
    key = 'eft1__' + name
    if key not in z:
        raise AdapterBlocked(f'B missing {name} at MS, neither canonical nor eft1__ trajectory')
    mu = _numerical_array(z.get('mu_eft1_gev'),'B/mu_eft1_gev',real=True)
    if mu.ndim != 1 or len(mu)<2 or not np.all(np.diff(mu)<0) or not np.isclose(mu[-1],ms,rtol=1e-10,atol=1e-8):
        raise AdapterBlocked('B trajectory must descend and terminate at the common MS')
    traj = _numerical_array(z[key],key)
    if traj.shape[0] != len(mu):
        raise AdapterBlocked(f'B/{key}: time axis differs from EFT1 scales')
    return traj[-1]


def _sm_from_b(meta: dict, z: dict, ms: float) -> tuple[dict, dict]:
    sm = meta.get('sm_couplings_at_ms')
    if isinstance(sm, dict) and all(k in sm for k in ('gY','g2','g3','lambdaH')):
        constants = {}
        for k in ('gY','g2','g3','lambdaH'):
            val = _numerical_array(sm[k],f'B/{k}',shape=(),real=True)
            constants[k] = float(val)
    else:
        constants={}
        for k in ('gY','g2','g3','lambdaH'):
            v = _numerical_array(_end_of_eft1(z,name=k,ms=ms),f'B/{k} at MS',shape=(),real=True)
            constants[k]=float(v)
    matrices={}
    for k in ('ye','yu','yd'):
        src=z.get(k)
        if src is None:
            src=_end_of_eft1(z,name=k,ms=ms)
        matrices[k]=_numerical_array(src,'B/'+k,shape=(3,3))
    return constants, matrices


def _permute_full_scalar_axes(meta:dict, n:int, z:dict) -> tuple[np.ndarray, np.ndarray | None]:
    """Require an actual Higgs-first scalar component map, not positional renaming."""
    q=_numerical_array(z.get('quartic_real_full_ms'),'B/full scalar quartic',shape=(4+n,)*4,real=True)
    channels = z.get('quartic_real_channels_ms')
    if channels is not None:
        channels=_numerical_array(channels,'B/quartic channels',real=True)
        if channels.ndim != 5 or channels.shape[1:] != (4+n,)*4:
            raise AdapterBlocked('B full quartic channels have incorrect shape')
    if meta.get('format')==B_FORMAT:
        return q,channels
    full_labels=meta.get('scalar_full_real_labels')
    if isinstance(full_labels,list) and len(full_labels)==n+4:
        targets=['H[1,R]','H[1,I]','H[2,R]','H[2,I]']
        try:
            hidx=[full_labels.index(x) for x in targets]
        except ValueError as exc:
            raise AdapterBlocked('B full quartic labels lack standard Higgs R/I components') from exc
        order=hidx+[i for i in range(n+4) if i not in hidx]
        if len(set(full_labels))!=n+4:
            raise AdapterBlocked('B full scalar labels are duplicated')
        for ax in range(4):q=np.take(q,order,axis=ax)
        if channels is not None:
            for ax in range(1,5):channels=np.take(channels,order,axis=ax)
        # The heavy basis of the mass matrix must be in the same relative order.
        listed=[full_labels[i] for i in order[4:]]
        if listed != meta.get('scalar_background_real_labels'):
            raise AdapterBlocked('B full quartic and mass/background heavy-scalar label ordering mismatch')
    elif meta.get('quartic_scalar_order') != 'H_RI_interleaved_then_heavy_background_v1':
        raise AdapterBlocked('B full real quartic lacks explicit Higgs-first axis convention')
    return q,channels


def _source_background(meta: dict, z: dict, model: str, identity: dict) -> tuple[dict,dict]:
    fmt=meta.get('format',meta.get('schema'))
    if fmt not in (B_FORMAT,'t3-parallel-B-background-v1'):
        raise AdapterBlocked(f'Unrecognised B format {fmt!r}')
    if meta.get('quartic_potential_convention') != Q_CONVENTION or meta.get('quartic_inventory_status') != 'full_tensor':
        raise AdapterBlocked('B full Matchete/RGBeta-to-V4/4! quartic tensor is not declared complete')
    n=n_real(model)
    q,channels=_permute_full_scalar_axes(meta,n,z)
    m=z.get('mass_squared_gev2_ms')
    if m is None:
        # Native B's mS1Sq/mS2Sq are NOT silently turned into a mass Hessian.
        # This path is allowed only when the owner explicitly attests to the
        # source quadratic convention and the canonical mass frame.
        if meta.get('mass_hessian_constructor') != 'explicit_t3_native_quadratic_v1':
            raise AdapterBlocked('B lacks full scalar mass Hessian; native mS* beta variables are not sufficient without an attested mass-basis conversion')
        if meta.get('background_scalar_basis_id') != model+':canonical_interleaved':
            raise AdapterBlocked('Native quadratic constructor requires the canonical-interleaved mass frame')
        from .models import from_parameters
        s1=complex(_end_of_eft1(z,name='mS1Sq',ms=identity['mu_s_gev']))
        s2=complex(_end_of_eft1(z,name='mS2Sq',ms=identity['mu_s_gev']))
        if max(abs(s1.imag),abs(s2.imag))>1e-10:
            raise AdapterBlocked('Real scalar mass-square parameters must be real')
        mu=(complex(_end_of_eft1(z,name='muT3',ms=identity['mu_s_gev']))
            if model in ('B_-1','C_-1','A_-2','A_+0','D_-2','D_+0','E_-2','E_+0') else None)
        m=from_parameters(model,s1.real,s2.real,mu_sq=mu).mass2_canonical_gev2
    elif np.asarray(m).shape == (4+n,4+n):
        if not isinstance(meta.get('scalar_full_real_labels'),list):
            raise AdapterBlocked('Full H+S Hessian requires explicit full scalar ordering')
        original=meta['scalar_full_real_labels']
        hidx=[original.index(x) for x in ('H[1,R]','H[1,I]','H[2,R]','H[2,I]')]
        heavy=[i for i in range(n+4) if i not in hidx]
        m=_numerical_array(m,'B/H+S Hessian',real=True)
        if np.max(abs(m[np.ix_(hidx,heavy)]))>1e-9*max(1.,np.max(abs(m))):
            raise AdapterBlocked('B Hessian mixes Z2-even Higgs with Z2-odd heavy scalars')
        m=m[np.ix_(heavy,heavy)]
    m=_numerical_array(m,'B/heavy scalar Hessian',shape=(n,n),real=True)
    R=_numerical_array(z.get('R_wilson_from_background'), 'B/Wilson-to-background rotation',shape=(n,n),real=True)
    if not np.allclose(R@R.T,np.eye(n),rtol=1e-10,atol=1e-10):
        raise AdapterBlocked('B claimed scalar rotation is not orthogonal')
    b_labels=meta.get('scalar_background_real_labels')
    if not isinstance(b_labels,list) or len(b_labels)!=n or len(set(b_labels))!=n:
        raise AdapterBlocked('B requires complete ordered native heavy real-field labels')
    if not isinstance(meta.get('background_scalar_basis_id'),str) or not meta['background_scalar_basis_id']:
        raise AdapterBlocked('B background mass/quartic basis identifier absent')
    if not isinstance(meta.get('wilson_scalar_basis_id'),str) or not meta['wilson_scalar_basis_id']:
        raise AdapterBlocked('B rotation source basis identifier absent')
    if meta.get('single_common_scalar_threshold') is not True:
        raise AdapterBlocked('B must explicitly integrate BOTH scalars at ONE MS')
    if meta.get('scalar_mass_origin') not in ('unbroken_electroweak','EWSB_resummed','synthetic_algebra_only'):
        raise AdapterBlocked('B scalar mass origin not identified')
    sm, mats=_sm_from_b(meta,z,identity['mu_s_gev'])
    data={'mass_squared_gev2_ms':m, 'R_wilson_from_background':R,
          'quartic_real_full_ms':q, **mats}
    if channels is not None:data['quartic_real_channels_ms']=channels
    if 'scalar_gauge_generators_full' in z:
        data['scalar_gauge_generators_full']=z['scalar_gauge_generators_full']
    out={**identity,'format':B_FORMAT,'tensor_npz':'background.npz',
         'quartic_potential_convention':Q_CONVENTION,
         'quartic_inventory_status':'full_tensor',
         'wilson_scalar_basis_id':meta['wilson_scalar_basis_id'],
         'background_scalar_basis_id':meta['background_scalar_basis_id'],
         'scalar_background_real_labels':b_labels,
         'single_common_scalar_threshold':True,
         'scalar_mass_origin':meta['scalar_mass_origin'],
         'sm_couplings_at_ms':sm,
         'quartic_channel_names':meta.get('quartic_channel_names', []),
         'scalar_operator_inventory_certified':meta.get('scalar_operator_inventory_certified',False),
         'native_matchete_rgbeta_normalization_certified':meta.get('native_matchete_rgbeta_normalization_certified',False),
         'source_B_format':fmt}
    return out,data


def _save_pair(path: Path, meta: dict, arrays: dict) -> Path:
    path.parent.mkdir(parents=True,exist_ok=True)
    npz=path.with_name(meta['tensor_npz'])
    np.savez_compressed(npz,**arrays)
    meta['tensor_npz_sha256']=_sha(npz)
    path.write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return path


def prepare_pair(model: str, wilson_json: Path, background_json: Path,
                 destination: Path, *, allow_synthetic: bool=False) -> PreparedPair:
    """Adapt both producers without regenerating any fermion or scalar EFT running."""
    model=normalize_model(model)
    wm,wz,wh=_read_input(Path(wilson_json),nominal_npz='wilson_at_MS.npz')
    bm,bz,bh=_read_input(Path(background_json),nominal_npz='background.npz')
    wi=_identity_fields(wm,model,label='A')
    bi=_identity_fields(bm,model,label='B')
    _same_input(wi,bi)
    if wi['synthetic_fixture'] and not allow_synthetic:
        raise AdapterBlocked('Synthetic input forbidden for RoughNumericalBenchmark')
    wout, wa = _source_wilson(wm,wz,model,wi)
    bout, ba = _source_background(bm,bz,model,bi)
    if wout['scalar_basis_id'] != bout['wilson_scalar_basis_id']:
        raise AdapterBlocked('B rotation does not identify the A real scalar source basis')
    if wout.get('scalar_source_real_labels') is None or len(wout['scalar_source_real_labels'])!=n_real(model):
        raise AdapterBlocked('A source scalar labels missing or incomplete')
    for key in ('matching_context_id','matched_seed_id','lagrangian_source_id','mu_s_gev','mu_f_gev'):
        if wout.get(key) != bout.get(key):
            raise AdapterBlocked(f'Post-conversion A/B provenance mismatch {key}')
    outdir=Path(destination)
    a=_save_pair(outdir/'wilson.json',wout,wa)
    b=_save_pair(outdir/'background.json',bout,ba)
    prov={'A':wh,'B':bh,'C_input_A_manifest':str(a),'C_input_B_manifest':str(b),
          'adapted_A_format':wm.get('format'),
          'adapted_B_format':bm.get('format',bm.get('schema')),
          'synthetic_fixture':wi['synthetic_fixture'],
          'no_native_coupling_guess':True,
          'no_recomputed_EFT1_log':True}
    return PreparedPair(a,b,prov)
