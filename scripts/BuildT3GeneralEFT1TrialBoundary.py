"""Build an A(-4) UNCERTIFIED general-EFT1 master-tensor trial input.

Intended path: scripts/BuildT3GeneralEFT1TrialBoundary.py

Reads the actual saved UV Yukawa benchmark, runs the native UV segment to
muF, recomputes the Matchete-factorised LLSS source and evaluates the native
HHSS scalar quartic tensor on the EFT1 muF boundary. Embeds both tensors into
GeneralEFT1C5Mixing's left-Weyl/real-scalar array conventions.

Trial scope: LLSS boundary + HHSS quartics + full SU(2)L/U(1)Y gauge
representations of the *six SM lepton-doublet Weyl components* and H,S1,S2.
SM singlet Weyl fermions, SM Yukawa vertices, non-HHSS scalar quartics and
other dimension-five operators are NOT included. Thus all seven algebraic
terms are evaluated, but this is NOT the complete general EFT1 RGE.
The Matchete GammaCC -> Weyl coefficient factor is provisionally 1, chosen
explicitly, not derived. No physical-certification claims are made.

PowerShell:
 python -X utf8 -m scripts.BuildT3GeneralEFT1TrialBoundary --model A_-4 \
   --source-run thesis16_merged_v1
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from itertools import permutations
from pathlib import Path

import numpy as np

from RGE.running.intermediate.GeneralEFT1C5Mixing import (
    FORMAT, CERTIFICATES, load_boundary,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/studies/general_scalar_extended/general_eft1_trial_boundaries'


def encode(z):
    z = complex(z)
    return {'real': float(z.real), 'imag': float(z.imag)}


def dense_components(tensor, nfirst, *, cutoff=0.0):
    """Store only one upper-triangle representative for each independent pair."""
    a = np.asarray(tensor, dtype=complex)
    assert a.ndim == 4 and a.shape[:2] == (nfirst, nfirst)
    out = []
    for i in range(nfirst):
        for j in range(i, nfirst):
            for k in range(a.shape[2]):
                for l in range(k, a.shape[3]):
                    x = a[i, j, k, l]
                    if abs(x) > cutoff:
                        out.append({'indices': [i+1, j+1, k+1, l+1],
                                    'value': encode(x)})
    return out


def quartic_hhss(q):
    """Full symmetric lambda on H(4)+S(n), with HHSS sector supplied."""
    q = np.asarray(q)
    if q.ndim != 4 or q.shape[:2] != (4, 4) or q.shape[2] != q.shape[3]:
        raise ValueError('Expected HHSS quartic [4,4,n,n]')
    if not np.isfinite(q).all() or np.max(np.abs(np.imag(q))) > 1e-12:
        raise ValueError('HHSS quartic must be real and finite')
    n = q.shape[2]
    lam = np.zeros((4+n,)*4, dtype=float)
    used = set()
    # Each of the 24 permutations belongs to one original H,H,S,S tuple.
    for A in range(4):
        for B in range(A, 4):
            for a in range(n):
                for b in range(a, n):
                    val = float(np.real(q[A, B, a, b]))
                    if val == 0:
                        continue
                    idx = (A, B, 4+a, 4+b)
                    for key in set(permutations(idx)):
                        if key in used and not np.isclose(lam[key], val, atol=1e-11, rtol=1e-9):
                            raise ValueError('HHSS channel does not form symmetric quartic')
                        lam[key] = val
                        used.add(key)
    rows = []
    for i in range(4+n):
        for j in range(i, 4+n):
            for k in range(j, 4+n):
                for l in range(k, 4+n):
                    val = lam[i, j, k, l]
                    if val:
                        rows.append({'indices':[i+1,j+1,k+1,l+1],'value':float(val)})
    return rows


def gauge_sector(state, n, *, d1=1, d2=3, y1=-2, y2=-1):
    """Exact H,S1,S2 real-field gauge generators from existing source conventions."""
    from RGE.general.GaugeGenerators import (
        su2_complex_generators, real_scalar_generator,
    )
    from scipy.linalg import block_diag
    import sympy as sp

    dims = (2, d1, d2)
    ys = (Fraction(1,2), Fraction(y1), Fraction(y2))
    if 4 + 2*(d1+d2) != 4+n:
        raise ValueError('Scalar basis dimensions disagree')
    g2, gy = float(state['g2']), float(state['gY'])
    if not np.isfinite([g2, gy]).all() or min(g2,gy) <= 0:
        raise ValueError('Invalid native EFT1 SU(2)/U(1) gauge couplings')
    sig = [np.asarray(x, complex) for x in su2_complex_generators(2)]
    result = []
    for name, g, ngen in (('SU2L',g2,3),('U1Y',gy,1)):
        gens = []
        for a in range(ngen):
            f = np.kron(np.eye(3), sig[a]) if name == 'SU2L' else -0.5*np.eye(6)
            if name == 'SU2L':
                blocks = [np.asarray(real_scalar_generator(su2_complex_generators(d)[a]),complex)
                          for d in dims]
            else:
                blocks = [np.asarray(real_scalar_generator(sp.Rational(y)*sp.eye(d)),complex)
                          for d,y in zip(dims,ys)]
            s = block_diag(*blocks)
            if s.shape != (n+4,n+4):
                raise ValueError('Incorrect embedded scalar generator dimension')
            gens.append({'fermion': {'real':f.real.tolist(),'imag':f.imag.tolist()},
                         'scalar': {'real':s.real.tolist(),'imag':s.imag.tolist()}})
        result.append({'name':name,'coupling':g,'generators':gens})
    return result


def build(args):
    if args.model != 'A_-4':
        raise ValueError('Only A_-4 is implemented; do not reuse the field basis for other models')
    from scripts import RunT3GeneralSequentialFitAblation as original
    from scripts import RunT3GeneralAll16FromReference as general
    from Numerical.orchestration.GeneralRunningBenchmark import initial_state, _state_values
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.running.GeneralT3DynamicFermionWilson import nonmixed_boundary
    from Numerical.running.GeneralT3SequentialLeadingLog import scalar_tensors
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest, load_source
    from Numerical.lo_certification.couplings.extra_hhss_geometry import model_extra_geometry

    saved = general.OUTPUT / args.source_run / original.tag(args.model)
    if not saved.is_dir():
        raise FileNotFoundError(f'Missing fitted benchmark: {saved}')
    settings = original.load(general.OUTPUT / args.source_run / 'settings.json')
    canonical_id = settings.get('canonical_run_id','full16_general_v1')
    betauv, _, _, _, _, _ = general.nonmixed_inputs(args.model, canonical_id,
                                                     generate_general_c5=False)
    overrides_path = saved / 'general_uv_overrides.json'
    overrides = original.load(overrides_path)
    _, start = initial_state(betauv, mu_uv=original.SCALES[0],
                             mu_f=original.SCALES[1], mu_s=original.SCALES[2],
                             overrides=overrides)
    sol = run_uv_segment(start, betauv, original.SCALES[1], rtol=1e-7, atol=1e-10)
    at_f = sol.final_state
    tree, tree_info = nonmixed_boundary(args.model, at_f, repo_root=original.ROOT,
                                        mu_f=original.SCALES[1])
    eft_f = build_intermediate_scalar_boundary(at_f, matching_scale_gev=original.SCALES[1])
    state = _state_values(eft_f)
    release = load_release_manifest(args.release_dir)
    ref = load_source(release[args.model], expected_model=args.model)
    channels = original._load_channels(args.model, ref, v4=args.v4_channel_dir,
                                       seeds=args.quartic_seed_dir)
    registry = original.load(args.registry)
    geometry = model_extra_geometry(args.model, registry['models'][args.model],
                                    include_components=True)
    snap = original.snapshot(args.model, state, saved/'fitted_config.json')
    tensor = scalar_tensors(args.model, tree, channels, snap, registry, geometry=geometry)
    raw = np.asarray(tree.by_source, dtype=complex)
    if raw.ndim != 7 or raw.shape[1:5] != (3,3,2,2):
        raise ValueError(f'Unexpected source-by-source tensor shape: {raw.shape}')
    n = raw.shape[-1]
    if raw.shape[-2] != n or n != 8 or tensor.quartic.shape != (4,4,n,n):
        raise ValueError('A_-4 requires 8 real heavy scalars and matching HHSS basis')
    if tuple(tree.labels) != tuple(['NewScalar1']*2+['NewScalar2']*6):
        raise ValueError('Heavy field ordering differs from intended A_-4 basis')
    scale = max(1e-40, float(np.max(np.abs(raw))))
    if not np.allclose(raw, raw.transpose(0,2,1,4,3,5,6), rtol=1e-9, atol=scale*1e-9):
        raise ValueError('LLSS source violates combined Weyl-component exchange symmetry')
    if not np.allclose(raw, raw.swapaxes(5,6), rtol=1e-9, atol=scale*1e-9):
        raise ValueError('LLSS source violates scalar-index exchange symmetry')

    operators = []
    factor = complex(args.source_factor)
    for index, src in zip(tree.source_indices, raw):
        c = np.zeros((6,6,4+n,4+n),complex)
        for p in range(3):
            for q in range(3):
                for i in range(2):
                    for j in range(2):
                        c[2*p+i,2*q+j,4:,4:] = factor * src[p,q,i,j]
        if not np.allclose(c,c.swapaxes(0,1),atol=scale*1e-9,rtol=1e-9):
            raise ValueError(f'F_tree source {index}: Weyl pair is not symmetric')
        if not np.allclose(c,c.swapaxes(2,3),atol=scale*1e-9,rtol=1e-9):
            raise ValueError(f'F_tree source {index}: scalar pair is not symmetric')
        entries = dense_components(c,6)
        if entries:
            operators.append({'source_id':f'F_tree_{index}', 'components':entries})
    if not operators:
        raise ValueError('No nonzero F-tree Wilson operators')
    payload = {
        'format':FORMAT, 'authority':'synthetic', 'model':args.model,
        'scheme':'MSbar', 'mu_F_gev':float(original.SCALES[1]),
        'mu_S_gev':float(original.SCALES[2]),
        'n_fermions':6, 'n_real_scalars':4+n,
        'lepton_components':[[1,2],[3,4],[5,6]],
        'higgs_real_components':[1,2,3,4],
        'scalar_field_blocks':[
            {'name':'H','role':'Higgs','indices':[1,2,3,4]},
            {'name':'S1','role':'heavy','indices':[5,6]},
            {'name':'S2','role':'heavy','indices':list(range(7,13))}],
        'master_lagrangian_convention':'C_ijab psi_i psi_j phi_a phi_b + hc',
        'weinberg_lagrangian_convention':'(1/2) C5_ij (L_i.H)(L_j.H) + hc',
        'certificates':{name:False for name in CERTIFICATES},
        'operators':operators,
        'quartic_components':quartic_hhss(tensor.quartic),
        # No invented SM Yukawa inputs. The truncated Weyl space contains
        # only the three lepton doublets, so all singlet Yukawas are omitted.
        'yukawa_components':[],
        'gauge_sectors':gauge_sector(state,n),
        'trial_provenance':{
            'status':'UNCERTIFIED_TRUNCATED_FIELD_CONTENT',
            'source_run':args.source_run, 'canonical_run_id':canonical_id,
            'muF_state_from_actual_uv_running':True,
            'source_boundary_sha256':tree.hash,
            'fitted_uv_overrides_sha256':sha256(overrides_path.read_bytes()).hexdigest(),
            'source_factor_assumed_real':args.source_factor,
            'fermion_source_diagnostics':tree_info,
            'included':'LLSS sources, native HHSS quartic directions, SU2/U1 lepton+scalar gauge generators',
            'omitted':'Charged-lepton and quark singlet Weyl states, SM Yukawas, HHHH and SSSS quartics, other dim-5 Wilson classes, C5(muF) hard matching',
            'not_a_complete_general_eft1_evaluation':True,
            'master_rge_gauge_covariance_not_certified':True,
            'tensor_convention_conversion_not_certified':True,
        },
    }
    # Real solver structural validation, with certification intentionally off.
    load_boundary(payload,allow_synthetic=True)
    target = args.output or (OUT/'A_m4_master_boundary.json')
    target = target.resolve()
    if target.exists() and not args.overwrite:
        raise FileExistsError(f'File exists (use --overwrite): {target}')
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Boundary:',target)
    print('F-tree sources:',len(operators),'HHSS components:',len(payload['quartic_components']))
    print('Gauge sectors: SU2L, U1Y; omitted SM Yukawa/singlet sector')
    print('UNCERTIFIED truncated seven-term trial only; not physical EFT1 running')


def main():
    from scripts import RunT3GeneralSequentialFitAblation as original
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model',default='A_-4')
    ap.add_argument('--source-run',default='thesis16_merged_v1')
    ap.add_argument('--source-factor',type=float,default=1.0,
                    help='Explicit assumed raw Matchete-to-Weyl multiplier; NOT certified')
    ap.add_argument('--release-dir',type=Path,default=original.DEFAULT_RELEASE)
    ap.add_argument('--v4-channel-dir',type=Path,default=original.DEFAULT_V4)
    ap.add_argument('--quartic-seed-dir',type=Path,default=original.DEFAULT_SEEDS)
    ap.add_argument('--registry',type=Path,default=original.DEFAULT_REGISTRY)
    ap.add_argument('--output',type=Path,default=None)
    ap.add_argument('--overwrite',action='store_true')
    build(ap.parse_args())
    return 0


if __name__=='__main__':
    raise SystemExit(main())
