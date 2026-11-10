#!/usr/bin/env python3
"""A(-4) three-flavour LLSS -> full 45x45x12x12 EFT1 Wilson RG bridge.

Uses source-exact Matchete T3MixCG and previously solved Batch 1 RG data.
No Mathematica calls or UV re-matching. Operator and spinor normalization
certificates remain missing, so outputs are FORMAL and not physical C5.

Fixed-order mode computes all seven beta terms at mu_F, times ln(mu_S/mu_F).
Full-ode mode integrates the SAME master tensor with the previously solved,
running renormalisable EFT1 couplings over the entire F -> S interval.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.running.T3A4NativeWilsonBridge import (
    NativeWilsonBridgeError, derive_triplet_map, formal_tree_master,
    source_validation, three_flavour_wbar,
)

A4 = ROOT/'output/parallel_t3/full_general_final_one_batch/thesis_a4_fresh_f_v3/A_-4'
B1 = ROOT/'output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1/A_-4'
OUT = ROOT/'output/parallel_t3/a4_native_wilson_rg_v1/A_-4'
LOOP = 1./(16*math.pi**2)


def read_json(p: Path) -> dict:
    ans = json.loads(p.read_text(encoding='utf-8-sig'))
    if not isinstance(ans, dict):
        raise NativeWilsonBridgeError('Expected JSON object: '+str(p))
    return ans


def _source_file(path: Path, records: dict, label: str) -> Path:
    if not path.is_file():
        raise NativeWilsonBridgeError('Missing verified source: '+str(path))
    try:
        source_label = str(path.resolve().relative_to(ROOT))
    except ValueError:
        source_label = str(path.resolve())
    records[label] = {'project_path':source_label,
                      'sha256':sha256(path.read_bytes()).hexdigest()}
    return path


def build_running_context(batch1: Path, mu_f: float, mu_s: float):
    """Reuse Batch 1 solved EFT1 RG trajectory, never rerun RGBeta.

    Produces full Yukawa, scalar-quartic and gauge-generator tensors in the
    existing Aminus4EFT1Couplings *reference RGBeta* coordinate basis.
    Matchete-to-RGBeta scalar quartic coupling factors are not certified.
    """
    from Numerical.parallel_t3.background.general_batch1_api import NativeBatch1Trajectory
    from Numerical.running.Aminus4EFT1Couplings import (
        BETA_NAMES, F_LABELS, S_LABELS, gauge_generators, tensor_point,
    )
    from Numerical.running.T3EFT1WilsonODETransport import Aminus4Trajectory
    native = NativeBatch1Trajectory.load(batch1)
    if native.manifest.get('model') != 'A_-4':
        raise NativeWilsonBridgeError('Wrong native Batch 1 model')
    names = tuple(native.manifest['coupling_order']['eft1'])
    if set(names) != BETA_NAMES:
        raise NativeWilsonBridgeError('Expected all 17 A(-4) running EFT1 couplings')
    with np.load(batch1/'trajectories.npz',allow_pickle=False) as data:
        mu = np.asarray(data['eft1__mu_gev'],float)
        points = [tensor_point(float(t),
                  {key:data['eft1__'+key][i] for key in names})
                  for i,t in enumerate(mu)]
    if len(mu)<3 or not np.isfinite(mu).all():
        raise NativeWilsonBridgeError('Missing full RG trajectory')
    if not np.isclose(mu[0],mu_f,rtol=1e-9) or not np.isclose(mu[-1],mu_s,rtol=1e-9):
        raise NativeWilsonBridgeError('Batch 1 F/S scale endpoints inconsistent')
    order = np.argsort(mu)
    gens = gauge_generators()
    gauge = {kind: tuple(zip(gens[kind]['fermion'],gens[kind]['scalar']))
             for kind in ('SU2','U1','SU3')}
    trajectory = Aminus4Trajectory(
        np.log(mu[order]), mu[order],
        np.stack([points[i].Y for i in order]),
        np.stack([points[i].lambda4 for i in order]),
        F_LABELS, S_LABELS,
        tuple({k:float(points[i].couplings[k]) for k in ('g2','gY','g3')}
              for i in order), gauge)
    return trajectory


def _c5_json(c: np.ndarray) -> dict:
    return {'real':np.asarray(c).real.tolist(),
            'imag':np.asarray(c).imag.tolist()}


def _norm(c: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(c)))


def run(a4: Path, b1: Path, *, mode: str = 'fixed-order',
        rtol: float = 3e-6, max_step_log: float = 1.0
        ) -> tuple[dict, np.ndarray]:
    if mode not in ('fixed-order','full-ode'):
        raise NativeWilsonBridgeError('Unknown Wilson transport mode')
    inputs = {}
    seed = read_json(_source_file(a4/'data/eft1_after_F_wilson_seed.json',
                                  inputs, 'native_CG_and_tree_Wilson_seed'))
    text = _source_file(a4/'A_EFT1_tree_contracted.txt',
                        inputs, 'fresh_EFT1_tree').read_text(encoding='utf-8-sig')
    status = read_json(_source_file(a4/'native_two_threshold_status.json',
                                    inputs,'fresh_native_status'))
    lag = read_json(_source_file(a4/'A_joint_scalar_lagrangian_audit.json',
                                 inputs,'native_hermiticity'))
    first = read_json(_source_file(b1/'fermion_matching_input.json',
                                   inputs,'batch1_fermion_endpoint'))
    _source_file(b1/'manifest.json',inputs,'batch1_rg_manifest')
    _source_file(b1/'trajectories.npz',inputs,'batch1_existing_rg_trajectory')
    source_validation(seed,text,status,lag)
    U = derive_triplet_map(seed)
    Wbar, mass = three_flavour_wbar(first)
    # Unit FORMAL GammaCC factor, explicitly NOT a certified convention.
    tree = formal_tree_master(Wbar,U)
    mf = float(first['scale_GeV'])
    ms = 1e3
    tr = build_running_context(b1,mf,ms)
    from RGE.running.intermediate.GeneralEFT1C5Mixing import (
        calculate_master_beta, project_c5, TERMS,
    )
    at_f = tr.context(mf,mf,ms)
    rg_terms = calculate_master_beta(at_f,tree)
    if set(rg_terms) != set(TERMS) or any(not np.isfinite(x).all()
                                          for x in rg_terms.values()):
        raise NativeWilsonBridgeError('Incomplete/nonfinite seven-term master beta')
    projected = {}
    for name in TERMS:
        mat, _, _ = project_c5(rg_terms[name],at_f)
        projected[name] = mat
    total_beta = sum(rg_terms.values())
    beta_c5, beta_ward, beta_flavour = project_c5(total_beta,at_f)
    if not np.allclose(beta_c5,sum(projected.values()),rtol=1e-12,atol=1e-26):
        raise NativeWilsonBridgeError('Seven-term beta projection does not add up')
    if mode == 'fixed-order':
        log_interval = math.log(ms/mf)
        final = tree + LOOP*log_interval*total_beta
        ode = {'method':'fixed-order-one-loop-leading-log',
               'running_renormalizable_couplings_inside_beta':False,
               'fixed_order_higher_order_terms_omitted':True,
               'ln_muS_over_muF':log_interval}
    else:
        from Numerical.running.T3EFT1WilsonODETransport import integrate_master_tensor
        final, ode = integrate_master_tensor(
            tree,lambda mu:tr.context(mu,mf,ms),mf,ms,
            rtol=rtol,atol=1e-20,max_step_log=max_step_log)
        ode['method']='full-master-numerical-ODE'
        ode['running_renormalizable_couplings_inside_beta']=True
        ode['resums_higher_order_Wilson_leading_logs']=True
    ctx = tr.context(ms,mf,ms)
    c5, ward, sym = project_c5(final,ctx)
    zero_tree_c5, _, _ = project_c5(tree,at_f)
    if np.max(np.abs(zero_tree_c5)) > 1e-27:
        raise NativeWilsonBridgeError('Tree LLSS unexpectedly contains LLHH')
    report = {
        'format':'t3-a4-source-checked-tree-master-Wilson-rg-v1',
        'status':'FormalTreeMasterRGComputed_NotPhysicalC5',
        'model':'A_-4',
        'integrated_heavy_fields_in_order':[['F'],['S1','S2']],
        'renormalizable_eft1_rg_reused':True,
        'fermion_flavour_count':3,
        'master_tensor_dimensions':list(tree.shape),
        'source_files':inputs,
        'cg_conversion':{
            'source_CG_shape':[2,2,3],
            'spherical_normalization_factor':U.spherical_cg_scale,
            'cartesian_from_spherical_real':U.cartesian_from_spherical.real.tolist(),
            'cartesian_from_spherical_imag':U.cartesian_from_spherical.imag.tolist(),
            'real_triplet_map':U.real_s2_cartesian_from_spherical.tolist(),
            'CG_reconstruction_residual':U.max_reconstruction_residual,
            'triplet_unitarity_residual':U.max_unitarity_residual,
            'real_scalar_orthogonality_residual':U.max_real_orthogonality_residual,
            'exact_source_CG_to_reference_spherical_coordinate_map':True,
            'native_Matchete_to_RGBeta_full_quartic_dictionary_certified':False,
        },
        'heavy_mass_and_F_tree_flavour_audit':mass,
        'wilson_solver':ode,
        'master_tree_Wilson_frobenius_GeVinv':_norm(tree),
        'master_at_muS_frobenius_GeVinv':_norm(final),
        'per_seven_master_beta_channel_C5_Frobenius':{
            k:_norm(v) for k,v in projected.items()},
        'beta_C5_muF_16pi2_formal_GeVinv':_c5_json(beta_c5),
        'beta_C5_muF_ward_absolute_residual':float(beta_ward),
        'beta_C5_muF_flavour_symmetry_residual':float(beta_flavour),
        'C5_LLHH_muS_formal_GeVinv':_c5_json(c5),
        'C5_LLHH_muS_formal_frobenius_GeVinv':_norm(c5),
        'C5_LLHH_muS_formal_ward_absolute_residual':float(ward),
        'C5_LLHH_muS_formal_flavour_symmetry_residual':float(sym),
        'gamma_cc_to_weyl_factor_and_sign_certified':False,
        'Matchete_to_RGBeta_coupling_dictionary_certified':False,
        'all_F_1loop_counterterm_operator_sources_integrated':False,
        'full_eft1_dimension_five_basis_closure_certified':False,
        'finite_MSbar_matching_verified':False,
        'joint_scalar_matching_inherited':False,
        'physical_C5_computed':False,
        'physical_neutrino_masses_computed':False,
        'warning':('The CG component transformation is algebraic, but this is '
                   'only a formal master-tensor RG source. The unknown signed '
                   'GammaCC/Weyl factor and independent scalar-coupling basis '
                   'map forbid physical Wilson/C5 interpretation.'),
    }
    return report, final


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--a4-dir',type=Path,default=A4)
    p.add_argument('--batch1-dir',type=Path,default=B1)
    p.add_argument('--out',type=Path,default=OUT)
    p.add_argument('--mode',choices=('fixed-order','full-ode'),default='fixed-order')
    p.add_argument('--rtol',type=float,default=3e-6)
    p.add_argument('--max-step-log',type=float,default=1.)
    args=p.parse_args(argv)
    try:
        if args.out.exists():
            raise NativeWilsonBridgeError('Refusing overwrite: '+str(args.out))
        report, final=run(args.a4_dir,args.batch1_dir,mode=args.mode,
                          rtol=args.rtol,max_step_log=args.max_step_log)
        args.out.mkdir(parents=True,exist_ok=False)
        path=args.out/'A4_formal_tree_master_RG.json'
        path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',
                        encoding='utf-8')
        np.savez_compressed(args.out/'A4_formal_master_tensor_muS.npz',
                            C_LLphiPhi_muS=final)
        print('A(-4) exact native CG -> RGBeta triplet coordinates: VERIFIED')
        print('A(-4) one-loop master tensor transport:',args.mode,'COMPUTED')
        print('Formal LLHH Frobenius [GeV^-1]:',
              f"{report['C5_LLHH_muS_formal_frobenius_GeVinv']:.9e}")
        print('Physical C5 / full MSbar matching: NOT COMPUTED')
        print('Output:',path)
        return 0
    except (NativeWilsonBridgeError, ValueError, KeyError, FileNotFoundError,
            OSError, RuntimeError) as e:
        print('BLOCKED:',e,file=sys.stderr)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
