#!/usr/bin/env python3
"""B/C fixed-scale, leading one-loop masses from existing native Matchete A.

No production-module changes. Requires archived Matchete InputForm sources.
Does not certify UV-to-100-GeV RG evolution or an independent UV amplitude.
"""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import numpy as np


def as_complex(z):
    if isinstance(z, (int, float, complex, np.number)):
        return complex(z)
    if isinstance(z, (list, tuple)) and len(z) == 2:
        return complex(float(z[0]), float(z[1]))
    raise ValueError('Complex values must be real numbers or [real, imag]')


def matrix(z, shape):
    if isinstance(z, dict) and 'real' in z and 'imag' in z:
        z = (np.asarray(z['real'],dtype=float)+1j*np.asarray(z['imag'],dtype=float)).tolist()
    a = np.asarray([[as_complex(v) for v in row] for row in z], dtype=complex)
    if a.shape != shape or not np.isfinite(a).all():
        raise ValueError(f'Expected finite matrix of shape {shape}, received {a.shape}')
    return a


def pairs(a):
    return [[[float(z.real), float(z.imag)] for z in row] for row in a]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--input', type=Path, help='Explicit JSON with models B_m1 and C_m1')
    parser.add_argument('--archived-reference', action='store_true', help='Use the archived frozen v18 reference, not user physical inputs')
    parser.add_argument('--output', type=Path, default=Path('Reports/validation/bc_fixed_order_masses.json'))
    args = parser.parse_args()
    root = args.root.resolve()
    sys.path.insert(0, str(root))
    from importlib.util import spec_from_file_location, module_from_spec
    src = root / 'scripts/DirectVsSequentialOneLoop.experimental.py'
    if not src.is_file():
        raise FileNotFoundError('Missing existing Workstream C v22 source: '+str(src))
    spec = spec_from_file_location('t3_v22_direct', src)
    v22 = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(v22)
    if args.input and args.archived_reference:
        raise ValueError('Choose either --input or --archived-reference')
    if args.input:
        payload = json.loads(args.input.read_text(encoding='utf-8-sig'))
        configs = payload['models']
    elif args.archived_reference:
        v18_path = root / 'scripts/EvaluateT3BCThreeFermionTriangle.experimental.py'
        if not v18_path.is_file():
            raise FileNotFoundError('Missing v18 frozen reference helper: '+str(v18_path))
        v18spec = spec_from_file_location('t3_v18_ref', v18_path)
        v18 = module_from_spec(v18spec)
        assert v18spec.loader is not None
        v18spec.loader.exec_module(v18)
        configs = {}
        def enc(a):
            return [[float(z.real),float(z.imag)] for z in a] if np.asarray(a).ndim==1 else [[[float(z.real),float(z.imag)] for z in row] for row in a]
        for stem in v22.MODELS:
            params = v22.read_json(root/v22.BASE/f'{stem}_final_parameters.json')
            reference = v18.load_frozen_reference(root,stem)
            configs[stem] = {
                'MF_GeV':reference['mf'].tolist(),
                'scalar_eigenmasses_GeV':params['heavy_masses_gev'],
                'scalar_U':params['unitary_U'],
                'y1':enc(reference['y1']), 'y2':enc(reference['y2']),
                'quartics':{n:([float(reference['quartics'][n].real),float(reference['quartics'][n].imag)] if n in reference['quartics'] else 0) for n in v22.Q_NAMES},
                'vev_GeV':reference['vev_gev'],
            }
    else:
        raise ValueError('Supply --input for explicit same-scale inputs or --archived-reference for existing frozen benchmark')
    output = {'format':'t3-bc-native-direct-one-loop-fixed-scale-mass-v1',
              'source':str(src.relative_to(root)),
              'source_sha256':sha256(src.read_bytes()).hexdigest(),
              'scope':'Direct fixed-scale one-loop Matchete A, 3 diagonal positive Majorana masses, two mixed scalar masses',
              'full_rg_evolved_prediction_certified':False,
              'independent_non_matchete_uv_derivation_performed':False,
              'parameter_source':'archived_frozen_v18_reference' if args.archived_reference else 'user_explicit_common_scale_input',
              'models':{}}
    for stem, (model, _) in v22.MODELS.items():
        if stem not in configs:
            raise ValueError('Missing input model '+stem)
        cfg = configs[stem]
        F = np.asarray(cfg['MF_GeV'],dtype=float)
        scalar = np.asarray(cfg['scalar_eigenmasses_GeV'],dtype=float)
        if F.shape != (3,) or scalar.shape != (2,) or not np.isfinite(F).all() or not np.isfinite(scalar).all():
            raise ValueError('MF_GeV must have 3 finite masses; scalar_eigenmasses_GeV [MX,MY]')
        X,Y=map(float, scalar)
        if not (np.all(F>Y) and Y>X>0):
            raise ValueError('Require MF_i > MY > MX > 0 for v22 native direct evaluator')
        U=matrix(cfg['scalar_U'], (2,2))
        if np.linalg.norm(U.conj().T@U-np.eye(2)) > 1e-10:
            raise ValueError('Non-unitary scalar_U')
        y1=matrix(cfg['y1'],(3,3)); y2=matrix(cfg['y2'],(3,3))
        q={key:as_complex(cfg['quartics'][key]) for key in v22.Q_NAMES}
        vev=float(cfg['vev_GeV'])
        if not np.isfinite(vev) or vev <= 0:
            raise ValueError('vev_GeV must be positive and finite')
        direct_path=root/v22.BASE/'one_loop_probe'/f'{stem}_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt'
        report_path=root/v22.BASE/'one_loop_probe'/f'{stem}_final_one_loop_probe.json'
        if not direct_path.is_file():
            raise FileNotFoundError('Missing native amplitude '+str(direct_path))
        # Check original Matchete report and archived scope/topology gates.
        v22.source_gates(root,stem)
        compiled,denom,_,_=v22.compile_source(direct_path.read_text(encoding='utf-8-sig'), model)
        C5=np.zeros((3,3),complex)
        per_mass=[]
        for i,m in enumerate(F):
            a11,a12,a22=v22.direct_ordered(compiled,denom,float(m),X,Y,U,q)
            K=np.array([[2*a11,a12],[a12,2*a22]],complex)
            # Source conventions: barred Yukawa couplings in ordered A.
            Ycol=np.column_stack((y1[:,i].conj(),y2[:,i].conj()))
            C5+=Ycol@K@Ycol.T
            per_mass.append({'MF_GeV':float(m),'ordered_A_coefficients_GeV_inverse':[[float(z.real),float(z.imag)] for z in (a11,a12,a22)]})
        sym=float(np.linalg.norm(C5-C5.T)/max(np.linalg.norm(C5),1e-250))
        if sym>1e-12:
            raise ValueError(f'{model} symmetric C5 check failed: {sym}')
        mnu=-0.5*vev**2*C5*1e9
        masses=np.sort(np.linalg.svd(mnu,compute_uv=False))
        output['models'][stem]={
            'model':model,'status':'FIXED_SCALE_ONE_LOOP_COMPUTED_NOT_RG_CERTIFIED',
            'native_source_sha256':sha256(direct_path.read_bytes()).hexdigest(),
            'source_report_sha256':sha256(report_path.read_bytes()).hexdigest(),
            'C5_GeV_inverse_complex_pairs':pairs(C5),
            'mnu_eV_complex_pairs':pairs(mnu),
            'masses_eV_ascending':[float(z) for z in masses],
            'C5_transpose_symmetry_residual':sym,
            'MF_contributions':per_mass,
            'full_physical_rg_evolved_mass_certified':False,
        }
    out=args.output if args.output.is_absolute() else root/args.output
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(output,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    for key,v in output['models'].items():
        print(key, 'fixed-scale leading-one-loop masses [eV]:',v['masses_eV_ascending'])
    print('Saved:',out)
    print('NOTE: source is native Matchete; running high-scale-to-low-scale result is NOT certified')


if __name__=='__main__':
    main()
