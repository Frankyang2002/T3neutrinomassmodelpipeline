"""Read-only formal-source versus canonical LLSS Ward audit (A_-4 retained trials)."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical


def variation(c, f, s):
    return (np.einsum('ki,kjab->ijab',f,c,optimize=True)
            +np.einsum('kj,ikab->ijab',f,c,optimize=True)
            +np.einsum('ca,ijcb->ijab',s,c,optimize=True)
            +np.einsum('cb,ijac->ijab',s,c,optimize=True))


def relative(x, c):
    return float(np.linalg.norm(x)/max(np.linalg.norm(c),1e-100))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--capture',type=Path,default=Path('output/parallel_t3/all16_ward_capture01'))
    args=ap.parse_args()
    root=Path.cwd()
    reports=sorted(args.capture.glob('fit16_*_ward_diagnostic.json'))
    if not reports: raise FileNotFoundError('No trial diagnostic manifests')
    for report in reports:
        token=report.name.removesuffix('_ward_diagnostic.json')
        folder=root/'output/parallel_t3/rough_inputs'/token/'a/A_m4'
        case=load_canonical(root,folder,'A_-4',c5_boundary='require',assume_formal_identity=False)
        refs=case.source['file_sha256']
        input_path=Path(refs.get('canonical_input_npz',refs.get('boundary_npz'))['path'])
        with np.load(input_path,allow_pickle=False) as z:
            if 'C_LLSS_formal_muF' not in z or 'heavy_source_from_canonical' not in z:
                raise ValueError('Formal source and explicit rotation required')
            raw=z['C_LLSS_formal_muF'].astype(complex)
            R=z['heavy_source_from_canonical'].astype(float)
        basis=case.basis
        fidx=[i for pair in basis.leptons for i in pair]
        nh=basis.ns-4
        if raw.shape!=(3,3,2,2,nh,nh) or R.shape!=(nh,nh): raise ValueError('Unexpected source dimensions')
        orth=float(np.linalg.norm(R@R.T-np.eye(nh)))
        # Compare source and canonical LLSS up to the unknown overall signed spinor factor.
        src=raw.transpose(0,2,1,3,4,5).reshape(6,6,nh,nh)
        can=case.initial[np.ix_(fidx,fidx,range(4,basis.ns),range(4,basis.ns))]
        transformed=np.einsum('ijab,aA,bB->ijAB',src,R,R,optimize=True)
        factor=np.vdot(transformed,can)/max(np.vdot(transformed,transformed).real,1e-100)
        map_err=relative(can-factor*transformed,can)
        g=basis.generators()['SU2']
        results=[]
        for k,(ff,ss) in enumerate(g):
            lf=ff[np.ix_(fidx,fidx)]
            sc=ss[4:,4:]
            # phi_source = R phi_canonical; generators transform accordingly.
            sf=R@sc@R.T
            dw_source=variation(src,lf,sf)
            dw_canon=variation(can,lf,sc)
            transformed_ward=np.einsum('ijab,aA,bB->ijAB',dw_source,R,R,optimize=True)*factor
            equiv=relative(dw_canon-transformed_ward,dw_canon)
            results.append({'generator':f'SU2_T{k+1}',
                            'source_ward_relative':relative(dw_source,src),
                            'canonical_ward_relative':relative(dw_canon,can),
                            'source_to_canonical_ward_covariance_error':equiv})
        output={'status':'READ_ONLY_NOT_A_REPAIR','trial':token,'raw_source':str(input_path),
                'raw_shape':list(raw.shape),'rotation_orthogonality_error':orth,
                'canonical_vs_rotated_source_relative_error':map_err,
                'fitted_global_spinor_factor':{'real':float(factor.real),'imag':float(factor.imag)},
                'ward':results,
                'interpretation':('Source Ward uses conjugated heavy scalar generators R T R^T. '
                                  'Small map/covariance errors with order-one source Ward residual '
                                  'locate the failure in the formal tensor relative to the declared generators, '
                                  'but do not independently certify those generators or matching.')}
        path=args.capture/f'{token}_formal_source_map_ward.json'
        path.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
        print('\n',token,'rotation orthogonality:',orth,'source-to-canonical error:',map_err)
        print('fitted factor:',factor)
        for x in results: print(x)
        print('report:',path)

if __name__=='__main__': main()
