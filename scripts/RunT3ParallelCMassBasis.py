"""Export C-owned scalar quadratic mass-basis interface for Workstream B."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.parallel_t3.threshold.models import from_parameters,normalize_model,n_real,scalar_labels
from Numerical.parallel_t3.threshold.matching import mass_eigenbasis


def main(argv=None):
    p=argparse.ArgumentParser(description='Source-convention scalar mass Hessians and native field maps')
    p.add_argument('--model',required=True)
    p.add_argument('--mS1sq',type=float)
    p.add_argument('--mS2sq',type=float)
    p.add_argument('--mu-real',type=float)
    p.add_argument('--mu-imag',type=float)
    p.add_argument('--background-npz',type=Path,
                  help='Alternative: diagonalize the B-supplied arbitrary real mass matrix')
    p.add_argument('--mass-key',default='mass_squared_gev2_ms')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(argv)
    model=normalize_model(a.model)
    if a.background_npz:
        with np.load(a.background_npz,allow_pickle=False) as data:
            M=np.array(data[a.mass_key])
        n=n_real(model)
        if M.shape!=(n,n):raise ValueError('Background real scalar mass dimension invalid')
        R=np.eye(n)
        extra={'origin':'explicit_B_mass_matrix',
               'gauge_invariance_certified':False,'native_rotation_unsupplied':True}
        native_labels=scalar_labels(model)
        kind='explicit_full_matrix'
    else:
        if a.mS1sq is None or a.mS2sq is None:
            p.error('Supply both mS1sq and mS2sq for analytic native mass construction')
        mu=None
        if a.mu_real is not None or a.mu_imag is not None:
            mu=complex(a.mu_real or 0.,a.mu_imag or 0.)
        ans=from_parameters(model,a.mS1sq,a.mS2sq,mu)
        M=ans.mass2_canonical_gev2
        R=ans.R_canonical_from_native
        extra=ans.extra
        native_labels=list(ans.native_labels)
        kind=ans.kind
    eig,U=mass_eigenbasis(M)
    a.output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(a.output/'scalar_mass_basis.npz',
                  mass_squared_canonical_gev2=M,
                  R_canonical_from_native=R,
                  mass_squared_native_gev2=R.T@M@R,
                  U_canonical_from_eigen=U,
                  eigenvalues_sq_gev2=eig)
    report={'format':'t3-parallel-c-scalar-mass-export-v1','model':model,
            'kind':kind,'original_basis':'canonical_interleaved_S1_then_S2',
            'canonical_labels':scalar_labels(model),
            'native_labels':native_labels,
            'field_map':'phi_canonical = R_canonical_from_native phi_native',
            'diagonalisation':'phi_canonical=U_canonical_from_eigen phi_mass',
            'all_scalars_remain_active_until_one_MS':True,
            'masses2_gev2':eig.tolist(),'masses_gev':np.sqrt(eig).tolist(),
            'extra':extra,
            'native_CG_numerical_alignment_certified':False,
            'source':'Inspected original public source at GitHub main a944e003 (read-only)'}
    (a.output/'scalar_mass_basis.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(a.output/'scalar_mass_basis.json')
    return 0

if __name__=='__main__':sys.exit(main())
