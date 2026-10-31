"""Independent two-component Majorana flavour-action consistency witness.

Tests the coefficient transformation under chi_old=U chi_new, with
mass M_new=U.T M_old U and y_new=y_old U. The conjugate chiral
component transforms with U*, but CConj(F)=F for a Majorana field.
This does NOT certify Matchete/RGBeta coupling normalization or C5.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

MODELS = ('B_-1', 'C_-1')

def rel(a,b):
    return float(np.linalg.norm(np.asarray(a)-np.asarray(b))/max(np.linalg.norm(np.asarray(b)),1e-30))

def unitary(rng):
    a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    q,r=np.linalg.qr(a)
    return q @ np.diag(np.conj(np.diag(r))/np.abs(np.diag(r)))

def one_trial(seed):
    rng=np.random.default_rng(seed)
    U=unitary(rng)
    chi=rng.normal(size=3)+1j*rng.normal(size=3)
    ell=rng.normal(size=3)+1j*rng.normal(size=3)
    y1=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    y2=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    A=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    M=A+A.T+np.diag([5.,7.,9.])
    old=U@chi
    mass_old=old.T@M@old
    mass_new=chi.T@(U.T@M@U)@chi
    # bilinear coefficient tests (light and scalar factors are spectators)
    a1=ell.T@y1@old; a2=ell.T@y2@old
    b1=ell.T@(y1@U)@chi; b2=ell.T@(y2@U)@chi
    wrong1=ell.T@(y1@U.conj())@chi
    # For a 4-component Majorana field (chi, epsilon chi*) both chiralities
    # rotate differently, while CConj(F) = F and its left component is chi.
    conjugate_old=U.conj()@chi.conj()
    conjugate_new=(U@chi).conj()
    chiral_wrong=rel(y1@U.conj(),y1@U)
    return {'unitarity':rel(U.conj().T@U,np.eye(3)),
            'majorana_mass':rel(mass_old,mass_new),
            'y1_majorana_left':rel(a1,b1),
            'y2_majorana_left':rel(a2,b2),
            'right_chiral_conjugate':rel(conjugate_old,conjugate_new),
            'wrong_y1_conjugate_U_action':rel(wrong1,a1),
            'wrong_y1_matrix_relative':chiral_wrong}

def run(trials=64,base_seed=2311,native=None):
    if trials<2: raise ValueError('at least two complex trials required')
    records=[one_trial(base_seed+i) for i in range(trials)]
    keys=('unitarity','majorana_mass','y1_majorana_left','y2_majorana_left','right_chiral_conjugate')
    worst={k:max(x[k] for x in records) for k in keys}
    if any(v>2e-12 for v in worst.values()):raise AssertionError(f'action covariance failed {worst}')
    controls={k:min(x[k] for x in records) for k in ('wrong_y1_conjugate_U_action','wrong_y1_matrix_relative')}
    if controls['wrong_y1_matrix_relative']<1e-3:raise AssertionError('negative control not discriminating')
    n={'status':'NativeEvidenceNotSupplied'}
    if native is not None:
        p=Path(native)
        if not p.is_file():raise FileNotFoundError(str(p))
        data=json.loads(p.read_text(encoding='utf-8-sig'))
        models=data.get('models',{})
        statuses={k:models.get(k.replace('-','m'),{}).get('passed',False) for k in MODELS}
        # B_m1 / C_m1 in uploaded native JSON
        if not all(statuses.values()):raise ValueError('native report lacks B/C passed checks')
        checks={k:models[k.replace('-','m')]['checks'] for k in MODELS}
        if not all(all(isinstance(v,bool) and v for v in row.values()) and len(row)>=6 for row in checks.values()):
            raise ValueError('native checks incomplete')
        n={'status':'NativeAbstractChiralProbePassed','models':statuses,'checks':checks,
           'concrete_CC_orientation_measured':bool(data.get('concrete_CC_orientation_measured',False)),
           'cross_package_certified':bool(data.get('rgbeta_cross_package_dictionary_certified',False))}
        if n['cross_package_certified']:raise ValueError('unexpected native certification')
    return {'format':'t3-bc-majorana-action-flavour-batch-v1',
            'status':'MajoranaWeylActionCovariancePassedNativeDictionaryPending',
            'n_trials':trials,'seed_start':base_seed,'max_covariance_residuals':worst,
            'negative_control_min_residuals':controls,
            'native':n,
            'verified_local_identity':'For self-conjugate F, chi_old=U chi_new implies y1_new=y1_old U and y2_new=y2_old U in a common left-Weyl action convention',
            'matchete_rgbeta_absolute_dictionary_certified':False,
            'production_yukawa_map_modified':False,
            'physical_C5_certified':False,
            'next_gate':'Derive native Matchete-to-RGBeta indexed coupling/CG/action normalization before editing running one-scale triangle'}
