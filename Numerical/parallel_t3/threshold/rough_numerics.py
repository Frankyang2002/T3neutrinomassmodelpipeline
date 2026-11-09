"""Rough Stream-3 finalisation; delegates matching/RGE to existing C pipeline.

Only low-energy Takagi/PMNS diagnostics and benchmark bookkeeping are new.
No hard one-loop term is manufactured; existing C explicitly assigns zero to
independently unavailable terms as a benchmark approximation.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import json
import math
import shutil
import numpy as np

from .pipeline import _load_pair, B_FORMAT, run_case
from .rough_adapter import PreparedPair

STATUS='RoughNumericalBenchmark'
OMITTED=(
    'Independent finite one-loop scalar threshold diagrams beyond full-tensor two-scalar bubble',
    'Finite scalar wavefunction and field-normalisation/EOM matching',
    'Independently unavailable fermion-hard and other dimension-five hard matching',
    'Two-loop and higher scalar/SMEFT corrections',
    'Pole-neutrino masses and electroweak matching corrections',
)


def _digest(path:Path)->str:
    return sha256(path.read_bytes()).hexdigest()


def _json(path:Path, obj:dict)->None:
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def _serialized_complex(x:np.ndarray)->dict:
    return {'real':np.asarray(x).real.tolist(), 'imag':np.asarray(x).imag.tolist()}


def neutrino_observables(m_gev:np.ndarray, ye_ms:np.ndarray) -> tuple[dict,dict]:
    """Sorted Takagi eigenmasses, and PMNS if mass and Ye bases are nondegenerate.

    SM one-loop beta(Ye) = a Ye + (3/2) Ye Ye† Ye. It preserves Ye's left
    singular vectors from MS to 173 GeV. Therefore Ye(MS) is enough to define
    the charged-lepton basis of an already-evolved m_nu(173 GeV).

    Sorting by absolute mass is NOT a physical normal/inverted relabelling.
    """
    m=np.asarray(m_gev,dtype=np.complex128)
    y=np.asarray(ye_ms,dtype=np.complex128)
    if m.shape!=(3,3) or y.shape!=(3,3) or not np.isfinite(m).all() or not np.isfinite(y).all():
        raise ValueError('Nonfinite 3x3 mass/Yukawa matrices')
    if not np.allclose(m,m.T,rtol=1e-9,atol=1e-18):
        raise ValueError('Majorana mass matrix must be symmetric, NOT Hermitian')
    ms=np.sort(np.linalg.svd(m,compute_uv=False))*1e9
    diffs=[float(ms[1]**2-ms[0]**2),float(ms[2]**2-ms[0]**2),
           float(ms[2]**2-ms[1]**2)]
    values={'mass_order':'ascending_Takagi_not_physical_ordering_assignment',
            'masses_eV':ms.tolist(), 'dm21_sq_eV2':diffs[0],
            'dm31_sq_eV2':diffs[1], 'dm32_sq_eV2':diffs[2],
            'sum_masses_eV':float(sum(ms)),
            'm_beta_beta_eV':float(abs(m[0,0])*1e9) if np.linalg.norm(y-np.diag(np.diag(y))) <= 1e-12*max(float(np.linalg.norm(y)),1e-40) else None,
            'm_beta_beta_note':'Only flavour m_ee if charged leptons are already diagonal; otherwise unreported',
            'mixing_defined':False,'mixing_blocker':None,
            'sin2_theta12':None, 'sin2_theta13':None,'sin2_theta23':None,
            'theta12_deg':None,'theta13_deg':None,'theta23_deg':None,
            'J_CP':None,'delta_CP_deg':None}
    matrices={}
    # Singular squared masses in GeV²; no precision-losing subtraction of
    # almost-zero entries in the mass matrix itself.
    e2,Ue=np.linalg.eigh(y@y.conj().T)
    n2,Un=np.linalg.eigh(m.conj().T@m)
    e=np.sqrt(np.maximum(e2,0))
    n=np.sqrt(np.maximum(n2,0))
    scale=max(float(np.max(n)),1e-60)
    escale=max(float(np.max(e)),1e-60)
    if (min(np.diff(n)) < 1e-8*scale or min(np.diff(e)) < 1e-8*escale):
        values['mixing_blocker']='Degenerate/near-degenerate neutrino or charged-lepton singular vectors'
        return values,matrices
    D=Un.T@m@Un
    residual=float(np.linalg.norm(D-np.diag(np.diag(D)))/max(np.linalg.norm(m),1e-100))
    if residual>1e-7:
        values['mixing_blocker']=f'Takagi numerical diagonalisation residual {residual:.3e}'
        return values,matrices
    for j in range(3):
        if abs(D[j,j])>0:
            Un[:,j]*=np.exp(-.5j*np.angle(D[j,j]))
    U=Ue.conj().T@Un
    s13=abs(U[0,2])**2
    c13=max(0.,1-s13)
    if c13<1e-12:
        values['mixing_blocker']='theta13=90 degrees, conventional theta12 undefined'
        return values,matrices
    s12=float(np.clip(abs(U[0,1])**2/c13,0,1))
    s23=float(np.clip(abs(U[1,2])**2/c13,0,1))
    s13=float(np.clip(s13,0,1))
    angles=[math.degrees(math.asin(math.sqrt(s))) for s in (s12,s13,s23)]
    J=float(np.imag(U[0,0]*U[1,1]*np.conj(U[0,1])*np.conj(U[1,0])))
    s12r,s13r,s23r=np.sqrt([s12,s13,s23])
    c12r,c13r,c23r=np.sqrt([1-s12,1-s13,1-s23])
    factor=s12r*c12r*s23r*c23r*s13r
    delta=None
    if factor>1e-10:
        sin_delta=J/(factor*c13r**2)
        cos_delta=(abs(U[1,0])**2-s12*c23r**2-c12r**2*s23*s13)/(2*factor)
        delta=float(math.degrees(math.atan2(sin_delta,cos_delta))%360)
    values.update(mixing_defined=True,mixing_blocker=None,
                  sin2_theta12=s12,sin2_theta13=s13,sin2_theta23=s23,
                  theta12_deg=angles[0],theta13_deg=angles[1],theta23_deg=angles[2],
                  J_CP=J,delta_CP_deg=delta,
                  Takagi_relative_residual=residual,
                  PMNS_ordering='columns_in_sorted_Takagi_mass_order')
    matrices={'Takagi_U':Un,'charged_lepton_left_U':Ue,'PMNS_U':U}
    return values,matrices


def run_rough_case(model: str, pair:PreparedPair, output_dir:Path, *,
                   mu_low:float=173.,vev:float=246.22,overwrite:bool=False,
                   allow_synthetic:bool=False,
                   allow_projection_residual:bool=False) -> dict:
    """Existing Workstream C common-MS calculation + reproducible observables."""
    output_dir=Path(output_dir)
    bmeta,barray,_,_=_load_pair(pair.background,B_FORMAT)
    if bmeta.get('synthetic_fixture',False) and not allow_synthetic:
        raise ValueError('Synthetic C fixture must not enter the genuine rough-results tree')
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError(f'Existing Stream3 result; use --overwrite: {output_dir}')
    stage=output_dir.with_name(output_dir.name+'.partial')
    if stage.exists():shutil.rmtree(stage)
    stage.mkdir(parents=True)
    try:
        # The only physical matching/RGE calculation is the existing C runner.
        result=run_case(wilson_manifest=pair.wilson,background_manifest=pair.background,
                        output_dir=stage,mode='rge_improved',mu_low=mu_low,vev=vev,
                        match_only=False,overwrite=False,
                        allow_projection_residual=allow_projection_residual)
        if result['SMEFT_running']['success'] is not True:
            raise RuntimeError('Existing SM+Weinberg numerical RGE did not return success')
        with np.load(stage/'C5_smeft_low.npz',allow_pickle=False) as d:
            low={name:np.asarray(d[name]).copy() for name in d.files}
        obs,matrix=neutrino_observables(low['m_nu_GeV'],barray['ye'])
        if not np.allclose(low['mass_singular_values_eV'],obs['masses_eV'],rtol=1e-8,atol=1e-13):
            raise RuntimeError('Independent sorted masses disagree with original C SVD')
        # The observable basis diagnostics are appended without altering C5.
        np.savez_compressed(stage/'C5_smeft_low.npz',**low,**matrix)
        obs.update(format='t3-rough-stream3-observables-v1',
                   status='SyntheticFixtureSmokeTest' if bmeta['synthetic_fixture'] else STATUS,
                   model=model,mu_low_gev=float(mu_low),vev_gev=float(vev),
                   mass_convention='m_nu=-(v^2/2)*C5_low_GeV^-1',
                   charged_lepton_rotation='Left eigenvectors of ye(MS) ye(MS)^dagger; invariant at one-loop in SM',
                   m_nu_GeV=_serialized_complex(low['m_nu_GeV']))
        _json(stage/'observables.json',obs)
        result['repository_main_inspected_readonly']='25eb5e13d09dd2f138176d4e36c4fc65edfa4bb6'
        result['underlying_C_status']=result['status']
        result['status']=obs['status']
        result['classification']=obs['status']
        result['input_producer_sha256']=pair.provenance
        result['benchmark_approximation']={
            'mode':'rge_improved',
            'other_independently_unavailable_finite_terms_assumed_zero':True,
            'terms_omitted':list(OMITTED),
            'one_common_MF_all_fermions':True,'one_common_MS_all_scalars':True,
            'C5_EFT1_already_evolved_and_retained':True,
            'direct_LLSS_to_LLHH_log_readded':False,
            'not_independently_certified':True,
            'allow_nonweinberg_projection_residual':bool(allow_projection_residual)}
        result['observables_file']='observables.json'
        result['scalar_matching_status']='FullTensorPairBubbleOnly'
        result['output_npz_sha256']={filename:_digest(stage/filename) for filename in
            ('C5_smeft_at_MS.npz','C5_smeft_low.npz','scalar_mass_basis.npz')}
        _json(stage/'manifest.json',result)
        if output_dir.exists():
            shutil.rmtree(output_dir)
        stage.rename(output_dir)
        return {'model':model,'status':obs['status'],'masses_eV':obs['masses_eV'],
                'dm21_sq_eV2':obs['dm21_sq_eV2'],'dm31_sq_eV2':obs['dm31_sq_eV2'],
                'dm32_sq_eV2':obs['dm32_sq_eV2'],
                'mixing_defined':obs['mixing_defined'],
                'sin2_theta12':obs['sin2_theta12'],
                'sin2_theta13':obs['sin2_theta13'],
                'sin2_theta23':obs['sin2_theta23'],
                'mu_f_gev':result['mu_f_gev'],'mu_s_gev':result['mu_s_gev'],
                'mu_low_gev':result['mu_low_gev'],'output_manifest':str(output_dir/'manifest.json'),
                'scalar_mass_kind':result['scalar_mass_kind'],
                'scalar_eigenstates_integrated_out':result['scalar_eigenstates_integrated_out'],
                'omitted_terms':list(OMITTED)}
    except Exception:
        if stage.exists():shutil.rmtree(stage)
        raise
