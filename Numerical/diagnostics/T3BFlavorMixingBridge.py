"""T3-B full-flavour OFF-DIAGONAL mixing study on the actual numerical background.

Calculates a matched 3x3 C12(MF) and evolves the audited flavour-blind
C12<->C5 mixing only, retaining all three light- and heavy-flavour indices.
Not the complete 3-generation anomalous dimension and NOT scalar matching.
Does not alter fixed-order C5 or generate physical neutrino predictions.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.integrate import solve_ivp

LOOP = 16 * np.pi**2
FWD = 4 * np.sqrt(3) / 3
REV = 8 * np.sqrt(3) / 3


def matrix_from_json(p):
    a = np.asarray(p['real'], float) + 1j * np.asarray(p['imag'], float)
    if a.shape != (3, 3) or not np.all(np.isfinite(a)):
        raise ValueError('Expected finite 3x3 complex matrix')
    return a


def matrix_to_json(a):
    a = np.asarray(a, complex)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist(),
            'frobenius_norm': float(np.linalg.norm(a))}


def mixed_boundary(y1, y2, masses):
    """C12=(conj(y1) diag(1/M) conj(y2)^T + 1<->2)/2."""
    y1, y2 = np.asarray(y1, complex), np.asarray(y2, complex)
    masses = np.asarray(masses, float)
    if y1.shape != (3, 3) or y2.shape != (3, 3) or masses.shape != (3,):
        raise ValueError('This benchmark requires (3,3) Yukawas and three heavy masses')
    if not np.all(np.isfinite(masses)) or np.any(masses <= 0):
        raise ValueError('Heavy masses must be finite positive values')
    first = (np.conjugate(y1) / masses) @ np.conjugate(y2).T
    c12 = (first + first.T) / 2
    if not np.allclose(c12, c12.T, atol=1e-15, rtol=1e-12):
        raise AssertionError('Boundary is not symmetric')
    return c12


def offdiagonal_evolution(samples, c12, *, rtol=1e-10, atol=1e-13):
    mus = np.asarray([float(p['mu_gev']) for p in samples], float)
    if mus.size < 2 or np.any(~np.isfinite(mus)) or np.any(mus <= 0) or not np.all(np.diff(mus)<0):
        raise ValueError('Require ordered downward MF->MS EFT1 samples')
    lam = np.asarray([complex(*p['lambdaT3']) for p in samples])
    if not np.all(np.isfinite(lam)):
        raise ValueError('Nonfinite lambdaT3 samples')
    t = np.log(mus)[::-1]
    l = lam[::-1]
    def coupling(tnow):
        return np.interp(tnow,t,l.real)+1j*np.interp(tnow,t,l.imag)
    c12 = np.asarray(c12,complex)
    if c12.shape != (3,3):raise ValueError('C12 must be 3x3')
    initial = np.concatenate((c12.ravel(),np.zeros(9,complex)))
    def rhs(tnow, v, *, reverse=True):
        s,k = v[:9],v[9:]
        lamnow = coupling(tnow)
        return np.concatenate(((REV*np.conjugate(lamnow)*k/LOOP if reverse else np.zeros(9,complex)), FWD*lamnow*s/LOOP))
    t_eval = np.log(mus)
    coupled = solve_ivp(rhs,(t_eval[0],t_eval[-1]),initial,rtol=rtol,atol=atol,method='DOP853',t_eval=t_eval)
    frozen = solve_ivp(lambda x,y:rhs(x,y,reverse=False),(t_eval[0],t_eval[-1]),initial,rtol=rtol,atol=atol,method='DOP853',t_eval=t_eval)
    if not (coupled.success and frozen.success):raise RuntimeError('Wilson integration failed')
    s_final=coupled.y[:9,-1].reshape(3,3)
    k_final=coupled.y[9:,-1].reshape(3,3)
    k_forward=frozen.y[9:,-1].reshape(3,3)
    # Validate the leading-log (frozen MF coupling) coefficient separately.
    leading_log=(FWD*lam[0]/LOOP)*np.log(mus[-1]/mus[0])*c12
    norm=float(np.linalg.norm(k_forward))
    return {'c12_at_mf':matrix_to_json(c12),
            'c12_at_ms_offdiagonal_only':matrix_to_json(s_final),
            'c5_at_ms_coupled_offdiagonal_only':matrix_to_json(k_final),
            'c5_at_ms_forward_only_running_lambda':matrix_to_json(k_forward),
            'c5_at_ms_frozen_mf_leading_log':matrix_to_json(leading_log),
            'relative_feedback_on_offdiagonal_c5':float(np.linalg.norm(k_final-k_forward)/norm) if norm else None,
            'coupled_minus_forward_only':matrix_to_json(k_final-k_forward),
            'c12_transpose_residual':float(np.linalg.norm(s_final-s_final.T)),
            'c5_transpose_residual':float(np.linalg.norm(k_final-k_final.T)),
            'nfev_coupled':int(coupled.nfev)}


def run(config, uv_rgbeta, control):
    from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
    from Numerical.fitting.ScanCLI import build_uv_state_from_config
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.FinalC5TrajectoryAdapter import _majorana_takagi_mass_basis
    cfg = json.loads(Path(config).read_text(encoding='utf-8'))
    source = json.loads(Path(control).read_text(encoding='utf-8'))
    if source.get('status')!='FixedOrderControlComplete' or source.get('model')!='T3-B alpha=-1':
        raise ValueError('Require verified T3-B alpha=-1 control file')
    raw=dict(cfg);raw.setdefault('scan',{'bindings':{}})
    state=build_uv_state_from_config(SimpleNamespace(raw=raw),{}).validated()
    r=state.representation
    if (r.d_s1,r.d_s2,r.d_f,r.alpha,r.shared_scalar)!=(2,2,1,-1,False):
        raise ValueError('Only T3-B alpha=-1 is supported')
    mf=float(source['mu_f_gev'])
    uv=run_uv_segment(state,load_rgbeta_payload(uv_rgbeta),mf,
                      rtol=float(cfg.get('numerical',{}).get('rtol',1e-8)),
                      atol=float(cfg.get('numerical',{}).get('atol',1e-11)))
    basis=_majorana_takagi_mass_basis(uv.final_state,
                                     symmetry_atol=1e-10,residual_rtol=1e-9,
                                     residual_atol=1e-10)
    c12=mixed_boundary(basis.y1,basis.y2,basis.masses_gev)
    samples=source['intermediate_running']['samples']
    if not np.isclose(samples[0]['mu_gev'],mf,rtol=1e-10):
        raise ValueError('UV boundary scale and control background do not match')
    if not np.isclose(samples[-1]['mu_gev'],float(source['mu_s_gev']),rtol=1e-10):
        raise ValueError('Scalar threshold and background do not match')
    results=offdiagonal_evolution(samples,c12)
    return {'status':'OffDiagonalFullFlavorDiagnosticComplete','model':'T3-B alpha=-1',
            'scope':'3x3 full-flavour matched boundary with verified off-diagonal mixing only',
            'heavy_masses_gev':basis.masses_gev.tolist(),
            'heavy_takagi_residual':float(basis.takagi_residual),
            'mu_f_gev':mf,'mu_s_gev':float(source['mu_s_gev']),
            'fixed_order_c5_at_ms':source['fixed_order_c5_at_ms'],
            **results,
            'physical_c5_prediction_available':False,
            'neutrino_mass_prediction_available':False,
            'limitations':['No full-flavour diagonal/Wilson anomalous dimensions',
                           'No additional induced flavour operator closure checked',
                           'Zero C5(MF) boundary used for this isolated mixing response',
                           'No scalar threshold hard contribution or subtraction',
                           'Cannot add this result to fixed-order C5 without double counting']}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--uv-rgbeta',type=Path,required=True)
    p.add_argument('--control',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=run(a.config,a.uv_rgbeta,a.control)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['status','heavy_masses_gev',
              'relative_feedback_on_offdiagonal_c5','c12_transpose_residual',
              'c5_transpose_residual','physical_c5_prediction_available']},indent=2))

if __name__=='__main__':main()
