"""Read-only, conditional flavour sensitivity of the *fermion* LLSS boundary.

Not a physical finite matching coefficient. Uses the repository's genuine joint
scalar bubble and Weinberg projection; a frozen-background coupled C_LLSS/C5
IVP is compared to the repository's fixed-order leading-log implementation.
Does not overwrite production data or infer a GammaCC->Weyl factor.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

LOOP=16*math.pi**2
REL=(.01,.05,.10)

def symmetric(C):
    return (C+C.transpose(1,0,3,2,4,5)+C.swapaxes(4,5)+C.transpose(1,0,3,2,5,4))/4

def directions(C):
    """Flavour-only transformations preserve every original gauge contraction."""
    n=float(np.linalg.norm(C))
    if not n>0: raise ValueError('Zero baseline LLSS norm')
    aligned=C.copy()
    out={'aligned':aligned}
    matrices=(np.diag([1.,-1.,0.]),np.array([[0.,1.,0.],[1.,0.,0.],[0.,0.,0.]]),
              np.diag([0.,1.,-1.]))
    for i,H in enumerate(matrices):
        D=np.einsum('pr,rqijab->pqijab',H,C,optimize=True)
        D+=np.einsum('qr,prijab->pqijab',H,C,optimize=True)
        D=symmetric(D)
        D-=np.vdot(C,D)/np.vdot(C,C)*C
        D=symmetric(D)
        if np.linalg.norm(D)>1e-10*n:
            out[f'orthogonal_flavour_{i+1}']=D
    for label,D in list(out.items()):
        out[label]=D*(n/np.linalg.norm(D))
    return out

def observables(M):
    """Majorana Takagi singular values and PMNS moduli in the charged-lepton basis.

    Angles are conditional on an assumed diagonal Ye. Mass ordering uses
    ascending singular values; degeneracies are flagged as ill-conditioned.
    """
    if M.shape!=(3,3):raise ValueError('C5 must be 3x3')
    if np.linalg.norm(M-M.T)>1e-7*max(np.linalg.norm(M),1e-40):
        raise ValueError('C5 is not complex symmetric')
    U,s,Vh=np.linalg.svd(M)
    ix=np.argsort(s); s=s[ix]; U=U[:,ix]
    # singular-vector moduli are sufficient for mixing angles
    t13=float(np.arcsin(np.clip(abs(U[0,2]),0,1)))
    t12=float(np.arctan2(abs(U[0,1]),abs(U[0,0])))
    t23=float(np.arctan2(abs(U[1,2]),abs(U[2,2])))
    return {'relative_masses':(s/s[-1]).tolist() if s[-1] else [0]*3,
            'mass_ratio_m2_m3':float(s[1]/s[2]) if s[2] else None,
            'mass_ratio_m1_m3':float(s[0]/s[2]) if s[2] else None,
            'theta12_deg':math.degrees(t12),'theta13_deg':math.degrees(t13),
            'theta23_deg':math.degrees(t23),'near_degeneracy':bool(min(np.diff(s))<1e-6*max(s[-1],1e-60))}

def evaluate(C, q, m2, mu_f, mu_s, *, rtol=1e-10):
    from Numerical.joint_scalar_matching.large_batch import _response, _beta_project
    from Numerical.joint_scalar_matching.core import scalar_bubble_matching,project_weinberg
    beta,br=_beta_project(C,q,strict=False)
    t0,t1=np.log(mu_f),np.log(mu_s)
    # Explicit coupled initial-value system for constant C_LLSS and sourced C5.
    # This is the fixed-order truncation with frozen renormalisable background,
    # NOT a claim to the full running-background seven-term coupled IVP.
    y0=np.r_[C.reshape(-1),np.zeros(9,dtype=complex)]
    def rhs(t,y):
        ll=y[:C.size].reshape(C.shape)
        b,_=_beta_project(ll,q,strict=False)
        return np.r_[np.zeros(C.size,dtype=complex),(b/LOOP).reshape(-1)]
    sol=solve_ivp(rhs,(t0,t1),y0,method='DOP853',rtol=rtol,atol=1e-25)
    if not sol.success:raise RuntimeError(sol.message)
    running=sol.y[C.size:,-1].reshape(3,3)
    fixed=_response(C,q,m2,mu_f,mu_s)
    bubble=fixed['bubble']
    total=running+bubble
    res=float(np.linalg.norm(running-fixed['running'])/max(np.linalg.norm(fixed['running']),1e-35))
    if res>2e-7:raise AssertionError('IVP disagrees with native repository LL leading log')
    if max(br,fixed['bubble_residual'])>1e-7:
        raise ValueError('Full perturbed tensor violates Weinberg Ward geometry')
    return {'matrix':total,'running':running,'bubble':bubble,'beta':beta,'ward_max':max(br,fixed['bubble_residual']),
            'native_ivp_relative_residual':res,'ivp_nfev':int(sol.nfev)}

def encode(M):return {'real':np.asarray(M).real.tolist(),'imag':np.asarray(M).imag.tolist()}

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path('.'))
    p.add_argument('--model',choices=('A_-4','D_-2'),required=True)
    p.add_argument('--boundary',type=Path,required=True,help='source-released LLSS boundary JSON')
    p.add_argument('--background',type=Path,required=True,help='joint-scalar background JSON with adjacent NPZ')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--assume-diagonal-ye',action='store_true',help='enable conditional mixing-angle extraction')
    a=p.parse_args(argv)
    root=a.repo.resolve(); out=a.output.resolve()
    allowed=root/'output'/'audits'/'worker2_fermion_matching'
    if not out.is_relative_to(allowed.resolve()):p.error('output must be under output/audits/worker2_fermion_matching')
    if out.exists():p.error('refusing to overwrite existing output directory')
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.large_batch import load_background,_rotate
    source=load_source(a.boundary,expected_model=a.model)
    bg=load_background(a.background,source)
    if not math.isclose(source.mu_f_gev,1e5) or not math.isclose(bg.mu_s,1e3):
        raise ValueError('Expected muF=1e5, muS=1e3 GeV')
    C=_rotate(source.raw,bg.U)
    if np.linalg.norm(C-symmetric(C))>1e-8*max(np.linalg.norm(C),1e-30):
        raise ValueError('Input source has unexpected exchange symmetries')
    base=evaluate(C,bg.quartic,bg.mass2,bg.mu_f,bg.mu_s)
    base_norm=float(np.linalg.norm(base['matrix']))
    rows=[]; matrices={'baseline':encode(base['matrix'])}
    for name,D in directions(C).items():
        for fraction in REL:
            for sign in (-1,1):
                key=f'{name}_{sign:+d}_{int(fraction*100)}pct'
                result=evaluate(C+sign*fraction*D,bg.quartic,bg.mass2,bg.mu_f,bg.mu_s)
                change=result['matrix']-base['matrix']
                row={'scenario':key,'category':'HYPOTHETICAL_NOT_PHYSICAL_MATCHING',
                     'relative_boundary_norm':fraction,'direction':name,'sign':sign,
                     'flavour_boundary_overlap_real':float(np.vdot(C,D).real/np.vdot(C,C).real),
                     'delta_c5_relative_frobenius':float(np.linalg.norm(change)/max(base_norm,1e-35)),
                     'c5_shape_distance':float(np.linalg.norm(result['matrix']/max(np.linalg.norm(result['matrix']),1e-35)-base['matrix']/max(base_norm,1e-35))),
                     'ward_residual':result['ward_max'], 'native_ivp_residual':result['native_ivp_relative_residual'],
                     'running_relative_norm':float(np.linalg.norm(result['running'])/max(np.linalg.norm(result['matrix']),1e-35)),
                     'scalar_bubble_relative_norm':float(np.linalg.norm(result['bubble'])/max(np.linalg.norm(result['matrix']),1e-35))}
                if a.assume_diagonal_ye:
                    ob=observables(result['matrix']); bobs=observables(base['matrix'])
                    row.update(ob)
                    for angle in ('theta12_deg','theta13_deg','theta23_deg'):
                        row['delta_'+angle]=ob[angle]-bobs[angle]
                    row['delta_mass_ratio_m2_m3']=ob['mass_ratio_m2_m3']-bobs['mass_ratio_m2_m3']
                rows.append(row); matrices[key]=encode(result['matrix'])
    out.mkdir(parents=True)
    with (out/'scenarios.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'matrices.json').write_text(json.dumps(matrices,indent=2)+'\n')
    manifest={'status':'CONDITIONAL_FIXED_BACKGROUND_COUPLED_IVP_NOT_PHYSICAL_HARD_MATCHING',
              'model':a.model,'source_sha256':hashlib.sha256(a.boundary.read_bytes()).hexdigest(),
              'background_sha256':hashlib.sha256(a.background.read_bytes()).hexdigest(),
              'muF_gev':bg.mu_f,'muS_gev':bg.mu_s,'baseline_c5':encode(base['matrix']),
              'baseline_observables_diagonal_Ye_only':observables(base['matrix']) if a.assume_diagonal_ye else None,
              'baseline_norm':base_norm,'scenario_count':len(rows),
              'max_shape_distance':max(x['c5_shape_distance'] for x in rows),
              'max_relative_c5_norm_change':max(x['delta_c5_relative_frobenius'] for x in rows),
              'rge_stage_scope':{'EFT1':'frozen background LLSS -> C5 source integrated with solve_ivp',
                                 'joint_scalar':'actual repository simultaneous scalar bubble',
                                 'SMEFT':'not propagated, missing gauge-compatible full C5 low-energy inputs',
                                 'UV':'unmodified'},
              'limitations':['No UV-hard - EFT1-soft subtraction/counterterms',
              'No physical signed GammaCC-to-Weyl normalization','Assumed fixed quartics and scalar masses',
              'Not the full running-background seven-term coupled IVP',
              'Potential flavour dependence in independent hard operator channels not bounded',
              'No physical absolute neutrino mass inference']}
    (out/'summary.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'output':str(out),'scenarios':len(rows),'max_shape_distance':manifest['max_shape_distance']},indent=2))

if __name__=='__main__':main()
