"""Formal F-hard native LLHH coefficient, evaluated over three running Majorana states.

Only original source constant and logarithmic sectors are used.  The pole
sector is independently retained as an RG diagnostic.  No physical EFT
counterterm or Wilson handoff is inferred.
"""
from __future__ import annotations
import math
import numpy as np

LOOP = 16.0 * math.pi**2
SECTORS = ('pole', 'log', 'constant')

def require(test, text):
    if not test:
        raise ValueError('F-hard numeric: '+text)

def normalise_inputs(masses, y1, y2, U, quartics):
    m=np.asarray(masses,float)
    a=np.asarray(y1,complex); b=np.asarray(y2,complex)
    u=np.asarray(U,complex)
    require(m.shape==(3,) and np.isfinite(m).all() and np.min(m)>0, 'invalid Majorana masses')
    require(a.shape==b.shape==(3,3) and np.isfinite(a).all() and np.isfinite(b).all(), 'invalid Yukawas')
    require(u.shape==(2,2) and np.isfinite(u).all() and np.linalg.norm(u.conj().T@u-np.eye(2))<1e-10, 'invalid scalar rotation')
    require(set(quartics)=={'lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1'},'wrong source quartics')
    require(all(np.isfinite(complex(z).real) and np.isfinite(complex(z).imag) for z in quartics.values()),'nonfinite quartic')
    return m,a,b,u

def evaluate_sector(ledger, *, U, quartics, masses, y1, y2):
    """Source-owned ledger, PL-only Warsaw A and conventional C5=A+A.T.

    Evaluate each Majorana eigenstate independently. No inadvertent PR addition.
    """
    from Numerical.experimental.T3BCNativeYLLHHFlavorKernel import evaluate_spurion_kernel
    require(ledger.get('model') in ('B_-1','C_-1'),'wrong ledger model')
    require(ledger.get('monomials') and len(ledger['monomials'])>0,'empty native sector')
    m,a,b,u=normalise_inputs(masses,y1,y2,U,quartics)
    total=np.zeros((3,3),complex)
    individual=[]
    hc=[]
    for j,mj in enumerate(m):
        matrices=evaluate_spurion_kernel(ledger,U=u,quartics=quartics,y1=a[:,j],
                    y2=b[:,j],MF=float(mj),hbar=1/LOOP)
        pl=np.asarray(matrices['PL'],complex)
        pr=np.asarray(matrices['PR'],complex)
        require(pl.shape==pr.shape==(3,3) and np.isfinite(pl).all() and np.isfinite(pr).all(),
                'nonfinite or malformed source matrices')
        scale=max(np.linalg.norm(pl),np.linalg.norm(pr),1e-35)
        err=float(np.linalg.norm(pr-pl.conj())/scale)
        require(err<2e-10,'native PL/PR conjugacy broken')
        c5=pl+pl.T
        require(np.linalg.norm(c5-c5.T)<1e-10*max(np.linalg.norm(c5),1e-35),
                'symmetric Warsaw coefficient failed')
        total+=c5
        individual.append(c5)
        hc.append(err)
    return total,np.asarray(individual),max(hc)

def evaluate_three(ledgers,*,masses,y1,y2,U,quartics,mu_f):
    """Formal source = constant + log(mu_f^2/M_r^2)*log coefficient.

    The pole is used for derivative diagnostics only. There is no physical
    subtraction prescription implied by the formula.
    """
    m,a,b,u=normalise_inputs(masses,y1,y2,U,quartics)
    require(np.isfinite(mu_f) and mu_f>0,'invalid matching scale')
    require(set(ledgers)==set(SECTORS),'incomplete F-hard source sectors')
    results={}
    for name in SECTORS:
        x,by,hc=evaluate_sector(ledgers[name],U=u,quartics=quartics,
                       masses=m,y1=a,y2=b)
        results[name]=(x,by,hc)
    pole,logs,constant=(results[n][0] for n in SECTORS)
    require(np.linalg.norm(pole-logs)<=2e-10*max(np.linalg.norm(pole),np.linalg.norm(logs),1e-35),
            'pole/log coefficient parity failed')
    per_log=np.log((mu_f/m)**2)
    formal=constant+np.einsum('r,rij->ij',per_log,results['log'][1])
    derivative=2*logs
    h=1e-5
    plus=constant+np.einsum('r,rij->ij',np.log((mu_f*math.exp(h)/m)**2),results['log'][1])
    minus=constant+np.einsum('r,rij->ij',np.log((mu_f*math.exp(-h)/m)**2),results['log'][1])
    err=float(np.linalg.norm((plus-minus)/(2*h)-derivative)/max(np.linalg.norm(derivative),1e-35))
    require(err<2e-9,'log derivative does not equal twice pole coefficient')
    return {
        'sector_matrix':{k:v[0] for k,v in results.items()},
        'per_heavy_sector':{k:v[1] for k,v in results.items()},
        'formal_at_mu_f':formal,
        'matching_scale_log_per_heavy':per_log,
        'derivative_dln_mu_formal':derivative,
        'finite_difference_log_derivative_relative_residual':err,
        'pole_log_matrix_relative_residual':float(np.linalg.norm(pole-logs)/max(np.linalg.norm(pole),1e-35)),
        'max_PR_PL_conjugacy_relative_residual':max(v[2] for v in results.values()),
        'renormalized_boundary':False,
        'physical_C5':None,
    }

def scale_scan(ledgers,*,masses,y1,y2,U,quartics,mu_f,ratios=(0.5,1.0,2.0)):
    raw=evaluate_three(ledgers,masses=masses,y1=y1,y2=y2,U=U,quartics=quartics,mu_f=mu_f)
    const=raw['sector_matrix']['constant'];per=raw['per_heavy_sector']['log']
    m=np.asarray(masses,float)
    out=[]
    for ratio in ratios:
        require(np.isfinite(ratio) and ratio>0,'invalid scale ratio')
        val=const+np.einsum('r,rij->ij',np.log((ratio*mu_f/m)**2),per)
        out.append({'mu_over_mu_f':float(ratio),'scale_gev':float(ratio*mu_f),
                    'matrix':val,'frobenius_norm':float(np.linalg.norm(val))})
    require(np.linalg.norm(out[1]['matrix']-raw['formal_at_mu_f'])<1e-18*max(1.,np.linalg.norm(raw['formal_at_mu_f'])),
            'central scan not equal formal evaluation')
    return raw,out
