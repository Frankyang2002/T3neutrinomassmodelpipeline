"""Provisional model-independent numerical evaluation of exported full-general T3 C5.

The symbolic expression determines Yukawa monomials and non-Yukawa couplings.
Only diagonal real fermion masses and the current Matchete/RGBeta mapping
assumptions are supported. No claim of certified full-EFT1 matching is made.
"""
from __future__ import annotations
import argparse
import json
import math
import re
import subprocess
import tempfile
from pathlib import Path
import numpy as np
from Numerical.running.SMWeinbergEvolution import SMWeinbergInitialConditions, evolve_sm_weinberg
from physics.NeutrinoMass import numerical_neutrino_mass_matrix
from physics.NeutrinoObservables import calculate_neutrino_observables

YUK = {'y1', 'y2', 'y1Conj', 'y2Conj'}
COUPLING_PATTERN = re.compile(r'Coupling\[([A-Za-z][A-Za-z0-9]*)')
MIXED = {'B_-1', 'C_-1'}

def c(v):
    if isinstance(v, dict) and set(v) == {'real','imag'}:
        return complex(float(v['real']), float(v['imag']))
    return complex(v)

def encoded(v):
    z = c(v)
    return {'real':float(z.real), 'imag':float(z.imag)}

def matrix(v):
    a = np.array([[c(z) for z in row] for row in v], dtype=complex)
    if a.shape != (3,3):
        raise ValueError(f'Expected 3x3 matrix; got {a.shape}')
    return a

def matrix_json(a):
    return [[encoded(z) for z in row] for row in a]

def mass_params(model, eft, kappa):
    m1sq, m2sq = float(eft['mS1Sq']), float(eft['mS2Sq'])
    if min(m1sq,m2sq) <= 0: raise ValueError('Nonpositive scalar mass square')
    vals = {'MS1':math.sqrt(m1sq),'MS2':math.sqrt(m2sq)}
    if model in MIXED: raise ValueError('Mixed scalar S1/S2 mass basis not implemented')
    if kappa is not None:
        if abs(kappa) == 0: raise ValueError('Zero splitting requires analytic degenerate limit')
        if model.split('_')[1] == '-2':
            sq, species = m2sq, '2'
        elif model.split('_')[1] == '+0':
            sq, species = m1sq, '1'
        else:
            raise ValueError('Self-bilinear supplied for model without neutral scalar')
        if sq - 2*abs(kappa) <= 0: raise ValueError('Tachyonic split scalar')
        vals.update({f'T3{species}PlusMass':math.sqrt(sq+2*abs(kappa)),
                     f'T3{species}MinusMass':math.sqrt(sq-2*abs(kappa)),
                     'T3NeutralPhase':encoded(np.exp(-0.5j*np.angle(kappa)))})
    return vals

def kernel_parameters(symbols, model, uv, eft, fmass, *, quartic_override=None,
                      kappa_override=None, strict_mapping=False):
    scalar = mass_params(model, eft, kappa_override)
    result = {'MF':float(fmass), **scalar}
    mapping_assumptions = []
    for name in sorted(symbols-YUK):
        if name in result: continue
        if name in ('MS1','MS2','T3NeutralPhase','T31PlusMass','T31MinusMass',
                    'T32PlusMass','T32MinusMass'):
            raise ValueError(f'Unprovided scalar mass-basis variable {name}')
        if name == 'lambdaT3':
            val = eft.get('lambdaT3',uv.get('lambdaT3'))
        elif name.startswith('lambdaT3ExtraQ'):
            candidate = name.replace('ExtraQ','GeneralQ',1)
            if quartic_override is not None:
                val = quartic_override
                mapping_assumptions.append(f'{name}: explicit benchmark override')
            elif candidate in eft:
                if strict_mapping: raise ValueError(f'Uncertified map {name} -> {candidate}')
                val = eft[candidate]
                mapping_assumptions.append(f'UNCERTIFIED {name} -> {candidate}')
            else: raise ValueError(f'Missing extra quartic {name}: tried {candidate}')
        elif name in eft and not isinstance(eft[name],list):
            val = eft[name]
        elif name in uv and not isinstance(uv[name],list):
            val = uv[name]
        else: raise ValueError(f'Unmapped Matchete parameter {name}')
        if val is None: raise ValueError(f'Missing {name}')
        result[name] = encoded(val) if c(val).imag else float(c(val).real)
    return result, mapping_assumptions

def run(args):
    data = json.loads(args.running_json.read_text(encoding='utf-8-sig'))
    model = data['model']
    if model in MIXED: raise ValueError(f'{model} requires mixed-scalar mass diagonalization')
    coefficient_text = args.coefficient.read_text(encoding='utf-8-sig')
    symbols = set(COUPLING_PATTERN.findall(coefficient_text))
    if not ({'y1','y2'} <= symbols):
        raise ValueError('Coefficient lacks required y1/y2 flavour factors; inspect expression')
    uv,eft = data['uv_at_fermion_threshold'],data['eft1_at_scalar_threshold']
    fmass_matrix = matrix(uv['MF'])
    if not np.allclose(fmass_matrix,np.diag(np.diag(fmass_matrix)),rtol=1e-10,atol=1e-7):
        raise ValueError('Only diagonal fermion masses supported')
    fmass = np.diag(fmass_matrix)
    if np.any(abs(fmass.imag)>1e-7) or np.any(fmass.real<=0):
        raise ValueError('Heavy fermion masses must be positive real')
    missing_y = (symbols & YUK)-set(uv)
    if missing_y: raise ValueError(f'Missing Yukawa arrays: {missing_y}')
    # Only models with an actual matching phase or split-mass variable need kappa.
    has_split = any(s.startswith(('T31PlusMass','T32PlusMass')) for s in symbols)
    kappa = complex(args.kappa_real,args.kappa_imag) if args.kappa_real is not None else None
    if has_split and kappa is None:
        # The proposed naming relation is not yet certified.
        if not args.use_running_bilinear:
            raise ValueError('Split-mass coefficient: provide --use-running-bilinear or --kappa-real')
        val = eft.get('muT3GeneralS2Sq' if model.endswith('-2') else 'muT3GeneralS1Sq')
        if val is None: raise ValueError('No running scalar bilinear available')
        kappa = c(val)
    quartic_override = (complex(args.quartic_real,args.quartic_imag)
                        if args.quartic_real is not None else None)
    points = []
    assumptions = []
    for i,m in enumerate(fmass.real):
        params, note = kernel_parameters(symbols,model,uv,eft,m,
                     quartic_override=quartic_override,kappa_override=kappa,
                     strict_mapping=args.strict_mapping)
        points.append({'index':i,'couplings':params})
        assumptions.extend(note)
    with tempfile.TemporaryDirectory() as td:
        inp,out = Path(td)/'points.json',Path(td)/'kernels.json'
        inp.write_text(json.dumps({'points':points,'hbar_value':args.hbar_value}),encoding='utf-8')
        command = [args.wolframscript,'-file',str(args.evaluator_wl),
                   str(args.coefficient),str(inp),str(out)]
        p = subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        if p.returncode:
            raise RuntimeError(f'Wolfram return code {p.returncode}:\n{p.stdout[-7000:]}')
        result = json.loads(out.read_text(encoding='utf-8'))
    kernels = result['results']
    ys = {y:matrix(uv[y]) for y in symbols & YUK}
    c5 = np.zeros((3,3),dtype=complex)
    channel_info = {}
    for item in kernels:
        if not item['ok']: raise RuntimeError(f'Unresolved kernel: {item}')
        factors = item['factors']
        if len(factors) != 2: raise RuntimeError(f'Unexpected Yukawa factors: {factors}')
        def yvector(factor, idx):
            match = re.fullmatch(r'T3YUK(y1Conj|y2Conj|y1|y2)(Bar|Plain)',factor)
            if not match: raise ValueError(f'Unknown Yukawa marker {factor}')
            arr = ys[match[1]][:,idx]
            return np.conjugate(arr) if match[2]=='Bar' else arr
        i = int(item['index'])
        kernel = complex(item['real'],item['imag'])
        c5 += kernel*np.outer(yvector(factors[0],i),yvector(factors[1],i))
        channel_info.setdefault(item['channel'],[]).append(encoded(kernel))
    asym = float(np.linalg.norm(c5-c5.T)/max(np.linalg.norm(c5),1e-60))
    c5 = (c5+c5.T)/2
    sm = SMWeinbergInitialConditions(gY=float(eft['gY']),g2=float(eft['g2']),
           g3=float(eft['g3']),lambdaH=float(eft['lambdaH']),
           ye=matrix(eft['ye']),yu=matrix(eft['yu']),yd=matrix(eft['yd']),K=c5)
    mu_s=float(data['scales_gev']['S'])
    status='Success'
    errors=[]
    try:
        evolved=evolve_sm_weinberg(sm,mu_s,args.mu_low).K
    except (RuntimeError,ValueError,FloatingPointError) as exc:
        errors.append(str(exc)); evolved=c5
        status='LowEnergyRGEFailed_ThresholdOnly'
    nu=numerical_neutrino_mass_matrix(evolved)
    obs=calculate_neutrino_observables(nu,ordering=args.ordering)
    outdata={'model':model,'status':'ProvisionalDirectMatching_NotFullEFT1',
             'running_input':str(args.running_json),'coefficient':str(args.coefficient),
             'low_energy_rge_status':status,'low_energy_rge_errors':errors,
             'observables_scale_gev':args.mu_low if status=='Success' else mu_s,
             'neutrino_masses_ev':obs.masses_ev.tolist(),
             'delta_m21_sq_ev2':obs.delta_m21_sq_ev2,
             'delta_m31_sq_ev2':obs.delta_m31_sq_ev2,
             'pmns_abs':obs.pmns_abs.tolist(),'takagi_residual':obs.takagi_residual,
             'flavour_unsymmetrized_relative_asymmetry':asym,
             'C5_at_scalar_threshold':matrix_json(c5),
             'C5_at_low_scale':matrix_json(evolved) if status=='Success' else None,
             'numeric_kernels':channel_info,'coupling_mapping_assumptions':sorted(set(assumptions)),
             'mass_basis_kappa':encoded(kappa) if kappa is not None else None,
             'not_certified':['Matchete-RGBeta coupling dictionary','scalar mass basis normalization',
             'EFT1 higher-dimensional Wilson running','flavour-index normalization']}
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(outdata,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(model,':',status,'masses [eV]:',obs.masses_ev.tolist(),flush=True)
    return outdata

def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--running-json',type=Path,required=True)
    p.add_argument('--coefficient',type=Path,required=True)
    p.add_argument('--output-json',type=Path,required=True)
    p.add_argument('--evaluator-wl',type=Path,default=Path('Lagrangian/EvaluateGeneralC5Channels.wl'))
    p.add_argument('--wolframscript',default='wolframscript')
    p.add_argument('--kappa-real',type=float)
    p.add_argument('--kappa-imag',type=float,default=0.)
    p.add_argument('--use-running-bilinear',action='store_true')
    p.add_argument('--quartic-real',type=float)
    p.add_argument('--quartic-imag',type=float,default=0.)
    p.add_argument('--strict-mapping',action='store_true')
    p.add_argument('--hbar-value',type=float,default=1.)
    p.add_argument('--mu-low',type=float,default=91.1876)
    p.add_argument('--ordering',default='AUTO',choices=['AUTO','NO','IO'])
    return p

def main():
    run(parser().parse_args());return 0

if __name__=='__main__':raise SystemExit(main())
