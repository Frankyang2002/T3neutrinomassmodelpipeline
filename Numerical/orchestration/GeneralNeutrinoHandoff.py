"""Continue an existing general-T3 running JSON through the existing SM+Weinberg
solver, using an explicitly supplied scalar-threshold C5.

This is a numerical handoff, NOT a general-T3 matching calculation. The C5 file
must be independently computed; never assume restricted T3 matching applies
to additional interactions.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from Numerical.running.SMWeinbergEvolution import (
    SMWeinbergInitialConditions, evolve_sm_weinberg,
)
from physics.NeutrinoTrajectory import (
    neutrino_mass_from_c5, charged_lepton_mass_basis_matrix,
)
from physics.NeutrinoObservables import calculate_neutrino_observables

def decode(x):
    if isinstance(x, dict) and set(x)=={"real","imag"}:
        return complex(x["real"],x["imag"])
    if isinstance(x,list):
        return [decode(v) for v in x]
    return x

def encode(x):
    if isinstance(x,np.ndarray):
        return encode(x.tolist())
    if isinstance(x,(list,tuple)):
        return [encode(v) for v in x]
    if isinstance(x,complex):
        return {"real":float(x.real),"imag":float(x.imag)}
    if isinstance(x,np.generic):
        return encode(x.item())
    if isinstance(x,dict):
        return {str(k):encode(v) for k,v in x.items()}
    return x

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--running-json",type=Path,required=True)
    ap.add_argument("--c5-json",type=Path,required=True,
        help="JSON object with a 'C5' symmetric 3x3 complex matrix in GeV^-1")
    ap.add_argument("--mu-low",type=float,default=91.1876)
    ap.add_argument("--vev",type=float,default=246.22)
    ap.add_argument("--ordering",choices=("AUTO","NO","IO"),default="AUTO")
    ap.add_argument("--output-json",type=Path,required=True)
    a=ap.parse_args(argv)
    running=json.loads(a.running_json.read_text(encoding="utf-8-sig"))
    if running.get("status")!="NumericalRenormalisableRunningOnly":
        raise ValueError("Expected successful general UV+EFT1 running JSON")
    c5data=json.loads(a.c5_json.read_text(encoding="utf-8-sig"))
    c5=np.asarray(decode(c5data["C5"]),dtype=complex)
    if c5.shape!=(3,3) or not np.all(np.isfinite(c5)) or not np.allclose(c5,c5.T):
        raise ValueError("C5 must be finite symmetric 3x3")
    sm=running["eft1_at_scalar_threshold"]
    mu_s=float(running["scales_gev"]["S"])
    ini=SMWeinbergInitialConditions(
        gY=float(sm["gY"]),g2=float(sm["g2"]),g3=float(sm["g3"]),
        lambdaH=float(sm["lambdaH"]),
        ye=np.asarray(decode(sm["ye"]),dtype=complex),
        yu=np.asarray(decode(sm["yu"]),dtype=complex),
        yd=np.asarray(decode(sm["yd"]),dtype=complex),K=c5,
    )
    evolved=evolve_sm_weinberg(ini,mu_s,a.mu_low)
    m=charged_lepton_mass_basis_matrix(
        neutrino_mass_from_c5(evolved.K,vev_gev=a.vev),evolved.ye)
    obs=calculate_neutrino_observables(m,ordering=a.ordering)
    result={
        "model":running["model"],
        "status":"LowEnergyObservablesFromExternallySuppliedC5",
        "general_interaction_matching_performed":False,
        "c5_source":str(a.c5_json),
        "mu_scalar_gev":mu_s,"mu_low_gev":a.mu_low,
        "C5_at_low_gev_inverse":encode(evolved.K),
        "mass_matrix_gev":encode(m),
        "observables":encode(vars(obs)),
    }
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(f"{result['model']}: SM+Weinberg and observables completed -> {a.output_json}")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
