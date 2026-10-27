"""Independent U(1)xSU(2)xZ2 enumeration of ordinary-T3 scalar bilinears.

This is an algebraic inventory. It does not take the archived Stage20 manifest
as a definition and does not mistake a removable mass-matrix off-diagonal
entry for a redundant interaction in the full flavour/Yukawa/quartic theory.
"""
from __future__ import annotations

import json
from pathlib import Path

from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
    REP, quadratic_kind, load_manifest,
)

ROOT=Path(__file__).resolve().parents[1]
MODELS={
    'A':(-4,-2,0),'B':(-3,-1,1),'C':(-3,-1,1),
    'D':(-2,0,2),'E':(-4,-2,0,2),
}


def predict(d1:int,d2:int,alpha:int) -> list[str]:
    # Each non-SM scalar is odd under Z2; two such fields are even.
    # SU(2) singlet in R_j ⊗ R_j exists for any j; the symmetric bosonic
    # square contains it only for integer j (odd SU2 dimension).
    ops=[]
    if alpha == 0 and d1%2==1:
        ops.append('self_s1')
    if alpha == -2 and d2%2==1:
        ops.append('self_s2')
    # Two distinct irreps can form an SU(2) singlet only if d1=d2.
    # Hypercharge Y1+Y2=(alpha+1) must vanish for S1 S2.
    if alpha == -1 and d1==d2:
        ops.append('mixed')
    return ops


def audit(*,manifest_path:Path|None=None,output:Path|None=None) -> dict:
    manifest=load_manifest(manifest_path or
        ROOT/'docs/general_scalar_trial_operator_manifest.json')
    found={}
    for family,alphas in MODELS.items():
        d1,d2,_=REP[family]
        for alpha in alphas:
            name=f'{family}_{alpha:+d}'
            predicted=predict(d1,d2,alpha)
            recorded_kind=quadratic_kind(manifest,name)
            if predicted != ([] if recorded_kind=='none' else [recorded_kind]):
                raise ValueError(f'Independent quadratic count disagrees with source for {name}')
            found[name]={
                'su2_dimensions':[d1,d2],
                'alpha':alpha,
                'hypercharge_s1':alpha/2,
                'hypercharge_s2':(alpha+2)/2,
                'holomorphic_mass_operators':predicted,
                'independent_complex_mass_parameters_in_declared_basis':len(predicted),
                'physically_redundant_operator_proven':False,
                'phase_removal_or_mass_diagonalization_possible': bool(predicted),
                'simultaneous_yukawa_quartic_basis_change_required':bool(predicted),
            }
    if set(found)!=set(manifest['models']) or \
       sum(x['independent_complex_mass_parameters_in_declared_basis'] for x in found.values())!=8:
        raise ValueError('Expected precisely eight quadratic entries in sixteen models')
    report={
        'status':'IndependentGaugeInvariantQuadraticInventoryVerified',
        'models':found,
        'models_checked':16,
        'quadratic_slots':8,
        'neutral_self_bilinear_slots':6,
        'opposite_hypercharge_doublet_mixing_slots':2,
        'field_convention':'complex S1 and S2, exact Z2 odd; SU2 dimensions 1,2,3',
        'claim_of_physical_basis_redundancy':False,
        'matchete_rgbeta_coupling_map_certified':False,
        'mass_threshold_matching_completed':False,
    }
    if output is not None:
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return report


if __name__=='__main__':
    output=ROOT/'output/general_t3_quadratic_eft1/quadratic_basis_audit.json'
    result=audit(output=output)
    print(f"Independent quadratic basis: {result['quadratic_slots']} slots / "
          f"{result['models_checked']} models; physical map UNCERTIFIED")
    print(output)
