"""Batch-1 A(-4) signed-spinor obstacle audit. No unsupported certification.

Requires the previously installed a_minus4_boundary_verify.py (from Batch 1)
and its native F-only Matchete seed, complete EFT, flavour and tensor sources.
First reruns all of its actual source-derived numerical checks. Then imports a
LIVE Matchete UV-action/abstract-CC probe result and performs independent
complex 4x4 spinor and exact Grassmann contractions.

This is NOT a connected native UV scattering-amplitude measurement: native
Matchete's abstract GammaCC rules alone admit both C and -C. Neither choice
is silently selected, and no physical Wilson coefficient is exported.

Run:
  python -X utf8 -m Numerical.lo_certification.spinor.a_minus4_native_phase_audit \
    --root . --native output/general_observable/a_minus4_boundary_batch1/native_spinor_action.json
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np

from . import a_minus4_boundary_verify as boundary
from .spinors import G, wedge

FORMAT = 't3-a-minus4-signed-spinor-ambiguity-audit-v1'
NATIVE_FORMAT = 't3-a-minus4-native-spinor-action-probe-v1'
NATIVE_STATUS = 'NativeAminus4UVActionAndAbstractCCVerified_PhaseUnmeasured'
NATIVE_DEFAULT = Path('output/general_observable/a_minus4_boundary_batch1/native_spinor_action.json')
OUT_DEFAULT = Path('output/general_observable/a_minus4_boundary_batch1/native_spinor_audit.json')
ACTIONS = ('native_validated_UV_action','only_F_heavy','three_heavy_flavours',
           'both_native_topology_yukawas_present',
           'y1_native_action_contains_charge_conjugation',
           'heavy_free_action_has_MF','heavy_free_action_has_NewFermion')
SPINOR = ('native_CC_square_minus_one','native_CC_transpose_antisymmetric',
          'native_charge_conjugated_bilinear','native_chiral_projector_flip')


class NativeBlocked(ValueError):
    pass


def require(test: bool, msg: str) -> None:
    if not test:
        raise NativeBlocked(msg)


def load_native(path: Path) -> dict:
    require(path.is_file(),f'Missing live native probe result: {path}')
    obj=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj,dict) and obj.get('format')==NATIVE_FORMAT
            and obj.get('model')=='A_-4' and obj.get('status')==NATIVE_STATUS,
            'Live Matchete A(-4) action report missing, failed, or wrong schema')
    for category,keys in (('checks_action',ACTIONS),('checks_abstract_spinor',SPINOR)):
        d=obj.get(category)
        require(isinstance(d,dict) and all(d.get(k) is True for k in keys),
                'Missing/failed native checks: '+category)
    require(obj.get('native_heavy_fields')==['F']
            and obj.get('native_light_fields')==['S1','S2'],
            'Live action did not preserve F-only threshold hierarchy')
    for key,coupling in (('native_yukawa_1_terms_inputform','y1'),
                         ('native_yukawa_2_terms_inputform','y2')):
        terms=obj.get(key)
        require(isinstance(terms,list) and bool(terms)
                and all(isinstance(s,str) for s in terms)
                and any(f'Coupling[{coupling},' in s for s in terms),
                'Live native Yukawa expression absent: '+coupling)
    require('GammaCC' in ''.join(obj['native_yukawa_1_terms_inputform'])
            and all(token in obj.get('native_heavy_free_inputform','')
                    for token in ('MF','NewFermion')),
            'Live native heavy action/charge conjugation is incomplete')
    require('GammaCC' in obj.get('native_spinor_bilinear_explicit_inputform',''),
            'Live native CC bilinear not captured')
    require(obj.get('native_concrete_4x4_cc_matrix') is None
            and obj.get('native_4x4_CC_to_Weyl_epsilon_ratio') is None
            and obj.get('native_connected_UV_amplitude_calculated') is False
            and obj.get('native_source_matched_contact_amplitude_calculated') is False
            and obj.get('full_F_matching_performed_by_this_probe') is False
            and obj.get('physical_certification') is False,
            'Unexpected physical certification or numerical CC claim in abstract report')
    return obj


def numeric_spinor_convention() -> dict:
    """Two independent nonzero complex spinors: C and C^-1 give opposite PL."""
    zero=np.zeros((2,2),complex)
    unit=np.eye(2,dtype=complex)
    s2=np.array([[0,-1j],[1j,0]],complex)
    gamma0=np.block([[zero,unit],[unit,zero]])
    gamma2=np.block([[zero,s2],[-s2,zero]])
    gamma5=np.block([[-unit,zero],[zero,unit]])
    pl=(np.eye(4,dtype=complex)-gamma5)/2
    eps=1j*s2
    c=1j*gamma2@gamma0
    assert np.allclose(c.T,-c,rtol=0,atol=0)
    assert np.allclose(c@c,-np.eye(4),rtol=0,atol=0)
    p=np.array([1+2j,3-1j],complex)
    q=np.array([-2+0.4j,1+5j],complex)
    p4=np.r_[p,0j,0j]
    q4=np.r_[q,0j,0j]
    reference=q.T@eps@p
    chosen=q4.T@c@pl@p4
    opposite=q4.T@np.linalg.inv(c)@pl@p4
    require(abs(reference)>1 and abs(chosen-reference)<1e-12 and
            abs(opposite+reference)<1e-12,
            'Explicit Weyl matrix/spinor consistency failed')
    # A single signed complex source component with 2 distinct flavours;
    # no diagonal Yukawa simplifications or real-CG phase replacements.
    y1=0.8+0.3j
    y2=-0.2+0.7j
    mass=1900.0
    g=3**(-0.25)
    s1=(0.3-0.6j)/np.sqrt(2)
    s2_value=(-0.9+1.1j)/np.sqrt(2)
    source=np.conj(y1)*np.conj(y2)*g*s1*s2_value/mass
    a_plus=source*chosen
    a_minus=source*opposite
    require(abs(a_plus)>1e-6 and abs(a_plus+a_minus)<1e-14,
            'Candidate spinor phases did not produce distinct amplitudes')
    def packed(z):
        return {'real':float(z.real),'imag':float(z.imag)}
    return {
        'metric':'+---','Weyl_epsilon_12':1,
        'fixed_reference_C':'i gamma^2 gamma^0',
        'fixed_reference_C_inverse':'-i gamma^2 gamma^0',
        'reference_weyl_bilinear_q_eps_p':packed(reference),
        'C_PL_to_Weyl_ratio':packed(chosen/reference),
        'C_inverse_PL_to_Weyl_ratio':packed(opposite/reference),
        'sample_source_coefficient_gev_inv':packed(source),
        'sample_C_amplitude':packed(a_plus),
        'sample_C_inverse_amplitude':packed(a_minus),
        'native_Matchete_sign_measured':False,
        'why_not_measured':'Native GammaCC abstract algebra does not assign either numeric 4x4 matrix to its symbol.'}


def exact_dirac_saddle() -> dict:
    """Exterior-algebra verification of a DECLARED UV Dirac source action."""
    m=Fraction(13)
    a=(G.var('a0'),G.var('a1'))
    b=(G.var('b0'),G.var('b1'))
    # Declared UV L = -M (u.v) + (u.A) + (v.B), with eps_12=+1.
    # Integrating u,v by their Grassmann EOM gives + (A.B)/M.
    u=tuple(-x*(Fraction(1)/m) for x in b)
    v=tuple(x*(Fraction(1)/m) for x in a)
    uv=-m*wedge(u,v)+wedge(u,a)+wedge(v,b)
    eft=wedge(a,b)*(Fraction(1)/m)
    require(uv==eft and eft!=G.scalar(0),
            'Exact signed Dirac heavy-current elimination failed')
    return {'declared_UV_action':'-M u.eps.v + u.eps.A + v.eps.B',
            'result':'A.eps.B/M', 'mass':'13',
            'matched_Grassmann_coefficient':'1/M',
            'eom_equality':True,
            'native_Aminus4_UV_vertex_identification_established':False}


def evaluate(root:Path, native_path:Path, *,check_running:bool=True) -> dict:
    root=root.resolve()
    native_path=native_path if native_path.is_absolute() else root/native_path
    live=load_native(native_path)
    # Reconstruct from the underlying native tree, CGs, numeric flavour source,
    # and independent existing tensors -- NEVER trust a saved status field.
    conditional=boundary.verify_from_sources(root,check_running=check_running)
    require(conditional.get('physical_certification') is False
            and conditional.get('physical_C_ijab_at_MF') is None,
            'Existing source verifier has changed its physical contract')
    witness=numeric_spinor_convention()
    saddle=exact_dirac_saddle()
    return {
        'format':FORMAT,'model':'A_-4',
        'status':'NATIVE_ACTION_AND_CONDITIONAL_WILSON_VERIFIED_ABSOLUTE_PHASE_OPEN',
        'root':str(root),
        'native_probe_sha256':hashlib.sha256(native_path.read_bytes()).hexdigest(),
        'independent_source_file_sha256':conditional['sources'],
        'mu_F_gev':conditional['numeric_flavour_checks'].get('mu_F_gev',None),
        'native_action_checks':live['checks_action'],
        'native_abstract_spinor_checks':live['checks_abstract_spinor'],
        'spinor_convention_witness':witness,
        'declared_UV_dirac_Grassmann_saddle':saddle,
        'flavour_complex_3x3_boundary_recomputed':conditional['numeric_flavour_checks']['UV_evolved_3x3_flavour_recomputed'],
        'conditional_tensor_independent_nonzero':len(conditional['independent_conditional_tensor_entries']),
        'conditional_tensor_full_shape':conditional['basis_contract']['full_tensor_shape'],
        'physical_C_ijab_at_MF':None,
        'native_GammaCC_to_Weyl_factor':None,
        'certified':False,
        'blocking_probes':[
            'Native installed GammaCC needs a concrete 4x4 matrix/orientation in the same external Weyl convention.',
            'A connected UV F-exchange amplitude and the Matchete PL EFT contact amplitude must be compared with signed complex external spinors and source Yukawa vertices.',
            'Matchete triplet SU(2) component axes must be equated to the master RGE generators.',
            'The master three-family fermion-basis constructor must be checked before claiming a 45x45x12x12 physical contract.'],
        'scope':'Batch 1 source/action and phase only; no full Wilson ODE or coupling RGEs.'}


def main(argv:list[str]|None=None)->int:
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--root',type=Path,default=Path.cwd())
    cli.add_argument('--native',type=Path,default=NATIVE_DEFAULT)
    cli.add_argument('--out',type=Path,default=OUT_DEFAULT)
    cli.add_argument('--skip-running-check',action='store_true',
                     help='Fixture testing only; real physics should retain the evolved UV comparison')
    args=cli.parse_args(argv)
    try:
        obj=evaluate(args.root,args.native,check_running=not args.skip_running_check)
    except (NativeBlocked,boundary.BoundaryBlocked,ImportError,ValueError,KeyError,OSError) as exc:
        print('A(-4) native phase audit BLOCKED:',exc)
        return 2
    out=args.out if args.out.is_absolute() else args.root.resolve()/args.out
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Native A(-4) action + source audit:',obj['status'])
    print('Native absolute GammaCC phase: UNMEASURED; physical C_ijab: null')
    print('Independent conditional tensor entries:',obj['conditional_tensor_independent_nonzero'])
    print('Report:',out)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
