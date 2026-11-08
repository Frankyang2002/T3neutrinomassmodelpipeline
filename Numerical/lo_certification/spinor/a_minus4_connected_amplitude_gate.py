"""A(-4) Batch-1, source-derived CG + signed Dirac-exchange gate.

This is a NEW verification tool; it does not replace the previous conditional
boundary or modify the RGE. It uses the actual A(-4) F-only Matchete source,
registered *UV* T3Y1CG/T3Y2CG tensors, native action and independent existing
conditional boundary verifier. The contact tensor is read separately from the
F-only Wilson seed, not assumed equal to a product of UV Yukawa CGs.

For a DECLARED Weyl embedding, perform a connected zero-momentum Dirac-source
saddle, with explicit complex external spinors/Yukawas and exact Grassmann
sign test. Compare its signed coefficient to the native F-only PL contact
for BOTH charge-conjugation representatives C and -C. They agree for both,
so the installed Matchete GammaCC -> master-Weyl absolute phase is *not*
measured by this test. Fail closed, never export a certified C tensor.

Run from the repository root:
  python -X utf8 -m Numerical.lo_certification.spinor.a_minus4_connected_amplitude_gate \
      --root . \
      --native output/general_observable/a_minus4_boundary_batch1/native_spinor_action.json \
      --audit output/general_observable/a_minus4_boundary_batch1/native_spinor_audit.json \
      --out output/general_observable/a_minus4_boundary_batch1/connected_amplitude_gate.json
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import subprocess
from typing import Any

import numpy as np

from . import a_minus4_boundary_verify as boundary
from . import a_minus4_native_phase_audit as phase
from .spinors import G, wedge

FORMAT = 't3-a-minus4-connected-component-two-cc-sign-gate-v1'
UV_CG_PATH = Path('output/general_observable/lo_direct_one_loop/source_cg_all16/'
                  'T3_dS1_1_dS2_3_dF_2_alpha_m4/uv_cg_registry.json')
MATCHED_SEED = boundary.RAW / 'eft1_after_F_wilson_seed.json'
DEFAULT_NATIVE = phase.NATIVE_DEFAULT
DEFAULT_AUDIT = phase.OUT_DEFAULT
DEFAULT_OUT = Path('output/general_observable/a_minus4_boundary_batch1/connected_amplitude_gate.json')


class GateBlocked(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise GateBlocked(message)


def digest(path: Path) -> str:
    require(path.is_file(), 'Missing original source: '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    require(path.is_file(), 'Missing original JSON: '+str(path))
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), 'Expected JSON object: '+str(path))
    return obj


def native_cg_sparse(item: dict, expected_name: str, shape: tuple[int, ...]) -> np.ndarray:
    """Decode the current native Matchete CSR InputForm, without interpreting reports."""
    require(item.get('Name') == expected_name, 'Missing expected native CG '+expected_name)
    source = str(item.get('TensorInputForm', ''))
    # 1-based column coordinates; pointer array partitions by FIRST index.
    regex = (r'SparseArray\[Automatic,\s*\{(?P<dims>[0-9, ]+)\},\s*0,\s*'
             r'\{1,\s*\{\{(?P<ptrs>[0-9, ]+)\},\s*\{(?P<columns>.*?)\}\},\s*'
             r'\{(?P<values>[^{}]*)\}\}\]')
    match = re.fullmatch(regex, source, re.S)
    require(match is not None, expected_name+': unknown Matchete sparse encoding')
    dims = tuple(int(z.strip()) for z in match['dims'].split(','))
    require(dims == shape, expected_name+': wrong native CG shape')
    ptrs = [int(z.strip()) for z in match['ptrs'].split(',')]
    cols = [tuple(int(z.strip()) for z in vec.split(','))
            for vec in re.findall(r'\{([^{}]+)\}', match['columns'])]
    vals = [boundary.signed_real_expr(z.strip()) for z in match['values'].split(',')]
    require(len(ptrs) == shape[0]+1 and ptrs[0] == 0 and
            ptrs[-1] == len(cols) == len(vals) and
            all(x <= y for x, y in zip(ptrs, ptrs[1:])) and
            all(len(v) == len(shape)-1 for v in cols),
            expected_name+': inconsistent sparse pointers/indices/values')
    tensor = np.zeros(shape, dtype=complex)
    occupied: set[tuple[int, ...]] = set()
    for row in range(shape[0]):
        for pos in range(ptrs[row], ptrs[row+1]):
            key = (row,)*(1) + tuple(v-1 for v in cols[pos])
            require(len(key) == len(shape) and
                    all(0 <= ix < dim for ix, dim in zip(key, shape)) and
                    key not in occupied,
                    expected_name+': invalid or duplicate sparse coordinates')
            occupied.add(key)
            tensor[key] = vals[pos]
    return tensor


def validate_actions(native: dict) -> None:
    """Check literal native terms needed for the DECLARED Weyl UV reduction."""
    y1 = native['native_yukawa_1_terms_inputform']
    y2 = native['native_yukawa_2_terms_inputform']
    require(len(y1) == len(y2) == 2, 'Unexpected native Yukawa/HC inventory')
    x1 = next((x for x in y1 if 'Bar[Coupling[y1,' in x), None)
    x2 = next((x for x in y2 if 'Bar[Coupling[y2,' in x), None)
    require(isinstance(x1, str) and isinstance(x2, str), 'Missing holomorphic HC Yukawas')
    c1 = boundary.compact(x1)
    c2 = boundary.compact(x2)
    for token in ('Bar[Field[NewScalar1,Scalar,{},{}]]',
                  'CG[Bar[T3Y1CG]', 'DiracProduct[GammaCC,Proj[-1]]',
                  'NCM[Transp[Field[l,Fermion,',
                  'Field[NewFermion,Fermion,'):
        require(token in c1, 'y1 HC action changed: '+token)
    for token in ('Field[NewScalar2,Scalar,',
                  'CG[Bar[T3Y2CG]', 'DiracProduct[Proj[-1]]',
                  'NCM[Bar[Field[NewFermion,Fermion,',
                  'Field[l,Fermion,'):
        require(token in c2, 'y2 HC action changed: '+token)
    mf = boundary.compact(native['native_heavy_free_inputform'])
    require(mf.startswith('-(Coupling[MF,') and
            'NCM[Bar[Field[NewFermion,Fermion,' in mf and
            'Field[NewFermion,Fermion,' in mf,
            'Heavy Dirac mass/source sign no longer matches probe hypothesis')
    require(native.get('native_concrete_4x4_cc_matrix') is None and
            native.get('native_4x4_CC_to_Weyl_epsilon_ratio') is None and
            native.get('native_connected_UV_amplitude_calculated') is False,
            'Native result unexpectedly claims a concrete C sign')


def native_cg_identity(uv: dict, seed: dict) -> tuple[dict[str, np.ndarray], dict]:
    """Compare UV Yukawa CG contraction with INDEPENDENTLY matched contact CG."""
    require(uv.get('status') == 'Success' and
            uv.get('metadata', {}).get('dS1') == 1 and
            uv.get('metadata', {}).get('dS2') == 3 and
            uv.get('metadata', {}).get('dF') == 2 and
            uv.get('metadata', {}).get('alpha') == -4,
            'UV native CG registry is not A(-4)')
    registry = {row['Name']: row for row in uv['CGRegistry']}
    require(set(registry) == {'T3Y1CG','T3Y2CG','T3MixCG'},
            'UV CG registry completeness changed')
    expected_reps = {
        'T3Y1CG':'{SU2L[fund],T3BSMd2}',
        'T3Y2CG':'{SU2L[fund],Bar[T3BSMd2],T3BSMd3}',
        'T3MixCG':'{Bar[SU2L[fund]],Bar[SU2L[fund]],T3BSMd3}',
    }
    for name, reps in expected_reps.items():
        require(boundary.compact(registry[name].get('RepsInputForm','')) == reps,
                'Native SU2 representation/order changed: '+name)
    g1 = native_cg_sparse(registry['T3Y1CG'], 'T3Y1CG', (2,2))
    g2 = native_cg_sparse(registry['T3Y2CG'], 'T3Y2CG', (2,2,3))
    gm_uv = native_cg_sparse(registry['T3MixCG'], 'T3MixCG', (2,2,3))
    gm_tree, proof = boundary.parse_native_cg(seed['CGRegistry'])
    require(boundary.compact(registry['T3MixCG']['TensorInputForm']) ==
            boundary.compact(seed['CGRegistry'][0]['TensorInputForm']),
            'Native UV mixing CG and matching-exported CG differ')
    require(np.allclose(g1, np.array([[0,1],[-1,0]], dtype=complex), atol=1e-14, rtol=0),
            'Native y1 singlet epsilon convention changed')
    # y1 HC contains G1*; y2 HC contains G2* and uses the Dirac mass delta.
    # This is a component source equality, NOT a guessed CG trace factor.
    from_uv = np.einsum('ah,bhc->abc', np.conjugate(g1), np.conjugate(g2))
    cgerr = float(np.max(np.abs(from_uv-gm_tree)))
    uv_err = float(np.max(np.abs(gm_uv-gm_tree)))
    require(cgerr < 1e-13 and uv_err < 1e-13,
            'Source Yukawa SU(2) contraction does not reproduce native F-only CG')
    return {'G1':g1,'G2':g2,'Gmatch':gm_tree}, {
        'y1_CG_dimensions':[2,2], 'y2_CG_dimensions':[2,2,3],
        'matched_CG_dimensions':[2,2,3],
        'UV_Y1_Y2_CG_to_Matchete_F_only_max_abs_residual':cgerr,
        'UV_registry_to_matching_registry_max_abs_residual':uv_err,
        'native_CG1_nonzero':int(np.count_nonzero(g1)),
        'native_CG2_nonzero':int(np.count_nonzero(g2)),
        'matched_nonzero':int(np.count_nonzero(gm_tree)),
        'matched_symmetric_SU2_pair':proof['cg_symmetric_pair_residual'] < 1e-13,
        'native_to_master_triplet_generator_basis_certified':False}


def exact_two_current_saddle() -> dict:
    """Zero-momentum connected Dirac exchange; exact GRASSMANN sign and 1/M."""
    mass = Fraction(17)
    q = (G.var('q0'), G.var('q1'))
    p = (G.var('p0'), G.var('p1'))
    for phase_sign in (+1,-1):
        A, B = q, tuple(phase_sign*z for z in p)
        # L = -M eta.eps.chi + eta.eps.A + chi.eps.B.
        # Eta = -B/M, chi = A/M, so L_on-shell = +A.eps.B/M.
        eta = tuple(-z*(Fraction(1)/mass) for z in B)
        chi = tuple(z*(Fraction(1)/mass) for z in A)
        uv = -mass*wedge(eta,chi) + wedge(eta,A) + wedge(chi,B)
        eff = wedge(A,B)*(Fraction(1)/mass)
        require(uv == eff and eff != G.scalar(0),
                'Connected exact Grassmann EOM fails (including mass sign)')
        require(wedge(A,B) == wedge(B,A),
                'Incorrect fermion crossing sign in Grassmann contraction')
    return {'UV_holomorphic_dirac_saddle':
            '-M eta.epsilon.chi + eta.epsilon.A + chi.epsilon.B',
            'sources':'A = y2* conj(G2) S2 ell_q; B = (CC sign) y1* conj(G1) S1* ell_p',
            'Weyl_field_convention':'epsilon_12=+1; independent anticommuting fields',
            'phase_signs_tested':[1,-1],
            'UV_EFT_saddle_equal_for_both_signs':True,
            'fixed_action_relative_dirac_tree_coefficient':'1/M',
            'mapping_native_Dirac_action_to_declared_Weyl_components_certified':False}


def sample_connected_amplitudes(cg: dict[str,np.ndarray]) -> dict:
    """Complex, offdiagonal 3-flavour, nontrivial C-phase witness per sign."""
    rng = np.random.default_rng(20261108)
    y1 = (rng.normal(size=(3,3)) + 1j*rng.normal(size=(3,3)))/5
    y2 = (rng.normal(size=(3,3)) + 1j*rng.normal(size=(3,3)))/6
    masses = np.array([113.,171.,229.])
    q = np.array([1.1+0.7j, -0.4+1.9j], complex)
    p = np.array([-0.8+1.6j, 1.3-0.5j], complex)
    s1bar = 0.4-1.5j
    s2 = -0.7+0.9j
    eps=np.array([[0,1],[-1,0]], complex)
    z=np.zeros((2,2), complex)
    unit=np.eye(2, dtype=complex)
    sigma2=np.array([[0,-1j],[1j,0]], complex)
    gamma0=np.block([[z,unit],[unit,z]])
    gamma2=np.block([[z,sigma2],[-sigma2,z]])
    gamma5=np.block([[-unit,z],[z,unit]])
    pl=(np.eye(4)-gamma5)/2
    c0=1j*gamma2@gamma0
    require(np.allclose(c0[:2,:2],eps,atol=0,rtol=0) and
            np.allclose(c0@c0,-np.eye(4),atol=0,rtol=0),
            'Declared 4-component C sign convention inconsistent')
    p4 = np.r_[p, np.zeros(2,complex)]
    q4 = np.r_[q, np.zeros(2,complex)]
    (a,b,c) = (0,0,0)  # native Gmatch[1,1,1] != 0
    def wedge_num(x,y):
        return x[0]*y[1] - x[1]*y[0]
    rows = []
    for c_sign in (+1,-1):
        C=c_sign*c0
        uv=0j
        # INDEPENDENT UV calculation: sum over three Dirac mass eigenstates
        # and both heavy SU(2) components, using the TWO source Yukawa CGs.
        for r,m in enumerate(masses):
            for h in range(2):
                A = np.conjugate(y2[1,r])*np.conjugate(cg['G2'][b,h,c])*s2*q
                B = c_sign*np.conjugate(y1[0,r])*np.conjugate(cg['G1'][a,h])*s1bar*p
                uv += wedge_num(A,B)/m
        # INDEPENDENT EFT source: contracted matching-generated CG and Wpq.
        w = boundary.dense_flavour(y1,y2,masses)[0,1]
        eft = w*cg['Gmatch'][a,b,c]*s1bar*s2*(q4.T@C@pl@p4)
        rel = float(abs(uv-eft)/max(abs(eft),1e-300))
        require(rel<5e-13 and abs(eft)>1e-7,
                'Signed zero-momentum UV/EFT component does not match')
        rows.append({'C_choice':'+' if c_sign==1 else '-',
                     'C_definition':('+' if c_sign==1 else '-')+'i gamma^2 gamma^0',
                     'UV_zero_momentum_connected':{'real':float(uv.real),'imag':float(uv.imag)},
                     'EFT_Matchete_PL_contact':{'real':float(eft.real),'imag':float(eft.imag)},
                     'relative_amplitude_residual':rel,
                     'UV_over_EFT_real':float((uv/eft).real),
                     'UV_over_EFT_imag':float((uv/eft).imag),
                     'native_GammaCC_concrete_choice_verified':False})
    require(abs(complex(**{'real':rows[0]['UV_zero_momentum_connected']['real'],
                           'imag':rows[0]['UV_zero_momentum_connected']['imag']}) +
                complex(**{'real':rows[1]['UV_zero_momentum_connected']['real'],
                           'imag':rows[1]['UV_zero_momentum_connected']['imag']}))<1e-12,
            'C sign flip did not flip full complex amplitude')
    return {'flavour_matrices':'random seeded, independent full complex 3x3 (offdiagonal)',
            'random_seed':20261108,
            'mass_eigenvalues_gev':masses.tolist(),
            'selected_external':{'lepton_p':1,'lepton_q':2,'SU2_p':1,'SU2_q':1,
                                 'triplet_component':1},
            'spinors':'two independent complex LH spinors (declared basis)',
            'candidate_amplitudes':rows,
            'both_signs_give_UV_over_EFT_equal_one':True,
            'native_phase_discriminated':False,
            'physical_master_C_tensor_certified':False}


def verify_evidence(root: Path, native_path: Path, audit_path: Path) -> dict:
    """Trust source content and computations; verify hashes, not saved PASS flags."""
    root = root.resolve()
    native_path = native_path if native_path.is_absolute() else root/native_path
    audit_path = audit_path if audit_path.is_absolute() else root/audit_path
    live = phase.load_native(native_path)
    validate_actions(live)
    audit = load(audit_path)
    require(audit.get('format') == phase.FORMAT and
            audit.get('status') == 'NATIVE_ACTION_AND_CONDITIONAL_WILSON_VERIFIED_ABSOLUTE_PHASE_OPEN' and
            audit.get('model') == 'A_-4' and
            audit.get('native_probe_sha256') == digest(native_path) and
            audit.get('certified') is False and
            audit.get('physical_C_ijab_at_MF') is None,
            'Uploaded native/audit provenance conflict or false physical certification')
    conditional = boundary.verify_from_sources(root, check_running=True)
    require(set(audit.get('independent_source_file_sha256',{})) ==
            set(conditional['sources']),
            'Missing or extra native audit source provenance records')
    require(audit.get('conditional_tensor_full_shape') == [45,45,12,12] and
            abs(float(audit.get('mu_F_gev',0)) -
                float(conditional['numeric_flavour_checks']['mu_F_gev'])) < 1e-9,
            'Native audit tensor shape or fermion matching scale changed')
    for tag, row in audit['independent_source_file_sha256'].items():
        require(row.get('sha256') == conditional['sources'][tag]['sha256'],
                'Native audit hash differs from current original source: '+tag)
    require(audit.get('conditional_tensor_independent_nonzero') ==
            len(conditional['independent_conditional_tensor_entries']) == 60,
            'Local conditional entries differ from uploaded native audit')
    uv_path = root/UV_CG_PATH
    seed_path = root/MATCHED_SEED
    uv, seed = load(uv_path), load(seed_path)
    cg, cg_report = native_cg_identity(uv,seed)
    signed = sample_connected_amplitudes(cg)
    grassmann = exact_two_current_saddle()
    return {
        'format':FORMAT,'model':'A_-4',
        'status':'NATIVE_CG_AND_DECLARED_WEYL_DIRAC_EXCHANGE_VERIFIED_BOTH_CC_SIGNS_REMAIN',
        'source_revision':conditional['git_revision'],
        'conditional_source_sha256':conditional['sources'],
        'mu_F_gev':conditional['numeric_flavour_checks']['mu_F_gev'],
        'three_heavy_fermion_mass_eigenstates_gev':
            conditional['numeric_flavour_checks']['three_mass_eigenstates_gev'],
        'source_files_sha256':{
            'native_action':digest(native_path),'native_audit':digest(audit_path),
            'UV_registered_CG':digest(uv_path),'Matchete_F_only_seed':digest(seed_path)},
        'source_term_contract':{'UV_yukawas':['y1* S1* (L^T C PL F)',
                                                 'y2* S2 (Fbar PL L)'],
                                'heavy_mass':'-MF Fbar F',
                                'holomorphic_EFT1':'y1* y2* / MF S1* S2 (Lq^T GammaCC PL Lp)',
                                'native_y1_and_y2_flavour_order':'y1[p,r] then y2[q,r]',
                                'fermion_bilinear_order':'q left / p right'},
        'source_native_SU2_CG_contraction':cg_report,
        'declared_Weyl_exact_Grassmann_saddle':grassmann,
        'declared_Weyl_complex_connected_component':signed,
        'meaning':'UV/EFT relative sign and tree 1/M agree in BOTH admissible C sign conventions.',
        'native_Cconj_to_master_Weyl_factor':None,
        'physical_C_ijab_at_MF':None,
        'physical_certification':False,
        'remaining_gates':[
            'An externally fixed native GammaCC matrix/charge-conjugation convention relative to the master Weyl fields, or an explicit, justified field rephasing dictionary.',
            'Native Matchete triplet generator axis transformation relative to master RGE scalar generators.',
            'Current true 3-generation master Weyl ordering (the 45-Weyl extension is presently a proposal).'],
        'no_RGE_or_scalar_threshold_changes':True,
    }


def main(argv: list[str]|None=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path.cwd())
    p.add_argument('--native',type=Path,default=DEFAULT_NATIVE)
    p.add_argument('--audit',type=Path,default=DEFAULT_AUDIT)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    args=p.parse_args(argv)
    try:
        result=verify_evidence(args.root,args.native,args.audit)
    except (GateBlocked,phase.NativeBlocked,boundary.BoundaryBlocked,ValueError,
            FileNotFoundError,KeyError,OSError,ImportError) as exc:
        print('A(-4) connected amplitude gate BLOCKED:',exc)
        return 2
    target=args.out if args.out.is_absolute() else args.root.resolve()/args.out
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Native UV CG1 x conj(CG2) -> matched CG: PASS')
    print('Declared Weyl Dirac UV/EFT zero-momentum amplitude: PASS for C and -C')
    print('Physical GammaCC/master Weyl mapping: BLOCKED (both signs equivalent)')
    print('Report:',target)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
