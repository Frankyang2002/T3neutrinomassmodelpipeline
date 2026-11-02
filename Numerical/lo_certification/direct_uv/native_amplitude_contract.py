"""Workstream 4: read-only UV source-operator / fermion-amplitude contract.

This checks the CURRENT local Matchete/RGBeta source forms, constructs exact
Dirac-doublet and Majorana-singlet Gaussian mass-insertion witnesses, and
identifies the physical factors which source text does not fix.

It NEVER certifies a numerical 4-component->Weyl map. In particular, it does
not equate RGBeta beta-function normalisations with Matchete Feynman rules.
There is no modification to the existing engine, benchmark, EFT, or matrices.
"""
from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
from typing import Any


CLASSES = {
    'A': (1, 3, 2, (-4, -2, 0)),
    'B': (2, 2, 1, (-3, -1, 1)),
    'C': (2, 2, 3, (-3, -1, 1)),
    'D': (3, 1, 2, (-2, 0, 2)),
    'E': (3, 3, 2, (-4, -2, 0, 2)),
}
SOURCES = {
    'model': 'Lagrangian/model/T3ModelCatalog.wl',
    'fields': 'Lagrangian/model/T3Fields.wl',
    'cg': 'Lagrangian/model/SU2Invariants.wl',
    'topology': 'Lagrangian/interactions/T3Topology.wl',
    'scalar_potential': 'Lagrangian/interactions/ScalarPotential.wl',
    'rgbeta': 'RGE/running/rgbeta/T3RGBetaModel.wl',
    'input_builder': 'scripts/BuildT3DirectOneLoopInputs.py',
    'core': 'Numerical/direct_one_loop/core.py',
}
# These are literal *executed-source structural* guards, not a search in comments.
# Source changes cause a BLOCKED report until the physical contract is reviewed.
GUARDS = {
    'model': (
        '"A"->{1,3,2}', '"B"->{2,2,1}', '"C"->{2,2,3}',
        '"D"->{3,1,2}', '"E"->{3,3,2}',
        '"Scalar1"->alpha/2', '"Scalar2"->(alpha+2)/2',
        '"Fermion"->(alpha+1)/2',
    ),
    'fields': (
        'DefineCoupling[y1,Indices->{Flavor,NFlavor},SelfConjugate->False]',
        'DefineCoupling[y2,Indices->{Flavor,NFlavor},SelfConjugate->False]',
        'DefineCoupling[lambdaT3,SelfConjugate->False]',
        'SelfConjugate->True', 'Mass->T3MassSpec["F",',
    ),
    'cg': (
        'T3Y1CG,{2,dF,d1},{True,True,False},{},{1}',
        'T3Y2CG,{2,dF,d2},{True,False,True},{},{1}',
        'T3MixCG,{2,2,d1,d2},{False,False,False,True},{1,2},{1,2}',
    ),
    'topology': (
        'useCConj=which===1;', 'scalarBar=which===2;',
        'cg=If[which===1,T3Y1CG,T3Y2CG]',
        'If[useCConj,CConj[NewFermion[r]],NewFermion[r]]',
        'If[useCConj,CConj[NewFermion[j,r]],NewFermion[j,r]]',
        'y[p,r]NCM[Bar[l[i,p]],f]scalarCG[cg,labels]',
        'lambdaT3[]H[i]H[j]s1s2CG[T3MixCG,labels]',
        'T3ExperimentalGeneralFermionsQ[]',
        'BuildT3ExtraYukawaCandidates',
    ),
    'scalar_potential': (
        'coupling[]Bar[H[i]]H[i]ScalarNorm[which,d,j]',
        'BuildHiggsPortalInvariantCandidates',
        'coupling[]fieldProductIf[cg===None,1,CG[cg,labels]]',
    ),
    'rgbeta': (
        'secondHeavy=If[yF===0,F,Fc]',
        'AddYukawa[y1,{Bar@S1,l,F}',
        'AddYukawa[y2,{S2,l,secondHeavy}',
        'AddFermionMass[MF,{F,F}',
        'AddFermionMass[MF,{F,Fc}',
        'dF===2,(eps[T3RGBetaRep[dF],#1,#2]&)',
        'dF===3,(del[T3RGBetaRep[dF],#1,#2]&)',
    ),
    'input_builder': (
        "'name':'y1_HC','conjugate':encoded(y1dag)",
        "'name':'y2_HC','field':encoded(y2plain)",
        'factor=cg[\'T3Y1CG\'][0,f,s].conjugate()',
        'factor=cg[\'T3Y2CG\'][0,f,s].conjugate()',
        'mm[nc+t*df+w,r*df+w]=mf[r,t]',
    ),
    'core': (
        "_complex(data['fermion_mass_gev']", "_complex(data['scalar_A_gev2']",
        "_complex(data['scalar_B_gev2']", 'b0_finite(',
    ),
}


class ContractError(ValueError):
    pass


def _require(condition: bool, msg: str) -> None:
    if not condition:
        raise ContractError(msg)


def _wl_without_comments(data: str) -> str:
    """Remove Mathematica nested comments and whitespace, preserving strings."""
    out: list[str] = []
    level, in_string, i = 0, False, 0
    while i < len(data):
        if not in_string and data.startswith('(*', i):
            level += 1; i += 2; continue
        if level and data.startswith('*)', i):
            level -= 1; i += 2; continue
        if level:
            i += 1; continue
        c = data[i]
        if c == '"' and (i == 0 or data[i-1] != '\\'):
            in_string = not in_string
        if not c.isspace():
            out.append(c)
        i += 1
    _require(level == 0 and not in_string, 'Malformed Mathematica comment/string')
    return ''.join(out)


def _python_without_comments(data: str) -> str:
    """Tokenise Python; do not allow a commented-out guard to count as code."""
    import io
    import tokenize
    tokens = tokenize.generate_tokens(io.StringIO(data).readline)
    parts = []
    for tok in tokens:
        if tok.type in (tokenize.COMMENT, tokenize.NL, tokenize.NEWLINE,
                        tokenize.INDENT, tokenize.DEDENT, tokenize.ENDMARKER):
            continue
        parts.append(tok.string)
    return ''.join(parts)


def source_audit(root: Path) -> dict[str, Any]:
    details = {}
    errors = []
    for name, rel in SOURCES.items():
        p = root / rel
        if not p.is_file():
            errors.append(f'{rel}: missing source file')
            continue
        content = p.read_text(encoding='utf-8-sig')
        compact = (_python_without_comments(content) if p.suffix == '.py'
                   else _wl_without_comments(content))
        missing = [expr for expr in GUARDS[name] if expr not in compact]
        if missing:
            errors.append(f'{rel}: required source form changed: {missing}')
        details[name] = {'path': rel, 'sha256': sha256(p.read_bytes()).hexdigest(),
                         'source_forms_satisfied': not missing,
                         'missing_forms': missing}
    return {'status': 'source_forms_verified' if not errors else 'blocked_source_forms',
            'all_required_sources_present': len(details) == len(SOURCES),
            'source_files': details, 'errors': errors}


def epsilon_doublet_mass_insertion(mass: int = 7) -> dict[str, Any]:
    """Explicit source {F,Fc} epsilon contraction in UNBARRED fund axes.

    In the ordered left-Weyl basis (F0,F1,Fc0,Fc1), the symmetric M in
    -psi.T M psi/2 contains the upper-right block mass*epsilon and its
    transpose. This is not the post-rotation identity pairing used by the
    existing numerical input builder.
    """
    _require(mass > 0, 'Mass must be positive')
    m = Fraction(mass)
    mat = [[Fraction(0) for _ in range(4)] for _ in range(4)]
    eps = [[0, 1], [-1, 0]]
    for i in range(2):
        for j in range(2):
            mat[i][2+j] = m*eps[i][j]
            mat[2+j][i] = m*eps[i][j]
    square = [[sum(mat[i][k]*mat[k][j] for k in range(4))
               for j in range(4)] for i in range(4)]
    _require(all(square[i][j] == (m*m if i == j else 0)
                 for i in range(4) for j in range(4)),
             'Epsilon Dirac mass square not scalar identity')
    inv = [[mat[i][j]/(m*m) for j in range(4)] for i in range(4)]
    # No ad-hoc x2: the source-pair Gaussian cross term gives J1.T M^-1 J2.
    _require(inv[0][3] == 1/m and inv[1][2] == -1/m,
             'Off-diagonal SU2 orientation/mass insertion wrong')
    # Negative control: delta is NOT invariant for 2 x 2 unbarred fundamentals.
    # For T3=diag(1/2,-1/2), T3^T delta + delta T3 != 0.
    # T3^T*I + I*T3 = diag(+1,-1), whereas T3^T*eps+eps*T3=0.
    T3 = [Fraction(1, 2), Fraction(-1, 2)]
    delta_residual = [[(T3[i]+T3[j]) * (1 if i == j else 0)
                       for j in range(2)] for i in range(2)]
    eps_residual = [[(T3[i]+T3[j]) * eps[i][j]
                     for j in range(2)] for i in range(2)]
    _require(any(v != 0 for row in delta_residual for v in row) and
             all(v == 0 for row in eps_residual for v in row),
             'Doublet epsilon versus delta generator Ward check failed')
    return {'model': 'A_+0', 'basis': ['F0', 'F1', 'Fc0', 'Fc1'],
            'unbarred_fundamental_epsilon': eps,
            'symmetric_mass_matrix_in_units_MF':
                [[int(x/m) for x in row] for row in mat],
            'fermion_mass_square_MF2_identity': True,
            'inverse_mixed_component_plus': str(inv[0][3]),
            'inverse_mixed_component_minus': str(inv[1][2]),
            'coefficient_of_mixed_sources_has_no_extra_factor_two': True,
            'native_GammaCC_and_RGBeta_epsilon_numeric_alignment': 'UNRESOLVED'}


def majorana_singlet_mass_insertion(mass: int = 7) -> dict[str, Any]:
    """B_-1: -MF chi chi/2, J=g1*x+g2*y; L_eff=J**2/(2MF).

    The *distinct* x*y term has coefficient g1*g2/MF, not 2*g1*g2/MF.
    Only a source-normalised MF may be substituted for the numerical mass.
    """
    m = Fraction(mass)
    _require(m > 0, 'Mass must be positive')
    lhs = Fraction(1, 2) * (2*Fraction(1, m))
    _require(lhs == Fraction(1, m), 'Majorana Gaussian cross factor failed')
    _require(lhs != 2*Fraction(1, m), 'Majorana overcount negative control')
    return {'model': 'B_-1', 'basis': ['F'],
            'quadratic_action': '-MF/2 * F F + h.c.',
            'source_action': '-F*(J1+J2) + h.c.',
            'mixed_LLSS_coefficient_in_units_y1y2': str(lhs),
            'coefficient_of_mixed_sources_has_no_extra_factor_two': True,
            'native_Matchete_MF_half_factor_and_phase': 'UNRESOLVED'}


def model_dictionary() -> list[dict[str, Any]]:
    rows = []
    for family, (ds1, ds2, df, alphas) in CLASSES.items():
        for alpha in alphas:
            majorana = df % 2 == 1 and alpha == -1
            extra = ('y1Conj' if alpha == 0 else
                     'y2Conj' if alpha == -2 else None)
            if family not in ('A', 'D', 'E'):
                extra = None
            rows.append({
                'model': f'{family}_{alpha:+d}',
                'SU2_dimensions': {'S1': ds1, 'S2': ds2, 'F': df},
                'hypercharges_Y': {'S1': alpha/2, 'S2': (alpha+2)/2,
                                   'F': (alpha+1)/2},
                'fermion_type': 'Majorana' if majorana else 'Dirac',
                'Matchete_original_y1': 'Bar[L] CConj[F] S1 * T3Y1CG + h.c.',
                'Matchete_original_y2': 'Bar[L] F Bar[S2] * T3Y2CG + h.c.',
                'RGBeta_y1': '{Bar[S1], l, F}',
                'RGBeta_y2': '{S2,l,F}' if majorana else '{S2,l,Fc}',
                'RGBeta_mass': '{F,F}' if majorana else '{F,Fc}',
                'RGBeta_SU2_mass': 'epsilon' if df == 2 else
                                   ('delta' if df == 3 else 'singlet'),
                'input_builder_vertex1': 'y1_HC in conjugate scalar channel',
                'input_builder_vertex2': 'y2_HC in field scalar channel',
                'extra_symmetry_allowed_yukawa': extra,
                'extra_yukawa_restricted_benchmark_value':
                    'assumed_zero' if extra else 'not_applicable',
                'source_operator_form_checked': True,
                'physical_component_normalisation_certified': False,
            })
    _require(len(rows) == 16, 'Incorrect number of frozen models')
    _require(sum(x['fermion_type'] == 'Majorana' for x in rows) == 2,
             'Incorrect Majorana count')
    _require(sum(x['extra_symmetry_allowed_yukawa'] is not None for x in rows) == 6,
             'Incorrect conjugate-scalar Yukawa count')
    return rows


def report(root: Path) -> dict[str, Any]:
    source = source_audit(root)
    if source['status'] != 'source_forms_verified':
        raise ContractError('; '.join(source['errors']))
    return {
        'format': 't3-ws4-native-amplitude-source-contract-v1',
        'status': 'source_operator_and_mass_insertion_verified_not_physical',
        'source_audit': source,
        'native_amplitude_witnesses': {
            'A_+0_Dirac_doublet': epsilon_doublet_mass_insertion(),
            'B_-1_Majorana_singlet': majorana_singlet_mass_insertion(),
        },
        'models': model_dictionary(),
        'model_count': 16,
        'physical_source_normalisation_certified_count': 0,
        'native_inputs_missing_for_physical_certification': [
            'Installed Matchete GammaCC numerical phase in common left-Weyl field basis',
            'RGBeta native epsilon and triplet component axes including relative phase',
            'Full y1, y2 and MF signed map under a single heavy-field rephasing',
            'Native scalar potential-to-Hessian portal and T3 mixing coefficient map',
            'Nonzero independent extra BSM operators (zero under frozen restriction)',
        ],
        'frozen_restricted_benchmark_assumptions': {
            'general_yukawas': 'y1Conj/y2Conj set to zero in six allowed models',
            'additional_quadratics': 'set to zero unless explicitly present in frozen input',
            'nonminimal_scalar_invariants': 'not included unless frozen component inputs supply them',
            'representative_nonminimal_couplings': [
                'lambdaHHdagS2S2', 'lambdaHHdagS1barS1bar',
                'lambdaHHdagS1barS2barCross', 'lambdaH1Adj',
                'lambdaH2Adj', 'lambdaS1Adj', 'lambdaS2Adj',
            ],
            'these_zeros_are_not_symmetry_predictions': True,
        },
        'formal_EFT_bridge': {
            'action': 'L = -psi.T MF psi/2 - psi.T J + h.c.',
            'tree_level_only': 'L_eff = +J.T MF^-1 J/2 + h.c.',
            'all_BSM_fermions_eliminated': True,
            'both_BSM_scalars_retained_in_EFT1': True,
            'no_MSbar_hard_one_loop_matching': True,
            'mu0_gev': 10000000,
            'muF_gev': 100000,
            'muS_gev': 1000,
            'same_scale_EFT1_matching_certified': False,
        },
        'numerical_UV_matrices_recomputed': False,
        'physical_masses_claimed': False,
    }


def audit_frozen_inputs(root: Path) -> dict[str, Any]:
    """Check all 16 *existing* frozen component tensors, without evaluating Mnu.

    The checks are necessary source-structure tests, not a proof of the four-
    component spinor phase or of scalar potential normalization. In particular,
    the Dirac input is already in the conjugate-index basis, where its mass
    off-diagonal SU2 block is identity rather than unbarred epsilon.
    """
    import numpy as np
    directory = root / 'output/general_observable/lo_direct_one_loop/provisional_inputs'

    def z(raw: Any, shape: tuple[int, ...], label: str) -> Any:
        if isinstance(raw, dict) and set(raw) == {'real', 'imag'}:
            v = np.asarray(raw['real'], float) + 1j * np.asarray(raw['imag'], float)
        else:
            v = np.asarray(raw, complex)
        _require(v.shape == shape, f'{label}: wrong shape {v.shape}, expected {shape}')
        _require(np.isfinite(v.real).all() and np.isfinite(v.imag).all(),
                 f'{label}: non-finite components')
        return v

    results = []
    for model in model_dictionary():
        label = model['model']
        path = directory / f'{label}.json'
        _require(path.is_file(), f'Missing frozen component input {path}')
        data = json.loads(path.read_text(encoding='utf-8-sig'))
        _require(data.get('model') == label, f'{label}: signed model label mismatch')
        d1, d2, df = (model['SU2_dimensions'][x] for x in ('S1', 'S2', 'F'))
        majorana = model['fermion_type'] == 'Majorana'
        n = 3*df*(1 if majorana else 2)
        ns = d1+d2
        m = z(data['fermion_mass_gev'], (n,n), f'{label} MF')
        A = z(data['scalar_A_gev2'], (ns,ns), f'{label} A')
        B = z(data['scalar_B_gev2'], (ns,ns), f'{label} B')
        scale = max(float(np.linalg.norm(m)), 1.)
        _require(np.linalg.norm(m-m.T) < 1e-10*scale,
                 f'{label}: left-Weyl mass is not symmetric')
        _require(np.linalg.norm(A-A.conj().T) < 1e-10*max(np.linalg.norm(A), 1.),
                 f'{label}: Hermitian scalar A is invalid')
        _require(np.linalg.norm(B-B.T) < 1e-10*max(np.linalg.norm(B), 1.),
                 f'{label}: scalar B is not symmetric')
        _require(np.linalg.norm(B) == 0.,
                 f'{label}: frozen restricted benchmark requires no B quadratic')
        if not majorana:
            half = n//2
            _require(np.linalg.norm(m[:half,:half]) < 1e-10*scale and
                     np.linalg.norm(m[half:,half:]) < 1e-10*scale,
                     f'{label}: incorrect Dirac left-Weyl block basis')
        vertices = data['vertices']
        _require(len(vertices) == 2 and
                 [v.get('name') for v in vertices] == ['y1_HC','y2_HC'],
                 f'{label}: unexpected restricted Yukawa vertex basis')
        _require('conjugate' in vertices[0] and 'field' not in vertices[0] and
                 'field' in vertices[1] and 'conjugate' not in vertices[1],
                 f'{label}: source scalar conjugation/orientation changed')
        v1 = z(vertices[0]['conjugate'], (3,n,ns), f'{label} y1')
        v2 = z(vertices[1]['field'], (3,n,ns), f'{label} y2')
        _require(np.linalg.norm(v1[:,:,d1:]) == 0. and
                 np.linalg.norm(v2[:,:,:d1]) == 0.,
                 f'{label}: scalar source families are swapped')
        if not majorana:
            half = n//2
            _require(np.linalg.norm(v1[:,:half,:]) == 0. and
                     np.linalg.norm(v2[:,half:,:]) == 0.,
                     f'{label}: wrong heavy-Dirac chirality vertex sector')
        results.append({'model':label,
                        'source_tensor_structure_verified':True,
                        'provisional_input_sha256':sha256(path.read_bytes()).hexdigest(),
                        'complex_fermion_states':n,
                        'complex_scalar_states':ns,
                        'native_phase_and_normalisation_verified':False})
    return {'format':'t3-ws4-existing-provisional-input-shape-audit-v1',
            'status':'restricted_component_structure_verified_not_physical',
            'models_checked':len(results), 'results':results,
            'UV_mass_matrices_recomputed':False,
            'physical_source_mapping_certified':False}
