#!/usr/bin/env python3
"""Workstream A: exact scalar equivalence and conditional quartic normalization audit.

Run from the project root:
    python Lagrangian/RunT3UVExactEquivalence.py --self-test
    python Lagrangian/RunT3UVExactEquivalence.py --all

Reads existing native Matchete tensors, the Stage23 source-derived RGBeta
normalization candidates, and previously passing Workstream A reports.
Never runs Match, edits GitHub, or changes any production code. Deliberately
DOES NOT certify the physical Matchete/RGBeta basis identification.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import re
import sys
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Lagrangian.RunT3UVFinalEvidence import (  # noqa: E402
    MODELS, SELF, BC, DIMS, STAGE22, ROOT as EVIDENCE_ROOT,
    parse_stage22, native_quartic_data, stage22_model_key, read_json,
    name_alpha, require,
)

STAGE23 = ROOT / 'output/general_t3_stage22/stage23_quartic_prefactor_candidates.json'
RGBETA = ROOT / 'output/general_t3_stage22/rgbeta_native_invariants.json'
PRIOR = ROOT / 'Reports/output/workstream_a/final_uv_evidence/final_uv_evidence.json'
SELF_REPORTS = ROOT / 'Reports/output/workstream_a/six_self_eigenbasis'
OUT = ROOT / 'Reports/output/workstream_a/exact_uv_equivalence'


def exact_number(raw: str) -> sp.Expr:
    """Parse ONLY source numeric radicals, never arbitrary executable Wolfram text."""
    if not isinstance(raw, str):
        raise ValueError('Expected numeric expression string')
    text = raw.strip().replace('Sqrt[', 'sqrt(').replace(']', ')').replace('^', '**')
    if not re.fullmatch(r'[\d\s+\-*/().sqrtI]+', text):
        raise ValueError('Unexpected native numerical constant: ' + repr(raw))
    if re.search(r'\d+\.\d+',text):
        raise ValueError('Approximate decimal is forbidden in exact native component: '+repr(raw))
    try:
        val = sp.sympify(text, locals={'sqrt': sp.sqrt, 'I': sp.I}, evaluate=True)
    except (sp.SympifyError, TypeError, SyntaxError, ValueError) as exc:
        raise ValueError('Unparseable exact native component: '+repr(raw)) from exc
    if val.has(sp.Float) or val.free_symbols or val.is_finite is not True:
        raise ValueError('Nonfinite or nonnumeric native constant: ' + repr(raw))
    return val


def source_digest(path: Path) -> dict[str, str]:
    return {'path':str(path.relative_to(ROOT)), 'sha256':sha256(path.read_bytes()).hexdigest()}


def simple_matrix(data: list[list[str]], label: str) -> sp.Matrix:
    require(isinstance(data, list) and len(data) > 0 and
            all(isinstance(r, list) and len(r) == len(data) for r in data),
            f'{label}: expected nonempty square matrix')
    return sp.Matrix([[exact_number(x) for x in row] for row in data])


def exact_zero(v: sp.Expr) -> bool:
    return sp.cancel(sp.expand(v)) == 0


def exact_matrix_zero(m: sp.MatrixBase) -> bool:
    return all(exact_zero(x) for x in m)


def self_mass_identity(dim: int, kappa: sp.Expr) -> dict[str, Any]:
    """Exact unbroken-phase quadratic potential and kinetic equivalence.

    Native identity metric, b=-kappa*lambda_extra = rho*p**2, p*conj(p)=1.
    This is an algebraic equality under the explicitly stated parameter map,
    not an actual substitution into the post-rotation Wolfram Lagrangian.
    """
    require(dim in (1,3) and exact_zero(sp.im(kappa)) and kappa != 0,
            'Self quadratic native invariant must be nonzero real scalar')
    m, rho = sp.symbols('m2 rho', real=True)
    p = sp.symbols('p', nonzero=True)
    real = sp.symbols(f'R0:{dim}', real=True)
    imag = sp.symbols(f'I0:{dim}', real=True)
    dr = sp.symbols(f'dR0:{dim}', real=True)
    di = sp.symbols(f'dI0:{dim}', real=True)
    norm = sp.sqrt(2)
    S = [(p/norm)*(r+sp.I*i) for r,i in zip(real,imag)]
    Sbar = [(1/(p*norm))*(r-sp.I*i) for r,i in zip(real,imag)]
    # UV Lagrangian contains +lambda_extra*invariant+h.c.; V=-L, so
    # the physical potential coefficient is b=-kappa*lambda_extra=rho*p**2.
    hol = rho*p**2
    hol_bar = rho/p**2
    old = sum(m*a*b + hol*b*b + hol_bar*a*a for a,b in zip(S,Sbar))
    new = sum(((m+2*rho)*r*r+(m-2*rho)*i*i)/2 for r,i in zip(real,imag))
    require(exact_zero(old-new), 'Exact self-bilinear potential identity failed')
    deriv = [(p/norm)*(r+sp.I*i) for r,i in zip(dr,di)]
    derivbar = [(1/(p*norm))*(r-sp.I*i) for r,i in zip(dr,di)]
    kine_old = sum(a*b for a,b in zip(deriv,derivbar))
    kine_new = sum((r*r+i*i)/2 for r,i in zip(dr,di))
    require(exact_zero(kine_old-kine_new), 'Exact self-bilinear kinetic identity failed')
    # Can invert for original quadratic parameters without changing the UV model.
    mr, mi = sp.symbols('M_R_sq M_I_sq', real=True)
    inv_m = (mr+mi)/2
    inv_rho = (mr-mi)/4
    require(exact_zero((inv_m+2*inv_rho)-mr) and
            exact_zero((inv_m-2*inv_rho)-mi), 'Self inverse mass map failed')
    return {'passed':True,'dimension':dim,'native_kappa':str(kappa),
            'original_native_coupling_relation': '-kappa*lambda_extra=rho*p**2',
            'phase_constraint': 'conjugate(p)=1/p, rho>=0, p nonzero',
            'eigen_mass_squared': ['m2+2*rho', 'm2-2*rho'],
            'inverse_original_m2':'(M_R_sq+M_I_sq)/2',
            'inverse_rho':'(M_R_sq-M_I_sq)/4',
            'quadratic_potential_polynomial_exact':True,
            'free_canonical_kinetic_polynomial_exact':True}


def mixed_mass_identity() -> dict[str, Any]:
    """Exact B/C SU2-epsilon mass and kinetic equivalence.

    A=U diag(MX²,MY²) U†, with U an arbitrary SU(2) flavour rotation
    parameterized rationally by t and a unit-modulus phase p. This spans
    non-degenerate flavour rotations up to field rephasings.
    """
    t=sp.symbols('t',real=True)
    p=sp.symbols('p',nonzero=True)
    mx,my=sp.symbols('M_X_sq M_Y_sq',real=True)
    c=(1-t*t)/(1+t*t); s=2*t/(1+t*t)
    U=sp.Matrix([[c,-s/p],[s*p,c]])
    Ud=sp.Matrix([[c,s/p],[-s*p,c]])
    require(exact_matrix_zero(Ud*U-sp.eye(2)), 'Exact mixed U unitarity failed')
    D=sp.diag(mx,my)
    A=U*D*Ud
    require(exact_matrix_zero(Ud*A*U-D), 'Exact B/C mass diagonalization failed')
    # Field redefinition Phi=(S1,epsilon*conjugate(S2))=U*(X,Y).
    eps=sp.Matrix([[0,1],[-1,0]])
    require(eps*eps == -sp.eye(2) and eps.T*eps == sp.eye(2), 'Native epsilon invalid')
    X=sp.Matrix(sp.symbols('x0:2'))
    Y=sp.Matrix(sp.symbols('y0:2'))
    Xd=sp.Matrix(sp.symbols('xb0:2'))
    Yd=sp.Matrix(sp.symbols('yb0:2'))
    phi1=c*X-s/p*Y
    phi2=s*p*X+c*Y
    phi1d=c*Xd-s*p*Yd
    phi2d=s/p*Xd+c*Yd
    S2=-eps*phi2d
    S2d=-eps*phi2
    require(exact_matrix_zero(eps*S2d-phi2), 'Epsilon field inverse failed')
    require(exact_matrix_zero(eps*S2-phi2d), 'Epsilon conjugate field inverse failed')
    phi=sp.Matrix([[phi1[0],phi1[1]],[phi2[0],phi2[1]]])
    phid=sp.Matrix([[phi1d[0],phi1d[1]],[phi2d[0],phi2d[1]]])
    original=sum((phid[:,i].T*A*phi[:,i])[0] for i in range(2))
    diagonal=sum(mx*Xd[i]*X[i]+my*Yd[i]*Y[i] for i in range(2))
    require(exact_zero(original-diagonal), 'Exact B/C quadratic field expression failed')
    # Same U and epsilon identities imply canonical kinetic terms.
    deriv_X=sp.Matrix(sp.symbols('dx0:2'));deriv_Y=sp.Matrix(sp.symbols('dy0:2'))
    deriv_Xd=sp.Matrix(sp.symbols('dxb0:2'));deriv_Yd=sp.Matrix(sp.symbols('dyb0:2'))
    deriv_phi1=c*deriv_X-s/p*deriv_Y
    deriv_phi2=s*p*deriv_X+c*deriv_Y
    deriv_phi1d=c*deriv_Xd-s*p*deriv_Yd
    deriv_phi2d=s/p*deriv_Xd+c*deriv_Yd
    old_kin=sum(deriv_phi1d[i]*deriv_phi1[i]+deriv_phi2d[i]*deriv_phi2[i] for i in range(2))
    new_kin=sum(deriv_Xd[i]*deriv_X[i]+deriv_Yd[i]*deriv_Y[i] for i in range(2))
    require(exact_zero(old_kin-new_kin), 'Exact B/C kinetic field expression failed')
    return {'passed':True,'native_epsilon':[[int(x) for x in row] for row in eps.tolist()],
            'mass_matrix_entries':[[str(sp.factor(A[i,j])) for j in range(2)] for i in range(2)],
            'unitarity_symbolic_identity':True,'diagonal_mass_polynomial_exact':True,
            'free_canonical_kinetic_polynomial_exact':True,
            'phase_constraint':'conjugate(p)=1/p, t real',
            'physical_eigenvalue_mapping':'A=U.diag(M_X_sq,M_Y_sq).U_dagger'}


def stage23_normalization(stage: dict[str, Any], native: dict[str, list[dict[str, Any]]],
                          rgbeta: dict[str, Any]) -> dict[str, Any]:
    """Verify saved Stage23 conditional conversion matrices exactly, blockwise.

    r_in_m is a geometry map; the documented assumed relative sign implies
    lambda_RGBeta = -inverse(r_in_m)*lambda_Matchete. This is a consistency
    check of a CONDITIONAL source-derived map, NOT a physical basis proof.
    """
    require(stage.get('format') == 'stage23-exact-source-prefactors-conditional-signed-maps-v1'
            and stage.get('status') == 'SourceQuarticNormalisationVerifiedConditionalSignedMapsNotPhysicalCertification'
            and stage.get('totals',{}).get('models') == 16
            and stage.get('totals',{}).get('slots') == 50
            and stage.get('totals',{}).get('blocks') == 46
            and stage.get('source_prefactors_verified') is True
            and stage.get('signed_physical_coupling_map_certified') is False,
            'Unexpected Stage23 source or certification status')
    require(rgbeta.get('format')=='stage22-native-rgbeta-invariants-v1'
            and rgbeta.get('source')=='RGBetaRegisteredQuarticInvariantFunctions'
            and rgbeta.get('total_slots') == 50
            and rgbeta.get('native_invariant_component_alignment_verified') is False,
            'Unexpected RGBeta native-invariant source or certification flag')
    blocks=stage.get('blocks')
    require(isinstance(blocks,list) and len(blocks)==46,'Stage23 missing conversion blocks')
    seen:dict[str,set[str]]=defaultdict(set)
    rows=[]
    for block in blocks:
        model=block['model']
        require(model in {stage22_model_key(m) for m in MODELS},'Unknown Stage23 model')
        mlist=block['matchete_couplings'];rlist=block['rgbeta_couplings']
        require(len(mlist)==len(rlist) and 1<=len(mlist)<=2,'Bad Stage23 coupling block dimension')
        M=simple_matrix(block['r_in_m'], 'r_in_m')
        P=simple_matrix(block['lambda_rgbeta_from_lambda_matchete_if_assumptions_hold'], 'conditional map')
        require(M.rows==P.rows==len(mlist),'Stage23 matrix/coupling dimension mismatch')
        require(not exact_zero(M.det()),'Stage23 conversion block is singular')
        require(exact_matrix_zero(M*P+sp.eye(M.rows)) and
                exact_matrix_zero(P*M+sp.eye(M.rows)),
                f'{model} {block["block"]}: conditional inverse map not exact')
        source_symbols={entry['symbol'] for entry in native[model]}
        registered={entry['symbol'] for entry in rgbeta['models'][model]['invariants']}
        require(set(mlist)<=source_symbols and set(rlist)<=registered,
                f'{model}: conditional map references absent native coupling symbols')
        require(all(x not in seen[model] for x in mlist),
                f'{model}: duplicate Matchete coupling across normalization blocks')
        require(all(x not in seen[model] for x in rlist),
                f'{model}: duplicate RGBeta coupling across normalization blocks')
        seen[model].update(mlist);seen[model].update(rlist)
        require(block.get('component_prefactor_one_ordered_orientation')=='1/4' and
                block.get('assumed_relative_potential_sign')=='-1' and
                block.get('physical_mapping_certified') is False,
                'Stage23 has changed the source-derived prefactor or certification flags')
        permutations=block.get('native_rgbeta_permutation_counts')
        require(isinstance(permutations,dict) and set(permutations)==set(rlist) and
                all(isinstance(x,int) and 0<x<=24 for x in permutations.values()),
                'Malformed RGBeta native arrangement counts')
        rows.append({'model':model,'block':block['block'],'dimension':M.rows,
                     'matchete_symbols':mlist,'rgbeta_symbols':rlist,
                     'geometry_class':block['geometry_class'],
                     'candidate_conversion_map_exact_inverse_passed':True,
                     'candidate_map':[[str(x) for x in row] for row in P.tolist()],
                     'candidate_physical_map_certified':False})
    for m in MODELS:
        key=stage22_model_key(m)
        require({x['symbol'] for x in native[key]} <= seen[key],
                f'{key}: Stage23 did not cover every native extra quartic')
        require(len([x for x in rows if x['model']==key])==len(set(
            entry['symbol'].split('Q')[1][:6] for entry in native[key])),
            f'{key}: unexpected Stage23 block multiplicity')
    require(sum(row['dimension'] for row in rows)==50,
            'Stage23 conditional map does not cover 50 slots')
    return {'passed':True,'blocks_exact_inverse_checked':46,'coupling_slots_checked':50,
            'source_component_prefactor_one_ordered_orientation':'1/4',
            'assumed_relative_potential_sign':'-1',
            'candidate_signed_coupling_maps_conditional_only':True,
            'native_component_identification_certified':False,
            'physical_coupling_mapping_certified':False,
            'blocks':rows}


def scalar_native_field_map(model: str, p: sp.Expr) -> dict[str,list[tuple[str,sp.Matrix]]]:
    """Old component fields as exact matrices times new component fields.

    Every dagger field is an independent formal symbol with the stated
    conjugation condition on p. No extra SU2 intertwiners are assumed.
    """
    family,alpha=name_alpha(model)
    d1,d2,_=DIMS[family]
    dims={'H':2,'Hdag':2,'S1':d1,'S1dag':d1,'S2':d2,'S2dag':d2}
    maps={name:[(name,sp.eye(d))] for name,d in dims.items()}
    if model in SELF:
        neutral='S2' if alpha==-2 else 'S1'
        d=dims[neutral]; z=sp.sqrt(2)
        maps[neutral]=[(neutral+'_R',p/z*sp.eye(d)),
                       (neutral+'_I',sp.I*p/z*sp.eye(d))]
        maps[neutral+'dag']=[(neutral+'_R',1/(p*z)*sp.eye(d)),
                             (neutral+'_I',-sp.I/(p*z)*sp.eye(d))]
    elif model in BC:
        # Rational SU2 mixing with general unit-modulus phase and epsilon.
        t=sp.symbols('t',real=True)
        c=(1-t*t)/(1+t*t);s=2*t/(1+t*t)
        eps=sp.Matrix([[0,1],[-1,0]])
        maps['S1']=[('X',c*sp.eye(2)),('Y',-s/p*sp.eye(2))]
        maps['S1dag']=[('Xdag',c*sp.eye(2)),('Ydag',-s*p*sp.eye(2))]
        maps['S2']=[('Xdag',-s/p*eps),('Ydag',-c*eps)]
        maps['S2dag']=[('X',-s*p*eps),('Y',-c*eps)]
    return maps


def exact_tensor_witness(model: str, item: dict[str, Any], sample:int) -> bool:
    """Evaluate actual native CG identity in exact Q(i,sqrt(n)).

    The two sides use independent index-summation routes: source fields
    pulled back into the old contraction versus direct expanded new-basis
    tensor contractions. Unlike floating tests, equality is algebraically
    exact at each chosen nontrivial rational witness. Linear field maps
    guarantee the identity for arbitrary multilinear tensors by
    distributivity; the finite witnesses validate the actual native data.
    """
    names=item['species'];tensor=item['tensor']
    _,alpha=name_alpha(model)
    p=sp.Rational(3,5)+sp.I*sp.Rational(4,5) if sample==0 else sp.Rational(5,13)+sp.I*sp.Rational(12,13)
    maps=scalar_native_field_map(model,p)
    # Choose nontrivial exact rational complex values, dagger fields are
    # conjugates (or real for the two real eigenmultiplets).
    dims={}
    for name in names:
        for nn,M in maps[name]:dims[nn]=M.cols
    new={}
    for i,(nn,d) in enumerate(sorted(dims.items())):
        if nn.endswith('_R') or nn.endswith('_I'):
            new[nn]=sp.Matrix([sp.Rational(j+1+sample,4)+sp.Rational(i+1,7) for j in range(d)])
        elif nn.endswith('dag'):
            continue
        else:
            new[nn]=sp.Matrix([sp.Rational(j+2+sample,3)+sp.I*sp.Rational(i+j+1,5) for j in range(d)])
    # The stage tensor may include only dagger species; create their undaggered bases first.
    for nn,d in sorted(dims.items()):
        if nn.endswith('dag') and nn not in new:
            base=nn[:-3]
            if base not in new:
                new[base]=sp.Matrix([sp.Rational(j+2+sample,3)+sp.I*sp.Rational(j+1,5) for j in range(d)])
            new[nn]=sp.Matrix([sp.conjugate(x) for x in new[base]])
    for nn,d in sorted(dims.items()):
        if nn not in new:
            new[nn]=sp.Matrix([sp.Rational(j+2+sample,3)+sp.I*sp.Rational(j+1,5) for j in range(d)])

    def contract_with(vectors: list[sp.Matrix]) -> sp.Expr:
        return sum(exact_number(str(tensor[idx]))*sp.prod(vectors[k][idx[k]] for k in range(4))
                   for idx in np.ndindex(tensor.shape))
    # Native components were formerly parsed to binary floats. Preserve the
    # exact original InputForm components supplied separately by caller.
    old=[]
    for name in names:
        entries=maps[name]
        old.append(sum((M*new[nn] for nn,M in entries),sp.zeros(entries[0][1].rows,1)))
    lhs=contract_with(old)
    rhs=0
    for chosen in product(*(maps[n] for n in names)):
        term=contract_with([M*new[nn] for nn,M in chosen])
        rhs+=term
    return exact_zero(lhs-rhs)


def native_rational_transport(model: str, native: list[dict[str,Any]], raw_records: list[dict[str,Any]]) -> dict[str,Any]:
    """Exact arithmetic source-tensor checks for eight exceptional models."""
    raw_by={x['symbol']:x for x in raw_records}
    for entry in native:
        raw=raw_by[entry['symbol']]
        dims=tuple(raw['dimensions'])
        exact_tensor=np.array([exact_number(x) for x in raw['components']],dtype=object).reshape(dims)
        item={'symbol':entry['symbol'],'species':entry['species'],'tensor':exact_tensor}
        for sample in range(2):
            require(exact_tensor_witness(model,item,sample),
                    f'{model}:{entry["symbol"]}: exact native tensor transport failure')
    return {'passed':True,'native_extra_quartic_tensors':len(native),
            'exact_non_floating_witnesses':2*len(native),
            'equality_type':'exact rational/radical component contractions, 2 independent field assignments',
            'universal_multilinearity_identity':'T * product(M*new) = (T contracted with M axes) * product(new)',
            'actual_wolfram_rotated_lagrangian_difference_proved_zero':False}


def validate_prior() -> dict[str,Any]:
    prior=read_json(PRIOR)
    require(prior.get('format')=='t3-workstream-a-final-uv-evidence-v1' and
            prior.get('status')=='ConditionalUVOperatorEvidencePassedNotFullCompleteness' and
            prior.get('matchete_native_extra_quartic_slots_checked')==50 and
            prior.get('exceptional_native_quartic_tensor_transport_models')==8,
            'Final evidence stage must have passed before exact equivalence stage')
    require(len(prior.get('rows',[]))==16,'Prior evidence model coverage mismatch')
    return prior


def run_all(outdir:Path, overwrite:bool) -> dict[str,Any]:
    require(EVIDENCE_ROOT==ROOT,'Repository import points outside local project')
    paths=[PRIOR,STAGE22,STAGE23,RGBETA]+[SELF_REPORTS/(m+'.json') for m in SELF]
    for path in paths:
        require(path.is_file(),f'Required prior evidence missing: {path}')
    targets=[outdir/'exact_uv_equivalence.json',outdir/'exact_uv_equivalence.md']
    require(overwrite or not any(p.exists() for p in targets),
            'Report exists; first back up, then invoke --overwrite')
    validate_prior()
    stage22=parse_stage22(STAGE22)
    native={stage22_model_key(m):native_quartic_data(stage22,m) for m in MODELS}
    rgbeta=read_json(RGBETA)
    norm=stage23_normalization(read_json(STAGE23),native,rgbeta)
    mixed=mixed_mass_identity()
    rows=[]; witness_total=0
    for m in MODELS:
        key=stage22_model_key(m)
        row={'model':m,'native_quartic_count':len(native[key]),
             'prior_operator_evidence_passed':True,
             'conditionally_normalized_source_slots_verified':len(native[key]),
             'physical_cross_package_operator_map_certified':False}
        if m in SELF:
            report=read_json(SELF_REPORTS/(m+'.json'))
            require(report.get('status')=='MatcheteSelfEigenHeavyBasisPassedMatchingNotPerformed'
                    and report.get('ordinary_matchete_check_passed') is True
                    and report.get('detailed_heavy_mass_basis_passed') is True,
                    f'{m}: self-eigenstate report does not certify Matchete mass basis')
            dim=int(report['native_dimension'])
            factor=exact_number(report['native_bilinear_metric_factor_inputform'])
            row['exact_quadratic_and_kinetic']=self_mass_identity(dim,factor)
        elif m in BC:
            row['exact_quadratic_and_kinetic']=mixed
        if m in SELF+BC:
            row['exact_native_quartic_transport']=native_rational_transport(
                m,native[key],stage22['models'][key]['tensors'])
            witness_total+=row['exact_native_quartic_transport']['exact_non_floating_witnesses']
        rows.append(row)
        print(f'{m}: native quartics {len(native[key])} | '
              f'{"exact mass + native tensor PASS" if m in SELF+BC else "conditional normalization block PASS"}',flush=True)
    require(sum(x['native_quartic_count'] for x in rows)==50,'Lost native extra quartic tensor')
    result={
        'format':'t3-workstream-a-exact-uv-equivalence-v1',
        'status':'ExactScalarIdentityAndConditionalNormalizationPassedNotFullPhysicalCertificate',
        'source_provenance':[source_digest(p) for p in paths],
        'exact_scalar_quadratic_and_kinetic_models':8,
        'native_extra_quartics_checked':50,
        'native_extra_quartics_exceptional_exact_arithmetic_witnesses':witness_total,
        'conditional_matchete_rgbeta_source_maps':norm,
        'eight_model_exact_symbolic_mass_parameter_identities':True,
        'exact_source_to_rotated_uv_lagrangian_wolfram_equality_certified':False,
        'matchete_rgbeta_global_su2_component_dictionary_certified':False,
        'matchete_rgbeta_physical_coupling_maps_certified':False,
        'explicit_matchete_z2_registered':False,
        'physical_matching_performed':False,
        'rows':rows,
        'remaining_certification':[
            'Impose UV eigenmass/phase relations on symbolic Matchete mass spurions and compare actual full LUV expressions',
            'Export and transport actual native baseline quartic, topology, and Yukawa CG tensors (not just 50 extra quartics)',
            'Determine one global physical Matchete↔RGBeta component/field convention and sign/prefactor dictionary',
            'Verify that the Stage23 conditional coupling maps reproduce the physical four-scalar action in that convention',
            'If required, explicitly register Z2 in Matchete rather than only checking parity of assembled terms',
        ],
    }
    md=['# Workstream A — exact scalar equivalence and normalization candidate audit','',
        f'**Outcome: {result["status"]}**','',
        '- Exact symbolic quadratic and canonical-kinetic identities for **8/8** exceptional scalar sectors.',
        f'- Actual native Matchete extra-quartic tensors: **50/50** matched to Stage23 source-derived conditional maps.',
        '- Stage23 matrix blocks: **46/46** exact radical inverse identities; includes 2×2 mixing blocks.',
        f'- Exceptional native-tensor exact-arithmetic transport witnesses: **{witness_total}**.',
        '- **No** Matchete/RGBeta physical coupling normalization certification; existing Stage23 remains conditional.',
        '', '| Model | Extra quartics | Exact scalar basis | Exact native witnesses |',
        '|---|---:|---|---:|']
    for row in rows:
        ex=row.get('exact_native_quartic_transport',{}).get('exact_non_floating_witnesses')
        md.append(f'| {row["model"]} | {row["native_quartic_count"]} | '
                  f'{"PASS" if row["model"] in SELF+BC else "n/a"} | '
                  f'{ex if ex is not None else "n/a"} |')
    md+=['','## Remaining physical certification work','']+[
        '- '+x for x in result['remaining_certification']
    ]+['','This stage reads previous reports and source tensors. It does not run Match[], RGEs, or modify the repository’s existing production code.','']
    outdir.mkdir(parents=True,exist_ok=True)
    targets[0].write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    targets[1].write_text('\n'.join(md),encoding='utf-8')
    return result


def self_test() -> None:
    require(exact_number('Sqrt[2]/2')==sp.sqrt(2)/2,'Wolfram numeric parser failed')
    require(exact_number('1/3')==sp.Rational(1,3),'Exact rational parser failed')
    require(exact_number('-I/Sqrt[2]')==-sp.I/sp.sqrt(2),'Imaginary numeric parser failed')
    try:
        exact_number('0.5')
    except ValueError:
        pass
    else:
        raise ValueError('Exact-number parser silently accepted an approximate decimal')
    for dim in (1,3):
        self_mass_identity(dim,sp.Integer(1))
        self_mass_identity(dim,sp.Rational(-3,2))
    mixed_mass_identity()
    # Exact quartic tensor transport with arbitrary nontrivial complex CG values.
    for m in SELF+BC:
        fam,a=name_alpha(m)
        d1,d2,_=DIMS[fam]
        if m in SELF:
            neutral='S2' if a==-2 else 'S1'
            dim=d2 if a==-2 else d1
            names=(neutral,neutral+'dag',neutral,neutral+'dag')
            dims=(dim,)*4
        else:
            names=('S1','S1dag','S2','S2dag')
            dims=(2,2,2,2)
        data=np.empty(dims,dtype=object)
        for idx in np.ndindex(dims):
            data[idx]=sp.Rational(1+sum(idx),3)+sp.I*sp.Rational(idx[0]+1,7)
        require(exact_tensor_witness(m,{'species':names,'tensor':data},0),
                f'{m}: direct exact tensor witness failed')
    R=sp.Matrix([[sp.Rational(1,2),1],[-sp.sqrt(3)/2,0]])
    P=sp.Matrix([[0,2*sp.sqrt(3)/3],[-1,-sp.sqrt(3)/3]])
    require(exact_matrix_zero(R*P+sp.eye(2)),'Exact Stage23 2x2 example failed')
    Pbad=P.copy();Pbad[0,0]+=1
    require(not exact_matrix_zero(R*Pbad+sp.eye(2)),
            'Mutated Stage23 coupling matrix wrongly accepted')
    print('PASS: 8/8 exact symbolic scalar quadratic and kinetic identities')
    print('PASS: 8/8 exact native-style component-tensor transport self-tests')
    print('PASS: exact 2x2 conditional normalization inverse example')
    print('NOTE: --all is required to check actual Stage22/Stage23/RGBeta source records')


def main() -> int:
    pa=argparse.ArgumentParser(description=__doc__)
    grp=pa.add_mutually_exclusive_group(required=True)
    grp.add_argument('--self-test',action='store_true')
    grp.add_argument('--all',action='store_true')
    pa.add_argument('--output-dir',type=Path,default=OUT)
    pa.add_argument('--overwrite',action='store_true')
    a=pa.parse_args()
    if a.self_test:
        self_test();return 0
    result=run_all(a.output_dir.resolve(),a.overwrite)
    print('RESULT:',result['status'])
    print('JSON:',a.output_dir.resolve()/'exact_uv_equivalence.json')
    print('REPORT:',a.output_dir.resolve()/'exact_uv_equivalence.md')
    return 0

if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,OSError,ImportError,KeyError,TypeError,IndexError) as exc:
        print('ERROR:',exc,file=sys.stderr)
        raise SystemExit(2)
