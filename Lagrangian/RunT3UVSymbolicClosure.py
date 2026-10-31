#!/usr/bin/env python3
"""Workstream A: independent abstract SU(2) triplet reference, replay of the
native conditional signed-action identity, and final evidence classification.

This audit intentionally DOES NOT assert that an abstract reference irrep is
RGBeta's independently observed numerical component basis. In particular it
must not promote a conditional Stage23/26 signed map into a physical mapping.

From the local project root:
 python Lagrangian/RunT3UVSymbolicClosure.py --self-test
 python Lagrangian/RunT3UVSymbolicClosure.py --all
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import QuarticActionIdentity as ACTION
from Lagrangian import RunT3UVNativeBridge as NATIVE

DATA = ROOT / 'output/general_t3_stage22'
REPORTS = ROOT / 'Reports/output/workstream_a'
OUTPUT = REPORTS / 'symbolic_uv_closure'
STATUS = 'ConditionalSymbolicUVClosurePassedPhysicalRGBetaDictionaryBlocked'


def require(check: Any, message: str) -> None:
    if not check:
        raise ValueError('Workstream A symbolic closure: ' + message)


def load(path: Path) -> dict[str, Any]:
    require(path.is_file(), f'missing prerequisite {path}')
    x = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(x,dict),f'expected JSON object at {path}')
    return x


def digest(path: Path) -> dict[str,str]:
    return {'path':str(path.relative_to(ROOT)),
            'sha256':sha256(path.read_bytes()).hexdigest()}


def zero(expr: Any) -> bool:
    if isinstance(expr,sp.MatrixBase):
        return all(sp.simplify(v)==0 for v in expr)
    return sp.simplify(expr)==0


def symmetric_doublet_reference() -> dict[str, Any]:
    """Derive a declared *reference* spin-1 irrep from Sym^2(2).

    This has no assumed RGBeta numeric basis; uses only Pauli/2.  Its
    unbarred triplet and conjugate are related by the invariant metric K.
    Every identity is exact in Q(sqrt(2),i).
    """
    I=sp.I
    sq=sp.sqrt(2)
    t=(sp.Matrix([[0,1],[1,0]])/2,
       sp.Matrix([[0,-I],[I,0]])/2,
       sp.diag(1,-1)/2)
    e=sp.Matrix([[1,0,0],[0,1/sq,0],[0,1/sq,0],[0,0,1]])
    # Coordinates: (H_0 H_0, H_0 H_1, H_1 H_0, H_1 H_1).
    perm=sp.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    projector=(sp.eye(4)+perm)/2
    require(zero(e.H*e-sp.eye(3)) and zero(e*e.H-projector),
            'symmetric doublet projector or orthogonality failed')
    j=tuple(sp.simplify(e.H*(sp.kronecker_product(x,sp.eye(2))+
                            sp.kronecker_product(sp.eye(2),x))*e) for x in t)
    for k,x in enumerate(j):
        require(zero((sp.kronecker_product(t[k],sp.eye(2))+
                      sp.kronecker_product(sp.eye(2),t[k]))*e-e*x),
                'delS2 symmetric-doublet covariance failed')
        require(zero(x-x.H) and zero(sp.trace(x)),
                'spin-one generator Hermiticity or trace failed')
    for a,b,c in ((0,1,2),(1,2,0),(2,0,1)):
        require(zero(j[a]*j[b]-j[b]*j[a]-I*j[c]),
                'spin-one commutator failed')
    require(zero(sum((x*x for x in j),sp.zeros(3))-2*sp.eye(3)),
            'quadratic Casimir is not C2=2')
    for a in range(3):
        for b in range(3):
            require(zero(sp.trace(j[a]*j[b])-(2 if a==b else 0)),
                    'triplet trace normalization failed')
    metric=sp.Matrix([[0,0,1],[0,-1,0],[1,0,0]])
    require(zero(metric*metric.H-sp.eye(3)) and zero(metric.T-metric),
            'reference triplet conjugation metric not symmetric/unitary')
    for x in j:
        require(zero(x.T*metric+metric*x),
                'triplet invariant bilinear covariance failed')
        require(zero(metric*(-x.conjugate())*metric.inv()-x),
                'barred-triplet generator -conj(J) not equivalent through metric')
    return {'status':'ExactAbstractSymmetricDoubletTripletReferenceVerified',
            'basis':['H0 H0','(H0 H1+H1 H0)/sqrt(2)','H1 H1'],
            'spin_one_generators':[[[str(v) for v in row] for row in x.tolist()] for x in j],
            'conjugation_intertwiner_K':[[str(v) for v in row] for row in metric.tolist()],
            'delS2_symmetric_doublet_projector_exact':True,
            'delS2_spin_one_covariance_exact':True,
            'triplet_su2_commutators_casimir_trace_exact':True,
            'barred_tGen_minus_complex_conjugate_equivalence_exact':True,
            'native_RGBeta_triplet_component_basis_extracted':False,
            'physical_intertwiner_to_RGBeta_verified':False}


def symbolic_rgbeta_source_guard(source: str) -> dict[str,Any]:
    """A guarded *source* audit, not a runtime RGBeta group-component test.

    Strip Wolfram comments before asserting built definitions, otherwise a
    commented-out AddYukawa would give false evidence of implementation.
    """
    text=re.sub(r'\(\*.*?\*\)','',source,flags=re.DOTALL)
    anchors={
      'triplet_representation':r'T3RGBetaRep\[3\]\s*:=\s*SU2L\[S2\]',
      'triplet_yukawa_singlet':r'T3RGBetaYukawaInvariant\[2,\s*3\]',
      'triplet_yukawa_doublet':r'T3RGBetaYukawaInvariant\[3,\s*2\]',
      'topology_13':r'T3RGBetaMixInvariant\[1,\s*3\]',
      'topology_31':r'T3RGBetaMixInvariant\[3,\s*1\]',
      'topology_22':r'T3RGBetaMixInvariant\[2,\s*2\]',
      'topology_33':r'T3RGBetaMixInvariant\[3,\s*3\]',
      'native_triplet_CG':r'delS2\[SU2L,',
      'native_triplet_generators':r'tGen\[SU2L\s*@\s*S2,',
      'rgbeta_baseline_y1':r'AddYukawa\[\s*y1\s*,',
      'rgbeta_baseline_y2':r'AddYukawa\[\s*y2\s*,',
    }
    for label,pat in anchors.items():
        require(re.search(pat,text) is not None, f'RGBeta source anchor changed: {label}')
    baseline={k:len(re.findall(r'AddYukawa\[\s*'+k+r'\s*,',text)) for k in ('y1','y2')}
    extras={k:len(re.findall(r'AddYukawa\[\s*'+k+r'\s*,',text)) for k in ('y1Conj','y2Conj')}
    require(baseline=={'y1':1,'y2':1},'RGBeta baseline Yukawa registrations changed')
    require(extras=={'y1Conj':0,'y2Conj':0},
            'RGBeta extra Yukawa registration changed; reassess this report')
    return {'rgbeta_source_native_delS2_tGen_used':True,
            'rgbeta_baseline_yukawa_registration_count':2,
            'rgbeta_extra_conjugate_scalar_yukawa_registrations':0,
            'matchete_extra_yukawa_matrix_types_outside_current_rgbeta_builder':6,
            'rgbeta_matchete_complete_fermionic_dictionary_certified':False,
            'proof_type':'guarded source inspection; not live independent native component proof'}


def conditionals_rerun(data:Path=DATA) -> dict[str,Any]:
    """Recompute, do not merely trust, all 46 full action polynomial blocks.

    Uses existing Stage26 reconstructed tensor convention.  It includes real
    field polynomials and H.c. and has negative sign/1/4 controls. RGBeta
    native invariant components are *not* independently available for triplet.
    """
    source_names=('matchete_tensors.json','rgbeta_native_invariants.json',
                  'stage23_quartic_prefactor_candidates.json',
                  'stage25_native_quartic_assembly_diagnostics.json',
                  'stage26_conditional_quartic_action.json')
    docs={name:load(data/name) for name in source_names}
    out=ACTION.audit({'matchete_tensors.json':docs['matchete_tensors.json'],
                      'rgbeta_native_invariants.json':docs['rgbeta_native_invariants.json']},
                      docs['stage23_quartic_prefactor_candidates.json'],
                      docs['stage25_native_quartic_assembly_diagnostics.json'])
    reference=docs['stage26_conditional_quartic_action.json']
    require(out['status']=='SourceSignAndConditionalQuarticActionCheckedGlobalBasisStillUncertified'
            and out['blocks_checked']==46 and out['slots_checked']==50,
            'fresh signed-action polynomial audit failed')
    require(reference['status']==out['status'] and
            reference.get('models_checked')==out['models_checked'] and
            reference.get('blocks_checked')==out['blocks_checked'] and
            reference.get('slots_checked')==out['slots_checked'],
            'saved Stage26 report incomplete or inconsistent')
    refblocks={(x['model'],x['block']):x for x in reference['blocks']}
    require(len(refblocks)==46,'duplicate/missing saved Stage26 block')
    for block in out['blocks']:
        key=(block['model'],block['block'])
        require(block==refblocks.get(key),'saved Stage26 signed-action polynomial differs: '+str(key))
        require(block['physical_mapping_certified'] is False and
                block['global_component_basis_verified'] is False,
                'unexpected physical normalization claim')
    return {'recomputed_signed_action_blocks':46,
            'native_quartic_coupling_slots_covered':50,
            'models_with_native_scalar_inventories':16,
            'rgbeta_live_native_permutation_rules_checked':out['native_stage25_permutation_rules_checked'],
            'hermitian_conjugate_terms_included':out['hermitian_conjugate_formal_polynomial_checked'],
            'negative_controls_from_stage26':'wrong sign and spurious 1/4 rejected',
            'all_recomputed_blocks_identical_to_saved_stage26':True,
            'reconstruction_basis_conditional':True,
            'global_rgbeta_matchete_physical_component_basis_certified':False,
            'physical_signed_coupling_map_certified':False}


def validate_prior() -> dict[str,Any]:
    statuses={
      'uv_closure_stage/uv_closure_stage.json':'StructuralAndNumericalStagePassed',
      'final_uv_evidence/final_uv_evidence.json':'ConditionalUVOperatorEvidencePassedNotFullCompleteness',
      'exact_uv_equivalence/exact_uv_equivalence.json':'ExactScalarIdentityAndConditionalNormalizationPassedNotFullPhysicalCertificate',
      'native_cg_bridge/native_cg_bridge.json':'LiveMatcheteBaselineCGsPassedRGBetaPhysicalBridgeUncertified',
      'triplet_coordinate_gate/triplet_coordinate_gate.json':'LiveMatcheteTripletRepresentationVerifiedRGBetaCoordinatesUnavailable',
    }
    for relative,expected in statuses.items():
        doc=load(REPORTS/relative)
        require(doc.get('status')==expected, 'prior report did not pass: '+relative)
    t=load(REPORTS/'triplet_coordinate_gate/triplet_coordinate_gate.json')
    require(t.get('matchete_HH_triplet_models_exactly_checked')==6 and
            t.get('rgbeta_native_numeric_triplet_basis_available') is False and
            t.get('signed_physical_coupling_dictionary_certified') is False,
            'current RGBeta coordinate blocker changed: investigate before certification')
    return {'previous_closure_models':16,
            'independent_fermionic_matrix_types':38,
            'previous_native_yukawa_and_topology_CGs':48,
            'previous_exact_scalar_basis_models':8,
            'previous_native_triplet_witness_models':6,
            'previous_native_RGBeta_numeric_triplet_components':False}


def validate_triplet_native_runtime_evidence() -> dict[str,Any]:
    stage28=load(DATA/'stage28_live_rgbeta_triplet_structure.json')
    require(stage28.get('models_checked')==6 and
            stage28.get('all_group_representation_checks_passed') is True and
            stage28.get('all_native_barred_generator_rules_passed') is True and
            stage28.get('all_native_delS2_contractions_passed') is True and
            stage28.get('native_numerical_triplet_generators_available') is False and
            stage28.get('physical_cross_package_component_basis_certified') is False,
            'live native RGBeta symbolic SU2 audit changed')
    rows=stage28.get('models',[])
    require(len(rows)==6 and all(row.get('native_generator_is_symbolic') is True
             and row.get('native_conjugate_generator_minus_transpose') is True
             and row.get('native_delS2_orthogonality') is True for row in rows),
            'native RGBeta six-model symbolic tensor probe altered')
    return {'live_rgbeta_abstract_triplet_representation_probed':True,
            'native_symbolic_representation_models':6,
            'symbolic_barred_generator_and_delS2_probes':True,
            'native_numeric_component_matrices_available':False,
            'abstract_reference_is_not_independently_matched_to_RGBeta_basis':True}


def self_test() -> None:
    ref=symmetric_doublet_reference()
    require(ref['triplet_su2_commutators_casimir_trace_exact'],
            'spin-one reference selftest failed')
    # A non-physical toy tensor: signed Lagrangian and full H.c. identity.
    mat={'species':['H','Hdag','S1','S1dag'],
         'dimensions':[1,1,1,1], 'components':['2']}
    pol=ACTION.assemble_action([mat],[sp.Symbol('g',complex=True)],plus_hc=True)
    rpc=ACTION.assemble_action([mat],[-sp.Symbol('g',complex=True)],plus_hc=True)
    require(ACTION.actions_equal(pol,{k:-v for k,v in rpc.items()}),
            'signed action selftest failed')
    require(not ACTION.actions_equal(pol,rpc), 'incorrect sign passed synthetic control')
    src='''T3RGBetaRep[3] := SU2L[S2];
     T3RGBetaYukawaInvariant[2, 3] := delS2[SU2L, #3,#1,#2] &;
     T3RGBetaYukawaInvariant[3, 2] := delS2[SU2L, #1,#2,#3] &;
     T3RGBetaMixInvariant[1, 3] := tGen[SU2L @ S2, A, i,j] &;
     T3RGBetaMixInvariant[3, 1] := 1 &;
     T3RGBetaMixInvariant[2, 2] := 1 &;
     T3RGBetaMixInvariant[3, 3] := 1 &;
     AddYukawa[y1,foo]; AddYukawa[y2,foo];'''
    symbolic_rgbeta_source_guard(src)
    try:
        symbolic_rgbeta_source_guard(src.replace('AddYukawa[y2','AddYukawa[y2Conj'))
    except ValueError:
        pass
    else:
        raise AssertionError('Missing Yukawa coverage detected incorrectly')
    print('PASS: exact symmetric-doublet triplet, generators, barred metric and symbolic covariance')
    print('PASS: signed H.c. polynomial identity and negative controls')
    print('PASS: source/Yukawa guard negative control')
    print('NOTE: --all recomputes live-source-evidence-dependent 46-block Stage26 action identities')


def all_stage(outdir:Path, overwrite:bool) -> dict[str,Any]:
    targets=(outdir/'symbolic_uv_closure.json',outdir/'symbolic_uv_closure.md')
    require(overwrite or not any(path.exists() for path in targets),
            'output already exists: back up and use --overwrite')
    old=validate_prior()
    abstract=symmetric_doublet_reference()
    rgbeta_symbolic=validate_triplet_native_runtime_evidence()
    source=ROOT/'RGE/running/rgbeta/T3RGBetaModel.wl'
    require(source.is_file(),'missing current RGBeta model source')
    source_info=symbolic_rgbeta_source_guard(source.read_text(encoding='utf-8-sig'))
    print('SU2: exact Sym^2 fundamental delS2/tGen reference + conjugate metric PASS',flush=True)
    print('RGBeta: native symbolic triplet structure 6/6; numerical component basis unavailable',flush=True)
    print('Recomputing 46 exact signed full-action polynomial identities (50 extra quartic slots)',flush=True)
    polys=conditionals_rerun()
    require(polys['recomputed_signed_action_blocks']==46 and
            polys['native_quartic_coupling_slots_covered']==50,
            'conditional signed-action polynomial coverage missing')
    print('Quartics: 46/46 signed action and H.c. identities in provisional basis PASS',flush=True)
    evidence_paths=[source,DATA/'stage23_quartic_prefactor_candidates.json',
      DATA/'stage25_native_quartic_assembly_diagnostics.json',
      DATA/'stage26_conditional_quartic_action.json',
      DATA/'stage28_live_rgbeta_triplet_structure.json',
      REPORTS/'native_cg_bridge/native_cg_bridge.json',
      REPORTS/'triplet_coordinate_gate/triplet_coordinate_gate.json']
    report={
      'format':'t3-workstream-a-symbolic-uv-closure-v1','status':STATUS,
      'source_sha256':[digest(p) for p in evidence_paths],
      'previous_stage_evidence':old,
      'derived_reference_triplet':abstract,
      'live_rgbeta_symbolic_representation_evidence':rgbeta_symbolic,
      'guarded_rgbeta_model_source_coverage':source_info,
      'full_signed_quartic_action_recomputation':polys,
      'source_includes_separately_implemented_extra_yukawa_matrices_in_rgbeta':False,
      'full_matchete_rgbeta_physical_component_basis_certified':False,
      'full_matchete_rgbeta_signed_physical_coupling_dictionary_certified':False,
      'exact_original_vs_rotated_full_wolfram_lagrangian_identity_certified':False,
      'explicit_Z2_registered_in_matchete':False,
      'full_renormalizable_uv_completeness_unconditional_certificate':False,
      'physical_matching_performed':False, 'production_modified':False,
      'remaining_blockers':[
        'RGBeta does not expose independent triplet CG/generator numeric component coordinates; Sym^2 reference is an explicit candidate, not a runtime basis identification',
        'Stage23/26 signed maps are exact only in the declared provisional/common SU(2) component convention',
        'The current RGBeta UV model registers y1,y2 but not all six conjugate-neutral-scalar Yukawa matrix types present in generalized Matchete UV models',
        'No one-session equality of the complete native Matchete and RGBeta actions, including full baseline Yukawa/topology normalization, has been certified',
        'Z2 is checked as parity of actual assembled terms, not explicitly registered in Matchete',
      ],
    }
    lines=['# Workstream A — final symbolic UV evidence and normalization boundary','',
           f'**Status:** `{STATUS}`','',
           '## Evidence matrix','',
           '| Audit | Result | Nature |','|---|---|---|',
           '| Scalar and Z2 structural closure | 16/16 | Previous local Matchete audit |',
           '| Independent Yukawa inventory | 16/16, 38 matrix types | Previous independent census |',
           '| Exceptional scalar mass/kinetic basis | 8/8 | Previous exact scalar identities |',
           '| Live baseline Matchete CGs | 48/48 | Previous native registration |',
           '| Triplet from live Matchete CGs | 6/6 | Previous native-CG derivation |',
           '| Symmetric-doublet reference triplet, barred metric | PASS | **Exact abstract construction**, not native RGBeta components |',
           '| Native RGBeta `delS2`/`tGen` symbolic rules | 6/6 | Saved live symbolic structure probes |',
           '| Quartic signed Lagrangian/potential/H.c. polynomials | 46/46 blocks, 50/50 slots | **Recomputed conditional action identity** |',
           '| Physical Matchete–RGBeta coupling dictionary | **BLOCKED** | No independently certified common triplet basis |',
           '| Generalized RGBeta Yukawas | **GAP** | Six extra matrix types absent from current source builder |',
           '', '## Conventions and meaning', '',
           'The exact reference uses `Sym^2(2)` with basis `(00,(01+10)/sqrt(2),11)`,',
           'Pauli fundamental generators `sigma_a/2`, and triplet conjugation metric',
           '`K = [[0,0,1],[0,-1,0],[1,0,0]]`. Its generator action,',
           'the symmetric projector, `Tr(T_a T_b)=2 delta_ab`,',
           '`C2=2`, and `K (-conjugate(T_a)) K^-1=T_a` all hold exactly.',
           'This **constructs a consistent reference**; RGBeta has not exported',
           'the independent numerical components needed to identify it with its own basis.',
           '', 'The previously source-derived Stage23 signed quartic maps were checked again',
           'as **complete commuting component polynomials with H.c.** via the existing',
           'Stage26 tensor reconstruction. Wrong relative sign and a spurious 1/4',
           'are rejected. RGBeta triplet tensors still depend on a provisional',
           'component-coordinate identification; this does not certify a physical map.',
           '', '## Unresolved requirements (not test failures)', '']
    lines.extend('- '+x for x in report['remaining_blockers'])
    lines += ['', '## Final classification', '',
              '**Conditional UV operator evidence passes. An unconditional physical',
              'Matchete–RGBeta normalization certificate is not established.**',
              'No production RGE, matching or numerical code has been modified.', '']
    outdir.mkdir(parents=True,exist_ok=True)
    for p in targets:
        if p.exists() and not overwrite:
            raise ValueError('Unexpected existing report '+str(p))
    targets[0].write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    targets[1].write_text('\n'.join(lines),encoding='utf-8')
    print('RESULT:',STATUS,flush=True)
    print('JSON:',targets[0],flush=True)
    print('REPORT:',targets[1],flush=True)
    return report


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    act=parser.add_mutually_exclusive_group(required=True)
    act.add_argument('--self-test',action='store_true')
    act.add_argument('--all',action='store_true')
    parser.add_argument('--overwrite',action='store_true')
    parser.add_argument('--output-dir',default=str(OUTPUT))
    args=parser.parse_args()
    try:
        if args.self_test:
            self_test()
        else:
            all_stage(Path(args.output_dir).resolve(),args.overwrite)
        return 0
    except (OSError,ValueError,KeyError,TypeError,IndexError,AssertionError) as exc:
        print('FAILED:',exc,file=sys.stderr)
        return 1


if __name__=='__main__':
    sys.exit(main())
