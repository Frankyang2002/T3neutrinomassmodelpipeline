"""Auditable all-16 T3 operator/contribution ledger, fail closed.

This inventory separates CODE IMPLEMENTATION from PHYSICAL CERTIFICATION.
GitHub read-only snapshot: 2c29e5ba2c0350b7e9be292e02ceb8bc7c52c0d3,
2026-11-02. A checked-in module's presence does not certify its source mapping.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha1
from pathlib import Path
from typing import Any

MODEL_ALPHAS={'A':(-4,-2,0),'B':(-3,-1,1),'C':(-3,-1,1),
              'D':(-2,0,2),'E':(-4,-2,0,2)}
ALL_MODELS=tuple(f'{c}_{x:+d}' if x>=0 else f'{c}_{x}'
                 for c,vs in MODEL_ALPHAS.items() for x in vs)
GIT_COMMIT='2c29e5ba2c0350b7e9be292e02ceb8bc7c52c0d3'
SEQ='UV -> ALL heavy BSM fermions together at muF -> EFT1(H,S1,S2,SM fermions) -> BOTH heavy scalar multiplets together at muS -> SMEFT'
STATUSES={'INCLUDED','ZERO_BY_SYMMETRY','HIGHER_ORDER','MISSING','UNRESOLVED'}
SOURCES={
    'tree':'scripts/RunGeneralT3FermionFirstEFT1Tree.py',
    'tree_report':'Numerical/joint_scalar_matching/fermion_boundary_bridge.py',
    'seven':'RGE/running/intermediate/GeneralEFT1C5Mixing.py',
    'seven_source':'RGE/running/intermediate/SourceBackedEFT1Batch.py',
    'bubble':'Numerical/joint_scalar_matching/core.py',
    'bubble_adapter':'Numerical/joint_scalar_matching/adapter.py',
    'direct':'Numerical/direct_one_loop/core.py',
    'hard_raw':'scripts/RunGeneralT3EFT1HardPoleAudit.py',
    'hard_matchete':'scripts/RunGeneralT3FermionFirstEFT1OneLoop.py',
    'matched_scalar_bridge':'Numerical/joint_scalar_matching/fermion_boundary_bridge.py',
    'uv':'Lagrangian/interactions/T3Topology.wl',
    'uv_builder':'Lagrangian/model/LagrangianBuilder.wl',
    'uv_catalog':'Lagrangian/model/T3ModelCatalog.wl',
    'quartic':'Lagrangian/interactions/GeneralScalarOperators.wl',
}
# Git blob IDs fetched via the GitHub connector on 2026-11-02. Only paths
# actually verified are pinned; no fabricated hash for the rest.
BLOBS={
 'scripts/RunGeneralT3FermionFirstEFT1Tree.py':'44d9d1dcd1dbf3d1a04b2d10680e60487a4004fb',
 'scripts/RunGeneralT3FermionFirstEFT1OneLoop.py':'3df9f231de5df211d3bc30e0df472291f8e8815a',
 'scripts/RunGeneralT3EFT1HardPoleAudit.py':'5a1267f07c5aa7dc4d6f13412573374200917af9',
 'Numerical/joint_scalar_matching/core.py':'b2afc8c71b0adc0d0783995e0ea65513e622b712',
 'Numerical/joint_scalar_matching/adapter.py':'7a8653772c9e8725744599aefc216695b5d34944',
 'Numerical/joint_scalar_matching/fermion_boundary_bridge.py':'6929de2aabbd2e69b382db239f558add06059e69',
 'Numerical/direct_one_loop/core.py':'abbf261fbf2ca2561482d3e53b84452d1654610e',
 'Lagrangian/interactions/T3Topology.wl':'a45a099a8b766cbad0114a0559cb0bcb80895649',
 'Lagrangian/interactions/ScalarPotential.wl':'58f5bbef2688f9e87a1b8992bca11076ee254437',
 'Lagrangian/interactions/GeneralScalarOperators.wl':'7b5cda73830be1717eb543c6c5d1d99a3dec1dd3',
 'Lagrangian/model/LagrangianBuilder.wl':'98d549cda503b751380dc10c4a56212e5312ccd8',
 'Lagrangian/model/T3ModelCatalog.wl':'9a6c47e10fc83d9480f08353682151180c9fbb4a',
 'Lagrangian/model/T3Fields.wl':'65467ede6ef06a7743e3ca2878debf777dc94507',
}

# Status means "this contribution is represented by code or a mathematical
# proof"; the independent physical_certified field is ALWAYS false unless
# determined from explicitly sourced certificates. We never auto-promote.
ENTRIES=(
 ('EXTRA_NEUTRAL_SCALAR_CONJUGATE_YUKAWA','UNRESOLVED','uv','Additional independent Yukawa exists ONLY if scalar S1 or S2 has Y=0; construction is opt-in via T3_EXPERIMENTAL_GENERAL_FERMIONS, so its contribution to LLSS and one-loop mass must be audited.'),
 ('GENERAL_SCALAR_PORTAL_CHANNEL_CLOSURE','UNRESOLVED','quartic','Extra renormalizable scalar quartics are gated behind environment flags; no justified default zero for symmetry-allowed HHSS couplings.'),
 ('GENERAL_HEAVY_FERMION_FLAVOUR_MASS','UNRESOLVED','uv_builder','Builder explicitly does not certify general off-diagonal flavour mass matrices.'),
 ('TREE_LLSS','INCLUDED','tree','All heavy fermions eliminated at a common muF; formal LLSS tree boundary exists, signed master-Weyl normalization unresolved.'),
 ('TREE_LLHH','ZERO_BY_SYMMETRY','tree','Exact inert Z2 forbids a lepton-Higgs-BSM-fermion Yukawa; no direct Weinberg tree seesaw in the specified T3 UV.'),
 ('ODD_SCALAR_DIM5_LLHS','ZERO_BY_SYMMETRY','uv','One Z2-odd scalar with SM-even external leptons/Higgs is Z2-odd.'),
 ('TREE_LLSS_ALL_SCALAR_PAIRS','UNRESOLVED','tree_report','All saved PL sources retained, but independent completeness of every symmetry-allowed LLSS scalar/SU2 structure is not certified.'),
 ('F_THRESHOLD_HARD_LLHH','MISSING','hard_raw','One-loop hard Weinberg term is allowed, generally nonzero; B_-1/C_-1 raw hard poles/logs exist but no fully subtracted, signed, projected C5 for all 16.'),
 ('F_THRESHOLD_HARD_SCALAR_TRIANGLE_DIM5','UNRESOLVED','bubble','Conditional per-fermion mass eigenstate universal hard kernel added by this batch; physical per-species LLSS tensors/normalization not supplied.'),
 ('F_THRESHOLD_OTHER_HARD_TOPOLOGIES','UNRESOLVED','hard_matchete','Exhaustive UV one-loop operator/counterterm reduction including other Yukawa and gauge topologies not certified.'),
 ('F_THRESHOLD_ONE_LOOP_LLSS','HIGHER_ORDER','tree','One-loop correction to LLSS needs another scalar loop to produce C5 when inert scalars have no VEV: two loops in neutrino mass.'),
 ('EFT1_SEVEN_TERM_C5_MIXING','INCLUDED','seven','Seven-term one-loop master Wilson source implemented; numerical full source/coupling dictionary and signed normalization remain conditional.'),
 ('EFT1_SCALAR_C5_LOG','INCLUDED','seven','A scalar-pair one-loop C_LLSS*lambda contribution; only once in EFT1 logarithmic evolution.'),
 ('EFT1_LLSS_ONE_LOOP_SELF_RUNNING_IN_BUBBLE','HIGHER_ORDER','seven','A one-loop beta correction to tree LLSS inserted into an already one-loop bubble is two loops.'),
 ('JOINT_SCALAR_FINITE_BUBBLE','INCLUDED','bubble','All eigenstate pairs, both scalar multiplets integrated together at one muS; finite B0bar and complex flavour preserved.'),
 ('JOINT_SCALAR_OTHER_OPERATOR_INSERTIONS','UNRESOLVED','bubble','Finite matching of all possible independent EFT1 dimension-five sources and field/counterterm effects not fully certified.'),
 ('DIRECT_UV_ONE_LOOP_P0','INCLUDED','direct','Independent unexpanded one-loop Weyl/real-scalar mass eigenstate kernel; source vertex mapping provisional.'),
 ('ONE_LOOP_HIGHER_DIMENSION_POWER_TERMS','UNRESOLVED','direct','Order mS^2/MF^2 times one-loop is NOT higher-loop; exact triangle kernel available conditionally, but full higher-dimensional EFT basis uncomputed.'),
 ('SMEFT_RUNNING_ON_ONE_LOOP_C5','HIGHER_ORDER','seven','One-loop SMEFT running applied to an O(hbar) Weinberg boundary changes fixed-order mass at O(hbar^2).'),
 ('DOUBLE_COUNT_SCALAR_BUBBLE_LOG','UNRESOLVED','bubble','The EFT1 scalar beta must cancel explicit muS derivative of the joint bubble; certification requires common normalization and all RGE terms.'),
 ('GAUGE_FLAVOUR_COUNTERTERM_WARD','UNRESOLVED','hard_matchete','Wavefunction and scheme/counterterm cancellations are not certified for all sixteen physical inputs.'),
)


def _git_blob(path: Path) -> str | None:
    if not path.is_file(): return None
    content=path.read_bytes()
    return sha1(f'blob {len(content)}\0'.encode()+content).hexdigest()


def source_audit(root: Path) -> dict[str,dict[str,Any]]:
    result={}
    for label,path in SOURCES.items():
        expect=BLOBS.get(path)
        actual=_git_blob(Path(root)/path)
        result[label]={'path':path,'present_locally':actual is not None,
                       'local_git_blob_sha':actual,'github_snapshot_blob_sha':expect,
                       'snapshot_agreement':None if expect is None or actual is None else expect==actual}
    return result


def build_ledger(root: Path | None = None) -> dict[str,Any]:
    audit=source_audit(root) if root is not None else {}
    rows=[]
    for model in ALL_MODELS:
        alpha=int(model.split('_')[1])
        ys1=alpha/2
        ys2=(alpha+2)/2
        for code,status,source,reason in ENTRIES:
            if code=='EXTRA_NEUTRAL_SCALAR_CONJUGATE_YUKAWA' and ys1 != 0 and ys2 != 0:
                status='ZERO_BY_SYMMETRY'
                reason='Hypercharge forbids the opposite-complex-scalar extra Yukawa for this model: neither scalar has Y=0.'
            rec={'model':model,'contribution':code,'status':status,
                 'source_path':SOURCES[source],
                 'quantum_numbers':{'Y_S1':ys1,'Y_S2':ys2,'Y_F':(alpha+1)/2},
                 'basis':'full complex 3x3; both scalar multiplets',
                 'source_supported_reason':reason,'physical_certified':False,
                 'implementation_state': ('conditional_numerical_hard_kernel_new' if code=='F_THRESHOLD_HARD_SCALAR_TRIANGLE_DIM5'
                     else 'source_algorithm_present' if status=='INCLUDED'
                     else 'not_newly_implemented'),
                 'external_blockers':(['signed GammaCC-to-Weyl normalization','full UV/EFT1 operator closure']
                                      if status in ('UNRESOLVED','MISSING') else [])}
            if source in audit: rec['local_source_audit']=audit[source]
            # Only the two models have a verified hard-LLHH *raw factorisation*
            # path. Never infer that physical renormalized C5 is obtained.
            if code=='F_THRESHOLD_HARD_LLHH':
                rec['raw_hard_pole_log_factorization_available']=model in ('B_-1','C_-1')
            rows.append(rec)
    counts=Counter(x['status'] for x in rows)
    return {'format':'t3-leading-one-loop-all16-contribution-ledger-v1',
            'repository_commit_readonly_snapshot':GIT_COMMIT,
            'eft_sequence':SEQ,'models':list(ALL_MODELS),'n_models':len(ALL_MODELS),
            'n_contributions_per_model':len(ENTRIES),'n_rows':len(rows),
            'status_counts':dict(sorted(counts.items())),
            'fully_physics_complete':False,'physical_C5_certified':False,
            'different_loop_order_not_inferred_from_small_coefficient':True,
            'notes':['Included means algorithmically present, never independently physical-certified.',
                     'Unresolved does not imply zero.',
                     'Exact one-loop mS^2/MF^2-suppressed terms are retained as separately tagged triangle remainder, not called two loops.',
                     'A full one-loop neutrino mass is not the same as a dimension-five leading-power EFT truncation.'],
            'sources':audit,'rows':rows}
