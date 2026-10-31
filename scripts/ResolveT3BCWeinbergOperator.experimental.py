#!/usr/bin/env python3
"""B/C native LLHH versus Warsaw operator: installed-source-anchored, fail-closed.

This proves the standard chiral/Grassmann spinor relation in an explicit
representation and audits the actual B/C InputForm orientations.  It does not
assume an installed Matchete model nor infer a physical threshold Wilson.
An optional live Matchete source/operator probe provides an additional check.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from collections import Counter
import sys
import sympy as sp

class AuditError(ValueError):
    pass


def spinor_algebra() -> dict:
    """Exact chiral gamma matrices and antisymmetric Grassmann bivectors."""
    I = sp.I
    z = sp.zeros(2)
    eye = sp.eye(2)
    sigma2 = sp.Matrix([[0, -I], [I, 0]])
    gamma0 = sp.BlockMatrix([[z, eye], [eye, z]]).as_explicit()
    gamma2 = sp.BlockMatrix([[z, sigma2], [-sigma2, z]]).as_explicit()
    C = sp.simplify(I * gamma2 * gamma0)
    pl = sp.diag(1, 1, 0, 0)
    pr = sp.diag(0, 0, 1, 1)
    eye4 = sp.eye(4)
    assert (C.T + C) == sp.zeros(4)
    assert (C*C + eye4) == sp.zeros(4)
    assert (C.conjugate().T + C) == sp.zeros(4)
    assert gamma0*C*gamma0 == -C
    assert gamma0*pl*gamma0 == pr
    assert C*pl == pl*C
    assert C*pr == pr*C

    # Grassmann generators: left-flavour p/r and the conjugate generators.
    # A polynomial term is indexed by its sorted wedge basis, with canonical
    # sign from sorting; ordinary commuting multiplication would get this wrong.
    def product(gens: tuple[str, ...], coef=1):
        if len(set(gens)) != len(gens):
            return {}, 0
        inversions = sum(gens[a] > gens[b]
                         for a in range(len(gens))
                         for b in range(a+1, len(gens)))
        return tuple(sorted(gens)), sp.Integer(coef)*(-1 if inversions%2 else 1)

    def bilinear(u, matrix, v):
        out = {}
        for a in range(2):
            for b in range(2):
                val = sp.simplify(matrix[a,b])
                if val == 0: continue
                key,sgn = product((u[a],v[b]))
                if not key: continue
                out[key] = sp.simplify(out.get(key,0)+val*sgn)
        return {k:v for k,v in out.items() if v != 0}

    eps = C[:2,:2]
    cpr = C[2:,2:]
    fields = {fl: [f"{fl}_{a}" for a in (1,2)] for fl in ("p","r","P","R")}
    W_pr = bilinear(fields['p'],eps,fields['r'])
    W_rp = bilinear(fields['r'],eps,fields['p'])
    if W_pr != W_rp or not W_pr:
        raise AuditError('Majorana bilinear is not flavour symmetric')

    # Hermitian conjugation is an anti-automorphism (XY)^dag = Y^dag X^dag.
    # Uppercase generators are conjugates of lowercase left Weyl components.
    conj={"p":"P","r":"R"}
    def dagger(poly):
        out={}
        for basis,v in poly.items():
            # Reverse ordered Grassmann operators before mapping flavours.
            reversed_fields=tuple(conj[g[0]]+g[1:] for g in reversed(basis))
            key,sgn=product(reversed_fields)
            if not key:continue
            out[key]=sp.simplify(out.get(key,0)+sp.conjugate(v)*sgn)
        return {k:v for k,v in out.items() if v!=0}

    # Native PL = l_p^T C PL l_r; native PR = bar(l_p) C PR bar(l_r)^T.
    # Thus native PR(p,r) = [Warsaw(r,p)]^dag = [Warsaw(p,r)]^dag.
    native_pr = bilinear(fields['P'],cpr,fields['R'])
    if native_pr != dagger(W_pr) or not native_pr:
        raise AuditError('Native PR fails Warsaw Hermitian-conjugate sign test')
    if bilinear(fields['P'],-cpr,fields['R']) == dagger(W_pr):
        raise AuditError('Incorrect opposite sign was not rejected')
    if bilinear(fields['p'],eps,fields['r']) == {}:
        raise AuditError('Unexpected vanishing Majorana operator')

    # Hermitian gamma identity: [l_p^T C PL l_r]^dag =
    # bar(l_r) C PR bar(l_p)^T with no extra minus sign.
    hc_matrix=sp.simplify(gamma0*pl*C.conjugate().T*gamma0)
    if hc_matrix != C*pr:
        raise AuditError('Charge-conjugation Hermitian-conjugate matrix identity failed')
    return {
        'chiral_C_matrix': [[str(C[i,j]) for j in range(4)] for i in range(4)],
        'C_transpose_antisymmetric': True,
        'C_squared_minus_identity': True,
        'C_hermitian_conjugate_minus_C': True,
        'PL_warsaw_equals_native_PL_grassmann': True,
        'PR_native_equals_warsaw_Hc_grassmann': True,
        'relative_spinor_sign': '+1',
        'wrong_opposite_sign_rejected': True,
        'symmetric_flavour_bilinear': True,
    }


def evidence_check(path: Path) -> dict:
    if not path.is_file():
        raise AuditError(f'Installed Matchete source evidence not found: {path}')
    data=json.loads(path.read_text(encoding='utf-8-sig'))
    if (data.get('status')!='InstalledGammaCCSourceLocated' or
        not data.get('source_only') or data.get('files_scanned',0)<1):
        raise AuditError('Unsupported installed-source probe metadata')
    rows=data.get('source_occurrences',[])
    def snippets(path_part, symbol):
        return '\n'.join(sn['excerpt'] for x in rows
          if path_part.lower() in x.get('file','').lower() and x.get('symbol')==symbol
          for sn in x.get('source_snippets',[]))
    dirac=snippets('DiracAlgebra.m','CConj')
    cc=snippets('DiracAlgebra.m','GammaCC')
    proj=snippets('DiracAlgebra.m','Proj')
    warsaw=snippets('SMEFT_Warsaw.m','CConj')
    mass=snippets('Definitions.m','CConj')
    anchors={
      'CConj_field_C_barTranspose': 'DiracProduct@ GammaCC \\[CenterDot] Transp@Bar@f' in dirac,
      'CConj_bar_field_transpose_C': 'Transp@Bar@f \\[CenterDot] DiracProduct@GammaCC' in dirac,
      'C_squared_negative': 'GammaCC, GammaCC, b___]:= - DiracProduct' in cc,
      'PL_makes_minus_one_projector': 'PL= DiracProduct@ Proj@ -1' in proj,
      'PR_makes_plus_one_projector': 'PR= DiracProduct@ Proj@ 1' in proj,
      'Warsaw_Lagd5_uses_BarCConj_without_half':
          'cllHH[p,r] Bar@eps[SU2L][i,j]Bar@eps[SU2L][k,m] H[i] H[k] Bar@CConj@l[j,p]**l[m,r]' in warsaw,
      'Majorana_mass_reference_has_explicit_half':
          'm/2 Bar[CConj[field[ind1]]]' in mass,
    }
    if not all(anchors.values()):
        raise AuditError('Missing/tampered installed source identity: '+str([k for k,v in anchors.items() if not v]))
    return {'installed_source_status':data['status'], 'source_files_scanned':data['files_scanned'],
            'anchored_definitions':anchors,
            'source_root_in_evidence':data.get('source_root')}


def native_audit(root:Path,stem:str)->dict:
    if str(root) not in sys.path: sys.path.insert(0, str(root))
    from Numerical.experimental.T3BCMatcheteFiniteTerms import split_inputform
    path=(root/'output'/'general_t3_mass_eigen_matchete'/'eft1_scalar_Y_only_canonical_probe'/
           f'{stem}_final_y_only_canonical_match_raw_scalar_llhh.inputform.txt')
    if not path.is_file(): raise AuditError(f'Native source absent: {path}')
    src=path.read_text(encoding='utf-8-sig').strip()
    parts=split_inputform(src,stage='Y',sector='LLHH')
    counts={}
    for name,items in [('pole',parts.pole_terms),('log',parts.log_terms),('finite',parts.constant_terms)]:
        c=Counter()
        for idx,t in enumerate(items):
            if 'Proj[-1]' in t and 'Proj[1]' not in t:
                chir='PL'
                needed=['NCM[Transp[Field[l, Fermion,','DiracProduct[GammaCC, Proj[-1]]',
                        'CG[Bar[eps[SU2L]]','Field[H, Scalar']
                if ('Bar[Field[H, Scalar' in t or 'NCM[Bar[Field[l' in t):
                    raise AuditError(f'{stem} {name} {idx}: inconsistent direct Warsaw orientation')
            elif 'Proj[1]' in t and 'Proj[-1]' not in t:
                chir='PR'
                needed=['NCM[Bar[Field[l, Fermion,','DiracProduct[GammaCC, Proj[1]]',
                        'CG[eps[SU2L]','Bar[Field[H, Scalar']
                if ('CG[Bar[eps[SU2L]]' in t or 'NCM[Transp[Field[l' in t):
                    raise AuditError(f'{stem} {name} {idx}: inconsistent H.c. Warsaw orientation')
            else:
                raise AuditError(f'{stem} {name} {idx}: unknown projector chirality')
            if not all(x in t for x in needed):
                raise AuditError(f'{stem} {name} {idx}: wrong spinor/H/epsilon inventory')
            if t.count('GammaCC')!=1 or t.count('Field[l, Fermion')!=2:
                raise AuditError(f'{stem} {name} {idx}: unexpected lepton structure')
            if t.count('Field[H, Scalar')!=2:
                raise AuditError(f'{stem} {name} {idx}: unexpected Higgs structure')
            c[chir]+=1
        counts[name]=dict(c)
    if any(counts[name] != {'PL':exp,'PR':exp} for name,exp in [('pole',25),('log',25),('finite',16)]):
        raise AuditError(f'{stem}: wrong known native Y/LLHH term counts: {counts}')
    return {'source':str(path.relative_to(root)),'terms':counts,
            'Warsaw_chirality_and_H_orientations_verified':True}


def live_probe(path:Path) -> dict:
    data=json.loads(path.read_text(encoding='utf-8-sig'))
    if data.get('format')!='t3-bc-warsaw-direct-operator-identity-probe-v11':
        raise AuditError('Unexpected live Matchete probe schema')
    flag=data.get('direct_equivalence_status')
    if flag!='ExactPlusEquality' or data.get('negative_sign_equality'):
        raise AuditError('Installed Matchete failed direct positive-sign operator equivalence: '+str(flag))
    return {'direct_equivalence_status':flag,'probe_source':str(path)}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path('.'))
    ap.add_argument('--installed-source',type=Path,default=Path('output/audits/t3_bc_y_gammacc_convention_v1/installed_source.json'))
    ap.add_argument('--matchete-identity',type=Path)
    ap.add_argument('--output',type=Path,default=Path('output/audits/t3_bc_y_warsaw_direct_identity_v11/summary.json'))
    ap.add_argument('--self-test',action='store_true')
    ns=ap.parse_args()
    algebra=spinor_algebra()
    if ns.self_test:
        print('PASS: exact chiral GammaCC, Grassmann flavour, and native PL/PR H.c. sign')
        print('PASS: wrong sign rejected; normalized-Warsaw factor 2 follows from symmetric flavour')
        print('PASS: no full EFT matching or physical C5 asserted')
        return
    root=ns.root.resolve()
    evidence_path=ns.installed_source if ns.installed_source.is_absolute() else root/ns.installed_source
    ev=evidence_check(evidence_path)
    native={stem:native_audit(root,stem) for stem in ('B_m1','C_m1')}
    live=None
    if ns.matchete_identity is not None:
        identity=ns.matchete_identity if ns.matchete_identity.is_absolute() else root/ns.matchete_identity
        live=live_probe(identity)
    report={'format':'t3-bc-warsaw-spinor-orientation-v11',
            'status':'InstalledSourceAndNativeSpinorAlgebraVerified' if live is None else 'InstalledMatcheteDirectIdentityAndSpinorAlgebraVerified',
            'spinor_algebra':algebra, 'installed_source':ev,'native':native,
            'native_PL_represents_Warsaw_Q_in_tested_spinor_convention':True,
            'native_PR_represents_Warsaw_Q_Hc_in_tested_spinor_convention':True,
            'native_spinor_relative_sign':'+1',
            'Warsaw_coefficient_to_half_normalized_C5':'C5_ij = A_ij + A_ji, for +1 native/Warsaw mapping',
            'mass_convention':'m_nu = -v^2/2 * C5',
            'live_Matchete_operator_identity':live,
            'absolute_installed_Matchete_operator_factor_certified':live is not None,
            'full_EFT1_EFT2_counterterms_certified':False,
            'fermion_first_EOMEFT_vs_canonical_equivalence_certified':False,
            'physical_C5_certified':False,
            'physical_C5':None}
    out=ns.output if ns.output.is_absolute() else root/ns.output
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    for stem in ('B_m1','C_m1'):
        print(f'{stem}: native Y LLHH 50 pole + 50 log + 32 finite; PL/PR Warsaw orientations PASS')
    print(f'Installed source anchors: {len(ev["anchored_definitions"])} PASS; explicit Weyl/Grassmann sign +1')
    print('Live Matchete exact equality:', live['direct_equivalence_status'] if live else 'NOT RUN')
    print('Physical C5: NOT CERTIFIED; result:',out)

if __name__=='__main__':
    try: main()
    except (AuditError, ImportError, ValueError, AssertionError, OSError) as e:
        raise SystemExit(f'FAIL: {e}')
