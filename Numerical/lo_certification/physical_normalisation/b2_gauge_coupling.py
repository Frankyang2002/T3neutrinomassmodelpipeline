"""B2: read-only, conditional gauge-tensor cross-package audit for restricted T3.

Do not use this module as a Feynman-rule conversion. RGBeta's actual numerical
Clebsch components and the four-spinor-to-Weyl action remain unmeasured here.
All source/CG data are read at runtime, never silently substituted.
"""
from __future__ import annotations
import hashlib
import json
import re
import subprocess
from itertools import product
from pathlib import Path
from typing import Any
import sympy as sp

# Checked public GitHub main 0b210f18 on 2026-11-03. These are *Git blob*
# digests, NOT SHA256. They intentionally block silent source drift.
HEAD = '0b210f18f6a2cfb10fff03059e2ec4a73d2195a4'
SOURCE_BLOBS = {
    'RGE/running/rgbeta/T3RGBetaModel.wl': '7af7a4ddf99a7e38210543a21b3b1cda203eaf02',
    'Lagrangian/model/SU2Invariants.wl': 'dbe0cd490204dec363fcf4278f079b58c8ddcefc',
    'Lagrangian/interactions/T3Topology.wl': 'a45a099a8b766cbad0114a0559cb0bcb80895649',
    'Lagrangian/interactions/ScalarPotential.wl': '58f5bbef2688f9e87a1b8992bca11076ee254437',
    'Lagrangian/model/T3Fields.wl': '65467ede6ef06a7743e3ca2878debf777dc94507',
    'scripts/BuildT3DirectOneLoopInputs.py': '976a5c6303c637b7d1e45942a8dfa06f7e780ac7',
}
CLASSES = {
    'A': (1, 3, 2, (-4, -2, 0)),
    'B': (2, 2, 1, (-3, -1, 1)),
    'C': (2, 2, 3, (-3, -1, 1)),
    'D': (3, 1, 2, (-2, 0, 2)),
    'E': (3, 3, 2, (-4, -2, 0, 2)),
}
FLAGS = {'T3Y1CG': (True, True, False),
         'T3Y2CG': (True, False, True),
         'T3MixCG': (False, False, False, True)}
AUDITS = Path('output/general_observable/lo_direct_one_loop/source_audits')

class AuditBlocked(ValueError):
    pass


def require(condition: bool, explanation: str) -> None:
    if not condition:
        raise AuditBlocked(explanation)


def git_blob(data: bytes) -> str:
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()


def compact_wl(s: str) -> str:
    """Remove nested Mathematica comments, then whitespace: a guard, not parser."""
    chunks, depth, i = [], 0, 0
    while i < len(s):
        if s.startswith('(*', i):
            depth += 1; i += 2
        elif s.startswith('*)', i):
            require(depth > 0, 'stray WL comment terminator')
            depth -= 1; i += 2
        else:
            if depth == 0 and not s[i].isspace():
                chunks.append(s[i])
            i += 1
    require(depth == 0, 'unterminated WL comment')
    return ''.join(chunks)


def git_result(root: Path, *args: str) -> tuple[int | None, str]:
    """Read-only Git operation. Never mutate the index, checkout, or worktree.

    Return (None, '') when Git cannot run. This is deliberately not a
    substitute for a matched source hash.
    """
    try:
        result = subprocess.run(
            ['git', '-C', str(root), *args],
            capture_output=True, text=True, timeout=8, check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None, ''
    return result.returncode, result.stdout.strip()


def canonical_git_blob(root: Path, rel: str, expected: str, head: str | None) -> dict[str, Any]:
    """Inspect tracked source through Git's own clean/smudge conversion.

    On Windows, Git's HEAD blob can be LF while the working-tree source is
    CRLF. `git hash-object --path=...` applies the *current path-specific Git
    clean filter*, unlike a raw SHA1 over the worktree bytes. To trust that
    canonical value require the known commit, unchanged HEAD tree/index,
    and no working-tree or staged changes. A failed/unknown check blocks.

    Never use `git status --short` alone as evidence. Never replace native
    source evidence with a normalised string or heuristic LF conversion.
    """
    result: dict[str, Any] = {
        'git_canonical_checked': False, 'git_canonical_matches_reference': False,
        'git_path_is_clean_against_HEAD': False,
    }
    if head != HEAD:
        result['git_canonical_reason'] = 'local HEAD differs from audited reference'
        return result
    tree_rc, tree = git_result(root, 'ls-tree', 'HEAD', '--', rel)
    index_rc, index = git_result(root, 'ls-files', '--stage', '--', rel)
    canonical_rc, canonical = git_result(root, 'hash-object', f'--path={rel}', '--', rel)
    diff_rc, _ = git_result(root, 'diff', '--quiet', 'HEAD', '--', rel)
    result['git_canonical_checked'] = True
    tree_blob = tree.split()[2] if tree_rc == 0 and len(tree.split()) >= 4 else None
    index_line = index.splitlines()
    index_blob = None
    if index_rc == 0 and len(index_line) == 1:
        fields = index_line[0].split()
        if len(fields) >= 3 and fields[2] == '0':
            index_blob = fields[1]
    result.update({
        'git_HEAD_tree_blob_sha': tree_blob,
        'git_index_blob_sha': index_blob,
        'git_canonical_blob_sha': canonical if canonical_rc == 0 else None,
        'git_path_is_clean_against_HEAD': diff_rc == 0,
    })
    ok = (tree_blob == expected and index_blob == expected
          and canonical_rc == 0 and canonical == expected and diff_rc == 0)
    result['git_canonical_matches_reference'] = bool(ok)
    if not ok:
        result['git_canonical_reason'] = (
            'HEAD tree, index, canonical worktree content or Git diff differs; '
            'no source-provenance promotion'
        )
    return result


def source_inventory(root: Path) -> dict[str, Any]:
    # Direct byte equality is strongest. Canonical Git equality is also valid
    # if path-specific clean filters account for Windows checkout bytes.
    # `git -C some/subdirectory` normally walks upwards to an ancestor repo.
    # That ancestor is NOT evidence for the requested project root. Likewise,
    # GIT_DIR/GIT_WORK_TREE in the environment can redirect Git. Accept a HEAD
    # only if Git identifies this exact root as the worktree top level.
    top_rc, top = git_result(root, 'rev-parse', '--show-toplevel')
    top_path = Path(top).resolve() if top_rc == 0 and top else None
    root_path = root.resolve()
    top_matches_root = bool(top_path is not None and top_path == root_path)
    rev_rc, rev = git_result(root, 'rev-parse', 'HEAD') if top_matches_root else (None, '')
    head = rev if rev_rc == 0 and re.fullmatch(r'[0-9a-f]{40}', rev) else None
    report = {}
    for rel, expected in SOURCE_BLOBS.items():
        path = root / rel
        if not path.is_file():
            report[rel] = {'state': 'missing', 'expected_blob_sha': expected}
            continue
        payload = path.read_bytes()
        actual = git_blob(payload)
        record = {
            'state': 'reference_match' if actual == expected else 'drift_or_unreviewed',
            'raw_worktree_git_blob_sha': actual,
            'git_blob_sha': actual,  # keep previous output schema
            'expected_blob_sha': expected,
            'sha256': hashlib.sha256(payload).hexdigest(),
            'verification': 'exact_worktree_bytes' if actual == expected else 'not_verified',
        }
        if actual != expected:
            canonical = canonical_git_blob(root, rel, expected, head)
            record.update(canonical)
            if canonical['git_canonical_matches_reference']:
                record['state'] = 'reference_match'
                record['verification'] = 'path_specific_git_clean_filter'
        report[rel] = record
    if all(v['state'] == 'reference_match' for v in report.values()):
        files = {rel: compact_wl((root / rel).read_text(encoding='utf-8-sig'))
                 for rel in SOURCE_BLOBS if rel.endswith('.wl')}
        rgb = files['RGE/running/rgbeta/T3RGBetaModel.wl']
        topology = files['Lagrangian/interactions/T3Topology.wl']
        inv = files['Lagrangian/model/SU2Invariants.wl']
        scalar = files['Lagrangian/interactions/ScalarPotential.wl']
        for needle in ('AddYukawa[y1,{Bar@S1,l,F}', 'AddYukawa[y2,{S2,l,secondHeavy}',
                       'secondHeavy=If[yF===0,F,Fc]',
                       'AddFermionMass[MF,{F,Fc}', 'AddFermionMass[MF,{F,F}',
                       'dF===2,(eps[T3RGBetaRep[dF],#1,#2]&)',
                       'T3RGBetaMixInvariant[3,3]',
                       'AddQuartic[lambdaT3,{H,H,S1,Bar@S2}',
                       'AddQuartic[lambdaH1,{Bar@H,H,Bar@S1,S1}',
                       'AddQuartic[lambdaH2,{Bar@H,H,Bar@S2,S2}'):
            require(needle in rgb, 'RGBeta source contract changed: ' + needle)
        for needle in ('useCConj=which===1', 'scalarBar=which===2',
                       'y[p,r]NCM[Bar[l[i,p]],f]scalarCG[cg,labels]',
                       'PlusHc[lambdaT3[]H[i]H[j]s1s2CG[T3MixCG,labels]]'):
            require(needle in topology, 'Matchete source contract changed: ' + needle)
        require('InvariantTensors[SU[2],algebraReps' in inv and
                'DefineCG[cgName,cgReps,First[tensors]]' in inv,
                'Matchete CG source contract changed')
        require('1/2PlusHc[coupling[]fieldProductIf[cg===None,1,CG[cg,labels]]]' in scalar,
                'Matchete portal Hermitian convention changed')
    return {'reference_commit': HEAD, 'current_git_head': head,
            'current_git_head_measured': head is not None,
            'git_repository_root': str(top_path) if top_path is not None else None,
            'git_repository_root_matches_project': top_matches_root,
            'source_files': report,
            'source_reference_matches': all(v['state'] == 'reference_match' for v in report.values()),
            'source_contract_guard_is_full_parser': False}


def exact(src: Any) -> sp.Expr:
    text = str(src).replace('^', '**')
    require(len(text) < 160 and bool(re.fullmatch(r'[0-9IisqrtSQRt\-+*/(). _]*', text)),
            'non-numeric/untrusted native component: ' + text[:80])
    try:
        val = sp.sympify(text, locals={'I': sp.I, 'sqrt': sp.sqrt, 'Sqrt': sp.sqrt}, evaluate=True)
    except (ValueError, TypeError, SyntaxError) as exc:
        raise AuditBlocked('cannot parse exact native component') from exc
    require(val.is_number is True and val.is_finite is True,
            'native component is not exact finite numeric')
    return sp.simplify(val)


def source_rotations(generator: dict) -> dict[int, sp.Matrix]:
    require(generator.get('status') in ('generator_algebra_verified_CG_basis_unverified',
                                       'generator_algebra_verified', 'all_representations_verified'),
            'unsupported generator audit status')
    u = {1: sp.eye(1)}
    for row in generator['representations']:
        require(row['status'] == 'generator_algebra_verified', 'generator witness not verified')
        d = int(row['dimension'])
        require(d in (2, 3) and d not in u, 'invalid generator dimension')
        matrix = sp.Matrix([[exact(v) for v in items]
                            for items in row['U_columns_are_weight_eigenvectors']])
        require(matrix.shape == (d, d) and sp.simplify(matrix.H*matrix) == sp.eye(d),
                'source generator rotation not unitary')
        cartan = sp.Matrix([[exact(z) for z in line] for line in row['T3_in_original_matchete_basis']])
        target = sp.diag(*(sp.Rational(d-1,2)-i for i in range(d)))
        require(sp.simplify(matrix.H*cartan*matrix) == target,
                'wrong triplet/doublet weight-axis ordering or rotation')
        u[d] = matrix
    require(set(u) == {1, 2, 3}, 'generator rotations incomplete')
    return u


def native_tensor(vertex: dict, flags: tuple[bool, ...], u: dict[int, sp.Matrix]) -> sp.MutableDenseNDimArray:
    dims = tuple(int(d) for d in vertex['dimensions'])
    require(len(dims) == len(flags), 'native CG fields/axes mismatch')
    require(all(d in u for d in dims), 'unsupported native dimension')
    tensor = sp.MutableDenseNDimArray.zeros(*dims)
    for x in vertex['components']:
        idx = tuple(int(k)-1 for k in x['indices'])
        require(len(idx) == len(dims) and all(0 <= i < d for i,d in zip(idx,dims)),
                'invalid native CG index')
        value = exact(x['cg'])
        for target in product(*(range(d) for d in dims)):
            term = value
            for axis,(m,n) in enumerate(zip(idx,target)):
                term *= sp.conjugate(u[dims[axis]][m,n]) if flags[axis] else u[dims[axis]][m,n]
            tensor[target] += term
    return sp.MutableDenseNDimArray([sp.simplify(tensor[ix]) for ix in product(*(range(d) for d in dims))], dims)


def charge(d: int, y: sp.Rational, ix: int) -> sp.Rational:
    return sp.Rational(d-1,2)-ix+y


def mass_metrics(u: dict[int, sp.Matrix]) -> dict[str, Any]:
    eps = sp.Matrix([[0,1],[-1,0]])
    Pauli = [sp.Matrix([[0,1],[1,0]])/2,
             sp.Matrix([[0,-sp.I],[sp.I,0]])/2,
             sp.diag(sp.Rational(1,2),-sp.Rational(1,2))]
    require(all(t.T*eps+eps*t == sp.zeros(2) for t in Pauli),
            'epsilon SU2 Ward check failed')
    require(any(t.T*sp.eye(2)+sp.eye(2)*t != sp.zeros(2) for t in Pauli),
            'delta doublet negative control was accepted')
    # Cartesian adjoint chosen such that (T_a)_{bc}=-i epsilon_{abc}.
    tc = [sp.Matrix(3,3,lambda b,c:-sp.I*sp.LeviCivita(a,b,c)) for a in range(3)]
    g = sp.simplify(u[3].T*u[3]);tw = [sp.simplify(u[3].H*t*u[3]) for t in tc]
    require(all(t.T*g+g*t == sp.zeros(3) for t in tw), 'triplet weight metric Ward failed')
    # Weight-basis field rephasings change g by congruence and are valid.
    return {'epsilon_fundamental_convention': [[str(z) for z in r] for r in eps.tolist()],
            'triplet_metric_in_source_weight_basis': [[str(z) for z in r] for r in g.tolist()],
            'triplet_metric_is_a_convention_dependent_ward_witness': True,
            'native_RGBeta_epsilon_orientation_measured': False,
            'native_RGBeta_triplet_del_axes_measured': False}


def models() -> list[dict]:
    return [{'family': name, 'alpha': a, 'dS1': d1, 'dS2': d2, 'dF': df,
             'majorana': bool(a == -1 and df in (1,3)),
             'YF': str(sp.Rational(a+1,2)),
             'YS1': str(sp.Rational(a,2)), 'YS2': str(sp.Rational(a+2,2)),
             'label': f'{name}_{"m"+str(-a) if a<0 else "p"+str(a)}'}
            for name,(d1,d2,df,alphas) in CLASSES.items() for a in alphas]


def charge_channels(m: dict) -> list[dict]:
    d1,d2,df,a = (m[k] for k in ('dS1','dS2','dF','alpha'))
    qf = [charge(df,sp.Rational(a+1,2),i) for i in range(df)]
    q1 = [charge(d1,sp.Rational(a,2),i) for i in range(d1)]
    q2 = [charge(d2,sp.Rational(a+2,2),i) for i in range(d2)]
    return [{'Q':str(q), 'f':i+1,'s1':j+1,'s2':k+1}
            for i,q in enumerate(qf) for j,v in enumerate(q1) for k,w in enumerate(q2)
            if q == v == w]


def tensor_channels(m: dict, case: dict, u: dict[int,sp.Matrix]) -> list[dict]:
    """CG diagnostic: not a spinor-normalised loop coefficient."""
    d1,d2,df,a = (m[k] for k in ('dS1','dS2','dF','alpha'))
    t = {key:native_tensor(case['vertices'][key],FLAGS[key],u) for key in FLAGS}
    out = []
    for channel in charge_channels(m):
        f,s1,s2 = (channel[k]-1 for k in ('f','s1','s2'))
        y1,y2,mix = (t['T3Y1CG'][0,f,s1], t['T3Y2CG'][0,f,s2],
                     t['T3MixCG'][1,1,s1,s2])
        coeff = sp.simplify(sp.conjugate(y1)*sp.conjugate(y2)*mix)
        if all(sp.simplify(z) != 0 for z in (y1,y2,mix)):
            out.append({**channel, 'y1_CG':str(y1), 'y2_CG':str(y2),
                        'H0H0_mixing_CG':str(mix),
                        'conditional_CG_product': str(coeff),
                        'physical_amplitude_coefficient': None})
    return out


PORTAL_PRIOR = Path('output/general_observable/lo_certification_direct_uv/native_portal_dictionary_02/native_portal_operator_dictionary.json')


def prior_portal_dictionary(root: Path) -> dict[str, Any]:
    """Import ALL prior Hermitian nullspaces, only with current source-byte evidence."""
    p = root/PORTAL_PRIOR
    if not p.is_file():
        return {'status': 'prior_report_not_available', 'portal_dictionary_by_class': None,
                'native_source_hash_verified': False}
    try:
        obj = json.loads(p.read_text(encoding='utf-8-sig'))
        require(obj.get('format') == 't3-ws4-registered-native-portal-singlet-dictionary-v2'
                and obj.get('exact_operator_embedding_verified') is True
                and obj.get('physical_source_to_Weyl_model_certified') == 0
                and obj.get('native_RGBeta_quantum_action_prefactors_independently_certified') is False,
                'prior portal proof metadata/physical-certification flags invalid')
        listed = obj['source_evidence']['source_sha256']
        require(isinstance(listed,dict) and len(listed)>=6, 'prior portal source hashes incomplete')
        for rel,expected in listed.items():
            actual_path=root/Path(rel.replace('\\','/'))
            require(actual_path.is_file(), f'prior portal source missing: {rel}')
            require(hashlib.sha256(actual_path.read_bytes()).hexdigest() == expected,
                    f'prior portal proof stale against current file {rel}')
        byclass = obj['portal_dictionary_by_class']
        require(set(byclass)==set(CLASSES) and all(set(v)=={'S1','S2'} for v in byclass.values()),
                'prior portal 10-sector basis incomplete')
        nonunique=[]
        for cls,sub in byclass.items():
            for label,rec in sub.items():
                n=rec['native_invariant_count']
                rank=rec['real_Hermitian_operator_rank']
                ns=rec['nullspace_basis_real_coupling_shifts']
                require(rec['symbolic_identity_residual_zero'] is True and
                        len(ns)==n-rank and len(rec['coefficients_native_lambdaHiInvN_per_RGBeta_lambdaHi'])==n,
                        f'{cls}/{label} nullspace/identity evidence malformed')
                require(rec['unique_coupling_dictionary'] is (len(ns)==0),
                        f'{cls}/{label} falsely declared unique')
                require(all(len(vec)==n for vec in ns),
                        f'{cls}/{label} nullspace vector dimension invalid')
                if ns: nonunique.append(f'{cls}/{label}')
        require(set(nonunique)=={'A/S2','D/S1','E/S1','E/S2'},
                'expected four nonunique portal sectors altered')
        return {'status':'guarded_prior_exact_operator_identity_loaded',
                'prior_report_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
                'native_source_hash_verified':True,
                'portal_dictionary_by_class':byclass,
                'nonunique_sectors':nonunique,
                'physical_action_prefactor_certified':False}
    except (KeyError,TypeError,ValueError,OSError,json.JSONDecodeError,AuditBlocked) as exc:
        return {'status':'prior_report_blocked_or_stale',
                'reason':str(exc), 'native_source_hash_verified':False,
                'portal_dictionary_by_class':None}


def compute(root: Path) -> tuple[dict,dict]:
    src=source_inventory(root)
    base = {'format':'b2-gauge-coupling-conditional-v1',
            'git_reference':HEAD,'source_inventory':src,
            'physical_crosspackage_complete':False,'physical_source_certified':False,
            'physics_action_normalisation_certified':False,
            'missing_native_measurements':[
               'RGBeta eps, del, delS2 numerical component tensor in registered gauge axes and normalization',
               'RGBeta full potential/action coefficient, not internal AddYukawa/AddFermionMass prefactors',
               'Matchete↔RGBeta signed action-level mixing and portal relative factors',
               'B1 four-spinor→Weyl GammaCC and heavy-mass phase eta'],
            'source_ops':{'matchete_y1':'PlusHc[y1 Bar[L] CConj[F] S1 CG] (schematic)',
                          'matchete_y2':'PlusHc[y2 Bar[L] F Bar[S2] CG] (schematic)',
                          'matchete_mix':'PlusHc[lambdaT3 H H S1 Bar[S2] CG], no 1/2',
                          'rgbeta_y1':['Bar[S1]','l','F'],
                          'rgbeta_y2_dirac':['S2','l','Fc'],
                          'rgbeta_y2_majorana':['S2','l','F'],
                          'rgbeta_mass_dirac':['F','Fc'],
                          'rgbeta_mass_majorana':['F','F'],
                          'rgbeta_mix':['H','H','S1','Bar[S2]']},
            'spurion_invariant_conditional':'y1*y2*conjugate(MF)*conjugate(lambdaT3)',
            'portal_singlet_sectors':10,
            'portal_hermitian_nonunique_sectors':['A/S2','D/S1','E/S1','E/S2'],
            'portal_nullspaces_currently_loaded':False,
            'cross_package_yukawa_ratios':None,
            'cross_package_mass_ratio':None,'cross_package_mixing_ratio':None}
    prior = prior_portal_dictionary(root)
    base['prior_portal_evidence'] = prior
    base['portal_nullspaces_currently_loaded'] = prior['native_source_hash_verified']
    records = models()
    bycase = {}
    channels_verified = False
    if src['source_reference_matches']:
        try:
            gen=json.loads((root/AUDITS/'matchete_generator_analysis.json').read_text(encoding='utf-8'))
            comp=json.loads((root/AUDITS/'component_support_audit.json').read_text(encoding='utf-8'))
            cov=json.loads((root/AUDITS/'cg_covariance_audit.json').read_text(encoding='utf-8'))
            require(cov.get('status')=='cg_generator_covariance_verified' and cov.get('count')==16,
                    'covariance audit incomplete')
            require(len(comp['cases']) == 16, 'component support audit incomplete')
            u=source_rotations(gen)
            base['gauge_metrics']=mass_metrics(u)
            for case in comp['cases']:
                info=case['model'];key=tuple(int(info[k]) for k in ('dS1','dS2','dF','alpha'))
                require(key not in bycase, 'duplicate native CG case')
                bycase[key]=case
            require(set(bycase) == {(m['dS1'],m['dS2'],m['dF'],m['alpha']) for m in records},
                    'native audit model set differs')
            channels_verified=True
        except (OSError, KeyError, TypeError, json.JSONDecodeError, AuditBlocked) as exc:
            base['native_audit_blocker']=str(exc)
    else:
        base['native_audit_blocker']='Local source missing/drift relative to inspected main; no native-CG conclusions'
    rows=[]
    for m in records:
        key=tuple(m[k] for k in ('dS1','dS2','dF','alpha'))
        row={**m, 'charge_compatible_scalar_channels':charge_channels(m),
             'operator_charge_conservation_exact':True,
             'conditional_native_CG_channels':None,
             'independent_numeric_RGBeta_tensor_checked':False,
             'physical_source_certified':False}
        if channels_verified:
            row['conditional_native_CG_channels']=tensor_channels(m,bycase[key],u)
        rows.append(row)
    base['native_cg_snapshot_loaded'] = channels_verified
    base['native_cg_source_provenance_cryptographically_bound_to_current_execution'] = False
    base['all16_charge_skeletons_count']=len(rows)
    base['native_CG_diagnostic_models']=sum(r['conditional_native_CG_channels'] is not None for r in rows)
    return base, {'format':'b2-all16-conditional-CG-channels-v1',
                  'model_count':len(rows),'models':rows,
                  'charge_skeletons_are_analytic_only':True,
                  'conditional_CG_product_is_a_physical_loop_coefficient':False,
                  'native_physical_source_certified_models':0}


def portal_singlet_embedding(raw_tensors: list[sp.Matrix], d: int) -> dict:
    """Real-coupling linear solve after Matchete 1/2 PlusHc projection.

    Each raw tensor has rows (Hbar,H), columns (Sbar,S). Returns the ENTIRE
    affine solution, not a unique map when the Hermitian operators are dependent.
    """
    require(1 <= len(raw_tensors) <= 2, 'unexpected portal basis dimension')
    n = len(raw_tensors)
    require(all(t.shape == (4,d*d) for t in raw_tensors), 'portal tensor axes')
    def col(m):
        return sp.Matrix(list(m)).reshape(4*d*d,1)
    raw = sp.Matrix.hstack(*(col(t) for t in raw_tensors))
    require(raw.rank() == n, 'raw portal CGs linearly dependent')
    projected = [sp.Matrix(4,d*d,lambda ij,ab:
                 sp.simplify((t[ij,ab] +
                   sp.conjugate(t[(ij%2)*2+ij//2, (ab%d)*d+ab//d]))/2))
                 for t in raw_tensors]
    full=sp.Matrix.hstack(*(col(t) for t in projected))
    target=sp.Matrix(4,d*d,lambda ij,ab:
                     int(ij//2==ij%2 and ab//d==ab%d))
    rhs=col(target)
    r=full.applyfunc(sp.re).col_join(full.applyfunc(sp.im))
    b=rhs.applyfunc(sp.re).col_join(rhs.applyfunc(sp.im))
    require(r.rank()>0 and r.row_join(b).rank()==r.rank(),
            'RGBeta portal singlet outside Matchete Hermitian span')
    pivots=list(r.rref()[1]); vals,par=r[:,pivots].gauss_jordan_solve(b)
    require(len(par)==0, 'pivot coefficients not uniquely solved')
    sol=sp.zeros(n,1)
    for j,v in zip(pivots,vals): sol[j]=sp.simplify(v)
    ns=r.nullspace()
    require(full*sol==rhs and all(r*v==sp.zeros(r.rows,1) for v in ns),
            'Hermitian portal reconstruction failed')
    return {'unique':len(ns)==0,'rank':r.rank(),'native_CG_count':n,
            'canonical_pivot_coefficients':[str(z) for z in sol],
            'nullspace_basis':[[str(z) for z in v] for v in ns],
            'operator_equality_exact':True,
            'physical_RGBeta_action_normalisation_certified':False}


def scalar_A_mixing_offdiag(lam: sp.Expr, v:sp.Expr, cg:sp.Expr) -> tuple[sp.Expr,sp.Expr]:
    """Source-builder signed potential convention, NOT native RGBeta conversion."""
    left=sp.simplify(-lam*v**2*cg/2)
    return (left,sp.conjugate(left))


def spurion_phase_ward(majorana: bool) -> dict:
    """Coupling rephasing charges, given the explicit *conditional* B1 action."""
    y1=sp.Matrix([-1,-1,-1,0]) if majorana else sp.Matrix([0,-1,-1,0])
    y2=sp.Matrix([-1,0,0,1])
    mass=sp.Matrix([-2,0,0,0]) if majorana else sp.Matrix([-1,-1,0,0])
    mixing=sp.Matrix([0,0,-1,1])
    # Majorana uses a (no independent b); set b=0, y1 F*S1.
    total=y1+y2-mass-mixing
    if majorana:
        # remove unused 'b' by setting y1 = (-1,0,-1,0)
        y1=sp.Matrix([-1,0,-1,0])
        total=y1+y2-mass-mixing
    require(total==sp.zeros(4,1),'phase invariance failed')
    require(y1+y2-mass+mixing != sp.zeros(4,1),
            'wrong non-conjugated lambda unexpectedly invariant')
    return {'charges_y1':[int(z) for z in y1], 'charges_y2':[int(z) for z in y2],
            'charges_MF':[int(z) for z in mass],
            'charges_lambdaT3':[int(z) for z in mixing],
            'spurion_charges':[int(z) for z in total],
            'spurion':'y1*y2*conjugate(MF)*conjugate(lambdaT3)',
            'conditional_B1_action_only':True}
