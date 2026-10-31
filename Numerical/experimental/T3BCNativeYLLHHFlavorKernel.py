"""Strict source-level neutral-component reduction of native Y->LLHH constants.

Each Matchete LLHH InputForm *term* is treated as a commutative product
of opaque factors until its exact field, gauge epsilon, flavour Yukawa,
spurion and spinor contents have been audited. The SU(2) neutral slice is
H_2=nonzero, nu_L=l_1, epsilon_12=+1. This is a COMPONENT coefficient,
NOT the normalized Weinberg Wilson coefficient. No Matchete hbar convention,
Majorana combinatorial factor, field-redefinition finite threshold,
MS-bar counterterm, or numerical quartic map is assumed.

No execution or sympification of input expressions; unrecognized factors abort.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
import re

from Numerical.experimental.T3BCPoleLogConsistency import polynomial, terms, ConsistencyError


class KernelError(ValueError):
    pass

SU2 = re.compile(r'Index\[(d\$\$\d+), SU2L\[fund\]\]')
LEPTON = re.compile(r'Field\[l, Fermion, \{Index\[(d\$\$\d+), SU2L\[fund\]\], Index\[(d\$\$\d+), Flavor\]\}')
YUKAWA = re.compile(r'^Coupling\[(y[12]), \{Index\[(d\$\$\d+), Flavor\], Index\[(d\$\$\d+), NFlavor\]\}, 0\]$')
COUPLING = re.compile(r'^Coupling\[([A-Za-z_0-9]+), \{\}, 0\]$')
EXPECTED_Q = frozenset(('lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1'))


def _parse_coupling(atom: str) -> tuple[str, bool] | None:
    isbar=atom.startswith('Bar[') and atom.endswith(']')
    inner=atom[4:-1] if isbar else atom
    m=COUPLING.fullmatch(inner)
    return (m.group(1),isbar) if m else None


def _index_from_field(atom: str, field: str) -> str | None:
    if atom.startswith(f'Field[{field}, ') or atom.startswith(f'Bar[Field[{field}, '):
        indices=SU2.findall(atom)
        if len(indices)!=1:
            raise KernelError(f'Unexpected {field} SU(2) index count')
        return indices[0]
    return None


def _epsilon_value(i: int, j: int) -> int:
    if (i,j)==(1,2):return 1
    if (i,j)==(2,1):return -1
    if i==j:return 0
    raise KernelError('SU(2) value is not fundamental 1 or 2')


def _record_from_term(term: str, model: str):
    if model not in ('B_-1','C_-1'):
        raise KernelError('Unsupported model')
    tab={}
    try:
        poly=polynomial((term,),table=tab)
    except ConsistencyError as exc:
        raise KernelError('Unrecognized Matchete input monomial') from exc
    if len(poly)!=1:
        raise KernelError('Expected one nonzero native monomial')
    (atoms,coeff),=poly.items()
    reverse={v:k for k,v in tab.items()}
    assigned={}
    hfields=[]; spinors=[]; epsilons=[]
    y={};u={};quartics=[];hbar=0;mf=0;sqrt3=0
    for opaque,exp in atoms:
        atom=reverse.get(opaque,opaque)
        if atom=='hbar':
            hbar+=exp;continue
        if atom=='Sqrt[3]':
            sqrt3+=exp;continue
        if atom.startswith('NCM['):
            if exp!=1:raise KernelError('Nonlinear NCM')
            spinors.append(atom);continue
        if atom.startswith('CG['):
            if exp!=1:raise KernelError('Nonlinear CG')
            if not (atom.startswith('CG[eps[SU2L], {') or
                    atom.startswith('CG[Bar[eps[SU2L]], {')):
                raise KernelError('Unrecognized SU(2) Clebsch')
            indices=SU2.findall(atom)
            if len(indices)!=2:raise KernelError('Non-binary epsilon')
            epsilons.append(indices)
            continue
        idx=_index_from_field(atom,'H')
        if idx is not None:
            if exp!=1:raise KernelError('Nonlinear Higgs field')
            hfields.append((idx,atom.startswith('Bar[')))
            if idx in assigned and assigned[idx]!=2:
                raise KernelError('Contradictory Higgs and lepton index')
            assigned[idx]=2
            continue
        barred=atom.startswith('Bar[') and atom.endswith(']')
        inner=atom[4:-1] if barred else atom
        ym=YUKAWA.fullmatch(inner)
        if ym:
            if exp!=1:raise KernelError('Nonlinear Yukawa')
            name,flavour,heavy=ym.groups()
            if heavy!='d$$1':raise KernelError('Unexpected heavy flavour slot')
            if flavour not in ('d$$1','d$$2') or flavour in y:
                raise KernelError('Unexpected lepton flavour slot')
            y[flavour]=(name,barred)
            continue
        coup=_parse_coupling(atom)
        if coup:
            name,isbar=coup
            if name=='MF':
                if isbar:raise KernelError('Conjugated MF not supported')
                mf+=exp
            elif name.startswith('lambda'):
                if exp!=1:raise KernelError('Nonlinear quartic')
                quartics.append((name,isbar))
            elif name.startswith('T3MEU'):
                if exp<1:raise KernelError('Inverse spurion')
                u[(name,isbar)]=u.get((name,isbar),0)+exp
            else:
                raise KernelError('Unexpected coupling '+name)
            continue
        raise KernelError('Unknown monomial factor: '+atom[:100])

    if len(spinors)!=1 or len(hfields)!=2 or len(epsilons)!=2 or len(quartics)!=1:
        raise KernelError('Field / epsilon / quartic inventory mismatch')
    if len(y)!=2 or sum(u.values())!=4 or mf!=-1 or hbar!=1:
        raise KernelError('Yukawa / mixing / mass / hbar inventory mismatch')
    if sqrt3!=(-1 if model=='C_-1' else 0):
        raise KernelError('Unexpected group normalization')
    nc=spinors[0]
    if nc.count('GammaCC')!=1 or nc.count('Field[l, Fermion')!=2:
        raise KernelError('Unexpected Majorana bilinear')
    lm=LEPTON.findall(nc)
    if len(lm)!=2 or {fl for _,fl in lm}!={'d$$1','d$$2'}:
        raise KernelError('Unexpected Majorana flavour indices')
    for gauge,flavour in lm:
        if gauge in assigned and assigned[gauge]!=1:
            raise KernelError('Contradictory lepton / Higgs index')
        assigned[gauge]=1
    if set(y)!={'d$$1','d$$2'}:
        raise KernelError('Incomplete Yukawa indices')
    p=re.findall(r'Proj\[(-?\d+)\]',nc)
    if len(p)!=1 or p[0] not in ('1','-1'):
        raise KernelError('Unrecognized Majorana chirality')
    chirality='PR' if p[0]=='1' else 'PL'
    hbarred={isbar for _,isbar in hfields}
    if len(hbarred)!=1 or (chirality=='PR')!=next(iter(hbarred)):
        raise KernelError('Higgs / Majorana chirality mismatch')
    if any(isbar!=(chirality=='PL') for _,isbar in y.values()):
        raise KernelError('Yukawa conjugation mismatch')
    sign=1
    for a,b in epsilons:
        if a not in assigned or b not in assigned:
            raise KernelError('Unassigned SU(2) index')
        sign*=_epsilon_value(assigned[a],assigned[b])
    if sign==0:raise KernelError('Unexpected zero neutral epsilon contraction')
    coeff*=sign
    qname,qbar=quartics[0]
    if qname not in EXPECTED_Q:
        raise KernelError('Unexpected finite Higgs quartic')
    record={
        'chirality':chirality,
        'quartic':qname,'quartic_conjugated':qbar,
        'yukawa_i':y['d$$1'][0], 'yukawa_j':y['d$$2'][0],
        'yukawas_conjugated':y['d$$1'][1],
        'mixing':tuple(sorted((name,bar,exp) for (name,bar),exp in u.items())),
        'sqrt3_exponent':sqrt3,
    }
    return record,coeff


def reduce_finite_y_llhh(source: str,model: str) -> dict:
    """Return exact neutral component spurion monomials, not a Wilson C5.

    The two chirality sectors are never summed: doing so would double count
    the Hermitian conjugate of the Majorana operator.
    """
    try:
        original=terms(source)
    except ConsistencyError as exc:
        raise KernelError('Invalid finite InputForm sum') from exc
    if len(original)!=32:
        raise KernelError('Expected exactly 32 audited Y/LLHH finite terms')
    out=defaultdict(Fraction)
    for raw in original:
        key,value=_record_from_term(raw,model)
        canonical=(key['chirality'],key['quartic'],key['quartic_conjugated'],
                   key['yukawa_i'],key['yukawa_j'],key['yukawas_conjugated'],
                   key['mixing'],key['sqrt3_exponent'])
        out[canonical]+=value
    out={key:val for key,val in out.items() if val}
    if not out or {k[0] for k in out}!={'PR','PL'} or {k[1] for k in out}!=EXPECTED_Q:
        raise KernelError('Missing chirality or quartic group after reduction')
    entries=[]
    for (ch,qn,qc,yi,yj,yc,mix,sqrt3),val in sorted(out.items()):
        entries.append({'chirality':ch,'quartic':qn,'quartic_conjugated':qc,
                        'yukawa_i':yi,'yukawa_j':yj,
                        'yukawas_conjugated':yc,
                        'mixing':[{'name':n,'conjugated':b,'exponent':e} for n,b,e in mix],
                        'sqrt3_exponent':sqrt3,
                        'coefficient_numerator':val.numerator,
                        'coefficient_denominator':val.denominator})
    chirality_counts={ch:sum(1 for e in entries if e['chirality']==ch) for ch in ('PR','PL')}
    groups={ch:{q:sum(1 for e in entries if e['chirality']==ch and e['quartic']==q)
                for q in sorted(EXPECTED_Q)} for ch in ('PR','PL')}
    return {
      'model':model, 'source_additive_terms':len(original),
      'neutral_SU2_slice':'H_component_2__lepton_component_1__epsilon12_plus_1',
      'original_fermion_bilinear_and_chirality_preserved':True,
      'one_common_factor':'hbar/Coupling[MF,{},0]',
      'chirality_monomial_counts':chirality_counts,
      'chirality_quartic_group_counts':groups,
      'monomials':entries,
      'MSbar_counterterms_certified':False,
      'hbar_to_loop_factor_conversion_certified':False,
      'Majorana_Wilson_operator_normalization_certified':False,
      'native_Y_finite_C5_physical_certified':False,
    }


def evaluate_spurion_kernel(record:dict, *, U, quartics, y1, y2, MF, hbar):
    """Evaluate formal neutral-component matrices at explicitly supplied inputs.

    Returns independent PL/PR complex 3x3 coefficient matrices, with the
    common hbar/MF included. Not a neutrino-mass or Weinberg C5 evaluation.
    """
    import numpy as np
    if record.get('model') not in ('B_-1','C_-1') or len(record.get('monomials',[]))<1:
        raise KernelError('Unrecognized kernel report')
    U=np.asarray(U,dtype=complex)
    y={'y1':np.asarray(y1,dtype=complex),'y2':np.asarray(y2,dtype=complex)}
    if U.shape!=(2,2) or any(v.shape!=(3,) for v in y.values()) or not np.isfinite(U).all():
        raise KernelError('Invalid unitary mixing or 3-flavour Yukawa arrays')
    if np.linalg.norm(U.conj().T@U-np.eye(2))>1e-9:
        raise KernelError('Mixing U must be unitary')
    if not all(np.isfinite(v).all() for v in y.values()):
        raise KernelError('Invalid Yukawas')
    if not np.isfinite(MF) or not np.isfinite(hbar) or MF<=0:
        raise KernelError('Invalid MF/hbar')
    mats={'PR':np.zeros((3,3),complex),'PL':np.zeros((3,3),complex)}
    for mon in record['monomials']:
        name=mon['quartic']
        if name not in quartics:raise KernelError('Missing quartic '+name)
        lam=complex(quartics[name])
        if mon['quartic_conjugated']:lam=lam.conjugate()
        weight=(mon['coefficient_numerator']/mon['coefficient_denominator'])*\
            (3**(mon['sqrt3_exponent']/2))*lam*(hbar/MF)
        for fac in mon['mixing']:
            s=fac['name']
            if not re.fullmatch(r'T3MEU[12][12]',s):
                raise KernelError('Unknown mixing spurion')
            k=U[int(s[-2])-1,int(s[-1])-1]
            if fac['conjugated']:k=k.conjugate()
            weight*=k**fac['exponent']
        left=y[mon['yukawa_i']]; right=y[mon['yukawa_j']]
        if mon['yukawas_conjugated']:
            left=left.conj();right=right.conj()
        mats[mon['chirality']]+=weight*np.outer(left,right)
    return mats


def check_exact_chirality_conjugacy(record:dict) -> bool:
    """Exact matching of PL/PR monomials under formal complex conjugation.

    This checks the *scalar/flavour prefactor* only; it does not establish
    charge-conjugation matrix conventions for the bilinear or C5 normalization.
    """
    mon=record.get('monomials',[])
    def key(x):
        return (x['chirality'],x['quartic'],x['quartic_conjugated'],
                x['yukawa_i'],x['yukawa_j'],x['yukawas_conjugated'],
                tuple(sorted((v['name'],v['conjugated'],v['exponent'])
                             for v in x['mixing'])),x['sqrt3_exponent'])
    bykey={key(x):Fraction(x['coefficient_numerator'],x['coefficient_denominator'])
           for x in mon}
    if len(bykey)!=len(mon):
        raise KernelError('Repeated canonical monomial in report')
    for x in mon:
        other=('PL' if x['chirality']=='PR' else 'PR',
               x['quartic'],not x['quartic_conjugated'],
               x['yukawa_i'],x['yukawa_j'],not x['yukawas_conjugated'],
               tuple(sorted((v['name'],not v['conjugated'],v['exponent'])
                            for v in x['mixing'])),x['sqrt3_exponent'])
        if bykey.get(other)!=Fraction(x['coefficient_numerator'],x['coefficient_denominator']):
            raise KernelError('PL/PR scalar prefactors not exact conjugate pairs')
    return True
