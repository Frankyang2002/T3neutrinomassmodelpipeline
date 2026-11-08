"""A(-4) EFT1 renormalisable tensors in the explicit *RGBeta* invariant frame.

This does NOT certify a Matchete-to-RGBeta scalar-CG/flavour-phase map.
The reference SU(2) triplet is the highest-weight basis (m=+1,0,-1).
Coupling signs use V4=+sum(lambda*RGBeta invariant). Any map from the
Matchete Lagrangian coefficients must be established algebraically first.

The UV F field and its y1/y2 Yukawas are deliberately absent after the
all-fermions-first threshold. Source input is the existing 17-key EFT1 RGBeta
report_betas together with a complete 3-generation numerical state.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path
from typing import Mapping
import json
import math

import numpy as np

MODEL = 'A_-4'
REVISION_CHECKED = '12128de602ac6bf109dbad50dbac8446ed629434'
BETA_CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'
BETA_NAMES = frozenset((
    'gY', 'g2', 'g3', 'yu', 'yd', 'ye', 'mS1Sq', 'mS2Sq',
    'lambdaH', 'lambdaS1', 'lambdaS2', 'lambdaH1', 'lambdaH2',
    'lambda12', 'lambdaT3', 'lambdaH2Adj', 'lambdaS2Adj',
))
SCALAR_COUPLINGS = (
    'lambdaH', 'lambdaS1', 'lambdaS2', 'lambdaS2Adj',
    'lambdaH1', 'lambdaH2', 'lambdaH2Adj', 'lambda12', 'lambdaT3',
)
SIGN_CONVENTION = (
    'V4 = +(lambdaH/2)(HdagH)^2 + (lambdaS1/2)|S1|^4 '
    '+ (lambdaS2/2)(S2dagS2)^2 '
    '+ (lambdaS2Adj/2)sum_A(S2dag T_A S2)^2 '
    '+ lambdaH1(HdagH)|S1|^2 '
    '+ lambdaH2(HdagH)(S2dagS2) '
    '+ lambdaH2Adj sum_A(Hdag t_A H)(S2dag T_A S2) '
    '+ lambda12|S1|^2(S2dagS2) '
    '+ [lambdaT3 (HH)_triplet S1 dot S2dag + h.c.]; '
    'H,S1,S2 complex components=(Re+i Im)/sqrt(2). '
    'Signs are RGBeta-potential reference definitions, NOT a certified '
    'Matchete coupling conversion.'
)


class EFT1BasisError(ValueError):
    """Explicit blocker for a missing physical or numerical source input."""


def scalar_labels() -> tuple[str, ...]:
    return tuple(f'{field}[{m},{ri}]' for field, dim in (
        ('H', 2), ('S1', 1), ('S2', 3)) for m in range(1, dim+1)
        for ri in ('R', 'I'))


def fermion_labels() -> tuple[str, ...]:
    result=[]
    for species in ('Q', 'uc', 'dc', 'L', 'ec'):
        for gen in range(1, 4):
            if species == 'Q':
                result += [f'Q[{gen},{col},{w}]' for col in range(1,4)
                           for w in range(1,3)]
            elif species in ('uc', 'dc'):
                result += [f'{species}[{gen},{col}]' for col in range(1,4)]
            elif species == 'L':
                result += [f'L[{gen},{w}]' for w in range(1,3)]
            else:
                result += [f'ec[{gen}]']
    return tuple(result)


F_LABELS = fermion_labels()
S_LABELS = scalar_labels()
FI = {s:i for i,s in enumerate(F_LABELS)}
SI = {s:i for i,s in enumerate(S_LABELS)}
assert len(F_LABELS)==45 and len(S_LABELS)==12


def index_map(external_labels: list[str] | tuple[str, ...],
              *, sector: str, require_complete: bool = False) -> np.ndarray:
    """External slot -> adapter 0-based slot; never infer Matchete order.

    For Wilson boundaries, use the returned arrays to place *each* external
    coordinate explicitly. This function deliberately does not pad C with zero.
    """
    if sector not in ('fermion','scalar'):
        raise EFT1BasisError('sector must be fermion or scalar')
    expected = FI if sector == 'fermion' else SI
    if len(set(external_labels)) != len(external_labels):
        raise EFT1BasisError('Duplicate external basis label')
    absent = set(external_labels)-set(expected)
    if absent:
        raise EFT1BasisError(f'Unknown {sector} labels: {sorted(absent)}')
    if require_complete and set(external_labels) != set(expected):
        raise EFT1BasisError('External basis is not the complete EFT1 basis')
    return np.asarray([expected[k] for k in external_labels],dtype=np.int64)


def su2_generators(dim: int) -> tuple[np.ndarray, ...]:
    spin=(dim-1)/2
    ms=[spin-i for i in range(dim)]
    plus=np.zeros((dim,dim),complex)
    for col,m in enumerate(ms):
        if m+1 in ms:
            plus[ms.index(m+1),col]=math.sqrt((spin-m)*(spin+m+1))
    minus=plus.conj().T
    return ((plus+minus)/2, (plus-minus)/(2j),
            np.diag(ms).astype(complex))


def su3_generators() -> tuple[np.ndarray, ...]:
    r=math.sqrt(3)
    g=[
        [[0,1,0],[1,0,0],[0,0,0]],
        [[0,-1j,0],[1j,0,0],[0,0,0]],
        [[1,0,0],[0,-1,0],[0,0,0]],
        [[0,0,1],[0,0,0],[1,0,0]],
        [[0,0,-1j],[0,0,0],[1j,0,0]],
        [[0,0,0],[0,0,1],[0,1,0]],
        [[0,0,0],[0,0,-1j],[0,1j,0]],
        [[1/r,0,0],[0,1/r,0],[0,0,-2/r]],
    ]
    return tuple(np.asarray(t,complex)/2 for t in g)


def real_scalar_generator(t: np.ndarray) -> np.ndarray:
    """Same interleaved convention as RGE/general/GaugeGenerators.py."""
    d=t.shape[0]; out=np.zeros((2*d,2*d),complex)
    for p in range(d):
        for q in range(d):
            a,b=t[p,q].real,t[p,q].imag
            out[2*p:2*p+2,2*q:2*q+2]=1j*np.array([[b,a],[-a,b]])
    return out


def _embed(out: np.ndarray, small: np.ndarray, ids: list[int]):
    out[np.ix_(ids,ids)]=small


def gauge_generators() -> dict[str, dict[str, tuple[np.ndarray,...]]]:
    """All 3+1+8 generators for 45 Weyl and 12 real components."""
    ff={'SU2':np.zeros((3,45,45),complex),
        'U1':np.zeros((1,45,45),complex),
        'SU3':np.zeros((8,45,45),complex)}
    s2=su2_generators(2)
    c3=su3_generators()
    q={'Q':1/6, 'uc':-2/3, 'dc':1/3, 'L':-1/2, 'ec':1}
    for field in ('Q','uc','dc','L','ec'):
        for gen in range(1,4):
            if field in ('Q','L'):
                if field=='Q':
                    for color in range(1,4):
                        ix=[FI[f'Q[{gen},{color},{w}]'] for w in range(1,3)]
                        for a in range(3):_embed(ff['SU2'][a],s2[a],ix)
                else:
                    ix=[FI[f'L[{gen},{w}]'] for w in range(1,3)]
                    for a in range(3):_embed(ff['SU2'][a],s2[a],ix)
            if field in ('Q','uc','dc'):
                nweak=2 if field=='Q' else 1
                for w in range(1,nweak+1):
                    ix=[FI[f'Q[{gen},{c},{w}]'] if field=='Q'
                        else FI[f'{field}[{gen},{c}]'] for c in range(1,4)]
                    for a in range(8):
                        _embed(ff['SU3'][a],c3[a] if field=='Q'
                               else -c3[a].conj(),ix)
    for i,s in enumerate(F_LABELS):
        ff['U1'][0,i,i]=q[s.split('[')[0]]
    ss={'SU2':np.zeros((3,12,12),complex),
        'U1':np.zeros((1,12,12),complex),
        'SU3':np.zeros((8,12,12),complex)}
    for name,d,y,offset in (('H',2,.5,0),('S1',1,-2.,4),('S2',3,-1.,6)):
        ids=list(range(offset,offset+2*d))
        for a,t in enumerate(su2_generators(d)):
            _embed(ss['SU2'][a],real_scalar_generator(t),ids)
        _embed(ss['U1'][0],real_scalar_generator(y*np.eye(d)),ids)
    return {k:{'fermion':tuple(ff[k]),'scalar':tuple(ss[k])}
            for k in ('SU3','SU2','U1')}


def _complex_array(value, shape, name):
    if isinstance(value,dict) and set(value)=={'real','imag'}:
        arr=np.asarray(value['real'],float)+1j*np.asarray(value['imag'],float)
    elif (isinstance(value,(list,tuple)) and len(value)==shape[0]
          and any(isinstance(item,dict) for row in value for item in row)):
        # Native RGBeta saved running files encode individual matrix entries
        # as {real,imag}, unlike SMNumericalState's combined 3x3 arrays.
        if any(len(row)!=shape[1] or any(not isinstance(item,dict)
                  or set(item)!={'real','imag'} for item in row) for row in value):
            raise EFT1BasisError(f'{name} invalid native RGBeta complex-entry matrix')
        arr=np.asarray([[complex(entry['real'],entry['imag'])
                         for entry in row] for row in value],complex)
    else: arr=np.asarray(value,complex)
    if arr.shape!=shape or not np.isfinite(arr).all():
        raise EFT1BasisError(f'{name} requires finite {shape} complex matrix')
    return arr


def _complex_scalar(value,name):
    if isinstance(value,dict) and set(value)=={'real','imag'}:
        z=complex(value['real'],value['imag'])
    else:z=complex(value)
    if not np.isfinite(z):raise EFT1BasisError(f'{name} nonfinite')
    return z


def _point(point: Mapping) -> dict:
    """Validate all inputs; even zero benchmark couplings MUST be supplied."""
    names=BETA_NAMES - {'mS1Sq','mS2Sq'}
    missing=names-set(point)
    if missing:raise EFT1BasisError(f'Missing EFT1 couplings: {sorted(missing)}')
    unknown=set(point)-BETA_NAMES
    if unknown:raise EFT1BasisError(f'Unexpected A(-4) couplings: {sorted(unknown)}')
    ans={}
    for name in BETA_NAMES:
        if name not in point:
            continue  # scalar mass parameters are optional for tensor evaluation
        v=point[name]
        if name in ('yu','yd','ye'):ans[name]=_complex_array(v,(3,3),name)
        elif name=='lambdaT3':ans[name]=_complex_scalar(v,name)
        else:
            z=_complex_scalar(v,name)
            if abs(z.imag)>1e-12:
                raise EFT1BasisError(f'{name}: required real')
            ans[name]=z.real
    return ans


def state_point(state) -> dict:
    """Native T3IntermediateScalarState -> current RGBeta coordinates."""
    state=state.validated()
    rep=state.representation
    if (rep.d_s1,rep.d_s2,rep.d_f,rep.alpha,rep.shared_scalar)!=(1,3,2,-4,False):
        raise EFT1BasisError('Expected ordinary A(-4) EFT1 state')
    if state.general_couplings:
        raise EFT1BasisError('A(-4) has no additional quartic invariant')
    names=('gY','g2','g3','lambdaH','yu','yd','ye')
    ans={k:getattr(state.sm,k) for k in names}
    for name in BETA_NAMES-set(ans):
        ans[name]=getattr(state,name)
    return _point(ans)


def yukawa_tensor(point: Mapping) -> np.ndarray:
    """L_Y=-(1/2)Y_ija psi_i psi_j phi_a+h.c.; Y_ija=Y_jia.

    Matrix ordering is (left-handed doublet generation, singlet generation).
    Quark colours are contracted in the 3 x anti-3 representation; the up
    doublets contract through epsilon_12=+1. Down and charged leptons use H*.
    No F, y1, y2, or scalar Yukawas survive in this EFT1 theory.
    """
    p=_point(point);Y=np.zeros((45,45,12),complex)
    def add(i,j,weak,hcomponent,v,conjugated):
        k=2*hcomponent
        Y[i,j,k]+=v/math.sqrt(2)
        Y[i,j,k+1]+=(-1j if conjugated else 1j)*v/math.sqrt(2)
        Y[j,i,k]=Y[i,j,k]
        Y[j,i,k+1]=Y[i,j,k+1]
    for g,h in product(range(1,4),repeat=2):
        for c in range(1,4):
            for w in range(1,3):
                for hcomp in range(1,3):
                    eps=(1 if (w,hcomp)==(1,2) else
                         -1 if (w,hcomp)==(2,1) else 0)
                    if eps:
                        add(FI[f'Q[{g},{c},{w}]'],FI[f'uc[{h},{c}]'],w,
                            hcomp-1,eps*p['yu'][g-1,h-1],False)
                    if w==hcomp:
                        add(FI[f'Q[{g},{c},{w}]'],FI[f'dc[{h},{c}]'],w,
                            hcomp-1,p['yd'][g-1,h-1],True)
        for w in range(1,3):
            add(FI[f'L[{g},{w}]'],FI[f'ec[{h}]'],w,w-1,
                p['ye'][g-1,h-1],True)
    return Y


def _real_hessian(t: np.ndarray) -> np.ndarray:
    """Hessian K of z^dag T z for z=(R+i I)/sqrt(2)."""
    d=t.shape[0]; K=np.zeros((2*d,2*d),float)
    for p in range(d):
        for q in range(d):
            a,b=t[p,q].real,t[p,q].imag
            K[2*p:2*p+2,2*q:2*q+2]=[[a,-b],[b,a]]
    return K


def _block_hessian(t: np.ndarray, off: int) -> np.ndarray:
    K=np.zeros((12,12),float);local=_real_hessian(t)
    K[off:off+local.shape[0],off:off+local.shape[1]]=local
    return K


def _quartic_product(A:np.ndarray,B:np.ndarray) -> np.ndarray:
    """Exact fourth derivative of (phi.A.phi)(phi.B.phi)/4."""
    result=np.zeros((12,12,12,12),float)
    for perm in ((0,1,2,3),(0,2,1,3),(0,2,3,1)):
        x=np.einsum('ij,kl->ijkl',A,B,optimize=True)
        x=x.transpose(perm)
        y=np.einsum('ij,kl->ijkl',B,A,optimize=True)
        y=y.transpose(perm)
        result+=x+y
    return result


def _mix_tensor() -> np.ndarray:
    """Reference highest-weight symmetric 2x2->3 CG, with unit CG norm."""
    cg=np.zeros((2,2,3),float)
    cg[0,0,0]=cg[1,1,2]=1
    cg[0,1,1]=cg[1,0,1]=1/math.sqrt(2)
    out=np.zeros((12,12,12,12),complex)
    for i,j,m in product(range(2),range(2),range(3)):
        if not cg[i,j,m]:continue
        for bits in product(range(2),repeat=4):
            ids=(2*i+bits[0],2*j+bits[1],4+bits[2],6+2*m+bits[3])
            v=cg[i,j,m]/4*(1j)**sum(bits[:3])*(-1j)**bits[3]
            for perm in permutations(range(4)):
                out[tuple(ids[k] for k in perm)]+=v
    return out


def quartic_columns() -> dict[str,np.ndarray]:
    """Reference RGBeta invariant coefficient columns, in real-scalar basis.

    No source-native Matchete CG normalization is inferred. The 4th derivative
    uses z=(R+iI)/sqrt(2), and includes all 24 assignments of each monomial.
    """
    idh=_block_hessian(np.eye(2),0)
    ids1=_block_hessian(np.eye(1),4)
    ids2=_block_hessian(np.eye(3),6)
    hsu=[_block_hessian(t,0) for t in su2_generators(2)]
    tsu=[_block_hessian(t,6) for t in su2_generators(3)]
    mix=_mix_tensor()
    return {
        'lambdaH':.5*_quartic_product(idh,idh),
        'lambdaS1':.5*_quartic_product(ids1,ids1),
        'lambdaS2':.5*_quartic_product(ids2,ids2),
        'lambdaS2Adj':sum((.5*_quartic_product(t,t) for t in tsu),
                          np.zeros((12,)*4)),
        'lambdaH1':_quartic_product(idh,ids1),
        'lambdaH2':_quartic_product(idh,ids2),
        'lambdaH2Adj':sum((_quartic_product(a,b) for a,b in zip(hsu,tsu)),
                          np.zeros((12,)*4)),
        'lambda12':_quartic_product(ids1,ids2),
        'lambdaT3':mix,
    }


def field_content_gauge_beta_coefficients() -> dict[str, str]:
    """Independent 1-loop b factors from the full active representations.

    For Weyl fermions and complex scalars:
      b = -11*C2(G)/3 + 2*sum_f(T_f)/3 + sum_s(T_s)/3.
    U(1) uses Y^2 with Q=T3+Y; no SU(5) rescaling.
    Includes all three generations and SU(3) colour dimensions.
    """
    F=Fraction
    fermions=(('Q',3,3,2,F(1,6)),('uc',3,3,1,F(-2,3)),
              ('dc',3,3,1,F(1,3)),('L',3,1,2,F(-1,2)),
              ('ec',3,1,1,F(1,1)))
    scalars=(('H',1,2,F(1,2)),('S1',1,1,F(-2,1)),
             ('S2',1,3,F(-1,1)))
    def isospin_index(d):return F(d*(d*d-1),12)
    by=F(2,3)*sum((n*dc*dw*y*y for _,n,dc,dw,y in fermions),F(0)) \
       +F(1,3)*sum((dc*dw*y*y for _,dc,dw,y in scalars),F(0))
    b2=F(-22,3)+F(2,3)*sum((n*dc*isospin_index(dw)
                  for _,n,dc,dw,_ in fermions),F(0)) \
       +F(1,3)*sum((dc*isospin_index(dw)
                  for _,dc,dw,_ in scalars),F(0))
    b3=F(-11)+F(2,3)*sum((n*dw*F(1,2) if dc==3 else F(0)
                  for _,n,dc,dw,_ in fermions),F(0))
    return {'gY':str(by),'g2':str(b2),'g3':str(b3)}


def one_loop_scalar_higgs_coefficients(columns: Mapping | None = None) -> dict[str, float]:
    """Independent phi^4 one-loop contraction projected onto lambdaH.

    Uses 16*pi^2 beta_lambda[abcd] |_(lambda^2) =
        lambda[abef]lambda[efcd] + lambda[acef]lambda[efbd]
      + lambda[adef]lambda[efbc].
    Valid because the lambdaH real-component column has value 3 at (0,0,0,0).
    This is a *diagnostic* for native RGBeta invariant normalization.
    """
    col = quartic_columns() if columns is None else columns
    result = {}
    for name, source in col.items():
        unit = (source + source.conjugate()).real if name == 'lambdaT3' else source.real
        beta = (np.einsum('abef,efcd->abcd', unit, unit, optimize=True)
                + np.einsum('acef,efbd->abcd', unit, unit, optimize=True)
                + np.einsum('adef,efbc->abcd', unit, unit, optimize=True))
        result[name] = float(beta[0, 0, 0, 0] / col['lambdaH'][0,0,0,0])
    return result


def source_normalization_audit(beta_payload: Mapping, *, columns: Mapping | None = None) -> dict:
    """Compare the ACTUAL native beta_lambdaH to independent scalar loops.

    Both the original individual-BetaTerm export and corrected global
    BetaTerm[Quartic,1] export are supported. No old quartic coefficient is
    hard-coded, fitted or substituted. This audit concerns a declared RGBeta
    reference frame, NOT Matchete-to-Python physical basis certification.
    """
    validate_beta_payload(beta_payload)
    gauge_betas = {'gY':'(55*gY^3)/6', 'g2':'(-5*g2^3)/2', 'g3':'-7*g3^3'}
    if field_content_gauge_beta_coefficients() != {'gY':'55/6','g2':'-5/2','g3':'-7'}:
        raise EFT1BasisError('Active-field gauge beta calculation inconsistent')
    if any(beta_payload['report_betas'][k].replace(' ','') != v
           for k,v in gauge_betas.items()):
        raise EFT1BasisError('EFT1 gauge beta functions disagree with active field content')
    from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
    expression = beta_payload['report_betas']['lambdaH']
    reference = one_loop_scalar_higgs_coefficients(columns)
    probe_names = ('lambdaH','lambdaH1','lambdaH2','lambdaH2Adj','lambdaT3')
    env = {k: 0. for k in BETA_NAMES}
    env.update({'yu':np.zeros((3,3),complex),
                'yd':np.zeros((3,3),complex),
                'ye':np.zeros((3,3),complex),
                'lambdaT3':0j, 'mS1Sq':1.0e6, 'mS2Sq':1.0e6})
    native = {}
    for name in probe_names:
        test = dict(env)
        test[name] = 1.+0j if name == 'lambdaT3' else 1.
        result = complex(evaluate_inputform_expression(expression,test))
        if not np.isfinite(result.real) or not np.isfinite(result.imag) or abs(result.imag)>1e-10:
            raise EFT1BasisError(f'Non-real/nonfinite native scalar beta coefficient: {name}')
        native[name] = float(result.real)
    mismatches = {k:{'reference':reference[k],'native':native[k]}
                  for k in probe_names if abs(reference[k]-native[k])>1e-10}
    return {
        'source':'actual native RGBeta beta_lambdaH evaluated on unit scalar couplings; independent real quartic loop',
        'native_beta_projection_kind': beta_payload.get('metadata',{}).get(
            'GeneralScalarQuarticBetaBasis','IndividualBetaTermLegacy'),
        'independent_gauge_beta_coefficients':field_content_gauge_beta_coefficients(),
        'reference_coefficients':reference,
        'native_coefficients':native,
        'unresolved_normalization_mismatches':mismatches,
        'reference_scalar_beta_lambdaH_agrees':not mismatches,
        'scalar_basis_physically_certified':False,
        'matchete_to_RGBeta_invariant_mapping_certified':False,
    }


def quartic_tensor(point: Mapping, *, columns: Mapping|None=None)->np.ndarray:
    p=_point(point)
    col=quartic_columns() if columns is None else columns
    # c*M + c* M_h.c. is used ONLY for the one complex topology quartic.
    result=np.zeros((12,)*4,complex)
    for name in SCALAR_COUPLINGS:
        if name=='lambdaT3':
            result+=p[name]*col[name]+p[name].conjugate()*col[name].conjugate()
        else:result+=p[name]*col[name]
    if not np.isfinite(result).all() or np.max(np.abs(result.imag))>1e-10:
        raise EFT1BasisError('Real-scalar lambda tensor must be finite and real')
    result=result.real
    if any(np.max(np.abs(result-result.transpose(perm)))>5e-11
           for perm in ((1,0,2,3),(0,2,1,3),(0,1,3,2))):
        raise EFT1BasisError('lambda not fully symmetric')
    return result


def source_contract(root:Path, payload:Mapping|None=None)->dict:
    """Check executable source and raw basis inventories, not completion flags."""
    root=Path(root)
    cat=(root/'Lagrangian/model/T3ModelCatalog.wl').read_text(encoding='utf8')
    t3=(root/'RGE/running/rgbeta/T3RGBetaModel.wl').read_text(encoding='utf8')
    other=(root/'Lagrangian/interactions/T3Topology.wl').read_text(encoding='utf8')
    v=(root/'Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl').read_text(encoding='utf8')
    reg_path=root/'RGE/general/generated/general_rgbeta_quartic_registry.json'
    reg=json.loads(reg_path.read_text(encoding='utf8'))['models']['A_-4']
    tests=(('model class','"A" -> {1, 3, 2}',cat),
           ('hypercharge','"Scalar1" -> alpha/2',cat),
           ('extra fermion trigger','T3GeneralExtraYukawaCount',other),
           ('quartic manifest','"1_3_-4" -> {}',v),
           ('EFT1 builder','T3RGBetaBuildEFT1[',t3),
           ('triplet self quartic','lambdaS2Adj,',t3),
           ('triplet portal','lambdaH2Adj,',t3),
           ('topology quartic','T3RGBetaMixInvariant[1, 3]',t3))
    for label,needle,content in tests:
        if needle not in content:
            raise EFT1BasisError(f'Source contract changed: {label}')
    if 'Dim[gen] = 3' not in t3 or 'AddGaugeGroup[g3, SU3c, SU[3]]' not in t3:
        raise EFT1BasisError('SM three-generation/colour source was changed')
    eft1_block=t3.split('T3RGBetaBuildEFT1[',1)[1].split('T3RGBetaEFT1OneLoopBetas[]',1)[0]
    for name in ('T3RGBetaAddSM[]', 'T3RGBetaAddBasicQuartics[',
                 'T3RGBetaAddMixingQuartic[', 'T3RGBetaAddRGClosedQuartics['):
        if name not in eft1_block:
            raise EFT1BasisError(f'Active EFT1 native builder missing {name}')
    for illegal in ('T3RGBetaAddHeavyFermion[','T3RGBetaAddT3Yukawas[',
                    'T3RGBetaAddGeneralFermionYukawas['):
        if illegal in eft1_block:
            raise EFT1BasisError('An integrated-out F Yukawa was retained in EFT1')
    runner=(root/'Numerical/running/IntermediateScalarRunner.py').read_text(encoding='utf8')
    evaluator=(root/'Numerical/running/IntermediateScalarRGBetaEvaluator.py').read_text(encoding='utf8')
    if 'def run_intermediate_scalar_segment(' not in runner or 'solve_ivp(' not in runner:
        raise EFT1BasisError('Native EFT1 numerical integration interface was changed')
    if 'FERMION_PARAMETERS' not in evaluator or 'validate_beta_keys' not in evaluator:
        raise EFT1BasisError('Native EFT1 beta/decoupling verification was changed')
    if (reg['dimensions'],reg['alpha'],reg['total_additional_quartics'],
        reg['quadratics_not_implemented']) != ([1,3,2],-4,0,0):
        raise EFT1BasisError('A(-4) raw registry differs from source inventory')
    if reg['new_rgbeta_symbols'] or reg['preexisting_rgbeta_symbols']:
        raise EFT1BasisError('A(-4) unexpected extra invariant registry')
    if payload is not None:validate_beta_payload(payload)
    sources=('Lagrangian/model/T3ModelCatalog.wl',
             'Lagrangian/interactions/T3Topology.wl',
             'Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl',
             'RGE/running/rgbeta/T3RGBetaModel.wl',
             'RGE/general/generated/general_rgbeta_quartic_registry.json',
             'Numerical/running/IntermediateScalarRunner.py',
             'Numerical/running/IntermediateScalarRGBetaEvaluator.py')
    native_path=root/'output/full16_symbolic/T3_A_alpha_m4/data/eft1_after_F_scalar_quartic_seed.json'
    native_audit=None
    if native_path.exists():
        native_audit=audit_native_seed(json.loads(native_path.read_text(encoding='utf-8-sig')))
        native_audit['native_archived_seed_sha256']=sha256(native_path.read_bytes()).hexdigest()
        # This archived Matchete seed may be restricted: record the inventory,
        # never treat absence there as absence from the active general theory.
    return {'model':MODEL,'source_sha256':{
        n:sha256((root/n).read_bytes()).hexdigest() for n in sources},
        'quartic_extra_invariants':0,'quadratic_extra_invariants':0,
        'additional_active_nonSM_yukawas':0,
        'native_archived_seed_audit':native_audit,
        'basis':'RGBeta reference, NOT Matchete CG',
        'matchete_rgbeta_normalisation_certified':False}


def validate_beta_payload(payload:Mapping)->None:
    meta=payload.get('metadata',{})
    required={'SharedScalar':False,'dS1':1,'dS2':3,'dF':2,'alpha':-4,
              'IntegratedField':'F', 'ActiveBSMFields':['S1','S2'],
              'ReportBetaConvention':BETA_CONVENTION}
    for k,v in required.items():
        if meta.get(k)!=v:
            raise EFT1BasisError(f'Wrong A(-4) EFT1 RGBeta metadata {k}')
    betas=payload.get('report_betas')
    if not isinstance(betas,dict) or set(betas)!=BETA_NAMES:
        raise EFT1BasisError(f'Incomplete EFT1 beta names: expected {sorted(BETA_NAMES)}')
    if any(not isinstance(s,str) or not s.strip() for s in betas.values()):
        raise EFT1BasisError('Empty beta function')
    forbidden=('y1','y2','y1Conj','y2Conj','MF')
    if set(forbidden)&set(betas):
        raise EFT1BasisError('F was not decoupled')
    if payload.get('status')!='Success':
        raise EFT1BasisError('Native RGBeta export failed')


def audit_native_seed(seed: Mapping)->dict:
    """Expose native term names WITHOUT claiming a basis map."""
    if seed.get('IntegratedFields')!=['F'] or seed.get('ActiveHeavyFields')!=['S1','S2']:
        raise EFT1BasisError('Native seed does not describe requested EFT1')
    terms=seed.get('ScalarQuarticTerms')
    if not isinstance(terms,list) or len(terms)!=seed.get('ScalarQuarticTermCount'):
        raise EFT1BasisError('Malformed native scalar seed')
    import re
    names=sorted({n for t in terms for n in
            re.findall(r'Coupling\[([^,\]]+)',t.get('TermInputForm',''))})
    return {'native_term_count':len(terms),'native_coupling_symbols':names,
            'native_source_tensor_map_verified':False,
            'note':'Native CG coefficients and exact basis change require algebra, not a positional alias.'}


@dataclass(frozen=True)
class EFT1TensorPoint:
    mu_gev: float
    fermion_basis: tuple[str,...]
    real_scalar_basis: tuple[str,...]
    Y: np.ndarray
    lambda4: np.ndarray
    gauge: dict
    couplings: dict

    def to_machine(self)->dict:
        return {'model':MODEL, 'mu_gev':self.mu_gev,
                'fermion_basis':list(self.fermion_basis),
                'real_scalar_basis':list(self.real_scalar_basis),
                'yukawa_operator':'L=-(1/2) Y_ija psi_i psi_j phi_a + h.c.',
                'quartic_operator':'V4=lambda_abcd phi_a phi_b phi_c phi_d / 4!',
                'source_invariant_convention':SIGN_CONVENTION,
                'yukawa_flavour_matrix_axis_order':'(Q or L generation, singlet generation)',
                'g1_convention':'gY: Q=T3+Y, not SU(5) rescaled',
                'fermion_symmetric':True,'quartic_fully_symmetric':True,
        'representation_order':{
            'Q':{'generations':3,'SU3':'3','SU2':2,'Y':'1/6'},
            'uc':{'generations':3,'SU3':'bar3','SU2':1,'Y':'-2/3'},
            'dc':{'generations':3,'SU3':'bar3','SU2':1,'Y':'1/3'},
            'L':{'generations':3,'SU3':'1','SU2':2,'Y':'-1/2'},
            'ec':{'generations':3,'SU3':'1','SU2':1,'Y':'1'},
            'H':{'SU3':'1','SU2':2,'Y':'1/2'},
            'S1':{'SU3':'1','SU2':1,'Y':'-2'},
            'S2':{'SU3':'1','SU2':3,'Y':'-1'}},
                'matchete_scalar_and_Wilson_basis_map_certified':False,
                'gauge_couplings':{k:float(self.couplings[k])
                    for k in ('gY','g2','g3')},
                'array_format':'NPZ; zero-based indices; complex NumPy arrays',
                'gauge_generator_array_shapes':{k:{s:list(np.shape(a[s]))
                    for s in ('fermion','scalar')} for k,a in self.gauge.items()}}


def tensor_point(mu_gev:float,point: Mapping,
                 *,columns:Mapping|None=None)->EFT1TensorPoint:
    p=_point(point)
    if not np.isfinite(mu_gev) or mu_gev<=0:
        raise EFT1BasisError('Invalid scale')
    return EFT1TensorPoint(float(mu_gev),F_LABELS,S_LABELS,
        yukawa_tensor(p),quartic_tensor(p,columns=columns),
        gauge_generators(),p)


def integrate_eft1(root:Path,initial_state,beta_file:Path,mu_s:float,
                   *,n_samples:int=5,rtol:float=1e-8,atol:float=1e-11):
    """Call existing 17-beta RGBeta EFT1 numerical solver; do NOT invent betas."""
    beta_file=Path(beta_file)
    payload=json.loads(beta_file.read_text(encoding='utf-8-sig'))
    validate_beta_payload(payload)
    source_contract(root,payload)
    state_point(initial_state)
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
    if n_samples<3:raise EFT1BasisError('Require at least 3 sampled scales')
    scales=np.geomspace(initial_state.mu_gev,mu_s,int(n_samples))
    out=run_intermediate_scalar_segment(initial_state,payload,mu_s,
        save_scales_gev=scales,rtol=rtol,atol=atol)
    if not out.solver_success or out.n_points<3:
        raise EFT1BasisError('EFT1 numerical running did not return >=3 scales')
    return out


def point_on_grid(running, mu_gev:float):
    """Piecewise LINEAR in t=ln(mu) on packed ODE state (no extrapolation)."""
    from Numerical.running.IntermediateScalarStateVector import unpack_intermediate_scalar_state
    grid=np.asarray(running.mu_gev,float)
    lo,hi=np.min(grid),np.max(grid)
    if not lo<=mu_gev<=hi:
        raise EFT1BasisError('Requested scale outside solved EFT1 interval')
    if not np.isfinite(mu_gev) or mu_gev<=0:
        raise EFT1BasisError('Invalid requested scale')
    t=np.log(grid)
    order=np.argsort(t)
    packed=np.asarray([np.interp(np.log(mu_gev),t[order],line[order])
                       for line in running.y],float)
    state=unpack_intermediate_scalar_state(packed,running.layout,mu_gev=mu_gev)
    return tensor_point(mu_gev,state_point(state))


def write_sampled_npz(running,folder:Path,*,source_root:Path,
                      source_beta:Path)->dict:
    """Full sampled tensors; generator matrices constant throughout EFT1."""
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    cols=quartic_columns()
    pts=[tensor_point(float(mu),state_point(running.state_at_index(i)),
                      columns=cols)
         for i,mu in enumerate(running.mu_gev)]
    matrices={'mu_gev':np.asarray(running.mu_gev),
              'Y':np.stack([p.Y for p in pts]),
              'lambda4':np.stack([p.lambda4 for p in pts])}
    for kind,sec in pts[0].gauge.items():
        for rep in ('fermion','scalar'):
            matrices[f'generator_{kind}_{rep}']=np.stack(sec[rep])
    np.savez_compressed(folder/'A_m4_eft1_tensors.npz',**matrices)
    info=pts[0].to_machine()
    info['sampled_scales_gev']=list(map(float,running.mu_gev))
    info['interpolation']='piecewise linear in log(mu) for packed native numerical state, then reconstruct tensors; no extrapolation'
    info['rge_payload_sha256']=sha256(Path(source_beta).read_bytes()).hexdigest()
    info['native_beta_normalization_audit']=source_normalization_audit(
        json.loads(Path(source_beta).read_text(encoding='utf-8-sig')),columns=cols)
    info['running_couplings']=[{
        name: ({'real':val.real.tolist(),'imag':val.imag.tolist()}
               if isinstance(val,np.ndarray) else
               {'real':val.real,'imag':val.imag}
               if isinstance(val,complex) else float(val))
        for name,val in p.couplings.items()} for p in pts]
    info['running_scalar_mass_squared_gev2']=[
        {'mS1Sq':p.couplings['mS1Sq'],'mS2Sq':p.couplings['mS2Sq']} for p in pts]
    info['both_heavy_scalars_active_until_mu_s']=True
    info['scalar_pole_mass_matching_certified']=False
    info['provenance']=source_contract(source_root)
    info['physical_full_Wilson_running_ready']=False
    info['why_not_fully_certified']='Matchete-to-RGBeta native invariant and Wilson-index mapping not certified'
    (folder/'A_m4_eft1_tensor_basis.json').write_text(json.dumps(info,indent=2)+'\n',encoding='utf8')
    return info
