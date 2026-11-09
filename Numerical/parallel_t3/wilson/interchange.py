"""Strict B->A numeric contract; never map unlabelled native tensors by alias.

Format t3-parallel-a-input-v1:
 background.npz: mu_gev[N] (any order), Y[N,nf,nf,ns],
   lambda4[N,ns,ns,ns,ns], g2[N], gY[N], g3[N].
 boundary.npz: C_full[nf,nf,ns,ns] OR
   C_LLSS_formal_muF[3,3,2,2,ns-4,ns-4], and optional C5_hard[3,3].
 manifest JSONs: model, basis_id, fermion_labels, scalar_labels,
   MF_gev, MS_gev, npz_sha256, source_sha256{repo_relpath:SHA256},
   and boundary 'operator_basis'. Native boundary also needs an explicit
   orthogonal 'heavy_source_from_canonical' in boundary.npz and numeric
   gamma_cc_to_weyl_signed_factor in manifest. Only user-selected
   --assume-formal-identity mode supplies formal identity/factor=1.

The generator matrices are generated in canonical RGBeta reference basis and
can be checked against optional native generator_* arrays from Workstream B.
No Matchete-to-RGBeta coupling map is silently inferred.
"""
from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass
import hashlib, json
import numpy as np
from .basis import ModelBasis
from .coordinates import embed_c5,embed_formal_llss,symmetry_residual

FORMAT='t3-parallel-a-input-v1'
BASIS_ID='T3_RGBeta_reference_Weyl45_real_interleaved_v1'

class BlockedInput(ValueError):
    """Actionable physics/provenance blocker rather than silent tensor substitution."""


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with open(path,'rb') as fh:
        for chunk in iter(lambda:fh.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()


def _json(path:Path)->dict:
    data=json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(data,dict): raise BlockedInput('Blocked_ManifestNotObject')
    return data


def _source_check(manifest:dict, root:Path, npz_file:Path)->None:
    if manifest.get('npz_sha256')!=sha256(npz_file):
        raise BlockedInput('Blocked_NpzSHA256Mismatch: '+str(npz_file))
    refs=manifest.get('source_sha256')
    if not isinstance(refs,dict) or not refs:
        raise BlockedInput('Blocked_SourceHashMissing: B must supply source_sha256')
    for name,digest in refs.items():
        rel=Path(name)
        if rel.is_absolute() or '..' in rel.parts:
            raise BlockedInput('Blocked_UnsafeSourcePath')
        path=root/rel
        if not path.is_file():raise BlockedInput('Blocked_SourceMissing: '+str(rel))
        if sha256(path)!=digest:raise BlockedInput('Blocked_SourceHashMismatch: '+str(rel))


def _base_check(meta:dict,basis:ModelBasis)->None:
    if meta.get('model')!=basis.model:raise BlockedInput('Blocked_ModelMismatch')
    if meta.get('basis_id')!=BASIS_ID:raise BlockedInput('Blocked_BasisIdentifierMismatch')
    if tuple(meta.get('fermion_labels',()))!=basis.fermion_labels:
        raise BlockedInput('Blocked_WeylCoordinateMismatch: explicit unitary map required')
    if tuple(meta.get('scalar_labels',()))!=basis.scalar_labels:
        raise BlockedInput('Blocked_RealScalarCoordinateMismatch: explicit orthogonal map required')
    if meta.get('full_SM_three_generation_Weyl_basis') is not True:
        raise BlockedInput('Blocked_IncompleteWeylBasis')


def _portal_preflight(root:Path,basis:ModelBasis)->None:
    if basis.model!='A_-4':return
    path=root/'Lagrangian/interactions/ScalarPotential.wl'
    if not path.is_file():
        raise BlockedInput('Blocked_Aminus4ScalarPotentialSourceMissing')
    content=path.read_text(encoding='utf8')
    if not all(x in content for x in ('T3HiggsPortalHermitianPhase',
                                       'phase coupling[] fieldProduct',
                                       'phase === I', 'PlusHc')):
        raise BlockedInput('Blocked_Aminus4HermitianPortalRegression')


@dataclass(frozen=True)
class Background:
    basis:ModelBasis
    t:np.ndarray
    y:np.ndarray
    quartic:np.ndarray
    gauge_couplings:np.ndarray  # per row: g2,gY,g3
    generators:dict
    provenance:dict

    def at(self,t:float):
        if t<self.t[0]-1e-10 or t>self.t[-1]+1e-10:
            raise BlockedInput('Blocked_BackgroundExtrapolation')
        t=np.clip(t,self.t[0],self.t[-1])
        ix=min(int(np.searchsorted(self.t,t,side='right'))-1,len(self.t)-2)
        ix=max(ix,0)
        w=(t-self.t[ix])/(self.t[ix+1]-self.t[ix])
        y=(1-w)*self.y[ix]+w*self.y[ix+1]
        lam=(1-w)*self.quartic[ix]+w*self.quartic[ix+1]
        gs=(1-w)*self.gauge_couplings[ix]+w*self.gauge_couplings[ix+1]
        gauges=tuple((float(g),self.generators[k]) for g,k in zip(gs,('SU2','U1','SU3')))
        return y,lam,gauges


def load_background(root:Path,basis:ModelBasis,npz_path:Path,json_path:Path,
                    *, synthetic:bool=False)->Background:
    meta=_json(json_path)
    _base_check(meta,basis)
    if meta.get('MF_gev') is None or meta.get('MS_gev') is None:
        raise BlockedInput('Blocked_ThresholdMetadataMissing')
    mf,ms=float(meta['MF_gev']),float(meta['MS_gev'])
    if not 0<ms<mf:raise BlockedInput('Blocked_ThresholdOrderingMFgtMS')
    if not synthetic:_source_check(meta,root,npz_path)
    _portal_preflight(root,basis)
    with np.load(npz_path,allow_pickle=False) as d:
        req=('mu_gev','Y','lambda4','g2','gY','g3')
        if not all(k in d for k in req):
            raise BlockedInput('Blocked_BackgroundKeysMissing: '+repr(req))
        mu=np.array(d['mu_gev'],float)
        y=np.array(d['Y'],complex)
        lam=np.array(d['lambda4'],float)
        gauge=np.column_stack((d['g2'],d['gY'],d['g3'])).astype(float)
        gen=basis.generators()
        for k,n in [('SU2',3),('U1',1),('SU3',8)]:
            fk=f'generator_{k}_fermion';sk=f'generator_{k}_scalar'
            if (fk in d)!=(sk in d):
                raise BlockedInput('Blocked_HalfGeneratorPairProvided: '+k)
            if fk in d:
                fs,ss=np.asarray(d[fk]),np.asarray(d[sk]);expected=np.stack([t for t,_ in gen[k]])
                want=np.stack([t for _,t in gen[k]])
                if fs.shape!=(n,basis.nf,basis.nf) or ss.shape!=(n,basis.ns,basis.ns):
                    raise BlockedInput('Blocked_GeneratorShapeMismatch: '+k)
                if not np.allclose(fs,expected,rtol=1e-10,atol=1e-10) or not np.allclose(ss,want,rtol=1e-10,atol=1e-10):
                    raise BlockedInput('Blocked_GeneratorBasisMismatch: '+k)
    n=len(mu)
    if n<2 or y.shape!=(n,basis.nf,basis.nf,basis.ns) or lam.shape!=(n,)+(basis.ns,)*4 or gauge.shape!=(n,3):
        raise BlockedInput('Blocked_BackgroundShapes')
    if not all(np.isfinite(x).all() for x in (mu,y,lam,gauge)) or np.any(mu<=0):
        raise BlockedInput('Blocked_NonfiniteBackground')
    ix=np.argsort(mu);mu,y,lam,gauge=mu[ix],y[ix],lam[ix],gauge[ix]
    if np.any(np.diff(mu)<=0) or abs(np.log(mu[0]/ms))>1e-8 or abs(np.log(mu[-1]/mf))>1e-8:
        raise BlockedInput('Blocked_BackgroundDoesNotSpanExactThresholds')
    if not np.allclose(y,y.swapaxes(1,2),atol=1e-11,rtol=1e-10):
        raise BlockedInput('Blocked_YukawaExchange')
    if any(not np.allclose(lam,lam.swapaxes(1,axis),atol=1e-11,rtol=1e-10) for axis in (2,3,4)):
        raise BlockedInput('Blocked_QuarticPermutation')
    return Background(basis,np.log(mu),y,lam,gauge,gen,
                      {'background_npz_sha256':sha256(npz_path),'background_manifest':str(json_path),
                       'background_source_hashes':meta.get('source_sha256',{}),
                       'synthetic':synthetic})


def load_boundary(root:Path,basis:ModelBasis,npz_path:Path,json_path:Path,
                  *, assume_formal_identity:bool=False,synthetic:bool=False):
    meta=_json(json_path)
    _base_check(meta,basis)
    if meta.get('operator_basis')!='psi2phi2_no_factorial':
        raise BlockedInput('Blocked_WilsonOperatorNormalizationBasis')
    if not synthetic:_source_check(meta,root,npz_path)
    if basis.model=='A_-4' and any(x in str(meta.get('native_source_path','')) for x in (
        'fermion_first_all16_v1','eft1_stageB_saved_a_m4')):
        raise BlockedInput('Blocked_StaleAminus4FormalBoundary')
    if basis.model=='A_-4' and not synthetic:
        origin=str(meta.get('native_source_path','')).replace('\\','/')
        if (not meta.get('corrected_Aminus4_portal_evidence',False) or not origin):
            raise BlockedInput('Blocked_Aminus4CorrectedSourceProvenanceMissing')
        rel=Path(origin)
        if (rel.is_absolute() or '..' in rel.parts or
                not origin.startswith('output/general_observable/fermion_first_hermitian_portal_v1/A_m4')):
            raise BlockedInput('Blocked_Aminus4CorrectedSourcePathInvalid')
        if not any(name.replace('\\','/').startswith(origin.rstrip('/')+'/')
                   for name in meta.get('source_sha256',{})):
            raise BlockedInput('Blocked_Aminus4CorrectedSourceHashMissing')
    with np.load(npz_path,allow_pickle=False) as d:
        if 'C_full' in d:
            C=np.asarray(d['C_full'],np.complex128).copy()
            if meta.get('C_full_is_canonical') is not True:
                raise BlockedInput('Blocked_UncertifiedCanonicalFullTensor')
            origin='canonical_full_psi2phi2'
        elif 'C_LLSS_formal_muF' in d:
            raw=np.asarray(d['C_LLSS_formal_muF'],np.complex128)
            if 'heavy_source_from_canonical' in d:
                R=np.asarray(d['heavy_source_from_canonical'],float)
                basis_status='EXPLICIT_ORTHOGONAL_INPUT'
            elif assume_formal_identity:
                R=np.eye(basis.ns-4)
                basis_status='IDENTITY_ASSUMED_FOR_QUALITATIVE_TEST_ONLY'
            else:
                raise BlockedInput('Blocked_MatcheteRGBetaScalarMapMissing')
            fac=meta.get('gamma_cc_to_weyl_signed_factor')
            if fac is None:
                if assume_formal_identity: fac=1.
                else:raise BlockedInput('Blocked_GammaCCWeylSignedFactorMissing')
            C=embed_formal_llss(raw,basis,R,complex(fac))
            origin='formal_LLSS_'+basis_status
        else:
            raise BlockedInput('Blocked_MissingWilsonSourceCfullOrFormalLLSS')
        hard=np.asarray(d['C5_hard'],np.complex128) if 'C5_hard' in d else np.zeros((3,3),complex)
    if C.shape!=(basis.nf,basis.nf,basis.ns,basis.ns) or not np.isfinite(C).all():
        raise BlockedInput('Blocked_FullCGeometryOrFinite')
    if symmetry_residual(C)>max(1e-24,1e-10*np.max(np.abs(C))):
        raise BlockedInput('Blocked_BoundaryPairSymmetry')
    C=C+embed_c5(hard,basis)
    return C, {'origin':origin,'C5_hard_supplied':bool(np.any(hard)),
               'GammaCC_Weyl_factor_certified':False,
               'Matchete_RGBeta_map_certified':False,
               'boundary_npz_sha256':sha256(npz_path),
               'boundary_source_hashes':meta.get('source_sha256',{}),
               'boundary_manifest':str(json_path),'synthetic':synthetic}
