"""Exactly reconstruct selected RGBeta scalar invariants as named real quartic tensors.

No additional scalar interactions are guessed: every nonzero channel is taken
from Workstream B's SAME saved EFT1 state and the source-invariant unit tensor
already used to construct the completed all-16 benchmark. Native quartic
coordinates outside that restricted benchmark remain excluded explicitly.

At MF: decompose the existing heavy-fermion hard triangle by HHSS source
quartic channel. At MS: pass the named quartic tensors to Workstream C so
its existing simultaneous two-scalar bubble is resolved by channel.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np

from Numerical.parallel_t3.background.api import Background
from Numerical.parallel_t3.background.contracts import tag
from Numerical.parallel_t3.rough_inputs import builder
from Numerical.parallel_t3.wilson.source_audited_all16 import source_invariant_units
from Numerical.parallel_t3.wilson.fermion_hard_triangle_batch1 import (
    convert_hard_projector_normalisation)
from Numerical.lo_certification.completeness.hard_matching import evaluate_hard_scalar_channel

B_RUN='B_all16_qualitative_v1'


def named_quartics(root:Path,model:str,scale:str,reference:np.ndarray):
    """Return (names, individual tensors) for exactly the retained benchmark.

    scale 'MF'/'MS' selects saved B's native EFT1 endpoint at that scale.
    All omitted native couplings remain omitted, as in Batch 1; this is a
    decomposition of the *existing* conditional potential, not a completion.
    """
    if scale not in ('MF','MS'):
        raise ValueError('Expected MF or MS')
    folder=Path(root)/'output/parallel_t3/background'/B_RUN/tag(model)
    bg=Background.load(folder)
    pos=0 if scale=='MF' else -1
    units=source_invariant_units(model)
    expected=np.asarray(reference,float)
    names=[]; arrays=[]
    omitted={}
    for key in bg.manifest['beta_symbol_order_eft1']:
        if not key.startswith('lambda'):
            continue
        z=np.asarray(bg.arrays['eft1__'+key][pos]).item()
        if key not in units:
            if abs(z)>1e-15:
                omitted[key]={'real':float(np.real(z)),'imag':float(np.imag(z))}
            continue
        q,unused=builder.assemble_quartic(units,{key:z},reduced=True)
        if unused:
            raise ValueError('A supposedly mapped quartic was omitted: '+key)
        if np.max(abs(q))>0:
            names.append(key)
            arrays.append(q)
    if not arrays:
        raise ValueError('No saved invariant scalar quartic channels')
    split=np.stack(arrays)
    if split.shape[1:]!=expected.shape:
        raise ValueError('Named quartic tensor shape differs from reference')
    deviation=float(np.linalg.norm(np.sum(split,axis=0)-expected)/max(np.linalg.norm(expected),1e-100))
    if deviation>1e-10:
        raise ValueError(f'Named RGBeta quartic components do not reproduce benchmark tensor: {deviation:g}')
    if len(names)!=len(set(names)):
        raise ValueError('Duplicate quartic channel names')
    return tuple(names),split,{'scale':scale,'n_channels':len(names),
                              'reassembly_relative':deviation,
                              'omitted_nonzero_native_couplings':omitted,
                              'native_quartic_completion_claimed':False}


def split_existing_f_hard(root:Path,model:str,batch1:Path,names:tuple[str,...],
                          quartics_at_mf:np.ndarray)->tuple[dict,dict]:
    """Reproduce, then resolve Batch-1 hard coefficient by named HHSS quartics.

    Uses Batch-1's already-source-resolved per-heavy-fermion LLSS tensors and
    the same masses, same source scalar frame and same triangle kernel.
    """
    from Numerical.parallel_t3.wilson.source_audited_all16 import hash_file
    from Numerical.parallel_t3.wilson.basis import build_basis
    project=Path(root)
    hp=Path(batch1)/'hard_fermion'/tag(model)/'fermion_hard_triangle.npz'
    manifest=hp.with_name('fermion_hard_triangle.json')
    if not hp.is_file() or not manifest.is_file():
        raise FileNotFoundError('Batch-1 source-resolved triangle output required: '+str(hp))
    import json
    meta=json.loads(manifest.read_text(encoding='utf-8-sig'))
    if meta.get('model')!=model or meta.get('hard_npz_sha256')!=hash_file(hp):
        raise ValueError('Source-resolved heavy hard triangle provenance mismatch')
    with np.load(hp,allow_pickle=False) as file:
        per=np.asarray(file['C_LLSS_per_fermion_canonical'],complex)
        masses=np.asarray(file['heavy_masses_gev'],float)
        m2=np.asarray(file['mass_squared_heavy_MS'],float)
        qtotal=np.asarray(file['HHSS_quartic_MF'],float)
        saved=np.asarray(file['C5_hard_F_triangle_MF'],complex)
    allq=np.asarray(quartics_at_mf,float)
    if allq.ndim!=5 or len(allq)!=len(names):
        raise ValueError('Missing MF quartic channel decomposition')
    # All quartic channels are (k,ns,ns,ns,ns).
    qch=allq[:, :4, :4, 4:, 4:]
    if qch.shape!=(len(names),4,4,len(m2),len(m2)):
        raise ValueError('Incorrect MF HHSS quartic shapes')
    if not np.allclose(np.sum(qch,axis=0),qtotal,rtol=1e-10,atol=1e-12):
        raise ValueError('MF channel quartics fail source hard-triangle reassembly')
    if per.shape!=(3,3,3,2,2,len(m2),len(m2)):
        raise ValueError('Source-resolved per-heavy-fermion geometry changed')
    hard=evaluate_hard_scalar_channel(
        per,masses,qtotal,m2,100_000.,1_000.,
        C_total_muF=np.sum(per,axis=0),
        lambda_channels=qch,channel_names=names)
    c5=convert_hard_projector_normalisation(hard['matrices']['hard_muF'])
    relative=float(np.linalg.norm(c5-saved)/max(np.linalg.norm(saved),1e-100))
    if relative>1e-9:
        raise ValueError('Hard MF quartic source decomposition differs from Batch 1')
    pieces={key:convert_hard_projector_normalisation(
        hard['quartic_channels'][key]['matrices']['hard_muF']) for key in names}
    if np.linalg.norm(sum(pieces.values(),np.zeros_like(c5))-c5)>max(1e-22,1e-8*np.linalg.norm(c5)):
        raise ValueError('Hard triangle quartic channel C5 terms do not close')
    return pieces,{'n_channels':len(names),'relative_reproduction_batch1':relative,
                   'source_triangle_npz_sha256':hash_file(hp),
                   'all_three_fermions_included':True,
                   'source_normalisation_conditional':True}
