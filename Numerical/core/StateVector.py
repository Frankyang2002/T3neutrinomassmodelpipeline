"""Canonical UV T3 real ODE vector packing, including optional general RGBeta couplings."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Literal
import numpy as np
from Numerical.core.State import SMNumericalState, SharedT3UVState, T3Representation, T3UVState
from Numerical.core.GeneralCouplingLayout import spec_for

BlockKind = Literal['real_scalar','complex_scalar','complex_matrix']
@dataclass(frozen=True)
class VectorBlock:
    name: str
    kind: BlockKind
    start: int
    stop: int
    shape: tuple[int,...] = ()
    @property
    def size(self) -> int: return self.stop-self.start

@dataclass(frozen=True)
class StateVectorLayout:
    representation: T3Representation
    reference_mu_gev: float
    shared_scalar: bool
    blocks: tuple[VectorBlock,...]
    @property
    def size(self) -> int: return self.blocks[-1].stop if self.blocks else 0
    def block(self, name: str) -> VectorBlock:
        matching = [b for b in self.blocks if b.name == name]
        if len(matching) != 1:
            raise KeyError(f'Expected exactly one vector block named {name!r}, found {len(matching)}.')
        return matching[0]
    @property
    def names(self) -> tuple[str,...]: return tuple(b.name for b in self.blocks)

class _VectorBuilder:
    def __init__(self) -> None:
        self.values: list[np.ndarray] = []
        self.blocks: list[VectorBlock] = []
        self.offset = 0
    def _append(self,name: str,kind: BlockKind,data: np.ndarray,shape: tuple[int,...]) -> None:
        if any(b.name == name for b in self.blocks):
            raise ValueError(f'Duplicate ODE state name: {name}')
        stop=self.offset+int(data.size)
        self.values.append(data)
        self.blocks.append(VectorBlock(name,kind,self.offset,stop,shape))
        self.offset=stop
    def add_real_scalar(self,name: str,value: float) -> None:
        self._append(name,'real_scalar',np.asarray([float(value)]),())
    def add_complex_scalar(self,name: str,value: complex) -> None:
        z=complex(value)
        self._append(name,'complex_scalar',np.asarray([z.real,z.imag]),())
    def add_complex_matrix(self,name: str,value: np.ndarray) -> None:
        m=np.asarray(value,dtype=complex)
        self._append(name,'complex_matrix',np.concatenate([m.real.reshape(-1),m.imag.reshape(-1)]),m.shape)
    def finish(self,*,representation: T3Representation,reference_mu_gev: float,shared_scalar: bool) -> tuple[np.ndarray,StateVectorLayout]:
        vector=np.concatenate(self.values) if self.values else np.empty(0,dtype=float)
        return vector,StateVectorLayout(representation,float(reference_mu_gev),shared_scalar,tuple(self.blocks))

def _pack_sm(builder: _VectorBuilder, sm: SMNumericalState) -> None:
    for name in ('gY','g2','g3'): builder.add_real_scalar(name,getattr(sm,name))
    for name in ('yu','yd','ye'): builder.add_complex_matrix(name,getattr(sm,name))
    builder.add_real_scalar('lambdaH',sm.lambdaH)

def pack_uv_state(state: T3UVState | SharedT3UVState) -> tuple[np.ndarray,StateVectorLayout]:
    validated=state.validated()
    b=_VectorBuilder()
    _pack_sm(b,validated.sm)
    if isinstance(validated,T3UVState):
        for name in ('y1','y2','MF'): b.add_complex_matrix(name,getattr(validated,name))
        for name in ('mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12'):
            b.add_real_scalar(name,getattr(validated,name))
        b.add_complex_scalar('lambdaT3',validated.lambdaT3)
        for name in ('lambdaH1Adj','lambdaH2Adj','lambdaS1Adj','lambdaS2Adj','lambda12Adj','lambda12Cross'):
            value=getattr(validated,name)
            if value is not None: b.add_real_scalar(name,value)
        for name in ('lambdaHHdagS2S2','lambdaHHdagS1barS1bar','lambdaS1bar2S2bar2',
                     'lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar','lambdaHHdagS1barS2barCross'):
            value=getattr(validated,name)
            if value is not None: b.add_complex_scalar(name,value)
        for name,value in sorted(validated.general_couplings.items()):
            kind=spec_for(name).kind
            if kind=='real': b.add_real_scalar(name,value)
            elif kind=='complex': b.add_complex_scalar(name,value)
            elif kind=='matrix': b.add_complex_matrix(name,value)
            else: raise ValueError(f'Unrecognised general coupling type: {name}')
        return b.finish(representation=validated.representation,reference_mu_gev=validated.mu_gev,shared_scalar=False)
    b.add_complex_matrix('h',validated.h)
    b.add_complex_matrix('MF',validated.MF)
    for name in ('mSSq','lambdaS','lambda3','lambda4','lambda5'):
        b.add_real_scalar(name,getattr(validated,name))
    return b.finish(representation=validated.representation,reference_mu_gev=validated.mu_gev,shared_scalar=True)

def _read_real_scalar(vector: np.ndarray,block: VectorBlock) -> float:
    d=vector[block.start:block.stop]
    if block.kind!='real_scalar' or d.size!=1: raise ValueError(f'Invalid real-scalar block for {block.name}.')
    return float(d[0])

def _read_complex_scalar(vector: np.ndarray,block: VectorBlock) -> complex:
    d=vector[block.start:block.stop]
    if block.kind!='complex_scalar' or d.size!=2: raise ValueError(f'Invalid complex-scalar block for {block.name}.')
    return complex(float(d[0]),float(d[1]))

def _read_complex_matrix(vector: np.ndarray,block: VectorBlock) -> np.ndarray:
    if block.kind!='complex_matrix' or len(block.shape)!=2:
        raise ValueError(f'Invalid complex-matrix block for {block.name}.')
    n=int(np.prod(block.shape));d=vector[block.start:block.stop]
    if d.size!=2*n: raise ValueError(f'Block {block.name} has {d.size} entries; expected {2*n}.')
    return d[:n].reshape(block.shape)+1j*d[n:].reshape(block.shape)

def _check_vector(vector: np.ndarray,layout: StateVectorLayout) -> np.ndarray:
    v=np.asarray(vector,dtype=float)
    if v.ndim!=1: raise ValueError('State vector must be one-dimensional.')
    if v.size!=layout.size: raise ValueError(f'State vector has length {v.size}; layout expects {layout.size}.')
    if not np.all(np.isfinite(v)): raise ValueError('State vector contains non-finite entries.')
    return v

def unpack_uv_state(vector: np.ndarray,layout: StateVectorLayout,*,mu_gev: float | None=None) -> T3UVState | SharedT3UVState:
    v=_check_vector(vector,layout)
    blocks={b.name:b for b in layout.blocks}
    def real(name: str) -> float: return _read_real_scalar(v,blocks[name])
    def cscalar(name: str) -> complex: return _read_complex_scalar(v,blocks[name])
    def cmatrix(name: str) -> np.ndarray: return _read_complex_matrix(v,blocks[name])
    sm=SMNumericalState(gY=real('gY'),g2=real('g2'),g3=real('g3'),lambdaH=real('lambdaH'),
                        yu=cmatrix('yu'),yd=cmatrix('yd'),ye=cmatrix('ye'))
    scale=layout.reference_mu_gev if mu_gev is None else float(mu_gev)
    if layout.shared_scalar:
        return SharedT3UVState(mu_gev=scale,representation=layout.representation,sm=sm,
            h=cmatrix('h'),MF=cmatrix('MF'),mSSq=real('mSSq'),lambdaS=real('lambdaS'),
            lambda3=real('lambda3'),lambda4=real('lambda4'),lambda5=real('lambda5')).validated()
    optional_real_names=('lambdaH1Adj','lambdaH2Adj','lambdaS1Adj','lambdaS2Adj','lambda12Adj','lambda12Cross')
    optional_complex_names=('lambdaHHdagS2S2','lambdaHHdagS1barS1bar','lambdaS1bar2S2bar2',
                            'lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar','lambdaHHdagS1barS2barCross')
    optional_real={name:real(name) if name in blocks else None for name in optional_real_names}
    optional_complex={name:cscalar(name) if name in blocks else None for name in optional_complex_names}
    base_names={'gY','g2','g3','yu','yd','ye','lambdaH','y1','y2','MF','mS1Sq','mS2Sq',
                'lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12','lambdaT3'}
    base_names.update(optional_real_names);base_names.update(optional_complex_names)
    extra={}
    for name,block in blocks.items():
        if name not in base_names:
            if block.kind=='real_scalar':extra[name]=real(name)
            elif block.kind=='complex_scalar':extra[name]=cscalar(name)
            elif block.kind=='complex_matrix':extra[name]=cmatrix(name)
            else:raise ValueError(f'Unknown ODE block kind for {name}')
    return T3UVState(mu_gev=scale,representation=layout.representation,sm=sm,
        y1=cmatrix('y1'),y2=cmatrix('y2'),MF=cmatrix('MF'),mS1Sq=real('mS1Sq'),mS2Sq=real('mS2Sq'),
        lambdaS1=real('lambdaS1'),lambdaS2=real('lambdaS2'),lambdaH1=real('lambdaH1'),
        lambdaH2=real('lambdaH2'),lambda12=real('lambda12'),lambdaT3=cscalar('lambdaT3'),
        **optional_real,**optional_complex,general_couplings=extra).validated()
