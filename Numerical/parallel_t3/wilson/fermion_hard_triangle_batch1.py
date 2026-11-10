"""Batch 1: source-resolved all-heavy-fermion hard LLHH triangle boundary.

The established source and matching implementations are reused without changes:
  B: same-run UV(MF) trajectory, Matchete PL LLSS sources and scalar CGs;
  A: explicit Matchete -> canonical real scalar R;
  hard_matching.py: finite MSbar heavy-fermion/scalar triangle local kernel.

Output is one topology, not a complete one-loop UV->EFT1 matching. The
GammaCC->Weyl and quartic normalisations retain the baseline assumptions.
No scalar is integrated out at MF, and no finite scalar term is inserted twice.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Mapping

import numpy as np

from Numerical.parallel_t3.background.api import Background, load_formal_boundary
from Numerical.parallel_t3.background.contracts import tag
from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.coordinates import embed_c5
from Numerical.parallel_t3.wilson.source_audited_all16 import source_scalar_map, hash_file
from Numerical.lo_certification.completeness.hard_matching import (
    evaluate_hard_scalar_channel, weinberg_tensor)
from Numerical.parallel_t3.threshold.matching import master_weinberg_template

SCHEMA = 't3-source-resolved-fermion-hard-triangle-batch1-v1'
B_RUN = 'B_all16_qualitative_v1'
MF = 100_000.0
MS = 1_000.0
MIXED = frozenset(('B_-1', 'C_-1'))


def _load_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as file:
        return {name: np.asarray(file[name]).copy() for name in file.files}


def _json(path: Path) -> dict:
    result = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(result, dict):
        raise ValueError(f'Expected object {path}')
    return result


def _sha_array(x: np.ndarray) -> str:
    a = np.ascontiguousarray(x)
    return hashlib.sha256(a.tobytes()).hexdigest()


def complex_matrix(raw: Mapping) -> np.ndarray:
    a = np.asarray(raw['real'], float) + 1j * np.asarray(raw['imag'], float)
    if not np.isfinite(a).all():
        raise ValueError('Nonfinite source complex matrix')
    return a


def _uv_at_mf(bg: Background, model: str):
    """Reconstruct exactly B's saved UV endpoint, never refit or rerun RGBeta."""
    raw = bg.arrays
    scales = np.asarray(raw['mu_uv_gev'], float)
    if not np.isclose(scales[-1], MF, atol=1e-6, rtol=1e-11):
        raise ValueError('B UV trajectory does not end at the common MF')
    required = ('MF', 'y1', 'y2')
    for key in required:
        if 'uv__' + key not in raw:
            raise ValueError('B UV endpoint lacks '+key)
    yextra = {
        key: np.asarray(raw['uv__'+key][-1], complex)
        for key in ('y1Conj','y2Conj') if 'uv__'+key in raw
    }
    return SimpleNamespace(
        mu_gev=MF,
        MF=np.asarray(raw['uv__MF'][-1], complex),
        y1=np.asarray(raw['uv__y1'][-1], complex),
        y2=np.asarray(raw['uv__y2'][-1], complex),
        general_couplings=yextra,
        representation=SimpleNamespace(majorana_fermion=model in MIXED),
    )


def _source_rows(project: Path, model: str, source_path: Path,
                 mixing_U: np.ndarray | None):
    """Same source-operator parser and CG registry as B's formal boundary."""
    if model in MIXED:
        from Numerical.running.GeneralT3MixedExactFermionWilson import (
            read_native, normalise_mixed_term, insert_current_mixing)
        source, pl, cg, _ordered = read_native(model,project)
        if source_path.resolve() != (project/'output/general_observable/ordered_wilson_tensor_v4'/tag(model)/'source_flavour.json').resolve():
            raise ValueError('B/C native source path disagrees with B manifest')
        if mixing_U is None:
            raise ValueError('B/C native X/Y source requires same-run MF mixing U')
        rows=[]
        for saved in pl:
            row=normalise_mixed_term(saved['source_term_inputform'])
            row['source_index']=int(saved['source_index'])
            expanded, names=insert_current_mixing(row['operator_template_unit_flavor_inputform'],mixing_U)
            if set(names)!=set(saved['spurion_names']):
                raise ValueError('Mixed native U spurion inventory changed')
            if row['factors']!=saved['factors']:
                raise ValueError('Mixed Yukawa inventory changed')
            row['operator_template_unit_flavor_inputform']=expanded
            rows.append(row)
        return rows,cg['CGRegistry']
    from scripts.RunT3GeneralWilsonBoundary import _validate_seed
    doc=_json(source_path)
    return _validate_seed(doc,model),doc['CGRegistry']


def _per_fermion_boundary(model: str, rows: list[dict], cg_registry: list,
                          masses: np.ndarray, yukawas: Mapping[str,np.ndarray]
                          ) -> tuple[np.ndarray, np.ndarray]:
    """One Matchete tree term *per mass eigenstate*, using exact heavy flavour.

    Restrict the already heavy-mass-rotated Yukawa matrices to column k; each
    term retains its exact source CG, all source signs, spinor symmetrisation,
    and 1/M_k. The sum over k reconstructs the saved formal tree coefficient.
    """
    from scripts.RunT3GeneralWilsonBoundary import flavour_matrix, encode_matrix
    from Numerical.running.GeneralT3OrderedTensor import project_ordered_term
    from Numerical.running.GeneralT3DynamicFermionWilson import _dense_operator
    masses=np.asarray(masses,float)
    if masses.shape!=(3,) or np.any(masses<=0):
        raise ValueError('Three positive heavy masses required')
    all_components=[]
    index_reference=None
    for k in range(3):
        masked={}
        for name, matrix in yukawas.items():
            y=np.asarray(matrix,complex)
            if y.shape!=(3,3):
                raise ValueError('Unexpected heavy-mass-basis Yukawa geometry')
            single=np.zeros_like(y)
            single[:,k]=y[:,k]
            masked[name]=single
        operators=[]
        current=[]
        for saved in rows:
            row=dict(saved)
            row['flavour_matrix_axes_external_indices']=[f['external_index'] for f in row['factors']]
            row['ordered_flavour_matrix_gev_inv']=encode_matrix(flavour_matrix(row,masses,masked))
            operators.append(project_ordered_term(row,cg_registry,model))
            current.append(row)
        by,indices,_audit=_dense_operator(operators,current,model)
        if index_reference is None:
            index_reference=tuple(indices)
        elif tuple(indices)!=index_reference:
            raise ValueError('Heavy eigenstate resolved Matchete indices changed')
        all_components.append(by.sum(axis=0))
    return np.stack(all_components),np.asarray(index_reference,dtype=np.int64)


def _mass_matrix_ms(project: Path, input_run_id: str, model: str) -> np.ndarray:
    p=project/'output/parallel_t3/rough_inputs'/input_run_id/'c'/tag(model)
    return np.asarray(_load_npz(p/'background.npz')['mass_squared_gev2_ms'],float)


def _quartic_mf(project: Path, input_run_id: str, model: str) -> np.ndarray:
    p=project/'output/parallel_t3/rough_inputs'/input_run_id/'a'/tag(model)
    d=_load_npz(p/'background.npz')
    scales=np.asarray(d['mu_gev'],float)
    indices=np.flatnonzero(np.isclose(scales,MF,atol=1e-6,rtol=1e-11))
    if len(indices)!=1:
        raise ValueError('Missing unique EF1 quartic at common MF')
    return np.asarray(d['lambda4'][int(indices[0])],float)


def _basis_transform(raw: np.ndarray, R: np.ndarray, gamma_factor: complex) -> np.ndarray:
    raw=np.asarray(raw,complex)
    r=np.asarray(R,float)
    if r.shape!=(raw.shape[-1],raw.shape[-1]) or not np.allclose(r@r.T,np.eye(len(r)),atol=1e-10):
        raise ValueError('Invalid source-to-canonical hard matching basis')
    return complex(gamma_factor)*np.einsum('kfgijab,aA,bB->kfgijAB',raw,r,r,optimize=True)


def _relative(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-100))


def convert_hard_projector_normalisation(value: np.ndarray) -> np.ndarray:
    """Hard module uses amplitude tensor 2*T_master; convert to C full-sum."""
    t_hard=weinberg_tensor()
    t_master=master_weinberg_template()
    if not np.allclose(t_hard,2*t_master,rtol=1e-13,atol=1e-14):
        raise RuntimeError('Hard and downstream Weinberg coefficient templates changed')
    result=2*np.asarray(value,complex)
    if result.shape!=(3,3) or not np.isfinite(result).all():
        raise ValueError('Invalid hard matching 3x3 matrix')
    return (result+result.T)/2


def evaluate_source_resolved_fhard(project: Path, model: str, *, input_run_id: str,
                                   gamma_factor: float=1.0,
                                   hard_projection_tol: float=1e-7) -> tuple[dict,dict]:
    """Return arrays and provenance for *hard triangle topology only*.

    Source data are read ONLY from Workstream B (with source hash gates).
    MF quartic is the constructed RG-improved conditional tensor at MF.
    The corrected MS heavy scalar mass matrix is used for optional full-
    triangle diagnostics but does NOT enter the local hard_muF coefficient.
    """
    project=Path(project).resolve()
    from Numerical.matching.GeneralT3FlavourC5 import prepare_heavy_flavour
    folder=project/'output/parallel_t3/background'/B_RUN/tag(model)
    bg=Background.load(folder)
    formal, fmeta=load_formal_boundary(folder)
    original=np.asarray(formal['C_LLSS_formal'],complex)
    if fmeta.get('model')!=model or fmeta.get('active_BSM_fields_below_MF')!=['S1','S2']:
        raise ValueError('B source does not integrate all fermions while retaining both scalars')
    rel=str(fmeta['source_matchete_seed']).replace('\\','/')
    if Path(rel).is_absolute() or '..' in Path(rel).parts:
        raise ValueError('Unsafe Matchete source path')
    seed_path=project/rel
    if hash_file(seed_path)!=fmeta.get('source_matchete_seed_sha256'):
        raise ValueError('Saved Workstream B Matchete source SHA changed')
    R,evidence=source_scalar_map(project,model,seed_path)
    U=None
    if model in MIXED:
        U=complex_matrix(fmeta['mixed_native_diagnostics']['mixing_U_muS'])
    rows,cg=_source_rows(project,model,seed_path,U)
    needed={f['name'] for row in rows for f in row['factors']}
    uv=_uv_at_mf(bg,model)
    masses,yukawas,heavy_diag=prepare_heavy_flavour(uv,needed|{'y1','y2'})
    saved_masses=np.asarray(formal['heavy_mass_eigenvalues_gev'],float)
    if not np.allclose(masses,saved_masses,rtol=2e-8,atol=5e-5):
        raise ValueError('Current mass-basis decomposition does not match saved B heavy eigenmasses')
    per_native,indices=_per_fermion_boundary(model,rows,cg,masses,yukawas)
    per_total=per_native.sum(axis=0)
    match_residual=_relative(per_total,original)
    if match_residual>2e-7:
        raise ValueError(f'Resolved heavy eigenstates fail B formal boundary reproduction: {match_residual:g}')
    c=_basis_transform(per_native,R,gamma_factor)
    mass=_mass_matrix_ms(project,input_run_id,model)
    qfull=_quartic_mf(project,input_run_id,model)
    b=build_basis(model)
    if mass.shape!=(b.ns-4,b.ns-4) or qfull.shape!=(b.ns,)*4:
        raise ValueError('Invalid canonical quartic or corrected heavy scalar Hessian')
    quartic=qfull[:4,:4,4:,4:]
    hard=evaluate_hard_scalar_channel(c,masses,quartic,mass,MF,MS,
                                     projection_tolerance=hard_projection_tol,
                                     C_total_muF=c.sum(axis=0))
    c5=convert_hard_projector_normalisation(hard['matrices']['hard_muF'])
    if not np.isfinite(c5).all() or not hard['projection_passed']:
        raise ValueError('F hard triangle did not project to a pure Weinberg tensor: '+repr(hard['residuals']))
    # Conversion check against the full EFT1 embedding used downstream.
    initial=embed_c5(c5,b)
    if not np.isfinite(initial).all():
        raise ValueError('Nonfinite hard C5 canonical embedding')
    metadata={
        'schema':SCHEMA,'model':model,'status':'HardFScalarQuarticTriangleOnlyConditional',
        'mf_gev':MF,'ms_gev':MS,'gamma_cc_to_weyl_benchmark_factor':float(gamma_factor),
        'hard_amplitude_template_to_master_C5_factor':2.0,
        'tree_formal_B_reproduction_relative':match_residual,
        'source_operator_indices':indices.tolist(),
        'three_heavy_fermion_masses_gev':masses.tolist(),
        'Matchete_source_path':rel,'Matchete_source_sha256':hash_file(seed_path),
        'B_background_npz_sha256':hash_file(folder/'background.npz'),
        'B_fermion_boundary_npz_sha256':hash_file(folder/'fermion_boundary_formal.npz'),
        'source_map':evidence,'heavy_basis_diagnostics':heavy_diag,
        'hard_one_loop_projection_residuals':hard['residuals'],
        'hard_triangle_one_loop_not_complete_F_matching':True,
        'one_common_F_decoupling_scale':True,
        'scalars_active_below_MF':['S1','S2'],
        'both_scalars_out_at_one_common_MS':True,
        'scalar_pair_eigenvalues_gev2':np.linalg.eigvalsh(mass).tolist(),
        'mass_hierarchy_max_scalar_over_min_fermion':hard['mass_hierarchy_max_scalar_over_min_fermion'],
        'hard_input_C_LLSS_per_fermion_sha256':_sha_array(c),
        'MF_quartic_HHSS_sha256':_sha_array(quartic),
        'source_spinor_factor_certified':False,
        'complete_UV_operator_matching_certified':False,
        'physical_C5_certified':False,
        'matching_scope':'one fermion + two scalars with HHSS insertion; no other hard topologies or counterterms',
    }
    arrays={'C5_hard_F_triangle_MF':c5,
            'C5_hard_F_triangle_native_projector':hard['matrices']['hard_muF'],
            'C5_triangle_soft_diagnostic':convert_hard_projector_normalisation(hard['matrices']['scalar_threshold_muS']),
            'C5_triangle_eft1_scalar_log_diagnostic':convert_hard_projector_normalisation(hard['matrices']['eft1_scalar_log']),
            'C5_full_UV_triangle_diagnostic':convert_hard_projector_normalisation(hard['matrices']['full_triangle']),
            'heavy_masses_gev':np.asarray(masses,float),
            'source_indices':indices,
            'C_LLSS_per_fermion_source_native':per_native,
            'C_LLSS_per_fermion_canonical':c,
            'HHSS_quartic_MF':quartic,
            'mass_squared_heavy_MS':mass}
    return arrays,metadata


def install_hard_boundary(project: Path, model: str, *, input_run_id: str,
                          output: Path, gamma_factor: float=1.0) -> dict:
    """Replace only the isolated batch's zero C5 seed, preserving LLSS verbatim."""
    from Numerical.parallel_t3.wilson.coordinates import project_c5
    project=Path(project).resolve()
    out=Path(output)
    arrays,meta=evaluate_source_resolved_fhard(project,model,
        input_run_id=input_run_id,gamma_factor=gamma_factor)
    path=project/'output/parallel_t3/rough_inputs'/input_run_id/'a'/tag(model)
    npz=path/'boundary.npz';manifest=path/'boundary.json'
    d=_load_npz(npz)
    old=np.asarray(d['C5_hard'],complex)
    if not np.array_equal(old,np.zeros((3,3),complex)):
        raise ValueError('Hard triangle injection requires original explicit zero C5 benchmark')
    c5=arrays['C5_hard_F_triangle_MF']
    d['C5_hard']=c5
    # The provisional tensor has already been updated to the correct LLSS
    # frame, but is ignored by the boundary loader. Store it consistently.
    if 'C_full_provisional' in d:
        basis=build_basis(model)
        original=np.asarray(d['C_full_provisional'],complex)
        c0,_=project_c5(original,basis)
        if np.linalg.norm(c0)>1e-25:
            raise ValueError('C_full_provisional already contains an MF hard contribution')
        d['C_full_provisional']=original+embed_c5(c5,basis)
    npz_tmp=npz.with_suffix('.partial.npz')
    np.savez_compressed(npz_tmp,**d)
    npz_tmp.replace(npz)
    bmeta=_json(manifest)
    bmeta['npz_sha256']=hash_file(npz)
    bmeta['C5_hard_origin']='computed_conditional_F_scalar_quartic_triangle_only'
    bmeta['source_audited_hard_triangle_schema']=SCHEMA
    bmeta['complete_one_loop_fermion_hard_matching']=False
    manifest.write_text(json.dumps(bmeta,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
    out.mkdir(parents=True,exist_ok=True)
    out_npz=out/'fermion_hard_triangle.npz'
    tmp=out/'fermion_hard_triangle.partial.npz'
    np.savez_compressed(tmp,**arrays)
    tmp.replace(out_npz)
    meta['hard_npz_sha256']=hash_file(out_npz)
    meta['updated_boundary_npz_sha256']=hash_file(npz)
    meta_path=out/'fermion_hard_triangle.json'
    meta_path.write_text(json.dumps(meta,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
    return {'status':meta['status'],'MF_hard_C5_frobenius':float(np.linalg.norm(c5)),
            'C5_hard_MF_real':c5.real.tolist(),
            'C5_hard_MF_imag':c5.imag.tolist(),
            'fermion_hard_manifest':str(meta_path.relative_to(project)),
            'fermion_hard_npz':str(out_npz.relative_to(project)),
            'B_source_reproduction_relative':meta['tree_formal_B_reproduction_relative'],
            'hard_projection_residuals':meta['hard_one_loop_projection_residuals']}
