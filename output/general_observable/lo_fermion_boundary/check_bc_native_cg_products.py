"""B/C native Yukawa CG contraction and mixed-source tensor-basis audit.

Project path: output/general_observable/lo_fermion_boundary/check_bc_native_cg_products.py

Input is live Matchete PhysicalSU2InvariantBasis output, plus original Matchete
matched CG registry, pinned source records, and 16-model formal release.
No fit replaces mismatched tensors. Epsilon/bridge orientations are explicit.
The C_-1 tensor decomposition is NOT an independent native Feynman amplitude.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np

import bc_majorana_exchange_audit as bc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PROBES = HERE / 'normalization_probes'
EXPECTED = {'B_-1':1,'C_-1':3}


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array(record, name, dims):
    if not isinstance(record, dict) or record.get('name') != name or record.get('dimensions') != list(dims):
        raise ValueError(f'{name}: name/dimensions mismatch')
    if not record.get('representations_inputform'):
        raise ValueError(f'{name}: missing Matchete representation orientation')
    real, imag = record.get('real_flat'), record.get('imag_flat')
    if not isinstance(real,list) or not isinstance(imag,list) or len(real) != math.prod(dims) or len(imag) != math.prod(dims):
        raise ValueError(f'{name}: malformed components')
    a = np.array(real,dtype=float)+1j*np.array(imag,dtype=float)
    if not np.all(np.isfinite(a)) or np.max(abs(a))==0:
        raise ValueError(f'{name}: nonfinite/zero tensors')
    return a.reshape(dims)


def cg_registry(root, model):
    import re
    folder=root/'output/general_observable/ordered_wilson_tensor_v4'/model.replace('-','m').replace('+','p')
    path=folder/'mixed_native_tree_llss_cg_registry.json'
    data=read(path)
    if data.get('format')!='t3-mixed-matchete-native-cg-registry-v1' or data.get('model')!=model:
        raise ValueError(f'{model}: incompatible native matched CG registry')
    # Parse exact SparseArray InputForm as *data*. The numeric component lists
    # are generated from a fixed sparse-array grammar, never eval/exec.
    # Entries use Mathematica 1-based coordinates. The structure is sourced
    # by the existing Matchete matching run, not invented tensor components.
    out={}
    for entry in data['CGRegistry']:
        form=entry['TensorInputForm']
        m=re.fullmatch(r'SparseArray\[Automatic, \{([0-9, ]+)\}, 0, \{1, \{\{[0-9, ]+\}, \{(.+)\}\}, \{(.+)\}\}\]', form)
        if not m:
            raise ValueError(f'{model}: unsupported Matchete sparse-array serialization for {entry["Name"]}')
        dims=tuple(int(s.strip()) for s in m.group(1).split(','))
        coord_matches=re.findall(r'\{([0-9]+(?:, [0-9]+)*)\}',m.group(2))
        if not coord_matches:
            raise ValueError('Empty sparse-array coordinate inventory')
        coords=[tuple(int(k.strip()) for k in c.split(',')) for c in coord_matches]
        values=m.group(3).split(', ')
        if len(coords)!=len(values):
            raise ValueError(f'{entry["Name"]}: sparse components count mismatch')
        arr=np.zeros(dims,dtype=complex)
        for ix, val in zip(coords, values):
            if len(ix)!=len(dims) or any(not 1<=k<=d for k,d in zip(ix,dims)):
                raise ValueError('Sparse coordinate outside dimensions')
            v=val.strip()
            # The certified registries contain only 1, -1, 1/2.
            # Reject unfamiliar symbolic values rather than approximating.
            if v not in ('1','-1','1/2','-1/2'):
                raise ValueError(f'Unsupported native CG entry {v!r}')
            from fractions import Fraction
            arr[tuple(k-1 for k in ix)]=float(Fraction(v))
        out[entry['Name']]=arr
    return out,sha(path)


def relative_error(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-20))


def project(t, basis):
    a=np.column_stack([x.reshape(-1) for x in basis])
    coeff,_,rank,_=np.linalg.lstsq(a,t.reshape(-1),rcond=None)
    if rank != len(basis):
        raise ValueError('Dependent native invariant tensor basis')
    residue=np.linalg.norm(a@coeff-t.reshape(-1))/max(np.linalg.norm(t),1e-20)
    return {'real':[float(x.real) for x in coeff],
            'imag':[float(x.imag) for x in coeff],
            'relative_residual':float(residue), 'basis_rank':int(rank)}


def contraction(y1,y2):
    # The triplet's self-conjugate internal representation carries a delta
    # contraction in the native real basis. Scalars' doublet indices stay open.
    if y1.shape!=y2.shape or y1.ndim!=3 or y1.shape[0]!=2 or y1.shape[-1]!=2:
        raise ValueError('Unexpected Yukawa SU2 index structure')
    return np.einsum('ima,jmb->ijab',np.conjugate(y1),np.conjugate(y2))


def analyze(model, rec, registry):
    d=EXPECTED[model]
    if rec.get('fermion_su2_dimension')!=d:
        raise ValueError(f'{model}: unexpected heavy representation')
    dims=(2,d,2) if d==3 else (2,2)
    y1=array(rec['native_yukawa1_cg'],'T3MENativeY1CG',dims)
    y2=array(rec['native_yukawa2_cg'],'T3MENativeY2CG',dims)
    expected_alpha='-1' if model=='B_-1' else '1'
    if rec.get('old_y2_epsilon_transport_alpha_inputform')!=expected_alpha or rec.get('native_y1_to_new_basis_alpha_inputform')!='1':
        raise ValueError(f'{model}: Yukawa transport coefficient changed')
    if rec.get('uv_y2_explicit_prefactor_inputform')!= ('1' if d==1 else '-1'):
        raise ValueError(f'{model}: missing signed y2 UV factor')
    if d==1:
        # The Matchete singlet B native Yukawa CGs are archived exactly.
        for name,y in [('T3MENativeY1CG',y1),('T3MENativeY2CG',y2)]:
            if name not in registry or relative_error(y,registry[name])>1e-12:
                raise ValueError(f'{model}: native Yukawa CG differs from matched-source registry')
        k=np.einsum('ia,jb->ijab',y1.conj(),y2.conj())
        return {'model':model,'status':'B_YukawaCGMatchesMatchedRegistry',
                'native_yukawa_cg_equal_to_matched_source':True,
                'mixed_channel_tensor_norm':float(np.linalg.norm(k)),
                'native_y2_uv_scalar_transport_sign':1,
                'complete_feynman_amplitude_certified':False}
    q=registry['T3ExtraCGQ020020I1'].conj()
    b=registry['T3MEBridgeBB']
    eps=registry['eps[SU2L]'].conj()
    if q.shape!=(2,2,2,2) or b.shape!=(2,2) or eps.shape!=(2,2):
        raise ValueError('C matched SU2 invariant dimensions wrong')
    # Source 23/24 actually uses BridgeBB[scalar2,scalar1] and
    # Bar[eps[SU2L]][lepton1,lepton2]. Keep this exact order.
    e=np.einsum('ij,ba->ijab',eps,b)
    base=[q,e]
    if not np.allclose(y1,y2,atol=1e-12,rtol=1e-12):
        raise ValueError('C native Y1/Y2 tensor conventions differ unexpectedly')
    kval={}
    for label,aa,bb in [('11',y1,y1),('12',y1,y2),('22',y2,y2)]:
        k=contraction(aa,bb)
        kval[label]=project(k,base)
        if kval[label]['relative_residual']>1e-11:
            raise ValueError(f'C {label}: native Yukawa product outside matched Q/bridge invariant span')
        expected=[1/math.sqrt(3),math.sqrt(3)/2]
        observed=kval[label]['real']
        if max(abs(observed[n]-expected[n]) for n in (0,1))>1e-10 or max(abs(x) for x in kval[label]['imag'])>1e-10:
            raise ValueError(f'C {label}: signed Q/bridge normalization differs from native source expectation: {observed}')
    return {'model':model,'status':'C_YukawaProductsMatchSignedQAndBridgeFactors',
            'basis':['Bar[T3ExtraCGQ020020I1]','Bar[eps[SU2L]](i,j)*T3MEBridgeBB(b,a)'],
            'native_yukawa_product_projections':kval,
            'native_y2_uv_scalar_transport_sign':-1,
            'all_three_yukawa_products_in_native_invariant_span':True,
            'Q_factor_relative_to_raw_yukawa_product':'1/sqrt(3)',
            'bridge_factor_relative_to_raw_yukawa_product':'sqrt(3)/2',
            'signed_two_tensor_identity_verified':True,
            'explicit_C23_C24_bridge_contribution_compared_to_full_UV_exchange':False,
            'complete_feynman_amplitude_certified':False}


def load_separate_reports(paths):
    """Join two fresh-kernel evidence records, rejecting stale/wrong models."""
    if set(paths) != set(EXPECTED):
        raise ValueError('Exactly B_-1 and C_-1 report paths are required')
    models = {}
    hashes = {}
    for model in EXPECTED:
        path = Path(paths[model])
        obj = read(path)
        if obj.get('format') != 't3-bc-independent-yukawa-cg-components-v1' or \
           obj.get('model') != model or \
           obj.get('status') != 'IndependentNativeYukawaCGComponentsExtracted_SingleModel':
            raise ValueError(f'{model}: incorrect model/format/single-kernel status')
        for flag in ('master_weyl_phase_certified',
                     'native_BC_full_exchange_amplitude_computed',
                     'previous_release_changed',
                     'native_cg_product_computed_in_python'):
            if obj.get(flag) is not False:
                raise ValueError(f'{model}: unexpected certification flag {flag}')
        if not isinstance(obj.get('record'), dict):
            raise ValueError(f'{model}: missing extracted Yukawa tensor record')
        models[model] = obj['record']
        hashes[model] = sha(path)
    joined_sha = hashlib.sha256(('B_-1:' + hashes['B_-1'] + '\n' +
                                 'C_-1:' + hashes['C_-1'] + '\n').encode('ascii')).hexdigest()
    return {'format': 't3-bc-independent-yukawa-cg-components-v1',
            'status': 'IndependentNativeYukawaCGComponentsExtracted',
            'models': models,
            'master_weyl_phase_certified': False,
            'native_BC_full_exchange_amplitude_computed': False,
            'previous_release_changed': False}, joined_sha, hashes


def run(path=None, check_only=False, model_paths=None):
    if model_paths is None:
        raw=read(path)
        report_hash=sha(path)
        input_hashes=None
    else:
        raw, report_hash, input_hashes=load_separate_reports(model_paths)
    if raw.get('format')!='t3-bc-independent-yukawa-cg-components-v1' or raw.get('status')!='IndependentNativeYukawaCGComponentsExtracted':
        raise ValueError('Native extractor report missing/blocked')
    if raw.get('master_weyl_phase_certified') is not False or raw.get('previous_release_changed') is not False or raw.get('native_BC_full_exchange_amplitude_computed') is not False:
        raise ValueError('Input claims physical certification without supporting amplitude')
    if set(raw.get('models',{}))!=set(EXPECTED):
        raise ValueError('Missing or extra native B/C models')
    prior=bc.source_audit(ROOT)
    models={}
    for model in EXPECTED:
        reg,reg_sha=cg_registry(ROOT,model)
        models[model]=analyze(model,raw['models'][model],reg)
        models[model]['matched_cg_registry_sha256']=reg_sha
        models[model]['source_operators_checked']=prior['models'][model]['pl_source_terms_checked']
    result={'format':'t3-bc-independent-cg-product-audit-v1',
            'status':'BC_NativeCGProductDecompositionVerified_FullAmplitudeAndWeylPhaseOpen',
            'source_report_sha256':report_hash,
            'formal_release_summary_sha256':prior['release_summary_sha256'],
            'models':models,
            'native_BC_full_feynman_rule_exchange_computed':False,
            'C_bridge_source_amplitude_coefficient_certified':False,
            'matchete_gamma_cc_to_master_weyl_phase_certified':False,
            'full_physical_16_model_normalization_certified':False,
            'previous_release_changed':False}
    if input_hashes is not None:
        result['source_individual_report_sha256'] = input_hashes
    print('PASS: B native Yukawa CGs match Matchete matched registry')
    print('PASS: C native Yukawa CG products reproduce Q/sqrt(3) + sqrt(3)/2 bridge exactly')
    print('C projections (Q, bridge):',models['C_-1']['native_yukawa_product_projections'])
    print('OPEN: C source 23/24 full exchange coefficient; native Majorana amplitudes; master Weyl phase')
    if check_only:
        print('CHECK ONLY: no output written')
        return result
    encoded=json.dumps(result,sort_keys=True,indent=2)+'\n'
    out=PROBES/('bc_native_cg_products_'+report_hash[:16]+'.json')
    PROBES.mkdir(parents=True,exist_ok=True)
    if out.exists():
        if out.read_text(encoding='utf8')!=encoded:
            raise ValueError('Refusing to replace an existing differently-valued report')
        print('UNCHANGED:',out)
    else:
        with out.open('x',encoding='utf8') as f:f.write(encoded)
        print('CREATED:',out)
    return result


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--report',type=Path,help='Legacy combined B/C report')
    group.add_argument('--report-b',type=Path,help='Fresh-kernel B_-1 report (also provide --report-c)')
    ap.add_argument('--report-c',type=Path,help='Fresh-kernel C_-1 report; requires --report-b')
    ap.add_argument('--check-only',action='store_true')
    args=ap.parse_args(argv)
    try:
        if args.report_b is not None and args.report_c is None:
            raise ValueError('--report-c required with --report-b')
        if args.report is not None and args.report_c is not None:
            raise ValueError('--report-c requires --report-b, not --report')
        inputs={'B_-1': args.report_b, 'C_-1': args.report_c} if args.report_b is not None else None
        run(args.report,args.check_only,inputs)
        return 0
    except (AssertionError,KeyError,ValueError,TypeError,OSError) as e:
        print(f'BLOCKED: {type(e).__name__}: {e}',file=sys.stderr)
        return 2

if __name__=='__main__':raise SystemExit(main())
