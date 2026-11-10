"""Source-cached F-only Matchete LLHH Laurent kernel inventory for all 16.

This is a NON-INJECTED source-coverage register. Its native kernels contain
pole, finite and log terms but were not UV-minus-EFT1-soft/counterterm matched.
They must not be added on top of the analytic hard triangle: that would risk
both double-counting and using an unrenormalized bare coefficient.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import json
import re

ROOT_REL=Path('output/general_observable/lo_certification_spinor/pole_log_e515d82d118284ba_bridge01')
DICT_REL=Path('output/parallel_t3/wilson/additional_quartic_dictionary_conditional_v3.json')
FORMAT='t3-source-native-LLHH-laurent-coverage-batch3-v1'
COUPLING=re.compile(r'Coupling\[([A-Za-z][A-Za-z0-9]*),')


def source_inventory(root:Path,model:str,slug:str,available_eft1_symbols:set[str]) -> dict:
    src=Path(root)/ROOT_REL/slug/'trios.json'
    data=json.loads(src.read_text(encoding='utf-8-sig'))
    if data.get('format')!='t3-f-first-source-resolved-pole-log-bridge-v1' or data.get('model')!=model:
        raise ValueError('Wrong native Matchete one-loop Laurent inventory: '+model)
    if (data.get('EFT1_soft_subtraction_performed') is not False or
        data.get('counterterm_subtraction_performed') is not False or
        data.get('physical_MSbar_hard_coefficient') is not None):
        raise ValueError('Unexpected native first-threshold matching provenance flags')
    records=data.get('records')
    if not isinstance(records,list) or len(records)!=data.get('exact_kernel_groups'):
        raise ValueError('Missing source-native Laurent groups')
    seen=set();pl=[];pr=[];names=set()
    for v in records:
        if v.get('chirality') not in ('PL','PR') or v.get('coefficient_factor')!='1 + 1/ϵ + Log[μbar2/MF^2]':
            raise ValueError('Unrecognized native Laurent kernel')
        kernel=v['kernel_inputform']
        digest=sha256(kernel.encode('utf8')).hexdigest()
        if digest!=v.get('kernel_sha256') or v.get('pole_residue_equals_constant_kernel') is not True or v.get('log_scale_coefficient_equals_pole_residue') is not True:
            raise ValueError('Native Laurent kernel identity/source hash failure')
        sidx=v.get('source_indices',{})
        if set(sidx)!={'constant','pole','scale_log'} or any(type(x)!=int for x in sidx.values()):
            raise ValueError('Incorrect source term provenance indexes')
        for idx in sidx.values():
            if idx in seen:raise ValueError('Duplicate native Laurent source term')
            seen.add(idx)
        if 'GammaCC' not in kernel or 'Coupling[MF,' not in kernel or 'hbar' not in kernel:
            raise ValueError('Not an F-only one-loop LLHH source kernel')
        present=set(COUPLING.findall(kernel))
        quartic={s for s in present if s.startswith('lambda')}
        names.update(quartic)
        record={'kernel_sha256':digest,'quartic_symbols':sorted(quartic),
                'mixing_spurions':sorted(s for s in present if s.startswith('T3MEU')),
                'source_indices':sidx,'native_scale_log':v.get('mass_log_inputform')}
        (pl if v['chirality']=='PL' else pr).append(record)
    if len(pl)!=len(pr) or len(pl)!=data['PL_groups'] or len(pr)!=data['PR_groups']:
        raise ValueError('Native PL/PR inventory incomplete')
    mapping={}
    dictpath=Path(root)/DICT_REL
    if dictpath.exists():
        d=json.loads(dictpath.read_text(encoding='utf-8-sig'))
        entry=d.get('models',{}).get(model,{})
        matnames=entry.get('matchete_symbols',[])
        rgbnames=entry.get('rgbeta_symbols',[])
        matrix=entry.get('map_rgbeta_from_matchete',{})
        real=matrix.get('real',[])
        imag=matrix.get('imag',[])
        if matnames:
            if (len(real)!=len(rgbnames) or len(imag)!=len(rgbnames) or
                any(len(row)!=len(matnames) for row in real+imag)):
                raise ValueError('Conditional source quartic dictionary matrix malformed')
            # Source convention: c_RGBeta = A @ c_Matchete.
            # Non-diagonal blocks MUST NOT be collapsed to zipped aliases.
            for j,mat in enumerate(matnames):
                mapping[mat]=[rgb for i,rgb in enumerate(rgbnames)
                    if abs(float(real[i][j]))+abs(float(imag[i][j]))>1e-10]
    for s in names:
        if s=='lambdaT3':mapping[s]=['lambdaT3']
    represented=sorted(s for s in names if s in mapping)
    missing=sorted(s for s in names if s not in mapping)
    return {'format':FORMAT,'model':model,'source_path':str(src.relative_to(root)).replace('\\','/'),
            'source_sha256':sha256(src.read_bytes()).hexdigest(),
            'num_native_PL_kernels':len(pl),'num_native_PR_kernels':len(pr),
            'complete_source_laurent_triplets':True,
            'native_quartics':sorted(names),
            'quartics_with_conditional_dictionary_record':represented,
            'unmapped_native_quartics':missing,
            'conditional_RGBeta_symbol_candidates':mapping,
            'candidate_RGBeta_symbols_available_at_MF':sorted(set(available_eft1_symbols).intersection(
                {x for targets in mapping.values() for x in targets})),
            'PL_source_records':pl,
            'requires_EFT1_soft_subtraction':True,
            'requires_counterterms':True,
            'native_GammaCC_Weyl_sign_unfixed':True,
            'raw_native_terms_added_to_C5':False,
            'complete_one_loop_F_matching_claimed':False}
