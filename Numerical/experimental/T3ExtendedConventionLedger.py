"""Source provenance and physical-normalization blockers for 6 extended T3 models.

Checks naming and scalar monomial correspondences; explicitly does NOT certify
native SU2 coordinate phases, fermion charge conjugation, Dirac mass phase,
Matchete-to-RGBeta overall coupling factors or flavour threshold matching.
"""
from __future__ import annotations
from pathlib import Path
import hashlib
import json
import re
from typing import Any

from Numerical.experimental import T3RGAssistedExtendedC5 as BR
from Numerical.experimental import T3ExtendedC5ExpressionReduction as C5

ROOT=Path(__file__).resolve().parents[2]
BASE=Path('Lagrangian/interactions/T3Topology.wl')
SCALAR=Path('Lagrangian/interactions/GeneralScalarOperators.wl')
RG=Path('RGE/running/rgbeta/T3RGBetaModel.wl')


def require(ok:bool, message:str)->None:
    if not ok:raise ValueError(message)


def audited_source(path:Path,words:list[str])->dict[str,Any]:
    require(path.is_file(),'Missing source '+str(path))
    data=path.read_text(encoding='utf-8-sig')
    for word in words:
        require(word in data, 'Changed source contract: '+str(path)+' expected '+word)
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def audit(root:Path=ROOT, *,native_source:Path|None=None)->dict[str,Any]:
    root=Path(root)
    sources={
        'Matchete_T3_topology':audited_source(root/BASE, ['lambdaT3[] H[i] H[j] s1 s2 CG[T3MixCG, labels]',
                                                    '"Expression" -> PlusHc[']),
        'Matchete_general_quartic':audited_source(root/SCALAR, ['"Expression" -> PlusHc[contribution]',
                                                           'contribution = coupling[] product If[cg === None, 1, CG[cg, labels]]']),
        'RGBeta_T3_model':audited_source(root/RG,['T3RGBetaAddMixingQuartic[','lambdaT3,',
                                                   '{H, H, S1, Bar @ S2}','T3RGBetaYukawaInvariant[']),
    }
    container=Path(native_source) if native_source is not None else root/'Reports/output/workstream_a/extended_c5_native/native'
    rows=[]
    for model,(matchete_quartic,rgbeta_quartic) in BR.EXTRA_QUARTIC.items():
        native=C5.verify_model(container,model)
        sourcefile=root/'RGE/general/generated'/f'{model[0]}_{(0 if model.endswith("p0") else -2):+d}_extra_quartics.wl'
        text=sourcefile.read_text(encoding='utf-8-sig') if sourcefile.is_file() else ''
        require(text!='','Missing extra-quartic RGBeta source '+str(sourcefile))
        suffix=model.split('_')[1]
        word='{Bar @ H, Bar @ H, S1, S2}' if suffix=='p0' else '{Bar @ H, Bar @ H, Bar @ S1, Bar @ S2}'
        require(re.search(r'AddQuartic\[\s*'+re.escape(rgbeta_quartic)+r'\s*,\s*'+re.escape(word)+r'\s*,\s*GroupInvariant\s*->',text,re.S) is not None,
                'Extra-quartic RGBeta field monomial changed for '+model)
        require(native['extra_quartic_matchete']==matchete_quartic,'Matchete quartic name mismatch '+model)
        rows.append({'model':model,'extra_yukawa':native['extra_yukawa'],
                     'Matchete_extra_quartic':matchete_quartic,'RGBeta_extra_quartic':rgbeta_quartic,
                     'RGBeta_field_word':word,
                     'RGBeta_source_sha256':hashlib.sha256(sourcefile.read_bytes()).hexdigest(),
                     'native_Matchete_export_digests':native['source_sha256'],
                     'native_exact_one_step_C5_reduced':True,
                     'field_word_correspondence_established':True,
                     'native_SU2_tensor_coordinate_ratio_certified':False,
                     'physical_complex_coupling_conversion_certified':False})
    prior=root/'Reports/output/workstream_a/yukawa_spinor_phase/a0_spinor_phase.json'
    phase_status='NotPresent'
    if prior.exists():
        p=json.loads(prior.read_text(encoding='utf-8-sig'))
        require(p.get('physical_cross_package_yukawa_dictionary_certified') is False,
                'Prior spinor phase evidence unexpectedly claims a physical dictionary')
        phase_status=str(p.get('status','Unknown'))
    return {
        'status':'SixModelNativeC5AndSourceFieldWordsCheckedPhysicalDictionaryBlocked',
        'sources':sources,'models':rows,'n_native_matchete_reductions':len(rows),
        'prior_A0_spinor_phase_status':phase_status,
        'candidate_relative_Lagrangian_sign':
            'Matchete topology/extra scalar interactions PlusHc(+coupling*monomial); RGBeta scalar sector has conventional minus sign. A sign candidate alone is not sufficient to transport field coordinate phases/CG normalization.',
        'required_remaining_native_certificates':[
            'Same-component Matchete/RGBeta Clebsch-Gordan tensor normalization and phase for each A/D/E scalar vertex',
            'Absolute Weyl/Dirac charge-conjugation and heavy-flavour mass-basis map for all four Yukawas',
            'S1/S2 triplet basis conversion for all A/D/E models, including conjugate representations',
            'Full EFT1 fermion/scalar sequential matching with new Yukawa and quartic operators',
            'Counterterm and matching-scale cancellation through the declared loop order'
        ],
        'physical_coupling_dictionary_certified':False,
        'hierarchical_threshold_matching_certified':False,
        'physical_neutrino_mass_prediction_certified':False,
        'production_files_changed':False}
