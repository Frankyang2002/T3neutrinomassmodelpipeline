"""Stage 15: source-checked B/C(-1) RGBeta-to-Matchete quartic convention.

RGBeta's published model-builder Lagrangian (Thomsen, EPJC 81 (2021) 408,
Eq. 2.3) contains -sum(lambda phi phi phi phi). Matchete's *project*
experimental extra-quartic builder instead adds PlusHc[lambda*fields*CG]
with a plus sign. The two field words are matched explicitly here, and the
previous independent SU(2) tensor comparison is required to be +1.

This establishes the source-level convention for the two named extra
quartics in a common scalar-field basis. It does NOT validate their
threshold-scale running, the complete matching scheme, or other couplings.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping

from Numerical.experimental.GeneralT3TensorMap import MAPPING, MODELS

RGBETA_MBL_REFERENCE = 'https://doi.org/10.1140/epjc/s10052-021-09142-4'
RGBETA_MBL_EQUATION = 'Equation (2.3): L contains - sum_couplings lambda_abcd phi_a phi_b phi_c phi_d'
MAPPING_FIELDS = {
    'lambdaT3ExtraQ020020I1': ('lambdaHHdagS2S2', ('Hdag', 'Hdag', 'S2', 'S2')),
    'lambdaT3ExtraQ020200I1': ('lambdaHHdagS1barS1bar', ('Hdag', 'Hdag', 'S1dag', 'S1dag')),
}


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def _read(path: Path) -> str:
    return Path(path).read_text(encoding='utf-8-sig')


def _check_matchete_source(text: str) -> None:
    checks = {
        'extra quartics are positive PlusHc[contribution]': r'"Expression"\s*->\s*PlusHc\[contribution\]',
        'contribution is coupling times field product and CG': r'contribution\s*=\s*coupling\[\]\s*product\s*If\[cg\s*===\s*None\s*,\s*1\s*,\s*CG\[cg\s*,\s*labels\]\]',
        'quartics use Q family': r'prefix\s*=\s*If\[degree\s*===\s*4\s*,\s*"Q"\s*,\s*"M"\]',
        'extra couplings explicitly complex': r'DefineCoupling\[coupling\s*,\s*SelfConjugate\s*->\s*False\]',
        'Hdag is a barred H field': r'"Hdag"\s*,\s*Bar\[H\[index\]\]',
        'S1dag is barred Scalar1': r'"S1dag"\s*,\s*ScalarPotentialFieldValue\[1\s*,\s*True',
        'S2 is unbarred Scalar2': r'"S2"\s*,\s*ScalarPotentialFieldValue\[2\s*,\s*False',
    }
    for label, pattern in checks.items():
        if re.search(pattern, text, re.S) is None:
            raise ValueError(f'Matchete source contract changed: {label}')
    # Avoid silently accepting a new prefactor or sign multiplying the interaction.
    if re.search(r'"Expression"\s*->\s*-\s*PlusHc', text):
        raise ValueError('Matchete extra-quartic sign changed')


def _check_rgbeta_source(text: str) -> None:
    for mat, (rgb, fields) in MAPPING_FIELDS.items():
        field_words = {
            'Hdag': 'Bar\\s*@\\s*H',
            'S1dag': 'Bar\\s*@\\s*S1',
            'S2': 'S2',
        }
        expression = r'\s*,\s*'.join(field_words[name] for name in fields)
        pattern = (r'AddQuartic\s*\[\s*' + re.escape(rgb) + r'\s*,\s*\{\s*'
                   + expression + r'\s*\}\s*,\s*GroupInvariant\s*->\s*\('
                   + r'\s*delS2\[SU2L,\s*a,\s*#1,\s*#2\]\s*\*'
                   + r'\s*delS2\[SU2L,\s*a,\s*#3,\s*#4\]\s*&\s*\)\s*,'
                   + r'\s*SelfConjugate\s*->\s*False')
        if re.search(pattern, text, re.S) is None:
            raise ValueError(f'RGBeta source contract changed for {mat} / {rgb}')


def _check_tensor_report(obj: Mapping[str, Any]) -> None:
    if obj.get('status') != 'DoubletGroupTensorRatiosVerified':
        raise ValueError('Stage 13 tensor report must be verified')
    if obj.get('full_lagrangian_coupling_identification_verified') is not False:
        raise ValueError('Expected the earlier independent tensor-only report')
    expected = {(model, mat): rgb for model in MODELS for mat, (rgb, _) in MAPPING_FIELDS.items()}
    received = {}
    for row in obj.get('mappings', []):
        key = (row.get('model'), row.get('matchete_coupling'))
        if key not in expected or key in received:
            raise ValueError('Unexpected or duplicate tensor mapping')
        if row.get('rgbeta_coupling') != expected[key]:
            raise ValueError('Tensor mapping has wrong RGBeta field structure')
        if float(row.get('rgbeta_over_matchete_tensor_factor', 0)) != 1.0:
            raise ValueError('Expected tensor ratio +1')
        received[key] = True
    if set(received) != set(expected):
        raise ValueError('Tensor report does not certify four B/C(-1) mappings')


def derive_lagrangian_convention(
    matchete_source: Path,
    rgbeta_source: Path,
    tensor_report: Path,
) -> dict[str, Any]:
    """Certify source-level relative sign, retaining explicit scope caveats."""
    mat_text, rgb_text = _read(matchete_source), _read(rgbeta_source)
    _check_matchete_source(mat_text)
    _check_rgbeta_source(rgb_text)
    tensor = json.loads(Path(tensor_report).read_text(encoding='utf-8-sig'))
    _check_tensor_report(tensor)
    pairs = [
        {
            'model': model,
            'matchete_coupling': mat,
            'rgbeta_coupling': rgb,
            'field_monomial': list(fields),
            'rgbeta_over_matchete_tensor_factor': 1,
            'matchete_in_terms_of_rgbeta': 'negative',
            'relative_sign': -1,
            'requires_complex_conjugation': False,
        }
        for model in MODELS
        for mat, (rgb, fields) in MAPPING_FIELDS.items()
    ]
    return {
        'status': 'SourceLagrangianConventionDerived',
        'scope': 'ordinary split-scalar B/C(alpha=-1), two C5-active quartics only',
        'matchete_lagrangian_term': '+lambda_mat * Q + H.c.',
        'rgbeta_model_builder_lagrangian_term': '-lambda_rgb * Q + H.c. with overall minus on the quartic sector',
        'reference': RGBETA_MBL_REFERENCE,
        'reference_equation': RGBETA_MBL_EQUATION,
        'matchete_source_sha256': _sha256(mat_text),
        'rgbeta_model_source_sha256': _sha256(rgb_text),
        'tensor_report_status': tensor['status'],
        'source_checks_passed': True,
        'results': pairs,
        'same_field_basis_assumed': True,
        'full_threshold_matching_verified': False,
        'physical_rg_improved_prediction': False,
        'disclaimer': (
            'The sign/conjugation map follows the two source Lagrangian conventions and independently '
            'verified SU(2) tensor ratios in the same named scalar-field basis. This is not a '
            'threshold-consistent coupled matching validation or a check of all T3 quartics.'
        ),
    }
