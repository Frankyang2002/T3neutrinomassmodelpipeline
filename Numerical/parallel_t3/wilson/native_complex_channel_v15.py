"""Read-only v15 report on native complex-quartic/conjugate-channel observations.

Unsuccessful probes are observations, not physical coefficients. In particular
RGBeta's numerical S2 frame and the Matchete non-portal coupling map remain
unverified even if all native symbolic channels are resolved.
"""
from __future__ import annotations
from collections import Counter
import math
from pathlib import Path
from .native_context_quartic_v14 import (
    audit_witness, read, file_hash, FORMAT, V14Blocked,
)

STATUS = 'NativeRGBetaComplexQuartic28SourceConjugateChannelsObserved_Conditional'
BLOCKED = 'NativeRGBetaComplexQuarticChannelsUnresolved_Conditional'


def audit_from_native(root: Path, native_path: Path) -> dict:
    root=root.resolve()
    native_path=native_path if native_path.is_absolute() else root/native_path
    native=read(native_path)
    summary=audit_witness(native,root,allow_blocked=True)
    summary_rows={}
    counts=Counter()
    for model,model_row in sorted(native['models'].items()):
        rows=[]
        for row in model_row['native_complex_conjugate_projector_probes']:
            coupling=row['symbol']
            obs=row['native_complex_conjugate_response']
            counts['attempted']+=1
            if not isinstance(obs,dict) or obs.get('ok') is not True:
                tag='native_conjugate_channel_unresolved'
                reasons=str(obs.get('reason','missing_native_response')) if isinstance(obs,dict) else 'malformed_native_response'
            else:
                a,b=obs.get('direct_channel'),obs.get('conjugate_channel')
                if (type(a) not in (float,int) or type(b) not in (float,int) or
                        not math.isfinite(a) or not math.isfinite(b)):
                    tag='non_numeric_native_channel';reasons='missing_or_nonfinite'
                elif abs(a-1.)<2e-8 and abs(b)<2e-8 and obs.get('unit_direct_zero_conjugate') is True:
                    tag='native_source_unit_conjugate_zero';reasons=None
                    counts['unit_direct_zero_conjugate']+=1
                else:
                    tag='native_complex_mixed_response';reasons='observed_not_unit_direct_zero_conjugate'
            rows.append({'symbol':coupling,'classification':tag,'reason':reasons,
                         'native_response':obs})
        summary_rows[model]=rows
    counts=dict(counts)
    fully_observed=(counts.get('unit_direct_zero_conjugate',0)==counts['attempted'])
    return {'format':'t3-parallel-a-native-rgbeta-complex-conjugate-channels-v15',
            'status':STATUS if fully_observed else BLOCKED,
            'smoke_only':summary['smoke_only'],
            'native_v14_input_sha256':file_hash(native_path),
            'source_sha256':native['source_sha256'],
            'counts':counts,'models':summary_rows,
            'native_complex_direct_conjugate_symbolic_verified':fully_observed,
            'native_RGBeta_numeric_S2_basis_verified':False,
            'Matchete_complex_quartic_coupling_dictionary_verified':False,
            'native_complex_action_physical_normalization_certified':False,
            'physical_Wilson_evolution_authorized':False,
            'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT'}
