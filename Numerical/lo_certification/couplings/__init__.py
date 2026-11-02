"""Fail-closed T3 EFT1 coupling and scalar-basis handoff; no physics aliases inferred."""
from .basis import (orthogonal_matrix, transport_mass2, transport_hhss,
                    transport_llss, relative_basis_rotation, self_bilinear_hessian)

__all__ = ['orthogonal_matrix', 'transport_mass2', 'transport_hhss',
           'transport_llss', 'relative_basis_rotation', 'self_bilinear_hessian']
