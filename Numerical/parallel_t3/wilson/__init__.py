"""Workstream A: conditional, source-checked full EFT1 Wilson transport.

No physical one-loop normalization or complete Matchete/RGBeta map is certified.
"""
from .basis import ModelBasis, ALL_MODELS, build_basis
from .coordinates import embed_c5, project_c5, embed_formal_llss
from .adm import beta, gauge_ward_residual, SYMMETRIC_TERMS
from .workflow import run_wilson, synthetic_case, run_from_files, BlockedInput

__all__ = ["ModelBasis", "ALL_MODELS", "build_basis", "embed_c5",
           "project_c5", "embed_formal_llss", "beta", "gauge_ward_residual",
           "SYMMETRIC_TERMS", "run_wilson", "synthetic_case", "run_from_files",
           "BlockedInput"]
