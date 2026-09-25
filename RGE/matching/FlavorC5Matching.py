"""Compatibility surface for the former mixed flavor-matching module.

New code should import one-generation flavor lifting from
``RGE.matching.WeinbergFlavorMatching`` and hierarchical final-C5 adaptation
from ``RGE.matching.FinalWeinbergAdapter``.
"""

from RGE.matching.FinalWeinbergAdapter import (
    is_final_weinberg_json,
    is_hierarchical_final_c5,
    load_final_weinberg_flavor_matrix,
    load_hierarchical_majorana_c5,
)
from RGE.matching.WeinbergFlavorMatching import (
    build_flavor_c5_from_matchete,
    build_flavor_c5_matrix,
    build_majorana_c5_flavor_matrix,
    extract_one_generation_loop_kernel,
    extract_t3_loop_kernel,
    match_c5_flavor_from_file,
    symbolic_t3_yukawas,
    write_flavor_matching_outputs,
)

__all__ = [
    "build_flavor_c5_from_matchete",
    "build_flavor_c5_matrix",
    "build_majorana_c5_flavor_matrix",
    "extract_one_generation_loop_kernel",
    "extract_t3_loop_kernel",
    "is_final_weinberg_json",
    "is_hierarchical_final_c5",
    "load_final_weinberg_flavor_matrix",
    "load_hierarchical_majorana_c5",
    "match_c5_flavor_from_file",
    "symbolic_t3_yukawas",
    "write_flavor_matching_outputs",
]
