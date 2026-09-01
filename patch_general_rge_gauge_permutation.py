from pathlib import Path

path = Path("GeneralWeinbergRGEGenerator.py")

if not path.exists():
    raise FileNotFoundError(path)

text = path.read_text(encoding="utf-8")

old = """    # In the real-scalar component implementation the gauge contraction is
    # already symmetric under the duplicated external-pair exchanges, so only
    # the two simultaneous exchanges are retained.
    for a_perm, b_perm, i_perm, j_perm in _simultaneous_pair_permutations(
        a, b, i, j
    ):
"""

new = """    # The gauge contraction carries the same four independent external-pair
    # swaps.  Unlike the Yukawa contraction above, it has no summed k index.
    for a_perm, b_perm, i_perm, j_perm in _pair_product_permutations(
        a, b, i, j
    ):
"""

if old not in text:
    raise RuntimeError(
        "Could not locate the separated gauge permutation block."
    )

text = text.replace(old, new, 1)

path.write_text(text, encoding="utf-8")

print("Patched:", path.resolve())

section = text[
    text.find("def mixed_yukawa_gauge_term"):
    text.find("def crossed_yukawa_term")
]

print(
    "four-swap loops in mixed term:",
    section.count("_pair_product_permutations("),
)
print(
    "simultaneous-swap loops in mixed term:",
    section.count("_simultaneous_pair_permutations("),
)
