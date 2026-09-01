from pathlib import Path

path = Path("GeneralWeinbergRGEGenerator.py")

if not path.exists():
    raise FileNotFoundError(path)

text = path.read_text(encoding="utf-8")

# Add a new helper without changing the existing permutation helper/comment.
anchor = """def _pair_product_permutations(a, b, i, j):
    \"\"\"Return the four independent swaps in sigma({a,b} x {i,j}).\"\"\"

    return tuple(
        (a_perm, b_perm, i_perm, j_perm)
        for a_perm, b_perm in _swap_pair(a, b)
        for i_perm, j_perm in _swap_pair(i, j)
    )


"""

addition = """def _simultaneous_pair_permutations(a, b, i, j):
    \"\"\"Return the two simultaneous pair exchanges used by Eq. (4.85).\"\"\"

    return (
        (a, b, i, j),
        (b, a, j, i),
    )


"""

if "_simultaneous_pair_permutations" not in text:
    if anchor not in text:
        raise RuntimeError("Could not locate _pair_product_permutations.")
    text = text.replace(anchor, anchor + addition, 1)

# Replace exactly the two Eq. (4.85) product-permutation uses.
old_loop = """for a_perm, b_perm, i_perm, j_perm in _pair_product_permutations(
        a, b, i, j
    ):"""

new_loop = """for a_perm, b_perm, i_perm, j_perm in _simultaneous_pair_permutations(
        a, b, i, j
    ):"""

count = text.count(old_loop)

if count != 2:
    raise RuntimeError(
        f"Expected exactly 2 product-permutation uses, found {count}."
    )

text = text.replace(old_loop, new_loop)

# Add the missing gauge-coupling factor to gamma_cf.
old_gamma = """        result += -3 * casimir[i - 1, j - 1]
"""

new_gamma = """        result += (
            -3
            * sector.coupling**2
            * casimir[i - 1, j - 1]
        )
"""

if old_gamma not in text:
    raise RuntimeError(
        "Could not locate fermion anomalous-dimension gauge term."
    )

text = text.replace(old_gamma, new_gamma, 1)

path.write_text(text, encoding="utf-8")

print("Patched:", path.resolve())
print(
    "simultaneous permutation uses:",
    text.count("_simultaneous_pair_permutations(") - 1,
)
print(
    "fermion gamma contains gauge coupling:",
    "* sector.coupling**2" in text,
)
