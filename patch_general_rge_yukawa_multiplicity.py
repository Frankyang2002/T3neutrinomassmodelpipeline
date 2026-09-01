from pathlib import Path

path = Path("GeneralWeinbergRGEGenerator.py")

if not path.exists():
    raise FileNotFoundError(path)

text = path.read_text(encoding="utf-8")

old = """    for a_perm, b_perm, i_perm, j_perm in _simultaneous_pair_permutations(
        a, b, i, j
    ):
        for k in range(1, nf + 1):
            for l in range(1, nf + 1):
                for c in range(1, ns + 1):
                    bracket = (
                        2
                        * y(j_perm, k, c)
                        * sp.conjugate(y(k, l, b_perm))
                        + y(j_perm, k, b_perm)
                        * sp.conjugate(y(k, l, c))
                    )

                    for sector in inputs.gauge_sectors:
                        for t_generator, theta in zip(
                            sector.fermion_generators,
                            sector.scalar_generators,
                        ):
                            bracket += (
                                4
                                * sector.coupling**2
                                * t_generator[l - 1, j_perm - 1]
                                * theta[b_perm - 1, c - 1]
                            )

                    result += (
                        bracket
                        * coefficient[i_perm, l, a_perm, c]
                    )

    return sp.simplify(result)
"""

new = """    # The Yukawa contractions require the four independent external-pair swaps.
    for a_perm, b_perm, i_perm, j_perm in _pair_product_permutations(
        a, b, i, j
    ):
        for k in range(1, nf + 1):
            for l in range(1, nf + 1):
                for c in range(1, ns + 1):
                    yukawa_bracket = (
                        2
                        * y(j_perm, k, c)
                        * sp.conjugate(y(k, l, b_perm))
                        + y(j_perm, k, b_perm)
                        * sp.conjugate(y(k, l, c))
                    )

                    result += (
                        yukawa_bracket
                        * coefficient[i_perm, l, a_perm, c]
                    )

    # In the real-scalar component implementation the gauge contraction is
    # already symmetric under the duplicated external-pair exchanges, so only
    # the two simultaneous exchanges are retained.
    for a_perm, b_perm, i_perm, j_perm in _simultaneous_pair_permutations(
        a, b, i, j
    ):
        for l in range(1, nf + 1):
            for c in range(1, ns + 1):
                gauge_bracket = sp.S.Zero

                for sector in inputs.gauge_sectors:
                    for t_generator, theta in zip(
                        sector.fermion_generators,
                        sector.scalar_generators,
                    ):
                        gauge_bracket += (
                            4
                            * sector.coupling**2
                            * t_generator[l - 1, j_perm - 1]
                            * theta[b_perm - 1, c - 1]
                        )

                result += (
                    gauge_bracket
                    * coefficient[i_perm, l, a_perm, c]
                )

    return sp.simplify(result)
"""

if old not in text:
    raise RuntimeError(
        "Could not locate the current mixed_yukawa_gauge_term loop body."
    )

text = text.replace(old, new, 1)
path.write_text(text, encoding="utf-8")

print("Patched:", path.resolve())
print(
    "mixed term has four-swap Yukawa loop:",
    "_pair_product_permutations(" in text[text.find("def mixed_yukawa_gauge_term"):text.find("def crossed_yukawa_term")]
)
print(
    "mixed term has two-swap gauge loop:",
    "_simultaneous_pair_permutations(" in text[text.find("def mixed_yukawa_gauge_term"):text.find("def crossed_yukawa_term")]
)
