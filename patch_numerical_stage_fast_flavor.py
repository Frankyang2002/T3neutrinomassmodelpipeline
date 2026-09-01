from pathlib import Path

path = Path("NumericalPipelineStage.py")

if not path.exists():
    raise FileNotFoundError(path)

text = path.read_text(encoding="utf-8")

old_import = "from FlavorMatchedC5 import flavor_match_from_c5_file\n"
new_import = (
    "from FlavorMatchedC5 import extract_t3_loop_kernel\n"
    "from MatchedEFTRGE import parse_matchete_c5\n"
)

if old_import not in text:
    raise RuntimeError("Could not locate FlavorMatchedC5 import.")

text = text.replace(old_import, new_import, 1)

start_marker = "def evaluate_symbolic_c5(\n"
end_marker = "\n\ndef run_numerical_pipeline_stage(\n"

start = text.find(start_marker)
end = text.find(end_marker)

if start == -1 or end == -1 or end <= start:
    raise RuntimeError("Could not locate evaluate_symbolic_c5 block.")

new_function = '''def evaluate_symbolic_c5(
    c5_path: Path,
    config: dict,
) -> np.ndarray:
    """Evaluate matched flavor C5 without building a symbolic 3x3 matrix."""

    c5_path = Path(c5_path)

    kappa_1g = parse_matchete_c5(
        c5_path.read_text(encoding="utf-8")
    )
    kernel = extract_t3_loop_kernel(kappa_1g)

    model = config["t3"]
    heavy_masses = np.asarray(model["MF"], dtype=float)

    y1 = _complex_matrix_from_json(
        model["y1_real"],
        model["y1_imag"],
    )
    y2 = _complex_matrix_from_json(
        model["y2_real"],
        model["y2_imag"],
    )

    if y1.shape != y2.shape:
        raise ValueError("y1 and y2 must have the same shape.")

    if y1.shape[0] != 3:
        raise ValueError(
            "T3 Yukawa matrices must have three lepton-flavor rows."
        )

    if y1.shape[1] != len(heavy_masses):
        raise ValueError(
            "The number of Yukawa columns must match the number of heavy masses."
        )

    base_substitutions = {
        sp.Symbol("MS1"): model["MS1"],
        sp.Symbol("MS2"): model["MS2"],
        sp.Symbol("lambdaT3"): model["lambdaT3"],
        sp.Symbol("hbar"): model.get(
            "hbar",
            1.0 / (16.0 * np.pi**2),
        ),
    }

    MF = sp.Symbol("MF")
    kernel_values = np.empty(len(heavy_masses), dtype=complex)

    for r, mass in enumerate(heavy_masses):
        kernel_r = kernel.subs(
            {
                **base_substitutions,
                MF: float(mass),
            }
        )
        kernel_values[r] = complex(sp.N(kernel_r, 18))

    K = np.zeros((3, 3), dtype=complex)

    for p in range(3):
        for q in range(p, 3):
            value = 0.0j

            for r in range(len(heavy_masses)):
                value += 0.5 * kernel_values[r] * (
                    np.conjugate(y1[p, r])
                    * np.conjugate(y2[q, r])
                    + np.conjugate(y2[p, r])
                    * np.conjugate(y1[q, r])
                )

            K[p, q] = value
            K[q, p] = value

    return K
'''

text = text[:start] + new_function + text[end:]

path.write_text(text, encoding="utf-8")

print("Patched:", path.resolve())
print("Fast numerical flavor evaluation enabled.")
