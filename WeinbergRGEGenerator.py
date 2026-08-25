"""Symbolic helpers for converting the general psi^2 phi^2 RGE to the
scotogenic-model conventions used in the project.

Matrix convention: matrix[row, column].
Real scalar basis for each complex multiplet: (R1, I1, R2, I2, ...).
"""

import sympy as sp


# -----------------------------------------------------------------------------
# Symbols and field conventions
# -----------------------------------------------------------------------------

HALF = sp.Rational(1, 2)

g1, g2, g3 = sp.symbols("g1 g2 g3")
kappa = sp.Symbol("kappa")
C = sp.IndexedBase("C")
K11, K22 = sp.symbols("K11 K22")
lambda1, lambda2, lambda3, lambda4, lambda5 = sp.symbols(
    "lambda1 lambda2 lambda3 lambda4 lambda5"
)
LAMBDA = sp.IndexedBase("lambda")

FIELDS = {
    "L": {
        "id": 1,
        "symbol": sp.Symbol("L"),
        "kind": "fermion",
        "rep": (1, 2, -HALF),
        "components": {"nu": 1, "eL": 2},
    },
    "eR": {
        "id": 3,
        "symbol": sp.Symbol("e_R"),
        "kind": "fermion",
        "rep": (1, 1, -1),
    },
    "eta": {
        "id": 4,
        "symbol": sp.Symbol("eta"),
        "kind": "complex_scalar",
        "rep": (1, 2, HALF),
    },
    "NR": {
        "id": 5,
        "symbol": sp.Symbol("N_R"),
        "kind": "fermion",
        "rep": (1, 1, 0),
    },
    "H": {
        "id": 6,
        "symbol": sp.Symbol("H"),
        "kind": "complex_scalar",
        "rep": (1, 2, HALF),
    },
}

NU = FIELDS["L"]["components"]["nu"]
EL = FIELDS["L"]["components"]["eL"]

# Global real-scalar basis used by lambda_abcd and C_ijab:
#   eta -> 1..4, H -> 5..8.
SCALAR_GLOBAL_RANGES = {
    "eta": range(1, 5),
    "H": range(5, 9),
}
SCALAR_OFFSETS = {"eta": 0, "H": 4}
TOTAL_REAL_SCALAR_DIMENSION = 8

# Neutral complex component written in the real basis:
#   eta^0 = (R_3 + i I_4)/sqrt(2) -> K22
#   H^0   = (R_7 + i I_8)/sqrt(2) -> K11
NEUTRAL_C_TO_K = {
    (3, 3): K22 / 2,
    (3, 4): sp.I * K22 / 2,
    (4, 3): sp.I * K22 / 2,
    (4, 4): -K22 / 2,
    (7, 7): K11 / 2,
    (7, 8): sp.I * K11 / 2,
    (8, 7): sp.I * K11 / 2,
    (8, 8): -K11 / 2,
}

DEFAULT_OUTPUT_COMPONENT = (NU, NU, 3, 3)
# Backwards-compatible name used by the original script.
output_component = DEFAULT_OUTPUT_COMPONENT


# -----------------------------------------------------------------------------
# Gauge generators
# -----------------------------------------------------------------------------

def su2_complex_generators(dimension):
    """Return Hermitian SU(2) generators (T1, T2, T3) for dimension 2j+1."""
    if not isinstance(dimension, int) or dimension < 1:
        raise ValueError("SU(2) representation dimension must be a positive integer.")

    j = sp.Rational(dimension - 1, 2)
    m_values = [j - position for position in range(dimension)]
    t_plus = sp.zeros(dimension)

    # T_+ |j,m> = sqrt((j-m)(j+m+1)) |j,m+1>.
    for column, m in enumerate(m_values):
        target_m = m + 1
        if target_m in m_values:
            row = m_values.index(target_m)
            t_plus[row, column] = sp.sqrt((j - m) * (j + m + 1))

    t_minus = t_plus.conjugate().T
    t1 = (t_plus + t_minus) / 2
    t2 = (t_plus - t_minus) / (2 * sp.I)
    t3 = sp.diag(*m_values)
    return tuple(sp.simplify(generator) for generator in (t1, t2, t3))


def real_scalar_generator(generator):
    """Convert a complex generator T to theta acting on interleaved (R,I) fields.

    With phi_complex = (R+iI)/sqrt(2), the convention is

        delta phi_real = i epsilon theta phi_real.

    For T = Re(T) + i Im(T), each complex matrix element contributes the block

        [[ i Im(T),  i Re(T)],
         [-i Re(T),  i Im(T)]].
    """
    dimension = generator.rows
    theta = sp.zeros(2 * dimension)

    for p in range(dimension):
        for q in range(dimension):
            entry = sp.expand_complex(generator[p, q])
            real_part = sp.re(entry)
            imaginary_part = sp.im(entry)
            r_p, i_p = 2 * p, 2 * p + 1
            r_q, i_q = 2 * q, 2 * q + 1

            theta[r_p, r_q] += sp.I * imaginary_part
            theta[r_p, i_q] += sp.I * real_part
            theta[i_p, r_q] += -sp.I * real_part
            theta[i_p, i_q] += sp.I * imaginary_part

    return sp.simplify(theta)


def su2_real_scalar_generators(dimension):
    """Return the three SU(2) generators in the real-scalar basis."""
    return tuple(
        real_scalar_generator(generator)
        for generator in su2_complex_generators(dimension)
    )


def u1_real_scalar_generator(dimension, hypercharge):
    """Return the U(1)_Y generator in the real-scalar basis."""
    return real_scalar_generator(hypercharge * sp.eye(dimension))


# -----------------------------------------------------------------------------
# Scalar quartic tensor lambda_abcd
# -----------------------------------------------------------------------------

MIXED_QUARTIC_COMPONENTS = {
    (1, 1, 5, 5): lambda3 + lambda4 + lambda5,
    (1, 1, 6, 6): lambda3 + lambda4 - lambda5,
    (1, 1, 7, 7): lambda3,
    (1, 1, 8, 8): lambda3,
    (2, 2, 5, 5): lambda3 + lambda4 - lambda5,
    (2, 2, 6, 6): lambda3 + lambda4 + lambda5,
    (2, 2, 7, 7): lambda3,
    (2, 2, 8, 8): lambda3,
    (3, 3, 5, 5): lambda3,
    (3, 3, 6, 6): lambda3,
    (3, 3, 7, 7): lambda3 + lambda4 + lambda5,
    (3, 3, 8, 8): lambda3 + lambda4 - lambda5,
    (4, 4, 5, 5): lambda3,
    (4, 4, 6, 6): lambda3,
    (4, 4, 7, 7): lambda3 + lambda4 - lambda5,
    (4, 4, 8, 8): lambda3 + lambda4 + lambda5,
    (1, 2, 5, 6): lambda5,
    (3, 4, 7, 8): lambda5,
    (1, 3, 5, 7): (lambda4 + lambda5) / 2,
    (1, 3, 6, 8): (lambda4 - lambda5) / 2,
    (1, 4, 5, 8): (lambda4 + lambda5) / 2,
    (1, 4, 6, 7): -(lambda4 - lambda5) / 2,
    (2, 3, 5, 8): -(lambda4 - lambda5) / 2,
    (2, 3, 6, 7): (lambda4 + lambda5) / 2,
    (2, 4, 5, 7): (lambda4 - lambda5) / 2,
    (2, 4, 6, 8): (lambda4 + lambda5) / 2,
}


def _validate_scalar_indices(indices):
    if len(indices) != 4:
        raise ValueError("Exactly four scalar indices are required.")
    if any(not isinstance(index, int) for index in indices):
        raise TypeError("Scalar indices must be Python integers.")
    if any(index < 1 or index > TOTAL_REAL_SCALAR_DIMENSION for index in indices):
        raise IndexError("Scalar indices must lie between 1 and 8.")


def _quartic_self_coupling(key, coupling):
    """Return lambda1/lambda2 contribution for four fields from one doublet."""
    multiplicities = sorted(
        (key.count(index) for index in set(key)),
        reverse=True,
    )
    if multiplicities == [4]:
        return 3 * coupling
    if multiplicities == [2, 2]:
        return coupling
    return sp.S.Zero


def scalar_quartic_component(*indices):
    """Return symmetric scalar coupling lambda_abcd in the global real basis.

    eta occupies indices 1..4 and H occupies 5..8. Since lambda_abcd is fully
    symmetric, the sorted index tuple uniquely identifies each component.
    """
    _validate_scalar_indices(indices)
    key = tuple(sorted(indices))

    if all(index in SCALAR_GLOBAL_RANGES["eta"] for index in key):
        return _quartic_self_coupling(key, lambda2)
    if all(index in SCALAR_GLOBAL_RANGES["H"] for index in key):
        return _quartic_self_coupling(key, lambda1)
    return sp.simplify(MIXED_QUARTIC_COMPONENTS.get(key, sp.S.Zero))


def lambda_substitutions(lambda_tensor=LAMBDA):
    """Map all lambda[a,b,c,d] components in the 1..8 basis to lambda1..lambda5."""
    return {
        lambda_tensor[a, b, c, d]: scalar_quartic_component(a, b, c, d)
        for a in range(1, 9)
        for b in range(1, 9)
        for c in range(1, 9)
        for d in range(1, 9)
    }


def convert_lambda_expression(expression, lambda_tensor=LAMBDA):
    """Replace symbolic lambda[a,b,c,d] components by lambda1,...,lambda5."""
    substitutions = lambda_substitutions(lambda_tensor=lambda_tensor)
    return sp.simplify(sp.expand(expression.xreplace(substitutions)))


def print_lambda_component(*indices):
    """Print one lambda_abcd component and return its value."""
    value = scalar_quartic_component(*indices)
    label = "".join(str(index) for index in indices)
    print(f"lambda_{label} =")
    sp.pprint(value, use_unicode=True)
    return value


# -----------------------------------------------------------------------------
# RGE contributions
# -----------------------------------------------------------------------------

RGE_CONTRIBUTIONS = []


def register_rge_contribution(function):
    """Register a contribution with signature (scalar_field, output_component, C)."""
    RGE_CONTRIBUTIONS.append(function)
    return function


def _complex_scalar_field(scalar_field):
    """Return validated scalar field metadata from FIELDS."""
    if scalar_field not in FIELDS:
        raise KeyError(
            f"Unknown field {scalar_field!r}. Available fields are: {tuple(FIELDS)}"
        )
    field = FIELDS[scalar_field]
    if field["kind"] != "complex_scalar":
        raise ValueError(f"{scalar_field!r} is not a complex scalar field.")
    return field


def _validate_output_component(output_component, real_dimension=None):
    if len(output_component) != 4:
        raise ValueError("output_component must have the form (i, j, a, b).")

    if real_dimension is not None:
        _, _, a, b = output_component
        if not 1 <= a <= real_dimension:
            raise IndexError(f"Scalar index a={a} is outside 1,...,{real_dimension}.")
        if not 1 <= b <= real_dimension:
            raise IndexError(f"Scalar index b={b} is outside 1,...,{real_dimension}.")


def _gauge_generator_contribution(
    theta,
    gauge_coupling,
    i,
    j,
    a_index,
    b_index,
    real_dimension,
    coefficient,
):
    """Compute sum_cd g^2(theta_ac theta_bd + theta_bc theta_ad) C_ijcd."""
    result = sp.S.Zero
    for c_index in range(real_dimension):
        for d_index in range(real_dimension):
            c = c_index + 1
            d = d_index + 1
            symmetrized_generators = (
                theta[a_index, c_index] * theta[b_index, d_index]
                + theta[b_index, c_index] * theta[a_index, d_index]
            )
            result += (
                gauge_coupling**2
                * symmetrized_generators
                * coefficient[i, j, c, d]
            )
    return result


# The preceding Yukawa term in the general RGE vanishes for the lepton-sector
# indices used here because the required Yukawa interaction is absent.

@register_rge_contribution
def first_gauge_component(scalar_field, output_component, coefficient=C):
    r"""Return the first gauge contribution to beta(C_ijab).

    Implements

      - sum_alpha g_alpha^2 sum_cd
        (theta^alpha_ac theta^alpha_bd + theta^alpha_bc theta^alpha_ad) C_ijcd,

    where the overall minus sign is the -1/2 x 2 factor in the general RGE.
    SU(3) is absent for the colour-singlet scalars considered here.
    """
    field = _complex_scalar_field(scalar_field)
    _, su2_dimension, hypercharge = field["rep"]
    real_dimension = 2 * su2_dimension
    _validate_output_component(output_component, real_dimension=real_dimension)

    i, j, a, b = output_component
    a_index, b_index = a - 1, b - 1
    result = sp.S.Zero

    for theta in su2_real_scalar_generators(su2_dimension):
        result += _gauge_generator_contribution(
            theta, g2, i, j, a_index, b_index, real_dimension, coefficient
        )

    theta_y = u1_real_scalar_generator(su2_dimension, hypercharge)
    result += _gauge_generator_contribution(
        theta_y, g1, i, j, a_index, b_index, real_dimension, coefficient
    )

    return sp.simplify(-result)


@register_rge_contribution
def first_scalar_component(scalar_field, output_component, coefficient=C):
    r"""Return sum_cd lambda_abcd C_ijcd for the eta/H scalar system.

    The requested a,b indices are local to the selected doublet (1..4), while
    lambda_abcd and C_ijcd use the global eta=1..4, H=5..8 basis.
    """
    if scalar_field not in SCALAR_OFFSETS:
        raise ValueError("The first scalar component is implemented for 'eta' and 'H'.")

    _validate_output_component(output_component)
    i, j, a_local, b_local = output_component
    if not 1 <= a_local <= 4 or not 1 <= b_local <= 4:
        raise IndexError("Local scalar indices a and b must lie between 1 and 4.")

    offset = SCALAR_OFFSETS[scalar_field]
    a_global = a_local + offset
    b_global = b_local + offset
    result = sp.S.Zero

    # Both eta and H can appear in the internal scalar indices c,d.
    for c_global in range(1, 9):
        for d_global in range(1, 9):
            quartic = scalar_quartic_component(
                a_global, b_global, c_global, d_global
            )
            if quartic != 0:
                result += quartic * coefficient[i, j, c_global, d_global]

    return sp.simplify(result)


# -----------------------------------------------------------------------------
# C-basis to K-basis conversion
# -----------------------------------------------------------------------------

def c_to_k_substitutions(
    scalar_field=None,
    fermion_indices=(NU, NU),
    coefficient=C,
):
    """Return the neutrino-neutral-scalar C -> (K11,K22) substitutions.

    ``scalar_field`` is retained for compatibility with existing callers; the
    mapping intentionally covers both eta and H sectors because RGE terms can
    mix them through the internal scalar indices.
    """
    i, j = fermion_indices
    return {
        coefficient[i, j, a, b]: NEUTRAL_C_TO_K.get((a, b), sp.S.Zero)
        for a in range(1, 9)
        for b in range(1, 9)
    }


def convert_c_expression_to_k(
    expression,
    scalar_field,
    output_component,
    coefficient=C,
):
    """Replace supported neutrino-neutrino C components by K11 and K22."""
    i, j, _, _ = output_component
    if (i, j) != (NU, NU):
        raise NotImplementedError(
            "The current C-to-K conversion is implemented for "
            "neutrino-neutrino components only."
        )

    substitutions = c_to_k_substitutions(
        scalar_field=scalar_field,
        fermion_indices=(i, j),
        coefficient=coefficient,
    )
    return sp.simplify(sp.expand(expression.xreplace(substitutions)))


def kappa_dot_from_c_components(c_dot_33, c_dot_34, c_dot_44):
    r"""Convert neutral real-component beta functions to dot(K).

      dot(K) = 1/2 dot(C_33) - i dot(C_34) - 1/2 dot(C_44),

    using scalar symmetry dot(C_34)=dot(C_43).
    """
    return sp.simplify(HALF * c_dot_33 - sp.I * c_dot_34 - HALF * c_dot_44)


# -----------------------------------------------------------------------------
# RGE aggregation and display
# -----------------------------------------------------------------------------

def component_name(component):
    """Convert (i,j,a,b), e.g. (NU,NU,3,3), to a readable label."""
    i, j, a, b = component
    names = {NU: "nu", EL: "eL"}
    return f"{names.get(i, str(i))}{names.get(j, str(j))}{a}{b}"


def calculate_full_rge(
    scalar_field,
    output_component,
    coefficient=C,
    show_individual_terms=False,
):
    """Evaluate and sum every function registered in RGE_CONTRIBUTIONS."""
    contributions = []
    for contribution_function in RGE_CONTRIBUTIONS:
        contribution = sp.simplify(
            contribution_function(
                scalar_field=scalar_field,
                output_component=output_component,
                coefficient=coefficient,
            )
        )
        contributions.append(contribution)

        if show_individual_terms:
            print(f"\n{contribution_function.__name__}:")
            sp.pprint(contribution, use_unicode=True)

    return sp.simplify(sum(contributions, sp.S.Zero))


def print_full_rge(
    scalar_field,
    output_component,
    coefficient=C,
    show_individual_terms=False,
    show_latex=True,
    convert_to_k=True,
):
    """Calculate and print the accumulated RGE in either the C or K basis."""
    result_c = calculate_full_rge(
        scalar_field=scalar_field,
        output_component=output_component,
        coefficient=coefficient,
        show_individual_terms=show_individual_terms,
    )

    if convert_to_k:
        result = convert_c_expression_to_k(
            expression=result_c,
            scalar_field=scalar_field,
            output_component=output_component,
            coefficient=coefficient,
        )
        basis = "K basis"
    else:
        result = result_c
        basis = "C basis"

    print("\n" + "=" * 72)
    print(f"Full RGE result for {scalar_field}: {component_name(output_component)}")
    print(f"Included contributions: {len(RGE_CONTRIBUTIONS)}")
    print(f"Printed in: {basis}")
    print("-" * 72)
    sp.pprint(result, use_unicode=True)

    if show_latex:
        print("-" * 72)
        print("LaTeX:")
        print(sp.latex(result))

    print("=" * 72)
    return result


if __name__ == "__main__":
    FULL_RGE_RESULT = print_full_rge(
        scalar_field="eta",
        output_component=output_component,
        show_individual_terms=True,
        show_latex=True,
    )
