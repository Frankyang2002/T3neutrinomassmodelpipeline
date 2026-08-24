from dataclasses import dataclass
import sympy as sp

# Matrix convention: matrix[row, column]
# We want to convert from the general RGE in: Anomalous Dimension of a General Effective Gauge Theory II: Fermionic Sector
# To the RGE convention in: Running of Radiative Neutrino Masses: The Scotogenic Model – REVISITED

# ---------------
#   Definitions
# ---------------
HALF = sp.Rational(1, 2) #sp.Rational is there to have exact values instead of floating values like 0.3333...
g1, g2, g3 = sp.symbols("g1 g2 g3")
kappa = sp.Symbol("kappa")
C = sp.IndexedBase("C")
K11, K22 = sp.symbols("K11 K22")
lambda1, lambda2, lambda3, lambda4, lambda5 = sp.symbols(
    "lambda1 lambda2 lambda3 lambda4 lambda5"
)
LAMBDA = sp.IndexedBase("lambda")

# We define our fields, where sp.Symbol allows us to have symbols
FIELDS = {
    "L": {
        "id": 1,
        "symbol": sp.Symbol("L"),
        "kind": "fermion",
        "rep": (1, 2, -HALF),
        "components": {
            "nu": 1,
            "eL": 2,
        },
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

# Specified real component like nunu33
output_component=(NU,NU,3,3)




# --------------
#   Functions
# --------------
def su2_complex_generators(dimension):
    """
    Construct the usual complex Hermitian SU(2) generators
    T1, T2 and T3 for an irreducible representation.

    The representation dimension is dimension = 2*j + 1, where j is the weak isospin.

    Eg:
    dimension = 1  -> singlet, j = 0
    dimension = 2  -> doublet, j = 1/2
    dimension = 3  -> triplet, j = 1
    """

    # The SU(2) representation dimension must be a positive integer.
    if not isinstance(dimension, int) or dimension < 1:
        raise ValueError(
            "SU(2) representation dimension must be a positive integer."
        )

    # dim = 2*j + 1
    j = sp.Rational(dimension - 1, 2)

    # Get m values of: j, j-1, ..., -j.
    # EG: Doublet: m_values =[1/2, -1/2]
    # EG: Triplet: m_values = [1, 0, -1]
    m_values = [j - position for position in range(dimension)]

    # sp.zeros(n) creates an n x n symbolic zero matrix.
    T_plus = sp.zeros(dimension)

    # T_plus |j,m> = sqrt((j-m)(j+m+1)) |j,m+1>. for our raising operator
    for column, m in enumerate(m_values):
        # The raising operator maps m to m+1.
        target_m = m + 1

        # If m+1 not in m_values, T_plus|j,m+1>=0
        if target_m in m_values:

            # Find the row corresponding to the state |j,m+1>.
            row = m_values.index(target_m)

            # Calculate the sqrt((j-m)(j+m+1)) coefficient.
            coefficient = sp.sqrt((j - m) * (j + m + 1))

            # T_plus[m+1,m]
            T_plus[row, column] = coefficient

    # The lowering operator is the Hermitian conjugate of T_plus
    T_minus = T_plus.conjugate().T

    # We get T1,T2,T3 generators
    T1 = (T_plus + T_minus) / 2
    T2 = (T_plus - T_minus) / (2 * sp.I)
    T3 = sp.diag(*m_values)

    # Simplify every matrix before returning it.
    return (
        sp.simplify(T1),
        sp.simplify(T2),
        sp.simplify(T3),
    )


def real_scalar_generator(T):
    """
    Convert a complex n x n generator T into a 2n x 2n
    generator theta acting on real scalar components. 
    The separation of fermion and scalar generators is due to this n vs 2n basis

    For example: 
    A complex scalar multiplet is written as 
      eta_k = (R_k + i I_k)/sqrt(2).

    The real scalar basis is ordered as
      (R1, I1, R2, I2, ..., Rn, In).

    We use the convention to match the 
      delta(phi) = i * epsilon * theta * phi,

    where phi is the vector of real scalar fields.
    """
    # Number of complex scalar components.
    n = T.rows

    # Each complex field gives two real fields, so the real representation has dimension 2n.
    theta = sp.zeros(2 * n)

    # Loop over every complex-generator entry T[p,q].
    for p in range(n):
        for q in range(n):
            # Separate T[p,q] into its real and imaginary parts:
            #     T[p,q] = real_part + i*imag_part.
            entry = sp.expand_complex(T[p, q])
            real_part = sp.re(entry)
            imag_part = sp.im(entry)

            # Positions of R_p and I_p in the interleaved real basis, note python starts with 0 for indexing
            # p = 0 -> R1 at 0, I1 at 1
            # p = 1 -> R2 at 2, I2 at 3
            Rp = 2 * p
            Ip = 2 * p + 1

            # Positions of R_q and I_q.
            Rq = 2 * q
            Iq = 2 * q + 1

            # These four entries follow from expanding
            #     delta eta = i epsilon T eta where eta = (R + iI)/sqrt(2).
            # The result is written in the convention
            #     delta phi = i epsilon theta phi.
            # Our block matrix of theta, with phi = (R, I), is:
            #     theta = (iI, iR | -iR, iI)
            theta[Rp, Rq] += sp.I * imag_part  
            theta[Rp, Iq] += sp.I * real_part
            theta[Ip, Rq] += -sp.I * real_part
            theta[Ip, Iq] += sp.I * imag_part

    # Simplify the final real-scalar generator.
    return sp.simplify(theta)


def su2_real_scalar_generators(dimension):
    """
    Construct the three SU(2) generators acting directly
    on the real scalar components.

    This combines two steps:
        1. Generate the usual complex generators T1, T2, T3.
        2. Convert each one into the real-scalar basis.
    """

    # Generate T1, T2 and T3 in the complex representation.
    complex_generators = su2_complex_generators(dimension)

    # Convert each complex generator into a 2n x 2n real-scalar generator.
    real_generators = tuple(
        real_scalar_generator(T)
        for T in complex_generators
    )

    return real_generators


def u1_real_scalar_generator(dimension, hypercharge):
    """
    Construct the U(1)_Y generator acting on real scalar fields.

    For a complex scalar multiplet with hypercharge Y, the complex generator is
       T_Y = Y * identity.
    """
    T_Y = hypercharge * sp.eye(dimension)

    # Convert the complex U(1) generator into the real basis.
    theta_Y = real_scalar_generator(T_Y)

    return theta_Y




# -------
#   Scalar Components
# -------

def scalar_quartic_component(*indices):
    """
    Return the scalar quartic coupling lambda_{abcd}.

    The tensor is completely symmetric, so only the multiset of indices
    matters. The convention is

        eta real components: 1, 2, 3, 4
        H real components:   5, 6, 7, 8

    Any permutation of a listed component gives the same result.
    Components not present in the supplied tables return zero.
    """
    if len(indices) != 4:
        raise ValueError("Exactly four scalar indices are required.")

    if any(not isinstance(index, int) for index in indices):
        raise TypeError("Scalar indices must be Python integers.")

    if any(index < 1 or index > 8 for index in indices):
        raise IndexError("Scalar indices must lie between 1 and 8.")

    key = tuple(sorted(indices))

    # Pure eta sector: indices 1,...,4.
    if all(1 <= index <= 4 for index in key):
        multiplicities = sorted(
            [key.count(index) for index in set(key)],
            reverse=True,
        )

        if multiplicities == [4]:
            return 3 * lambda2

        if multiplicities == [2, 2]:
            return lambda2

        return sp.S.Zero

    # Pure Higgs sector: indices 5,...,8.
    if all(5 <= index <= 8 for index in key):
        multiplicities = sorted(
            [key.count(index) for index in set(key)],
            reverse=True,
        )

        if multiplicities == [4]:
            return 3 * lambda1

        if multiplicities == [2, 2]:
            return lambda1

        return sp.S.Zero

    # Mixed eta-Higgs sector. Keys are sorted, so all permutations
    # are handled automatically.
    mixed_components = {
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

    return sp.simplify(mixed_components.get(key, sp.S.Zero))


def lambda_substitutions(lambda_tensor=LAMBDA):
    """
    Build substitutions for every lambda[a,b,c,d] component with
    indices from 1 to 8.

    Because lambda_{abcd} is completely symmetric, each permutation
    receives the same converted value.
    """
    substitutions = {}

    for a in range(1, 9):
        for b in range(1, 9):
            for c in range(1, 9):
                for d in range(1, 9):
                    substitutions[lambda_tensor[a, b, c, d]] = (
                        scalar_quartic_component(a, b, c, d)
                    )

    return substitutions


def convert_lambda_expression(
    expression,
    lambda_tensor=LAMBDA,
):
    """
    Replace symbolic lambda[a,b,c,d] components by lambda1,...,lambda5.
    """
    return sp.simplify(
        sp.expand(
            expression.xreplace(
                lambda_substitutions(lambda_tensor=lambda_tensor)
            )
        )
    )


def print_lambda_component(*indices):
    """
    Print one scalar quartic component and its converted value.
    """
    value = scalar_quartic_component(*indices)
    label = "".join(str(index) for index in indices)

    print(f"lambda_{label} =")
    sp.pprint(value, use_unicode=True)

    return value

# -------
#   RGE
# -------
# NOTE: Things are are always true
# i and j are left handed Lepton sector indices, a b are scalar sector indices (either from eta or higgs)
RGE_CONTRIBUTIONS = []


def register_rge_contribution(function):
    """
    Register an RGE contribution function.

    Registered functions must accept the common arguments

        scalar_field, output_component, coefficient=C

    and return a SymPy expression.
    """
    RGE_CONTRIBUTIONS.append(function)
    return function


# ---------------------------
#   First Yukawa Component
# ---------------------------
# $y_{i j d}y_{k l d}^{*}\left[C_{\psi^{2}\phi^{2}}\right]_{k l a b}$
# Due to i and j being from the lepton sector, we lack the yukawa interaction for them and thus no contribution here and this is not affected by the representations

# ---------------------------
#   First Gauge Component
# ---------------------------

@register_rge_contribution
def first_gauge_component(
    scalar_field,
    output_component,
    coefficient=C,
):
    """
    Calculate the first gauge contribution
        -1/2 * 2 * (a,b permutations)
        * sum_alpha g_alpha^2
        * theta^alpha_ac theta^alpha_bd
        * C_ijcd

    which is implemented as
        -sum_alpha g_alpha^2 sum_cd [
             theta^alpha_ac theta^alpha_bd
           + theta^alpha_bc theta^alpha_ad
        ] C_ijcd.
    For our permutation

    Inputs
    ----------
    scalar_field : str
        Name of the scalar multiplet in FIELDS, for example "eta" or "H".

    output_component : tuple
        Requested component (i, j, a, b), using one-based indices.
        like for example: (NU, NU, 3, 3) -> C_{nu nu 3 3}.

    coefficient
        Symbolic coefficient tensor. By default this is the IndexedBase C.

    Outputs
    -------
    sympy expression
        The requested first gauge contribution.
    """

    # Check that the requested field exists in case mistype
    if scalar_field not in FIELDS:
        raise KeyError(
            f"Unknown field {scalar_field!r}. "
            f"Available fields are: {tuple(FIELDS)}"
        )

    field = FIELDS[scalar_field]

    # This function is specifically for scalar generators.
    if field["kind"] != "complex_scalar":
        raise ValueError(
            f"{scalar_field!r} is not a complex scalar field."
        )

    if len(output_component) != 4:
        raise ValueError(
            "output_component has form (i, j, a, b)."
        )

    i, j, a, b = output_component

    su3_dimension, su2_dimension, hypercharge = field["rep"]
    real_dimension = 2 * su2_dimension

    # Make sure our output components works for the representations we have, so we do not have a 7 for like a doublet
    if not 1 <= a <= real_dimension:
        raise IndexError(
            f"Scalar index a={a} is outside 1,...,{real_dimension}."
        )

    if not 1 <= b <= real_dimension:
        raise IndexError(
            f"Scalar index b={b} is outside 1,...,{real_dimension}."
        )

    # Convert from physics matrix to python matrix where python starts from 0
    a_python = a - 1
    b_python = b - 1

    # Initialise as 0
    result = sp.S.Zero

    # SU(3) contribution not included 

    # SU(2) contribution
    su2_generators = su2_real_scalar_generators(su2_dimension)
    # Sum over the three SU(2) generators.
    for theta in su2_generators:
        # Sum over internal real-scalar indices c and d.
        for c_python in range(real_dimension):
            for d_python in range(real_dimension):

                # One-based labels used when displaying C[i,j,c,d].
                c = c_python + 1
                d = d_python + 1

                permutation_term = (
                    theta[a_python, c_python]
                    * theta[b_python, d_python]
                    +
                    theta[b_python, c_python]
                    * theta[a_python, d_python]
                )

                result += (
                    g2**2
                    * permutation_term
                    * coefficient[i, j, c, d]
                )

    # U(1) contribution
    theta_Y = u1_real_scalar_generator(su2_dimension, hypercharge)

    for c_python in range(real_dimension):
        for d_python in range(real_dimension):
            c = c_python + 1
            d = d_python + 1

            permutation_term = (
                theta_Y[a_python, c_python]
                * theta_Y[b_python, d_python]
                +
                theta_Y[b_python, c_python]
                * theta_Y[a_python, d_python]
            )

            result += (
                g1**2
                * permutation_term
                * coefficient[i, j, c, d]
            )

    # The overall factor -1 comes from -1/2 * 2.
    return sp.simplify(-result)


# ---------------------------
#   First Scalar Component
# ---------------------------

@register_rge_contribution
def first_scalar_component(
    scalar_field,
    output_component,
    coefficient=C,
):
    """
    Calculate the scalar-quartic contribution

        sum_{c,d=1}^{8} lambda_{a b c d} C_{i j c d}.

    Global scalar-index convention:
        eta components: 1, 2, 3, 4
        H components:   5, 6, 7, 8

    The output indices a,b are supplied in the local real basis
    (1,2,3,4) of the selected scalar field. They are converted to
    global indices before looking up lambda_{abcd}.
    """
    if scalar_field not in ("eta", "H"):
        raise ValueError(
            "The first scalar component is implemented for 'eta' and 'H'."
        )

    if len(output_component) != 4:
        raise ValueError(
            "output_component must have the form (i, j, a, b)."
        )

    i, j, a_local, b_local = output_component

    if not 1 <= a_local <= 4 or not 1 <= b_local <= 4:
        raise IndexError(
            "Local scalar indices a and b must lie between 1 and 4."
        )

    # eta occupies global indices 1,...,4.
    # H occupies global indices 5,...,8.
    offset = 0 if scalar_field == "eta" else 4
    a_global = a_local + offset
    b_global = b_local + offset

    result = sp.S.Zero

    # Both scalar multiplets may run in the internal c,d indices.
    for c_global in range(1, 9):
        for d_global in range(1, 9):
            quartic = scalar_quartic_component(
                a_global,
                b_global,
                c_global,
                d_global,
            )

            if quartic != 0:
                result += (
                    quartic
                    * coefficient[i, j, c_global, d_global]
                )

    return sp.simplify(result)


# ---------------------------
#   C to K conversion
# ---------------------------


def c_to_k_substitutions(
    scalar_field=None,
    fermion_indices=(NU, NU),
    coefficient=C,
):
    """
    Build substitutions for both neutral scalar sectors.

    Global real-scalar ordering:
        eta: 1,2,3,4, with neutral components 3,4 -> K22
        H:   5,6,7,8, with neutral components 7,8 -> K11

    The nonzero neutrino-neutrino components are

        C_33 = K22/2,   C_34 = C_43 = i*K22/2,
        C_44 = -K22/2,

        C_77 = K11/2,   C_78 = C_87 = i*K11/2,
        C_88 = -K11/2.

    All other C_{nu nu a b} components in the 1,...,8 scalar basis
    are set to zero.
    """
    i, j = fermion_indices
    substitutions = {}

    for a in range(1, 9):
        for b in range(1, 9):
            component = coefficient[i, j, a, b]

            # eta neutral sector
            if (a, b) == (3, 3):
                substitutions[component] = K22 / 2
            elif (a, b) in ((3, 4), (4, 3)):
                substitutions[component] = sp.I * K22 / 2
            elif (a, b) == (4, 4):
                substitutions[component] = -K22 / 2

            # Higgs neutral sector
            elif (a, b) == (7, 7):
                substitutions[component] = K11 / 2
            elif (a, b) in ((7, 8), (8, 7)):
                substitutions[component] = sp.I * K11 / 2
            elif (a, b) == (8, 8):
                substitutions[component] = -K11 / 2

            else:
                substitutions[component] = sp.S.Zero

    return substitutions


def convert_c_expression_to_k(
    expression,
    scalar_field,
    output_component,
    coefficient=C,
):
    """
    Replace C components in an expression by K11 or K22 components.

    The current mapping applies to neutrino-neutrino coefficients and the
    neutral real components 3 and 4 of an SU(2) scalar doublet.
    """
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


def kappa_dot_from_c_components(
    c_dot_33,
    c_dot_34,
    c_dot_44,
):
    """
    Convert neutral real-component beta functions into dot(K):

        dot(K) = 1/2 dot(C_33) - i dot(C_34) - 1/2 dot(C_44),

    using scalar symmetry dot(C_34) = dot(C_43).
    """
    return sp.simplify(
        sp.Rational(1, 2) * c_dot_33
        - sp.I * c_dot_34
        - sp.Rational(1, 2) * c_dot_44
    )


def component_name(component):
    """
    Convert an index tuple such as (NU, NU, 3, 3) into a readable label.
    """

    i, j, a, b = component

    fermion_component_names = {
        NU: "nu",
        EL: "eL",
    }

    i_name = fermion_component_names.get(i, str(i))
    j_name = fermion_component_names.get(j, str(j))

    return f"{i_name}{j_name}{a}{b}"



def calculate_full_rge(
    scalar_field,
    output_component,
    coefficient=C,
    show_individual_terms=False,
):
    """
    Evaluate every registered RGE contribution and return their full sum.

    To include a new contribution automatically, define it with the same
    arguments and place @register_rge_contribution directly above it.
    """
    total_result = sp.S.Zero

    for contribution_function in RGE_CONTRIBUTIONS:
        contribution = sp.simplify(
            contribution_function(
                scalar_field=scalar_field,
                output_component=output_component,
                coefficient=coefficient,
            )
        )

        if show_individual_terms:
            print(f"\n{contribution_function.__name__}:")
            sp.pprint(contribution, use_unicode=True)

        total_result += contribution

    return sp.simplify(total_result)

def print_full_rge(
    scalar_field,
    output_component,
    coefficient=C,
    show_individual_terms=False,
    show_latex=True,
    convert_to_k=True,
):
    """
    Calculate and print the complete accumulated RGE expression.

    When convert_to_k=True, all supported C components are replaced by
    K11 for H or K22 for eta before printing.
    """
    full_result_c = calculate_full_rge(
        scalar_field=scalar_field,
        output_component=output_component,
        coefficient=coefficient,
        show_individual_terms=show_individual_terms,
    )

    if convert_to_k:
        full_result = convert_c_expression_to_k(
            expression=full_result_c,
            scalar_field=scalar_field,
            output_component=output_component,
            coefficient=coefficient,
        )
        expression_type = "K basis"
    else:
        full_result = full_result_c
        expression_type = "C basis"

    label = component_name(output_component)

    print("\n" + "=" * 72)
    print(f"Full RGE result for {scalar_field}: {label}")
    print(f"Included contributions: {len(RGE_CONTRIBUTIONS)}")
    print(f"Printed in: {expression_type}")
    print("-" * 72)
    sp.pprint(full_result, use_unicode=True)

    if show_latex:
        print("-" * 72)
        print("LaTeX:")
        print(sp.latex(full_result))

    print("=" * 72)
    return full_result


# This block runs only when this file is executed directly.
# It does not run when the file is imported into another Python file.
if __name__ == "__main__":
    FULL_RGE_RESULT = print_full_rge(
        scalar_field="eta",
        output_component=output_component,
        show_individual_terms=True,
        show_latex=True,
    )