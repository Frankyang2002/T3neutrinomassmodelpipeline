from __future__ import annotations

"""Shared human-facing notation for generated T3 reports.

This module owns presentation-only symbol conventions. It does not define
matching normalizations, RGE coefficients, or model physics.
"""

import re


PAPER_SYMBOL_LATEX: dict[str, str] = {
    "gY": r"g_1",
    "g2": r"g_2",
    "g3": r"g_3",
    "yu": r"Y_u",
    "yd": r"Y_d",
    "ye": r"Y_e",
    "y1": r"y_1",
    "y2": r"y_2",
    "h": r"h",
    "MF": r"M_F",
    "mSSq": r"m_S^2",
    "mS1Sq": r"m_1^2",
    "mS2Sq": r"m_2^2",
    "lambdaH": r"\lambda_1",
    "lambdaS": r"\lambda_2",
    "lambda3": r"\lambda_3",
    "lambda4": r"\lambda_4",
    "lambda5": r"\lambda_5",
    "lambdaS1": r"\lambda_{S_1}^{(1)}",
    "lambdaS2": r"\lambda_{S_2}^{(1)}",
    "lambdaH1": r"\lambda_{HS_1}^{(1)}",
    "lambdaH2": r"\lambda_{HS_2}^{(1)}",
    "lambdaH1Adj": r"\lambda_{HS_1}^{(A)}",
    "lambdaH2Adj": r"\lambda_{HS_2}^{(A)}",
    "lambdaT3": r"\lambda_5",
    "lambda12": r"\lambda_{12}^{(1)}",
    "lambda12Adj": r"\lambda_{12}^{(A)}",
    "lambda12Cross": r"\lambda_{12}^{(\times)}",
    "lambdaS1Adj": r"\lambda_{S_1}^{(A)}",
    "lambdaS2Adj": r"\lambda_{S_2}^{(A)}",
    "lambdaHHdagS1barS1bar": r"\lambda_{7}",
    "lambdaHHdagS2S2": r"\lambda_{8}",
    "lambdaS1bar2S2bar2": r"\lambda_{9}",
    "lambdaS1barS2S2bar2": r"\lambda_{10}",
    "lambdaS1S1bar2S2bar": r"\lambda_{11}",
    "lambdaHHdagS1barS2barCross": r"\lambda_{12}",
}


# RGBeta's raw special-cross symbol historically renders with its explicit
# field-content name before the later paper-notation cleanup. Keep that
# presentation behavior explicit rather than hiding it in RGEComparison.
RGBETA_SYMBOL_OVERRIDES: dict[str, str] = {
    "lambdaHHdagS1barS2barCross": (
        r"\lambda_{H H^\dagger S_1^\dagger S_2^\dagger}^{\rm Cross}"
    ),
}


RENDERED_PAPER_REPLACEMENTS: dict[str, str] = {
    r"\lambda_{H H^\dagger S_1^\dagger S_2^\dagger}^{\rm Cross}": r"\lambda_{12}",
    r"\lambda_{S_1 S_1^{\dagger 2}S_2^\dagger}": r"\lambda_{11}",
    r"\lambda_{S_1^\dagger S_2 S_2^{\dagger 2}}": r"\lambda_{10}",
    r"\lambda_{S_1^{\dagger 2}S_2^{\dagger 2}}": r"\lambda_9",
    r"\lambda_{H^\dagger H^\dagger S_2 S_2}": r"\lambda_8",
    r"\lambda_{H^\dagger H^\dagger S_1^\dagger S_1^\dagger}": r"\lambda_7",
    r"\lambda_{12}^{\rm Cross}": r"\lambda_{12}^{(\times)}",
    r"\lambda_{12}^{\rm Adj}": r"\lambda_{12}^{(A)}",
    r"\lambda_{H2}^{\rm Adj}": r"\lambda_{HS_2}^{(A)}",
    r"\lambda_{H1}^{\rm Adj}": r"\lambda_{HS_1}^{(A)}",
    r"\lambda_{S_2}^{\rm Adj}": r"\lambda_{S_2}^{(A)}",
    r"\lambda_{S_1}^{\rm Adj}": r"\lambda_{S_1}^{(A)}",
    r"\lambda_{T3}": r"\lambda_5",
    r"\lambda_{12}": r"\lambda_{12}^{(1)}",
    r"\lambda_{H2}": r"\lambda_{HS_2}^{(1)}",
    r"\lambda_{H1}": r"\lambda_{HS_1}^{(1)}",
    r"\lambda_{S_2}": r"\lambda_{S_2}^{(1)}",
    r"\lambda_{S_1}": r"\lambda_{S_1}^{(1)}",
    r"\lambda_H": r"\lambda_1",
    r"g_Y": r"g_1",
    r"m_{S_1}^2": r"m_1^2",
    r"m_{S_2}^2": r"m_2^2",
    r"M_{S_1}": r"m_1",
    r"M_{S_2}": r"m_2",
}


def paper_symbol_latex(name: str) -> str:
    """Return the canonical human-facing symbol for one pipeline quantity."""
    return PAPER_SYMBOL_LATEX.get(str(name), str(name))


def rgbeta_symbol_latex(name: str) -> str:
    """Return the report symbol used while normalising raw RGBeta LaTeX."""
    key = str(name)
    return RGBETA_SYMBOL_OVERRIDES.get(key, paper_symbol_latex(key))


BETA_LATEX_OVERRIDES: dict[str, str] = {
    # Preserve exact historical GroupFactorReports spelling.
    "h": r"\beta_h",
    "mSSq": r"\beta_{m_S^{2}}",
    "mS1Sq": r"\beta_{m_1^{2}}",
    "mS2Sq": r"\beta_{m_2^{2}}",
}


def beta_symbol_latex(name: str) -> str:
    """Return the canonical human-facing beta-function label."""
    key = str(name)

    if key in BETA_LATEX_OVERRIDES:
        return BETA_LATEX_OVERRIDES[key]

    if key in PAPER_SYMBOL_LATEX:
        return rf"\beta_{{{PAPER_SYMBOL_LATEX[key]}}}"

    return rf"\beta_{{\mathrm{{{key}}}}}"


def paper_notation_key_lines(*, heading: str = "Notation key") -> list[str]:
    """Return the common symbol/index key used in human-facing reports."""
    return [
        rf"\section*{{{heading}}}",
        r"\begin{longtable}{@{}p{0.23\linewidth}p{0.69\linewidth}@{}}",
        r"\toprule",
        r"symbol & interaction / definition \\",
        r"\midrule",
        r"$g_1$ & $D_\mu\supset i g_1YB_\mu$ \\",
        r"$g_2$ & $D_\mu\supset i g_2T^AW_\mu^A$ \\",
        r"$g_3$ & $D_\mu\supset i g_3t^AG_\mu^A$ \\",
        r"$Y_u$ & $-\bar QY_u\widetilde Hu+\mathrm{h.c.}$ \\",
        r"$Y_d$ & $-\bar QY_dHd+\mathrm{h.c.}$ \\",
        r"$Y_e$ & $-\bar LY_eHe+\mathrm{h.c.}$ \\",
        r"$y_1$ & $LFS_1+\mathrm{h.c.}$ \\",
        r"$y_2$ & $LFS_2+\mathrm{h.c.}$ \\",
        r"$M_F$ & $-\bar FM_FF$ \\",
        r"$T$ & $\operatorname{Tr}(Y_e^\dagger Y_e+3Y_u^\dagger Y_u+3Y_d^\dagger Y_d)$ \\",
        r"$T_\nu^{(1)}$ & $\operatorname{Tr}(y_1^\dagger y_1)$ \\",
        r"$T_\nu^{(2)}$ & $\operatorname{Tr}(y_2^\dagger y_2)$ \\",
        r"$T_\nu$ & $T_\nu^{(1)}+T_\nu^{(2)}$ \\",
        r"$m_1^2$ & $m_1^2S_1^\dagger S_1$ \\",
        r"$m_2^2$ & $m_2^2S_2^\dagger S_2$ \\",
        r"$\lambda_1$ & $\frac12\lambda_1(H^\dagger H)^2$ \\",
        r"$\lambda_{S_1}^{(1)}$ & $\frac12\lambda_{S_1}^{(1)}(S_1^\dagger S_1)^2$ \\",
        r"$\lambda_{S_2}^{(1)}$ & $\frac12\lambda_{S_2}^{(1)}(S_2^\dagger S_2)^2$ \\",
        r"$\lambda_{HS_1}^{(1)}$ & $\lambda_{HS_1}^{(1)}(H^\dagger H)(S_1^\dagger S_1)$ \\",
        r"$\lambda_{HS_2}^{(1)}$ & $\lambda_{HS_2}^{(1)}(H^\dagger H)(S_2^\dagger S_2)$ \\",
        r"$\lambda_{HS_1}^{(A)}$ & $\lambda_{HS_1}^{(A)}(H^\dagger T^AH)(S_1^\dagger T^AS_1)$ \\",
        r"$\lambda_{HS_2}^{(A)}$ & $\lambda_{HS_2}^{(A)}(H^\dagger T^AH)(S_2^\dagger T^AS_2)$ \\",
        r"$\lambda_5$ & $\lambda_5HHS_1S_2^\dagger+\mathrm{h.c.}$ \\",
        r"$\lambda_{12}^{(1)}$ & $\lambda_{12}^{(1)}(S_1^\dagger S_1)(S_2^\dagger S_2)$ \\",
        r"$\lambda_{12}^{(A)}$ & $\lambda_{12}^{(A)}(S_1^\dagger T^AS_1)(S_2^\dagger T^AS_2)$ \\",
        r"$\lambda_{12}^{(\times)}$ & crossed independent $S_1$--$S_2$ contraction \\",
        r"$\lambda_{S_1}^{(A)}$ & additional independent $S_1$ self-contraction \\",
        r"$\lambda_{S_2}^{(A)}$ & additional independent $S_2$ self-contraction \\",
        r"$\alpha$ & $Y_{S_1}=\alpha/2,\ Y_F=(\alpha+1)/2,\ Y_{S_2}=(\alpha+2)/2$ \\",
        r"$C_2(d)$ & $(d^2-1)/4$ \\",
        r"$T(d)$ & $d(d^2-1)/12$ \\",
        r"$C_{LLS_1S_2}^{ij;AB}$ & $L_iL_jS_1^AS_2^B$ \\",
        r"$C_{ijab}$ & $L_iL_j\phi_a\phi_b$ \\",
        r"$C_5^{ij}$ & $(L_i^TCL_j)HH$ \\",
        r"$R_W$ & $16\pi^2\beta_{C_5}\supset R_W\lambda_5C_{LLS_1S_2}$ \\",
        r"\bottomrule",
        r"\end{longtable}",
    ]


def apply_paper_notation(latex: str) -> str:
    """Convert rendered implementation symbols to common report notation."""
    if not latex:
        return latex

    text = str(latex)

    for raw, pretty in RENDERED_PAPER_REPLACEMENTS.items():
        text = text.replace(raw, pretty)

    for raw in sorted(PAPER_SYMBOL_LATEX, key=len, reverse=True):
        pretty = PAPER_SYMBOL_LATEX[raw]
        text = text.replace(r"\text{" + raw + "}", "{" + pretty + "}")
        text = re.sub(
            rf"(?<![A-Za-z0-9]){re.escape(raw)}(?![A-Za-z0-9])",
            lambda _match, replacement=pretty: replacement,
            text,
        )

    return text
