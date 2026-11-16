"""Insert source-derived, representation-specific C12 factors into EXISTING HTML.

This module modifies the supplied HTML in place, preserving its other content.
No new group-factor HTML template, physics integration, or matching is run.

Example:
 python -m Reports.C12GroupFactorHTML \
   --html Reports/output/analytical/GroupFactors/GF_EFT_1_after_F.html \
   --data-dir output/single/T3_B_alpha_m1/data
"""
from __future__ import annotations

import argparse
from html import escape
import json
from pathlib import Path
import re

import sympy as sp
from sympy.printing.mathml import mathml

from RGE.group_factors.validation.C12RepresentationGroupFactors import (
    derive_from_saved_data,
)

BEGIN = "<!-- T3-C12-REPRESENTATION-GROUP-FACTORS:BEGIN -->"
END = "<!-- T3-C12-REPRESENTATION-GROUP-FACTORS:END -->"

LABELS = {
    "total": "Total projected C12 beta",
    "scalar_pair": "Scalar-pair gauge + quartic tensor action",
    "mixed_yukawa_gauge": "Mixed Yukawa / gauge tensor action",
    "scalar_anomalous_dimension": "Scalar anomalous dimensions",
    "fermion_anomalous_dimension": "Lepton anomalous dimensions",
    "yukawa_wavefunction": "Yukawa wavefunction tensor term",
    "crossed_yukawa": "Crossed Yukawa tensor term",
    "conjugate_coefficient_yukawa": "Conjugate-Wilson / Yukawa tensor term",
}


def _mathml(expression: str | None) -> str:
    if expression is None:
        return "<span>Not proportional on C12 seed support</span>"
    try:
        body = mathml(sp.sympify(expression), printer="presentation")
        return f'<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline">{body}</math>'
    except Exception:
        return f"<code>{escape(expression)}</code>"


def render_c12_section(result: dict, *, label: str = "") -> str:
    """Render one embeddable HTML section, not a replacement HTML document."""
    rep = result["representation"]
    descriptor = (
        f"dS1={rep['dS1']}, dS2={rep['dS2']}, "
        f"dF={rep['dF']}, alpha={rep['alpha']}"
    )
    details = []
    for key in ("total", "scalar_pair", "mixed_yukawa_gauge", "scalar_anomalous_dimension",
                "fermion_anomalous_dimension", "yukawa_wavefunction", "crossed_yukawa",
                "conjugate_coefficient_yukawa"):
        row = result["sectors"][key]
        coefficient = row["coefficient"]
        details.append(
            "<tr><th scope=\"row\">" + escape(LABELS[key]) + "</th>"
            "<td>" + _mathml(coefficient) + "</td>"
            "<td>" + ("yes" if row["collinear_on_c12_seed_support"] else "NO") + "</td></tr>"
        )
    validation = result["sectors"]["total"]["collinear_on_c12_seed_support"]
    paper = result.get("t3b_paper_basis")
    paper_row = ""
    if paper:
        paper_row = (
            '<h3>Paper quartic basis (T3-B only)</h3>'
            '<p><code>lambda12Inv1 = -lambda12/2 - 3*lambda12Adj/8</code>; '
            '<code>lambda12Inv2 = -sqrt(3)*lambda12/2 + sqrt(3)*lambda12Adj/8</code>.</p>'
            '<p>16π² β<sub>C12</sub> |<sub>diagonal, 1 generation</sub> = '
            + _mathml(paper["coefficient"]) + ' C<sub>12</sub>.</p>'
            '<p>Independent two-operator diagonal coefficient check: <strong>'
            + ("passed" if paper["matches_previous_two_operator_diagonal_coefficient"] else "FAILED")
            + '</strong>.</p>'
        )
    source = result.get("sources", {})
    provenance = "".join(
        f"<li>{escape(name)}: <code>{escape(str(path))}</code></li>"
        for name, path in source.items()
    )
    generated = result["generated_mixed_components_outside_c12_seed_support"]
    generated_note = (
        f"{len(generated)} mixed-yukawa beta components lie outside the C12 matching seed support. "
        "These include operator mixing; this seed-support projection is not a full operator-closure proof."
    )
    note = ("Verified proportionality for the total beta on the C12 matching-tensor support."
            if validation else
            "No unique diagonal C12 coefficient exists on this matching-tensor support; inspect noncollinear residuals.")
    return "\n".join([
        BEGIN,
        '<section id="c12-representation-group-factors" style="margin-top:2rem;padding-top:1rem;border-top:2px solid #888">',
        '<h2>Representation-specific C<sub>12</sub> group factors</h2>',
        f'<p><strong>{escape(label or "Matched model")}</strong>: {escape(descriptor)}. '
        f'Scalar-only F-first EFT, with all heavy fermions removed together.</p>',
        '<p>Derived from saved <code>C_ijab</code> matching tensors and the complete one-loop '
        '<code>psi² phi²</code> master RGE, isolating the mixed '
        '<code>conjugate(y1)*conjugate(y2)</code> matching channel. '
        'All listed coefficients multiply <code>C12</code> in <code>16*pi²*beta_C12</code>.</p>',
        f'<p>Independent matched tensor components checked: '
        f'<strong>{result["seed_support_component_count"]}</strong>. '
        f'{escape(note)}</p>',
        '<div style="overflow-x:auto"><table><thead><tr><th>Master RGE contribution</th>'
        '<th>Group-factor coefficient (RGBeta tensor basis)</th><th>Collinear?</th>'
        '</tr></thead><tbody>',
        *details,
        '</tbody></table></div>',
        paper_row,
        '<p>' + escape(generated_note) + '</p>',
        '<p><strong>Scope:</strong> This table concerns the diagonal action on the matched C12 channel. '
        'It does not derive the reverse C5 → C12 factor or establish full three-generation closure. '
        'The standard fixed-order coefficient and the exploratory resummed numerical result remain separate.</p>',
        '<details><summary>Exact extraction metadata / residuals</summary><pre>'
        + escape(json.dumps(result, indent=2)) + '</pre></details>',
        '<details><summary>Saved physics inputs</summary><ul>' + provenance + '</ul></details>',
        '</section>', END,
    ])


def replace_existing_section(html: str, section: str) -> str:
    """Idempotently update one section while preserving the rest of the HTML."""
    if (BEGIN in html) != (END in html):
        raise ValueError("Existing HTML contains only one C12 report marker; refusing to corrupt it")
    if BEGIN in html:
        start = html.index(BEGIN)
        stop = html.index(END, start) + len(END)
        return html[:start] + section + html[stop:]
    # Keep original styling, tables, figures, scripts and report navigation.
    # Append inside the existing main element where possible.
    candidate = re.search(r"</main\s*>", html, flags=re.IGNORECASE)
    if candidate is None:
        candidate = re.search(r"</body\s*>", html, flags=re.IGNORECASE)
    if candidate is None:
        raise ValueError("Target is not an HTML report with a closing </main> or </body> tag")
    return html[:candidate.start()] + "\n" + section + "\n" + html[candidate.start():]


def insert_c12_into_existing_html(*, html_path: Path, data_dir: Path, label: str = "") -> dict:
    html_path = Path(html_path)
    if not html_path.is_file():
        raise FileNotFoundError(f"Existing report HTML is required: {html_path}")
    result = derive_from_saved_data(Path(data_dir))
    section = render_c12_section(result, label=label)
    old = html_path.read_text(encoding="utf-8-sig")
    new = replace_existing_section(old, section)
    if new != old:
        # Atomic replacement, no extra HTML file and no touching the PDF/TeX.
        tmp = html_path.with_name(html_path.name + ".c12tmp")
        tmp.write_text(new, encoding="utf-8")
        tmp.replace(html_path)
    return {
        "status": "RepresentationSpecificC12FactorsInserted",
        "html": str(html_path),
        "representation": result["representation"],
        "seed_component_count": result["seed_support_component_count"],
        "total_beta_collinear": result["sectors"]["total"]["collinear_on_c12_seed_support"],
        "total_beta_group_factor": result["sectors"]["total"]["coefficient"],
        "paper_basis": result.get("t3b_paper_basis"),
    }



def insert_many_c12_into_existing_html(*, html_path: Path, model_data: list[tuple[str, Path]]) -> list[dict]:
    """Update one existing combined-model HTML with independently derived factors.

    Never mix coefficients from different representations. Missing or invalid
    source physics is an error rather than a reason to substitute T3-B values.
    """
    if not model_data:
        return []
    html_path = Path(html_path)
    if not html_path.is_file():
        raise FileNotFoundError(f"Existing report HTML is required: {html_path}")
    blocks, summaries = [], []
    seen = set()
    for name, directory in model_data:
        directory = Path(directory)
        if directory.resolve() in seen:
            continue
        seen.add(directory.resolve())
        data = derive_from_saved_data(directory)
        subsection = render_c12_section(data, label=name)
        subsection = subsection.replace(BEGIN, '').replace(END, '')
        # The combined report has many sections; IDs must not be duplicated.
        subsection = subsection.replace('id="c12-representation-group-factors"',
                                        'class="c12-representation-group-factors"')
        blocks.append(subsection)
        summaries.append({
            'name': name, 'data_dir': str(directory),
            'representation': data['representation'],
            'total_beta_collinear': data['sectors']['total']['collinear_on_c12_seed_support'],
            'total_beta_group_factor': data['sectors']['total']['coefficient'],
        })
    content = BEGIN + '\n' + '\n'.join(blocks) + '\n' + END
    old = html_path.read_text(encoding='utf-8-sig')
    new = replace_existing_section(old, content)
    if old != new:
        temporary = html_path.with_name(html_path.name + '.c12tmp')
        temporary.write_text(new, encoding='utf-8')
        temporary.replace(html_path)
    return summaries

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", type=Path, required=True, help="EXISTING report HTML to update in place")
    parser.add_argument("--data-dir", type=Path, required=True, help="Saved model data containing the Wilson CG seed and one-loop beta")
    parser.add_argument("--label", default="", help="Display-only model name")
    args = parser.parse_args(argv)
    result = insert_c12_into_existing_html(html_path=args.html, data_dir=args.data_dir, label=args.label)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
