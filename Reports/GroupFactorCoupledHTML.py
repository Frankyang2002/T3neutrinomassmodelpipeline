"""Offline HTML companions for stage-aware GroupFactors reports.

Presentation-only: no RGE is calculated here. For the ordinary T3-B
(dS1,dS2,dF,alpha)=(2,2,1,-1) F-first EFT, display the retained
C12 <-> ordered-C5 one-loop matrix equation that has been checked against
the project's tensor RGE and numerical benchmark. Other representations
are deliberately not assigned the T3-B mixing coefficients.

CLI to rebuild saved GroupFactors HTML without rerunning Matchete/RGBeta:
    python -m Reports.GroupFactorCoupledHTML \
        --report-root Reports/output/analytical --t3b
"""
from __future__ import annotations

import argparse
from html import escape
from pathlib import Path
from typing import Iterable, Sequence


T3B_REPRESENTATION = (2, 2, 1, -1, False)


def _is_t3b(record: object) -> bool:
    actual = (
        int(getattr(record, "d_s1")),
        int(getattr(record, "d_s2")),
        int(getattr(record, "d_f")),
        int(getattr(record, "alpha")),
        bool(getattr(record, "shared_scalar", False)),
    )
    return actual == T3B_REPRESENTATION


def _is_f_first_scalar_eft(record: object, stage_label: str) -> bool:
    """Use actual serialized EFT content; a file name alone is insufficient."""
    summary = getattr(record, "summary", {})
    for stage in summary.get("EFTStages", []) or []:
        if str(stage.get("Label", "")).strip() != stage_label:
            continue
        return (
            set(stage.get("IntegratedFields", []) or []) == {"F"}
            and set(stage.get("ActiveHeavyFields", []) or []) == {"S1", "S2"}
        )
    return False


def _is_final_eft(record: object, stage_label: str) -> bool:
    for stage in getattr(record, "summary", {}).get("EFTStages", []) or []:
        if str(stage.get("Label", "")).strip() == stage_label:
            return not (stage.get("ActiveHeavyFields", []) or [])
    return False


def _head(title: str) -> str:
    return """<!doctype html>
<html lang="en">
<head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>""" + escape(title) + """</title>
<style>
:root {color-scheme:light dark; font-family: system-ui, 'Segoe UI', Arial, sans-serif;}
body {margin:0; background:#f6f7fa; color:#1f2937; line-height:1.6;}
main {max-width:1050px; margin:38px auto; padding:32px 42px; background:white;
      border:1px solid #dbe1e8; border-radius:12px;}
h1 {margin:0 0 8px; font-size:1.8rem;} h2 {margin-top:31px; border-bottom:1px solid #e5e7eb;
      padding-bottom:7px; font-size:1.3rem;} h3 {font-size:1rem;}
.subtitle,.muted {color:#526174;} .equation {background:#f7f8fb; border-left:3px solid #506990;
      padding:16px 19px; margin:15px 0; overflow-x:auto; font-family:'Cambria Math', Georgia,serif;
      font-size:1.11rem; white-space:nowrap;}
table {border-collapse:collapse; width:100%; font-size:.96rem;}
th,td {border-bottom:1px solid #dfe3e9; padding:10px 12px; text-align:left; vertical-align:top;}
th {background:#f2f4f8;} code,pre {font-family:Consolas,monospace;}
pre {white-space:pre-wrap; background:#f3f4f6; border-radius:6px; padding:14px; overflow-wrap:anywhere;}
a {color:#1959a7;} .links a {margin-right:18px;}
.label {font-size:.76rem; letter-spacing:.06em; text-transform:uppercase; font-weight:700; color:#526174;}
@media (prefers-color-scheme:dark) {body{background:#0d1320;color:#d9e3f3;}main{background:#151e2d;
 border-color:#38465b;}th{background:#243247;}td,th{border-color:#39475b;}.equation{background:#1d293a;}
 pre{background:#243247;}.subtitle,.muted,.label{color:#b4c1d2;}a{color:#91c2ff;}}
@media print {body{background:white;}main{border:0;margin:0;max-width:none;} }
</style></head><body><main>
"""


def _summary(title: str, pdf_name: str, tex_name: str) -> str:
    return (
        f'<h1>{escape(title)}</h1>'
        '<p class="subtitle">Stage-aware group factors and Wilson-operator mixing.</p>'
        '<p class="links">'
        f'<a href="{escape(pdf_name, quote=True)}">Full group-factors PDF</a>'
        f'<a href="{escape(tex_name, quote=True)}">LaTeX source</a>'
        '</p>'
        '<p class="muted">The linked PDF contains the complete representation and '
        'renormalisable group-factor tables. The equations below summarise the '
        'Wilson mixing sector where verified.</p>'
    )


def _t3b_coupled_section(names: Sequence[str]) -> str:
    shown = ', '.join(escape(x) for x in names)
    return """
<section id="c12-c5-rge">
<h2>Coupled dimension-five RGEs: T3-B (F-first EFT)</h2>
<p><span class="label">Scope: T3-B, &alpha; = &minus;1</span></p>
<p>Applicable records: """ + (shown or "T3-B (explicit CLI selection)") + """.</p>
<p>All three heavy fermions are integrated out together at &mu;<sub>F</sub>.
The remaining theory is SM + S<sub>1</sub> + S<sub>2</sub>.
The flavour-symmetric matrices C<sub>12</sub> and C<sub>5</sub>
are evolved simultaneously in the <strong>retained two-operator sector</strong>.</p>
<h3>One-loop anomalous-dimension system</h3>
<div class="equation">
16&pi;<sup>2</sup> dC<sub>12</sub>/dln&mu; =
(&minus;3g<sub>2</sub><sup>2</sup> + 2&lambda;<sub>12</sub> &minus; &frac12;&lambda;<sub>12Adj</sub>) C<sub>12</sub>
+ &frac12;(EC<sub>12</sub> + C<sub>12</sub>E<sup>T</sup>)
+ (8&radic;3/3)&lambda;<sub>T3</sub><sup>*</sup> C<sub>5</sub>
</div>
<div class="equation">
16&pi;<sup>2</sup> dC<sub>5</sub>/dln&mu; =
(&minus;3g<sub>2</sub><sup>2</sup> + 2&lambda;<sub>H</sub> + 2T) C<sub>5</sub>
&minus; (3/2)(EC<sub>5</sub> + C<sub>5</sub>E<sup>T</sup>)
+ (4&radic;3/3)&lambda;<sub>T3</sub> C<sub>12</sub>
</div>
<p>Here E = Y<sub>e</sub>Y<sub>e</sub><sup>&dagger;</sup> and
T = Tr(E + 3Y<sub>u</sub>Y<sub>u</sub><sup>&dagger;</sup>
+ 3Y<sub>d</sub>Y<sub>d</sub><sup>&dagger;</sup>).
In this intermediate-EFT system, C<sub>5</sub> denotes the
<strong>ordered Wilson-tensor convention</strong>; the physical Majorana
coefficient uses the matching layer's factor-of-two conversion.</p>
<h3>Group factors in 16&pi;<sup>2</sup>&beta;</h3>
<table>
<thead><tr><th>Source structure</th><th>&beta;<sub>C12</sub></th><th>&beta;<sub>C5</sub></th></tr></thead>
<tbody>
<tr><td>SU(2)<sub>L</sub> gauge</td><td>&minus;3g<sub>2</sub><sup>2</sup>C<sub>12</sub></td>
<td>&minus;3g<sub>2</sub><sup>2</sup>C<sub>5</sub></td></tr>
<tr><td>Scalar self/mixed quartics</td><td>(2&lambda;<sub>12</sub> &minus; &frac12;&lambda;<sub>12Adj</sub>)C<sub>12</sub></td>
<td>2&lambda;<sub>H</sub>C<sub>5</sub></td></tr>
<tr><td>Yukawa traces</td><td>&mdash;</td><td>2TC<sub>5</sub></td></tr>
<tr><td>Lepton flavour (left + transpose right)</td>
<td>+&frac12;(EC<sub>12</sub> + C<sub>12</sub>E<sup>T</sup>)</td>
<td>&minus;(3/2)(EC<sub>5</sub> + C<sub>5</sub>E<sup>T</sup>)</td></tr>
<tr><td><strong>C<sub>5</sub> &rarr; C<sub>12</sub></strong></td>
<td>+(8&radic;3/3)&lambda;<sub>T3</sub><sup>*</sup>C<sub>5</sub></td><td>&mdash;</td></tr>
<tr><td><strong>C<sub>12</sub> &rarr; C<sub>5</sub></strong></td><td>&mdash;</td>
<td>+(4&radic;3/3)&lambda;<sub>T3</sub>C<sub>12</sub></td></tr>
</tbody></table>
<details><summary>Copyable LaTeX equations</summary>
<pre>16\\pi^2 \\frac{dC_{12}}{d\\ln\\mu}=
(-3g_2^2+2\\lambda_{12}-\\tfrac12\\lambda_{12Adj})C_{12}
+\\tfrac12(EC_{12}+C_{12}E^T)
+\\frac{8\\sqrt3}{3}\\lambda_{T3}^* C_5

16\\pi^2 \\frac{dC_5}{d\\ln\\mu}=
(-3g_2^2+2\\lambda_H+2T)C_5
-\\tfrac32(EC_5+C_5E^T)
+\\frac{4\\sqrt3}{3}\\lambda_{T3} C_{12}</pre></details>
<h3>Boundary, power counting and validation</h3>
<p>The common matching point is &mu;<sub>F</sub>. Individual heavy masses M<sub>r</sub>
appear in matched coefficients and ln(&mu;<sub>F</sub>/M<sub>r</sub>),
while both Wilson matrices evolve down to the scalar threshold &mu;<sub>S</sub>.
Fixed-order production contains the O((16&pi;<sup>2</sup>)<sup>&minus;1</sup>)
forward contribution. The coupled solver is an <strong>opt-in leading-log
resummation</strong>; its first corrections and reverse feedback are higher-order
Wilson-running effects, not a complete NLO matching result.</p>
<p><strong>Evidence:</strong> scalar basis and operator mixing from
<code>RGE/general/WilsonTensorRGE.py</code>; charged-lepton flavour action from
<code>RGE/general/AnomalousDimensions.py</code>; diagonal Weinberg matrix action
from <code>RGE/running/weinberg/WeinbergRGE.py</code>.
The numerical benchmark's retained two-operator sector has passed one-generation
checks, but complete three-generation operator closure and higher-order
matching are not established.</p>
</section>
"""


def _final_c5_section() -> str:
    return """
<section id="c5-rge"><h2>Final SMEFT: Weinberg operator running</h2>
<p>Below the scalar threshold S<sub>1</sub>, S<sub>2</sub> and F are absent,
so C<sub>12</sub> &harr; C<sub>5</sub> mixing no longer occurs.
Only C<sub>5</sub> runs:</p>
<div class="equation">16&pi;<sup>2</sup> dC<sub>5</sub>/dln&mu; =
(2&lambda;<sub>H</sub> &minus; 3g<sub>2</sub><sup>2</sup> + 2T)C<sub>5</sub>
&minus; (3/2)(EC<sub>5</sub>+C<sub>5</sub>E<sup>T</sup>)</div>
<p>E = Y<sub>e</sub>Y<sub>e</sub><sup>&dagger;</sup>,
T = Tr(E + 3Y<sub>u</sub>Y<sub>u</sub><sup>&dagger;</sup>
+ 3Y<sub>d</sub>Y<sub>d</sub><sup>&dagger;</sup>).
Source: <code>RGE/running/weinberg/WeinbergRGE.py</code>.</p></section>
"""


def render_group_factor_html(
    *, stage_label: str, pdf_name: str, tex_name: str,
    t3b_names: Sequence[str] = (), final_eft: bool = False,
) -> str:
    """Build a complete offline-readable HTML report from verified scope."""
    title = f"Group factors — {stage_label}"
    sections = [_head(title), _summary(title, pdf_name, tex_name)]
    if stage_label == "UV":
        sections.append('<h2>UV field content</h2><p>Heavy fermions are still active. '
                        'The scalar-only EFT Wilson RGE does not apply above the F threshold.</p>')
    elif t3b_names:
        sections.append(_t3b_coupled_section(t3b_names))
    elif final_eft:
        sections.append(_final_c5_section())
    else:
        sections.append('<h2>Wilson RGE scope</h2><p>A T3-B-specific '
                        'coupled C<sub>12</sub>&harr;C<sub>5</sub> equation is not displayed '
                        'because the available records do not establish the required '
                        'F-first T3-B intermediate EFT. See the linked full group-factors report.</p>')
    sections.append('</main></body></html>\n')
    return '\n'.join(sections)


def write_group_factor_html_companions(
    records: Iterable[object], tex_paths: Iterable[Path], *,
    report_root: Path | None = None,
    explicit_t3b: bool = False,
    representation_data_dir: Path | None = None,
    preserve_existing_html: bool = False,
) -> list[Path]:
    """Write one companion HTML per actual stage TeX (no physics rerun)."""
    known = list(records)
    outputs: list[Path] = []
    for tex_file in tex_paths:
        tex = Path(tex_file)
        if not tex.is_file():
            raise FileNotFoundError(f"GroupFactors TeX does not exist: {tex}")
        if not tex.stem.startswith('GF_'):
            raise ValueError(f"Unexpected GroupFactors filename: {tex.name}")
        label = tex.stem[len('GF_'):]
        t3b = [str(getattr(record, 'name', 'T3-B')) for record in known
                if _is_t3b(record) and _is_f_first_scalar_eft(record, label)]
        # Explicit CLI mode operates on an already-generated F-first TeX;
        # it never guesses model identity from a filename alone.
        if explicit_t3b and label in ('EFT_1_after_F', 'EFT_1_after_F_with_C12') and not t3b:
            t3b = ['T3-B alpha=-1 (explicit CLI selection)']
        stage_records = [
            r for r in known
            if any(str(s.get('Label', '')).strip() == label
                   for s in getattr(r, 'summary', {}).get('EFTStages', []) or [])
        ]
        # Do not infer full decoupling from a filename: the actual EFT content
        # is given by the serialized stage records.
        final = bool(stage_records) and all(_is_final_eft(r, label)
                                            for r in stage_records)
        html_path = tex.with_suffix('.html')
        if not preserve_existing_html or not html_path.exists():
            html_path.write_text(
                render_group_factor_html(stage_label=label,
                                         pdf_name=tex.with_suffix('.pdf').name,
                                         tex_name=tex.name,
                                         t3b_names=t3b, final_eft=final),
                encoding='utf-8',
            )
        # Representation-specific coefficients are derived independently from
        # each model's saved tensor seed and master beta (never from a T3-B
        # universal constant). Insert into this same existing HTML report.
        from Reports.C12GroupFactorHTML import insert_many_c12_into_existing_html
        sources = []
        for r in known:
            if not _is_f_first_scalar_eft(r, label):
                continue
            if getattr(r, 'output_dir', None) is None:
                continue
            directory = Path(r.output_dir) / 'data'
            if all((directory / filename).is_file() for filename in (
                'eft1_wilson_rge.json', 'eft1_after_F_wilson_seed.json',
                'eft1_rgbeta_rge.json',
            )):
                sources.append((str(getattr(r, 'name', 'T3 model')), directory))
        if representation_data_dir is not None and t3b:
            # Explicit local T3-B run (no reconstructed RunRecord available).
            # Saved matching data may not exist in mock/unit or old runs;
            # the existing HTML remains valid and unchanged in that case.
            candidate = Path(representation_data_dir)
            if all((candidate / filename).is_file() for filename in (
                'eft1_wilson_rge.json', 'eft1_after_F_wilson_seed.json',
                'eft1_rgbeta_rge.json',
            )):
                sources = [(', '.join(t3b), candidate)]
        if sources:
            insert_many_c12_into_existing_html(html_path=html_path, model_data=sources)
        outputs.append(html_path)
    return outputs


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report-root', type=Path, default=Path('Reports/output/analytical'))
    parser.add_argument('--t3b', action='store_true',
                        help='Explicitly include verified T3-B alpha=-1 RGE in EFT_1_after_F')
    parser.add_argument('--data-dir', type=Path, help='Saved EFT1 RGE data for exact C12 representation projection')
    parser.add_argument('--preserve-existing-html', action='store_true', help='Append/update C12 section without overwriting existing HTML')
    args = parser.parse_args(argv)
    group_dir = args.report_root / 'GroupFactors'
    tex_paths = sorted(group_dir.glob('GF_*.tex'))
    if not tex_paths:
        parser.error(f'No group-factor TeX files under {group_dir}')
    output = write_group_factor_html_companions(
        [], tex_paths, report_root=args.report_root, explicit_t3b=args.t3b,
        representation_data_dir=args.data_dir,
        preserve_existing_html=args.preserve_existing_html,
    )
    for html in output:
        print(f'GroupFactors HTML: {html}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
