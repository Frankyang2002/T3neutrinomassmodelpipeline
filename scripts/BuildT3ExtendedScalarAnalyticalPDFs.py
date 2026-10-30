"""Read-only PDF/TeX analytical addenda for the general-scalar T3 extension.

Sources (all must exist in the LOCAL T3 repository):
  RGE/general/generated/general_rgbeta_quartic_registry.json
  docs/general_scalar_trial_operator_manifest.json
  output/general_t3_quadratic_eft1/<model>_{uv,eft1}_quadratic_live_rgbeta.json

Output (new, does NOT overwrite the historical reports):
  Reports/output/analytical/extended_scalar/{Lagrangian,GroupFactors,RGE}/

These are *candidate-basis analytical audits*, not a certification of the
RGBeta/Matchete normalisation or sequential EFT matching.

Run from repository root:
    python -m scripts.BuildT3ExtendedScalarAnalyticalPDFs --check
    python -m scripts.BuildT3ExtendedScalarAnalyticalPDFs
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any

import sympy as sp
from sympy.parsing.mathematica import parse_mathematica
from sympy.printing.latex import LatexPrinter

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
OPERATORS = ROOT / 'docs/general_scalar_trial_operator_manifest.json'
BETA_DIR = ROOT / 'output/general_t3_quadratic_eft1'
OUTPUT = ROOT / 'Reports/output/analytical/extended_scalar'
MODEL_ORDER = (
    'A_-4','A_-2','A_+0','B_-3','B_-1','B_+1',
    'C_-3','C_-1','C_+1','D_-2','D_+0','D_+2',
    'E_-4','E_-2','E_+0','E_+2',
)
BETA_FORMAT = 'general-t3-quadratic-rgbeta-experimental-v1'
BETA_STATUS = 'SuccessCandidateMassBeta'
BETA_CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'
MASS_KEYS = ('muT3GeneralS1Sq','muT3GeneralS2Sq','muT3GeneralS1S2')


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('Extended scalar PDF reports: ' + message)


def read_json(path: Path) -> dict:
    require(path.is_file(), 'missing source ' + str(path))
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data, dict), 'expected JSON object in ' + str(path))
    return data


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tag(model: str) -> str:
    return model.replace('+', 'p').replace('-', 'm')


def tex_escape(text: Any) -> str:
    return (str(text).replace('\\', r'\textbackslash{}')
            .replace('&', r'\&').replace('%', r'\%')
            .replace('$', r'\$').replace('#', r'\#')
            .replace('_', r'\_').replace('{', r'\{')
            .replace('}', r'\}').replace('~', r'\textasciitilde{}')
            .replace('^', r'\textasciicircum{}'))


def model_tex(model: str) -> str:
    cls, alpha = model.split('_')
    return rf'\mathrm{{T3-{cls}}}\,(\alpha={int(alpha)})'


def field_content(counts: list[int]) -> str:
    require(isinstance(counts,list) and len(counts)==6 and sum(counts) in (2,4),
            'bad scalar operator field counts')
    bases=(r'H',r'H^\dagger',r'S_1',r'S_1^\dagger',r'S_2',r'S_2^\dagger')
    parts=[]
    for i,n in enumerate(counts):
        if n:
            # Parentheses prevent a double-superscript on conjugated fields.
            parts.append((rf'\left({bases[i]}\right)^{{{n}}}' if i%2 and n>1
                          else bases[i] + (rf'^{{{n}}}' if n>1 else '')))
    return r'\, '.join(parts) or '1'


# The numbers 7--12 are already used by the original analytical-report
# notation for the structurally reused B/C quartics.  Never repurpose them.
# Extra *new* quartics are assigned 13+ globally by their RGBeta source names,
# so the same source operator receives the same display label in every model.
REUSED_PAPER_NAMES = {
    'lambdaHHdagS1barS1bar': r'\lambda_{7}',
    'lambdaHHdagS2S2': r'\lambda_{8}',
    'lambdaS1bar2S2bar2': r'\lambda_{9}',
    'lambdaS1barS2S2bar2': r'\lambda_{10}',
    'lambdaS1S1bar2S2bar': r'\lambda_{11}',
    'lambdaHHdagS1barS2barCross': r'\lambda_{12}',
}
BASE_PAPER_NAMES = {
    'lambdaT3': r'\lambda_{5}',
    'lambdaH': r'\lambda_{1}',
    'lambdaS1': r'\lambda_{S_1}^{(1)}',
    'lambdaS2': r'\lambda_{S_2}^{(1)}',
    'lambdaS1Adj': r'\lambda_{S_1}^{(A)}',
    'lambdaS2Adj': r'\lambda_{S_2}^{(A)}',
    'lambdaH1': r'\lambda_{HS_1}^{(1)}',
    'lambdaH2': r'\lambda_{HS_2}^{(1)}',
    'lambdaH1Adj': r'\lambda_{HS_1}^{(A)}',
    'lambdaH2Adj': r'\lambda_{HS_2}^{(A)}',
    'lambda12': r'\lambda_{12}^{(1)}',
    'lambda12Adj': r'\lambda_{12}^{(A)}',
    'lambda12Cross': r'\lambda_{12}^{(\times)}',
    'mS1Sq': r'm_1^2',
    'mS2Sq': r'm_2^2',
    'muT3GeneralS1Sq': r'\mu_{S_1}^2',
    'muT3GeneralS2Sq': r'\mu_{S_2}^2',
    'muT3GeneralS1S2': r'\mu_{12}^2',
    'gY': r'g_1', 'g2': r'g_2', 'g3': r'g_3',
    'y1': r'y_1', 'y2': r'y_2',
    'yu': r'Y_u', 'yd': r'Y_d', 'ye': r'Y_e',
    'MF': r'M_F',
}


def notation_map(reg: dict) -> dict[str, str]:
    """Single, injective paper notation for every additional scalar symbol.

    Source identifiers survive in manifest.json, not in typeset report pages.
    Explicitly exclude lambda5 and the previously numbered lambda7--lambda12.
    """
    names = sorted({item['name'] for m in MODEL_ORDER
                    for item in reg['models'][m]['new_rgbeta_symbols']})
    require(all(re.fullmatch(r'lambdaT3GeneralQ[0-9]{6}I[0-9]+', n) for n in names),
            'unrecognised new quartic symbol; update notation explicitly')
    extras = {name: rf'\lambda_{{{i}}}' for i, name in enumerate(names, start=13)}
    result = dict(BASE_PAPER_NAMES, **REUSED_PAPER_NAMES, **extras)
    used = set()
    for m in MODEL_ORDER:
        e = reg['models'][m]
        for item in e['new_rgbeta_symbols']:
            used.add(item['name'])
        for item in e['preexisting_rgbeta_symbols']:
            name = item['existing_name']
            require(name in REUSED_PAPER_NAMES,
                    f'{m}: no established paper symbol for reused {name}')
            used.add(name)
    require(len({result[n] for n in used}) == len(used),
            'two distinct extra quartics were assigned the same display label')
    return result


def pretty_symbol(name: str, notation: dict[str, str]) -> str:
    # Unmapped lambda names are a hard failure, rather than raw implementation
    # names leaking into a thesis PDF.
    require(name in notation, 'unmapped scalar notation symbol '+name)
    return notation[name]


def literal(text: Any) -> str:
    """Verbatim-safe ASCII string used only for the invariant tensor basis."""
    return str(text).encode('ascii','backslashreplace').decode('ascii')


def conjugated_symbol_tex(symbol: str) -> str:
    """Render the conjugate of a paper-notation symbol as valid LaTeX.

    A symbol such as ``\\mu_{S_2}^2`` or ``\\lambda_{HS_1}^{(1)}``
    already contains a superscript. Appending ``^{*}`` directly creates
    TeX's ``Double superscript`` error. Parenthesise the entire symbol in
    that case; for plain symbols retain the traditional ``\\lambda_i^*``.
    This is display-only and does not modify the RGBeta expressions.
    """
    if '^' in symbol:
        return r'\left(' + symbol + r'\right)^{*}'
    return symbol + r'^{*}'


def listing(text: Any, notation: dict[str,str] | None = None) -> str:
    """Preserve every native beta coefficient while typesetting readable symbols.

    listings/mathescape allows TeX math inside $...$; token substitution does
    not simplify, reorder, or change the Mathematica InputForm expression.
    Both the full source and exact alias dictionary are retained in JSON.
    """
    source = literal(text).replace('\\end{lstlisting}', 'END_LISTING')
    if notation is not None:
        # Recognise Mathematica Bar[singleSymbol] before replacing individual
        # names so conjugation prints as lambda_i^* or Y^*.
        def conj(match: re.Match[str]) -> str:
            name=match.group(1)
            return '$'+conjugated_symbol_tex(notation[name])+'$' if name in notation else match.group(0)
        source = re.sub(r'\bBar\[([A-Za-z][A-Za-z0-9]*)\]',conj,source)
        def rename(match: re.Match[str]) -> str:
            name=match.group(0)
            if name in notation:
                return '$'+notation[name]+'$'
            require(not name.startswith(('lambda', 'muT3General')),
                    'unmapped source symbol in typeset RGE: '+name)
            return name
        # Never replace inside inserted math-mode fragments.
        pieces=source.split('$')
        for i in range(0,len(pieces),2):
            pieces[i]=re.sub(r'\b[A-Za-z][A-Za-z0-9]*\b',rename,pieces[i])
        source='$'.join(pieces)
    return ('\\begin{lstlisting}[mathescape=true]\n'+source+
            '\n\\end{lstlisting}\n')



class _T3LatexPrinter(LatexPrinter):
    """Mathematica InputForm -> mathematical LaTeX, using report notation."""

    def __init__(self, notation: dict[str, str]):
        super().__init__({'mul_symbol': ' '})
        self.notation = notation

    def _print_Symbol(self, expr):
        name = str(expr)
        if name in self.notation:
            return self.notation[name]
        if name.startswith(('lambda', 'muT3General')):
            raise ValueError('Unmapped RGBeta coupling: ' + name)
        return super()._print_Symbol(expr)

    def _print_Function(self, expr):
        name = expr.func.__name__
        args = list(expr.args)
        if name == 'Bar' and len(args) == 1:
            inner = self._print(args[0])
            return conjugated_symbol_tex(inner) if args[0].is_Symbol else r'\overline{\left(' + inner + r'\right)}'
        if name == 'Trans' and len(args) == 1:
            return r'\left(' + self._print(args[0]) + r'\right)^{T}'
        if name == 'Tr' and len(args) == 1:
            return r'\operatorname{Tr}\!\left[' + self._print(args[0]) + r'\right]'
        if name == 'Dot':
            return r'\,'.join(self._print(a) for a in args)
        if name == 'Conjugate' and len(args) == 1:
            return r'\overline{' + self._print(args[0]) + '}'
        return super()._print_Function(expr)

def math_expression(text: str, notation: dict[str, str]) -> sp.Expr:
    try:
        expr = parse_mathematica(text)
    except Exception as exc:
        raise ValueError('Cannot parse RGBeta InputForm expression: '+text[:180]) from exc
    # The source is algebraically parsed, not copied as monospace text.
    # Reject unsupported function names instead of printing Mathematica syntax.
    valid = {'Bar', 'Trans', 'Tr', 'Dot', 'Conjugate'}
    unknown = {f.func.__name__ for f in expr.atoms(sp.Function) if f.func.__name__ not in valid}
    require(not unknown, 'Unsupported symbolic functions in beta: '+repr(unknown))
    return expr

def beta_equation(label: str, input_form: str, notation: dict[str,str]) -> str:
    """A page-breakable, source-exact symbolic equation with one additive term per row."""
    printer = _T3LatexPrinter(notation)
    expr = math_expression(input_form, notation)
    # Add is associative/commutative: formatting may reorder additive terms but
    # never changes numerical coefficients or mathematical contents.
    terms = list(expr.as_ordered_terms()) if expr.is_Add else [expr]
    rows = []
    for index, term in enumerate(terms):
        is_negative = term.could_extract_minus_sign()
        absolute = -term if is_negative else term
        sign = '-' if is_negative else ('+' if index else '')
        term_tex = printer.doprint(absolute)
        content = (r'\beta_{' + label + r'}^{(1)} = ' if index == 0 else r'\phantom{\beta_{' + label + r'}^{(1)}}\;')
        content += (sign + r'\;' if sign else '') + term_tex
        rows.append(r'\noindent\adjustbox{max width=\linewidth}{$\displaystyle ' + content + r'$}\par')
    return '\n'.join(rows) + '\n'

def notation_key(notation: dict[str,str], reg: dict) -> list[str]:
    """Typeset notation only; raw source-symbol map is saved in manifest.json."""
    out=[r'\section*{Scalar coupling notation}',
         r'The original topology interaction remains $\lambda_5$. Original '
         r'scalar and portal interactions retain the notation of the earlier '
         r'analytical reports. Additional quartic structures use $\lambda_7$--'
         r'$\lambda_{12}$ where these symbols already exist, followed by '
         r'$\lambda_{13},\lambda_{14},\ldots$ for newly named operators. '
         r'All labels refer to RGBeta potential coefficients in the saved source '
         r'convention; the Matchete--RGBeta mapping is not fully certified.',
         r'\begin{longtable}{@{}p{0.21\linewidth}p{0.46\linewidth}p{0.23\linewidth}@{}}',
         r'\toprule',r'Symbol & Scalar field structure & Models \\',
         r'\midrule',r'\endhead']
    used={}
    for m in MODEL_ORDER:
        e=reg['models'][m]
        for item in e['new_rgbeta_symbols']:
            key=item['name'];counts=item['counts']
            if key in used:
                require(used[key]['counts']==counts, 'inconsistent field content for '+key)
                used[key]['models'].append(m)
            else: used[key]={'counts':counts,'models':[m]}
        for item in e['preexisting_rgbeta_symbols']:
            key=item['existing_name'];counts=item['counts']
            if key in used:
                require(used[key]['counts']==counts,'inconsistent reused field content for '+key)
                used[key]['models'].append(m)
            else: used[key]={'counts':counts,'models':[m]}
    for key in sorted(used,key=lambda n:(int(re.search(r'\d+',notation[n]).group()),n)):
        e=used[key]
        models=', '.join(m.replace('_','') for m in e['models'])
        out.append(rf'${notation[key]}$ & ${field_content(e["counts"])}$ & {tex_escape(models)} \\')
    out.extend([r'\bottomrule',r'\end{longtable}',
      r'For distinct SU(2) contractions with the same field content, numbered '
      r'superscripts/index labels remain independent. These symbols denote '
      r'model-dependent invariant tensors, not universal Clebsch coefficients.'])
    return out


def document(title: str, overview: str, content: list[str]) -> str:
    header = [
        r'\documentclass[9pt,a4paper]{article}',
        r'\usepackage[a4paper,margin=17mm]{geometry}',
        r'\usepackage[T1]{fontenc}',
        r'\usepackage{amsmath,amssymb,booktabs,longtable,array,hyperref,listings,xcolor,adjustbox}',
        r'\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue}',
        r'\lstset{basicstyle=\ttfamily\footnotesize,breaklines=true,breakatwhitespace=false,columns=fullflexible,keepspaces=true,frame=single,framesep=4pt}',
        r'\setlength{\parindent}{0pt}',
        r'\setlength{\parskip}{5pt}',
        r'\sloppy',
        r'\begin{document}',
        rf'\section*{{{tex_escape(title)}}}',
        r'\textbf{Provenance:} saved 16-model operator registry and candidate one-loop RGBeta UV/EFT1 exports. '
        r'Source beta expressions are algebraically parsed from saved InputForm and typeset in paper coupling notation. No implicit Matchete rescaling.',
        r'\textbf{Status:} candidate tensor basis. General scalar RGBeta--Matchete mapping, '
        r'quadratic mass-basis matching and full sequential Wilson matching are not certified.',
        overview,
        r'\bigskip',
    ]
    return '\n'.join(header + content + [r'\end{document}', ''])


def validate_registry(reg: dict, manifest: dict) -> None:
    require(set(reg.get('models',{})) == set(MODEL_ORDER), 'registry must contain exactly 16 models')
    totals = reg.get('totals',{})
    require(totals.get('matchete_extra_quartics')==50 and
            totals.get('rgbeta_new_vertices')==38 and
            totals.get('rgbeta_preexisting_structural')==12 and
            totals.get('quadratic_couplings_blocked')==8,
            'unexpected general scalar registry totals')
    require(set(manifest.get('models',{}))==set(MODEL_ORDER), 'operator manifest differs from 16 models')
    require(sum(int(reg['models'][m]['total_additional_quartics']) for m in MODEL_ORDER)==50,
            'registry extra quartic count inconsistent')
    for model in MODEL_ORDER:
        e=reg['models'][model]
        names=[x['name'] for x in e['new_rgbeta_symbols']]
        reused=[x['existing_name'] for x in e['preexisting_rgbeta_symbols']]
        require(len(names)+len(reused)==e['total_additional_quartics'] and
                len(set(names+reused))==len(names)+len(reused) and 'lambdaT3' not in names+reused,
                'invalid new/reused split for '+model)
        require(manifest['models'][model]['extra_quartic_complex_couplings']==e['total_additional_quartics'],
                'matchete operator count mismatch '+model)


def beta_data(model: str, stage: str, reg: dict, beta_dir: Path) -> tuple[dict, Path]:
    assert stage in ('UV','EFT1')
    path = beta_dir / f'{tag(model)}_{stage.lower()}_quadratic_live_rgbeta.json'
    p = read_json(path)
    e = reg['models'][model]
    meta = p.get('metadata',{})
    require(p.get('format')==BETA_FORMAT and p.get('status')==BETA_STATUS,
            f'{model}/{stage}: unexpected candidate beta export status')
    require(meta.get('GeneralScalarModel')==model and meta.get('GeneralScalarStage')==stage,
            f'{model}/{stage}: beta export model/stage differs')
    require(meta.get('GeneralScalarAdditionalQuartics')==e['total_additional_quartics'] and
            meta.get('ReportBetaConvention')==BETA_CONVENTION,
            f'{model}/{stage}: source beta inventory/convention differs')
    b=p.get('report_betas')
    require(isinstance(b,dict) and b and set(b)==set(p.get('beta_symbols',[])),
            f'{model}/{stage}: missing/incomplete beta dictionary')
    require(all(isinstance(s,str) and s.strip() and not any(bad in s for bad in
            ('BetaTerm[','Missing[','$Failed','$Aborted')) for s in b.values()),
            f'{model}/{stage}: unresolved beta expression')
    new=[x['name'] for x in e['new_rgbeta_symbols']]
    reused=[x['existing_name'] for x in e['preexisting_rgbeta_symbols']]
    require(set(new+reused+['lambdaT3']).issubset(b),f'{model}/{stage}: essential scalar beta absent')
    require(p.get('physical_mapping_certified') is False,
            f'{model}/{stage}: unexpectedly changed certification metadata')
    return p,path


def inventory_table(reg: dict) -> list[str]:
    out=[r'\begin{longtable}{@{}lrrrrrr@{}}',r'\toprule',
         r'Model & $d_{S_1}$ & $d_{S_2}$ & $d_F$ & $\alpha$ & extra $\lambda$ & extra $m^2$\\',
         r'\midrule',r'\endfirsthead',r'\toprule',
         r'Model & $d_{S_1}$ & $d_{S_2}$ & $d_F$ & $\alpha$ & extra $\lambda$ & extra $m^2$\\',
         r'\midrule',r'\endhead']
    for m in MODEL_ORDER:
        x=reg['models'][m];a,b,c=x['dimensions']
        out.append(rf'${model_tex(m)}$ & {a} & {b} & {c} & {x["alpha"]} & {x["total_additional_quartics"]} & {x["quadratics_not_implemented"]} \\')
    out += [r'\bottomrule',r'\end{longtable}']
    return out


def model_quartic_sections(reg: dict, manifest: dict, *, include_invariants: bool, notation: dict[str,str]) -> list[str]:
    blocks=[]
    for m in MODEL_ORDER:
        x=reg['models'][m]
        blocks.append(rf'\subsection*{{$ {model_tex(m)} $}}')
        blocks.append(r'\begin{longtable}{@{}p{0.42\linewidth}p{0.36\linewidth}p{0.13\linewidth}@{}}')
        blocks += [r'\toprule',r'Coupling / class & field content & count / basis\\',r'\midrule',r'\endhead']
        for e in x['new_rgbeta_symbols']:
            blocks.append(rf'${pretty_symbol(e["name"],notation)}$ (new) & ${field_content(e["counts"])}$ & {e.get("matchete_count","?")} \\')
        for e in x['preexisting_rgbeta_symbols']:
            blocks.append(rf'${pretty_symbol(e["existing_name"],notation)}$ (reused) & ${field_content(e["counts"])}$ & 1 candidate coordinate \\')
        if x['total_additional_quartics']==0:
            blocks.append(r'\multicolumn{3}{l}{No additional scalar quartic operators in this assignment.}\\')
        blocks += [r'\bottomrule',r'\end{longtable}']
        q=manifest['models'][m].get('quadratic_structures_not_enabled',[])
        if q:
            blocks.append(r'\textbf{Additional quadratic operators:}')
            blocks.append(', '.join(r'$'+field_content(item['counts'])+r'$' for item in q) + r' (candidate, with Hermitian conjugates).')
        if include_invariants:
            for e in x['new_rgbeta_symbols']:
                blocks += [rf'\textbf{{$ {pretty_symbol(e["name"],notation)} $}}: candidate SU(2) contraction in registry notation:',
                           listing(e.get('invariant','(not supplied)'),notation)]
    return blocks


def lagrangian_content(reg: dict, manifest: dict, notation: dict[str,str]) -> list[str]:
    out=[r'\section*{Definition and coverage}',
         r'The extended BSM scalar sector is $V_{\mathrm{scalar}}=V_{\mathrm{restricted}}+'
         r'\Delta V_4+\Delta V_2$. The restricted couplings, including $\lambda_5$, are '
         r'\emph{retained}. $\Delta V_4$ consists of 50 model-dependent complex '
         r'quartic structures (38 newly named RGBeta vertices plus 12 structurally reused '
         r'coordinates, counting model occurrences). $\Delta V_2$ catalogues eight '
         r'possible holomorphic bilinears.',
         r'\[\Delta V_4=\sum_{a\in\mathcal{Q}}\left(\lambda_a\mathcal{I}_a+'
         r'\lambda_a^*\mathcal{I}_a^\dagger\right),\qquad '
         r'\Delta V_2=\sum_{b\in\mathcal{B}}\left(\mu_b^2\mathcal{J}_b+'
         r'(\mu_b^2)^*\mathcal{J}_b^\dagger\right).\]',
         r'The monomials and invariant tensors below come from the generated registry. '
         r'They do \emph{not} establish a common Matchete/RGBeta interaction coefficient normalization.',
         r'\section*{16-model scalar operator inventory}',*inventory_table(reg),
         *notation_key(notation,reg),r'\section*{Additional operators by model}',*model_quartic_sections(reg,manifest,include_invariants=False,notation=notation)]
    return out


def groupfactor_content(reg: dict, manifest: dict, betas: dict, notation: dict[str,str]) -> list[str]:
    out=[r'\section*{What can be reported reliably}',
         r'The original $SU(2)_L\times U(1)_Y$ representations have not changed. '
         r'The existing gauge and Yukawa group-factor tables therefore remain the '
         r'reference for the original interactions at the same loop order.',
         r'\[C_2(d)=\frac{d^2-1}{4},\qquad T(d)=\frac{d(d^2-1)}{12}.\]',
         r'For the extra scalar operators, the registry supplies candidate invariant '
         r'tensor contractions and multiplicities. These are listed below. '
         r'\textbf{They are not independent, certified Matchete-to-RGBeta group factors.} '
         r'The one-loop coefficients occurring in the exported beta equations are '
         r'source-calculated candidate coefficients, not re-derived in this report.',
         *notation_key(notation,reg),r'\section*{Operator and invariant tensor inventory}',*inventory_table(reg),
         *model_quartic_sections(reg,manifest,include_invariants=True,notation=notation),
         r'\section*{How the extra scalar terms enter topology-quartic running}',
         r'The following are the complete saved candidate UV expressions for '
         r'$16\pi^2\,d\lambda_{5}/d\ln\mu$, typeset in mathematical notation. '
         r'Use these as an audit of source coefficients; this is not an independent '
         r'group-factor proof.']
    for m in MODEL_ORDER:
        b=betas[(m,'UV')]['report_betas']
        out += [rf'\subsection*{{$ {model_tex(m)} $}}', beta_equation(notation['lambdaT3'], b['lambdaT3'], notation)]
    return out


def rge_content(reg: dict, betas: dict, notation: dict[str,str], *, scope: str) -> list[str]:
    out=[r'\section*{Conventions and completeness}',
         r'All expressions below are parsed from saved native RGBeta InputForm and typeset in mathematical notation '
         r'without changing field convention or polynomial normalisation:',
         r'\[16\pi^2\frac{d X}{d\ln\mu}=\beta_X^{(1)}.\]',
         r'\textbf{UV}: SM+$F$+$S_1$+$S_2$. '
         r'\textbf{EFT1}: fermion integrated out, both scalars retained. '
         r'These are the saved candidate renormalisable beta functions; '
         r'no Weinberg or $C_{12}$ Wilson evolution is inferred from them.',
         r'\textbf{Scope}: '+('all exported renormalisable beta variables.' if scope=='all' else
         r'all scalar quartics, scalar mass squares and holomorphic bilinears, including the original $\lambda_5$.')]
    out += notation_key(notation,reg)
    for stage in ('UV','EFT1'):
        out += [r'\clearpage',rf'\section*{{{stage} candidate one-loop beta functions}}']
        for m in MODEL_ORDER:
            out += [rf'\subsection*{{$ {model_tex(m)} $}}']
            p=betas[(m,stage)]
            keys=list(p['beta_symbols'])
            if scope=='scalar':
                keys=[k for k in keys if k.startswith('lambda') or k in ('mS1Sq','mS2Sq')+MASS_KEYS]
            require(keys, f'{m}/{stage}: no beta symbols selected')
            for key in keys:
                require(key in notation, f'{m}/{stage}: no paper symbol for {key}')
                out += [r'\medskip', beta_equation(notation[key], p['report_betas'][key], notation)]
    return out


def append_to_existing_report(category: str, *, appendix: list[str]) -> str:
    """Optional copy of the historical report with new appendix; never edits source."""
    names={'Lagrangian':'UV.tex','GroupFactors':'GF_UV.tex','RGE':'UV.tex'}
    require(category in names, 'unknown historical report category')
    old=ROOT/'Reports'/'output'/'analytical'/category/names[category]
    require(old.is_file(), 'historical TeX not available '+str(old))
    source=old.read_text(encoding='utf-8-sig')
    require(source.count(r'\end{document}')==1, 'unexpected historical TeX end')
    require(r'\begin{document}' in source, 'historical TeX preamble missing')
    if r'\usepackage{listings}' not in source:
        source=source.replace(r'\begin{document}',
            r'\usepackage{listings}'+'\n'+
            r'\lstset{basicstyle=\ttfamily\footnotesize,breaklines=true,breakatwhitespace=false,columns=fullflexible,keepspaces=true,frame=single,framesep=4pt}'+'\n'+
            r'\begin{document}',1)
    addition='\n'+r'\clearpage'+ '\n'+r'\section*{Extended scalar candidate-basis addendum}'+ '\n'+\
             r'\textbf{Uncertified RGBeta--Matchete normalisation and full sequential matching.}'+'\n'+\
             '\n'.join(appendix)+'\n'
    return source.replace(r'\end{document}',addition+r'\end{document}',1)


def compile_tex(source: Path, *, pdflatex: str) -> Path:
    command=[pdflatex,'-interaction=nonstopmode','-halt-on-error','-file-line-error',
             '-output-directory='+str(source.parent),source.name]
    for _ in range(2):
        p=subprocess.run(command,cwd=source.parent,capture_output=True,text=True,
                         encoding='utf-8',errors='replace')
        if p.returncode:
            raise RuntimeError(f'LaTeX compilation failed for {source}:\n'+
                               '\n'.join((p.stdout+'\n'+p.stderr).splitlines()[-36:]))
    output=source.with_suffix('.pdf')
    require(output.is_file() and output.stat().st_size>1000,
            'LaTeX reported success but PDF is absent '+str(output))
    return output


def run(argv: list[str] | None = None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--check',action='store_true',help='read-only completeness/provenance check')
    ap.add_argument('--rge-scope',choices=('scalar','all'),default='scalar',
                    help='all includes gauge/Yukawa/fermion beta functions as well')
    ap.add_argument('--output-root',type=Path,default=OUTPUT)
    ap.add_argument('--beta-dir',type=Path,default=BETA_DIR)
    ap.add_argument('--tex-only',action='store_true',help='write TeX but do not run pdflatex')
    ap.add_argument('--append-original',action='store_true',
                    help='copy the historical TeX first, then append extended scalar material; files stay separate')
    args=ap.parse_args(argv)
    reg=read_json(REGISTRY);manifest=read_json(OPERATORS)
    validate_registry(reg,manifest)
    notation=notation_map(reg)
    betas={};paths=[REGISTRY,OPERATORS]
    for m in MODEL_ORDER:
        for stage in ('UV','EFT1'):
            value,path=beta_data(m,stage,reg,args.beta_dir)
            betas[(m,stage)]=value
            paths.append(path)
    print(f'Verified {len(MODEL_ORDER)} models, 50 quartic slots, 8 quadratic slots, '
          f'{len(betas)} UV/EFT1 beta exports.',flush=True)
    for m in MODEL_ORDER:
        print(f'  {m:<5}: {reg["models"][m]["total_additional_quartics"]} quartics, '
              f'{reg["models"][m]["quadratics_not_implemented"]} quadratic',flush=True)
    if args.check:return 0
    exe=shutil.which('pdflatex')
    if not args.tex_only:require(bool(exe),'pdflatex missing; install MiKTeX or use --tex-only')
    dest=args.output_root.resolve()
    reports=[
        ('Lagrangian','UV_extended_scalar',
         'T3 UV Lagrangian: additional scalar operators',
         r'\textit{Model-dependent scalar-sector addendum; the original UV, '
         r'Yukawa and gauge terms are unchanged and documented separately.}',
         lagrangian_content(reg,manifest,notation)),
        ('GroupFactors','GF_UV_extended_scalar',
         'T3 UV group factors: additional scalar tensor structures',
         r'\textit{Operator tensor basis and candidate beta coefficients. '
         r'This does not claim an independently normalised Matchete/RGBeta group-factor map.}',
         groupfactor_content(reg,manifest,betas,notation)),
        ('RGE','RGE_UV_EFT1_extended_scalar',
         'T3 candidate UV and EFT1 scalar RGEs',
         r'\textit{Direct transcription of generated candidate one-loop RGBeta exports; '
         r'no new symbolic RGE computation is performed.}',
         rge_content(reg,betas,notation,scope=args.rge_scope)),
    ]
    created=[]
    for sub,stem,title,description,sections in reports:
        parent=dest/sub;parent.mkdir(parents=True,exist_ok=True)
        if args.append_original:
            stem += '_with_original'
            content=append_to_existing_report(sub,appendix=sections)
        else:
            content=document(title,description,sections)
        tex=parent/(stem+'.tex');tex.write_text(content,encoding='utf-8')
        print('TeX:',tex,flush=True)
        created.append(str(tex.relative_to(ROOT)) if tex.is_relative_to(ROOT) else str(tex))
        if not args.tex_only:
            pdf=compile_tex(tex,pdflatex=exe)
            print('PDF:',pdf,flush=True)
            created.append(str(pdf.relative_to(ROOT)) if pdf.is_relative_to(ROOT) else str(pdf))
    provenance={'status':'CandidateAnalyticalAddendumGenerated',
      'models':list(MODEL_ORDER),'quartic_slots':50,'quadratic_slots':8,
      'rge_scope':args.rge_scope,'rge_beta_stages':['UV','EFT1'],
      'full_general_scalar_physical_matching_certified':False,
      'extra_group_factors_independently_certified':False,
      'sources':{str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p):sha256(p) for p in paths},
      'created':created,
      'display_notation':notation,
      'notation_policy':'Existing original and B/C symbols retained; distinct new quartics numbered globally from 13.'}
    (dest/'manifest.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    return 0


if __name__=='__main__':
    raise SystemExit(run())
