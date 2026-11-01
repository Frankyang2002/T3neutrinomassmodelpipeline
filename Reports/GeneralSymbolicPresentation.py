"""Presentation-only symbol and TeXForm converters for the original reports."""
from __future__ import annotations
import json
from functools import lru_cache
from pathlib import Path
import re
from Reports.ReportGeneration import latex_escape_text
from scripts.RunT3GeneralSymbolicStudy import ROOT

# Reuse the already implemented saved-expression TeXForm functions.
from Reports.GeneralSymbolicAdapterLegacyConversions import matchete_terms, obtain_beta_latex

QUARTIC = re.compile(r'lambdaT3(?:General|Extra)Q([0-9]{6})I([0-9]+)')
MASS = re.compile(r'(?:mu|m)T3[A-Za-z0-9_]+')


def _registry():
    return json.loads((ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json').read_text(encoding='utf-8'))


@lru_cache(maxsize=1)
def _codes():
    """One deterministic label inventory for both Matchete Extra and RGBeta General.

    Read the existing registry, its Matchete basis aliases, the operator manifest,
    and saved model output. Source expressions are never changed. The files
    are only consulted to establish a stable shared presentation index.
    """
    codes: dict[str, set[int]] = {}

    def register(text):
        for match in QUARTIC.finditer(str(text)):
            codes.setdefault(match[1], set()).add(int(match[2]))

    registry = _registry()
    for model in registry['models'].values():
        for category in ('new_rgbeta_symbols', 'preexisting_rgbeta_symbols'):
            for item in model.get(category, []):
                if not isinstance(item, dict):
                    continue
                register(item.get('name', ''))
                # The Matchete alias can use I# for its invariant counter.
                alias = str(item.get('matchete_basis_name', ''))
                register(alias.replace('I#', 'I1'))
                for invariant in range(1, int(item.get('matchete_count', 1) or 1) + 1):
                    register(alias.replace('I#', 'I' + str(invariant)))

    manifest = ROOT / 'docs/general_scalar_trial_operator_manifest.json'
    if manifest.is_file():
        register(manifest.read_text(encoding='utf-8'))

    # Capture couplings present only in native Matchete output, including
    # aliases not listed as *new* RGBeta vertices. Bounded to saved model files.
    from scripts.RunT3GeneralSymbolicStudy import DEFAULT_DATA_ROOT
    data_root = DEFAULT_DATA_ROOT
    if data_root.is_dir():
        for run_dir in sorted(data_root.iterdir()):
            if not run_dir.is_dir():
                continue
            for model_dir in sorted(run_dir.iterdir()):
                if not model_dir.is_dir():
                    continue
                for file in sorted(model_dir.glob('*.json')) + sorted(model_dir.glob('*.txt')):
                    # Do not read enormous native output files solely to build
                    # notation labels. They are still handled by fallback below.
                    if file.stat().st_size <= 8_000_000:
                        register(file.read_text(encoding='utf-8', errors='replace'))

    return {code: (i, len(indices) > 1)
            for i, (code, indices) in enumerate(sorted(codes.items()), 1)}


def generated_symbol(name):
    m = QUARTIC.fullmatch(name)
    if m:
        code = _codes().get(m[1])
        if not code:
            # Never block report generation for a legitimate saved native
            # interaction absent from the current registry. This fallback is
            # clearly distinguished from the consecutively indexed labels.
            return r'\lambda_{X_{\mathrm{u' + str(int(m[1])) + '}}}^{(' + m[2] + ')}'
        return r'\lambda_{X_{' + str(code[0]) + '}}' + ('^{(' + m[2] + ')}' if code[1] else '')
    if name in ('y1Conj', 'y2Conj'):
        return r'y_{' + name[1] + r'}^{\mathrm{conj}}'
    if MASS.fullmatch(name):
        prefix = r'\mu' if name.startswith('mu') else 'm'
        tail = name[2:] if name.startswith('mu') else name[1:]
        return prefix + r'_{\mathrm{' + tail.replace('_', r'\,') + '}}'
    return None


def format_general_latex(tex):
    names = {m[0] for m in QUARTIC.finditer(tex)} | {m[0] for m in MASS.finditer(tex)}
    names |= {n for n in ('y1Conj', 'y2Conj') if n in tex}
    for name in sorted(names, key=len, reverse=True):
        display = generated_symbol(name)
        if not display: continue
        safe = '{' + display + '}'
        for wrapper in ('text', 'mathrm', 'operatorname', 'texttt'):
            tex = tex.replace('\\' + wrapper + '{' + name + '}', safe)
        tex = re.sub(r'(?<![A-Za-z0-9_])' + re.escape(name) + r'(?![A-Za-z0-9_])',
                     lambda _: safe, tex)
    return tex


def field_content(counts):
    fields = (r'H', r'H^\dagger', r'S_1', r'S_1^\dagger', r'S_2', r'S_2^\dagger')
    return r'\,'.join('(' + f + ')^{' + str(n) + '}' if n > 1 else f
                      for f, n in zip(fields, counts) if n) or '1'


# Names already defined by the original common notation key.
ORIGINAL_DEFINITIONS = frozenset({
    'gY','g2','g3','yu','yd','ye','y1','y2','MF',
    'mS1Sq','mS2Sq','lambdaH','lambdaS1','lambdaS2',
    'lambdaH1','lambdaH2','lambdaH1Adj','lambdaH2Adj',
    'lambdaT3','lambda12','lambda12Adj','lambda12Cross',
    'lambdaS1Adj','lambdaS2Adj',
})

# Exact field content is taken from the internal symbol names. These entries
# are schematic until an explicit invariant contraction is available.
NAMED_QUARTICS = {
    'lambdaHHdagS1barS1bar': (r'(H^\dagger)^2(S_1^\dagger)^2', r'\lambda_7'),
    'lambdaHHdagS2S2': (r'(H^\dagger)^2S_2^2', r'\lambda_8'),
    'lambdaS1bar2S2bar2': (r'(S_1^\dagger)^2(S_2^\dagger)^2', r'\lambda_9'),
    'lambdaS1barS2S2bar2': (r'S_1^\dagger S_2(S_2^\dagger)^2', r'\lambda_{10}'),
    'lambdaS1S1bar2S2bar': (r'S_1(S_1^\dagger)^2S_2^\dagger', r'\lambda_{11}'),
    'lambdaHHdagS1barS2barCross': (r'HH^\dagger S_1^\dagger S_2^\dagger', r'\lambda_{12}'),
}

def _notation_symbol(name):
    from Reports.ReportGeneration import PAPER_SYMBOL_LATEX
    return generated_symbol(name) or PAPER_SYMBOL_LATEX.get(name) or r'\mathrm{' + latex_escape_text(name) + '}'


def notation_rows(models, encountered=()):
    """Complete additional-coupling index for the actual selected study.

    Registry entries supply field content and contraction. Every remaining
    encountered beta key or Lagrangian coupling is included, but unknown
    operator definitions are *not* invented.
    """
    registry = _registry()['models']
    found = set(str(x) for x in encountered if x)
    entries = {}
    for model in models:
        for category in ('new_rgbeta_symbols', 'preexisting_rgbeta_symbols'):
            for entry in registry[model].get(category, []):
                if not isinstance(entry, dict):
                    continue
                name = str(entry.get('name',''))
                if not name: continue
                entries.setdefault(name, entry)
                # Include native Matchete aliases separately when they arise
                # in the saved Lagrangian expression.
                alias = str(entry.get('matchete_basis_name',''))
                if alias and 'I#' in alias:
                    for i in range(1, int(entry.get('matchete_count',1) or 1)+1):
                        entries.setdefault(alias.replace('I#', 'I'+str(i)), entry)
        # No guessed couplings for models absent from this report.
    found.update(entries)
    rows=[]
    for name in sorted(found):
        if name in ORIGINAL_DEFINITIONS:
            continue
        symbol = _notation_symbol(name)
        if name in entries:
            entry=entries[name]
            content=field_content(entry.get('counts', []))
            contraction=str(entry.get('invariant','')).strip()
            detail=r'\texttt{' + latex_escape_text(name) + '}'
            if contraction:
                detail += r'; contraction: \texttt{' + latex_escape_text(contraction) + '}'
        elif name in NAMED_QUARTICS:
            content, symbol = NAMED_QUARTICS[name]
            detail=r'\texttt{' + latex_escape_text(name) + r'}; invariant contraction model-dependent'
        elif name == 'y1Conj':
            content=r'L\,F\,S_1^\dagger'
            detail=r'\texttt{y1Conj}; conjugate-scalar Yukawa'
        elif name == 'y2Conj':
            content=r'L\,F\,S_2^\dagger'
            detail=r'\texttt{y2Conj}; conjugate-scalar Yukawa'
        elif MASS.fullmatch(name):
            content=r'\mathrm{scalar\ bilinear}'
            detail=r'\texttt{' + latex_escape_text(name) + r'}; see native Lagrangian for contraction'
        elif re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',name):
            content=r'\text{definition unavailable in registry}'
            detail=r'\texttt{' + latex_escape_text(name) + r'}; symbolic running parameter'
        else:
            continue
        rows.append((symbol,content,detail))
    return rows
