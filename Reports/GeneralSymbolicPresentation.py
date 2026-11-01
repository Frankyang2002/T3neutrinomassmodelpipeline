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


def notation_rows(models):
    registry = _registry()['models']
    seen = set(); rows = []
    for model in models:
        for category in ('new_rgbeta_symbols', 'preexisting_rgbeta_symbols'):
            for entry in registry[model].get(category, []):
                if not isinstance(entry, dict):continue
                name = entry.get('name', '')
                if name in seen or not name:continue
                seen.add(name)
                symbol = generated_symbol(name)
                content = field_content(entry.get('counts', []))
                detail = r'\texttt{' + latex_escape_text(name) + '}'
                rows.append((symbol, content, detail))
    for name in ('y1Conj', 'y2Conj'):
        rows.append((generated_symbol(name), r'L\,F\,S_i^\dagger', r'\texttt{' + name + '}'))
    return rows
