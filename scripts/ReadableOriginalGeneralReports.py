"""Presentation-only readability pass over ORIGINAL T3 generated .tex reports.

Preserves the original report tables. Writes a separate Readable/<stem>.tex/.pdf,
never overwrites the full native-source TeX or modifies any RGEs.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT_DIR = ROOT / 'Reports/output/full16_general_symbolic/GeneralSymbolic'


def _paren(text: str, start: int):
    depth = 0
    for i in range(start, len(text)):
        if text[i] == '(':
            depth += 1
        elif text[i] == ')':
            depth -= 1
            if depth == 0:
                return text[start+1:i], i+1
    return None


def _args(value: str):
    parts = []; depth = 0; begin = 0
    for i,ch in enumerate(value):
        if ch == '(' : depth += 1
        elif ch == ')': depth -= 1
        elif ch == ',' and depth == 0:
            parts.append(value[begin:i]); begin=i+1
    parts.append(value[begin:]); return parts


def _readable_field(name: str) -> str:
    return {'H':'H','S_1':'S_1','S_2':'S_2','F':'F','L':'L','Q':'Q'}.get(name,name)


def _operation(tag: str, inside: str) -> str:
    parts = _args(inside)
    head = parts[0].strip() if parts else ''
    if tag == 'Field':
        # Suppress representation and dummy-index wrappers in the *display* only.
        symbol = _readable_field(head)
        if 'Lorentz' in inside: return r'\partial_{\mu} '+symbol
        return symbol
    if tag == 'Coupling': return head
    if tag == 'Index': return r'\mathrm{i}'
    if tag == 'CG':
        # Keep a distinct labelled contraction, rather than discard CG information.
        key = re.search(r'(?:T3ExtraCG|T3GeneralCG)([QM]\d{6}I\d+)',head)
        if key:
            return r'\mathcal{I}_{\mathrm{'+key.group(1)+r'}}'
        key = re.search(r'T3HiggsPortal(\d)CGInv(\d+)',head)
        if key:
            return r'\mathcal{I}_{H'+key.group(1)+r'}^{('+key.group(2)+r')}'
        if 'T3' in head: return r'\mathcal{I}_{\mathrm{T3}}'
        return r'\mathcal{I}_{\mathrm{CG}}'
    if tag == 'Bar':
        return r'{'+head+r'}^{\dagger}'
    if tag == 'Transp': return r'{'+head+r'}^{T}'
    if tag == 'NCM': return ' '.join(parts)
    if tag == 'DiracProduct': return ' '.join(parts)
    if tag == 'Proj': return r'P_{'+(head or 'L')+'}'
    if tag == 'GammaM': return r'\gamma^{\mu}'
    return r'\operatorname{'+tag+r'}\left('+inside+r'\right)'


def compact_native_tex(text: str) -> str:
    """Render Matchete wrappers as symbols, retaining named CG contractions."""
    patt = re.compile(r'\\text\{(Field|Coupling|Index|CG|Bar|Transp|NCM|DiracProduct|Proj|GammaM)\}\(')
    # Innermost-first conversion: handles Field(...Index(...)) and Bar(Field(...)).
    for _ in range(35):
        matches=list(patt.finditer(text))
        if not matches: break
        changed = False
        for m in reversed(matches):
            parsed=_paren(text,m.end()-1)
            if parsed is None: continue
            inside,end=parsed
            # Avoid collapsing outer functions before nested wrappers are simplified.
            if patt.search(inside): continue
            new=_operation(m.group(1),inside)
            text=text[:m.start()]+new+text[end:]
            changed=True
        if not changed: break
    # Some Matchete dummy labels are awkward TeXForm tokens after conversion.
    text=re.sub(r'\\text\{d\\\$\\\$\d+\}',r'\\mathrm{i}',text)
    text=re.sub(r'(\\?\{?[A-Za-z][A-Za-z_0-9]*\}?\^\{\\dagger\})\^(\d+)',r'{\1}^{\2}',text)
    return text


def improve_latex(source: str, *, min_size: str = 'normalsize') -> str:
    """Keep original longtable layout but prohibit its tiny/scaled cells."""
    if min_size not in ('normalsize','small'):
        raise ValueError('min_size must be normalsize or small')
    s=source.replace(r'\documentclass[8pt]{article}',r'\documentclass[10pt,a3paper]{article}')
    s=s.replace(r'\tiny','\\'+min_size)
    s=s.replace(r'\scriptsize','\\'+min_size)
    s=compact_native_tex(s)
    # Disable geometric downscaling: the existing adjustbox otherwise reduces
    # even normalsize expressions to unreadably small glyphs. A3 helps fit.
    s=s.replace(r'\adjustbox{max width=\linewidth}{',r'\adjustbox{max width=1000pt}{')
    # Preserve original expressions as an unmodified TeX source beside readable copy.
    s=s.replace(r'\begin{document}', r'\begin{document}'+'\n'+r'\sloppy')
    return s


def run(report_dir: Path, *, min_size: str='normalsize') -> int:
    sources=[p for p in sorted(report_dir.rglob('*.tex')) if 'Readable' not in p.parts]
    if not sources:
        raise FileNotFoundError(f'No original .tex reports found under {report_dir}')
    latex=shutil.which('pdflatex')
    if not latex:
        raise RuntimeError('pdflatex is required to compile the original reports')
    failed=[]
    for src in sources:
        target=src.parent/'Readable'/src.name
        target.parent.mkdir(exist_ok=True)
        target.write_text(improve_latex(src.read_text(encoding='utf-8'), min_size=min_size),encoding='utf-8')
        command=[latex,'-interaction=nonstopmode','-halt-on-error',target.name]
        result=subprocess.run(command,cwd=target.parent,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,check=False)
        (target.parent/(src.stem+'.compile.log')).write_text(result.stdout,encoding='utf-8')
        exists=target.with_suffix('.pdf').is_file()
        print(f'{"OK" if result.returncode==0 and exists else "FAILED"}: {target.with_suffix(".pdf")}')
        if result.returncode!=0 or not exists: failed.append(str(target))
    if failed:
        print('Compilation failed for:',*failed,sep='\n - ')
        return 1
    return 0


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report-dir',type=Path,default=DEFAULT_REPORT_DIR)
    p.add_argument('--min-size',choices=['normalsize','small'],default='normalsize')
    a=p.parse_args()
    return run(a.report_dir,min_size=a.min_size)

if __name__=='__main__': raise SystemExit(main())
