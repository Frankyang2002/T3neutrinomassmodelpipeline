"""Self-contained mathematical notation for saved T3 analytical LaTeX reports.

The script does not claim undocumented loop-integral normalisations. It
inserts verified algebraic definitions and writes a transparent audit of
unresolved specialised symbols. No physics calculations are performed.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess

MATH_START = "% BEGIN T3 MATHEMATICAL DEFINITIONS (generated)"
MATH_END = "% END T3 MATHEMATICAL DEFINITIONS (generated)"

# Conventions taken from the report generator's notation and model definitions.
COMMON = r"""
\section*{Mathematical definitions and conventions}
\begin{align*}
Y_{S_1}&=\frac{\alpha}{2},&
Y_F&=\frac{\alpha+1}{2},&
Y_{S_2}&=\frac{\alpha+2}{2},\\
j(d)&=\frac{d-1}{2},&
C_2(d)&=\frac{d^2-1}{4},&
T(d)&=\frac{d(d^2-1)}{12},\\
D_\mu\Phi&=\left(\partial_\mu+i g_1Y_\Phi B_\mu+i g_2T^a_\Phi W^a_\mu\right)\Phi,\\
\widetilde H&=i\sigma^2H^*,&
\beta_X&=\mu\frac{dX}{d\mu},&
\beta_X&=\frac{\beta_X^{(1)}}{16\pi^2}+\cdots,\\
T_\nu^{(1)}&=\operatorname{Tr}(y_1^\dagger y_1),&
T_\nu^{(2)}&=\operatorname{Tr}(y_2^\dagger y_2),\\
T_\nu&=T_\nu^{(1)}+T_\nu^{(2)},&
T&=\operatorname{Tr}(Y_e^\dagger Y_e+3Y_u^\dagger Y_u+3Y_d^\dagger Y_d).
\end{align*}
\begin{align*}
\mathcal L_{Y,\mathrm{SM}}&\supset-\bar Q Y_u\widetilde H u-\bar QY_d H d-\bar L Y_eHe+\mathrm{h.c.},\\
V&\supset m_{S_1}^2S_1^\dagger S_1+m_{S_2}^2 S_2^\dagger S_2
+\lambda_1(H^\dagger H)^2\\
&\quad+\lambda_{HS_1}^{(1)}(H^\dagger H)(S_1^\dagger S_1)
+\lambda_{HS_2}^{(1)}(H^\dagger H)(S_2^\dagger S_2)\\
&\quad+\lambda_{12}^{(1)}(S_1^\dagger S_1)(S_2^\dagger S_2),\\
V_{\mathrm{adj}}&\supset \lambda_{HS_1}^{(A)}(H^\dagger T^aH)(S_1^\dagger T^aS_1)
+\lambda_{HS_2}^{(A)}(H^\dagger T^aH)(S_2^\dagger T^aS_2)\\
&\quad+\lambda_{12}^{(A)}(S_1^\dagger T^aS_1)(S_2^\dagger T^aS_2),\\
V_\times&\supset\lambda_{12}^{(\times)}(S_1^\dagger S_2)(S_2^\dagger S_1).
\end{align*}
"""

# These are algebraic definitions of the report's own notation, not
# identities claimed for all RGE conventions. The I3 function remains
# unresolved until the exact matching-library definition is established.
C5 = r"""
\subsection*{Weinberg-operator and flavour indices}
\begin{align*}
\mathcal L_{\mathrm{EFT}}&\supset\frac{1}{2}\,C_5^{ij}\,(L_iH)(L_jH)+\mathrm{h.c.},\\
C_5^{ij}&=C_5^{ji},\\
\mathcal Y^{ij}&=\sum_r\left(\overline{y_{1,i,r}}\,\overline{y_{2,j,r}}+
\overline{y_{1,j,r}}\,\overline{y_{2,i,r}}\right),\\
\epsilon&=\frac{4-d}{2}.
\end{align*}
\noindent\textbf{Unresolved matching convention:}
\[I_3(m_F^2,m_{S_1}^2,m_{S_2}^2)\quad\text{requires the exact matching-library definition.}\]
"""

GF = r"""
\subsection*{Group-theory contractions}
\begin{align*}
[T^a,T^b]&=i\epsilon^{abc}T^c,&
T^aT^a&=C_2(d)\,\mathbf 1_d,\\
\operatorname{Tr}(T^aT^b)&=T(d)\,\delta^{ab},&
\sum_{a=1}^{3}\delta^{aa}&=3.
\end{align*}
"""


def block_for(path: Path) -> str:
    lines = [COMMON]
    if path.name.lower() == 'c5.tex':
        lines.append(C5)
    if 'GroupFactors' in path.parts:
        lines.append(GF)
    return MATH_START + '\n' + '\n'.join(lines) + '\n' + MATH_END


def insert_definitions(tex: str, path: Path) -> str:
    # Each definition block is delimited and replaced, not duplicated.
    if MATH_START in tex:
        a=tex.index(MATH_START)
        b=tex.index(MATH_END,a)+len(MATH_END)
        tex=tex[:a]+tex[b:]
    # Keep existing symbolic derivations unchanged. Insert immediately
    # after the document begins, before the first analytical section.
    m=re.search(r'\\begin\{document\}',tex)
    if not m:
        raise ValueError(f'Missing LaTeX document start: {path}')
    offset=m.end()
    return tex[:offset]+'\n'+block_for(path)+'\n'+tex[offset:]


def discover_undefined(tex: str) -> list[str]:
    symbols=[]
    for sym, token in [('I3',r'I_3'),('loop_order_hbar',r'\hbar'),('Wilson_tensor',r'C_{ijab}')]:
        if token in tex:
            symbols.append(sym)
    return symbols


def apply(root: Path, *, compile_pdf: bool = False) -> dict:
    root=Path(root)
    paths=sorted(root.rglob('*.tex'))
    if not paths:
        raise FileNotFoundError(f'No .tex reports found in {root}')
    found=[]; unresolved={}; failed=[]
    for path in paths:
        original=path.read_text(encoding='utf-8')
        unresolved[str(path.relative_to(root))]=discover_undefined(original)
        updated=insert_definitions(original,path)
        backup=path.with_suffix('.tex.before_math_definitions.bak')
        if not backup.exists():
            backup.write_text(original,encoding='utf-8')
        path.write_text(updated,encoding='utf-8')
        found.append(str(path.relative_to(root)))
        if compile_pdf:
            exe=shutil.which('pdflatex')
            if not exe:
                raise RuntimeError('pdflatex not found')
            cmd=[exe,'-interaction=nonstopmode','-halt-on-error',path.name]
            output=subprocess.run(cmd,cwd=path.parent,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            if output.returncode:
                failed.append({'file':str(path),'tail':output.stdout[-2500:]})
    manifest={'status':'Complete' if not failed else 'CompilationFailed',
              'reports':found,'unresolved_symbols_by_report':unresolved,
              'compile_errors':failed,
              'important':'I3 matching-library convention not verified; no speculative definition inserted.'}
    (root/'mathematical_definitions_audit.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    return manifest


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('Reports/output/analytical'))
    parser.add_argument('--compile',action='store_true')
    args=parser.parse_args()
    result=apply(args.root,compile_pdf=args.compile)
    print(json.dumps({'status':result['status'],'reports':len(result['reports']),
      'unresolved':{k:v for k,v in result['unresolved_symbols_by_report'].items() if v},
      'compilation_failures':result['compile_errors']},indent=2))
    return 0 if result['status']=='Complete' else 2

if __name__=='__main__':
    raise SystemExit(main())
