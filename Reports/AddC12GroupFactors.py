"""Add the exact tensor-level C12 RGE prefactors to an existing EFT1 GF TeX.

Produces a new TeX/PDF alongside the source. It does NOT claim to provide
representation-projected CG coefficients: those require separate recoupling.
The terms are transcribed from RGE.general.WilsonTensorRGE, not a fitted model.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'Reports/output/analytical/GroupFactors/GF_EFT_1_after_F.tex'

SECTION = r'''
\clearpage
\section*{EFT1 $C_{12}$ / $C_{LLSS}$ tensor RGE: structural group factors}
\textbf{Scope.} The entries below are coefficients and index multiplicities
of the complete real-scalar $\psi^2\phi^2$ tensor beta function used by
\texttt{RGE.general.WilsonTensorRGE}. They are \emph{not} representation-specific
Clebsch--Gordan projections, which require the model's scalar basis and
projector. The Wilson tensor $C_{ijab}$ has lepton-flavour indices $i,j$
and real-scalar indices $a,b$.

\begin{longtable}{@{}p{0.29\linewidth}p{0.13\linewidth}p{0.49\linewidth}@{}}
\toprule
Structure & Factor & Tensor structure (indices summed) \\
\midrule
Yukawa wavefunction & $+1$ & $y_{ijd}y^*_{kld}C_{klab}$ \\
Scalar pair gauge & $-1$ & $\sum_{(a,b)}g_\alpha^2\theta^A_{ac}\theta^A_{bd}C_{ijcd}$ \\
Scalar pair quartic & $+\frac{1}{2}$ & $\sum_{(a,b)}\lambda_{abcd}C_{ijcd}$ \\
Mixed Yukawa (first) & $+2$ & $\sum_{(a,b)\times(i,j)}y_{jkc}y^*_{klb}C_{ilac}$ \\
Mixed Yukawa (second) & $+1$ & $\sum_{(a,b)\times(i,j)}y_{jkb}y^*_{klc}C_{ilac}$ \\
Mixed gauge & $+4$ & $\sum_{(a,b)\times(i,j)}g_\alpha^2(t^A)_{lj}(\theta^A)_{bc}C_{ilac}$ \\
Crossed Yukawa & $+1$ & $\sum_{(a,b),(i,j)\,\mathrm{simult.}}y_{jla}y_{ikc}C_{lkbc}$ \\
Conjugate Wilson (first) & $+2$ & $y_{ild}y_{jkd}C^*_{klab}$ \\
Conjugate Wilson (second) & $+1$ & $y_{ijd}y_{kld}C^*_{klab}$ \\
Scalar anomalous dimension & $+1$ & $\sum_{(a,b)}(\gamma_s)_{bd}C_{ijad}$ \\
Fermion anomalous dimension & $+1$ & $\sum_{(i,j)}(\gamma_f)_{jk}C_{ikab}$ \\
\bottomrule
\end{longtable}

\noindent In the scalar-pair rows, $\sum_{(a,b)}$ means both permutations of
external scalar indices, as implemented by \texttt{\_swap\_pair}. The mixed
Yukawa/gauge structures sum over the four independent scalar and flavour
pair exchanges. The crossed Yukawa structure includes the two simultaneous
exchanges. The scalar/fermion anomalous-dimension rows contribute only when
those tensors are provided to the RGE evaluator.

\textbf{Loop-order caveat.} Evolving $C_{12}$ at one loop and subsequently
inserting that correction in the one-loop $C_{12}\to C_5$ graph gives a
higher-order effect. The current fixed-order final-$C_5$ numerical adapter
keeps the direct one-loop $C_{12}\to C_5$ contribution, not this additional
self-running insertion.
'''


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,default=DEFAULT)
    ap.add_argument('--compile',action='store_true')
    args=ap.parse_args()
    original=args.input.read_text(encoding='utf-8')
    marker=r'\end{document}'
    if original.count(marker)!=1: raise ValueError('Expected one end-of-document marker')
    if 'EFT1 $C_{12}$ / $C_{LLSS}$ tensor RGE:' in original:
        raise ValueError('C12 section already exists; refusing duplicate')
    dest=args.input.with_name(args.input.stem+'_with_C12.tex')
    dest.write_text(original.replace(marker,SECTION+'\n'+marker),encoding='utf-8')
    print('Saved',dest)
    if args.compile:
        for _ in range(2):
            subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error',dest.name],cwd=dest.parent,check=True,stdout=subprocess.DEVNULL)
        print('PDF',dest.with_suffix('.pdf'))

if __name__=='__main__':
    main()
