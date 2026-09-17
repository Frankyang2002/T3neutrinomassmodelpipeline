# T3 Yukawa beta-function group factors

For the canonically normalised SU(2) invariant tensor
\[
C_{aAB}:\qquad \mathbf 2\otimes \mathbf{d_F}\to\mathbf{d_S},
\qquad d_S=d_F\pm1,
\]
define \(d_>=\max(d_F,d_S)\).  Orthogonality and Schur's lemma imply
\[
\sum_{a,A}C_{aAB}C^*_{aAB'}=\frac{d_>}{d_S}\delta_{BB'},
\]
\[
\sum_{A,B}C_{aAB}C^*_{a'AB}=\frac{d_>}{2}\delta_{aa'},
\]
\[
\sum_{a,B}C_{aAB}C^*_{aA'B}=\frac{d_>}{d_F}\delta_{AA'}.
\]

Hence define
\[
G_S=\frac{d_>}{d_S},\qquad
G_L=\frac{d_>}{2},\qquad
G_F=\frac{d_>}{d_F}.
\]

For \(y_1\), the matrix structures in the one-loop beta function have
\[
A_{\rm trace}=G_{S_1},
\]
\[
A_{\rm self}=\frac12\left(G_{L,1}+G_{F,1}\right)
=\frac{d_{>,1}(d_F+2)}{4d_F},
\]
and, because \(y_1\) and \(y_2\) share the lepton-doublet leg,
\[
A_{\rm cross}=\frac12G_{L,2}
=\frac{d_{>,2}}{4}.
\]
The charged-lepton contribution from the same \(L\)-leg wave-function
renormalisation is
\[
A_e=\frac12.
\]

For the five \(d\le3\) T3 classes this predicts

| model | \(A_{\rm trace}\) | \(A_{\rm self}\) | \(A_{\rm cross}\) | \(A_e\) |
|---|---:|---:|---:|---:|
| A | \(2\) | \(1\) | \(3/4\) | \(1/2\) |
| B | \(1\) | \(3/2\) | \(1/2\) | \(1/2\) |
| C | \(3/2\) | \(5/4\) | \(3/4\) | \(1/2\) |
| D | \(1\) | \(3/2\) | \(1/2\) | \(1/2\) |
| E | \(1\) | \(3/2\) | \(3/4\) | \(1/2\) |

Thus the non-gauge part of the generic \(y_1\) RGE is
\[
16\pi^2\beta_{y_1}\supset
G_{S_1}\operatorname{Tr}(y_1y_1^\dagger)y_1
+\frac12(G_{L,1}+G_{F,1})y_1y_1^\dagger y_1
+\frac12G_{L,2}y_2y_2^\dagger y_1
+\frac12Y_eY_e^\dagger y_1.
\]
