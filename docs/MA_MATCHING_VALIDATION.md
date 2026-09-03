# Step 3: Ma-to-Weinberg matching validation

## Scope

This audit validates the small-$\lambda_5$, common-threshold matching in
`match_scotogenic_c5` against both:

1. Merle and Platscher, JHEP 11 (2015) 148,
   [arXiv:1507.06314v3](https://arxiv.org/abs/1507.06314), Eqs. (5)--(6); and
2. the actual T3-B Matchete coefficient checked by `RegressionC5.wl`.

The project convention is

$$
m_\nu=-v^2C_5.
$$

## Paper expression

For nearly degenerate neutral inert scalars, the paper gives

$$
m_\nu=-\frac{v^2\lambda_5}{16\pi^2}
h^T\operatorname{diag}\!\left(\frac{f(M_i,m_0)}{M_i}\right)h,
$$

where

$$
f(M,m_0)=
\frac{M^4}{(M^2-m_0^2)^2}\ln\!\left(\frac{m_0^2}{M^2}\right)
+\frac{M^2}{M^2-m_0^2}.
$$

Therefore

$$
C_5^{\rm Ma}=\frac{\lambda_5}{16\pi^2}
h^T\operatorname{diag}\!\left(\frac{f(M_i,m_0)}{M_i}\right)h.
$$

## Matchete equal-scalar limit

After taking $M_{S1}=M_{S2}=m_0$, `RegressionC5.wl` verifies that the T3-B
loop kernel multiplying
$\hbar\,\lambda_{T3}\,\overline y_1\overline y_2$ is

$$
F_{\rm Matchete}(M,m_0)=
\frac{M}{m_0^2-M^2}
\left[
1-\frac{M^2}{m_0^2-M^2}\ln\!\left(\frac{m_0^2}{M^2}\right)
\right].
$$

Direct algebra gives

$$
F_{\rm Matchete}(M,m_0)=-\frac{f(M,m_0)}{M}.
$$

The bridge between the two implemented conventions is therefore

$$
\hbar=\frac{1}{16\pi^2},\qquad
\lambda_{T3}=-\lambda_5,\qquad
y_1=y_2=h^*.
$$

Since $\overline y_1=\overline y_2=h$, the two minus signs cancel and

$$
C_5^{\rm Matchete}=C_5^{\rm Ma}.
$$

## Fully degenerate limit

Writing $r=m_0^2/M^2$, the implemented weight is

$$
\frac{f}{M}=\frac{1}{M}
\left[\frac{\ln r}{(1-r)^2}+\frac{1}{1-r}\right].
$$

Its finite $r\to1$ limit is

$$
\frac{f}{M}\longrightarrow-\frac{1}{2M}.
$$

Consequently,

$$
C_5\longrightarrow
-\frac{\lambda_5}{32\pi^2M}h^Th.
$$

This agrees with the Matchete limit
$F_{\rm Matchete}\to+1/(2M)$ after applying $\lambda_{T3}=-\lambda_5$.
There is no additional factor of two.

## Flavour and basis conventions

For a complex symmetric Majorana matrix, the code uses a Takagi
factorisation

$$
U^TMU=D_M,
$$

and rotates the heavy-row Yukawa matrix as

$$
h_{\rm mass}=U^Th.
$$

The matched result is

$$
C_5=\frac{\lambda_5}{16\pi^2}
h_{\rm mass}^T\operatorname{diag}\!\left(\frac{f_i}{M_i}\right)
h_{\rm mass}.
$$

Numerical tests confirm that this expression:

1. is complex symmetric;
2. is invariant under unitary changes of the heavy-field basis; and
3. transforms as $C_5\to V^TC_5V$ under a lepton-doublet basis rotation
   $h\to hV$.

## Numerical stability

Near $r=1$, direct evaluation contains a cancellation between two singular
terms. The implementation uses

$$
\frac{f}{M}=\frac{1}{M}
\left[-\frac12+\frac{r-1}{3}-\frac{(r-1)^2}{4}+\cdots\right],
$$

which was checked on both sides of the degenerate point and across scalar- and
fermion-dominated mass regimes.

## Mass parameter used at matching

The code uses the unbroken-phase inert parameter `mEta2`. The physical
neutral-scalar average after electroweak breaking is

$$
m_0^2=m_\eta^2+v^2(\lambda_3+\lambda_4).
$$

Using `mEta2` is the appropriate leading common-threshold SMEFT matching
choice. Electroweak-suppressed threshold corrections of order $v^2/M^2$ are
not included and should not be confused with a sign or normalization error.

## Result

The implemented small-$\lambda_5$ Ma matching formula, sign, loop factor,
complex conjugation bridge, flavour orientation, and degenerate limit are
consistent with the T3-B Matchete result. No numerical formula correction was
required during Step 3; only the convention documentation was made explicit.
