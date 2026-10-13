# LLSS Fixed-Order Diagnostic

This note records the implementation convention for displaying the
dimension-five scalar-lepton Wilson coefficient between the fermion and scalar
thresholds.

## Implemented object

The existing component Wilson machinery in

```text
RGE/running/intermediate/ScalarOnlyWilsonFlow.py
```

implements the fixed-order leading-log transport

\[
C_{\rm LLSS}(\mu)
=
C_{\rm LLSS}^{(0)}
+
\frac{1}{16\pi^2}
\ln\!\left(\frac{\mu}{M_F}\right)
\beta^{(1)}[C_{\rm LLSS}^{(0)}]
+
\mathcal O(\hbar^2).
\]

`Numerical/diagnostics/IntermediateWilsonDiagnostics.py` samples this same
expression at requested intermediate scales. It does not define a new beta
function and does not solve an RG-improved Wilson ODE.

## Flavor limitation

The component Wilson transport is a symbolic one-generation/group-theory
object. The production numerical T3 calculation instead uses full complex
\(3\times3\) lepton-heavy Yukawa matrices.

The diagnostic therefore preserves symbolic component expressions. It does not
replace a Yukawa matrix by an arbitrary scalar, trace, matrix element, or norm
in order to manufacture a numerical \(C_{\rm LLSS}\) curve.

A numerical LLSS norm should only be added after a verified full-flavor lift of
the component Wilson RGE is implemented.

## Fixed-order power counting

The LLSS self-running correction is one-loop:

\[
\delta C_{\rm LLSS}^{(1)}=\mathcal O(\hbar).
\]

Converting that correction to the Weinberg operator through an additional
scalar loop gives

\[
\delta C_5
\sim
\hbar\,\delta C_{\rm LLSS}^{(1)}
=
\mathcal O(\hbar^2).
\]

It is therefore not inserted into the authoritative one-loop final \(C_5\).

The authoritative one-loop hierarchical result remains

\[
C_5(M_S)
=
C_{5,\rm hard}
+
\Delta C_{5,\rm direct}.
\]

The symbolic LLSS trajectory is diagnostic information only.

## Distinction from direct Weinberg running

The existing numerical intermediate plot

```text
intermediate_direct_weinberg_running.png
```

shows the full-flavor direct one-loop \(LLSS\to O_5\) contribution that actually
enters the authoritative final \(C_5\).

The symbolic LLSS fixed-order diagnostic is a different object and must not be
labelled as the physical Weinberg coefficient.
