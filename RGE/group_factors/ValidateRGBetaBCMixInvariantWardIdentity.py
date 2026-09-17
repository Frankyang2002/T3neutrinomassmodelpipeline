from __future__ import annotations

"""Ward-identity test for the T3-B/C RGBeta mixing invariant.

RGBeta defines, for dS1=dS2=2,

    I(h1,h2,s1,s2bar)
      = delS2[a,h1,h2] delS2[a,s1,s2bar]

and RefineGroupStructures gives the formal component tensor

    K = 1/2 (delta[h1,s2bar] delta[h2,s1]
           + delta[h1,s1]    delta[h2,s2bar]).

The physical field content is {H,H,S1,Bar[S2]}, so the infinitesimal SU(2)
action is

    +T on H,
    +T on H,
    +T on S1,
    -T^T on Bar[S2].

Because SU(2) doublets are pseudoreal, the formal RGBeta index on the S1 leg
must be mapped into the ordinary complex-component convention with epsilon.
This is exactly the same leg identified independently by the scalar-loop
closure diagnostic.

This script compares:
  1. the raw RGBeta-refined formal delta tensor;
  2. epsilon inserted on S1;
  3. epsilon inserted on Bar[S2];
  4. epsilon inserted on both scalar legs.

A correct ordinary-component embedding should satisfy the Ward identity
component-by-component.
"""

import sympy as sp


def generators():
    return (
        sp.Matrix([[0, 1], [1, 0]]) / 2,
        sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2,
        sp.Matrix([[1, 0], [0, -1]]) / 2,
    )


EPS = sp.Matrix([[0, 1], [-1, 0]])


def raw_tensor(h1, h2, s1, sb2):
    return sp.Rational(1, 2) * (
        int(h1 == sb2 and h2 == s1)
        + int(h1 == s1 and h2 == sb2)
    )


def eps_on_s1(h1, h2, s1, sb2):
    return sp.simplify(sum(
        raw_tensor(h1, h2, u, sb2) * EPS[u, s1]
        for u in range(2)
    ))


def eps_on_s2bar(h1, h2, s1, sb2):
    return sp.simplify(sum(
        raw_tensor(h1, h2, s1, u) * EPS[u, sb2]
        for u in range(2)
    ))


def eps_on_both(h1, h2, s1, sb2):
    return sp.simplify(sum(
        raw_tensor(h1, h2, u, v) * EPS[u, s1] * EPS[v, sb2]
        for u in range(2)
        for v in range(2)
    ))


def ward_component(tensor, T, h1, h2, s1, sb2):
    out = sp.S.Zero

    for p in range(2):
        out += T[h1, p] * tensor(p, h2, s1, sb2)
        out += T[h2, p] * tensor(h1, p, s1, sb2)
        out += T[s1, p] * tensor(h1, h2, p, sb2)

        # Bar[S2] transforms in the anti-fundamental:
        out -= T[p, sb2] * tensor(h1, h2, s1, p)

    return sp.simplify(out)


def residuals(tensor):
    rows = []
    for A, T in enumerate(generators(), start=1):
        nonzero = {}
        for h1 in range(2):
            for h2 in range(2):
                for s1 in range(2):
                    for sb2 in range(2):
                        r = ward_component(tensor, T, h1, h2, s1, sb2)
                        if r != 0:
                            nonzero[(h1, h2, s1, sb2)] = r
        rows.append((A, nonzero))
    return rows


def print_result(name, rows):
    print(name)
    total = 0
    for A, nonzero in rows:
        total += len(nonzero)
        print(f"  generator A={A}: nonzero Ward components = {len(nonzero)}")
        for key, value in list(nonzero.items())[:6]:
            print(f"    {key}: {value}")
        if len(nonzero) > 6:
            print("    ...")
    print(f"  total nonzero Ward components = {total}")
    print(f"  gauge invariant = {total == 0}")
    print()


def main():
    print("T3-B/C RGBeta mixing-invariant Ward-identity test")
    print("field ordering: H, H, S1, Bar[S2]")
    print()

    cases = [
        ("Raw RGBeta-refined tensor", raw_tensor),
        ("Epsilon on S1", eps_on_s1),
        ("Epsilon on Bar[S2]", eps_on_s2bar),
        ("Epsilon on both scalar legs", eps_on_both),
    ]

    totals = {}
    for name, tensor in cases:
        rows = residuals(tensor)
        print_result(name, rows)
        totals[name] = sum(len(x) for _, x in rows)

    ok = (
        totals["Raw RGBeta-refined tensor"] != 0
        and totals["Epsilon on S1"] == 0
        and totals["Epsilon on Bar[S2]"] != 0
    )

    if ok:
        print(
            "CONCLUSION: in an ordinary complex-component basis the RGBeta "
            "formal tensor becomes gauge invariant only when the SU(2) "
            "pseudoreality intertwiner is placed on the S1 leg. This matches "
            "the independent scalar-loop closure diagnostic."
        )
        return 0

    print("CONCLUSION: expected pseudoreal-map pattern was not obtained.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
