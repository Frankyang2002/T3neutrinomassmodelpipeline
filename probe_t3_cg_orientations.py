
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

ALG_VARIANTS = {
    "plain_plain": "{{1}, {3}, {2}}",
    "C2_plain4": "{CRep[{1}], {3}, {2}}",
    "plain2_C4": "{{1}, CRep[{3}], {2}}",
    "C2_C4": "{CRep[{1}], CRep[{3}], {2}}",
}

CG_VARIANTS = {
    "plain_plain": "{SU2L[fund], T3Probe4, T3Probe3}",
    "bar2_plain4": "{Bar[SU2L[fund]], T3Probe4, T3Probe3}",
    "plain2_bar4": "{SU2L[fund], Bar[T3Probe4], T3Probe3}",
    "bar2_bar4": "{Bar[SU2L[fund]], Bar[T3Probe4], T3Probe3}",
}

WL_TEMPLATE = r'''
load = UsingFrontEnd[Needs["Matchete`"]; True];
If[load =!= True, Exit[10]];

ResetAll[];
LSM = LoadModel["SM"];

DefineRepresentation[T3Probe3, SU2L, {2}, IndexAlphabet -> {"a","b","c","d"}];
DefineRepresentation[T3Probe4, SU2L, {3}, IndexAlphabet -> {"u","v","w","x"}];

tensors = InvariantTensors[SU[2], ALG_REPS];
If[!ListQ[tensors] || Length[tensors] == 0, Exit[11]];

UsingFrontEnd[
  DefineCG[T3ProbeCG, CG_REPS, First[tensors]];
];

Print["DEFINECG_SUCCESS"];
Exit[0];
'''

def run_variant(alg_name: str, alg_expr: str, cg_name: str, cg_expr: str) -> tuple[bool, str]:
    text = WL_TEMPLATE.replace("ALG_REPS", alg_expr).replace("CG_REPS", cg_expr)

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".wl",
        delete=False,
        encoding="utf-8",
    ) as f:
        f.write(text)
        path = Path(f.name)

    try:
        proc = subprocess.run(
            ["wolframscript", "-file", str(path)],
            capture_output=True,
            text=True,
            errors="replace",
        )
    finally:
        path.unlink(missing_ok=True)

    combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
    ok = proc.returncode == 0 and "DEFINECG_SUCCESS" in combined
    relevant = [
        line for line in combined.splitlines()
        if "DefineCG::" in line or "General::" in line or "$Aborted" in line
    ]
    return ok, "\n".join(relevant).strip()


def main() -> int:
    print("Testing Matchete 2 x 4 x 3 CG orientation conventions")
    print("=" * 72)

    passes: list[tuple[str, str]] = []

    for alg_name, alg_expr in ALG_VARIANTS.items():
        for cg_name, cg_expr in CG_VARIANTS.items():
            ok, detail = run_variant(alg_name, alg_expr, cg_name, cg_expr)
            status = "PASS" if ok else "FAIL"
            print(f"{status:4}  InvariantTensors={alg_name:12}  DefineCG={cg_name}")
            if detail and not ok:
                print("      " + detail.replace("\n", "\n      "))
            if ok:
                passes.append((alg_name, cg_name))

    print("\n" + "=" * 72)
    print("PASSING CONVENTIONS")
    print("=" * 72)

    if not passes:
        print("NONE")
        return 1

    for alg_name, cg_name in passes:
        print(f"InvariantTensors={alg_name}, DefineCG={cg_name}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
