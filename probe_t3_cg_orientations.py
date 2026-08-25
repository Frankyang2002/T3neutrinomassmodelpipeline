from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

# Matchete can represent conjugation either when constructing the algebraic
# invariant or when assigning that tensor to fields. Test both independently.
INVARIANT_VARIANTS = {
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

WOLFRAM_TEMPLATE = r'''
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


def run_variant(invariant_expr: str, cg_expr: str) -> tuple[bool, str]:
    """Test one InvariantTensors/DefineCG conjugation convention."""
    wolfram_code = (
        WOLFRAM_TEMPLATE.replace("ALG_REPS", invariant_expr)
        .replace("CG_REPS", cg_expr)
    )

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".wl", delete=False, encoding="utf-8"
    ) as handle:
        handle.write(wolfram_code)
        script_path = Path(handle.name)

    try:
        proc = subprocess.run(
            ["wolframscript", "-file", str(script_path)],
            capture_output=True,
            text=True,
            errors="replace",
            check=False,
        )
    finally:
        script_path.unlink(missing_ok=True)

    output = f"{proc.stdout or ''}\n{proc.stderr or ''}"
    success = proc.returncode == 0 and "DEFINECG_SUCCESS" in output
    diagnostics = [
        line
        for line in output.splitlines()
        if "DefineCG::" in line or "General::" in line or "$Aborted" in line
    ]
    return success, "\n".join(diagnostics).strip()


def main() -> int:
    print("Testing Matchete 2 x 4 x 3 CG orientation conventions")
    print("=" * 72)

    passing: list[tuple[str, str]] = []

    for invariant_name, invariant_expr in INVARIANT_VARIANTS.items():
        for cg_name, cg_expr in CG_VARIANTS.items():
            success, diagnostics = run_variant(invariant_expr, cg_expr)
            status = "PASS" if success else "FAIL"
            print(
                f"{status:4}  InvariantTensors={invariant_name:12}  "
                f"DefineCG={cg_name}"
            )
            if diagnostics and not success:
                print("      " + diagnostics.replace("\n", "\n      "))
            if success:
                passing.append((invariant_name, cg_name))

    print("\n" + "=" * 72)
    print("PASSING CONVENTIONS")
    print("=" * 72)

    if not passing:
        print("NONE")
        return 1

    for invariant_name, cg_name in passing:
        print(f"InvariantTensors={invariant_name}, DefineCG={cg_name}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
