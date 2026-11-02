"""Static regressions for the registered Matchete portal export guard.

These checks deliberately do NOT impersonate executing wolframscript.  The
InputForm fragments are taken from the previously successful native restricted
operator export for A_+0 and B_-1.
"""
from __future__ import annotations

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
EXPORT = ROOT / "scripts" / "RunLODirectUVNativePortalTensorExport.wl"


def inputform_contract(text: str, cg: str, coupling: str) -> bool:
    """Python equivalent of the two literal InputForm token-membership tests."""
    return f"CG[{cg}," in text and f"Coupling[{coupling}," in text


class NativePortalExportGuardTests(unittest.TestCase):
    def test_a_p0_registered_singlet_portal_is_present(self) -> None:
        real_native_fragment = (
            "CG[T3HiggsPortal1CGInv1, {Index[d$$1, SU2L[fund]], "
            "Bar[Index[d$$2, SU2L[fund]]]}]*"
            "Coupling[lambdaH1Inv1, {}, 0]*Field[H, Scalar, "
            "{Index[d$$2, SU2L[fund]]}, {}]"
        )
        self.assertTrue(inputform_contract(
            real_native_fragment, "T3HiggsPortal1CGInv1", "lambdaH1Inv1"))

    def test_b_m1_registered_doublet_portal_is_present(self) -> None:
        real_native_fragment = (
            "CG[T3HiggsPortal2CGInv2, {Index[d$$1, SU2L[fund]], "
            "Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd2], "
            "Bar[Index[d$$2, T3BSMd2]]}]*Coupling[lambdaH2Inv2, {}, 0]"
        )
        self.assertTrue(inputform_contract(
            real_native_fragment, "T3HiggsPortal2CGInv2", "lambdaH2Inv2"))

    def test_fail_closed_wrong_symbol_and_prefix_collision(self) -> None:
        fragment = "CG[T3HiggsPortal1CGInv10, {}]*Coupling[lambdaH1Inv10, {}, 0]"
        self.assertFalse(inputform_contract(
            fragment, "T3HiggsPortal1CGInv1", "lambdaH1Inv1"))
        self.assertFalse(inputform_contract(
            fragment, "T3HiggsPortal1CGInv10", "lambdaH1Inv1"))
        self.assertTrue(inputform_contract(
            fragment, "T3HiggsPortal1CGInv10", "lambdaH1Inv10"))

    def test_full_replacement_uses_inputform_not_held_local_symbol(self) -> None:
        script = EXPORT.read_text(encoding="utf-8")
        executable = script[script.index("WPPortalAssemblyPresentQ["):]
        self.assertIn('StringContainsQ[form,"CG["<>cgName<>","]', executable)
        self.assertIn('StringContainsQ[form,"Coupling["<>couplingName<>","]', executable)
        self.assertNotIn("If[FreeQ[expr,HoldPattern[CG[sym,_List]]]", executable)
        self.assertIn("GetCGTensor[sym]", executable)
        self.assertIn('"registered_in_assembled_interaction"->True', executable)
        self.assertIn('"T3_EXPERIMENTAL_GENERAL_FERMIONS"', script)
        self.assertIn('"T3_EXPERIMENTAL_GENERAL_QUADRATICS"', script)
        self.assertNotIn("Match[", executable)


if __name__ == "__main__":
    unittest.main()
