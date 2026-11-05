"""Regression tests for real-S2 D(-2) fermion-threshold boundary audit.

Fixtures test Python audit decisions; they do not simulate Wolfram matching.
"""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import sympy as sp

from scripts.RunGeneralDMinus2Stage1 import (
    audit_stage1, extra_quartic_coupling,
)

Q_NAMES = [
    "T3ExtraQ_0_0_0_0_0_4_1",
    "T3ExtraQ_0_0_0_0_1_3_1",
    "T3ExtraQ_0_0_1_1_0_2_1",
    "T3ExtraQ_0_2_0_1_0_1_1",
    "T3ExtraQ_1_1_0_0_0_2_1",
]
M_NAME = "lambdaT3ExtraM000002I1"
Q_SYMBOLS = [extra_quartic_coupling(name) for name in Q_NAMES]


def write_json(root: Path, relative: str, payload: dict) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def fixture(root: Path) -> None:
    write_json(root, "general_uv_inventory.json", {
        "Model": "D_-2", "BuildStatus": "Success", "HeavyFieldsAtFirstThreshold": ["F"],
        "AllowedInteractions": ["Yukawa1_GenericCG", "Yukawa2_GenericCG",
                                *Q_NAMES, "T3ExtraM_0_0_0_0_0_2_1",
                                "Yukawa2ConjugateScalar_GenericCG"],
        "AdditionalQuartics": Q_NAMES,
        "AdditionalQuadratics": ["T3ExtraM_0_0_0_0_0_2_1"],
        "AdditionalYukawas": ["Yukawa2ConjugateScalar_GenericCG"],
        "AdditionalQuarticsCount": 5,
        "AdditionalQuadraticsCount": 1,
        "AdditionalYukawasCount": 1,
        "SymmetryCompletenessCertified": False,
        "GeneralFermionFlavorMassCertified": False,
    })
    write_json(root, "scalar_mass_basis.json", {
        "Model": "D_-2", "Status": "RealSingletMassBasisConstructedNotNumericallyCertified",
        "OriginalComplexField": "NewScalar2", "RealFields": ["T3S2R", "T3S2I"],
        "PhysicalFieldDefinition": "S2=eta*(R+i*J)/sqrt(2)",
        "PhaseConstraint": "Abs[eta]=1; kappa/eta^2=Abs[kappa]",
        "HolomorphicCoupling": M_NAME,
        "LagrangianBilinearConvention": "+kappa*(S2dag)^2 + h.c.",
        "RealMassSquared": "MS2^2-2*Abs[kappa]",
        "ImagMassSquared": "MS2^2+2*Abs[kappa]",
        "MassPositivityCondition": "MS2^2>2*Abs[kappa]",
        "MatcheteRealMassParameters": ["T3S2RMass", "T3S2IMass"],
        "RemovedConjugateBilinearTermCount": 2,
        "GeneralQuarticCouplingNames": Q_SYMBOLS,
        "PhysicalMassRelationsImposedInMatchete": False,
        "ScalarThresholdPerformed": False,
    })
    write_json(root, "first_threshold_status.json", {
        "Status": "Success", "Model": "D_-2",
        "ThresholdPlan": [["F"]], "IntegratedFields": ["F"],
        "ActiveFields": ["S1", "S2R", "S2I"],
        "FermionAbsentFromEFT1Tree": True,
        "Scalar1PresentInEFT1Tree": True,
        "RealScalarRPresentInEFT1Tree": True,
        "RealScalarIPresentInEFT1Tree": True,
        "OriginalComplexS2AbsentFromEFT1Tree": True,
        "GeneralExtraQuarticsRetained": True,
        "ExtraYukawaRetainedInEFT1": True,
        "TreeWilsonTermCount": 2,
        "ScalarMassBasisDiagonalized": True,
        "PhysicalMassRelationsImposed": False,
        "ScalarThresholdPerformed": False,
        "IntermediateWilsonRunningPerformed": False,
        "StandardComplexScalarTensorAdapterCompatible": False,
    })
    write_json(root, "data/eft1_after_F_realbasis_wilson_inventory.json", {
        "Format": "D_m2_realbasis_wilson_inventory_not_tensor_adapter_input",
        "Status": "Success", "TreeWilsonTermCount": 2,
        "TreeWilsonTermsInputForm": [
            "GammaCC[...]*Field[T3S2R,...]*Field[T3S2R,...]",
            "GammaCC[...]*Field[T3S2I,...]*Field[NewScalar1,...]",
        ],
        "PhysicalScalarNames": ["NewScalar1", "T3S2R", "T3S2I"],
        "StandardComplexScalarTensorAdapterCompatible": False,
    })
    physical = ("Field[NewFermion, Fermion]*Field[NewScalar1, Scalar]*"
                "Field[T3S2R, Scalar]*Field[T3S2I, Scalar]*Coupling[y2Conj,{},0]*"
                + "*".join("Coupling[" + name + ",{},0]" for name in Q_SYMBOLS))
    tree = physical.replace("Field[NewFermion, Fermion]*", "") + "*GammaCC[...]"
    (root / "data").mkdir(exist_ok=True)
    (root / "data/physical_uv_realbasis_inputform.txt").write_text(physical, encoding="utf-8")
    (root / "data/eft1_after_F_realbasis_tree_inputform.txt").write_text(tree, encoding="utf-8")
    (root / "data/eft1_after_F_realbasis_loop_inputform.txt").write_text("hbar*...", encoding="utf-8")
    for name in (
        "eft1_after_F_realbasis_tree.wxf",
        "eft1_after_F_realbasis_loop.wxf",
        "eft1_after_F_realbasis_transition_tree.wxf",
        "eft1_after_F_realbasis_transition_full.wxf",
    ):
        (root / "data" / name).write_bytes(b"fixture")


class TestGeneralDMinus2FirstThreshold(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        fixture(self.root)

    def test_valid_boundary(self):
        result = audit_stage1(self.root)
        self.assertEqual(result["status"], "GeneralDMinus2RealMassBasisFermionThresholdVerified")
        self.assertEqual(result["first_matching"]["tree_wilson_terms"], 2)
        self.assertFalse(result["second_scalar_threshold_performed"])

    def test_quartic_symbol_mapping(self):
        self.assertEqual(extra_quartic_coupling("T3ExtraQ_0_0_0_0_0_4_1"),
                         "lambdaT3ExtraQ000004I1")
        self.assertEqual(extra_quartic_coupling("T3ExtraQ_0_2_0_1_0_1_1"),
                         "lambdaT3ExtraQ020101I1")

    def test_analytic_mass_split_sign(self):
        r, j, kappa, m2 = sp.symbols("r j kappa m2", real=True)
        s2 = (r + sp.I*j) / sp.sqrt(2)
        lag = sp.expand(-m2 * sp.conjugate(s2)*s2
                        + kappa*sp.conjugate(s2)**2 + kappa*s2**2)
        self.assertEqual(sp.simplify(lag.coeff(r, 2)), -m2/2 + kappa)
        self.assertEqual(sp.simplify(lag.coeff(j, 2)), -m2/2 - kappa)

    def test_reject_lost_quartic(self):
        path = self.root / "data/eft1_after_F_realbasis_tree_inputform.txt"
        path.write_text(path.read_text().replace(Q_SYMBOLS[0], "removed"))
        with self.assertRaisesRegex(ValueError, "General quartic lost"):
            audit_stage1(self.root)

    def test_reject_unrotated_complex_s2(self):
        path = self.root / "data/physical_uv_realbasis_inputform.txt"
        path.write_text(path.read_text() + "*Field[NewScalar2, Scalar]")
        with self.assertRaisesRegex(ValueError, "Original complex scalar"):
            audit_stage1(self.root)

    def test_reject_residual_off_diagonal_mass_term(self):
        path = self.root / "data/physical_uv_realbasis_inputform.txt"
        path.write_text(path.read_text() + "*" + M_NAME)
        with self.assertRaisesRegex(ValueError, "Off-diagonal"):
            audit_stage1(self.root)

    def test_reject_wrong_threshold(self):
        path = self.root / "first_threshold_status.json"
        x = json.loads(path.read_text())
        x["IntegratedFields"] = ["F", "S1"]
        write_json(self.root, "first_threshold_status.json", x)
        with self.assertRaisesRegex(ValueError, "matching stage"):
            audit_stage1(self.root)

    def test_reject_legacy_adapter_claim(self):
        path = self.root / "data/eft1_after_F_realbasis_wilson_inventory.json"
        x = json.loads(path.read_text())
        x["StandardComplexScalarTensorAdapterCompatible"] = True
        write_json(self.root, "data/eft1_after_F_realbasis_wilson_inventory.json", x)
        with self.assertRaisesRegex(ValueError, "legacy compatible"):
            audit_stage1(self.root)

    def test_reject_scalar_matching_claim(self):
        path = self.root / "first_threshold_status.json"
        x = json.loads(path.read_text())
        x["ScalarThresholdPerformed"] = True
        write_json(self.root, "first_threshold_status.json", x)
        with self.assertRaisesRegex(ValueError, "Unsupported completion"):
            audit_stage1(self.root)


if __name__ == "__main__":
    unittest.main()
