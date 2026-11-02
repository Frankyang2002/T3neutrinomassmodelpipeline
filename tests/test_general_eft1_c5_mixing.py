"""Independent arithmetic, normalization and fail-closed EFT1 source tests."""
from __future__ import annotations

from itertools import combinations_with_replacement
import copy
import math
import unittest

import numpy as np

from RGE.running.intermediate.GeneralEFT1C5Mixing import (
    BoundaryError, CERTIFICATES, MODELS, LOOP, _collinear,
    _weinberg_template, calculate_master_beta, load_boundary,
    project_c5, solve_boundary,
)
from scripts.BuildSyntheticGeneralEFT1Fixtures import fixture


class SchemaTests(unittest.TestCase):
    def test_16_synthetic_complex_and_separate(self):
        self.assertEqual(len(MODELS), 16)
        for model in MODELS:
            result = solve_boundary(fixture(model), allow_synthetic=True)
            self.assertEqual(result["status"], "SyntheticDiagnostic_NotPhysicalC5")
            self.assertEqual(len(result["operator_results"]), 2)
            self.assertEqual(len(result["summed_by_rge_term_beta_C5_16pi2"]), 7)
            mat = result["total_beta_C5_16pi2"]
            beta = np.array(mat["real"]) + 1j*np.array(mat["imag"])
            self.assertEqual(beta.shape, (3, 3))
            self.assertGreater(np.max(abs(beta.imag)), 0)
            np.testing.assert_allclose(beta, beta.T, rtol=0, atol=1e-26)
            log = result["total_delta_C5_log_gev_inverse"]
            actual = np.array(log["real"]) + 1j*np.array(log["imag"])
            np.testing.assert_allclose(actual, beta * math.log(2000/100000) * LOOP,
                                       rtol=2e-14, atol=1e-28)
            self.assertFalse(result["physics_completeness"]["complete_one_loop_neutrino_mass"])

    def test_scalar_quartic_contraction_has_correct_four_factor(self):
        raw = fixture("A_-4")
        raw["operators"] = [raw["operators"][0]]
        raw["quartic_components"] = [{"indices": [3, 3, 5, 7], "value": 0.125}]
        result = solve_boundary(raw, allow_synthetic=True)
        beta = result["total_beta_C5_16pi2"]
        # C_nu_nu,S1,S2 = c for both scalar permutations;
        # beta_C5 = 4*lambda*(2c) = 8 lambda c.
        initial = complex(2.3e-11, 2.1e-12)
        actual = complex(beta["real"][0][0], beta["imag"][0][0])
        self.assertAlmostEqual(abs(actual / (8 * 0.125 * initial) - 1), 0, places=13)

    def test_missing_and_legacy_certifications_fail_closed(self):
        raw = fixture("A_-4")
        with self.assertRaisesRegex(BoundaryError, "Synthetic inputs require"):
            load_boundary(raw)
        legacy = dict(raw, format="t3-ordered-source-flavour-and-real-scalar-tensor-v1")
        with self.assertRaisesRegex(BoundaryError, "legacy ordered"):
            load_boundary(legacy, allow_synthetic=True)
        raw["authority"] = "certified"
        with self.assertRaisesRegex(BoundaryError, "missing certifications"):
            load_boundary(raw)
        raw["certificates"] = {name: True for name in CERTIFICATES}
        with self.assertRaisesRegex(BoundaryError, r"SU\(2\)"):
            load_boundary(raw)

    def test_scale_and_symmetry_validation(self):
        raw = fixture("A_-4")
        raw["mu_S_gev"] = raw["mu_F_gev"]
        with self.assertRaisesRegex(BoundaryError, "ordering"):
            load_boundary(raw, allow_synthetic=True)
        raw = fixture("A_-4")
        # Same canonical tensor component specified with incompatible value.
        raw["operators"][0]["components"].append(
            {"indices": [1, 1, 7, 5], "value": -1.0})
        with self.assertRaisesRegex(BoundaryError, "conflicting"):
            load_boundary(raw, allow_synthetic=True)

    def test_template_normalization_C5_over_four(self):
        c = load_boundary(fixture("A_-4"), allow_synthetic=True)
        mat = np.array([[1+1j, 2-3j, 0.5], [2-3j, 0, -1j],
                        [0.5, -1j, 1.7]], dtype=complex) * 1e-12
        tensor = _weinberg_template(mat, c)
        recovered, ward, symmetry = project_c5(tensor, c)
        np.testing.assert_allclose(recovered, mat, rtol=5e-16, atol=1e-27)
        self.assertLess(ward, 1e-26)
        self.assertLess(symmetry, 1e-26)
        self.assertAlmostEqual(abs(tensor[0,0,2,2]/(mat[0,0]/4) - 1), 0, places=14)


    def test_chat1_adapter_refuses_uncertified_mixed_boundary(self):
        from RGE.running.intermediate.CertifiedFermionBoundaryAdapter import (
            adapt_chat1_boundary, CHAT1_FORMAT, CONTEXT_FORMAT,
        )
        base = fixture("A_-4")
        keys = ("model", "mu_F_gev", "n_fermions", "n_real_scalars",
                "lepton_components", "higgs_real_components", "scalar_field_blocks")
        chat1 = {"format": CHAT1_FORMAT, **{k: base[k] for k in keys},
                 "canonical_operators": base["operators"],
                 "certificates": {k: False for k in CERTIFICATES}}
        eft1 = {"format": CONTEXT_FORMAT, **{k: base[k] for k in keys},
                "mu_S_gev": base["mu_S_gev"],
                "yukawa_components": base["yukawa_components"],
                "quartic_components": base["quartic_components"],
                "gauge_sectors": base["gauge_sectors"],
                "certificates": {k: False for k in CERTIFICATES}}
        with self.assertRaisesRegex(BoundaryError, "Only Chat 1"):
            adapt_chat1_boundary(chat1, eft1)
        chat1["certificates"] = {k: True for k in CERTIFICATES}
        eft1["certificates"] = {k: True for k in CERTIFICATES}
        with self.assertRaisesRegex(BoundaryError, "SU\\(2\\)"):
            adapt_chat1_boundary(chat1, eft1)
        eft1["n_real_scalars"] += 2
        with self.assertRaisesRegex(BoundaryError, "disagrees"):
            adapt_chat1_boundary(chat1, eft1)

    def test_cli_exclusive_outputs(self):
        from tempfile import TemporaryDirectory
        from pathlib import Path
        import json
        from scripts.RunGeneralEFT1C5Mixing import main as run
        with TemporaryDirectory() as td:
            root = Path(td)
            inp = root / "input" / "A_m4" / "boundary.json"
            inp.parent.mkdir(parents=True)
            inp.write_text(json.dumps(fixture("A_-4")), encoding="utf-8")
            args = ["--boundary-root", str(root / "input"),
                    "--output-root", str(root / "result"),
                    "--model", "A_-4", "--allow-synthetic"]
            self.assertEqual(run(args), 0)
            output = root / "result" / "A_m4" / "c5_eft1_log.json"
            old = output.read_bytes()
            with self.assertRaises(FileExistsError):
                run(args)
            self.assertEqual(output.read_bytes(), old)


class IndependentMasterRGETests(unittest.TestCase):
    def test_all_seven_against_direct_component_sums(self):
        """Check exact indexed formula on nontrivial complex synthetic data."""
        raw = fixture("B_+1")
        raw["yukawa_components"] = [
            {"indices":[1,2,3], "value":{"real":0.11,"imag":0.07}},
            {"indices":[2,3,5], "value":{"real":-0.12,"imag":0.02}},
            {"indices":[1,3,6], "value":{"real":0.04,"imag":-0.08}},
        ]
        raw["gauge_sectors"] = [{"coupling":0.47,"generators":[{
            "fermion": np.diag([1,-1,0,0,0,0]).tolist(),
            "scalar": np.diag([1,-1,0,0,0,0,0,0,0,0,0,0]).tolist(),
        }]}]
        c = load_boundary(raw, allow_synthetic=True)
        wilson = c.operators[0][1]
        fast = calculate_master_beta(c, wilson)
        y, lam, gauge = c.yukawa, c.quartic, c.gauge
        gs, gf = _collinear(c)
        i, j, a, b = 0, 2, 2, 2
        # Independent loop-by-loop reference for the seven terms.
        terms = dict.fromkeys(fast, 0j)
        for k in range(c.nf):
            for l in range(c.nf):
                for d in range(c.ns):
                    terms["yukawa_wavefunction"] += (
                        y[i,j,d]*y[k,l,d].conjugate()*wilson[k,l,a,b])
                    terms["conjugate_coefficient_yukawa"] += (
                        2*y[i,l,d]*y[j,k,d]+y[i,j,d]*y[k,l,d]
                    ) * wilson[k,l,a,b].conjugate()
        for aa, bb in ((a,b), (b,a)):
            for d in range(c.ns):
                terms["scalar_anomalous_dimension"] += gs[bb,d]*wilson[i,j,aa,d]
            for cc in range(c.ns):
                for dd in range(c.ns):
                    piece = -lam[aa,bb,cc,dd]
                    for g, gens in gauge:
                        for _, theta in gens:
                            piece += 2*g*g*theta[aa,cc]*theta[bb,dd]
                    terms["scalar_pair"] += -0.5 * piece*wilson[i,j,cc,dd]
        for ii, jj in ((i,j), (j,i)):
            for k in range(c.nf):
                terms["fermion_anomalous_dimension"] += gf[jj,k]*wilson[ii,k,a,b]
        for aa,bb,ii,jj in ((a,b,i,j),(b,a,i,j),(a,b,j,i),(b,a,j,i)):
            for l in range(c.nf):
                for cc in range(c.ns):
                    gauge_piece = 0j
                    for g, gens in gauge:
                        for t, theta in gens:
                            gauge_piece += 4*g*g*t[l,jj]*theta[bb,cc]
                    terms["mixed_yukawa_gauge"] += gauge_piece*wilson[ii,l,aa,cc]
                    for k in range(c.nf):
                        terms["mixed_yukawa_gauge"] += (
                            2*y[jj,k,cc]*y[k,l,bb].conjugate()
                            + y[jj,k,bb]*y[k,l,cc].conjugate()
                        )*wilson[ii,l,aa,cc]
        for aa,bb,ii,jj in ((a,b,i,j),(b,a,j,i)):
            for k in range(c.nf):
                for l in range(c.nf):
                    for cc in range(c.ns):
                        terms["crossed_yukawa"] += (
                            y[jj,l,aa]*y[ii,k,cc]*wilson[l,k,bb,cc])
        for name in fast:
            np.testing.assert_allclose(fast[name][i,j,a,b],terms[name],rtol=4e-12,atol=2e-28,
                                       err_msg=name)

    def test_block_diagonal_gauge_no_LLSS_to_LLHH(self):
        raw=fixture("A_-4")
        raw["quartic_components"]=[]
        f=np.diag([0.5,-0.5]*3)
        # A block-diagonal scalar generator cannot convert S1S2 into HH.
        s=np.diag([0.5,-0.5,0.5,-0.5]+[0.2]*(raw["n_real_scalars"]-4))
        raw["gauge_sectors"]=[{"coupling":0.65,"generators":[
            {"fermion":f.tolist(),"scalar":s.tolist()}]}]
        result=solve_boundary(raw,allow_synthetic=True)
        self.assertEqual(max(abs(z) for arr in result["total_beta_C5_16pi2"].values()
                             for row in arr for z in row),0)


class ScotogenicLimitTests(unittest.TestCase):
    def test_holomorphic_doublet_quartic_SU2_ward_and_sign(self):
        """Exact Ma lambda5 doublet contraction in this master tensor convention."""
        import sympy as sp
        raw = fixture("B_+1")
        raw["n_real_scalars"] = 10  # eta doublet + inert spectator real fields
        raw["scalar_field_blocks"] = [
            {"name":"H", "role":"Higgs", "indices":[1,2,3,4]},
            {"name":"eta", "role":"heavy", "indices":[5,6,7,8]},
            {"name":"spectator", "role":"heavy", "indices":[9,10]},
        ]
        raw["quartic_components"] = []
        raw["operators"] = []
        # H=(H+ , H0), eta=(eta+,eta0), interleaved Re,Im.
        phi=sp.symbols("x0:8",real=True)
        H=[(phi[k]+sp.I*phi[k+1])/sp.sqrt(2) for k in (0,2)]
        S=[(phi[k]+sp.I*phi[k+1])/sp.sqrt(2) for k in (4,6)]
        z=sum(sp.conjugate(H[k])*S[k] for k in range(2))
        lambda5=sp.Rational(1,9)
        potential=sp.expand(lambda5*(z**2+sp.conjugate(z)**2)/2)
        for inds in combinations_with_replacement(range(8),4):
            deriv=sp.diff(potential,*(phi[k] for k in inds))
            if deriv != 0:
                raw["quartic_components"].append({"indices":[k+1 for k in inds],
                                                   "value":float(sp.re(deriv))})
        # C^(0) = kappa (L.eta)^2, therefore canonical C5_eta=2*kappa.
        base = load_boundary(dict(raw, operators=[{"source_id":"placeholder", "components":[
            {"indices":[1,1,7,7],"value":1e-11}]}]), allow_synthetic=True)
        kappa=1.7e-11+2.1e-12j
        k_mat=np.array([[kappa,0.23*kappa,0.41*kappa],
                        [0.23*kappa,0.8*kappa,0.19*kappa],
                        [0.41*kappa,0.19*kappa,1.2*kappa]],complex)
        from dataclasses import replace
        eta_context=replace(base,higgs=(4,5,6,7))
        full=_weinberg_template(2*k_mat,eta_context)
        sparse=[]
        for i in range(6):
            for j in range(i,6):
                for a in range(4,8):
                    for b in range(a,8):
                        z=full[i,j,a,b]
                        if abs(z)>1e-25:
                            sparse.append({"indices":[i+1,j+1,a+1,b+1],
                                           "value":{"real":z.real,"imag":z.imag}})
        raw["operators"]=[{"source_id":"scotogenic_LL_eta_eta", "components":sparse}]
        context=load_boundary(raw,allow_synthetic=True)
        betas=calculate_master_beta(context,context.operators[0][1])
        total=sum(betas.values(),np.zeros_like(context.operators[0][1]))
        got,ward,sym=project_c5(total,context)
        np.testing.assert_allclose(got,4*float(lambda5)*k_mat,rtol=2e-13,atol=1e-26)
        self.assertLess(ward,1e-24,"holomorphic SU2 mixing must be Weinberg shaped")
        self.assertLess(sym,1e-24)
        # A genuinely normalized analytical toy (not a T3 UV prediction) may
        # pass the same strict SU(2)/U(1) and Weinberg Ward gate that is
        # required of a future physical handoff.  This tests the certified
        # execution path without inventing Matchete normalization factors.
        sigma = [np.array([[0,1],[1,0]],complex)/2,
                 np.array([[0,-1j],[1j,0]],complex)/2,
                 np.diag([0.5,-0.5]).astype(complex)]
        def real_rep(t):
            out=np.zeros((4,4),complex)
            for p in range(2):
                for q in range(2):
                    x=t[p,q]
                    out[2*p,2*q] += 1j*x.imag
                    out[2*p,2*q+1] += 1j*x.real
                    out[2*p+1,2*q] -= 1j*x.real
                    out[2*p+1,2*q+1] += 1j*x.imag
            return out
        su2=[]
        for t in sigma:
            f=np.zeros((6,6),complex)
            for p in range(3): f[2*p:2*p+2,2*p:2*p+2]=t
            scalar=np.zeros((10,10),complex)
            scalar[:4,:4]=real_rep(t)
            scalar[4:8,4:8]=real_rep(t)
            su2.append({"fermion":{"real":f.real.tolist(),"imag":f.imag.tolist()},
                        "scalar":{"real":scalar.real.tolist(),"imag":scalar.imag.tolist()}})
        hyper=np.zeros((10,10),complex)
        for start in (0,4): hyper[start:start+4,start:start+4]=real_rep(np.eye(2)*0.5)
        u1=[{"fermion":(-0.5*np.eye(6)).tolist(),
             "scalar":{"real":hyper.real.tolist(),"imag":hyper.imag.tolist()}}]
        physical_toy=copy.deepcopy(raw)
        physical_toy["authority"]="certified"
        physical_toy["certificates"]={k:True for k in CERTIFICATES}
        physical_toy["gauge_sectors"]=[
            {"name":"SU2L","coupling":0.65,"generators":su2},
            {"name":"U1Y","coupling":0.36,"generators":u1},
        ]
        certified_run=solve_boundary(physical_toy)
        self.assertTrue(certified_run["physics_completeness"]["physical_eft1_source_verified"])
        self.assertFalse(certified_run["physics_completeness"]["complete_one_loop_neutrino_mass"])


if __name__=="__main__":
    unittest.main()
