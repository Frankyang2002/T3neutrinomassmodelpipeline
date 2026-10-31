#!/usr/bin/env python3
"""Workstream A: fermionic census, native extra-quartic geometry, 16-model evidence.

Run from project root:
  python Lagrangian/RunT3UVFinalEvidence.py --self-test
  python Lagrangian/RunT3UVFinalEvidence.py --all

This is read-only for every production source and historical report. No Match,
no RGE, no EFT, no GitHub mutation. The output explicitly keeps unresolved
Matchete/RGBeta identification, symbolic eigenmass constraints, and discrete
symmetry registration as NOT certified. Native extra quartic tensors are
trusted only as source data, and all direct numerical identities are tested.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path
import hashlib
import json
import math
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Lagrangian.RunT3UVClosureStage import (  # noqa: E402
    MODELS, SELF, BC, DIMS, ROOT as CLOSURE_ROOT, SPECIES,
    audit_operator_basis, name_alpha, scalar_operators, prior_eigen_evidence,
    complex_rand, read_json, require,
)
from Numerical.experimental.GeneralT3MixedScalarMassBasis import (  # noqa: E402
    MixedDoubletMassBasis, evaluate_transported_tensor,
)
from Numerical.experimental.T3SelfBilinearMassBasis import parse_matchete_constant  # noqa: E402

STAGE22 = ROOT / 'output/general_t3_stage22/matchete_tensors.json'
PRIOR = ROOT / 'Reports/output/workstream_a/uv_closure_stage/uv_closure_stage.json'
OUT = ROOT / 'Reports/output/workstream_a/final_uv_evidence'
SOURCE = ROOT / 'Lagrangian/interactions/GeneralScalarOperators.wl'
TOPOLOGY = ROOT / 'Lagrangian/interactions/T3Topology.wl'
CATALOG = ROOT / 'Lagrangian/model/T3ModelCatalog.wl'
FIELDS = ROOT / 'Lagrangian/model/T3Fields.wl'

# Only the two topologically chosen Yukawa chiral contractions in the
# repository are assumed. The exhaustive exercise tests all gauge orientations
# for SM L/eR/Q/uR/dR with the existing F and S1/S2 fields.
BASELINE = {('Fc', 'S1'), ('F', 'S2dag')}
SM = {
    'Lbar': ('lepton', 2, Fraction(1, 2)),
    'ebar': ('lepton', 1, Fraction(1, 1)),
    'Qbar': ('color_triplet', 2, Fraction(-1, 6)),
    'ubar': ('color_triplet', 1, Fraction(-2, 3)),
    'dbar': ('color_triplet', 1, Fraction(1, 3)),
}


def source_sha(path: Path) -> dict[str, str]:
    return {'path': str(path.relative_to(ROOT)),
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def trilinear_su2_singlets(d1: int, d2: int, d3: int) -> int:
    """SU(2): multiplicity of spin 0 in j1 x j2 x j3 (0 or 1)."""
    lo = abs(d1 - d2) + 1
    hi = d1 + d2 - 1
    return int(lo <= d3 <= hi and ((d3 - lo) % 2 == 0))


def fermion_census(model: str) -> dict[str, Any]:
    family, alpha = name_alpha(model)
    d1, d2, df = DIMS[family]
    y1 = Fraction(alpha, 2)
    y2 = Fraction(alpha + 2, 2)
    yf = Fraction(alpha + 1, 2)
    majorana = (yf == 0 and df % 2 == 1)
    species = {'S1': (d1, y1), 'S1dag': (d1, -y1),
               'S2': (d2, y2), 'S2dag': (d2, -y2)}
    raw: list[dict[str, Any]] = []
    rejected_by_color = []
    for sm_name, (color, ds, y_smbar) in SM.items():
        if color != 'lepton':
            # F and S are SU(3)c singlets, so one colored SM fermion cannot
            # give a color-singlet trilinear. Record rather than silently skip.
            rejected_by_color.append(sm_name)
            continue
        for fkind, fy in (('F', yf), ('Fc', -yf)):
            for s, (dim, sy) in species.items():
                if y_smbar + fy + sy == 0 and trilinear_su2_singlets(ds, df, dim):
                    raw.append({'sm': sm_name, 'fermion': fkind, 'scalar': s,
                                'su2_singlet_multiplicity': 1})
    canonical: set[tuple[str, str, str]] = set()
    for record in raw:
        f = 'F' if majorana else record['fermion']
        canonical.add((record['sm'], f, record['scalar']))
    expected = {(('Lbar', 'F' if majorana else f), s) for f, s in BASELINE}
    expected = {(sm, f, s) for (sm, f), s in expected}
    extra: set[tuple[str, str, str]] = set()
    if alpha == -2:
        extra.add(('Lbar', 'F', 'S2'))
    if alpha == 0:
        extra.add(('Lbar', 'Fc', 'S1dag'))
    missing = (expected | extra) - canonical
    unexpected = canonical - (expected | extra)
    passed = not missing and not unexpected and len(canonical) == 2 + (alpha in (-2, 0))
    require(passed, f'{model}: fermionic Yukawa census discrepancy: {missing}, {unexpected}')
    # F F H/Hdag: jF x jF has integer isospin; jH=1/2 cannot
    # form a singlet for any of the current integer dF multiplets.
    # All BSM-odd linear mixing with SM fermions is excluded by exact Z2.
    ff_h_su2 = trilinear_su2_singlets(df, df, 2)
    require(ff_h_su2 == 0, 'Unexpected FFH SU2 singlet')
    return {
        'model': model, 'passed': passed, 'fermion_type': 'Majorana' if majorana else 'Dirac',
        'canonical_yukawa_operators': [list(x) for x in sorted(canonical)],
        'number_of_yukawa_coupling_matrices': len(canonical),
        'additional_to_topology': len(extra),
        'raw_gauge_orientations_before_majorana_identification': len(raw),
        'majorana_identifies_F_and_Fc': majorana,
        'charged_SM_fermion_chirality_tested': 'ebar',
        'colored_SM_fermions_excluded_by_color': rejected_by_color,
        'ff_higgs_su2_singlets': ff_h_su2,
        'bsm_sm_mass_mixing_forbidden_by_z2': True,
        'assumptions': 'Exact Z2; BSM SU3 singlets; modeled Lbar/F/Fc chirality; no new SM fermion families; renormalizable operators',
        'full_lorentz_flavour_tensor_normalization_certified': False,
    }


def stage22_model_key(model: str) -> str:
    f, alpha = name_alpha(model)
    return f'{f}_{alpha:+d}'


def parse_stage22(source: Path) -> dict[str, Any]:
    obj = read_json(source)
    require(obj.get('format') == 'stage22-exact-tensors-v1' and
            obj.get('source') == 'MatcheteInvariantTensors' and
            obj.get('model_count') == 16 and obj.get('quartic_slot_count') == 50,
            'Unexpected native Matchete Stage22 source schema')
    models = obj.get('models')
    require(isinstance(models, dict), 'Native Stage22 model map missing')
    require({stage22_model_key(m) for m in MODELS} == set(models),
            'Stage22 set of models does not equal ordinary 16-model inventory')
    return obj


def native_quartic_data(stage: dict[str, Any], model: str) -> list[dict[str, Any]]:
    source = stage['models'][stage22_model_key(model)].get('tensors')
    require(isinstance(source, list), f'{model}: Matchete native source has no tensors')
    expected = scalar_operators(model)[4]
    from Lagrangian.RunT3UVClosureStage import BASELINE_KEYS, canonical  # noqa: E402
    extra = {k: n for k, n in expected.items() if k not in BASELINE_KEYS}
    require(len(source) == sum(extra.values()), f'{model}: source slot count wrong')
    observed_counts: Counter[tuple[int, ...]] = Counter()
    observed_symbols: set[str] = set()
    tensors = []
    for record in source:
        counts = tuple(record['counts'])
        species = tuple(record['species'])
        dims = tuple(record['dimensions'])
        symbol = record['symbol']
        require(len(counts) == 6 and len(species) == len(dims) == 4 and
                counts == canonical(counts) and counts in extra,
                f'{model}: native quartic has unexpected species/counts')
        require(Counter(species) == Counter({s: n for s, n in zip(SPECIES, counts) if n}),
                f'{model}: field labels do not match monomial')
        require(symbol.startswith('lambdaT3ExtraQ') and symbol not in observed_symbols,
                f'{model}: duplicate/native quartic symbol')
        observed_symbols.add(symbol)
        require(record.get('potential_prefactor') == '1',
                f'{model}: native quartic potential prefactor is not 1')
        d1, d2, _ = DIMS[name_alpha(model)[0]]
        dim_dict = {'H': 2, 'Hdag': 2, 'S1': d1, 'S1dag': d1,
                    'S2': d2, 'S2dag': d2}
        require(all(dim_dict[s] == n for s, n in zip(species, dims)),
                f'{model}: native SU2 tensor axis dimensions do not match field')
        raw = record.get('components')
        require(isinstance(raw, list) and len(raw) == math.prod(dims),
                f'{model}: native component length wrong')
        tensor = np.array([parse_matchete_constant(x) for x in raw],
                          dtype=complex).reshape(dims)
        require(np.max(np.abs(tensor)) > 0, f'{model}: zero native invariant')
        observed_counts[counts] += 1
        tensors.append({'symbol': symbol, 'counts': counts, 'species': species,
                        'tensor': tensor, 'prefactor': 1})
    require(dict(observed_counts) == extra, f'{model}: native multiplicity mismatch')
    # Independently test native tensor basis rank for each field content.
    by_structure: dict[tuple[int, ...], list[np.ndarray]] = defaultdict(list)
    for item in tensors:
        by_structure[item['counts']].append(item['tensor'].reshape(-1))
    for structure, rows in by_structure.items():
        rank = int(np.linalg.matrix_rank(np.stack(rows), tol=1e-11))
        require(rank == len(rows), f'{model}: native invariant tensors linearly dependent {structure}')
    return tensors


def rel_error(a: complex | float, b: complex | float) -> float:
    return float(abs(a - b) / max(1.0, abs(a), abs(b)))


def contract(tensor: np.ndarray, names: tuple[str, ...],
             fields: dict[str, np.ndarray]) -> complex:
    result = tensor
    for name in names:
        result = np.tensordot(result, fields[name], axes=(0, 0))
    return complex(result)


def self_native_transport(model: str, source: list[dict[str, Any]],
                          *, seed: int = 310, trials: int = 5) -> dict[str, Any]:
    """Test actual stored Matchete CG components under native Cartesian S=(R+iI)p/sqrt2.

    This establishes algebraic transport of every EXTRA quartic tensor, not
    equality of Matchete and RGBeta native bases or the actual symbolic CG
    substitution inside the post-rotation .wl driver.
    """
    family, alpha = name_alpha(model)
    d1, d2, _ = DIMS[family]
    neutral = 'S2' if alpha == -2 else 'S1'
    dim = d2 if neutral == 'S2' else d1
    require(dim in (1, 3), 'Native self-bilinear irrep must be dimension 1 or 3')
    rng = np.random.default_rng(seed + 11 * dim)
    maxerr = 0.0
    for _ in range(trials):
        r, im = rng.normal(size=dim), rng.normal(size=dim)
        phase = np.exp(0.5j * rng.uniform(-math.pi, math.pi))
        sg = phase * (r + 1j * im) / math.sqrt(2)
        old_other = complex_rand(rng, (d1 if neutral == 'S2' else d2,))
        h = complex_rand(rng, (2,))
        fields: dict[str, np.ndarray] = {'H': h, 'Hdag': h.conj(),
            neutral: sg, neutral+'dag': sg.conj(),
            neutral+'_R': r, neutral+'_I': im}
        other_label = 'S1' if neutral == 'S2' else 'S2'
        fields[other_label] = old_other
        fields[other_label+'dag'] = old_other.conj()
        for term in source:
            names, t = term['species'], term['tensor']
            left = contract(t, names, fields)
            options = []
            for name in names:
                if name in (neutral, neutral+'dag'):
                    dag = name.endswith('dag')
                    factor = phase.conjugate() if dag else phase
                    options.append(((neutral+'_R', factor / math.sqrt(2)),
                                    (neutral+'_I', (-1j if dag else 1j) * factor / math.sqrt(2))))
                else:
                    options.append(((name, 1.0),))
            right = 0j
            for choice in product(*options):
                right += np.prod([factor for _, factor in choice]) * contract(
                    t, tuple(k for k, _ in choice), fields)
            maxerr = max(maxerr, rel_error(left, right))
            # PlusHc[mu*T] in source: compare the actual real potential term.
            mu = complex(rng.normal(), rng.normal())
            before_hc = 2 * (mu * left).real
            after_hc = 2 * (mu * right).real
            maxerr = max(maxerr, rel_error(before_hc, after_hc))
    require(maxerr < 2e-11, f'{model}: native tensor field rotation mismatch {maxerr}')
    return {'passed': True, 'actual_native_extra_tensors_checked': len(source),
            'trials_per_tensor': trials, 'max_relative_residual': maxerr,
            'native_basis': 'Matchete Cartesian neutral scalar; exact native component tensors',
            'all_source_plus_hc_couplings_have_prefactor_one': True}


def bc_native_transport(model: str, source: list[dict[str, Any]],
                        *, trials: int = 5) -> dict[str, Any]:
    rng = np.random.default_rng(729)
    basis = MixedDoubletMassBasis.build(model.replace('_m1', '_-1'),
                                       950**2, 1100**2, 80**2 * np.exp(0.53j))
    maxerr = 0.0
    for term in source:
        require(term['tensor'].shape == (2, 2, 2, 2),
                'B/C native tensor must have four doublet axes')
        for _ in range(trials):
            h = complex_rand(rng, (2,))
            x = complex_rand(rng, (2,))
            y = complex_rand(rng, (2,))
            left, right = evaluate_transported_tensor(term['tensor'], term['species'],
                                                      basis.U, x, y, h)
            maxerr = max(maxerr, rel_error(left, right))
            mu = complex(rng.normal(), rng.normal())
            maxerr = max(maxerr, rel_error(2*(mu*left).real, 2*(mu*right).real))
    require(maxerr < 2e-11, f'{model}: native B/C tensor rotation mismatch {maxerr}')
    return {'passed': True, 'actual_native_extra_tensors_checked': len(source),
            'trials_per_tensor': trials, 'max_relative_residual': maxerr,
            'native_basis': 'Matchete doublets transported with epsilon and saved SU2 structure',
            'all_source_plus_hc_couplings_have_prefactor_one': True}


def check_symbolic_quadratic_identities() -> dict[str, Any]:
    """Closed-form, parameter-general proofs of the two 2x2 mass spectra.

    Self: b = - kappa * lambda_extra in V, complex S = p(R+iI)/sqrt(2),
          p^2 = b/|b|; this makes b (S*)^2 + b* S^2 equal |b|(R²-I²).
    BC: eigenvalues of [[m1²,mu],[mu*,m2²]] via trace and determinant.
    These do not impose the relations on symbolic Matchete mass spurions.
    """
    import sympy as sp
    m1, m2, a, b, x = sp.symbols('m1 m2 a b x', real=True)
    rho = sp.symbols('rho', positive=True, real=True)
    z = a + sp.I*b
    A = sp.Matrix([[m1, z], [sp.conjugate(z), m2]])
    detpoly = sp.expand((A - x * sp.eye(2)).det())
    expected = sp.expand((m1-x)*(m2-x) - (a*a+b*b))
    require(sp.simplify(detpoly - expected) == 0, 'B/C exact characteristic polynomial failed')
    delta = sp.sqrt((m1-m2)**2 + 4*(a*a+b*b))
    eig0 = (m1+m2-delta)/2
    eig1 = (m1+m2+delta)/2
    require(sp.simplify(detpoly.subs(x, eig0)) == 0 and
            sp.simplify(detpoly.subs(x, eig1)) == 0,
            'B/C eigenvalue identity failed')
    # Let p be a unit-modulus phase satisfying b/p^2 = |b|.
    # The equality is algebraic without selecting a complex branch.
    R, I, msq = sp.symbols('R I msq', real=True)
    orig = msq * (R**2+I**2)/2 + rho*(R**2-I**2)
    diagonal = ((msq+2*rho)*R**2 + (msq-2*rho)*I**2)/2
    require(sp.expand(orig-diagonal) == 0, 'Self exact quadratic expansion failed')
    return {'passed': True,
            'bc_characteristic_polynomial': str(expected),
            'bc_squared_mass_eigenvalues': [str(eig0), str(eig1)],
            'self_original_potential_identity': str(orig),
            'self_eigen_potential_identity': str(diagonal),
            'self_conditions': 'rho=abs(-kappa*lambda_extra), phase^2=(-kappa*lambda_extra)/rho; nonzero rho (rho=0 arbitrary phase)',
            'implemented_as_constraints_on_matchete_symbols': False}


def validate_closure() -> dict[str, dict[str, Any]]:
    data = read_json(PRIOR)
    require(data.get('format') == 't3-workstream-a-uv-closure-stage-v1' and
            data.get('status') == 'StructuralAndNumericalStagePassed' and
            data.get('native_z2_pass_count') == 16 and
            data.get('prior_eigenbasis_pass_count') == 8 and
            data.get('scalar_basis_totals', {}).get('quartic_slots') == 50 and
            data.get('scalar_basis_totals', {}).get('bilinear_slots') == 8,
            'Previous 16-model UV closure certificate did not pass')
    rows = data.get('rows', [])
    require(len(rows) == 16 and {x['model'] for x in rows} == set(MODELS),
            'Previous closure model list differs')
    for row in rows:
        require(row['independent_scalar_basis'].get('passed') is True and
                row['z2_parity_of_actual_unrotated_bsm_terms'].get('passed') is True and
                row['previous_heavy_eigenbasis'].get('passed') is True,
                f"Prior closure incomplete for {row['model']}")
    return {row['model']: row for row in rows}


def run_all(stage22_path: Path, outdir: Path, overwrite: bool) -> dict[str, Any]:
    require(CLOSURE_ROOT == ROOT, 'Using incorrect closure stage repository')
    sources = [STAGE22, PRIOR, SOURCE, TOPOLOGY, CATALOG, FIELDS]
    for p in sources[1:]:
        require(p.is_file(), f'Missing required source/evidence: {p}')
    require(stage22_path.is_file(), 'Missing Stage22 Matchete tensor source: '+str(stage22_path))
    target_json = outdir/'final_uv_evidence.json'
    target_md = outdir/'final_uv_evidence.md'
    require(overwrite or (not target_json.exists() and not target_md.exists()),
            'Outputs already exist; back up first, then use --overwrite')
    closure = validate_closure()
    stage = parse_stage22(stage22_path)
    basis, totals = audit_operator_basis()
    require(totals['all_models_match_manifest'], 'Independent scalar counter failed')
    exact = check_symbolic_quadratic_identities()
    evidence = []
    for row in basis:
        model = row['model']
        require(model in closure and row['passed'], f'{model}: previous source basis invalid')
        census = fermion_census(model)
        source = native_quartic_data(stage, model)
        original_z2 = closure[model]['z2_parity_of_actual_unrotated_bsm_terms'].get('data')
        require(isinstance(original_z2, dict),
                f'{model}: actual built-interaction inventory missing from previous Z2 report')
        require(original_z2.get('extra_yukawas_included') == census['additional_to_topology'] and
                original_z2.get('extra_quartics_included') == len(source),
                f'{model}: independent fermionic/native quartic census disagrees with actual Matchete builder')
        native: dict[str, Any] = {'passed': True, 'actual_native_extra_tensors_checked': len(source),
                                  'source_independent_invariant_ranks_verified': True,
                                  'all_source_plus_hc_couplings_have_prefactor_one': True}
        if model in SELF:
            native.update(self_native_transport(model, source))
        elif model in BC:
            native.update(bc_native_transport(model, source))
        # Re-evaluate current prior eigenreport, not just copied summary status.
        eigen = prior_eigen_evidence(model)
        require(eigen.get('passed') is True, f'{model}: previous eigenbasis evidence unavailable')
        evidence.append({'model': model,
                         'scalar_basis_passed': True,
                         'actual_native_z2_parity_passed': True,
                         'actual_builder_yukawa_quartic_inventory_crosschecked': True,
                         'fermionic_census': census,
                         'native_extra_quartics': native,
                         'heavy_eigenbasis_reverified': bool(eigen['passed']),
                         'physical_matching_performed': False,
                         'explicit_z2_registry_in_matchete': False,
                         'matchete_rgbeta_operator_normalization_certified': False})
        print(f"{model}: Yukawa matrices {census['number_of_yukawa_coupling_matrices']} | "
              f"native extra quartic tensors {len(source)} | "
              f"{'native transport PASS' if model in SELF + BC else 'native rank PASS'}", flush=True)
    total_native = sum(r['native_extra_quartics']['actual_native_extra_tensors_checked'] for r in evidence)
    require(total_native == 50, 'Native Matchete tensor source does not contain 50 slots')
    require(sum(r['fermionic_census']['number_of_yukawa_coupling_matrices'] for r in evidence) == 38,
            'Total Yukawa matrices not 2*16+6')
    status = 'ConditionalUVOperatorEvidencePassedNotFullCompleteness'  # no certificate inflation
    result = {
        'format': 't3-workstream-a-final-uv-evidence-v1',
        'status': status,
        'source_provenance': [source_sha(p) for p in sources[1:]] +
                             [{'path':str(stage22_path), 'sha256':hashlib.sha256(stage22_path.read_bytes()).hexdigest()}],
        'reused_prior_16_model_closure': True,
        'independent_fermionic_yukawa_count': 38,
        'additional_fermionic_yukawa_matrices': 6,
        'matchete_native_extra_quartic_slots_checked': total_native,
        'matchete_native_extra_quartic_invariant_rank_checked': True,
        'symbolic_quadratic_parameter_identities': exact,
        'exceptional_native_quartic_tensor_transport_models': len(SELF+BC),
        'unresolved_for_full_certificate': [
            'A symbolic Matchete equality tying the eigenmass and phase spurions to original UV couplings',
            'Native source tensor transport for all baseline quartics and the actual Yukawa/topology CG tensors, including conventions',
            'Matchete-to-RGBeta field and coupling dictionary and prefactor normalization',
            'Explicit Z2 symmetry registration/selection rules in Matchete rather than only direct parity of assembled terms',
            'General flavour/CP field-redefinition conventions for Dirac and Majorana fermion families',
        ],
        'workstream_a_complete_uv_physical_equivalence_certified': False,
        'physical_eft_matching_performed': False,
        'rows': evidence,
    }
    outdir.mkdir(parents=True, exist_ok=True)
    target_json.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    lines = ['# Workstream A — final UV evidence stage', '',
             f'**Outcome: {status}** (conditional evidence; not full UV completeness).', '',
             f'- Native Matchete extra-quartic CG tensors: **{total_native}/50** checked, including independent ranks.',
             '- Fermionic lepton Yukawa inventory: **38 matrices** across the 16 models (six beyond the two/model baseline).',
             '- Actual unrotated Z2 parity: **16/16** (from prior passing Matchete reports).',
             '- Prior heavy scalar eigenbases: **8/8** passing (reports re-read).',
             '- Exceptionally mixed/self mass parameter relations: exact algebraic identities and native extra quartic tensor transport.',
             '', '| Model | Fermion | Yukawa matrices | Extra native quartics | Native rank | Native rotated |',
             '|---|---|---:|---:|---|---|']
    for r in evidence:
        f=r['fermionic_census']; n=r['native_extra_quartics']
        lines.append(f"| {r['model']} | {f['fermion_type']} | {f['number_of_yukawa_coupling_matrices']} | "
                     f"{n['actual_native_extra_tensors_checked']} | PASS | "
                     f"{'PASS' if r['model'] in SELF+BC else 'n/a'} |")
    lines += ['', '## Outstanding for full UV physical-equivalence certification', '']
    lines.extend('- '+x for x in result['unresolved_for_full_certificate'])
    lines += ['', 'This does **not** run Match[], extend RGEs, or establish threshold matching. It does not modify the GitHub repository.', '']
    target_md.write_text('\n'.join(lines), encoding='utf-8')
    return result


def self_test() -> None:
    require(len(MODELS) == 16 and len(SELF) == 6 and len(BC) == 2,
            'Expected 16/6/2 models')
    all_census = [fermion_census(m) for m in MODELS]
    require(sum(x['number_of_yukawa_coupling_matrices'] for x in all_census) == 38,
            'Yukawa enumeration total wrong')
    require(sum(x['additional_to_topology'] for x in all_census) == 6,
            'Exactly six additional Yukawas expected')
    require(sum(x['fermion_type'] == 'Majorana' for x in all_census) == 2,
            'Exactly two Majorana cases expected')
    require(all(x['ff_higgs_su2_singlets'] == 0 for x in all_census),
            'FFH invariance proof failed')
    check_symbolic_quadratic_identities()
    require(trilinear_su2_singlets(2,2,1) == 1 and
            trilinear_su2_singlets(2,2,3) == 1 and
            trilinear_su2_singlets(2,2,2) == 0, 'SU2 tensor product wrong')
    # Actual tensor transport with nontrivial complex CG and arbitrary complex
    # Yukawa/quartic couplings; does not require Stage22 to be mounted.
    rng=np.random.default_rng(131)
    for model in SELF:
        fam, a = name_alpha(model)
        d1,d2,_=DIMS[fam]
        dim = d2 if a == -2 else d1
        species = ('S2','S2dag','S2','S2dag') if a == -2 else ('S1','S1dag','S1','S1dag')
        t=complex_rand(rng,(dim,dim,dim,dim))
        self_native_transport(model,[{'species':species,'tensor':t}],trials=4)
    for model in BC:
        t=complex_rand(rng,(2,2,2,2))
        bc_native_transport(model,[{'species':('S1','S1dag','S2','S2dag'),'tensor':t}],trials=4)
    # Native contraction independence: two distinct invariants must have rank 2.
    e1=np.diag([1.,2.,3.,4.]).astype(complex).reshape(2,2,2,2)
    e2=np.eye(4).reshape(2,2,2,2).astype(complex)
    require(np.linalg.matrix_rank(np.stack([e1.reshape(-1),e2.reshape(-1)]))==2,
            'Native rank test bad')
    print('PASS: independent fermionic census 16/16 (38 Yukawas; 2 Majorana points)')
    print('PASS: analytic self/mixed squared-mass identities')
    print('PASS: native-coordinate complex tensor transport on eight exceptional models')
    print('NOTE: native Stage22 data and saved user-generated reports are inspected by --all')


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    opts=parser.add_mutually_exclusive_group(required=True)
    opts.add_argument('--self-test', action='store_true')
    opts.add_argument('--all', action='store_true')
    parser.add_argument('--stage22', type=Path, default=STAGE22)
    parser.add_argument('--output-dir', type=Path, default=OUT)
    parser.add_argument('--overwrite', action='store_true')
    args=parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    report=run_all(args.stage22.resolve(),args.output_dir.resolve(),args.overwrite)
    print('RESULT:',report['status'])
    print('JSON:', args.output_dir.resolve()/'final_uv_evidence.json')
    print('REPORT:', args.output_dir.resolve()/'final_uv_evidence.md')
    return 0

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError, ImportError, KeyError, TypeError) as exc:
        print('ERROR:',exc,file=sys.stderr)
        raise SystemExit(2)
