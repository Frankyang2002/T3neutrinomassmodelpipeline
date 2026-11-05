"""Full-lepton-flavor direct LLSS -> Weinberg running for general D(-2).

The already-verified one-generation LLHH beta is lifted using the EXACT PL
Yukawa pairs in the physical real-scalar stage-1 Matchete seed.  No general
quartic or extra Yukawa is removed.  The output is the direct, fixed-order
one-loop contribution to C5; it does NOT claim a full flavor-dependent LLSS
self-RGE or the scalar hard threshold.

Power counting: C_LLSS(tree) ~ hbar^0.  Its direct mixing to C5 is hbar^1;
running C_LLSS at hbar^1 then matching via a scalar loop is hbar^2.

Run: python -X utf8 -m scripts.RunGeneralDMinus2FullFlavorRunning
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Any

import sympy as sp
from sympy.printing.mathematica import mathematica_code

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MATCHING = ROOT / 'output/general_sequential_pilot/D_m2_stage1'
YUKAWAS = {'y1', 'y2', 'y2Conj'}
EXPECTED_PAIRS = {('y1', 'y2'), ('y1', 'y2Conj')}
BARRED_YUKAWA = re.compile(
    r'Bar\[Coupling\[(?P<name>y1|y2Conj|y2),\s*'
    r'\{Index\[(?P<external>[^,\]]+),\s*Flavor\],\s*'
    r'Index\[(?P<internal>[^,\]]+),\s*NFlavor\]\},\s*0\]\]'
)


def load(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f'Missing required D(-2) output: {path}')
    result = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(result, dict):
        raise ValueError(f'Expected a JSON object: {path}')
    return result


def save(path: Path, obj: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def pl_pair_inventory(seed: dict[str, Any]) -> dict[str, Any]:
    """Independently audit matched Yukawa-flavor contractions before lifting."""
    if seed.get('Status') != 'Success' or seed.get('Model') != 'D_-2':
        raise ValueError('Not a successful general D(-2) real-basis seed')
    if seed.get('PhysicalScalarNames') != ['NewScalar1', 'T3S2R', 'T3S2I']:
        raise ValueError('Scalar-field basis differs from verified D(-2) EFT1')
    terms = seed.get('TreeWilsonTerms', [])
    if len(terms) != 8 or seed.get('TreeWilsonTermCount') != 8:
        raise ValueError('Expected exactly eight verified tree Wilson terms')
    left = [t for t in terms if t.get('Chirality') == 'PL']
    if len(left) != 4 or sum(t.get('Chirality') == 'PR' for t in terms) != 4:
        raise ValueError('Expected four PL and four PR tree terms')
    audit = []
    pairs = set()
    for term in left:
        raw = term['TermInputForm']
        matches = list(BARRED_YUKAWA.finditer(raw))
        if len(matches) != 2:
            raise ValueError(
                f"PL term {term['Index']}: expected two barred indexed Yukawas, "
                f'found {len(matches)}; inspect raw Matchete flavor orientation'
            )
        # Ensure that NO other Yukawa factor escapes the indexed parser.
        all_couplings = re.findall(r'(?:Bar\[)?Coupling\[(y1|y2|y2Conj)\s*,', raw)
        if len(all_couplings) != 2:
            raise ValueError(f"PL term {term['Index']}: unparsed Yukawa factor")
        names = [m.group('name') for m in matches]
        pair = tuple(sorted(names, key=lambda v: ('y1', 'y2', 'y2Conj').index(v)))
        if pair not in EXPECTED_PAIRS:
            raise ValueError(f"PL term {term['Index']}: unexpected Yukawa pair {pair}")
        internal = {m.group('internal') for m in matches}
        external = [m.group('external') for m in matches]
        if len(internal) != 1 or external[0] == external[1]:
            raise ValueError(f"PL term {term['Index']}: inconsistent heavy/external flavor indices")
        pairs.add(pair)
        audit.append({
            'index': term['Index'], 'pair': list(pair),
            'heavy_index_dummy': next(iter(internal)),
            'external_flavor_dummies': external,
        })
    if pairs != EXPECTED_PAIRS:
        raise ValueError(f'Missing full-flavour Yukawa pair(s): {EXPECTED_PAIRS-pairs}')
    return {'status': 'Success', 'pl_terms': len(left),
            'unique_yukawa_pairs': [list(p) for p in sorted(pairs)],
            'terms': audit, 'heavy_flavor_mass': 'common MF (not a general mass matrix)'}


def beta_pair_coefficients(beta_string: str) -> dict[tuple[str, str], sp.Expr]:
    """Get exact coefficients multiplying conjugate(y_a) conjugate(y_b)/MF.

    This is extracted from the independently calculated ONE-GENERATION beta,
    not imposed from a guessed group-theory normalisation.
    """
    MF = sp.Symbol('MF')
    y1, y2, y2c = sp.symbols('y1 y2 y2Conj')
    phase = sp.Symbol('T3S2Phase')
    sym = {'MF': MF, 'y1': y1, 'y2': y2, 'y2Conj': y2c,
           'T3S2Phase': phase, 'conjugate': sp.conjugate, 'sqrt': sp.sqrt,
           'I': sp.I}
    parsed = sp.sympify(beta_string, locals=sym)
    # Physical mass-basis constraint |eta| = 1. No assumption on its angle.
    parsed = parsed.subs(sp.conjugate(phase), 1/phase)
    parsed = sp.cancel(parsed)
    b1, b2, bc = sp.symbols('yb1 yb2 yb2Conj')
    parsed = sp.expand(MF * parsed).xreplace({
        sp.conjugate(y1): b1,
        sp.conjugate(y2): b2,
        sp.conjugate(y2c): bc,
    })
    poly = sp.Poly(parsed, b1, b2, bc, domain='EX')
    expected = {(1,1,0): ('y1','y2'), (1,0,1): ('y1','y2Conj')}
    if set(poly.monoms()) != set(expected):
        raise ValueError(f'Unexpected Yukawa monomials in one-generation beta: {poly.monoms()}')
    result = {expected[m]: sp.simplify(poly.coeff_monomial(m)) for m in expected}
    if any(v == 0 or v.has(b1, b2, bc, MF, phase) for v in result.values()):
        raise ValueError('Beta coefficients retain Yukawas, MF or neutral-scalar phase')
    allowed_couplings = {'lambdaT3', 'lambdaT3ExtraQ020101I1'}
    other_symbols = set().union(*(v.free_symbols for v in result.values()))
    if {x.name for x in other_symbols} - allowed_couplings:
        raise ValueError('Unexpected couplings in the independently exported D(-2) Weinberg beta')
    return result


def format_kernel(a: str, b: str) -> str:
    return (f'Sum_r (conjugate({a}[p,r])*conjugate({b}[q,r]) + '
            f'conjugate({b}[p,r])*conjugate({a}[q,r]))/2')


def wolfram_y(name: str, external: str) -> str:
    return f'Bar[Coupling[{name}, {{{external}, Index[r, NFlavor]}}, 0]]'


def wolfram_kernel(name: str, a: str, b: str) -> str:
    ap, bq = wolfram_y(a, 'p'), wolfram_y(b, 'q')
    bp, aq = wolfram_y(b, 'p'), wolfram_y(a, 'q')
    return (f'{name}[p_, q_] := Module[{{r = Unique["nfl$"]}},\n'
            f'  Sum[({ap}*{bq} + {bp}*{aq})/2, {{r, 1, NF}}]\n];')


def wl_coupling(value: sp.Expr) -> str:
    code = mathematica_code(value)
    for name in ('lambdaT3ExtraQ020101I1','lambdaT3'):
        code = re.sub(rf'\b{re.escape(name)}\b', f'Coupling[{name}, {{}}, 0]', code)
    # Sympy conjugate -> Mathematica Conjugate; Matchete uses Bar for couplings.
    code = re.sub(r'Conjugate\[(Coupling\[[^\n]*?, \{\}, 0\])\]', r'Bar[\1]', code)
    return code


def wolfram_insertion(coeff: dict[tuple[str,str], sp.Expr], mu_hi: str, mu_low: str) -> str:
    c12 = wl_coupling(coeff[('y1','y2')])
    c1c = wl_coupling(coeff[('y1','y2Conj')])
    # Coefficients for full-flavor K are /MF. Their one-generation reduction
    # K -> conjugate(ya)*conjugate(yb) exactly.
    if not re.fullmatch('[A-Za-z][A-Za-z0-9]*',mu_hi+mu_low):
        raise ValueError('Scale symbols must be simple Wolfram names')
    return f'''(* General D(-2) direct EFT1 LLSS -> LLHH contribution.
   Verified against all four Matchete PL seed terms and one-generation beta.
   Coefficient is beta_16pi2 * Log[muLow/muHigh] -- NO hbar factor here.
   The lower scalar matching must carry this pure-SM contribution SEPARATELY,
   not include it as a heavy operator in the new one-loop scalar Match.
   This file alone does NOT execute scalar matching.
*)
ClearAll[EFT1D2C12, EFT1D2C1Conj, EFT1D2WilsonRunningHeavyInsertion,
  EFT1D2WilsonRunningWeinbergCoefficient, EFT1D2WilsonRunningMetadata];
If[!KeyExistsQ[GetCouplings[], {mu_hi}],
  DefineCoupling[{mu_hi}, SelfConjugate -> True];
];
If[!KeyExistsQ[GetCouplings[], {mu_low}],
  DefineCoupling[{mu_low}, SelfConjugate -> True];
];
{wolfram_kernel('EFT1D2C12', 'y1', 'y2')}
{wolfram_kernel('EFT1D2C1Conj', 'y1', 'y2Conj')}
EFT1D2WilsonRunningHeavyInsertion = 0;
EFT1D2WilsonRunningWeinbergCoefficient[p_, q_] :=
  Log[Coupling[{mu_low}, {{}}, 0]/Coupling[{mu_hi}, {{}}, 0]] *
  (({c12})*EFT1D2C12[p, q] + ({c1c})*EFT1D2C1Conj[p, q]) /
  Coupling[MF, {{}}, 0];
EFT1D2WilsonRunningMetadata = <|
  "Status" -> "PreparedNotYetMatched",
  "Model" -> "D_-2",
  "DirectWeinbergOnlyAtOneLoop" -> True,
  "FullLeptonFlavor" -> True,
  "HeavyOperatorInsertion" -> 0,
  "HbarFactorIncluded" -> False,
  "PhysicalMassEigenvalueRelationsImposed" -> False,
  "ScalarHardThresholdPerformed" -> False
|>;
'''


def calculate(matching: Path, mu_hi: str = 'MF', mu_low: str = 'MuScalar') -> dict[str, Any]:
    source = matching / 'data/general_d_m2_eft1_rge/weinberg_mixing_beta.json'
    seed_path = matching / 'data/eft1_after_F_realbasis_wilson_tensor_seed.json'
    stage1 = matching / 'data/general_d_m2_stage1_summary.json'
    resultdir = matching / 'data/general_d_m2_full_flavor'
    summarypath = resultdir / 'pilot_summary.json'
    state: dict[str, Any] = {'status': 'Started', 'model': 'D_-2',
        'scope': 'fixed-one-loop full-lepton-flavour direct LLSS-to-Weinberg running',
        'full_flavor_direct_running_performed': False,
        'full_wilson_self_rge_performed': False,
        'scalar_hard_threshold_performed': False,
        'physical_scalar_mass_relations_imposed': False,
        'general_renormalisable_rgbeta_mapping_certified': False}
    save(summarypath, state)
    try:
        if load(stage1).get('status') != 'GeneralDMinus2RealMassBasisFermionThresholdVerified':
            raise ValueError('Physical D(-2) first threshold not verified')
        beta = load(source)
        val = beta.get('weinberg_subspace_validation', {})
        if (beta.get('status') != 'Success' or not val.get('matches_weinberg_subspace')
            or val.get('failure_count') != 0 or beta.get('checked_LLHH_components') != 30
            or beta.get('nonzero_LLHH_beta_components') != 10):
            raise ValueError('One-generation Weinberg beta has not passed the reference validation')
        pair_inventory = pl_pair_inventory(load(seed_path))
        coeffs = beta_pair_coefficients(val['beta_kappa_16pi2'])
        # Require exactly the two mixed seed pairs to be present in beta.
        if set(coeffs) != {tuple(pair) for pair in pair_inventory['unique_yukawa_pairs']}:
            raise ValueError('One-generation beta and Matchete pair inventories disagree')
        MF = sp.Symbol('MF', nonzero=True)
        y1,y2,yc=sp.symbols('y1 y2 y2Conj')
        ref = sp.sympify(val['beta_kappa_16pi2'], locals={
            'MF': MF,'y1': y1,'y2':y2,'y2Conj':yc,
            'T3S2Phase': sp.Symbol('T3S2Phase'), 'conjugate':sp.conjugate,
            'sqrt':sp.sqrt})
        ph = sp.Symbol('T3S2Phase')
        ref = ref.subs(sp.conjugate(ph), 1/ph)
        reduced = sum(coeffs[pair]*sp.conjugate({'y1': y1,'y2': y2,'y2Conj':yc}[pair[0]])*
            sp.conjugate({'y1': y1,'y2': y2,'y2Conj':yc}[pair[1]])/MF for pair in coeffs)
        if sp.simplify(ref-reduced) != 0:
            raise ValueError('Full flavor lift fails one-generation reduction')
        if mu_hi == mu_low:
            raise ValueError('Supply distinct symbolic scale names; equal-scale test is internal')
        beta_items = [
            {'pair': list(pair), 'coefficient': sp.sstr(value),
             'flavor_kernel': format_kernel(*pair),
             'beta_term': f'({sp.sstr(value)})/MF * K_{pair[0]}_{pair[1]}[p,q]'}
            for pair,value in coeffs.items()
        ]
        full_beta = ' + '.join(x['beta_term'] for x in beta_items)
        log_ratio = f'log({mu_low}/{mu_hi})'
        result = {
            'status': 'Success', 'model': 'D_-2',
            'source_one_generation_beta': val['beta_kappa_16pi2'],
            'phase_constraint_applied': 'T3S2Phase*conjugate(T3S2Phase)=1',
            'power_counting': '16*pi^2*d(kappa_W[p,q])/dlnmu = beta16[p,q]',
            'mass_convention': 'common MF for NFlavor; nondegenerate fermion masses not implemented',
            'lepton_flavor_indices': ['p','q'], 'internal_fermion_flavor': 'r',
            'independent_mixed_kernels': beta_items,
            'full_flavor_beta_tensor': full_beta,
            'mu_high': mu_hi, 'mu_low': mu_low, 'log_ratio': log_ratio,
            'direct_one_loop_correction': f'({full_beta})*({log_ratio})/(16*pi**2)',
            'one_generation_reduction_matches': True,
            'symmetric_under_external_flavor_exchange': True,
            'equal_scale_running_vanishes': True,
            'pl_seed_audit': pair_inventory,
            'real_scalar_hard_threshold_performed': False,
            'all_scalar_masses_equal_assumed': False,
            'physical_scalar_mass_relations_imposed': False,
        }
        save(resultdir/'full_flavor_direct_running.json', result)
        wlpath=resultdir/'full_flavor_weinberg_insertion.wl'
        wlpath.write_text(wolfram_insertion(coeffs,mu_hi,mu_low),encoding='utf-8')
        mass_info=load(matching/'scalar_mass_basis.json')
        if (mass_info.get('RealFields') != ['T3S2R','T3S2I']
            or mass_info.get('RealMassSquared') != 'MS2^2-2*Abs[kappa]'
            or mass_info.get('ImagMassSquared') != 'MS2^2+2*Abs[kappa]'):
            raise ValueError('Real scalar mass splitting differs from verified D(-2) first stage')
        wxfs = [
          'eft1_after_F_realbasis_tree.wxf',
          'eft1_after_F_realbasis_loop.wxf',
          'eft1_after_F_realbasis_transition_tree.wxf',
          'eft1_after_F_realbasis_transition_full.wxf',
        ]
        absent=[name for name in wxfs if not (matching/'data'/name).is_file()]
        if absent:
            raise FileNotFoundError('Missing exact physical EFT1 continuation WXF(s): '+', '.join(absent))
        preflight={
          'status':'FullFlavorDirectInsertionPrepared_ScalarHardThresholdNotRun',
          'input_eft1_continuation_files': [str(matching/'data'/name) for name in wxfs],
          'running_insertion_file':str(wlpath),
          'scalar_matching_integrate_out':['NewScalar1','T3S2R','T3S2I'],
          'scalar_masses_matchete':['MS1','T3S2RMass','T3S2IMass'],
          'real_mass_squared':'MS2^2-2*Abs[kappa]',
          'imag_mass_squared':'MS2^2+2*Abs[kappa]',
          'phase_constraint':'Abs[eta]=1; kappa/eta^2=Abs[kappa]',
          'physical_mass_relations_imposed':False,
          'scalar_hard_threshold_performed':False,
          'insertion_is_pure_SM_and_is_not_to_be_one_loop_matched':True,
          'insertion_coefficient_omits_hbar':True,
        }
        save(resultdir/'scalar_threshold_preflight.json',preflight)
        state.update(status='GeneralDMinus2FullFlavorDirectRunningPrepared',
                     full_flavor_direct_running_performed=True,
                     full_wilson_self_rge_performed=False,
                     scalar_hard_threshold_performed=False,
                     scalar_threshold_insertion_prepared=True,
                     scalar_threshold_preflight_path=str(resultdir/'scalar_threshold_preflight.json'),
                     full_flavor_running_path=str(resultdir/'full_flavor_direct_running.json'),
                     scalar_threshold_insertion_path=str(wlpath),
                     pl_seed_audit=pair_inventory)
    except Exception as exc:
        state.update(status='Failed',error=f'{type(exc).__name__}: {exc}')
    save(summarypath,state)
    return state


def main(argv: list[str] | None = None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matching-dir',type=Path,default=DEFAULT_MATCHING)
    parser.add_argument('--mu-high',default='MF')
    parser.add_argument('--mu-low',default='MuScalar')
    args=parser.parse_args(argv)
    result=calculate(args.matching_dir.resolve(),args.mu_high,args.mu_low)
    print(f"{result['status']} -> {args.matching_dir.resolve()/'data/general_d_m2_full_flavor/pilot_summary.json'}")
    if result['status']=='Failed':print(result.get('error'))
    return 0 if result['status']=='GeneralDMinus2FullFlavorDirectRunningPrepared' else 1


if __name__=='__main__':
    raise SystemExit(main())
