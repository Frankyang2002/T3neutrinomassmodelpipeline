"""Read-only, source-exact reduction of six native Matchete extended C5 exports.

The long InputForm expressions contain an *unindexed* MF introduced by
CanonicalWeinbergCoefficient as display metadata and an *indexed* MF[r].
Their equality is an EXPLICIT assumption made generation by generation.
It is not an independently certified cross-package mass dictionary.

No RGBeta -> Matchete coupling map, physical threshold matching, or full
scalar mass mixing is inferred by this module. Never use the result as a
certified neutrino-mass prediction.
"""
from __future__ import annotations

import ast
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping
import zipfile

import numpy as np
import sympy as sp
from scipy.integrate import quad

MODELS = {
    "A_p0": ("y1Conj", "lambdaT3ExtraQ021010I1", -2 / np.sqrt(3)),
    "A_m2": ("y2Conj", "lambdaT3ExtraQ020101I1", -2 / np.sqrt(3)),
    "D_p0": ("y1Conj", "lambdaT3ExtraQ021010I1", -2 / np.sqrt(3)),
    "D_m2": ("y2Conj", "lambdaT3ExtraQ020101I1", -2 / np.sqrt(3)),
    "E_p0": ("y1Conj", "lambdaT3ExtraQ021010I1", 2 * np.sqrt(2 / 3)),
    "E_m2": ("y2Conj", "lambdaT3ExtraQ020101I1", 2 * np.sqrt(2 / 3)),
}
MODEL_COEFFICIENT = {k: (-2 / sp.sqrt(3) if k[0] in 'AD' else 2 * sp.sqrt(sp.Rational(2, 3))) for k in MODELS}
SECTORS = ("full_C5", "zero_extra_C5", "extra_difference_C5")
M0, MR, MS1, MS2 = sp.symbols("M0 MR MS1 MS2", positive=True)
HBAR, K, L0, LE, YA, YB = sp.symbols("HBAR K L0 LE YA YB")
SYMBOLS = {str(v): v for v in (M0, MR, MS1, MS2, HBAR, K, L0, LE, YA, YB)}
YUKAWA_PATTERN = re.compile(r"Bar\[Coupling\[(y1|y2|y1Conj|y2Conj), \{Index\[d\$\$([12]), Flavor\], Index\[d\$\$1, NFlavor\]\}, 0\]\]")
QUARTIC_PATTERN = re.compile(r"Bar\[Coupling\[(lambdaT3ExtraQ[A-Za-z0-9]+), \{\}, 0\]\]")
INDEXED_MF = "Coupling[MF, {Index[d$$1, NFlavor]}, 0]"
GLOBAL_MF = "Coupling[MF, {}, 0]"


class C5SourceError(ValueError):
    pass


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise C5SourceError(reason)


def _safe_tree(node: ast.AST) -> sp.Expr:
    """Evaluate the restricted arithmetic grammar, NEVER exec/eval/sympify text."""
    if isinstance(node, ast.Expression):
        return _safe_tree(node.body)
    if isinstance(node, ast.Constant) and type(node.value) in (int, float):
        return sp.Rational(str(node.value))
    if isinstance(node, ast.Name) and node.id in SYMBOLS:
        return SYMBOLS[node.id]
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        value = _safe_tree(node.operand)
        return -value if isinstance(node.op, ast.USub) else value
    if isinstance(node, ast.BinOp):
        left, right = _safe_tree(node.left), _safe_tree(node.right)
        if isinstance(node.op, ast.Add): return left + right
        if isinstance(node.op, ast.Sub): return left - right
        if isinstance(node.op, ast.Mult): return left * right
        if isinstance(node.op, ast.Div): return left / right
        if isinstance(node.op, ast.Pow):
            require(right.is_number and abs(float(right)) <= 16 and right.is_integer,
                    "Unexpected exponent in native expression")
            return left ** right
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and not node.keywords:
        args = [_safe_tree(a) for a in node.args]
        if node.func.id == 'Sqrt' and len(args) == 1:
            return sp.sqrt(args[0])
        if node.func.id == 'Log' and len(args) == 1:
            return sp.log(args[0])
        if node.func.id == 'T3LoopI' and len(args) == 3:
            require(all(sp.simplify(a-b) == 0 for a, b in zip(args, (M0**2, MS1**2, MS2**2))),
                    "Canonical loop arguments have changed")
            return K
    raise C5SourceError("Unsupported native InputForm node: " + ast.dump(node)[:200])


def parse_native(text: str, model: str, sector: str) -> sp.Expr:
    """Strictly map the exact six-model InputForm vocabulary into symbols."""
    require(model in MODELS and sector in SECTORS, "Unsupported model/sector")
    extra, quartic, _ = MODELS[model]
    require(len(text) < 150000 and 1000 < len(text), "Unexpected native C5 expression length")
    seen = YUKAWA_PATTERN.findall(text)
    expected = (("y1", "1"), ("y2", "2")) if sector == 'zero_extra_C5' else (
        (extra, "1"), ("y2", "2")) if extra == 'y1Conj' else (("y1", "1"), (extra, "2"))
    if sector == 'full_C5':
        expected_counts = ({("y1", "1"):16,("y2", "2"):32,(extra, "1"):16}
                           if extra=='y1Conj' else {("y1", "1"):32,("y2", "2"):16,(extra, "2"):16})
    else:
        expected_counts = {expected[0]:16,expected[1]:16}
    from collections import Counter
    require(Counter(seen)==expected_counts, f"{model}/{sector}: Yukawa names/flavour slots changed: {Counter(seen)}")
    # All generic flavour products have the same ordered external (p,q) index structure.
    text = YUKAWA_PATTERN.sub(lambda m: 'YA' if m.group(2) == '1' else 'YB', text)
    text = QUARTIC_PATTERN.sub(lambda m: 'LE' if m.group(1) == quartic else 'INVALID_EXTRA_QUARTIC', text)
    if sector == 'zero_extra_C5':
        require('LE' not in text, "Zero sector contains extra quartic")
    else:
        require('LE' in text, "Missing requested extra quartic")
    text = text.replace('Coupling[lambdaT3, {}, 0]', 'L0')
    text = text.replace(INDEXED_MF, 'MR').replace(GLOBAL_MF, 'M0')
    text = text.replace('Coupling[MS1, {}, 0]', 'MS1').replace('Coupling[MS2, {}, 0]', 'MS2')
    text = re.sub(r'\bhbar\b', 'HBAR', text)
    text = text.replace('^','**').replace('[','(').replace(']',')')
    try:
        tree=ast.parse(text, mode='eval')
    except SyntaxError as exc:
        raise C5SourceError("Malformed native expression") from exc
    return _safe_tree(tree)


def _read_source(container: Path, name: str) -> bytes:
    if container.is_file() and container.suffix.lower() == '.zip':
        with zipfile.ZipFile(container) as z:
            require(name in z.namelist(), f"Missing archive entry: {name}")
            return z.read(name)
    path = container / name
    require(path.is_file(), f"Missing native export: {path}")
    return path.read_bytes()


def verify_model(container: Path, model: str) -> dict[str,Any]:
    require(model in MODELS, "Unsupported model")
    meta_name = f'{model}_h3.json'
    meta = json.loads(_read_source(container, meta_name).decode('utf-8-sig'))
    extra, quartic, _ = MODELS[model]
    require(isinstance(meta, dict) and meta.get('status')=='NativeOneStepExtendedWeinbergExtracted'
            and meta.get('model')==model and meta.get('heavy_flavours')==3
            and meta.get('additional_yukawa')==extra
            and meta.get('extraction_status')=='Success'
            and meta.get('matching_kind')=='full-one-step-one-loop-EFTOrder5'
            and meta.get('native_matchete_lagrangian_check_passed') is True
            and meta.get('extra_dependent_sector_diagnostic')=='ExplicitExtraYukawaDependence',
            f"{model}: native extraction provenance or status invalid")
    for flag in ('absolute_matchete_rgbeta_dictionary_certified','hierarchical_matching_certified',
                 'all_allowed_scalar_quadratics_included','physical_neutrino_mass_prediction',
                 'full_coefficient_flavor_index_structure_independently_verified',
                 'cross_package_numeric_substitution_authorized'):
        require(meta.get(flag) is False, f"{model}: unexpected certification field {flag}")
    require(meta.get('quartic_extensions_included') == (6 if model.startswith('E') else 5),
            f'{model}: changed native quartic inventory')
    coeff = MODEL_COEFFICIENT[model]
    expr = {}
    digests = {}
    for sector in SECTORS:
        name=f'{model}_h3.{sector}.inputform.txt'
        payload=_read_source(container,name)
        expr[sector]=parse_native(payload.decode('utf-8-sig'), model, sector)
        digests[sector]=hashlib.sha256(payload).hexdigest()
    # Compare all three independently exported expressions, BEFORE assigning M0=MR.
    require(sp.simplify(expr['full_C5']-expr['zero_extra_C5']-expr['extra_difference_C5'])==0,
            f'{model}: native full != zero + extra')
    base=coeff*HBAR*L0*YA*YB*MR*K
    extra_formula=coeff*HBAR*LE*YA*YB*MR*K
    require(sp.simplify(expr['zero_extra_C5'].subs(M0,MR) - base)==0,
            f'{model}: baseline kernel reduction failed')
    require(sp.simplify(expr['extra_difference_C5'].subs(M0,MR) - extra_formula)==0,
            f'{model}: extra kernel reduction failed')
    require(sp.simplify(expr['full_C5'].subs(M0,MR) - base-extra_formula)==0,
            f'{model}: full kernel reduction failed')
    return {
        'model':model, 'status':'ExactNativeOneStepOrderedCoefficientReduced',
        'extra_yukawa':extra,'extra_quartic_matchete':quartic,
        'common_group_prefactor_exact':str(coeff),
        'baseline_quartic_factor':'lambdaT3',
        'extra_quartic_factor':f'conjugate({quartic})',
        'baseline_flavour_product':'conj(y1[p,r])*conj(y2[q,r])',
        'extra_flavour_product':('conj(y1Conj[p,r])*conj(y2[q,r])' if extra=='y1Conj' else 'conj(y1[p,r])*conj(y2Conj[q,r])'),
        'loop_kernel':'T3LoopI(MF[r]^2,MS1^2,MS2^2)',
        'unindexed_MF_equal_generation_MF_assumed':True,
        'exported_coefficient_is_ordered_pq_not_independently_symmetric':True,
        'symmetric_C5_convention':'C5[p,q]=A5[p,q]+A5[q,p]',
        'exact_full_equals_zero_plus_extra_before_MF_identification':True,
        'exact_reduction_of_all_three_sectors_after_identification':True,
        'source_sha256':digests,
        'physical_rgbeta_matchete_dictionary_certified':False,
        'hierarchical_threshold_matching_certified':False,
        'physical_neutrino_mass_prediction_certified':False,
    }


def three_mass_loop(a: float,b: float,c: float) -> float:
    """The project's T3ThreeMassLoopExpanded(a,b,c), with squared GeV masses.

    Stable near degeneracies using a positive one-dimensional integral.
    I = integral_0^infinity t/((t+a)(t+b)(t+c)) dt.
    """
    vals=np.array([a,b,c],float)
    require(bool(np.isfinite(vals).all() and (vals>0).all()), 'Squared masses must be finite positive')
    scale=float(max(vals))
    x,y,z=vals/scale
    val,_=quad(lambda t:t/((t+x)*(t+y)*(t+z)),0,np.inf,epsabs=1e-12,epsrel=1e-11,limit=200)
    return float(val/scale)


def _scalar(value: Any,name: str) -> complex:
    require(isinstance(value,dict) and set(value)=={'real','imag'}, f'{name}: expected real/imag')
    z=complex(float(value['real']), float(value['imag']))
    require(np.isfinite(z.real) and np.isfinite(z.imag),f'{name}: nonfinite')
    return z


def _matrix(value: Any,name: str) -> np.ndarray:
    require(isinstance(value,dict) and set(value)=={'real','imag'},f'{name}: expected real/imag matrices')
    ar=np.asarray(value['real'],dtype=float);ai=np.asarray(value['imag'],dtype=float)
    require(ar.shape==(3,3) and ai.shape==(3,3) and np.isfinite(ar).all() and np.isfinite(ai).all(),
            f'{name}: expected finite complex 3x3 matrix')
    return ar+1j*ai


def _matrix_json(a: np.ndarray) -> dict[str,Any]:
    return {'real':np.asarray(a.real,dtype=float).tolist(),'imag':np.asarray(a.imag,dtype=float).tolist()}


def evaluate_conditional_c5(model: str, input_data: Mapping[str,Any]) -> dict[str,Any]:
    """Only explicitly Matchete-coordinate numerical inputs are accepted.

    Computes a one-step symmetric C5, NOT hierarchical matching or complete
    physical neutrino masses. Does not read RGBeta benchmark configurations.
    """
    require(model in MODELS, 'Unsupported model')
    require(input_data.get('model')==model and
            input_data.get('basis_declaration')=='explicit-matchete-source-couplings-no-rgbeta-transport',
            'Model/basis mismatch or Matchete coordinate declaration missing')
    p=input_data.get('parameters')
    require(isinstance(p,dict), 'Missing parameters')
    extra_name,quartic,_=MODELS[model]
    masses=np.asarray(p['MF_gev'],dtype=float)
    require(masses.shape==(3,) and np.isfinite(masses).all() and (masses>0).all(),
            'MF_gev must be three positive finite heavy masses')
    ms1=float(p['MS1_gev']);ms2=float(p['MS2_gev']);
    require(np.isfinite([ms1,ms2]).all() and min(ms1,ms2)>0,'Invalid scalar masses')
    y1=_matrix(p['y1'],'y1');y2=_matrix(p['y2'],'y2');yextra=_matrix(p[extra_name],extra_name)
    lam=_scalar(p['lambdaT3'],'lambdaT3');lam_extra=_scalar(p[quartic],quartic)
    loops=np.array([float(m)*three_mass_loop(float(m*m),ms1*ms1,ms2*ms2) for m in masses])
    pref=complex(MODEL_COEFFICIENT[model].evalf())/(16*np.pi**2)
    # Corresponding native expressions contain bar(y1), bar(y2), and bar(extra).
    def ordered(a:np.ndarray,b:np.ndarray, quartic_z:complex) -> np.ndarray:
        return pref*quartic_z*(np.conjugate(a)*loops[np.newaxis,:])@np.conjugate(b).T
    baseline_ordered=ordered(y1,y2,lam)
    extra_ordered=ordered(yextra if extra_name=='y1Conj' else y1,
                          y2 if extra_name=='y1Conj' else yextra,np.conjugate(lam_extra))
    base=baseline_ordered+baseline_ordered.T
    delta=extra_ordered+extra_ordered.T
    total=base+delta
    require(np.allclose(total,total.T,rtol=1e-14,atol=1e-25),'C5 symmetry failure')
    vev=float(p.get('vev_gev',246.22))
    require(np.isfinite(vev) and vev>0,'Invalid Higgs vev')
    mnu=-0.5*vev**2*total
    result={
        'status':'IllustrativeMatcheteCoordinateOneStepC5_NotPhysicalPrediction',
        'model':model,'basis_declaration':input_data['basis_declaration'],
        'masses_input_gev':masses.tolist(),'scalar_masses_input_gev':[ms1,ms2],
        'loop_weights_gev_inverse':loops.tolist(),
        'C5_baseline_gev_inverse':_matrix_json(base),'C5_extra_gev_inverse':_matrix_json(delta),
        'C5_total_gev_inverse':_matrix_json(total),
        'mnu_total_gev':_matrix_json(mnu),
        'max_abs_deltaC5_over_max_abs_baselineC5':float(np.max(np.abs(delta))/max(float(np.max(np.abs(base))),1e-300)),
        'absolute_matchete_rgbeta_dictionary_certified':False,
        'hierarchical_threshold_matching_certified':False,
        'complete_scalar_mass_mixing_included':False,
        'neutrino_mass_physical_prediction_certified':False,
    }
    # PMNS is basis-dependent. Only compute observables after supplying Ye.
    if 'ye' in p:
        from physics.NeutrinoTrajectory import charged_lepton_mass_basis_matrix
        from physics.NeutrinoObservables import calculate_neutrino_observables
        ye=_matrix(p['ye'],'ye')
        rotated=charged_lepton_mass_basis_matrix(mnu,ye)
        obs=calculate_neutrino_observables(rotated,ordering=str(p.get('ordering','NO')))
        require(obs.takagi_residual<1e-6,'Failed Takagi residual')
        result['conditional_observables']={
            'masses_ev':obs.masses_ev.tolist(),
            'delta_m21_sq_ev2':obs.delta_m21_sq_ev2,
            'delta_m31_sq_ev2':obs.delta_m31_sq_ev2,
            'pmns_abs':obs.pmns_abs.tolist(),
            'takagi_residual':obs.takagi_residual,
            'ordering':obs.ordering,
        }
    return result


def self_test():
    assert abs(three_mass_loop(2.,2.,2.)-0.25)<2e-13
    assert abs(three_mass_loop(4.,1.,9.)-0.1246970336020119)<2e-13
    try: _safe_tree(ast.parse("__import__('os').system('echo wrong')",mode='eval'))
    except C5SourceError: pass
    else: raise AssertionError('Unsafe AST accepted')
    for model in MODELS:
        rng=np.random.default_rng(100)
        vals={n:{'real':(rng.random((3,3))*0.003).tolist(),'imag':np.zeros((3,3)).tolist()}
              for n in ('y1','y2',MODELS[model][0])}
        vals.update({'lambdaT3':{'real':0.12,'imag':0.01},MODELS[model][1]:{'real':0.05,'imag':0.008},
                     'MF_gev':[2000,2400,3100],'MS1_gev':600,'MS2_gev':900})
        data={'model':model,'basis_declaration':'explicit-matchete-source-couplings-no-rgbeta-transport','parameters':vals}
        out=evaluate_conditional_c5(model,data)
        assert 'conditional_observables' not in out
        assert not out['neutrino_mass_physical_prediction_certified']
        total=np.array(out['C5_total_gev_inverse']['real'])+1j*np.array(out['C5_total_gev_inverse']['imag'])
        assert np.allclose(total,total.T)
        vals[MODELS[model][0]]={'real':np.zeros((3,3)).tolist(),'imag':np.zeros((3,3)).tolist()}
        zero=evaluate_conditional_c5(model,data)
        assert np.max(np.abs(np.array(zero['C5_extra_gev_inverse']['real'])))==0
    print('PASS: all six conditional 3x3 C5 models; zero-extra controls; stable loop integral; malicious AST rejected')
