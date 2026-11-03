"""Offline tests for the actual original three-colour T3 RGE template adapter."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'scripts' / 'BuildT3OriginalColouredRGE.py'
spec = importlib.util.spec_from_file_location('original_colour_adaptor', SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def sample_template() -> str:
    """Local rendering stub using the source-verified original template anchors."""
    return '''<!doctype html><html><head><title>__TITLE__</title></head><body>
<header>Numerical integration direction: UV → IR. Saved benchmark trajectories; no matching or RGE recalculation.</header>
<button id="zoomUV">Full T3</button><button id="zoomInt">Intermediate EFT</button>
<button id="zoomFinal">SM + Weinberg</button><select id="yScaleMode"></select>
<select id="quantity">__OPTIONS__</select>
<script>
const DATA=__DATA__;
const META=__META__;
const COLORS=['#2466ae','#c46c25','#278864'];
const regionColours={uv:'#eaf2ff',int:'#fff1cf',sm:'#e5f5e9'};
function draw(){
 const q='c5_fnorm',weinberg=true;
 const showParts=weinberg&&showComponents.checked;
 document.getElementById('weinbergLegend').hidden=!weinberg;
 componentControl.hidden=!weinberg;
 rectRegion(s.mu_uv,s.mu_fermion_threshold,regionColours.uv,'Full T3');
 rectRegion(s.mu_fermion_threshold,s.mu_scalar_threshold,regionColours.int,'Intermediate scalar EFT');
 rectRegion(s.mu_scalar_threshold,s.mu_low,regionColours.sm,'SM + Weinberg EFT');
 if(false){}else if(weinberg){caption.textContent='old c5 continuation';}else{
  caption.textContent=DATA.models.filter(m=>enabled[m]&&traces[m]).length+' model(s) shown. Curves are taken from saved diagnostics.';
 }
}
function defaultLimits(q){
 const s=getScales();
 if(isLLSS(q)&&s.mu_fermion_threshold>s.mu_scalar_threshold)return [s.mu_fermion_threshold,s.mu_scalar_threshold];
 return [s.mu_uv,s.mu_low];
}
</script></body></html>'''


def make_stages() -> dict:
    return {
        'uv': {'mu_gev': [1e7, 1e6, 1e5],
               'series': {'g2': [0.6,0.59,0.58], 'lambdaT3_abs':[0.1,0.11,0.12]}},
        'eft1': {'mu_gev': [1e5, 1e4, 1e3],
                 'series': {'g2':[0.58,0.56,0.55], 'lambdaT3_abs':[0.12,0.13,0.14]}},
        'smeft': {'mu_gev': [1e3, 5e2, 1e2],
                  'series': {'g2':[0.55,0.54,0.53], 'c5_fnorm':[1e-15,0.96e-15,0.9e-15],
                             'm1_ev':[0.005,0.0049,0.0048],
                             'delta_m21_sq_ev2':[7.9e-5,7.7e-5,7.5e-5]}},
    }


def test_joins_only_adjacent_real_solver_samples():
    curve = mod._joint_trace(make_stages(),'g2')
    assert curve['x'] == [1e7, 1e6, 1e5, 1e4, 1e3, 5e2, 1e2]
    assert curve['y'] == [0.6, 0.59, 0.58, 0.56, 0.55, 0.54, 0.53]
    assert mod._joint_trace(make_stages(),'c5_fnorm')['x'] == [1e3,5e2,1e2]
    broken=make_stages()
    broken['eft1']['series']['g2'][0]=0.45
    with pytest.raises(ValueError,match='threshold endpoint mismatch'):
        mod._joint_trace(broken,'g2')
    no_eft=make_stages(); del no_eft['eft1']['series']['g2']
    with pytest.raises(ValueError,match='uncomputed interval'):
        mod._joint_trace(no_eft,'g2')


def test_original_three_colours_and_controls_survive_adapter():
    data, meta = mod._assemble_curves({'A_-2': make_stages()})
    page=mod.actual_original_template_page(sample_template(),data,meta)
    assert "#2466ae','#c46c25" in page
    assert "uv:'#eaf2ff',int:'#fff1cf',sm:'#e5f5e9'" in page
    assert 'id="zoomUV"' in page and 'id="zoomInt"' in page and 'id="zoomFinal"' in page
    assert 'id="yScaleMode"' in page
    assert 'actual ODE samples' in page
    assert 'Yellow region: physical EFT1 C5 / LLSS Wilson evolution has NOT been computed' in page
    assert '__DATA__' not in page and '__META__' not in page
    assert '"c5_fnorm"' in page and '"m1_ev"' in page
    assert data['matched'] == {} and data['eft1_wilson_c5_running_certified'] is False
    assert data['quantities']['c5_fnorm']['A(-2)']['x'][0] == 1e3


def test_refuses_missing_actual_original_template():
    with pytest.raises(ValueError,match='not the original coloured'):
        mod.actual_original_template_page('<html>hand-made replacement</html>',{}, {})


def test_generates_original_filename_and_removes_parallel_dashboard(tmp_path,monkeypatch):
    run='fixture';folder=tmp_path/run; model_dir=folder/'A_m2';model_dir.mkdir(parents=True)
    (folder/'manifest.json').write_text(json.dumps({'interactive_rge_generated':True}))
    (model_dir/'centre.json').write_text(json.dumps({'model':'A_-2','status':'Success'}))
    sample=make_stages()
    (model_dir/'rge_trajectories.json').write_text(json.dumps({
        'format':mod.UV_EFT_FORMAT,'status':'Success','model':'A_-2',
        'uv':sample['uv'],'eft1':sample['eft1'],
    }))
    (folder/'interactive_RGE_with_SMEFT.html').write_text('previous wrong dashboard')
    template_file=tmp_path/'template.html';template_file.write_text(sample_template())
    def fake_smeft(centre, points):
        vals={k:[0.,0.,0.] for k in mod.SMEFT_KEYS}
        vals.update(sample['smeft']['series'])
        return {'model':'A_-2','format':mod.SMEFT_FORMAT,'status':'Success',
                'mu_gev':sample['smeft']['mu_gev'],'series':vals,
                'relative_endpoint_c5_residual':0.,'nfev':7}
    monkeypatch.setattr(mod,'sample_smeft',fake_smeft)
    report=mod.run(run,3,root=tmp_path,template_path=template_file)
    assert report['smeft_model_count']==1
    out=(folder/'interactive_RGE.html')
    assert out.exists()
    assert not (folder/'interactive_RGE_with_SMEFT.html').exists()
    assert '#eaf2ff' in out.read_text() and '#e5f5e9' in out.read_text()
    assert (model_dir/'smeft_trajectories.json').exists()
    assert (model_dir/'smeft_observables.csv').exists()
    assert (folder/'original_coloured_rge_status.json').exists()


def test_installed_original_template_contract():
    source=mod.TEMPLATE
    if not source.exists():
        pytest.skip('Real T3 repository template not mounted in offline test runtime')
    original=source.read_text(encoding='utf-8')
    data,meta=mod._assemble_curves({'A_-2': make_stages()})
    page=mod.actual_original_template_page(original,data,meta)
    assert "rectRegion(s.mu_scalar_threshold,s.mu_low,regionColours.sm,'SM + Weinberg EFT');" in page
    assert "const COLORS=['#2466ae','#c46c25'" in page
    assert 'id="zoomFinal"' in page


def certified_stub_module(monkeypatch, *, status=None, source_f=1.3e-15):
    """Mock only the existing external solver; do not claim physical calibration."""
    import sys
    import types
    module = types.ModuleType('RGE.running.intermediate.GeneralEFT1C5Mixing')
    module.CERTIFICATES = ('matchete_gamma_cc_to_weyl_factor_and_sign',
                           'canonical_C_ijab_operator_normalization')
    def fake_solver(raw, *, allow_synthetic, trajectory_points):
        assert allow_synthetic is False
        assert raw['authority'] == 'certified'
        def enc(value):
            return {'real':(np.eye(3)*value).tolist(), 'imag':np.zeros((3,3)).tolist()}
        scales=np.geomspace(1e5,1e3,trajectory_points)
        values=[source_f*np.log(m/1e5)/np.log(1e3/1e5) for m in scales]
        return {'format':mod.EFT1_SOURCE_FORMAT, 'model':raw['model'],
                'status':status or mod.EFT1_PHYSICAL_STATUS,
                'authority':'certified',
                'physics_completeness':{'physical_eft1_source_verified':True,
                                        'joint_two_scalar_finite_matching_included':False},
                'scales_gev':{'mu_F':1e5,'mu_S':1e3},
                'trajectory':[{'mu_gev':float(m),'delta_C5_log_gev_inverse':enc(v)}
                              for m,v in zip(scales,values)],
                'total_delta_C5_log_gev_inverse':enc(source_f)}
    module.solve_boundary = fake_solver
    monkeypatch.setitem(sys.modules,module.__name__,module)
    return module


def test_eft1_rejects_synthetic_and_missing_certifications(tmp_path,monkeypatch):
    ext=certified_stub_module(monkeypatch)
    source=tmp_path/'boundary.json'
    source.write_text(json.dumps({'model':'A_-2','authority':'synthetic','certificates':{}}))
    with pytest.raises(ValueError, match='Synthetic'):
        mod.certified_eft1_source('A_-2',source,3)
    source.write_text(json.dumps({'model':'A_-2','authority':'certified','certificates':{}}))
    with pytest.raises(ValueError, match='Missing independent'):
        mod.certified_eft1_source('A_-2',source,3)
    good={x:True for x in ext.CERTIFICATES}
    source.write_text(json.dumps({'model':'A_-2','authority':'certified','certificates':good}))
    result, stage = mod.certified_eft1_source('A_-2',source,3)
    assert result['status']==mod.EFT1_PHYSICAL_STATUS
    assert np.isclose(stage['series']['eft1_delta_c5_fnorm'][0],0)
    assert np.isclose(stage['series']['eft1_delta_c5_ee_abs'][-1],1.3e-15)
    assert len(stage['mu_gev'])==3


def test_eft1_source_is_separate_yellow_region_and_does_not_join_smeft():
    stages=make_stages()
    stages['eft1_c5']={
        'mu_gev':[1e5,1e4,1e3],
        'series':{'eft1_delta_c5_fnorm':[0,1e-15,2e-15]},
    }
    data,meta=mod._assemble_curves({'A_-2':stages})
    assert data['eft1_leading_log_source_only'] is True
    assert data['eft1_wilson_c5_running_certified'] is False
    assert data['quantities']['eft1_delta_c5_fnorm']['A(-2)']['x']==[1e5,1e4,1e3]
    assert data['quantities']['c5_fnorm']['A(-2)']['x']==[1e3,5e2,1e2]
    assert meta['eft1_delta_c5_fnorm']['eft1_only'] is True
    page=mod.actual_original_template_page(sample_template(),data,meta)
    assert 'externally certified one-loop LLSS-to-LLHH source' in page
    assert 'if(META[q]&&META[q].eft1_only)' in page
    assert 'id="zoomInt"' in page


def test_missing_certified_input_remains_unplotted_and_readiness_written(tmp_path,monkeypatch):
    run='fixture';folder=tmp_path/run; model_dir=folder/'A_m2';model_dir.mkdir(parents=True)
    (folder/'manifest.json').write_text(json.dumps({'interactive_rge_generated':True}))
    (model_dir/'centre.json').write_text(json.dumps({'model':'A_-2','status':'Success'}))
    sample=make_stages()
    (model_dir/'rge_trajectories.json').write_text(json.dumps({
        'format':mod.UV_EFT_FORMAT,'status':'Success','model':'A_-2',
        'uv':sample['uv'],'eft1':sample['eft1'],
    }))
    template_file=tmp_path/'template.html';template_file.write_text(sample_template())
    def fake_smeft(centre, points):
        vals={k:[0.,0.,0.] for k in mod.SMEFT_KEYS}
        vals.update(sample['smeft']['series'])
        return {'model':'A_-2','format':mod.SMEFT_FORMAT,'status':'Success',
                'mu_gev':sample['smeft']['mu_gev'],'series':vals,
                'relative_endpoint_c5_residual':0.,'nfev':7}
    monkeypatch.setattr(mod,'sample_smeft',fake_smeft)
    result=mod.run(run,3,root=tmp_path,template_path=template_file,
                   certified_root=tmp_path/'missing_certified')
    assert result['eft1_certified_source_count']==0
    readiness=json.loads((folder/'eft1_c5_readiness.json').read_text())
    assert readiness['models'][0]['status']=='MissingCertifiedBoundary'
    html=(folder/'interactive_RGE.html').read_text()
    assert '"eft1_delta_c5_fnorm"' not in html
    assert '"c5_fnorm"' in html


def test_cleanup_removes_only_named_obsolete_files(tmp_path):
    old=tmp_path/'scripts/BuildT3SharedA2ThreeStageRGE.py';old.parent.mkdir(parents=True)
    old.write_text('old dashboard')
    keep=tmp_path/'scripts/BuildT3OriginalColouredRGE.py';keep.write_text('KEEP')
    source=tmp_path/'output/rge_trajectories.json';source.parent.mkdir(parents=True)
    source.write_text('{"source":true}')
    stale=tmp_path/'output/interactive_RGE_with_SMEFT.html';stale.write_text('stale')
    removed=mod.remove_obsolete_dashboard_files(tmp_path, tmp_path/'output', clean_scripts=True)
    assert len(removed)==2
    assert not old.exists() and not stale.exists()
    assert keep.read_text()=='KEEP' and source.is_file()


def test_certified_source_is_plotted_only_in_original_yellow_stage(tmp_path, monkeypatch):
    ext=certified_stub_module(monkeypatch)
    run='fixture'; folder=tmp_path/run; model_dir=folder/'A_m2';model_dir.mkdir(parents=True)
    source_folder=tmp_path/'certified/A_m2';source_folder.mkdir(parents=True)
    (source_folder/'boundary.json').write_text(json.dumps({
        'model':'A_-2','authority':'certified',
        'certificates':{x:True for x in ext.CERTIFICATES},
    }))
    (folder/'manifest.json').write_text(json.dumps({'interactive_rge_generated':True}))
    (model_dir/'centre.json').write_text(json.dumps({'model':'A_-2','status':'Success'}))
    stages=make_stages()
    (model_dir/'rge_trajectories.json').write_text(json.dumps({
        'format':mod.UV_EFT_FORMAT,'status':'Success','model':'A_-2',
        'uv':stages['uv'],'eft1':stages['eft1'],
    }))
    def fake_smeft(centre, points):
        vals={key:[0.,0.,0.] for key in mod.SMEFT_KEYS}
        vals.update(stages['smeft']['series'])
        return {'model':'A_-2','format':mod.SMEFT_FORMAT,'status':'Success',
                'mu_gev':stages['smeft']['mu_gev'],'series':vals,
                'relative_endpoint_c5_residual':0.,'nfev':7}
    monkeypatch.setattr(mod,'sample_smeft',fake_smeft)
    template_file=tmp_path/'template.html';template_file.write_text(sample_template())
    result=mod.run(run,3,root=tmp_path,template_path=template_file,
                   certified_root=tmp_path/'certified')
    assert result['eft1_certified_source_count']==1
    assert (model_dir/'eft1_c5_source.json').is_file()
    readiness=json.loads((folder/'eft1_c5_readiness.json').read_text())
    assert readiness['models'][0]['status']=='CertifiedSourceComputed'
    page=(folder/'interactive_RGE.html').read_text()
    assert '"eft1_delta_c5_fnorm"' in page and '"c5_fnorm"' in page
    assert page.count('"eft1_delta_c5_ee_abs"')>=1
    assert 'if(META[q]&&META[q].eft1_only)' in page
