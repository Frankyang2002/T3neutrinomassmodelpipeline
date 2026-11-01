from pathlib import Path
import importlib.util
s=Path(__file__).resolve().parents[2]/'scripts/ReadableOriginalGeneralReports.py'
spec=importlib.util.spec_from_file_location('readable_original',s)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_nested_native():
    x=r'\text{Bar}(\text{Field}(S_1,\text{Scalar},\{\text{Index}(i,SU2)\},\{\}))'
    out=m.compact_native_tex(x)
    assert 'S_1' in out and 'dagger' in out and 'Field' not in out

def test_minimum_size():
    x=r'\documentclass[8pt]{article}'+'\n'+r'\tiny \scriptsize'
    y=m.improve_latex(x)
    assert r'\documentclass[10pt,a3paper]' in y and r'\tiny' not in y and r'\scriptsize' not in y
