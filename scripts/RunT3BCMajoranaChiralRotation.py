"""Run isolated B/C native chiral probe; do not modify production physics."""
import json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.experimental.T3BCMajoranaChiralRotation import check

def main():
    out=ROOT/'output/studies/t3_bc_majorana_chiral_rotation_v1';out.mkdir(parents=True,exist_ok=True)
    witness=check(); (out/'matrix_witness.json').write_text(json.dumps(witness,indent=2)+'\n')
    native=out/'native_matchete.json';native.unlink(missing_ok=True)
    env=os.environ.copy();env.update(T3_EXPERIMENTAL_GENERAL_FERMIONS='1',T3_EXPERIMENTAL_FLAVOR_MASSES='1',T3_EXPERIMENTAL_GENERAL_QUARTICS='1',T3_EXPERIMENTAL_GENERAL_QUADRATICS='0')
    cmd=['wolframscript','-file',str(ROOT/'Lagrangian/ProbeT3BCMajoranaNativeChiral.wl'),str(native)]
    try:
        proc=subprocess.run(cmd,cwd=ROOT,env=env,capture_output=True,text=True,errors='replace')
        (out/'native_matchete.log').write_text('COMMAND: '+' '.join(cmd)+'\nSTDOUT:\n'+proc.stdout+'\nSTDERR:\n'+proc.stderr,encoding='utf-8')
        code=proc.returncode
    except OSError as exc:
        (out/'native_matchete.log').write_text(str(exc));code=3
    summary={'status':'NativeAbstractProbePassedPhysicalDictionaryPending' if code==0 else 'NativeProbeBlockedOrFailed','native_exit_code':code,'matrix_witness':witness,'native_json_exists':native.is_file(),'physical_C5_certified':False,'production_files_modified':False}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(summary['status']);print('Output:',out)
    return 0 if code==0 else code
if __name__=='__main__':raise SystemExit(main())
