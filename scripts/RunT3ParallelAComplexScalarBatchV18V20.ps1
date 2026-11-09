# v18-v20: native Matchete CGs + native RGBeta complex-channel bridge + 178-slot census.
# Does not modify GitHub, installed packages, shared project sources or physical Wilson output.
param([switch]$VerifyOnly)
$ErrorActionPreference='Stop'
$root=Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $root
$python=Join-Path $root '.venv\Scripts\python.exe'
if(-not (Test-Path -LiteralPath $python)){throw "Missing project Python $python"}
function RequireExit([string]$label){ if($LASTEXITCODE -ne 0){throw "$label BLOCKED exit $LASTEXITCODE"} }
& $python -m unittest discover -s tests -p 'test_t3_parallel_a_complex_scalar_v18_v20.py' -v
RequireExit 'v18-v20 tests'
$base='output\parallel_t3\wilson'
$topo=Join-Path $base 'matchete_topology_native_v18.json'
$bc=Join-Path $base 'matchete_bc_extra_native_v18.json'
$flags=@('T3_EXPERIMENTAL_GENERAL_FERMIONS','T3_EXPERIMENTAL_GENERAL_QUADRATICS',
         'T3_EXPERIMENTAL_GENERAL_QUARTICS','T3_EXPERIMENTAL_FLAVOR_MASSES')
$previous=@{}
foreach($f in $flags){$previous[$f]=[Environment]::GetEnvironmentVariable($f,'Process')}
try{
    foreach($f in $flags){[Environment]::SetEnvironmentVariable($f,$null,'Process')}
    if(-not (Test-Path -LiteralPath $topo)){
        if($VerifyOnly){throw 'Missing native v18 topology output in VerifyOnly mode'}
        if(-not (Get-Command wolframscript -ErrorAction SilentlyContinue)){throw 'wolframscript unavailable'}
        & wolframscript -file scripts\RunT3ParallelAExportComplexMatcheteV18.wl $topo TOPOLOGY
        RequireExit 'v18 native Matchete all16 topology export'
    }else{Write-Host 'v18 topology witness exists; Python will revalidate source hashes'}
    [Environment]::SetEnvironmentVariable('T3_EXPERIMENTAL_GENERAL_QUARTICS','1','Process')
    if(-not (Test-Path -LiteralPath $bc)){
        if($VerifyOnly){throw 'Missing native v18 B/C extra CG output in VerifyOnly mode'}
        if(-not (Get-Command wolframscript -ErrorAction SilentlyContinue)){throw 'wolframscript unavailable'}
        & wolframscript -file scripts\RunT3ParallelAExportComplexMatcheteV18.wl $bc BCEXTRA
        RequireExit 'v18 native Matchete B/C extra CG export'
    }else{Write-Host 'v18 B/C witness exists; Python will revalidate source hashes'}
}finally{
    foreach($f in $flags){[Environment]::SetEnvironmentVariable($f,$previous[$f],'Process')}
    Write-Host 'Diagnostic optional flags restored in this process'
}
& $python -m scripts.RunT3ParallelAComplexBridgeV19 --project-root . --resume
RequireExit 'v19 complex source-conditioned Matchete/RGBeta bridge'
& $python -m scripts.RunT3ParallelAFullScalarCensusV20 --project-root . --resume
RequireExit 'v20 source-conditioned 178 quartic census'
Write-Host 'v18-v20 complete conditionally. NO physical Wilson evolution.'
