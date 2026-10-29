param(
  [string]$RunId = "overnight_six_v1",
  [int]$SobolPoints = 128,
  [int]$MaxNfev = 100,
  [switch]$AllowPoorA2,
  [switch]$Resume
)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Set-Location $root
$logRoot = Join-Path $root 'output\studies\general_scalar_extended\six_benchmarks'
New-Item -ItemType Directory -Path $logRoot -Force | Out-Null
$sessionLog = Join-Path $logRoot ($RunId + '_powershell.log')
$arguments = @('-u','-X','utf8','-m','scripts.RunT3GeneralScalarOvernight',
  '--run-id',$RunId,'--sobol-points',"$SobolPoints",'--max-nfev',"$MaxNfev")
if ($AllowPoorA2) { $arguments += '--allow-poor-a2' }
$arguments += '--resume' # Safe for both fresh and previously prepared runs; hashes are checked
"Overnight started $(Get-Date -Format o)" | Out-File -Encoding utf8 -Append $sessionLog
& python @arguments 2>&1 | Tee-Object -FilePath $sessionLog -Append
$code = $LASTEXITCODE
"Overnight finished exit=$code $(Get-Date -Format o)" | Out-File -Encoding utf8 -Append $sessionLog
exit $code
