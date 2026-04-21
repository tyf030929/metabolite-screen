$ErrorActionPreference = "Continue"
$workdir = "E:\基因组组装\网络药理学\11"
Set-Location $workdir

$lines = Get-Content -Path "$workdir\app.py" -Encoding UTF8

# Find the line indices (0-based) we want to manipulate
# Lines 409-415 (1-based) are the old Aura comment + duplicate MetaboLab block
# We want to keep the file structure but replace lines 409-415

Write-Output "Total lines: $($lines.Count)"
Write-Output "Lines around old CSS comment:"
for ($i = 406; $i -lt 420; $i++) {
    Write-Output "L$($i+1): $($lines[$i].Substring(0, [Math]::Min(80, $lines[$i].Length)))"
}
