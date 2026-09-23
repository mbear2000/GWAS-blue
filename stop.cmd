@echo off
setlocal
cd /d "%~dp0"
echo [GWAS] Stopping local GWAS console...

powershell.exe -NoProfile -ExecutionPolicy Bypass -Command ^
  "$root=(Resolve-Path '.').Path; " ^
  "$procs=Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -and ($_.CommandLine -like ('*' + $root + '*launcher.py*')) }; " ^
  "$procs | ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }"

timeout /t 1 /nobreak >nul

powershell.exe -NoProfile -ExecutionPolicy Bypass -Command ^
  "$root=(Resolve-Path '.').Path; " ^
  "$procs=Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -and ($_.CommandLine -like ('*' + $root + '*server.py*')) }; " ^
  "$procs | ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }"

timeout /t 1 /nobreak >nul

powershell.exe -NoProfile -ExecutionPolicy Bypass -Command ^
  "$root=(Resolve-Path '.').Path; " ^
  "$left=Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -and (($_.CommandLine -like ('*' + $root + '*launcher.py*')) -or ($_.CommandLine -like ('*' + $root + '*server.py*'))) }; " ^
  "if ($left) { Write-Host '[GWAS] Some local processes are still running:'; $left | Select-Object ProcessId,CommandLine | Format-Table -AutoSize } else { Write-Host '[GWAS] Local GWAS console stopped.' }"

endlocal
pause
