@echo off
setlocal
title Stop GWAS BLUE local service

echo ============================================================
echo Stop GWAS BLUE local service
echo Target: http://127.0.0.1:8766/
echo ============================================================
echo.
echo NOTE: this only stops local web monitoring/state.
echo It does NOT qdel remote PBS jobs and does NOT delete server results.
echo.

powershell -NoProfile -ExecutionPolicy Bypass -Command ^
  "$ErrorActionPreference='Continue';" ^
  "$conns=@(Get-NetTCPConnection -LocalPort 8766 -State Listen -ErrorAction SilentlyContinue);" ^
  "if($conns.Count -gt 0){foreach($c in $conns){$pidNum=$c.OwningProcess;$proc=Get-CimInstance Win32_Process -Filter ('ProcessId='+$pidNum) -ErrorAction SilentlyContinue;if($proc -and $proc.CommandLine -match 'server_blue\.py'){Write-Host ('Stopping server_blue.py, PID '+$pidNum+' ...');Stop-Process -Id $pidNum -Force -ErrorAction SilentlyContinue}else{Write-Host ('Port 8766 PID '+$pidNum+' is not server_blue.py; NOT stopped.')}}}else{Write-Host 'No process is listening on port 8766.'};" ^
  "Start-Sleep -Milliseconds 500;" ^
  "$runs=Join-Path (Get-Location) 'runs';$changed=0;" ^
  "if(Test-Path $runs){Get-ChildItem $runs -Directory | ForEach-Object {$state=Join-Path $_.FullName 'state.txt';if(Test-Path $state){try{$s=(Get-Content $state -Raw -ErrorAction Stop).Trim()}catch{$s=''};if($s -match '^(connecting|running)\|'){$msg='cancelled|用户已通过 stop_blue.cmd 中断本地网页监控；远程队列未自动 qdel';Set-Content $state $msg -Encoding UTF8;Set-Content (Join-Path $_.FullName 'cancelled.flag') '用户已通过 stop_blue.cmd 中断本地网页监控；远程队列/进程未自动 qdel。' -Encoding UTF8;Write-Host ($_.Name+' => '+$msg);$changed++}}}};" ^
  "if($changed -eq 0){Write-Host 'No local connecting/running states needed cleanup.'}else{Write-Host ('Cleaned '+$changed+' local running state(s).')};" ^
  "$left=@(Get-NetTCPConnection -LocalPort 8766 -State Listen -ErrorAction SilentlyContinue);if($left.Count -eq 0){Write-Host 'GWAS BLUE local service is stopped.'}else{Write-Host 'WARNING: port 8766 is still listening.'}"

echo.
pause
endlocal
