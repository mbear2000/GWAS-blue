@echo off
setlocal
title Stop GWAS BLUE local service

echo ============================================================
echo Stop GWAS BLUE local service
echo Target: http://127.0.0.1:8766/
echo ============================================================
echo.

powershell -NoProfile -ExecutionPolicy Bypass -Command ^
  "$ErrorActionPreference='Stop';" ^
  "$conns = @(Get-NetTCPConnection -LocalPort 8766 -State Listen -ErrorAction SilentlyContinue);" ^
  "if(-not $conns -or $conns.Count -eq 0){" ^
  "  Write-Host 'No process is listening on port 8766. GWAS BLUE is already stopped.';" ^
  "  exit 0" ^
  "}" ^
  "$stopped=$false;" ^
  "foreach($c in $conns){" ^
  "  $pidNum=$c.OwningProcess;" ^
  "  $proc=Get-CimInstance Win32_Process -Filter ('ProcessId='+$pidNum) -ErrorAction SilentlyContinue;" ^
  "  if($proc -and $proc.CommandLine -match 'server_blue\.py'){" ^
  "    Write-Host ('Stopping server_blue.py, PID '+$pidNum+' ...');" ^
  "    Stop-Process -Id $pidNum -Force;" ^
  "    $stopped=$true;" ^
  "  } else {" ^
  "    Write-Host ('Port 8766 is owned by PID '+$pidNum+', but it does not look like server_blue.py. It was NOT stopped.');" ^
  "  }" ^
  "}" ^
  "Start-Sleep -Milliseconds 500;" ^
  "$left=@(Get-NetTCPConnection -LocalPort 8766 -State Listen -ErrorAction SilentlyContinue);" ^
  "if($left.Count -eq 0){" ^
  "  Write-Host 'GWAS BLUE local service stopped successfully.';" ^
  "  exit 0" ^
  "} else {" ^
  "  Write-Host 'Port 8766 is still listening. Please check the process manually.';" ^
  "  exit 1" ^
  "}"

echo.
pause
endlocal
