@echo off
setlocal
cd /d "%~dp0"

set GWAS_WEB_PORT=8766

echo.
echo GWAS BLUE complete v9.2
echo http://127.0.0.1:8766/
echo.
echo Required original files:
echo   server.py
echo   workflow_v2.sh
echo   web\index.html
echo.
echo Keep this window open.
echo.

start "" powershell -NoProfile -WindowStyle Hidden -Command "Start-Sleep -Seconds 2; Start-Process 'http://127.0.0.1:8766/'"
python -B server_blue.py
pause
