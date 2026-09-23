@echo off
setlocal
cd /d "%~dp0"

set "PYTHONW=C:\Python314\pythonw.exe"

if not exist "%PYTHONW%" (
    echo [GWAS] Python not found:
    echo %PYTHONW%
    pause
    exit /b 1
)

if not exist "%~dp0launcher.py" (
    echo [GWAS] launcher.py not found:
    echo %~dp0launcher.py
    pause
    exit /b 1
)

echo [GWAS] Starting local GWAS console...
start "" "%PYTHONW%" "%~dp0launcher.py" --open

endlocal
