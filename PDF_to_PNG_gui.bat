@echo off
rem ============================================================
rem  PDF Pages to PNG - high-res  (GUI launcher)
rem  Double-click to open. Uses the local .venv Python if present,
rem  otherwise falls back to the system Python on PATH.
rem ============================================================
setlocal
set "SCRIPT_DIR=%~dp0"
cd /d "%SCRIPT_DIR%"

set "PY="
if exist "%SCRIPT_DIR%.venv\Scripts\python.exe" (
    set "PY=%SCRIPT_DIR%.venv\Scripts\python.exe"
) else (
    where python >nul 2>nul && set "PY=python"
)

if not defined PY (
    echo [ERROR] No Python found. Install Python or create the .venv.
    pause
    exit /b 1
)

echo Starting PDF Pages to PNG GUI...
"%PY%" "%SCRIPT_DIR%pdf_pages_to_png_gui.py"
if errorlevel 1 (
    echo.
    echo [ERROR] GUI exited with an error (see message above).
    pause
)
endlocal
