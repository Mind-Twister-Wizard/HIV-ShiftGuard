@echo off
setlocal EnableExtensions
cd /d "%~dp0"
title HIV-ShiftGuard - One Click Reproducible Research
color 1F
echo ============================================================
echo    HIV-SHIFTGUARD - DATA, ML, ABLATIONS, TABLES AND FIGURES
 echo    RESEARCH ONLY - NOT A HIV CURE OR CLINICAL TOOL
 echo ============================================================
 echo.
if exist ".venv\Scripts\python.exe" goto :have_venv
set "BASE_PY="
where py >nul 2>nul
if not errorlevel 1 (
    py -3.12 -c "import sys; assert sys.version_info >= (3,10)" >nul 2>nul
    if not errorlevel 1 set "BASE_PY=py -3.12"
    if not defined BASE_PY (
        py -3.11 -c "import sys; assert sys.version_info >= (3,10)" >nul 2>nul
        if not errorlevel 1 set "BASE_PY=py -3.11"
    )
    if not defined BASE_PY (
        py -3.13 -c "import sys; assert sys.version_info >= (3,10)" >nul 2>nul
        if not errorlevel 1 set "BASE_PY=py -3.13"
    )
)
if not defined BASE_PY (
    where python >nul 2>nul
    if not errorlevel 1 (
        python -c "import sys; assert sys.version_info >= (3,10)" >nul 2>nul
        if not errorlevel 1 set "BASE_PY=python"
    )
)
if not defined BASE_PY (
    echo Python 3.10 or newer not found. Attempting Python 3.12 installation using winget.
    where winget >nul 2>nul
    if errorlevel 1 goto :no_python
    winget install --exact --id Python.Python.3.12 --source winget --accept-package-agreements --accept-source-agreements
    if errorlevel 1 goto :no_python
    echo Please CLOSE and REOPEN this BAT file so the Python PATH updates.
    pause
    exit /b 0
)
echo Creating project-only virtual environment with %BASE_PY% ...
%BASE_PY% -m venv ".venv"
if errorlevel 1 goto :fail
:have_venv
set "PYTHONPATH=%~dp0src"
if not exist ".venv\Scripts\pip.exe" goto :fail
if not exist ".venv\.deps_installed" (
    echo Installing required packages automatically...
    ".venv\Scripts\python.exe" -m pip install --disable-pip-version-check --upgrade pip
    if errorlevel 1 goto :fail
    ".venv\Scripts\python.exe" -m pip install --disable-pip-version-check -r requirements.txt
    if errorlevel 1 goto :fail
    echo installed > ".venv\.deps_installed"
)
echo.
echo Testing code before running experiments...
".venv\Scripts\python.exe" -m unittest discover -s tests -v
if errorlevel 1 goto :fail
 echo.
echo Running or RESUMING real-data experiments...
".venv\Scripts\python.exe" -m hivshiftguard.pipeline --config config.json
if errorlevel 1 goto :fail
 echo.
echo Completed. Review results\latest_report.html
start "" "%~dp0results\latest_report.html"
echo Press any key to close.
pause >nul
exit /b 0
:no_python
 echo ERROR: Could not install Python automatically. Install Python 3.12 from https://www.python.org/downloads/ and rerun.
 goto :fail
:fail
 echo.
echo FAILED or INTERRUPTED. Nothing is lost; rerun RUN_ALL.bat to resume.
 echo Inspect results\run.log or check Internet connection / data\raw files.
pause
exit /b 1
