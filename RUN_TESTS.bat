@echo off
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
 echo First run RUN_ALL.bat to prepare the environment.
 pause
 exit /b 1
)
set "PYTHONPATH=%~dp0src"
".venv\Scripts\python.exe" -m unittest discover -s tests -v
pause
