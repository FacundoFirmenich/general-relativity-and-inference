@echo off
setlocal
cd /d "%~dp0"
set "GRAVITAS_MODE=%~1"
if not defined GRAVITAS_MODE set "GRAVITAS_MODE=all"
py -3 launch_gravitas.py "%GRAVITAS_MODE%"
if errorlevel 1 exit /b 1
endlocal
