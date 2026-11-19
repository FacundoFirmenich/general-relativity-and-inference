@echo off
cd /d "%~dp0"
where py >nul 2>nul
if %ERRORLEVEL% EQU 0 (
  py -3 launch_gravitas.py
) else (
  python launch_gravitas.py
)
if %ERRORLEVEL% NEQ 0 (
  echo Python 3.11 or newer is required. Browser alternative: the self-contained Colab notebook.
  pause
)
