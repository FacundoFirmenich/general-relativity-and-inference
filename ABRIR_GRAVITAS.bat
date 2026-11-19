@echo off
cd /d "%~dp0"
where py >nul 2>nul
if %ERRORLEVEL% EQU 0 (
  py -3 launch_gravitas.py
) else (
  python launch_gravitas.py
)
if %ERRORLEVEL% NEQ 0 (
  echo Se necesita Python 3.11 o superior. Alternativa sin instalacion: cuaderno Colab.
  pause
)
