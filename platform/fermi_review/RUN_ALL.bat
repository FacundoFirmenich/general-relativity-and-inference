@echo off
cd /d "%~dp0"
py -3 -m venv .venv
if errorlevel 1 exit /b 1
.venv\Scripts\python.exe -m pip install -r requirements.txt
if errorlevel 1 exit /b 1
.venv\Scripts\python.exe repaired\catalog_audit.py
if errorlevel 1 exit /b 1
.venv\Scripts\python.exe repaired\audit_originals.py
if errorlevel 1 exit /b 1
.venv\Scripts\python.exe repaired\audit_hunter.py
if errorlevel 1 exit /b 1
pause
