@echo off
REM GlycoLipid-Insight Project Startup Script
REM This script activates the 'chem_ai' environment and opens the project.

echo [INFO] Activating chem_ai environment...
call "C:\Users\Daniel Yang\anaconda3\Scripts\activate.bat" chem_ai

if errorlevel 1 (
    echo [ERROR] Failed to activate environment 'chem_ai'.
    pass
)

echo [INFO] Checking python...
where python

echo [INFO] Restoring environment...
pip install rdkit pandas tqdm openpyxl

echo [INFO] Opening VS Code...
code .
pause
