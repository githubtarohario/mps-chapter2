@echo off
rem ====================================================================
rem  make_pdf.bat
rem    技術解説書.md  ->  技術解説書.pdf
rem
rem    Usage: double-click this file.
rem
rem    Requires:
rem      - Python 3 with:  pip install markdown pypdf reportlab
rem      - Google Chrome (or Microsoft Edge)
rem
rem    The actual work is done by pdf_tools\build_pdf.ps1
rem ====================================================================
setlocal
cd /d "%~dp0"

powershell -NoProfile -ExecutionPolicy Bypass -File "pdf_tools\build_pdf.ps1" %*
set RC=%ERRORLEVEL%

echo.
if not "%RC%"=="0" (
    echo ******* FAILED *******
)
pause
exit /b %RC%
