@echo off
rem ====================================================================
rem  build.bat
rem    MPS simulation (mps.c) and the DirectX 11 viewer (mps_viewer.cpp)
rem
rem    Usage:  double-click this file, or run it from a command prompt.
rem            You do NOT need the "Developer Command Prompt" --
rem            this script locates Visual Studio by itself.
rem
rem    Output: mps.exe         ... the MPS solver
rem            mps_viewer.exe  ... the DirectX 11 animation viewer
rem ====================================================================
setlocal

cd /d "%~dp0"

rem --- Locate Visual Studio with vswhere.exe --------------------------
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" (
    echo [ERROR] vswhere.exe not found.
    echo         Visual Studio 2017 or later with the
    echo         "Desktop development with C++" workload is required.
    goto :failed
)

set "VSPATH="
for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -products * ^
    -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 ^
    -property installationPath`) do set "VSPATH=%%i"

if not defined VSPATH (
    echo [ERROR] No Visual Studio installation with the C++ toolset was found.
    goto :failed
)

set "VCVARS=%VSPATH%\VC\Auxiliary\Build\vcvars64.bat"
if not exist "%VCVARS%" (
    echo [ERROR] vcvars64.bat not found under "%VSPATH%".
    goto :failed
)

echo Using Visual Studio at: %VSPATH%
call "%VCVARS%" >nul
if errorlevel 1 goto :failed

rem --- Build the solver -----------------------------------------------
echo.
echo [1/2] Building mps.exe ...
cl /W3 /O2 /nologo mps.c /Fe:mps.exe /Fo:mps_c.obj
if errorlevel 1 goto :failed

rem --- Build the viewer -----------------------------------------------
echo.
echo [2/2] Building mps_viewer.exe ...
cl /W3 /O2 /EHsc /nologo mps_viewer.cpp /Fe:mps_viewer.exe /Fo:mps_viewer.obj
if errorlevel 1 goto :failed

rem --- Tidy up intermediate files -------------------------------------
if exist mps_c.obj      del mps_c.obj
if exist mps_viewer.obj del mps_viewer.obj

echo.
echo ==========================================================
echo  Build succeeded.
echo.
echo   1) Run the simulation   :  mps.exe
echo        (about 25 seconds; writes output_*.prof / particle_*.vtu)
echo   2) Watch the animation  :  mps_viewer.exe
echo ==========================================================
echo.
pause
exit /b 0

:failed
echo.
echo ******* BUILD FAILED *******
echo.
pause
exit /b 1
