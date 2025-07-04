@echo off
setlocal

REM === Set script directory safely ===
set "BASEDIR=%~dp0"
set "PYTHON_VERSION=3.12.2"
set "PYTHON_INSTALLER=python-%PYTHON_VERSION%-amd64.exe"
set "PYTHON_DIR=%BASEDIR%py"
set "VENV_DIR=%PYTHON_DIR%\venv"
set "APP_FILE=%BASEDIR%point_to_raster_app.py"

REM === Create directories ===
if not exist "%PYTHON_DIR%" mkdir "%PYTHON_DIR%"

REM === Check if Python exists ===
if not exist "%PYTHON_DIR%\python.exe" (
    echo [INFO] Python not found locally. Downloading and installing...

    REM Use PowerShell to download (quote carefully!)
    powershell -Command "Invoke-WebRequest -Uri 'https://www.python.org/ftp/python/%PYTHON_VERSION%/%PYTHON_INSTALLER%' -OutFile '%BASEDIR%%PYTHON_INSTALLER%'"

    REM Install Python silently
    "%BASEDIR%%PYTHON_INSTALLER%" /quiet InstallAllUsers=0 TargetDir="%PYTHON_DIR%" PrependPath=0 Include_launcher=0 Include_pip=1

    del "%BASEDIR%%PYTHON_INSTALLER%"
)

REM === Create virtual environment if it doesn't exist ===
if not exist "%VENV_DIR%\Scripts\python.exe" (
    echo [INFO] Creating virtual environment...
    "%PYTHON_DIR%\python.exe" -m venv "%VENV_DIR%"
)

REM === Activate virtual environment ===
call "%VENV_DIR%\Scripts\activate.bat"

REM === Upgrade pip and install required libraries ===
echo [INFO] Installing dependencies...
python -m pip install --upgrade pip
python -m pip install pandas numpy geopandas rasterio shapely pyproj scipy affine

REM === Run the GUI app ===
echo [INFO] Launching the application...
START "" python "%APP_FILE%"

endlocal
pause
