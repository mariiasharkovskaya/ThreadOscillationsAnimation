@echo off
setlocal enabledelayedexpansion

SET MINOR=
for /f "tokens=2 delims=." %%a in ('python --version') do (
    set MINOR=%%a
)

echo Compiling forcedOsc
g++ forcedOsc.cpp -o forcedOsc -I C:/msys64/mingw64/include -L/msys64/mingw64/lib -lgsl -lgslcblas -I "C:/Program Files/Python3!MINOR!/include/" -I C:/Users/Mary/AppData/Roaming/Python/Python3!MINOR!/site-packages/numpy/core/include/ -I "C:/opencv/build/include" -L "C:/opencv/build/" -L "C:/Program Files/Python3!MINOR!/libs" -lpython3!MINOR!
copy /Y forcedOsc.exe "C:\Users\Mary\Downloads\imagick_php_apache\webserver\Apache24\htdocs"
echo.
echo Running forcedOsc
forcedOsc.exe 1 1
@REM echo.
@REM echo.
@REM echo Compiling freeOsc
@REM g++ freeOsc.cpp -o freeOsc -I C:/msys64/mingw64/include -L/msys64/mingw64/lib -lgsl -lgslcblas -I "C:/Program Files/Python3!MINOR!/include/" -I C:/Users/Mary/AppData/Roaming/Python/Python3!MINOR!/site-packages/numpy/core/include/ -I "C:/opencv/build/include" -L "C:/opencv/build/" -L "C:/Program Files/Python3!MINOR!/libs" -lpython3!MINOR!
@REM echo.
@REM echo Running freeOsc
@REM freeOsc.exe
echo.
echo.
echo End
pause