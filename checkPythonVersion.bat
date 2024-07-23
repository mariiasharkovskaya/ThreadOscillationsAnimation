@echo off
setlocal enabledelayedexpansion

SET MINOR=
for /f "tokens=2 delims=." %%a in ('python --version') do (
    set MINOR=%%a
)

g++ forcedOsc.cpp -o mygif -I C:/msys64/mingw64/include -L/msys64/mingw64/lib -lgsl -lgslcblas -I "C:/Program Files/Python3!MINOR!/include/" -I C:/Users/Mary/AppData/Roaming/Python/Python3!MINOR!/site-packages/numpy/core/include/ -I "C:/opencv/build/include" -L "C:/opencv/build/" -L "C:/Program Files/Python3!MINOR!/libs" -lpython3!MINOR!
mygif.exe
