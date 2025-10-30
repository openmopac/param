echo off
call "C:\Program Files (x86)\IntelSWTools\parallel_studio_xe_2020.0.075\bin\psxevars.bat" intel64 > nul
cd /D M:\PARAM\Source Code\MOPAC and PARAM source code
m:\Utility\Make_Date_Stamp.exe
nmake -f Makefile_for_PARAM_Windows.txt 

