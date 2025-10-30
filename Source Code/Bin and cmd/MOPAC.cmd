if  %2x == x goto okay
echo off
echo *
echo *      The path or filename has been split into parts.
echo *      This typically happens if either the path to the folder or the filename contains a comma.
echo *      To correct this, make sure that neither the folder nor the filename contain a comma.
timeout 1000
:okay
start "Low cmd" /low /b "%PARAM_Start%\Source code\Bin and cmd\MOPAC2016.exe" %1
