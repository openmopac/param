REM
REM         Set up soft links to directories where JSmol is expected
REM
REM In this example, a HTML file in folder "PARAM\Level 1\Level 2\Level3\Level4\Level 5"
REM would be able to use the JSmol folder in "PARAM" via a soft link to jsmol in
REM "PARAM\Level 1\Level 2\Level3\"
REM
echo on
call home
mklink /D  "Level 1\Level 2\Level 3\jsmol" %PARAM_Start%\jsmol
pause
echo off