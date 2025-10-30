REM
REM  Construct a list of MOPAC data-sets, in order of the formula
REM
if exist all.txt del all.txt
dir /b *.mop > all.txt 2> NUL
call "%PARAM_Start%\Source code\Bin and cmd\files_program.exe"
if exist all.txt del all.txt


