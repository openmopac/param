echo off
REM
REM Build a list of solids, ordered according to the solid's position
REM in the periodic table.
REM
REM  Extend the path to allow all the local commands that will be used.
REM
PATH %PARAM_Start%\Source code\Bin and cmd;%PATH%
REM
REM Go to a folder that contains a complete set of ARC files for solids
REM
cd %PARAM_start%\Analysis\Analysis of solids\PM7
REM
REM Make a list of all the ARC files, unordered, in all.txt
REM
call ls.exe *.arc > all.txt
REM
REM Use this list to generate a sorted list, bits.txt
REM
call files_program.exe SOLID
REM
REM Build the file "Survey_of_Solids.txt" in the next-higher folder
REM
type  "Header for list of solids.txt"  bits.txt > ../Survey_of_Solids.txt
del bits.txt
del files.out
del all.txt

