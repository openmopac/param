echo off
REM
REM In order to keep the PARAM ZIP file as small as possible, many of the 
REM files that are used in parameter optimization were deliberately deleted
REM before the ZIP file was constructed.
REM
REM This command file will reconstruct all the files that were deleted.
REM
REM  Set up all the local commands that will be used in PARAM
REM
PATH %PARAM_Start%\Source code\Bin and cmd;%PATH%
REM
REM Run the statistical analyses
REM
cd %PARAM_start%/Analysis/Survey organic
echo Preparing Analysis for simple organic compounds
call Sets_of_elements.cmd
REM
REM Now move to the various reference data folders, and construct the local "files.txt",
REM rename it, then create maps of the various sets of data. 
REM
REM                                                 Normal compounds
REM
echo Building JSmol files for normal compounds
cd "%PARAM_start%/Reference data/Data normal"
call "files.cmd"
cd "%PARAM_start%/Reference data/"
copy Header_for_molecules.txt + files.txt  aa.txt
mv aa.txt files.txt
call call Set_Up_Map.exe "data normal.txt" "Data normal" 


