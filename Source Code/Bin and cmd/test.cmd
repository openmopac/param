REM
REM  Make a copy of the permanent files.
REM  If any permanent file has been deleted, restore it.
REM
set from=M:\PARAM\HTML Files
set to=M:\PARAM\HTML Permanent Files
FOR  %%p IN (PM7_PDB, PM7_10, PM7_3, PM7_PM7) DO (
CALL xcopy "M:\PARAM\HTML Files\%%p\*.txt" "M:\PARAM\HTML Permanent Files\%%p\"   /I /S /C /D /Y /H
)
REM
REM CALL xcopy "M:\PARAM\HTML Permanent Files\*" "M:\PARAM\HTML Files\"   /I /S /C /D /Y /H
pause