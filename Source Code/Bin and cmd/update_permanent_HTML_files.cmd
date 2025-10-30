REM
REM  Make a copy of the permanent files.
REM  If any permanent file has been deleted, restore it.
REM
set from=M:\PARAM\HTML Files
set to=M:\PARAM\HTML Permanent Files
CALL xcopy "%from%\Densit3.gif" "%to%\"   /I /S /C /D /Y /H
CALL xcopy "%from%\Densit4.gif" "%to%\"   /I /S /C /D /Y /H
CALL xcopy "%from%\Densities.html" "%to%\"   /I /S /C /D /Y /H
CALL xcopy "%from%\Heats_of_Formation.html" "%to%\"   /I /S /C /D /Y /H
CALL xcopy "%from%\Known faults in protein modeling.html" "%to%\"   /I /S /C /D /Y /H
CALL xcopy "%from%\Note_on_all_elements.html" "%to%\"   /I /S /C /D /Y /H
CALL xcopy "%from%\Notes_on_Proteins.html" "%to%\"   /I /S /C /D /Y /H
CALL xcopy "%from%\Periodic_table_solids.html" "%to%\"   /I /S /C /D /Y /H
CALL xcopy "%from%\Polarizabilities.html" "%to%\"   /I /S /C /D /Y /H
CALL xcopy "%from%\Accuracy of PM7 and PM6-D3H4.html" "%to%\"   /I /S /C /D /Y /H
CALL xcopy "%from%\PM7 and PM6-D3H4 Notes.html" "%to%\"   /I /S /C /D /Y /H
CALL xcopy "%from%\PM6-D3H4_PM6-D3H4\3CLpro inhibitor (4MDS).txt" "%to%\PM6-D3H4_PM6-D3H4\"   /I /S /C /D /Y /H
CALL xcopy "%from%\PM6-D3H4_PDB\3CLpro inhibitor (4MDS).txt" "%to%\PM6-D3H4_PDB\"   /I /S /C /D /Y /H
CALL xcopy "M:\PARAM\HTML Files\PM7_PDB\*.jpg" "M:\PARAM\HTML Permanent Files\PM7_PDB\"   /I /S /C /D /Y /H
CALL xcopy "M:\PARAM\HTML Files\PM7_PDB\*.cdxml" "M:\PARAM\HTML Permanent Files\PM7_PDB\"   /I /S /C /D /Y /H
CALL xcopy "M:\PARAM\HTML Files\PM7_PDB\Notes on*.html" "M:\PARAM\HTML Permanent Files\PM7_PDB\"   /I /S /C /D /Y /H
FOR  %%p IN (PM7_PDB, PM7_10, PM7_3, PM7_PM7, PM6-D3H4_PDB, PM6-D3H4_10, PM6-D3H4_3, PM6-D3H4_PM6-D3H4) DO (
CALL xcopy "M:\PARAM\HTML Files\%%p\*.txt" "M:\PARAM\HTML Permanent Files\%%p\"   /I /S /C /D /Y /H
)
REM
CALL xcopy "M:\PARAM\HTML Permanent Files\*" "M:\PARAM\HTML Files\"   /I /S /C /D /Y /H