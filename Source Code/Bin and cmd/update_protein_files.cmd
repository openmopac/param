echo on
REM
                   set from=M:\PARAM\Analysis\Proteins
                   set to=M:\PARAM\HTML Files
REM
REM   Copying PM6-D3H4 files
REM
FOR %%P IN (PM6-D3H4_PDB,PM6-D3H4_PM6-D3H4,PM6-D3H4_3,PM6-D3H4_10) DO CALL xcopy %from%\%%P\*.arc  "%to%\%%P\"  /I /S /C /D /Y
FOR %%P IN (PM6-D3H4_PDB,PM6-D3H4_PM6-D3H4,PM6-D3H4_3,PM6-D3H4_10) DO CALL xcopy %from%\%%P\*.pdb  "%to%\%%P\"  /I /S /C /D /Y
FOR %%P IN (PM6-D3H4_PDB,PM6-D3H4_PM6-D3H4,PM6-D3H4_3,PM6-D3H4_10) DO CALL xcopy %from%\%%P\*.html "%to%\%%P\"  /I /S /C /D /Y
#
# Do NOT update .txt files.  Edit them directly in ~/PARAM/HTML/PM7_PDB etc
#
# FOR %%P IN (PM6-D3H4_PDB,PM6-D3H4_PM6-D3H4,PM6-D3H4_3,PM6-D3H4_10) DO CALL xcopy %from%\%%P\*.txt  "%to%\%%P\"  /I /S /C /D /Y
REM
REM   Copying PM7 files
REM
REM
FOR %%P IN (PM7_PDB,PM7_PM7,PM7_3,PM7_10) DO CALL xcopy %from%\%%P\*.arc  "%to%\%%P\"  /I /S /C /D /Y
FOR %%P IN (PM7_PDB,PM7_PM7,PM7_3,PM7_10) DO CALL xcopy %from%\%%P\*.pdb  "%to%\%%P\"  /I /S /C /D /Y
FOR %%P IN (PM7_PDB,PM7_PM7,PM7_3,PM7_10) DO CALL xcopy %from%\%%P\*.html "%to%\%%P\"  /I /S /C /D /Y
# FOR %%P IN (PM7_PDB,PM7_PM7,PM7_3,PM7_10) DO CALL xcopy %from%\%%P\*.txt  "%to%\%%P\"  /I /S /C /D /Y

