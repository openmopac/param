echo on
REM
REM   Copying PM7 files
REM
                   set from_PM7=M:\PARAM\Analysis\Proteins
                   set from_PM6=M:\PARAM\Analysis\Proteins
                   set to_PM7=M:\PARAM\HTML Files
                   set to_PM6=M:\PARAM\HTML Files
REM
FOR %%P IN (PM7_PDB,PM7_PM7,PM7_3,PM7_10) DO CALL xcopy %from_PM7%\%%P\*.arc  "%to_PM7%\%%P\"  /I /S /C /D /Y
FOR %%P IN (PM7_PDB,PM7_PM7,PM7_3,PM7_10) DO CALL xcopy %from_PM7%\%%P\*.pdb  "%to_PM7%\%%P\"  /I /S /C /D /Y
FOR %%P IN (PM7_PDB,PM7_PM7,PM7_3,PM7_10) DO CALL xcopy %from_PM7%\%%P\*.html "%to_PM7%\%%P\"  /I /S /C /D /Y
FOR %%P IN (PM7_PDB,PM7_PM7,PM7_3,PM7_10) DO CALL xcopy %from_PM7%\%%P\*.txt  "%to_PM7%\%%P\"  /I /S /C /D /Y
REM FOR %%P IN (PM6_PDB,PM6_opt,PM6_3,PM6_10) DO CALL xcopy %from_PM6%\%%P\*.arc "%to_PM6%\%%P\"  /I /S /C /D /Y
REM FOR %%P IN (PM6_PDB,PM6_opt,PM6_3,PM6_10) DO CALL xcopy %from_PM6%\%%P\*.pdb "%to_PM6%\%%P\"  /I /S /C /D /Y
