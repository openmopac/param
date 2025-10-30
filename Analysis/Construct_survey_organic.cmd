REM
REM  Calculate properties of main-group elements for several methods
REM
cd "%PARAM_start%\Analysis\Survey organic" 
FOR  %%p IN (PM7, PM6-D3H4, PM6-D3H4X, PM6-DH+, PM6-DH2, PM6-DH2X, PM6, RM1) DO (
echo Starting %%p 
CALL "%PARAM_start%\Source Code\Bin and cmd\PARAM.exe" %%p
)
rm *.out
cd ../

