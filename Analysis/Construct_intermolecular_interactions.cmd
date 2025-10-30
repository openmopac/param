REM
REM  Calculate intermolecular interaction energies for several methods for several sets of systems
REM
ECHO OFF
FOR  %%p IN (S22, S66, L7, S12L, "Water dimers", "D3H4 Ionic", X40, "All intermolecular interactions") DO (
cd %%p
FOR  %%q IN (PM7, PM6-D3H4, PM6-D3H4X, PM6-DH+, PM6-DH2, PM6-DH2X, PM6) DO (
echo Starting %%p %%q
CALL "%PARAM_start%\Source Code\Bin and cmd\PARAM.exe" %%q
)
cd ../
)
REM rm */*.out
pause
