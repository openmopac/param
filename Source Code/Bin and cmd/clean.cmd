echo off
echo *
echo * Unwanted files from "Bin and cmd" deleted.
echo *
cd %PARAM_Start%/Source code/Bin and cmd
 del "*.lib" 2> NUL
 del "*exe.intermediate.manifest" 2> NUL
 del "*exe.embed.manifest.res" 2> NUL
 del "*exe.embed.manifest" 2> NUL
 del "*.exp" 2> NUL
 del "*.rc" 2> NUL
 del "*.pdb" 2> NUL
IF EXIST M:\utility  (
cp MOPAC2016.exe "M:\utility\Windows_Standalone_64_bit_MOPAC2016.exe"
)
