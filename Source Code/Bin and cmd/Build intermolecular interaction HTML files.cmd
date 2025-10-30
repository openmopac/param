echo on
REM
REM Build a list of non-covalent interaction pairs of molecules
REM
REM  Extend the path to allow all the local commands that will be used.
REM
PATH %PARAM_Start%\Source code\Bin and cmd;%PATH%
REM
REM Go to the folders that contains MOP data-sets for the non-covalently bonded pairs of molecules
REM and list the .mop files in numerical order within each folder.
REM
cd %PARAM_start%\Reference Data\Data intermolecular interactions\S22 Intermolecular          & call ls.exe *.mop >  "../all.txt"
cd %PARAM_start%\Reference Data\Data intermolecular interactions\S66 Intermolecular          & call ls.exe *.mop >> "../all.txt"
cd %PARAM_start%\Reference Data\Data intermolecular interactions\L7 Intermolecular           & call ls.exe *.mop >> "../all.txt"
cd %PARAM_start%\Reference Data\Data intermolecular interactions\Water dimers Intermolecular & call ls.exe *.mop >> "../all.txt"
cd %PARAM_start%\Reference Data\Data intermolecular interactions\S12L Intermolecular         & call ls.exe *.mop >> "../all.txt"
cd %PARAM_start%\Reference Data\Data intermolecular interactions\X40 Intermolecular          & call ls.exe *.mop >> "../all.txt"
cd %PARAM_start%\Reference Data\Data intermolecular interactions\D3H4 Ionic Intermolecular   & call ls.exe *.mop >> "../all.txt"
cd %PARAM_start%\Reference Data\Data intermolecular interactions\                            & move /Y all.txt "./all intermolecular"
REM
REM Use this list to generate an un-sorted list, bits.txt
REM
cd %PARAM_start%\Reference Data\Data intermolecular interactions\all intermolecular
call files_program.exe NOSORT
move /Y all.txt "../"
REM
REM                                                 Intermolecular Interaction compounds
REM
cd "%PARAM_start%/Reference data/Data Intermolecular interactions/"
copy "Header for all intermolecular.txt" + files.txt  aa.txt
mv aa.txt files.txt
"M:/param/source code/bin and cmd/set_up_map.exe" "All Intermolecular interactions.txt" "All Intermolecular"


