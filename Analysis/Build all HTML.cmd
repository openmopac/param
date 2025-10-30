REM 
REM  Batch command file to build all the files in HTML
REM
REM goto here
CALL "M:\PARAM\Source Code\Bin and cmd\update_permanent_HTML_files.cmd"
CALL "M:\PARAM\Analysis\Analyze_Co-Crystals.cmd"
CALL "M:\PARAM\Analysis\Analyze_intermolecular_interactions.cmd"
CALL "M:\PARAM\Analysis\Analysis of solids\Small sets of Solids.cmd
CALL "M:\PARAM\Analysis\Analysis of solids\Map all Solids.cmd
CALL "M:\PARAM\Source Code\Bin and cmd\build_protein_HTML_files.cmd
:here
CALL "M:\PARAM\Source Code\Bin and cmd\Build molecule HTML files.cmd"
CALL "M:\PARAM\Source Code\Bin and cmd\Build intermolecular interaction HTML files.cmd"
CALL "M:\PARAM\Analysis\Heats of Sublimation.cmd"
