echo on
del "%PARAM_start%\HTML Files\Statistics for Intermolecular interactions.html" 2> nul
call "%PARAM_start%\Source Code\Bin and cmd\stats" "%PARAM_start%\Analysis\S22\S22"
call "%PARAM_start%\Source Code\Bin and cmd\stats" "%PARAM_start%\Analysis\S12L\S12L"
call "%PARAM_start%\Source Code\Bin and cmd\stats" "%PARAM_start%\Analysis\S66\S66"
call "%PARAM_start%\Source Code\Bin and cmd\stats" "%PARAM_start%\Analysis\L7\L7"
call "%PARAM_start%\Source Code\Bin and cmd\stats" "%PARAM_start%\Analysis\Water dimers\Dimers of Water"
call "%PARAM_start%\Source Code\Bin and cmd\stats" "%PARAM_start%\Analysis\X40\X40"
call "%PARAM_start%\Source Code\Bin and cmd\stats" "%PARAM_start%\Analysis\D3H4 Ionic\D3H4 Ionic"
call "%PARAM_start%\Source Code\Bin and cmd\stats" "%PARAM_start%\Analysis\All intermolecular interactions\All"


