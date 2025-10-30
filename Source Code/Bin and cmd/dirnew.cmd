if %1x == x goto all
dir /o-d %1 | "%PARAM_Start%\Source Code\Bin and cmd\head.exe" -25
goto finish
:all
dir /o-d | "%PARAM_Start%\Source Code\Bin and cmd\head.exe" -25
:finish
