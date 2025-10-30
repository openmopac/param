echo off
if exist %1 goto is
if %2x == x goto heat
tail -6000f %1.out | grep  %2
 find /i %2 %1.out
:heat
tail -2000f %1.out 
:is
if %2x == x goto all
tail -6000f %1 | grep  %2
 find /i %2 %1
:all
tail -2000f %1 
 
