@REM Arguments:
@REM     - none
@REM
@REM Example:
@REM     CreateVCProjAuto
@REM
@REM Creates a Visual Studio project file (vcproj) for the Qt project file
@REM in the current directory.


@ECHO OFF
set CURDIR=%CD%
FOR %%A in (*.mcg) DO ECHO %CD%/%%A 
FOR %%A in (*.mcg) DO CALL MapleCodegenCompleter.exe %%A false

:end