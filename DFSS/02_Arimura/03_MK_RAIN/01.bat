@echo off
SET OSGEO4W_ROOT=C:\OSGeo4W
call %OSGEO4W_ROOT%\bin\o4w_env.bat

SET PYTHONHOME=%OSGEO4W_ROOT%\apps\Python39
PATH %PYTHONHOME%\Scripts;%PATH%

python3 "%~dp001_mk_fv.py" %*

