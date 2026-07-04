@echo off
REM scStudio launcher for Windows. Double-click this file.
REM Requires R installed and on PATH, with the scStudio package installed.
Rscript -e "scStudio::run_app()"
pause
