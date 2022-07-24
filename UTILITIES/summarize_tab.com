#!/bin/csh
setenv dir /c/Users/ACER/Downloads/PREDICTIONS
ls */*tab > TAB_FILES      # directory names
chmod 744 TAB_FILES
$dir/UTILITIES/summarize_tab.exe
$dir/UTILITIES/volume-additivity-densities-4.exe
echo " ----------------------------------------------------------------"
$dir/UTILITIES/atom_code-vol-additivity.exe
rm -f TAB_FILES
