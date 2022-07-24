
pwd > dir
./rename-dir.exe
mv -f fort.13 UTILITIES/make-files.com
mv -f fort.14 UTILITIES/summarize_tab.com
mv -f fort.15 UTILITIES/make_files.f90
mv -f fort.16 UTILITIES/volume-additivity-densities-4.f 
mv -f fort.17 UTILITIES/atom_code-vol-additivity.f
mv -f fort.18 UTILITIES/resort-summarize.com
rm -f dir
chmod 755 UTILITIES/*.com
exit

