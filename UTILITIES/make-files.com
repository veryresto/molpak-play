# shell script name = make-files.com
#
# shell script to make MOLPAK/WMIN/PMIN/DMAREL files for crystal structure/
#   density prediction
#
setenv dir PREDICTIONS
ls *molpak.xyz > MOLPAK.NAME
$dir/UTILITIES/make_files.exe
chmod 755 *
#rm MOLPAK.NAME
#end
