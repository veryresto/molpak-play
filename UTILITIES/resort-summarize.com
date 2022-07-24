#!/bin/csh
setenv dir /PREDICTIONS
goto JUMP
set f='AF BH BD BF AH AY AZ AV AB AQ AI CC DC AK CA AM FA FC DA DB DD DE CB AA AP BA BB CD CE AU AS'
set d='af bh bd bf ah ay az av ab aq ai cc dc ak ca am fa fc da db dd de cb aa ap ba bb cd ce au as'
set dir1 = ($d)
set fl  = ($f)
set n=1
while ($n != 32)
 if ( -d $dir1[$n] ) then
   cd $dir1[$n]
   echo "current job number: $n "
   echo "current directory : $dir1[$n]"
   mv *.tab $fl[$n].tab
   cd ..
 endif
@ n = ($n + 1)
end
JUMP:
#if (-e densities.list) goto NEXT
$dir/UTILITIES/summarize_tab.com > densities.list
NEXT:
$dir/UTILITIES/resort-summarize-tab.exe << EOF
densities.list
EOF
exit
