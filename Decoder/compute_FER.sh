#!/bin/bash
full_name="$1"
xpath=${full_name%/*} 
xbase=${full_name##*/}
xfext=${xbase##*.}
xpref=${xbase%.*}
#echo "path='${xpath}', pref='${xpref}', ext='${xfext}'"
fout=$xpref".fer"
#echo $fout
for ((i = 0 ; i < 700 ; i += 25));
#	do echo ./decoder $1 $i"E-02" 200 100000000 2	
	do ./decoder $1 $i"E-02" 200 100000000 2>> $fout
done
