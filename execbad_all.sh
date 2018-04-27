#!/bin/bash

for  dir in `ls .`
do
	if [[ -d $dir && -f $dir/POSCAR.relax ]]
	then
		cd $dir
		rm -rf badgen
		mkdir badgen
		cp str.out badgen/badgen
		cp ../BaderDict badgen/BaderDict
		cd badgen
		py_initbad -i atat -o badgen.inmad badgen
	    	calcmad badgen.inmad > bad.out
		echo $dir+"done!"
		cd ..
		cd ..
	else
	        echo "Failed to process:"+$dir+",no str.out file found in the directory or not a direcotory!"
	fi
done
