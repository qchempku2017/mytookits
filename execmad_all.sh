#!/bin/bash

for  dir in `ls .`
do
	if [[ -d $dir && -f $dir/CONTCAR.static ]]
	then
		cd $dir
		rm -rf madgen
		mkdir madgen
		cp CONTCAR.static madgen/madgen
		cd madgen
		py_conv -f madgen -i vasp -o w2k	
		py_initmad
	    	calcmad madgen.inmad > mad.out
		echo $dir+"done!"
		cd ..
		cd ..
	else
	        echo "Failed to process:"+$dir+",no CONTCAR file found in the directory or not a direcotory!"
	fi
done
