#!/bin/bash

for dir in `ls .`
do 
	if [[ -d $dir && -f  $dir/INCAR && -f $dir/POSCAR && -f $dir/POTCAR && -f $dir/KPOINTS ]]
	then
		cd $dir
		rm -rf energy
		mpirun -np 4 vasp_std
		if [[ -f OUTCAR ]]
		then
			cat OUTCAR|grep without|tail -1|awk '{print $4}' > energy
			echo $dir"done"
		else 
			echo $dir" VASP failed."
		fi
		cd ..
	fi
done
