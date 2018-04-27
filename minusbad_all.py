#!/usr/bin/env python

import os

for foldername in os.listdir('.'):
	if os.path.isdir(foldername):
                print "process: "+foldername
		os.chdir(foldername)
                if os.path.isdir('badgen'):
	        	f1 = open("badgen/badgen.outputmad",'r')
	        	f2 = open("energy",'r')
	        	f3 = open("energy-bad",'w')
	        	buf1 = f1.readline()
	            	buf2 = f2.readline()
	        	f1.close()
	        	f2.close()
	        	buf1 = [float(num) for num in buf1.split()]
	                buf2 = [float(num) for num in buf2.split()]
	                ene_bad = str(buf2[0]-buf1[1])
	                f3.write(ene_bad)
		        print "DFT energy without bader for "+foldername+":"+ene_bad
	                f3.close()
	        os.chdir("..")

