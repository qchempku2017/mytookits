#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
import numpy as np

print "Check rmsd.log and dat files in the directory,adjust rmsd and deviant criteria before use"
ffs_dfts_g = []
ffs_dfts_r = []
ffs_dfts_b = []
ffs_dfts_x = []
redlist = []
bluelist = []
xlist = []
if os.path.isfile('dftformation.dat'):
    f_lst = open('dftformation.dat','r')
    buf=[line.split() for line in f_lst.readlines()[1:]]
    buf=[[line[0],float(line[1]),float(line[2])] for line in buf]
    f_lst.close()
    for line in buf:
        if abs(line[2]-line[1])>0.10:
            bluelist.append(line[0])

if os.path.isfile('dft-ffformation.dat'):
    f_lst = open('dft-ffformation.dat','r')
    buf=[line.split() for line in f_lst.readlines()[1:]]
    buf=[[line[0],float(line[1]),float(line[2])] for line in buf]
    f_lst.close()
    for line in buf:
        if abs(line[2]-line[1])>0.10:
            xlist.append(line[0])

if os.path.isfile('rmsd.log'):
    f_lst = open('rmsd.log')
    buf=[line.split() for line in f_lst.readlines()]
    buf=[[line[0],float(line[1])] for line in buf]
    f_lst.close()
    for line in buf:
        if line[1]>0.25:
            redlist.append(line[0])

for name in os.listdir('.'):
    if os.path.isdir(name):
        os.chdir(name)
        if os.path.isfile('dftformation') and os.path.isfile('ffformation') and not(os.path.isfile('error')):
            f_ff = open('ffformation','r')
            buf1 = f_ff.readline()
            f_ff.close()
            buf1 = [float(num) for num in buf1.split()]
            f_dft = open('dftformation','r')
            buf2 = f_dft.readline()
            f_dft.close()
            buf2 = [float(num) for num in buf2.split()]
            ffs_dfts_g.append([buf1[0],buf2[0]])
            if name in redlist:
                ffs_dfts_r.append([buf1[0],buf2[0]])
            if name in bluelist:
                ffs_dfts_b.append([buf1[0],buf2[0]])
            if name in xlist:
                ffs_dfts_x.append([buf1[0],buf2[0]])
        os.chdir('..')

ffs_g = [ff_dft[0] for ff_dft in ffs_dfts_g]
dfts_g = [ff_dft[1] for ff_dft in ffs_dfts_g]
ffs_r = [ff_dft[0] for ff_dft in ffs_dfts_r]
dfts_r = [ff_dft[1] for ff_dft in ffs_dfts_r]
ffs_b = [ff_dft[0] for ff_dft in ffs_dfts_b]
dfts_b = [ff_dft[1] for ff_dft in ffs_dfts_b]
ffs_x = [ff_dft[0] for ff_dft in ffs_dfts_x]
dfts_x = [ff_dft[1] for ff_dft in ffs_dfts_x]
#print bluelist
#print redlist
#print xlist
#scale = (sum(ffs)/len(ffs))/(sum(dfts)/len(dfts))
#dfts = [dft*scale for dft in dfts]
plt.scatter(ffs_g,dfts_g,c='g',marker ='>')
plt.scatter(ffs_r,dfts_r,c='r',marker ='o',label='RMSD>0.25')
plt.scatter(ffs_b,dfts_b,c='b',marker ='D',label='DFT CE Outliers')
plt.scatter(ffs_x,dfts_x,c='k',marker ='x',label='DFT-FF CE Outliers')
x = np.arange(min(ffs_g),max(ffs_g),0.001)
plt.plot(x,x,c='k',label='y=x')
#xr = np.arange(min(ffs),max(ffs),(max(ffs)-min(ffs))/50)
#plt.plot(xr,xr,color='r')
plt.xlabel('Formation energy calculated with force field method/eV')
plt.ylabel('Formation energy calculated with PBESol/eV')
plt.legend()
plt.tight_layout()
plt.savefig('ff_DFT.png',dpi=300)
