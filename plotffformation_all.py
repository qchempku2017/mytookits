#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
import numpy as np

ffs_dfts = []
for name in os.listdir('.'):
    if os.path.isdir(name):
        os.chdir(name)
        if os.path.isfile('dftformation') and os.path.isfile('ffformation'):
            f_ff = open('ffformation','r')
            buf1 = f_ff.readline()
            f_ff.close()
            buf1 = [float(num) for num in buf1.split()]
            f_dft = open('dftformation','r')
            buf2 = f_dft.readline()
            f_dft.close()
            buf2 = [float(num) for num in buf2.split()]
            ffs_dfts.append([buf1[0],buf2[0]])
        os.chdir('..')
ffs = [ff_dft[0] for ff_dft in ffs_dfts]
dfts = [ff_dft[1] for ff_dft in ffs_dfts]
scale = (sum(ffs)/len(ffs))/(sum(dfts)/len(dfts))
#dfts = [dft*scale for dft in dfts]
plt.scatter(ffs,dfts,c='g',marker ='>')
#xr = np.arange(min(ffs),max(ffs),(max(ffs)-min(ffs))/50)
#plt.plot(xr,xr,color='r')
plt.xlabel('Formation energy calculated with force field method/eV')
plt.ylabel('Formation energy calculated with PBESol/eV')
plt.legend()
plt.tight_layout()
plt.savefig('ff_DFT.png',dpi=300)
