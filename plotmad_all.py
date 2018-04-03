#!/usr/bin/env python

import os
import matplotlib.pyplot as plt

index_mads_dfts = []
for name in os.listdir('.'):
    if os.path.isdir(name):
        os.chdir(name)
        if os.path.isdir('madgen'):
            f_mad = open('madgen/madgen.outputmad','r')
            buf1 = f_mad.readline()
            f_mad.close()
            buf1 = [float(num) for num in buf1.split()]
            f_dft = open('energy','r')
            buf2 = f_dft.readline()
            f_dft.close()
            buf2 = [float(num) for num in buf2.split()]
            index_mads_dfts.append([int(name),buf1[1],buf2[0]])
        os.chdir('..')
index_mads_dfts = sorted(index_mads_dfts,key =lambda t:t[0])
index = [index_mad_dft[0] for index_mad_dft in index_mads_dfts]
mad = [index_mad_dft[1] for index_mad_dft in index_mads_dfts]
dft = [index_mad_dft[2] for index_mad_dft in index_mads_dfts]
scale = (sum(mad)/len(mad))/(sum(dft)/len(dft))
dft = [dfti*scale for dfti in dft]
plt.plot(index,mad,color='b',label='Madelung Energy')
plt.plot(index,dft,color='r',label='Scaled DFT Energy')
plt.xlabel('Structure Number')
plt.ylabel('Energy/eV')
plt.title('Madelung/Scaled DFT Energy vs. Structure diagram')
plt.legend()
plt.tight_layout()
plt.savefig('Madelung_ScaDFT.png',dpi=300)
