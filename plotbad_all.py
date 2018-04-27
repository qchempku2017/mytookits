#!/usr/bin/env python

import os
import matplotlib.pyplot as plt

index_bads_dfts = []
for name in os.listdir('.'):
    if os.path.isdir(name):
        os.chdir(name)
        if os.path.isdir('badgen'):
            f_mad = open('badgen/badgen.outputmad','r')
            buf1 = f_mad.readline()
            f_mad.close()
            buf1 = [float(num) for num in buf1.split()]
            f_dft = open('energy','r')
            buf2 = f_dft.readline()
            f_dft.close()
            buf2 = [float(num) for num in buf2.split()]
            index_bads_dfts.append([int(name),buf1[1],buf2[0]])
        os.chdir('..')
index_bads_dfts = sorted(index_bads_dfts,key =lambda t:t[0])
index = [index_bad_dft[0] for index_bad_dft in index_bads_dfts]
bad = [index_bad_dft[1] for index_bad_dft in index_bads_dfts]
dft = [index_bad_dft[2] for index_bad_dft in index_bads_dfts]
scale = (sum(bad)/len(bad))/(sum(dft)/len(dft))
dft = [dfti*scale for dfti in dft]
plt.plot(index,bad,color='b',label='Bader Energy')
plt.plot(index,dft,color='r',label='Scaled DFT Energy,Scale '+str(scale))
plt.xlabel('Structure Number')
plt.ylabel('Energy/eV')
plt.title('Bader/Scaled DFT Energy vs. Structure diagram')
plt.legend()
plt.tight_layout()
plt.savefig('Bader_ScaDFT.png',dpi=300)
