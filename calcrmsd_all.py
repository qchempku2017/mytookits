#!/usr/bin/env python

import os
import numpy as np
import re
import math
import matplotlib.pyplot as plt

class aCIFStructure:
    "1,Currently only works for P1 space group.\
     2,Cannot read parameters inserted between loop_'s\
     3,CIF Files should be cut rigorously into 2 parts.\
     4.Loops are separated by empty lines. No emptylines should apprear between loops."
    def __init__(self):
        self.celldict = {}
        self.symGroup = 0
        self.loopdict = {}
        self.loopkeylist = []
        self.cellmat = np.matrix([])
        self.posmat = np.matrix([])
        self.filecutline = 0

    def Read_Cell_Sym(self,filebuf):
        pat_cell = re.compile(r'^_cell.*')
        pat_sym = re.compile(r'^_symmetry_Int_Table_number.*')
        pat_loop = re.compile(r'^loop_.*')
        for i in range(len(filebuf)):
            if pat_cell.match(filebuf[i]):
                buf = filebuf[i].split()
                buf = [buf[0],float(buf[1])]
                self.celldict[buf[0]]=buf[1]
            if pat_sym.match(filebuf[i]):
                buf = filebuf[i].split()
                self.symGroup = int(buf[1])
            if pat_loop.match(filebuf[i]):
                break
        self.filecutline = i

    def Read_Loopdict(self,filebuf):
        pat_key = re.compile(r'^_.*')
        pat_loop = re.compile(r'^loop_.*')
        pat_emptyline = re.compile(r'\n$')
        valuebuf = []
        loopkeylist = []
        #print len(filebuf)
        #print filebuf
        for i in range(self.filecutline,len(filebuf)):
            if pat_loop.match(filebuf[i]):
                continue
            else:
                if pat_emptyline.match(filebuf[i]) or i == len(filebuf)-1:
                    #print loopkeylist,valuebuf
                    for j in range(len(loopkeylist)):
                        self.loopdict[loopkeylist[j]] = [line[j] for line in valuebuf]
                    valuebuf = []
                    loopkeylist = []
                else:
                    if pat_key.match(filebuf[i]):
                        self.loopkeylist.append(filebuf[i])
                        loopkeylist.append(filebuf[i])
                    else:
                        valuebuf.append(filebuf[i].split())
        #print self.loopdict,self.loopkeylist

    def Conv_Celldict_To_Cellmat(self):
        a = self.celldict['_cell_length_a']
        b = self.celldict['_cell_length_b']
        c = self.celldict['_cell_length_c']
        alpha = self.celldict['_cell_angle_alpha']/180.0*math.pi
        beta = self.celldict['_cell_angle_beta']/180.0*math.pi
        gamma = self.celldict['_cell_angle_gamma']/180.0*math.pi
        vec_a = np.array([a,0,0])
        vec_b = np.array([math.cos(gamma),math.sin(gamma),0])*b
        cx = c*math.cos(beta)
        cy = c*(math.cos(alpha)-math.cos(beta)*math.cos(gamma))/math.sin(gamma)
        cz = math.sqrt(c**2 - cx**2 - cy**2)
        vec_c = np.array([cx,cy,cz])
        self.cellmat = np.matrix(np.vstack((vec_a,vec_b,vec_c)))
        print 'Cell parameters in vector form: ', self.cellmat

    def Conv_Loopdict_To_Posmat(self):
        pat_x = re.compile(r'^_atom_site_.*_x.*')
        pat_y = re.compile(r'^_atom_site_.*_y.*')
        pat_z = re.compile(r'^_atom_site_.*_z.*')
        pat_frac = re.compile(r'.*fract.*')
        is_frac= True
        pos_x = []
        pos_y = []
        pos_z = []
        for key in self.loopdict:
            if pat_x.match(key):
                pos_x = [float(x) for x in self.loopdict[key]]
                is_frac = is_frac and pat_frac.match(key)
            elif pat_y.match(key):
                pos_y = [float(y) for y in self.loopdict[key]]
                is_frac = is_frac and pat_frac.match(key)
            elif pat_z.match(key):
                pos_z = [float(z) for z in self.loopdict[key]]
                is_frac = is_frac and pat_frac.match(key)
        poslist = zip(pos_x,pos_y,pos_z)
        self.posmat = np.matrix(poslist)
        if is_frac:
            self.posmat = self.posmat * self.cellmat
        print "Atom coordinates converted!"

    def Read_From_File(self,filename):
        try:
            f = open(filename,'r')
        except e,IOError:
            print "Unable to open: "+filename
        buf = f.readlines()
        f.close()
        self.Read_Cell_Sym(buf)
        self.Read_Loopdict(buf)
        self.Conv_Celldict_To_Cellmat()
        self.Conv_Loopdict_To_Posmat()


def Calc_Rmsd(cif1,cif2):
    "In previous operations we have aligned the 0 point of the 2 structures to the same.\
     Indices of atoms are not changed during optimisation and cif convertion."
    rms = 0
    if cif1.posmat.shape == cif2.posmat.shape:
        natom = cif1.posmat.shape[0]
    else:
        print "Error: Atom number of compared structures don't match!"
    for i in range(natom):
        d = np.linalg.norm(cif1.posmat[i]-cif2.posmat[i])
        for j in range(3):
            d = min(d,np.linalg.norm(cif1.posmat[i]-cif2.posmat[i]- cif2.cellmat[j]),np.linalg.norm(cif1.posmat[i]-cif2.posmat[i]+cif2.cellmat[j]))
        rms = rms + d**2
    rmsd = math.sqrt(rms/natom)
    return rmsd

dir_list = os.listdir('.')
d_rmsd = []
for d in dir_list:
    if os.path.isdir(d):
        print '\nDoing directory:',d
        if os.path.isfile(d+'/error'):
            print "Error found in"+d
            continue
        os.chdir(d)
        os.system('py_conv -f CONTCAR -i vasp -o cif')
        dft_str = aCIFStructure()
        ff_str = aCIFStructure()
        dft_str.Read_From_File('CONTCAR.cif')
        ff_str.Read_From_File(d+'.cif')
        rmsd = Calc_Rmsd(ff_str,dft_str) #The second parameter is the reference
        d_rmsd.append((d,rmsd))
        print "RMSD value(use DFT structure as reference) of "+d+": ",rmsd,'\n'
        os.chdir('..')

d_rmsd = sorted(d_rmsd,key=lambda x:x[0]) 
names = [pair[0] for pair in d_rmsd]
rmsd = [pair[1] for pair in d_rmsd]
f_out = open('rmsd.log','w')
buf_out = []
for pair in d_rmsd:
    buf_out.append('{0}  {1:.8}'.format(pair[0],pair[1]))
buf_out='\n'.join(buf_out)+'\n'
f_out.write(buf_out)
f_out.close()
plt.bar(names,rmsd)
filtered_names = [pair[0] for pair in d_rmsd if pair[1]>0.8]
plt.xticks(filtered_names,fontsize=3,rotation=90)
plt.savefig('RMSD.png',dpi=600)




        
