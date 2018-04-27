#!/usr/bin/env python

import os
import numpy as np
import re
import shutil

class aStructure:
    def __init__(self,dirname,libname,ref,species):
        "ref - reference energy of pure compounds"\
        "species - species read from gulp lib file."
        self.scale = []
        self.cell = []
        self.atomnums = {}
        self.atompos = []
        self.atomnames = []
        self.dftref = {}
        self.ffref = {}
        self.dftformation = 0
        self.ffformation = 0
        self.minusformation = 0
        self.postype = ''
        self.species = species
        self.lib = libname
        self.dirname = dirname
        for element in ref:
            self.dftref[element[0]]=element[1]
            self.ffref[element[0]]=element[2]

    def ReadFromFile(self,filename='POSCAR'):
        "read POSCAR by default"
        try:
            fin = open(filename,'r')
        except IOError,e:
            print("Unable to open file: "+filename)
        buf = fin.readlines()
        self.scale = float(buf[1][0])
        self.cell = [[float(coord)*self.scale for coord in line.split()] for line in buf[2:5]]        
        self.postype = buf[6]
        self.atompos = [[line.split()[3]]+map(float,line.split()[0:3]) for line in buf[7:]]
        for atom in self.atompos:
            if atom[0] not in self.atomnames:
                self.atomnames.append(atom[0])
        fin.close()
        for i in range(len(self.atomnames)):
            self.atomnums[self.atomnames[i]] = int(buf[5].split()[i])

    def FormatizeCoord(self):
        coord_out_core = []
        coord_out_shel = []
        for atom in self.atompos:
            types_for_atom = []
            for specie in self.species:
                if specie[0] == atom[0]:
                    types_for_atom.append(specie[1])
            pat1 = re.compile(r'^(C|c).*')
            pat2 = re.compile(r'^(S|s).*')
            "Can be revised for bshels in the future."
            for aType in types_for_atom:
                if pat1.match(aType):
                    coord_out_core.append([atom[0],aType]+['{:.8f}'.format(num) for num in atom[1:]])
                if pat2.match(aType):
                    coord_out_shel.append([atom[0],aType]+['{:.8f}'.format(num) for num in atom[1:]])
        coord_out = coord_out_core + coord_out_shel
        coord_out = [' '.join(line) for line in coord_out]
        coord_out = '\n'.join(coord_out)
        return coord_out

    def FormatizeSpecies(self):
        "Warning: Make Sure that lib file is already modified so that atom names are the same with specie names."
        species_out = []
        for specie in self.species:
            species_out.append(specie+[specie[0]])
        species_out = [' '.join(line) for line in species_out]
        species_out = 'species\n'+'\n'.join(species_out)
        return species_out
            
    def ReadVasp(self):
        try:
            fin = open('energy')
        except IOError,e:
            print "Can't open DFT energy file."
        dftene = float(fin.readlines()[0])
        self.dftformation = dftene
        fin.close()
        for atom in self.dftref:
            self.dftformation = self.dftformation - self.atomnums[atom]*self.dftref[atom]
        print 'DFT energy for '+self.dirname+': ',dftene
        print 'DFT formation energy for '+self.dirname+': ',self.dftformation
        fout = open('dftformation','w')
        fout.write('{:.8f}'.format(self.dftformation))
        fout.close()

    def WriteToGulp(self,opt_condstring='1 1 1 1 1 1'):
        "condstring is usually set to 111111 when space = 1"
        "Warning: Remember to write all cores before all shels or the .gin file won't read!And do not forget the EOF\n mark!"
        fout = open(self.dirname+'.gin','w')
        fout.write('opti conp prop \ntitle \nend\nvectors\n')
        cell_out = [["%.8f" % coord for coord in vec] for vec in self.cell]
        cell_out = ['  '.join(vec) for vec in cell_out]
        cell_out = '\n'.join(cell_out)+'\n'
        print "Writing cell parameters:\n"+cell_out
        fout.write(cell_out)
        fout.write(opt_condstring+'\n')
        pattern = re.compile(r'^(C|c).*')
        coortype = 'Cart' if pattern.match(self.postype) else 'frac'
        print "Coordinate type:\n"+coortype
        fout.write(coortype+'\n')
        fout.write(self.FormatizeCoord()+'\nspace\n1\n')
        fout.write(self.FormatizeSpecies()+'\nlibrary '+self.lib+'\noutput cif '+dirname+'.cif\n')
        fout.close()

    def ExecGulp(self):
        try:
            os.system('mpirun -np 4 gulp < '+self.dirname+'.gin > '+self.dirname+'.gout')
        except e:
            print "Gulp execution error, please check your input files!"
        fin = open(self.dirname+'.gout','r')
        buf = fin.readlines()
        fin.close()
        pattern1 = re.compile(r'.*Final energy.*')
        pattern2 = re.compile(r'.*Final Gnorm.*')
        for line in buf:
            if pattern1.match(line):
                gulpene = float(line.split()[3])
            if pattern2.match(line):
                gnorm = float(line.split()[3])
        if gnorm > 0.01:
            print "Warning: Gnorm for this structure is too big, check the gulp calculation output!"
        self.ffformation = gulpene
        for atom in self.ffref:
            self.ffformation = self.ffformation - self.atomnums[atom]*self.ffref[atom]
        print 'Gulp FF energy for '+self.dirname+': ',gulpene
        print 'Gulp FF formation energy for '+self.dirname+': ',self.ffformation
        fout = open('ffformation','w')
        fout.write('{:.8f}'.format(self.ffformation))
        fout.close()

    def ExecVasp(self):
        "Attention: 1,not modified for cluster running on TH-tj; 2,Currently not useable."
        try:
            os.system('mpirun -np 4 vasp_std > vasp.out')
        except e:
            print "Vasp execution error, please check your input files!"
        fin = open('OUTCAR','r')
        buf = fin.readlines()
        fin.close()
        pattern1 = re.compile(r'.*Energy without entropy.*')
        for line in buf:
            if pattern1.match(line):
                vaspene = float(line.split()[3])
        self.dftformation = gulpene
        for atom in self.dftref:
            self.dftformation = self.dftformation - self.atomnums[atom]*self.dftref[atom]
        print 'DFT energy for '+self.dirname+': ',vaspene
        print 'DFT formation energy for '+self.dirname+': ',self.dftformation
        fout = open('dftformation','w')
        fout.write('{:.8f}'.format(self.dftformation))
        fout.close()

    def DFTminusFF_form(self):
        self.minusformation = self.dftformation - self.ffformation
        print 'DFT-FF formation energy for '+self.dirname+': ',self.minusformation
        fout = open('dft-ffformation','w')
        fout.write('{:.8f}'.format(self.minusformation))
        fout.close()

    def RunStruct(self,mode='readvasp'):
        print '\nRunning structure: '+self.dirname
        self.ReadFromFile()
        if mode == 'readvasp':
            self.ReadVasp()
        else:
            slef.ExecVasp()
        self.WriteToGulp()
        self.ExecGulp()
        self.DFTminusFF_form()

def ReadRef(filename='refstates'):
    "The first row is reference pure compund name by one of its element, the 2nd row of file 'refstates' is PBE energy, the 3rd is gulp energy.Do not leave a blank line!"
    fin = open(filename,'r')
    buf = fin.readlines()
    buf = [line.split() for line in buf]
    buf = map(list,zip(*buf))
    buf = [[element[0],float(element[1]),float(element[2])] for element in buf]
    fin.close()
    print "Reference states(eV): ",buf
    return buf

def ReadSpecies(filename):
    fin = open(filename,'r')
    buf = fin.readlines()
    pattern1 = re.compile(r'^species.*')
    pattern2 = re.compile(r'^\D*$')
    marker = 0
    index_begin = 0
    index_end = 0
    for line in buf:
        if marker == 2:
            break
        if pattern1.match(line) and marker == 0:
            index_begin = buf.index(line)
            marker = 1
            continue
        if pattern2.match(line):
            if marker == 1:
                index_end = buf.index(line)
                marker = 2
    spc = buf[index_begin+1:index_end]
    print "Species in system: \n"+''.join(spc)
    spc = [line.split()[0:2] for line in spc]
    return spc



print "Warning:\n"+\
"1,Currently only support ATAT generated POSCAR.relax type;\n"+\
"2,Make sure atoms are properly ordered according to atom counts before use!\n"+\
"3,Make sure your reference compound energy is divided to eV per cation!\n"+\
"4,Make Sure that lib file is already modified so that atom names are the same with specie names.\n"
ref = ReadRef()
libname = os.path.basename(os.getcwd())+'.lib'
print "Detected library: "+libname+'\n'
species = ReadSpecies(libname)
for dirname in os.listdir('.'):
    if os.path.isdir(dirname):
        os.chdir(dirname)
        shutil.copyfile('../'+libname,'./'+libname)
        struct = aStructure(dirname,libname,ref,species)
        struct.RunStruct()
        os.chdir('..')

