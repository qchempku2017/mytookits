#!/usr/bin/env python

"ValDict File Format:\
    List of atom names\
    corresponding atom number in POSCAR\
    Number of valence electrons read from POTCAR\
Please set up ValDict before running this script!"

class aBaderDict():
    def __init__(self):
        self.elec = []
        self.atoms = []
        self.atomsbuf = []
        self.elec_avg = []

    def Read_from_File(self):
        print 'Reading Bader analysis result.'
        try:
            f = open('ACF.dat','r')
        except IOError,e:
	    print "Can't open file ACF.dat!"
        buf = f.readlines()[2:-4]
        f.close()
        buf = [[float(num) for num in line.split()] for line in buf]
        self.elec = [line[4] for line in buf]
        print 'Electron count List:\n',self.elec
        print 'Reading Atom list from vasp file'
        try:
            f2 = open('ValDict','r')
        except IOError,e:
            print "Can't open file ValDict"
        buf2 = f2.readlines()
        self.atomsbuf = [buf2[0].split(),[int(num) for num in buf2[1].split()],[float(val) for val in buf2[2].split()]]
        self.atoms = [ [self.atomsbuf[0][i],self.atomsbuf[1][i],self.atomsbuf[2][i]] for i in range(len(self.atomsbuf[0])) ]
        print "Atom names and numbers:",self.atoms
        f2.close()

    "Take Average Over All Coordination Environments of the same Atom type."
    def AvgElec(self):
        index = 0
        for atomType in self.atoms:
            atomNum = atomType[1]
            self.elec_avg.append(atomType[2]-sum(self.elec[index:index+atomNum])/float(atomNum))
            index = index + atomNum

    def Write_to_File(self):
        f = open('BaderDict','w')
        elec_avg_strs = ["{:10.8f}".format(num) for num in self.elec_avg]
        f.write('  '.join(self.atomsbuf[0]))
        f.write('\n')
        f.write('  '.join(elec_avg_strs))

print "Warning: Please set up ValDict before running this script!"
Dict = aBaderDict()
Dict.Read_from_File()
Dict.AvgElec()
Dict.Write_to_File()

