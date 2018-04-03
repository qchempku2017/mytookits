#!/usr/bin/env python

import os
import numpy as np
import math
import cmath
import matplotlib.pyplot as plt

pi = math.pi

"Add more atoms to dict when needed."
atom_dict = { 'Na':1,'Cl':7,'Zn':12, 'S':6 ,'Ti':4 , 'O':6}

def read_elec():
    print 'Reading Bader analysis result.'
    if os.path.isfile('ACF.dat'):
        f = open('ACF.dat','r')
        buf = f.readlines()[2:-4]
        f.close()
        buf = [[float(num) for num in line.split()] for line in buf]
        elec = [line[4] for line in buf]
        print 'Electron count List:\n',elec
        return elec
    else:
        return False

def read_cell_atm_pos():
   print 'Reading structural information Warning: POSCAR must exist and is written in standard format'
   if os.path.isfile('POSCAR'):
       f = open('POSCAR','r')
       buf = f.readlines()
       f.close()
       buf_cell = buf[1:5]
       buf_atm = buf[5:7]
       buf_cell = [[float(num) for num in line.split()] for line in buf_cell]
       buf_atm = [buf_atm[0].split(),[int(num) for num in buf_atm[1].split()]]
       atm_num = sum(buf_atm[1])
       buf_pos = buf[8:8+atm_num]
       buf_pos = [[float(num) for num in line.split()[0:3]] for line in buf_pos]
       print '\nCell parameters:\n',buf_cell,'\nAtom type and numbers:\n',buf_atm,'\nAtom Positions:\n',buf_pos
       return buf_cell,buf_atm,buf_pos
   else:
       return False

def cell_invr(cell_para):
    V = np.linalg.det(cell_para) 
    if V != 0 :
        cell_para_inv = np.zeros((3,3))
        cell_para_inv[0]= 2*pi/V*np.cross(cell_para[1],cell_para[2])
        cell_para_inv[1]= 2*pi/V*np.cross(cell_para[2],cell_para[0])
        cell_para_inv[2]= 2*pi/V*np.cross(cell_para[0],cell_para[1])
        print 'Inverted cell:\n',cell_para_inv
        return np.matrix(cell_para_inv)
    else:
        print 'Original cell not invertible!'
        return False

def conv_cell_pos_tomat(buf_cell,buf_pos):
    cell_para = buf_cell[0][0]*np.matrix(buf_cell[1:])/0.52918
    'Vasp standard input is angstrom vs eV'
    atom_pos = np.matrix(buf_pos)*cell_para
    cell_para_inv = cell_invr(cell_para)
    return cell_para,cell_para_inv,atom_pos

def mod_chg(elec,buf_atm):
    print 'Calculating Bader charges.'
    "Bader electron counts of the same specie in different environments are taken average."
    chg_val = []
    avg_elec = []
    for atom in buf_atm[0]:
        chg_val.append(atom_dict[atom])
    index = 0
    for atom_num in buf_atm[1]:
        avg_elec.append(sum(elec[index:index+atom_num])/float(atom_num))
        index = index + atom_num
    chg_bad = []
    for i in range(len(chg_val)):
        chg_bad.append([chg_val[i]-avg_elec[i],buf_atm[1][i]])
        print 'Bader charge of atom ',buf_atm[0][i],': ',chg_val[i]-avg_elec[i]
    return chg_bad

def calcbad(cell_para,cell_para_inv,cell_lens,cell_V,atom_pos,chg_bad,G_cutoff,R_cutoff,sqrt_yita):
    "Use Ewald Summation to calculate electrostatic energies between bader charges. n_r is the real space convergence range in lattice vector unit."
    chg_lst = []
    for atom in chg_bad:
        for i in range(atom[1]):
            chg_lst.append(atom[0])
    chg_list = chg_lst
    chg_lst = np.matrix(chg_lst)
    print 'List of Charge:\n',chg_lst,'\nNumber of atoms:',chg_lst.size,'\nList of atom positions:\n',atom_pos
    '----------------------G and T list generation----------------------------'
    a = cell_lens[0]
    b = cell_lens[1]
    c = cell_lens[2]
    a_inv = cell_lens[3]
    b_inv = cell_lens[4]
    c_inv = cell_lens[5]
    print 'R_cutoff:',R_cutoff,',G_cutoff:',G_cutoff,',a:',a,',a_inv:',a_inv
    na = 1 + int(R_cutoff/a)
    nb = 1 + int(R_cutoff/b)
    nc = 1 + int(R_cutoff/c)
    print 'na:',na,',nb:',nb,',nc:',nc
    na_inv = 1 + int(G_cutoff/a_inv)
    nb_inv = 1 + int(G_cutoff/b_inv)
    nc_inv = 1 + int(G_cutoff/c_inv)
    print 'na_inv:',na_inv,',nb_inv:',nb_inv,',nc_inv',nc_inv
    T_lst = [i*cell_para[0]+j*cell_para[1]+k*cell_para[2] for i in range(-na,na+1) for j in range(nb,nb+1) for k in range(nc,nc+1)]
    G_lst = [i*cell_para_inv[0]+j*cell_para_inv[1]+k*cell_para_inv[2] for i in range(-na_inv,na_inv+1) for j in range(-nb_inv,nb_inv+1) for k in range(-nc_inv,nc_inv+1)]
    '--------------------Interaction matrix---------------------------'
    N = chg_lst.size
    D_mat = np.matrix(np.zeros([N,N]))
    for alpha in range(N):
        for beta in range(N):
            tao_a_b = atom_pos[alpha]-atom_pos[beta]
            sum_G = 0
            sum_T = 0
            delta = sqrt_yita/math.sqrt(pi) if alpha == beta else 0
            for G in G_lst:
                norm_G = np.linalg.norm(G)
                if norm_G !=0 :
                    GTao = -1j*G*tao_a_b.T
                    GTao = GTao[0,0]
                    sum_G += cmath.exp(GTao)*math.exp(-1.0*(norm_G**2)/(sqrt_yita**2))/norm_G**2
            sum_G = sum_G.real*4*pi/cell_V
            for T in T_lst:
                Ttao = T+tao_a_b
                norm_Ttao = np.linalg.norm(Ttao)
                if norm_Ttao != 0 :
                    sum_T += math.erfc(norm_Ttao*0.5*sqrt_yita)/norm_Ttao
            D_mat[alpha,beta] = sum_G - delta + sum_T
    print 'Interaction matrix D between atoms:\n',D_mat
    '------------------Total Charge correction-----------------'
    corr = 2*pi*sum(chg_list)**2/(cell_V*(sqrt_yita**2))
    print 'Total Charge correction(Hatree):',corr
    '------------------Final summation---------------------'
    W_mat = 0.5* chg_lst*D_mat*chg_lst.T
    W_N = W_mat[0,0]-corr
    print 'Electrostatic energy per cell(Hatree):',W_N,'\n'
    return W_N
 
def get_cell_lens(cell_para,cell_para_inv):
    a = np.linalg.norm(cell_para[0])
    b = np.linalg.norm(cell_para[1])
    c = np.linalg.norm(cell_para[2])
    a_inv = np.linalg.norm(cell_para_inv[0])
    b_inv = np.linalg.norm(cell_para_inv[1])
    c_inv = np.linalg.norm(cell_para_inv[2])
    return [a,b,c,a_inv,b_inv,c_inv]

N_cut = 50 #at least 2
eps = 0.00001
elec = read_elec()
buf_cell,buf_atm,buf_pos = read_cell_atm_pos()
if elec and buf_cell and buf_atm and buf_pos:
    cell_para,cell_para_inv,atom_pos = conv_cell_pos_tomat(buf_cell,buf_pos)
    chg_bad = mod_chg(elec,buf_atm)
    cell_V = np.linalg.det(cell_para)
    cell_lens = get_cell_lens(cell_para,cell_para_inv)
    a = cell_lens[0]
    b = cell_lens[1]
    c = cell_lens[2]
    a_inv = cell_lens[3]
    b_inv = cell_lens[4]
    c_inv = cell_lens[5]
    G_cutoff = 10 * max([a_inv,b_inv,c_inv])
    R_cutoff = 10 * max([a,b,c])
    buf_res = []
    for n_r in range(1,N_cut):
        sqrt_yita = min([G_cutoff,1/R_cutoff**2]) + abs(G_cutoff-1/R_cutoff**2)/N_cut*n_r
        print 'sqr_yita:',sqrt_yita
        bader_ene = calcbad(cell_para,cell_para_inv,cell_lens,cell_V,atom_pos,chg_bad,G_cutoff,R_cutoff,sqrt_yita)
        buf_res.append([sqrt_yita,bader_ene])
    for i in range(len(buf_res)):
        print buf_res[i][0],' Energy(eV):', 27.2107*buf_res[i][1]

    for i in range(1,len(buf_res)-1):
        if abs(buf_res[i][1]-buf_res[i+1][1])< eps and abs(buf_res[i][1]-buf_res[i-1][1])<0.0001:
            ind_left = i
            break

    for i in range(2,len(buf_res)):
        j = len(buf_res)-i
        if abs(buf_res[j][1]-buf_res[j+1][1])< eps and abs(buf_res[j][1]-buf_res[j-1][1])<0.0001:
            ind_right = j
            break
    ind = int((ind_left + ind_right)/2)
    print '\nEsitimated Bader interaction energy(eV):',buf_res[ind][1]*27.2107,'sqrt-yita',buf_res[ind][0]
    f = open('bader_ene','w')
    f.write(str(buf_res[ind][1]*27.2107))
    f.close()
    plt.plot([pair[0] for pair in buf_res],[pair[1]*27.2107 for pair in buf_res],color='b',label='Electrostatic Energy vs. sqrt-yita')
    plt.plot([pair[0] for pair in buf_res],[buf_res[ind][1]*27.2107 for pair in buf_res],color='r',label='E='+str(buf_res[ind][1]*27.2107)+'eV')
    plt.title('Electrostatic Energy vs. sqrt-yita')
    plt.xlabel('sqrt-yita')
    plt.ylabel('Electrostatic energy/eV')
    plt.legend()
    plt.show()
elif elec == 0:
    print 'Bader analysis file ACF.dat not found! Exiting.'
else:
    print 'Structural information file POSCAR not found! Exiting.'
