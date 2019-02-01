#!/usr/bin/python3

import numpy as np
import scipy.constants as consts
from scipy.constants import physical_constants as pc
from scipy.constants import pi
import sys

e_mc2 = pc["electron mass energy equivalent in MeV"][0] * 1e3 # keV now
print('electron rest mass energy %.4f keV' % e_mc2)
re = pc["classical electron radius"][0] * 1e2 # in centimeters
re2 = np.power(re,int(2))*1e24 # in barns now
print('classical electron area %.4f barns' % (pi*re2))

def PofE(en,E0 = 511.):
    k = E0/e_mc2 # in keV
    ka = np.power(E0,int(3))/np.power(e_mc2,int(2))
    b = en/(ka-k*en)
    result = (pi*re2)/ka * (2*b - np.power( b , int(2))) 
    result[np.where(result<0.)] = 0
    return result

def differentialcrosssection(th,E0 = 511.):
    # from http://xdb.lbl.gov/Section3/Sec_3-1.html (2 of 5) [2/14/2005 6:48:57 PM] Janos Kirz
    k=E0/e_mc2 # in keV
    return re2*(1+np.power(np.cos(th),int(2))) / (2.*np.power( ( 1 + k*(1-np.cos(th))) ,int(2)))


def recoilenergy(th,E0 = 511.):
    k=E0/e_mc2 # in keV
    return E0*(1-np.cos(th)) / (1 + k*(1-np.cos(th)) )

def x(th):
    return (1.-np.cos(th))

def xofen(en,E0 = 511.):
    k=E0/e_mc2 # in keV
    return en/(k* (E0 - en))

def yofen(en,E0 = 511.):
    k=E0/e_mc2 # in keV
    return 2. - en/(k* (E0 - en))

def main():
    if len(sys.argv)<2:
        print('syntax is :$ ./compton.py <dE=10[keV]> <E[in keV]=511.>')
    E = e_mc2
    dE = 25
    if len(sys.argv)>1:
        dE = float(sys.argv[1])
        if len(sys.argv)>2:
            E = float(sys.argv[2])

    th = np.arange(0,np.pi,np.pi/100)
    filename = "./data_fs/reference/compton/differentialenergy_compton.%ikeV.out" % int(E)
    headerstring = 'theta\tdifferential xsection [barns]\trecoil energy [keV]\tscattered photon energy [keV]'
    np.savetxt(filename,np.column_stack((th,differentialcrosssection(th,E),recoilenergy(th,E),E-recoilenergy(th,E))),fmt='%.4e',header=headerstring)

    en = np.arange(0,e_mc2,dE)
    filename = "./data_fs/reference/compton/differentialenergy_compton_recoil.%ikeV.out" % int(E)
    headerstring = 'recoil energies\tprobability [barns/keV]\tcorresponding compton photon energies [keV]'
    np.savetxt(filename,np.column_stack((en,PofE(en,E),E-en)),fmt='%.4e',header=headerstring)

    Eb = 50
    p1 = PofE(en,E)
    p1sum = np.sum(p1)
    inds = np.where(p1>0)
    nu1 = E-en
    out = np.copy(p1)
    '''
    secondary_photos = nu1-Eb
    p_secondary_photos = np.zeros(p1.shape,dtype=float)
    p_secondary_photos[:-Eb//dE] = np.copy(p1[Eb//dE:])
    print(secondary_photos)
    print(p1)
    print(p_secondary_photos) # need to get the secondary photos in ratio to the secondary compton events
    '''
    for i in range(len(p1)):
        if p1[i]>0:
            out += PofE(en,E0=nu1[i])*p1[i]/p1sum # need to multiply by the ratio of Compton to Photo... which is a function of en
            # note to self... the files to look at are in laptop$ $HOME/projects/2dtimetool_simulation/docs/Compton/
    filename = "./data_fs/reference/compton/compton_recoil_electrons.2ndOrder.%ikeV.out" % int(E)
    headerstring = 'recoil energies\tprobability [barns/keV] primary Compton\tprobability [barns/keV] including secondary Compton events still ignoring photoelectrons'
    np.savetxt(filename,np.column_stack((en,p1,out)),fmt='%.4e',header=headerstring)
    return


if __name__ == '__main__':
    main()