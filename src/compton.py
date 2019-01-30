#!/home/coffee/anaconda3/bin/python3

import numpy as np
import scipy.constants as consts
from scipy.constants import physical_constants as pc
import sys

(e_mc2,unit,err) = pc["electron mass energy equivalent in MeV"]
e_mc2 *= 1e3 # keV now
print(e_mc2)
(re,unit,err) = pc["classical electron radius"]
re *= 1e2 # in centimeters now

def differentialcrosssection(th,E0 = 511.):
    # from http://xdb.lbl.gov/Section3/Sec_3-1.html (2 of 5) [2/14/2005 6:48:57 PM] Janos Kirz
    k=E0/e_mc2 # in keV
    return np.power(re,int(2))*(1+np.power(np.cos(th),int(2))) / (2.*np.power( ( 1 + k*(1-np.cos(th))) ,int(2)))


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

def differentialcrosssectionofen(en,E0 = 511.):
    k=E0/e_mc2 # in keV
    return np.power(re,int(2))/2. * xofen(en,E0)*yofen(en,E0) / np.power( 1. + k*xofen(en,E0),int(2))

def Precoil(en,E0 = 511.):
    k=E0/e_mc2 # in keV
    return 2.*np.pi * differentialcrosssectionofen(en,E0) / (E0*k) * (1 + (k-1)*xofen(en,E0) + k*np.power(xofen(en,E0),int(2))) / (1.-np.power(k*xofen(en,E0),int(2))) 

def main():
    if len(sys.argv)<2:
        print('syntax is :$ ./compton.py E[in keV]')
        return
    th = np.arange(0,np.pi,np.pi/100)
    E = float(sys.argv[1])
    differentialenergy = differentialcrosssection(th,E)*recoilenergy(th,E)
    filename = "/home/coffee/projects/2dtimetool_simulation/data_fs/reference/compton/differentialenergy_compton.out"
    np.savetxt(filename,np.column_stack((th,differentialcrosssection(th,E),recoilenergy(th,E),differentialenergy)),fmt='%.4e')

    en = np.arange(0,511,5)
    filename = "/home/coffee/projects/2dtimetool_simulation/data_fs/reference/compton/differentialenergy_compton_recoil.out"
    np.savetxt(filename,np.column_stack((en,differentialcrosssectionofen(en,E),Precoil(en,E))),fmt='%.4e')
    return


if __name__ == '__main__':
    main()
