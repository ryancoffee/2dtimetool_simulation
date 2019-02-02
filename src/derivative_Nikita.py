#!/usr/bin/python3

import numpy as np
import sys

def main():
    for filename in sys.argv[1:]:
        data = np.loadtxt(filename,dtype = float)
        times=data[:,0]
        nv_holes=data[:,5]
        d_nv_holes=np.diff(data[:,5])/(times[1]-times[0])
        out = np.column_stack((times[1:],d_nv_holes))
        np.savetxt(filename + '.derivative',out,fmt='%.4f')
    return

if __name__ == '__main__':
    main()
