#!/usr/bin/python3

import numpy as np

def imgfile(dirstr,filehead,i):
    return dirstr + filehead + '.%i' % (i)

def outfile(dirstr,filehead,i,tailstr):
    return dirstr + filehead + '.%i.%s' % (i,tailstr)

def main():
    datadir = './data/raw/'
    calibfile=datadir + 'tfrecord_chirp-2400_617_interference.calibration'
    print(calibfile)
    calibmat = np.loadtxt(calibfile,dtype=int)
    filehead = 'tfrecord_chirp-2400_617_interference.out'
    imgnum = 10
    print(imgfile(datadir,filehead,imgnum))
    imgmat = np.loadtxt(imgfile(datadir,filehead,imgnum),dtype=int)
    print('imgmat.shape = ',imgmat.shape)
    print(np.matmul(calibmat,imgmat[0,:].T))
    procdir = './data/processed/'

    for i in range(150):
        imgmat = np.loadtxt(imgfile(datadir,filehead,i),dtype=int)
        f = open(imgfile(datadir,filehead,i),'r')
        line = f.readline()
        f.close()
        np.savetxt(outfile(procdir,filehead,i,'innerprod'),np.matmul(calibmat,imgmat.T),fmt='%i',header=line)

if __name__ == '__main__':
    main()
