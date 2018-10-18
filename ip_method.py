#!/usr/bin/python3

import numpy as np
import re as regexp
import subprocess
from scipy import sparse

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
	imgnum = 0
	#print(imgfile(datadir,filehead,imgnum))
	imgmat = np.loadtxt(imgfile(datadir,filehead,imgnum),dtype=int)
	(nfibers,npixels) = imgmat.shape
	#print('nfibers = ', nfibers)
	#print(np.matmul(calibmat,imgmat[0,:].T))
	procdir = './data/processed/'

	nimages=100
	OUT = np.zeros((nimages,2*imgmat.shape[0]+2),dtype=int)

	for i in range(nimages):
			imgmat = np.loadtxt(imgfile(datadir,filehead,i),dtype=int)
			result = np.matmul(calibmat,imgmat.T)
			f = open(imgfile(datadir,filehead,i),'r')
			line = f.readline().rstrip()
			f.close()
			string,delay = line.split('\t')
			d = float(delay)
			np.savetxt(outfile(procdir,filehead,i,'innerprod'),result,fmt='%i',header=line)
			inds = np.argmax(result,axis=0)
			maxs = np.max(result,axis=0)
			OUT[i,0] = i
			OUT[i,1] = d
                        OUT[i,2::2] = inds
                        OUT[i,3::2] = maxs
	outfname = procdir+filehead + '.innerprod_results'
	np.savetxt(outfname,OUT,fmt='%i')
        """
        OK, mask on maxs.  If where maxs > 2e8 then fit 
        OUT[:,2]... ugh too tired
        """


if __name__ == '__main__':
    main()
