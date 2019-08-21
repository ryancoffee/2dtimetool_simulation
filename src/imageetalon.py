#!/usr/bin/python3.5

import numpy as np
import sklearn as sk
import cv2
import sys
import re

def gauss(f,center,bwd):
    return np.exp(-np.power((f-center)/bwd,int(2)))


def main():
    if len(sys.argv)<1:
        print('Hey, syntax is: main imagefile fibermapfile')
        return
    m = re.search('(.+)raw/(.+).out.(\d+)',sys.argv[1])
    print('file:\t{}.out.{}'.format(m.group(2),m.group(3)))
    if m:
        basedir = m.group(1)
        filebase = m.group(2)
        imagnum = m.group(3)
        data = np.loadtxt(sys.argv[1],dtype=np.uint16)
        imgdata = np.zeros((data.shape[0]*10,data.shape[1]),dtype=float)
        imgdata[::10,:]=data
        #fibermap = np.loadtxt(sys.argv[2],dtype=float)
        IMGDATA = np.fft.fft2(imgdata)
        f0 = np.fft.fftfreq(imgdata.shape[0])
        f1 = np.fft.fftfreq(imgdata.shape[1])
        m0,m1 = np.meshgrid(gauss(f1,0,.2),gauss(f0,0,.075))
        s0,s1 = np.meshgrid(np.ones(f1.shape[0])*(np.abs(f1)>.0125),np.ones(f0.shape[0])) #*(np.abs(f0)>.00625)) # don't do it in the vertical dimension... too confusing
        result = np.fft.ifft2(IMGDATA * m0 * m1).real
        outfilename = '{}processed/{}.out.{}'.format(basedir,filebase,imagnum)
        np.savetxt(outfilename,result,fmt='%i')
        outfilename = '{}processed/{}.{}.jpg'.format(basedir,filebase,imagnum)
        cv2.imwrite(outfilename,result//2)
        result = np.fft.ifft2(IMGDATA * m0 * m1 * s0 * s1).real
        result -= np.min(result)
        result /= np.max(result)
        result *= (2**8) 
        smask = np.fft.ifft(np.sum(s0 * s1,axis=0))
        outfilename = '{}processed/{}.{}.schlieren.real.jpg'.format(basedir,filebase,imagnum)
        cv2.imwrite(outfilename,result.astype(np.uint))
        outfilename = '{}processed/{}.{}.schlieren.mask.dat'.format(basedir,filebase,imagnum)
        np.savetxt(outfilename,np.column_stack((f1,smask.real,smask.imag)),fmt='%.4e')
        ETALON = np.fft.fft(result,axis=1)
        m0,m1 = np.meshgrid(gauss(np.abs(f1),0.750,.25),np.ones(f0.shape))
        ETALON *= m0
        etalon = np.fft.ifft(ETALON,axis=1).real
        etalon -= np.mean(etalon[100:-100,:])
        etalon /= np.max(etalon[100:-100,:])
        etalon *= 2**7
        etalon += 2**7
        outfilename = '{}processed/{}.{}.schlieren.etalon.jpg'.format(basedir,filebase,imagnum)
        cv2.imwrite(outfilename,etalon.astype(np.uint))
        ETPOWER = np.power(np.abs(ETALON[:,:ETALON.shape[1]//2]),int(2))
        outfilename = '{}processed/{}.{}.schlieren.ETALON.dat'.format(basedir,filebase,imagnum)
        np.savetxt(outfilename,ETPOWER,fmt='%.4e')
        numerator = np.matmul(ETPOWER[:,100:],np.arange(100,ETPOWER.shape[1]).reshape(ETPOWER.shape[1]-100,1))
        denominator = np.sum( ETPOWER[:,100:] , axis=1).reshape(numerator.shape)
        inds = np.where(denominator>0)
        etalonfreq = np.zeros(numerator.shape)
        etalonfreq[inds] = numerator[inds]/denominator[inds]
        #etalonfreq = np.argmax(ETPOWER[:,100:],axis=1)
        outfilename = '{}processed/{}.{}.schlieren.etalonfreq.dat'.format(basedir,filebase,imagnum)
        np.savetxt(outfilename,etalonfreq,fmt='%.3f')

    return


if __name__ == '__main__':
    main()
