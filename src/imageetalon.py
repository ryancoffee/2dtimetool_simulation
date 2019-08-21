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
        s0,s1 = np.meshgrid(np.ones(f1.shape[0])*(np.abs(f1)>.0125),np.ones(f0.shape[0]))
        result = np.fft.ifft2(IMGDATA * m0 * m1).real
        outfilename = '{}processed/{}.out.{}'.format(basedir,filebase,imagnum)
        np.savetxt(outfilename,result,fmt='%i')
        outfilename = '{}processed/{}.{}.jpg'.format(basedir,filebase,imagnum)
        cv2.imwrite(outfilename,result//2)
        result = np.fft.ifft2(IMGDATA * m0 * m1 * s0 * s1)
        smask = np.fft.ifft2(s0 * s1)
        outfilename = '{}processed/{}.{}.schlieren.sqrabs.jpg'.format(basedir,filebase,imagnum)
        cv2.imwrite(outfilename,np.power(np.abs(result),int(2)).astype(np.uint)//8)
        outfilename = '{}processed/{}.{}.schlieren.mask.real.dat'.format(basedir,filebase,imagnum)
        np.savetxt(outfilename,smask.real.T*1e20)
        outfilename = '{}processed/{}.{}.schlieren.mask.imag.dat'.format(basedir,filebase,imagnum)
        np.savetxt(outfilename,smask.imag.T*1e20)
    return


if __name__ == '__main__':
    main()
