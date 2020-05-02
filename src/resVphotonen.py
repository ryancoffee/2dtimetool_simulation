#!/usr/bin/python3

import numpy as np
import sys
from scipy.optimize import leastsq
import random

def gauss(params,x):
    a = params[0]
    c = params[1]
    w = params[2]
    return a*np.exp(-1.*np.power((x-c)/w,int(2)))

def gaussresidual(params,x,data,eps_data):
    model = gauss(params,x)
    return (model - data)/eps_data

def getcenterofmass(data):
    sz = data.shape[1]
    x = np.arange(sz).astype(float).reshape((1,sz))
    num = np.inner(data,x)
    denom = np.sum(data,axis=1)
    return (num.reshape(denom.shape)/denom).astype(int)

def processfile(fname):
    data = np.loadtxt(fname)
    derivdata = np.loadtxt(fname + '.derivative')
    cominds = getcenterofmass(data)
    outmaxfits = []
    outminfits = []
    highlim = derivdata.shape[1]
    for i in range(derivdata.shape[0]):
        if cominds[i]<2:
            cominds[i] = data.shape[1]//2
        maxderivind = np.argmax(derivdata[i,:cominds[i]])
        cenind = maxderivind
        lowind = cenind-2
        highind = cenind+2
        maxval = derivdata[i,cenind]
        while lowind > 0 and derivdata[i,lowind] > maxval/4:
            lowind -= 1
        while highind < highlim-2 and derivdata[i,highind] > maxval/4:
            highind += 1
        params = (derivdata[i,maxderivind] , float(maxderivind) , (highind-lowind)/4.)
        xvals = np.arange(lowind,highind-1)
        eps_data = [1 for i in range(len(xvals))]
        yvals = derivdata[i,xvals]
        out,niters = leastsq(gaussresidual,params,args=(xvals,yvals,eps_data))
        if len(outmaxfits)<1:
            outmaxfits = out
        else:
            outmaxfits = np.row_stack((outmaxfits,out))
    for i in range(derivdata.shape[0]):
        minderivind = cominds[i] + np.argmin(derivdata[i,cominds[i]:])
        cenind = minderivind
        lowind = cenind-2
        highind = cenind+2
        minval = derivdata[i,cenind]
        while lowind > 0 and derivdata[i,lowind] < maxval/4:
            lowind -= 1
        while highind < highlim-2 and derivdata[i,highind] < maxval/4:
            highind += 1
        params = (derivdata[i,minderivind] , float(minderivind) , (highind-lowind)/4.)
        xvals = np.arange(lowind,highind-1)
        eps_data = [1 for i in range(len(xvals))]
        yvals = derivdata[i,xvals]
        out,niters = leastsq(gaussresidual,params,args=(xvals,yvals,eps_data))
        if len(outminfits)<1:
            outminfits = out
        else:
            outminfits = np.row_stack((outminfits,out))
    ofname = fname + '.maxminfits'
    headstring = 'max::amp\tmax::center\tmax::width\tmin::amp\tmin::center\tmin::width'
    np.savetxt(ofname,np.column_stack((outmaxfits,outminfits)),fmt='%.1f',header=headstring)
    print('Saved: %s'%(ofname))
    return

def main():
    if len(sys.argv)<2:
        print('SYntax: ./resVphotonen.py <filetoprocess -- list of as well>')
        return
    for fname in sys.argv[1:]:
        processfile(fname)
    return

if __name__ == '__main__':
    main()
