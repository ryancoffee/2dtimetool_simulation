#!/home/coffee/anaconda3/bin/python3

import numpy as np
from scipy.special import erf,erfc
from numpy.fft import fft,ifft,fftfreq

def diamond_cascade(x):
    wrise           = 1.19298 #         +/- 0.0371       (3.11%)
    y0              = 0.697986 #        +/- 0.004242     (0.6078%)
    slope           = 0.0147671 #       +/- 0.0003711    (2.513%)
    xfall           = 31.2418    #      +/- 0.06072      (0.1943%)
    wfall           = 6.72156     #     +/- 0.09024      (1.342%)
    result = 0.5*(1.+erf(x/wrise))*(y0+slope*x)*.5*(erfc((x-xfall)/wfall))
    return result

def integrate(x):
    _y = np.cumsum(x)
    return _y

def fourier_shift(x,t=0.,alpha=1e-4):
    _X=fft(x)
    _F=fftfreq(x.shape[0])
    _alpha = alpha if len(x)>int(1./alpha) else 1./(len(x))
    _denom = (_alpha+1j*_F)
    _denom = np.power(_denom,int(-1))
    _denom *=np.exp(-1j*2.*np.pi*_F*(x[0]))
    _X *=np.exp(-1j*2.*np.pi*_F*t)
    _X *= _denom
    return np.real(ifft(_X))

def fourier_shift_integrate(x,t=0.):
    X=fft(x)
    F=fftfreq(x.shape[0])
    X *=np.exp(-1j*2.*np.pi*F*t) # time shift
    denom = 1j*F # integrate
    denom[1:] = np.power(denom[1:],int(-1))
    denom[0]=0
    result = np.real(ifft(X*denom))
    result -= np.min(result)
    result /= np.max(result)
    return result

def fourier_shift_decay_integrate(x,t=0.,alpha=10e-3):
    X=fft(x)
    _alpha = alpha if len(x)>int(1./alpha) else 1./(len(x))
    F=fftfreq(x.shape[0])
    X *=np.exp(-1j*2.*np.pi*F*t) #time shift
    denom = 1j*F*(_alpha+1j*F)
    denom[1:] = np.power(denom[1:],int(-1))
    denom[0]=0
    result = np.real(ifft(X*denom))
    result -= np.min(result)
    result /= np.max(result)
    return result

def main():
    span = 1e4
    xvals = np.arange(-span//100,99*span//100,dtype=float)
    yvals = diamond_cascade(xvals)
    dvals = integrate(yvals)
    dvals *= np.sum(yvals)/np.sum(dvals)
    print(np.sum(yvals-dvals))
    out = np.column_stack((xvals,yvals,dvals))
    newyvals = yvals-dvals
    np.savetxt('./data_container/raw/testcascade.out',out,fmt='%.4f')
    out = np.column_stack((xvals,yvals,dvals,fourier_shift_integrate(yvals),fourier_shift_integrate(newyvals)))
    np.savetxt('./data_container/raw/testcascade.integrated.out',out,fmt='%.4f')


    return 0

if __name__ == '__main__':
    main()
