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

def fourier_shift(x,t=0.):
    X=fft(x)
    F=fftfreq(x.shape[0])
    X *=np.exp(-1j*2.*np.pi*F*t)
    return np.real(ifft(X))

def fourier_shift_integrate(x,t=0.):
    X=fft(x)
    F=fftfreq(x.shape[0])
    X *=np.exp(-1j*2.*np.pi*F*t)
    denom = 1j*F
    denom[1:] = np.power(denom[1:],int(-1))
    denom[0]=0
    result = np.real(ifft(X*denom))
    result -= np.min(result)
    result /= np.max(result)
    return result

def main():
    xvals = np.arange(-10,1000000)
    np.savetxt('testcascade.out',diamond_cascade(xvals),fmt='%.4f')
    yvals = fourier_shift(diamond_cascade(xvals),10)
    np.savetxt('testcascade.shifted.out',yvals,fmt='%.4f')
    yvals = fourier_shift_integrate(diamond_cascade(xvals),0)
    np.savetxt('testcascade.integrated.out',yvals,fmt='%.4f')


    return 0

if __name__ == '__main__':
    main()
