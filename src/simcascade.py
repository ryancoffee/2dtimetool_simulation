#!/home/coffee/anaconda3/bin/python3

import numpy as np
from scipy.special import erf,erfc
from numpy.fft import fft,ifft,fftfreq

def diamond_cascade_params(x,energy=9.5):
    # ultimately, these were had fit with Nikita's 2015 derivative curves... those are too slow by 50% or so he says
    y=np.zeros(x.shape)
    decay=np.zeros(x.shape)
    wfall = 0.77863 * energy
    xfall = 0.0916548 * np.power(float(energy),int(2)) + 2.56726 * energy
    quad =  0.014 * np.power(float(energy),int(-2))
    slope = .95 * np.power(float(energy),int(-2))
    y0 = .725 
    inds = np.where(x>0)
    y[inds] = (y0+slope*x[inds]+quad*np.power(x[inds],int(2)))*.5*(erfc((x[inds]-xfall)/wfall))
    decay[inds] = 0.5*erf((x[inds]-xfall)/wfall)
    D = np.sum(decay[inds])
    S = np.sum(y[inds])
    y[inds] -= decay[inds] * S/D
    return y
    
    

def diamond_cascade(x,energy=9):
    # default 9 keV case
    # wrise is really just the erf of the gaussian x-ray pulse duration
    wrise           = 1.19298 #         +/- 0.0371       (3.11%)
    y0              = 0.697986 #        +/- 0.004242     (0.6078%)
    slope           = 0.0147671 #       +/- 0.0003711    (2.513%)
    xfall           = 31.2418    #      +/- 0.06072      (0.1943%)
    wfall           = 6.72156     #     +/- 0.09024      (1.342%)
    quad=0;
    if (int(energy) == 24):
        #Final set of parameters            Asymptotic Standard Error
        #=======================            ==========================
        quad            = 2.26439e-05  #    +/- 1.894e-06    (8.364%)
        y0              = 0.673358     #    +/- 0.003428     (0.509%)
        slope           = 0.0016597    #    +/- 0.0001759    (10.6%)
        wrise           = 6.25389      #    +/- 0.1081       (1.729%)
        wfall           = 18.7944      #    +/- 0.1384       (0.7362%)
        xfall           = 114.331      #    +/- 0.1409       (0.1232%)
    if (int(energy) == 3):
        # this case is all hand fitted to g3data extracted and diffed, also to combination of SiO2@9keV for fast rising edge
        # 2015_Medvedev_ApplPhysB_XCASCADE.3keV.diamond.dat
# even better, using the fast rise of SiO2 derivative at 9 keV matches well the 3keV diamond derivative rise... more or less.  At least to get the slope term
        quad            = 0.
        y0              = 1.1
        slope           = 0.12
        wrise           = 1.
        wfall           = 12.
        xfall           = 7.
    result = 0.5*(1.+erf(x/wrise))*(y0+slope*x+quad*np.power(x,int(2)))*.5*(erfc((x-xfall)/wfall))
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
    y = diamond_cascade_params(xvals,energy=3)
    y2 = diamond_cascade_params(xvals,energy=5)
    y3 = diamond_cascade_params(xvals,energy=9)
    y4 = diamond_cascade_params(xvals,energy=15)
    y5 = diamond_cascade_params(xvals,energy=24)
    out = np.column_stack((xvals,y,y2,y3,y4,y5))
    np.savetxt('cascades_function.out',out,fmt='%.4f')
    out = np.column_stack((xvals,y5,fourier_shift_integrate(y5)))
    np.savetxt('./data_container/raw/testcascade.integrated.out',out,fmt='%.4f')

    return 0

if __name__ == '__main__':
    main()
