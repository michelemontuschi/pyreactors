import numpy as np
from . import varconvolve as varcon
import scipy

mp = 938272.08816
mn = 939565.42052
me = 510.99895000

def to_vis(E_nu, nan = False):
#     E_nu = np.linspace(1000,9000,1000)
    E_th = ((mn+me)**2-mp**2)/(2*mp)
    # E_nu_prime = E_nu - E_th
    E_nu_prime = E_nu
    # s = (E_nu+p_p)**2
    # s = s_cin(E_nu)
    s = 2*mp*E_nu_prime+mp**2
    E_nu_CM = (s-mp**2)/(2*np.sqrt(s))
    E_e_CM = (s-mn**2+me**2)/(2*np.sqrt(s))
#     print(s, (s-(mn-me)**2), (s-(mn+me)**2), (s-(mn-me)**2)*(s-(mn+me)**2))
    coeff = (s-(mn-me)**2)*(s-(mn+me)**2)
    if isinstance(E_nu, np.ndarray):
        mask = (coeff>0)
        p_e_CM = np.zeros(len(E_nu))
        p_e_CM[mask] = np.sqrt(coeff[mask])/(2*np.sqrt(s[mask]))
    else:
        if coeff <= 0:
            p_e_CM = 0
        else:
            p_e_CM = np.sqrt(coeff)/(2*np.sqrt(s))
    delta = (mn**2-mp**2-me**2)/(2*mp)
    E1 = E_nu - delta- 1/mp*E_nu_CM*(E_e_CM+p_e_CM)
    E2 = E_nu - delta- 1/mp*E_nu_CM*(E_e_CM-p_e_CM)
    E_e = np.mean(np.array([E1,E2]), axis = 0)
    d_e = (E2-E1)
    sum_E = E_e+me
    if isinstance(E_nu, np.ndarray) and nan:
        mask_nan = (E_nu<E_th)
        d_e[mask_nan] = 0
        sum_E[mask_nan] = np.nan
    elif not isinstance(E_nu, np.ndarray) and nan:
        if E_nu<E_th:
            sum_E = np.nan
            d_e = 0
    return sum_E, d_e

def gaussian(s):
    """
    Constructs a normalized discrete 1D gaussian kernel
    """

    size_grid = int(s*4)
    x= scipy.mgrid[-size_grid:size_grid+1]
    g = np.exp(-(x**2/float(s**2)/2.))
    return g / np.sum(g)

def sigma(E, a = 0.03, b = 0):
    """
    Creates a function that describes the kernel width
    """
    if isinstance(E, np.ndarray):
        E = to_vis(E)[0]
        E = E*1e-3 #from keV to MeV
        E[E<0] = 0.0001
    else:
        E = to_vis(E)[0]
        E = E*1e-3 #from keV to MeV
        if E<0: E=0.0001

    DeltaE_E = np.sqrt((a/np.sqrt(E))**2 + b**2)

    DeltaE = DeltaE_E*E

    k = 1
    std = DeltaE/k
    std = std*1e3 #from MeV to keV

    return std

def top_hat(s):
    size_grid = int(s*4)
    x= scipy.mgrid[-size_grid:size_grid+1]
    g = np.ones(x.shape)
    g[x<-float(s)/2] = 0
    g[x>float(s)/2] = 0
    norm = np.sum(g)
    if norm == 0:
        return g
    else:
        return g / norm

def width(E, nan = True):
    if isinstance(E, np.ndarray):
        w = to_vis(E, nan)[1]
        w[np.isnan(w)]=1e-1
        w[w<=0]=1e-1
    else:
        w = to_vis(E, nan)[1]
#         print(w)
        if w>=0:
            pass
        else:
            w = 1e-1
    return np.absolute(w)

def apply_resolution(E, spectrum, recoil = True):
    E_vis = to_vis(E, True)[0]
    if recoil:
        w = width(E)
        spectrum = varcon.varconvolve(E, spectrum, top_hat, w)
    spectrum = varcon.varconvolve(E, spectrum, gaussian, sigma)
    return E_vis, spectrum
