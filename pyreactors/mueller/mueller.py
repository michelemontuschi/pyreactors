#[T. Mueller et al. Phys. Rev. C 83, 054615 (2011).]
import numpy as np

a_U235 = dict(a1 = 3.217, a2 = -3.111, a3 = 1.395, a4 = -3.690e-1, a5 = 4.445e-2, a6 = -2.053e-3)
a_Pu239 = dict(a1 = 6.413, a2 = -7.432, a3 = 3.535, a4 = -8.820e-1, a5 = 1.025e-1, a6 = -4.550e-3)
a_U238 = dict(a1 = 4.833e-1, a2 = 1.927e-1, a3 = -1.283e-1, a4 = -6.762e-3, a5 = 2.233e-3, a6 = -1.536e-4)
a_Pu241 = dict(a1 = 3.251, a2 = -3.204, a3 = 1.428, a4 = -3.675e-1, a5 = 4.254e-2, a6 = -1.896e-3)

elements = dict(U235 = a_U235, Pu239 = a_Pu239, U238 = a_U238, Pu241 = a_Pu241)

def mueller(E, params): #Get Mueller Spectrum
    a_list = np.asarray(list(params.values()))
    y = np.zeros(len(E))
    for i, a in enumerate(a_list):
        temp = a * np.power(E, i)
        y = y + temp
   
    return np.exp(y)

def spectrum(E, element):
    return mueller(E, elements[element])
