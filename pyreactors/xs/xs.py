#Cross section of IBD calculation [Alessandro Strumia and Francesco Vissani, Phys. Lett. B 564 42 (2003); astro-ph/0302055.]
import numpy as np

def xs_IBD(E): # MeV
    E_nu = E
    if isinstance(E_nu, np.ndarray):
        pass
    else:
        E_nu = np.asarray([E_nu])
        
    Q = 1.2933 #MeV
    E_e = E_nu - Q #total energy of the positron
    dm_e = 0.511 #MeV
    E_e[E_e**2<dm_e**2] = dm_e
    p_e = np.sqrt(E_e**2 - dm_e**2) #Momentum of the positron
    dlog_E_nu= np.log(E_nu)
    exponent = -0.07056 + 0.02018*dlog_E_nu - 0.001953*dlog_E_nu**3
    sigmasection=1e-43*E_e*p_e*(E_nu**exponent)
    sigmasection[E_nu<1.806] = 0.0
    return sigmasection
