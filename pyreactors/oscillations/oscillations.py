import numpy as np

conv_factor = 1.26693 # See Lisi_JUNO_1st_Lecture_2021.pdf for (Dm_ij in ev2, x in m, E in MeV)

def getParams(hierarchy = 'NH'):
   
    params = dict()
    if hierarchy == 'NH':

        params['dm2_21'] = 7.54e-5
        params['Dm2_32'] = 0.0023969136#2.453e-3
        params['s2_12'] = 3.07e-1
        params['s2_13'] = 2.18e-2
        params['s2_23'] = 5.45e-1

        #dm2_21 = 7.34e-5
        #s2_12 = 0.304
        #s2_13 = 2.16e-2




        params['c2_12'] = 1 - params['s2_12']
        params['c2_13'] = 1 - params['s2_13']
        params['c2_23'] = 1 - params['s2_23']

    elif hierarchy == 'IH':
        
        params['dm2_21'] = 7.53e-5
        params['Dm2_32'] = 2.546e-3
        params['s2_12'] = 3.07e-1
        params['s2_13'] = 2.18e-2
        params['s2_23'] = 5.47e-1


        
        params['c2_12'] = 1 - params['s2_12']
        params['c2_13'] = 1 - params['s2_13']
        params['c2_23'] = 1 - params['s2_23']

    else: raise Exception('Choose between NH or IH')

    return params




######### Electron neutrino oscillation survival probability [F. Capozzi, E. Lisi, and A. Marrone Phys. Rev. D 89, 013001 – Published 9 January 2014] Eq. 39 ##########


def Pee3v_39(d, E_nu, hierarchy = 'NH', inMatter = False): # m, keV
    # Return electron neutrino oscillation survival probability in vacuum [F. Capozzi, E. Lisi, and A. Marrone Phys. Rev. D 89, 013001 – Published 9 January 2014] Eq. 39
    #
    # d = distance
    # E_nu =  Neutrino Energy axis
    # hierarchy = 'NH' if Normal Hierarchy, 'IH' if inverse hierarchy
    
    params = getParams(hierarchy)
    if hierarchy == 'NH':
        alpha = 1.
    elif hierarchy == 'IH':
        alpha = -1.

    eMev = E_nu # Neutrino Energuy in MeV
    
    dm2_21 = params['dm2_21']

    s2_12 = params['s2_12']
    c2_12 = params['c2_12']
    s2_13 = params['s2_13']
    c2_13 = params['c2_13']

    c4_13 = c2_13**2
    s4_13 = s2_13**2


    conv_factor = 1.26693
    d = 52500

    Delta_small = get_delta_small(dm2_21 , d , E_nu) # dm = dm2_21

    s2_2th12 = 4 * s2_12 * c2_12

    w = 1
    P2v = P2_vac(s2_2th12, Delta_small)

    if inMatter:
        #P2v = P2_mat(s2_2th12, eMev, dm2_21, conv_factor, d, Ne=1.3)
        w = get_damp_fact(d, E_nu, [], [], hierarchy)

    
    D_ee = get_Dee(d,eMev , hierarchy)
    phi, phi_2 = get_phi(d, eMev, hierarchy)

    prob = c4_13 * P2v + s4_13 + 2 * s2_13 * c2_13 * np.sqrt(P2v) * w * np.cos(2 * D_ee + alpha * phi)

    
    return prob



def Pee_39(d, E_nu, hierarchy = 'NH', inMatter=False): 
    # Return electron neutrino oscillation survival probability
    # d = distance
    # E_nu =  Neutrino Energy axis
    # hierarchy = 'NH' if Normal Hierarchy, 'IH' if inverse hierarchy
    
    return Pee3v_39(d, E_nu, hierarchy, inMatter)


def get_delta_small(dm2_21 , d , E_nu):
    Delta_small = conv_factor * dm2_21 * d / E_nu  # dm = dm2_21
    return Delta_small

def get_Deltam2(Dm2_32, dm2_21):
    Deltam2 = Dm2_32 + 0.5 * dm2_21
    return Deltam2
def get_Dm2_ee ( Deltam2 ,alpha,c2_12,s2_12,dm2_21):
    Dm2_ee = Deltam2 + (1 / 2) * alpha * (c2_12 - s2_12) * dm2_21
    return Dm2_ee

def get_phi(d, E_nu, hierarchy = 'NH'): # m, keV

    params = getParams(hierarchy)


    dm2_21 = params['dm2_21']
    s2_12 = params['s2_12']
    c2_12 = params['c2_12']

    Delta_small = get_delta_small(dm2_21 , d , E_nu)  # dm = dm2_21

    s2_2th12 = 4 * s2_12 * c2_12
    P2v_vac = P2_vac(s2_2th12, Delta_small)
    eta = np.sqrt(P2v_vac)

    phi = 2 * s2_12 * Delta_small * (1 - (np.sin(2 * Delta_small) / (2 * Delta_small * eta)))
    phi_2 = 2 * s2_12 * Delta_small * (1 - (np.sin(Delta_small) / (2*Delta_small * eta)))

    return phi, phi_2

def get_Dee(d, E_nu, hierarchy = 'NH'):

    params = getParams(hierarchy)
    if hierarchy == 'NH':
        alpha = 1.
    elif hierarchy == 'IH':
        alpha = -1.

    eMev = E_nu  # Neutrino Energy in MeV

    dm2_21 = params['dm2_21']
    Dm2_32 = params['Dm2_32']
    s2_12 = params['s2_12']
    c2_12 = params['c2_12']


    Deltam2 = get_Deltam2(Dm2_32, dm2_21)
    Dm2_ee= get_Dm2_ee(Deltam2, alpha, c2_12, s2_12, dm2_21)
    D_ee = conv_factor * Dm2_ee * d / (eMev)


    return D_ee

def P2_vac(s2_2th12, Delta_small):
    # Return electron 2-neutrino oscillation survival probability in vacuum [F. Capozzi, E. Lisi, and A. Marrone Phys. Rev. D 89, 013001 – Published 9 January 2014] Eq. 40
    #
    # s2_2th12 = 4*s2_12*c2_12

    P2v_vac = 1 - s2_2th12 * (np.sin(Delta_small)) ** 2
    return P2v_vac

def P2_mat(s2_2th12, eMev, dm2_21, conv_factor, d, Ne=1.3):

    # Return electron 2-neutrino oscillation survival probability in matter [F. Capozzi, E. Lisi, and A. Marrone Phys. Rev. D 89, 013001 – Published 9 January 2014] Eq. 56-57-59

    r = 1.526e-7 * Ne * eMev / dm2_21

    s_2th12 = np.sqrt(s2_2th12)
    c_2th12 = np.sqrt(1 - s2_2th12)

    #print(r * c_2th12)
    r_cos2th12=get_rcos_mat(s2_2th12, eMev, dm2_21, Ne)
    s_2th12_ma = (s_2th12 * (1 - r_cos2th12))
    s2_2th12_ma = s_2th12_ma**2

    dm2_21_ma = dm2_21 * (1 + r_cos2th12)
    Delta_small_ma = get_delta_small(dm2_21_ma, d, eMev)  # dm = dm2_21

    P2_mat = P2_vac(s2_2th12_ma, Delta_small_ma)
    #P2_mat = 1 - (s2_2th12_ma) * np.sin(Delta_small_ma) ** 2

    return P2_mat

def get_rcos_mat(s2_2th12, eMev, dm2_21, Ne=1.3):

    r = 1.526e-7 * Ne * eMev / dm2_21
    c_2th12 = np.sqrt(1 - s2_2th12)

    return r * c_2th12

def get_damp_fact(d, E_nu, Ln, wn, hierarchy = 'NH'):
    params = getParams(hierarchy)
    if hierarchy == 'NH':
        alpha = 1.
    elif hierarchy == 'IH':
        alpha = -1.

    eMev = E_nu  # Neutrino Energuy in MeV

    dm2_21 = params['dm2_21']
    Dm2_32 = params['Dm2_32']
    s2_12 = params['s2_12']
    c2_12 = params['c2_12']


    Deltam = Dm2_32 + 0.5 * dm2_21 #Dm = 1/2 * (Dm2_31 + Dm2_32)
    Dm_ee = Deltam + (1/2) * alpha * (c2_12 - s2_12) * dm2_21

    D_ee = conv_factor * Dm_ee * d / eMev

    #damp_fact = 0

    sum_wn_ln=2.16e-5

    for l_i, w_i in zip(Ln, wn):
        sum_wn_ln = sum_wn_ln + (w_i * (l_i / d)**2)

    damp_fact = (1 - 4 * (D_ee ** 2) * sum_wn_ln)

    return damp_fact


def ratioPhase(d, E_nu, hierarchy = 'NH'): # m, keV
    # Return electron neutrino oscillation survival probability in vacuum [F. Capozzi, E. Lisi, and A. Marrone Phys. Rev. D 89, 013001 – Published 9 January 2014] Eq. 39
    #
    # d = distance
    # E_nu =  Neutrino Energy axis
    # hierarchy = 'NH' if Normal Hierarchy, 'IH' if inverse hierarchy
    

    D_ee = get_Dee(d, E_nu, hierarchy = 'NH')
    phi, phi_2 = get_phi(d, E_nu, hierarchy = 'NH')
    ratio = phi / (2*D_ee)

    
    return phi,  D_ee, ratio, phi_2


