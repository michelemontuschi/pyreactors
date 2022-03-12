import numpy as np
from . import mueller, xs, fuels, oscillations
from . import distance as proj
import pandas as pd

days = dict(
Jan = 31,
Feb = 28, # in a common year and 29 days in leap years
Mar = 31,
Apr = 30,
May = 31,
Jun = 30,
Jul = 31,
Aug = 31,
Sep = 30,
Oct = 31,
Nov = 30,
Dec = 31 )

def readDB(filename, sheet = False): #Get Info from Reactor Database
    if sheet:
        data = pd.read_excel(filename, sheet)

    else:
        data = pd.read_excel(filename)

    fuels = np.asarray(['STD'] * data.shape[0])
    fuels[data['isMox'] == 1] = 'MOX'
    fuels[data['isCandu'] == 1] = 'CAN'
    fuels[(data['isCandu'] == 1) & (data['isMox'] ==1) ] = 'BU'
    data['Fuel'] = fuels #Add a column to data where is written if STD or MOX or CAN
    mean_load = 0
    tot_days = 0
    
    for month in days:
        mean_load = mean_load+data[month]*days[month]
        tot_days = tot_days+days[month]
        data['LF'] = mean_load/tot_days #Calculation of the load factor

    return data


def plotDB(data):
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.io.img_tiles as cimgt
    
    def fuel2col(fuel):
        if fuel == 'STD':
            return 'b'
        elif fuel == 'MOX':
            return 'r'
        elif fuel == 'CAN':
            return 'g'
        else:
            return 'k'

    my_proj = ccrs.PlateCarree()

    scale = 1.4
    fig = plt.figure(figsize=(scale*6.1637, scale*3))
    ax = plt.axes(projection=my_proj)
    ax.stock_img()
    ax.scatter(data['Lon'].to_numpy(), data['Lat'].to_numpy(), transform=ccrs.PlateCarree(), ec = 'k', s = 20, c = data['Fuel'].apply(fuel2col))

    fig.tight_layout()

    return

def mueller_spectra(E): #Convolve Muller spectrum with IBD cross section
    spe = dict()
    for el in mueller.elements:
        spe[el] = mueller.spectrum(E, el) #[#nu/(MeV*fission)]
    return spe

def cross_section(E): #Convolve Muller spectrum with IBD cross section
    c_section = xs.xs_IBD(E) # [cm^2]
    return c_section


def build_spectra(E): #Convolve Muller spectrum with IBD cross section
    spe = dict()
    for el in mueller.elements:

        c_section = xs.xs_IBD(E) # [cm^2]
        c_section = c_section * 1e-4  # [m^2]
        m_spectra = mueller.spectrum(E, el)     # [#nu/fission]
        spe[el] = c_section * m_spectra         # [(#nu*m^2)/fission]
    return spe

def spe_reactors(data, E, tau=3.1536e7):
    # Calculation of the spectrum as produced by every reactors
    # data = reactors data
    # E = Energy Axis
    # tau = 3.1536e7 #s = 1 year acquisition time

    spe_dict = dict()

    from tqdm.notebook import tqdm
    for i, row in tqdm(data.iterrows(), total=data.shape[0]):
        spe = spe_reactor(E, float(row['Pth [MW]']), row['LF'] / 100., row['Fuel'], tau)  # [#nu*m^2/10 keV]
        spe_dict[i] = spe

    return spe_dict

def spe_reactor(E, P, lf, fuel, tau = 3.1536e7):

    # Calculation of the spectrum as produced by one reactors
    # E = Energy Axis
    # P = Reactor Thermal Power in MW
    # lf = load factor
    # fuel = which fuel operate the reactor
    # tau = 3.1536e7 #s = 1 year acquisition time

    ###########################Norm factor due to reactor###########################

    Pn = P*1e6*6.241509e12 # Thermal Power from MW to MeV/s
    N_norm_reactor = tau*Pn*lf #[s*MeV/s] = [MeV]

    #print('-------> Constants_reactor')
    #print('-------> <Energy Axis>[MeV]: ', E.mean())
    #print('-------> P[MW]: ', P)
    #print('-------> Pn[MeV/s]: ', Pn)
    #print('-------> tau[s]: ', tau)
    #print('-------> lf: ', lf)
    #print('------->norm_reactor [MeV]: ', N_norm_reactor)

    spe = dict()

    # Get convolved Muller antineutrino spectra with IBD cross section for each element
    mu_spe = build_spectra(E) #[(#nu*m^2)/(10 keV*fission)]]

    p = fuels.p # Power fraction from fuels.py
    Q = fuels.Q # [MeV/fission] Energy released per fission from fuels.py

    p_fuels = p[fuel]

    for item in p_fuels.items():
        tot_spe = np.zeros(len(E))
        p_fuel = item[1]
        spe[item[0]] = dict()

        for el in mu_spe:
            norm = p_fuel[el]/Q[el]  # [fission/MeV]
            #print('-------> p_fuel[el]: ', p_fuel[el])
            #print('-------> Q[el]: ', Q[el])
            #print('-------> norm: ', norm)
            #print('-------> <mu_spe>: ', mu_spe[el].mean())
            temp_spe = mu_spe[el]*norm*N_norm_reactor  # [(#nu*m^2)/(10 keV*fission)]*[fission/MeV]*[MeV] = [#nu*m^2/10 keV]
            tot_spe = tot_spe+temp_spe
            spe[item[0]][el] = temp_spe
        spe[item[0]]['Tot'] = tot_spe
    return spe

def spectrum(data, E, det_lat, det_lon, det_depth, tau = 3.1536e7, hierarchy = 'NH', inMatter = False):
    # Calculation of spectrum for all reactors [Marica Baldoncini, et al. Phys. Rev. D 91, 065002] Eq. 15
    #   spe_dict=dict()
    #           |
    #           |
    #       [reactors]
    #           |
    #           |-->[no_burn_up]or[fuel compos. variation over time]
    #                               |
    #                               |-->[U235][U238][Pu239][Pu241]
    # data = reactors data
    # E = Energy Axis
    # det_lat, det_lon, det_depth = detector latitude, longitude and depth
    # tau = 3.1536e7 #s = 1 year acquisition time
    # hierarchy = 'NH' if Normal Hierarchy, 'IH' if inverse hierarchy
    # equation =  which equation to be used (from [F. Capozzi, E. Lisi, and A. Marrone Phys. Rev. D 89, 013001 – Published 9 January 2014])
    # inMatter = WIP
    # Calculation of the spectrum at Detector location due to different reactors
    # det_lat, det_lon, det_depth = detector latitude, longitude and depth
    # hierarchy = 'NH' if Normal Hierarchy, 'IH' if inverse hierarchy
    # all_spectra = if True return all the spectra due to different reactors

    spe_dict = dict()

    from tqdm.notebook import tqdm
    for i, row in tqdm(data.iterrows(), total=data.shape[0]):

        #Get spectrum for single reactor
        spe = spectrum_single(E, det_lat, det_lon, det_depth, float(row['Lat']), float(row['Lon']), float(row['Pth [MW]']), row['LF']/100., row['Fuel'], tau, hierarchy, inMatter)
        spe_dict[row['Name']]=spe

    return spe_dict


def spectrum_single(E, det_lat, det_lon, det_depth, rea_lat, rea_lon, P, lf, fuel, tau = 3.1536e7, hierarchy = 'NH', inMatter = False):
    
    # Calculation of spectrum for a single reactor [Marica Baldoncini, et al. Phys. Rev. D 91, 065002] Eq. 15
    #   spe = dict()
    #    |
    #    |-->[no_burn_up]or[fuel compos. variation over time]
    #                       |
    #                       |-->[U235][U238][Pu239][Pu241]
    
    # E = Energy Axis
    # det_lat, det_lon, det_depth = detector latitude, longitude and depth
    # rea_lat, rea_lon = reactor latitude and longitude
    # P = Reactor Thermal Power in MW
    # lf = load factor
    # fuel = which fuel operate the reactor
    # tau = 3.1536e7 #s = 1 year acquisition time
    # hierarchy = 'NH' if Normal Hierarchy, 'IH' if inverse hierarchy
    # equation =  which equation to be used (from [F. Capozzi, E. Lisi, and A. Marrone Phys. Rev. D 89, 013001 – Published 9 January 2014])
    # inMatter = WIP

    spe = dict()

    ###########################Detector constants from papers##########################

    #distance = 52474 # m [F. Capozzi, E. Lisi, and A. Marrone Phys. Rev. D 89, 013001 – Published 9 January 2014]
    eff = 1 #Efficiency
    targets = 1e32 #total number of free protons
    
    ############################Geometrical factor###########################
    distance = proj.dist([proj.earth_radius-det_depth, det_lon, det_lat], [proj.earth_radius, rea_lon, rea_lat]) #get distance between detector and reactor in m
    geom = 1/(4*np.pi*((distance)**2)) # get Geom factor of the antineutrino signal equation [m^-2]

    ###########################Power fraction from fuels.py###########################
    p = fuels.p
    p_fuels = p[fuel]

    ############################Survival probability###########################

    Pee = oscillations.Pee_39(distance, E, hierarchy, inMatter) #Get oscillation part of the antineutrino signal equation


    ###########################Norm factor due to detector###########################
    N_norm_detector = eff * targets * geom  # [#prot/m^2]

    #print('Constants')
    #print('Efficiency: ', eff)
    #print('Targets: ', targets)
    #print('Distance [m]: ', distance)
    #print('Geom_factor [m^-2]: ', geom)
    #print('norm_detector [prot*m^-2]: ', N_norm_detector)
    #print('Average Pee: ', Pee.mean())



    ###########################Spectra for each reactors###########################
    spe_react = spe_reactor(E, P, lf, fuel, tau)  # [#nu*m^2/10 keV] #


    for item in p_fuels.items(): # for loop for different fuel composition in time
        tot_spe = np.zeros(len(E))
        spe[item[0]]=dict()
        spe[item[0]]['dist'] = distance
        
        for el in mueller.elements:
            t_spe_ractor = spe_react[item[0]][el] #[#nu*m^2/10 keV] single element spectrum at reactor
            t_spe = t_spe_ractor*Pee*N_norm_detector # [#nu*m^2/10 keV] * [#prot/m^2] = [#nu*#prot/10 keV]  single element spectrum at detector
            tot_spe = tot_spe+t_spe
            spe[item[0]][el] = t_spe
        spe[item[0]]['Tot'] = tot_spe



    return spe





    

def get_Pee(E, distance, hierarchy = 'NH', inMatter = False):
    
    # Return electron neutrino oscillation survival probability in vacuum [F. Capozzi, E. Lisi, and A. Marrone Phys. Rev. D 89, 013001 – Published 9 January 2014]
    # E = Energy Axis
    # Distance = Distance between detectro and reactor in m
    # hierarchy = 'NH' if Normal Hierarchy, 'IH' if inverse hierarchy


    Pee = oscillations.Pee_39(distance, E, hierarchy, inMatter) #Get oscillation part of the antineutrino signal equation


    return Pee



def get_ratio(E, distance, hierarchy = 'NH', inMatter=False):
    
    # Return electron neutrino oscillation survival probability ratio of the phases in vacuum [F. Capozzi, E. Lisi, and A. Marrone Phys. Rev. D 89, 013001 – Published 9 January 2014]
    # E = Energy Axis
    # Distance = Distance between detectro and reactor in m
    # hierarchy = 'NH' if Normal Hierarchy, 'IH' if inverse hierarchy

    Pee = oscillations.ratioPhase(distance, E, hierarchy, inMatter=False) #Get oscillation part of the antineutrino signal equation

    return Pee



