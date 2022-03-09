E = np.arange(5, 10000, 10)
spe = build_spectra(E)
sign = dict()
for el in mueller.a:
    sign[el] = spe[el].sum()


def signal(det_lat, det_lon, det_depth, rea_lat, rea_lon, P, lf, fuel):
    En = P * 1e6 * tau * 6.241509e12  # MeV
    distance = proj.dist([proj.earth_radius - det_depth, det_lon, det_lat], [proj.earth_radius, rea_lon, rea_lat])
    geom = 1 / (4 * np.pi * ((distance * 1e2) ** 2))
    N_norm = targets * eff * En * geom * lf

    Pee = 0.55

    p = fuels.p
    Q = fuels.Q
    tot_sign = 0

    for el in spe:
        norm = p[fuel][el] / Q[el]
        temp_sign = sign[el] * Pee * norm * N_norm
        tot_sign = tot_sign + temp_sign
    return tot_sign


def build_map(data, resolution=1):
    lat_step = lon_step = resolution
    lat_min = -89.5  # [deg]
    lat_max = +89.5  # [deg]

    lon_min = -179.5  # [deg]
    lon_max = +180.0  # [deg]

    lat_list = np.arange(lat_min, lat_max + lat_step, lat_step)
    lon_list = np.arange(lon_min, lon_max, lon_step)

    lon_lat_grid = np.meshgrid(lat_list, lon_list)
    lon_lat_grid = np.array(lon_lat_grid).T
    map_shape = lon_lat_grid.shape[:2]
    pre_coords = lon_lat_grid
    lon_lat_grid = lon_lat_grid.reshape(-1, 2)
    lat_grid = lon_lat_grid[:, 0]
    lon_grid = lon_lat_grid[:, 1]

    from tqdm.notebook import tqdm
    FER = np.zeros(len(lon_grid))
    for i, row in tqdm(data.iterrows(), total=data.shape[0]):
        #     reactor = row['Kind']
        s = signal(lat_grid, lon_grid, 0, float(row['Lat']), float(row['Lon']), float(row['Pth [MW]']),
                   row['LF'] / 100., row['Fuel'])
        FER = FER + s

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.io.img_tiles as cimgt

    LER = dict(Lon=lon_grid, Lat=lat_grid, TNU=np.asarray(FER))
    # my_proj = ccrs.EqualEarth()
    my_proj = ccrs.PlateCarree()
    scale = 3
    fig = plt.figure(figsize=(scale * 6., scale * 2.5))
    ax = plt.axes(projection=my_proj)
    ax.stock_img()

    lat_mesh, lon_mesh = pre_coords[:, :, 0], pre_coords[:, :, 1]
    LER_mesh = np.asarray(FER).reshape(map_shape)

    from matplotlib.colors import LogNorm

    img_extent = (-179.5, 179.5, -89.5, 89.5)
    img = ax.imshow(LER_mesh[::-1][::1], extent=img_extent, transform=my_proj, alpha=0.7,
                    norm=LogNorm(vmin=1e0, vmax=1e4))

    # fig.tight_layout()
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.03), ncol=4, frameon = False)
    # ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.015), ncol=4, frameon = False)

    fig.colorbar(img, label='S [TNU]')

    fig.tight_layout()
    fig.savefig('map.svg')
    fig.savefig('map.png', dpi=300)


def readDB2(filename):
    head = ['Name', 'isMox', 'isCandu', 'Lat', 'Lon', 'Pth [MW]', 'Jan',
            'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov',
            'Dec']
    data = pd.read_csv(filename, sep='\t', names=head, header=0)

    fuels = np.asarray(['STD'] * data.shape[0])
    fuels[data['isMox'] == 1] = 'MOX'
    fuels[data['isCandu'] == 1] = 'CAN'
    data['Fuel'] = fuels
    mean_load = 0
    tot_days = 0
    for month in days:
        mean_load = mean_load + data[month] * days[month]
        tot_days = tot_days + days[month]

        data['LF'] = mean_load / tot_days

    return data


######### Electron neutrino oscillation survival probability in vacuum [F. Capozzi, E. Lisi, and A. Marrone Phys. Rev. D 89, 013001 – Published 9 January 2014] Eq. 37 ##########

def Pee3_vac_37(d, E_nu, hierarchy='NH'):  # m, keV
    # Return electron neutrino oscillation survival probability in vacuum [F. Capozzi, E. Lisi, and A. Marrone Phys. Rev. D 89, 013001 – Published 9 January 2014] Eq. 37
    #
    # d = distance
    # E_nu =  Neutrino Energy axis
    # hierarchy = 'NH' if Normal Hierarchy, 'IH' if inverse hierarchy

    params = getParams(hierarchy)
    if hierarchy == 'NH':
        alpha = 1.
    elif hierarchy == 'IH':
        alpha = -1.

    eMev = E_nu  # Neutrino Energuy in MeV

    Delta_small = conv_factor * params['dm2_21'] * d / eMev  # dm = dm2_21
    Deltam = np.absolute(params['Dm2_32'] + 0.5 * params['dm2_21'])  # Dm = 1/2 * (Dm2_31 + Dm2_32)
    Delta_big = conv_factor * Deltam * d / eMev

    s2_12 = params['s2_12']
    c2_12 = params['c2_12']
    s2_13 = params['s2_13']
    c2_13 = params['c2_13']

    c4_13 = c2_13 ** 2

    pezzo1 = 4 * c4_13 * s2_12 * c2_12 * (np.sin(Delta_small)) ** 2

    arg2 = alpha * Delta_big + Delta_small / 2
    pezzo2 = 4 * s2_13 * c2_13 * c2_12 * (np.sin(arg2)) ** 2

    arg3 = -alpha * Delta_big + Delta_small / 2
    pezzo3 = 4 * s2_13 * c2_13 * s2_12 * (np.sin(arg3)) ** 2

    prob = 1 - pezzo1 - pezzo2 - pezzo3

    return prob


def Pee_37(d, E_nu, hierarchy='NH'):
    # Return electron neutrino oscillation survival probability
    # d = distance
    # E_nu =  Neutrino Energy axis
    # hierarchy = 'NH' if Normal Hierarchy, 'IH' if inverse hierarchy

    return Pee3_vac_37(d, E_nu, hierarchy)



################################################ Other Stuff #############################################################################à

# cos_phi = (c2_12*np.cos(2*s2_12*Delta_small) + s2_12*np.cos(2*c2_12*Delta_small))/eta
# sin_phi = (c2_12*np.sin(2*s2_12*Delta_small) - s2_12*np.sin(2*c2_12*Delta_small))/eta

# def arctan (cos_phi, sin_phi):
# phi=[]
# for c, s, in zip(cos_phi, sin_phi):
# if c < 0:
# p=np.arctan(s/c) + np.pi
# else:
# p=np.arctan(s/c)
# phi.append(p)
# phi=np.array(phi)
# return phi

# phi = arctan(cos_phi, sin_phi)