import numpy as np

p_std = dict(no_burn_up = dict(x=0, U235 = 0.60, Pu239 = 0.27, U238 = 0.07, Pu241 = 0.06)) # from [F. Capozzi, E. Lisi, and A. Marrone Phys. Rev. D 89, 013001 – Published 9 January 2014]

p_candu = dict(no_burn_up = dict(x=0, U235 = 0.543, Pu239 = 0.411, U238 = 0.024, Pu241 = 0.022))

p_mox_clean = dict(U235 = 0.00, Pu239 = 0.708, U238 = 0.081, Pu241 = 0.212)

p_mox = dict()

for el in ['U235', 'Pu239', 'U238', 'Pu241']:
    p_mox[el] = np.around(0.7*p_std['no_burn_up'][el] + 0.3*p_mox_clean[el], 5)

p_mox['x'] = 0
p_mox = dict(n_time = p_mox)
    
import pandas as pd
df= pd.read_csv('./Papers/Graph/An2017.csv')
p_BU=df.to_dict('index')
p = dict(STD = p_std, MOX = p_mox, CAN = p_candu, BU = p_BU)

    

Q = dict(U235 = 202.36, Pu239 = 211.12, U238 = 205.99, Pu241 = 214.26) # MeV
Q['mean'] = np.asarray(list(Q.values())).mean()

dQ = dict(U235 = 0.26, Pu239 = 0.34, U238 = 0.52, Pu241 = 0.33)
dQ['mean'] = np.asarray(list(dQ.values())).mean()
