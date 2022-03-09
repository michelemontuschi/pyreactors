import numpy as np

earth_radius = 6371000 # [m]

def dist(point1, point2): # sphere
    r1, theta1, phi1 = point1 # radius, lon, lat
    r2, theta2, phi2 = point2 # radius, lon, lat

    x1 = r1 * np.cos(np.deg2rad(phi1)) * np.cos(np.deg2rad(theta1))
    y1 = r1 * np.cos(np.deg2rad(phi1)) * np.sin(np.deg2rad(theta1))
    z1 = r1 * np.sin(np.deg2rad(phi1))

    x2 = r2 * np.cos(np.deg2rad(phi2)) * np.cos(np.deg2rad(theta2))
    y2 = r2 * np.cos(np.deg2rad(phi2)) * np.sin(np.deg2rad(theta2))
    z2 = r2 * np.sin(np.deg2rad(phi2))

    dist = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    return dist

def dist_vec(r1, theta1, phi1, r2, theta2, phi2): # sphere
#     r1, theta1, phi1 = point1 # radius, lon, lat
#     r2, theta2, phi2 = point2 # radius, lon, lat

    x1 = r1 * np.cos(np.deg2rad(phi1)) * np.cos(np.deg2rad(theta1))
    y1 = r1 * np.cos(np.deg2rad(phi1)) * np.sin(np.deg2rad(theta1))
    z1 = r1 * np.sin(np.deg2rad(phi1))

    x2 = r2 * np.cos(np.deg2rad(phi2)) * np.cos(np.deg2rad(theta2))
    y2 = r2 * np.cos(np.deg2rad(phi2)) * np.sin(np.deg2rad(theta2))
    z2 = r2 * np.sin(np.deg2rad(phi2))

    dist = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    return dist

def distell(point1, point2): # sphere
    r1, theta1, phi1 = point1 # radius, lon, lat
    r2, theta2, phi2 = point2 # radius, lon, lat

    # x1 = r1 * np.cos(np.deg2rad(phi1)) * np.cos(np.deg2rad(theta1))
    # y1 = r1 * np.cos(np.deg2rad(phi1)) * np.sin(np.deg2rad(theta1))
    # z1 = r1 * np.sin(np.deg2rad(phi1))
    #
    # x2 = r2 * np.cos(np.deg2rad(phi2)) * np.cos(np.deg2rad(theta2))
    # y2 = r2 * np.cos(np.deg2rad(phi2)) * np.sin(np.deg2rad(theta2))
    # z2 = r2 * np.sin(np.deg2rad(phi2))
    #
    # dist = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    raise Exception('To implement')
    
    return dist
