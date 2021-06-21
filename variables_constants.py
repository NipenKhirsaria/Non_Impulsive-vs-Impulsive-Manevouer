import numpy as np


def constants():
    MU = 3.986*(10**14)          # Gravitational parameter for Earth (m^3/s^2)
    RE = 6378*(10**3)            # Radius of earth (m)
    GRAV_ACC = 9.807             # Gravitational acceleration (m/s^2)

    return MU,RE,GRAV_ACC

def variables():
    m_0 = 2000
    T = 10*(10**3) 
    Isp = 300 
    v_rp_0_BII = np.dot([480+6378,0,0],1000) 
    v_vp_0_BII = np.dot([0,7.7102,0],1000)  
    h_p_2 = 16000           # Perigee altitude of final target orbit (km)
    h_a_2 = 16000           # Apogee altitude of final target orbit (km)


    return m_0, T, Isp, h_p_2, h_a_2, v_rp_0_BII, v_vp_0_BII
