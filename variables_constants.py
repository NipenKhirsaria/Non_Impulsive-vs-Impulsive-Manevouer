import numpy as np

def constants():
    MU = 3.986*(10**14)   # Gravitational parameter for Earth (m^3/s^2)
    RE = 6378*(10**3)     # Radius of earth (m)
    GRAV_ACC = 9.807             # Gravitational acceleration (m/s^2)

    return MU,RE,GRAV_ACC

def variables():
    m_0 = 5 # 2000  #1112.4955
    T =  1400*(10**-3) # 10*(10**3)
    Isp = 8500 # 300
    #t_burn_1 = 261.1127
    #t_burn_2 = 118.88
    v_rp_0_BII = np.dot([200+6378,0,0],1000)    # np.dot([480+6378,0,0],1000)   #np.dot([-22141.56177,-3244.5289,0],1000)
    v_vp_0_BII =  np.dot([0,7.799,0],1000)    #  np.dot([0,7.7102,0],1000)     #[419.39,-2862.034,0]
    #h_p_2 = 16000
    #h_a_2 = 16000
    #Isp = 300
    #T = 10000

    return m_0, T, Isp, v_rp_0_BII, v_vp_0_BII

def for_impulsive():
    m_0 = 5   # 2000
    hp_0 = 200   # 480
    ha_0 = 250 # 800
    hp_2 = 345  # 16000
    ha_2 = 345  # 16000
    Isp = 8500  # 300

    return m_0,hp_0,ha_0,hp_2,ha_2,Isp

