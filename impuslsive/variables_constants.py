import numpy as np

def constants():
    MU = 3.986*(10**14)          # Gravitational parameter for Earth (m^3/s^2)
    RE = 6378*(10**3)            # Radius of earth (m)
    GRAV_ACC = 9.807             # Gravitational acceleration (m/s^2)

    return MU,RE,GRAV_ACC

def for_impulsive():
    m_0 = 5   # 2000             # Initial Mass of spacecraft
    hp_0 = 200   # 480           # Perigee altitude of initial orbit (km)
    ha_0 = 250 # 800             # Apogee altitude of initial orbit (km)
    hp_2 = 345  # 16000          # Perigee altitude of final target orbit (km)
    ha_2 = 345  # 16000          # Apogee altitude of final target orbit (km)
    Isp = 8500  # 300            # Specific Impulse (s)

    return m_0,hp_0,ha_0,hp_2,ha_2,Isp

