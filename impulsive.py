from math import *
import numpy as np


def impulsive(MU,RE,GRAV_ACC,m_0, Isp, h_p_2, h_a_2, v_rp_0_BII,a):

    # Perigee & Apogee of orbit 0 (initial parking orbit)
    r_p_0 = np.linalg.norm(v_rp_0_BII)
    r_a_0 = 2*a*1000 - r_p_0

    h_p_0 = (r_p_0 - RE)/(10**3)
    h_a_0 = (r_a_0 - RE)/(10**3)

    # Perigee & Apogee of orbit 1 (intermediate orbit)
    r_p_1 = r_p_0
    r_a_1 = RE + h_a_2*(10**3)

    # Perigee & Apogee of orbit 2 (final target orbit)
    r_p_2 = RE + h_p_2*(10**3)
    r_a_2 = r_a_1

    r_p = [r_p_0,r_p_1,r_p_2]
    r_a = [r_a_0,r_a_1,r_a_2]

    r_h = []
    r_vA = []
    r_vB = []
    for i in range(len(r_p)):
      r_h.append(sqrt(2*MU*r_a[i]*r_p[i]/(r_a[i]+r_p[i])))    # Angular momentum

      if i == 0:
          r_vA.append(r_h[i]/r_p[i])            # Velocity at A, Perigee of orbit 0 & 1
      elif i == 1:
          r_vA.append(r_h[i]/r_p[i])
          r_vB.append(r_h[i]/r_a[i])
      elif i == 2:
          r_vB.append(r_h[i]/r_a[i])            # Velocity at B, Apogee of orbit 1 & 2

    v_rp_0_BII = [r_p_0,0,0]
    v_ra_1_BII = [-r_a_1,0,0]
    v_ra_2_BII = v_ra_1_BII

    v_vp_0_BII = [0, r_vA[0], 0]
    v_va_1_BII = [0, -r_vB[0], 0]
    v_va_2_BII = [0, -r_vB[1], 0]

    delta_v_A = r_vA[1] - r_vA[0]
    delta_v_B = r_vB[1] - r_vB[0]

    delta_v_total = delta_v_A + delta_v_B    # Total velocity change required

    delta_m = m_0*(1-exp(-delta_v_total/(GRAV_ACC*Isp)))

    return h_p_0, h_a_0, delta_v_A,delta_v_B,delta_v_total,delta_m,v_rp_0_BII,\
           v_ra_1_BII,v_ra_2_BII,v_vp_0_BII,v_va_1_BII,v_va_2_BII




