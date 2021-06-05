from math import *

from variables_constants import *

m_0,h_p_0,h_a_0,h_p_2,h_a_2,Isp = for_impulsive()
mu,RE,g = constants()

def impulsive():

    # Perigee & Apogee of orbit 0 (initial parking orbit)
    r_p_0 = RE + h_p_0*(10**3)
    r_a_0 = RE + h_a_0*(10**3)

    # Perigee & Apogee of orbit 1 (intermediate orbit)
    r_p_1 = r_p_0
    r_a_1 = RE + h_a_2*(10**3)

    # Perigee & Apogee of orbit 2 (final target orbit)
    r_p_2 = RE + h_p_2*(10**3)
    r_a_2 = r_a_1

    r_p = [r_p_0,r_p_1,r_p_2]
    r_a = [r_a_0,r_a_1,r_a_2]

    h = []
    v_A = []
    v_B = []
    for i in range(len(r_p)):
      h.append(sqrt(2*mu*r_a[i]*r_p[i]/(r_a[i]+r_p[i])))    # Angular momentum

      if i == 0:
          v_A.append(h[i]/r_p[i])            # Velocity at A, Perigee of orbit 0 & 1
      elif i == 1:
          v_A.append(h[i]/r_p[i])
          v_B.append(h[i]/r_a[i])
      elif i == 2:
          v_B.append(h[i]/r_a[i])            # Velocity at B, Apogee of orbit 1 & 2

    r_p_0_vec = [r_p_0,0,0]
    r_a_1_vec = [-r_a_1,0,0]
    r_a_2_vec = r_a_1_vec

    v_p_0_vec = [0, v_A[0], 0]
    v_a_1_vec = [0, -v_B[0], 0]
    v_a_2_vec = [0, -v_B[1], 0]

    delta_v_A = v_A[1] - v_A[0]
    delta_v_B = v_B[1] - v_B[0]

    delta_v_total = delta_v_A + delta_v_B    # Total velocity change required

    delta_m = m_0*(1-exp(-delta_v_total/(g*Isp)))

    r_p_0_vec = [r_p_0,0,0]
    v_p_0_vec = [0,v_A[0],0]


    return delta_v_A,delta_v_B,delta_v_total,delta_m,r_p_0_vec,r_a_1_vec,r_a_2_vec,v_p_0_vec,v_a_1_vec,v_a_2_vec

#delta_v_A,delta_v_B,delta_v_total,delta_m,r_p_0_vec,r_a_1_vec,r_a_2_vec,v_p_0_vec,v_a_1_vec,v_a_2_vec = impulsive()

