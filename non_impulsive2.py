import numpy as np
from math import *

from impulsive import impulsive
from Lagrange_coefficients import Lagrange_coefficients
from ode_solver2 import end_of_burn_states2
from coe_from_sv import coe_from_sv
from variables_constants import *
from non_impulsive import non_impulsive


def non_impulsive2():
   MU,RE,GRAV_ACC = constants()
   m_0,T,Isp,v_rp_0_BII,v_vp_0_BII = variables()

   v_r_2_BII,v_v_2_BII,m_2,t_burn_2 = end_of_burn_states2()
   r_2 = np.linalg.norm(v_r_2_BII)
   v_2 = np.linalg.norm(v_v_2_BII)

   r_coe = coe_from_sv(v_r_2_BII,v_v_2_BII)
   e = r_coe[1]
   ta = r_coe[5]
   a = r_coe[6]

   if ta <= 180:
      delta_theta = 180 - ta
   else:
      delta_theta = 3*180 - ta

   v_ra_2_BII,v_va_2_BII = Lagrange_coefficients(v_r_2_BII,v_v_2_BII,delta_theta,MU)
   ra_2 = np.linalg.norm(v_ra_2_BII)
   va_2 = np.linalg.norm(v_va_2_BII)

   apse_rot_2 = acos(np.dot(v_ra_2_BII, [-1, 0, 0])/ra_2)   # apse line rotation

  # v_ra_1_BII, v_va_1_BII,m_1,v_1 = non_impulsive()

   print('\nStates after end of second burn\n')
   print('Spacecraft mass after',t_burn_2,'time is:',m_2,'kg')
   print('Position vector after',t_burn_2,'time is:',np.dot(v_r_2_BII,0.001),'km')
   print('Magnitude of Position vector is:',np.dot(r_2,0.001),'km')
   print('Velocity vector after',t_burn_2,'time is:',np.dot(v_v_2_BII,0.001),'km/s')
   print('Magnitude of Velocity vector is:',np.dot(v_2,0.001),'km/s','\n')

   print('After burn trajectory\n')
   print('[h(10^10) in m/s^2,  e , RAAN (deg), inclination (deg),w (rad), TA (deg), a(km)]')
   print(r_coe)
   print('Apse line rotation:',apse_rot_2*180/pi,'degree','\n')

   print('At Apogee of Orbit 2\n')
   print('Position vector at apogee of orbit 2 is:',np.dot(v_ra_2_BII,0.001),'km')
   print('Magnitude of Position vector is:',np.dot(ra_2,0.001),'km')
   print('Velocity vector at apogee of orbit 2 is:',np.dot(v_va_2_BII,0.001),'km/s')
   print('Magnitude of Velocity vector is:',np.dot(va_2,0.001),'km/s')
   print('Velocity change at B (non-impulsive) is:',v_2-np.linalg.norm(v_va_1_BII),'m/s')
   print('Total mass used:',m_0-m_2,'kg')
   print('Total velocity change:',v_2-np.linalg.norm(v_va_1_BII)+v_1-np.linalg.norm(v_vp_0_BII))

   return v_ra_2_BII,v_va_2_BII,m_2

v_ra_1_BII, v_va_1_BII,m_1,v_1 = non_impulsive()
v_ra_2_BII,v_va_2_BII,m_2 = non_impulsive2()







