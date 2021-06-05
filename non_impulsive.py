import numpy as np
from math import *

from impulsive import impulsive
from Lagrange_coefficients import Lagrange_coefficients
from ode_solver import end_of_burn_states
from coe_from_sv import coe_from_sv
from variables_constants import *

def non_impulsive():
   MU,RE,GRAV_ACC = constants()
   m_0,T,Isp,v_rp_0_BII,v_vp_0_BII = variables()

   v_r_1_BII,v_v_1_BII,m_1,t_burn_1 = end_of_burn_states()
   r_1 = np.linalg.norm(v_r_1_BII)
   v_1 = np.linalg.norm(v_v_1_BII)

   r_coe = coe_from_sv(v_r_1_BII,v_v_1_BII)
   e = r_coe[1]
   ta = r_coe[5]
   a = r_coe[6]

   if ta <= 180:
      delta_theta = 180 - ta
   else:
      delta_theta = 3*180 - ta

   v_ra_1_BII,v_va_1_BII = Lagrange_coefficients(v_r_1_BII,v_v_1_BII,delta_theta,MU)
   ra_1 = np.linalg.norm(v_ra_1_BII)
   va_1 = np.linalg.norm(v_va_1_BII)

   apse_rot_1  = acos(np.dot(v_ra_1_BII,[-1,0,0]) / ra_1)

   print('\nStates after end of first burn\n')
   print('Spacecraft mass after',t_burn_1,'time is:',m_1,'kg')
   print('Position vector after',t_burn_1,'time is:',np.dot(v_r_1_BII,0.001),'km')
   print('Magnitude of Position vector is:',np.dot(r_1,0.001),'km')
   print('Velocity vector after',t_burn_1,'time is:',np.dot(v_v_1_BII,0.001),'km/s')
   print('Magnitude of Velocity vector is:',np.dot(v_1,0.001),'km/s','\n')

   print('After burn trajectory\n')
   print('[h(10^10) in m/s^2,  e , RAAN (deg), inclination (deg),w (deg), TA (deg), a(km)]')
   print(r_coe)
   print('Apse line rotation:',apse_rot_1*180/pi,'degree','\n')

   print('At Apogee of Orbit 1\n')
   print('Position vector at apogee of orbit 1 is:',np.dot(v_ra_1_BII,0.001),'km')
   print('Magnitude of Position vector is:',np.dot(ra_1,0.001),'km')
   print('Velocity vector at apogee of orbit 1 is:',np.dot(v_va_1_BII,0.001),'km/s')
   print('Magnitude of Velocity vector is:',np.dot(va_1,0.001),'km/s')
   print('Velocity change at A (non-impulsive) is:',v_1-np.linalg.norm(v_vp_0_BII),'m/s')

   return v_ra_1_BII,v_va_1_BII,m_1,v_1






