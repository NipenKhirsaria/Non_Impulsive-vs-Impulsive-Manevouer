from math import *

from variables_constants import *
from ode_solver import end_of_burn_states
from coe_from_sv import coe_from_sv
from Lagrange_coefficients import Lagrange_coefficients
from impulsive import impulsive


MU, RE, GRAV_ACC = constants()
m_0, T, Isp,h_p_2, h_a_2, v_rp_0_BII, v_vp_0_BII = variables()
vp_0 = np.linalg.norm(v_vp_0_BII)


### Burn 1 ###

v_r_end_of_burn_1_BII, v_v_end_of_burn_1_BII, m_1, t_burn_1 = end_of_burn_states(MU, GRAV_ACC, m_0, T, Isp,
                                                                                 v_rp_0_BII, v_vp_0_BII)

r_1 = np.linalg.norm(v_r_end_of_burn_1_BII)
v_1 = np.linalg.norm(v_v_end_of_burn_1_BII)

print('\nStates after end of first burn')
print('Spacecraft mass after',t_burn_1,'time is:',m_1,'kg')
print('Position vector after',t_burn_1,'time is:',np.dot(v_r_end_of_burn_1_BII,0.001),'km')
print('Magnitude of Position vector is:',np.dot(r_1,0.001),'km')
print('Velocity vector after',t_burn_1,'time is:',np.dot(v_v_end_of_burn_1_BII,0.001),'km/s')
print('Magnitude of Velocity vector is:',np.dot(v_1,0.001),'km/s')


r_parameters_1 = coe_from_sv(MU, v_r_end_of_burn_1_BII, v_v_end_of_burn_1_BII)
ta_1 = r_parameters_1[5]

if ta_1 <= pi:
   delta_theta_1 = pi - ta_1
else:
   delta_theta_1 = 3*pi - ta_1

v_ra_1_BII, v_va_1_BII = Lagrange_coefficients(v_r_end_of_burn_1_BII, v_v_end_of_burn_1_BII, delta_theta_1, MU)
ra_1 = np.linalg.norm(v_ra_1_BII)
va_1 = np.linalg.norm(v_va_1_BII)

apse_rot_1 = acos(np.dot(v_ra_1_BII, [-1, 0, 0])/ra_1)

print('\nAfter burn trajectory')
print('[h(10^10) in m/s^2,  e , RAAN (rad), inclination (rad),w (rad), TA (rad), a(km)]')
print(r_parameters_1)
print('Apse line rotation:', apse_rot_1*180/pi, 'degree')

print('\nAt Apogee of Orbit 1')
print('Position vector at apogee of orbit 1 is:',np.dot(v_ra_1_BII,0.001),'km')
print('Magnitude of Position vector is:',np.dot(ra_1,0.001),'km')
print('Velocity vector at apogee of orbit 1 is:',np.dot(v_va_1_BII,0.001),'km/s')
print('Magnitude of Velocity vector is:',np.dot(va_1,0.001),'km/s')
print('Velocity change at A (non-impulsive) is:',v_1-vp_0,'m/s')



### Burn 2 ###

v_r_end_of_burn_2_BII, v_v_end_of_burn_2_BII, m_2, t_burn_2 = end_of_burn_states(MU, GRAV_ACC, m_1, T, Isp,
                                                                                 v_ra_1_BII,v_va_1_BII)

r_2 = np.linalg.norm(v_r_end_of_burn_2_BII)
v_2 = np.linalg.norm(v_v_end_of_burn_2_BII)

print('\nStates after end of first burn')
print('Spacecraft mass after',t_burn_2,'time is:',m_2,'kg')
print('Position vector after',t_burn_2,'time is:',np.dot(v_r_end_of_burn_2_BII,0.001),'km')
print('Magnitude of Position vector is:',np.dot(r_2,0.001),'km')
print('Velocity vector after',t_burn_2,'time is:',np.dot(v_v_end_of_burn_2_BII,0.001),'km/s')
print('Magnitude of Velocity vector is:',np.dot(v_2,0.001),'km/s')


r_parameters_2 = coe_from_sv(MU ,v_r_end_of_burn_2_BII, v_v_end_of_burn_2_BII)
ta_2 = r_parameters_2[5]

if ta_2 <= pi:
   delta_theta_2 = pi - ta_2
else:
   delta_theta_2 = 3*pi - ta_2

v_ra_2_BII, v_va_2_BII = Lagrange_coefficients(v_r_end_of_burn_2_BII, v_v_end_of_burn_2_BII, delta_theta_2, MU)
ra_2 = np.linalg.norm(v_ra_2_BII)
va_2 = np.linalg.norm(v_va_2_BII)

apse_rot_2 = acos(np.dot(v_ra_2_BII, [-1, 0, 0]) / ra_2)

print('\nAfter burn trajectory')
print('[h(10^10) in m/s^2,  e , RAAN (rad), inclination (rad),w (rad), TA (rad), a(km)]')
print(r_parameters_2)
print('Apse line rotation:', apse_rot_2*180/pi, 'degree')

print('\nAt Apogee of Orbit 1')
print('Position vector at apogee of orbit 1 is:',np.dot(v_ra_2_BII,0.001),'km')
print('Magnitude of Position vector is:',np.dot(ra_2,0.001),'km')
print('Velocity vector at apogee of orbit 1 is:',np.dot(v_va_2_BII,0.001),'km/s')
print('Magnitude of Velocity vector is:',np.dot(va_2,0.001),'km/s')
print('Velocity change at B (non-impulsive) is:',v_2-va_1,'m/s')

print('Total mass used:',m_0-m_2,'kg')
print('Total velocity change:',v_2-va_1+v_1-vp_0,'m/s')


### Impulsive ###

r_parameters_0 = coe_from_sv(MU,v_rp_0_BII,v_vp_0_BII)
a = r_parameters_0[6]

h_p_0, h_a_0, delta_v_A,delta_v_B,delta_v_total,delta_m,v_rp_0_BII, v_ra_1_BII,\
v_ra_2_BII,v_vp_0_BII,v_va_1_BII,v_va_2_BII =impulsive(MU,RE,GRAV_ACC,m_0, Isp, h_p_2, h_a_2, v_rp_0_BII,a)

r_parameters_1 = coe_from_sv(MU,v_ra_1_BII,v_va_1_BII)
r_parameters_2 = coe_from_sv(MU,v_ra_2_BII,v_va_2_BII)

print('\nImulsive Manevour')
print('Initial Perigee altitude:',h_p_0,'km')
print('Initial Apogee altitude:',h_a_0,'km')
print('Velocity change at A:',delta_v_A,'m/s')
print('Velocity change at B:',delta_v_B,'m/s')
print('Total Velocity change:',delta_v_total,'m/s')
print('Propellant mass used:',delta_m,'kg')

print('\nOrbital Parameters')
print('[h(10^10) in m^2/s,  e , RAAN (rad), inclination (rad),w (rad), TA (rad), a(km)]')
print(r_parameters_0)
print(r_parameters_1)
print(r_parameters_2)