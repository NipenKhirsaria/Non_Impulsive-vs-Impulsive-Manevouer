import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from math import *

from variables_constants import *
from non_impulsive import non_impulsive

MU,RE,GRAV_ACC = constants()
m_0,T,Isp,v_rp_0_BII,v_vp_0_BII = variables()

def model2(Z,t,T):
        x = Z[0]
        y = Z[1]
        z = Z[2]
        x_dot = Z[3]
        y_dot = Z[4]
        z_dot = Z[5]
        m = Z[6]

        dxdt = x_dot
        dydt = y_dot
        dzdt = z_dot
        dx_dotdt = -MU*x/(sqrt(x**2 + y**2 + z**2))**3 + T*x_dot/(m *(sqrt(x_dot**2 + y_dot**2 + z_dot**2)))
        dy_dotdt = -MU*y/(sqrt(x**2 + y**2 + z**2))**3 + T*y_dot/(m *(sqrt(x_dot**2 + y_dot**2 + z_dot**2)))
        dz_dotdt = -MU*z/(sqrt(x**2 + y**2 + z**2))**3 + T*z_dot/(m *(sqrt(x_dot**2 + y_dot**2 + z_dot**2)))
        dmdt = -T/(Isp*GRAV_ACC)

        return [dxdt,dydt,dzdt,dx_dotdt,dy_dotdt,dz_dotdt,dmdt]


v_ra_1_BII,v_va_1_BII,m_1,v_1 = non_impulsive()
Z0 = [v_ra_1_BII[0],v_ra_1_BII[1],v_ra_1_BII[2],v_va_1_BII[0],v_va_1_BII[1],v_va_1_BII[2],m_1]

t_burn_2 = float(input('Enter the burn time in sec: '))
t = np.linspace(0,t_burn_2)

Z = odeint(model2,Z0,t,args=(T,))

x = Z[:,0]
y = Z[:,1]
z = Z[:,2]
x_dot = Z[:,3]
y_dot = Z[:,4]
z_dot = Z[:,5]
m = Z[:,6]

v_r_2_BII = [x[-1],y[-1],z[-1]]
v_v_2_BII = [x_dot[-1],y_dot[-1],z_dot[-1]]
m_2 = m[-1]

def end_of_burn_states2():
    return v_r_2_BII,v_v_2_BII,m_2,t_burn_2

v_r_2_BII,v_v_2_BII,m_2,t_burn_2 = end_of_burn_states2()