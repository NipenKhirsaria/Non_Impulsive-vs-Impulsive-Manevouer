import numpy as np
from math import *

def Lagrange_coefficients(v_r_0_BII,v_v_0_BII,delta_theta,MU):
    v_h_I = np.linalg.norm(np.cross(v_r_0_BII,v_v_0_BII))
    vr_0 = np.dot(v_v_0_BII,v_r_0_BII)/np.linalg.norm(v_r_0_BII)
    r_0_norm = np.linalg.norm(v_r_0_BII)
    s = sin(delta_theta*pi/180)
    c = cos(delta_theta*pi/180)

    r = (v_h_I**2/MU)*(1/(1+(v_h_I**2/(MU*r_0_norm)-1)*c-v_h_I*vr_0*s/MU))

    f = 1 - MU*r*(1-c)/(v_h_I**2)
    g = r*r_0_norm*s/v_h_I

    g_dot = 1 - (MU*r_0_norm/(v_h_I**2))*(1-c)
    f_dot = (1/g)*(f*g_dot - 1)

    v_R_BII = np.dot(f,v_r_0_BII) + np.dot(g,v_v_0_BII)
    v_V_BII = np.dot(f_dot,v_r_0_BII) + np.dot(g_dot,v_v_0_BII)

    return v_R_BII,v_V_BII

