import numpy as np
from math import *


def coe_from_sv(MU,v_R_BII, v_V_BII):

    tol_e = 10**(-10)   # eccentricity limit to decide ellipse

    R = np.linalg.norm(v_R_BII)
    V = np.linalg.norm(v_V_BII)

    vr = np.dot(v_R_BII,v_V_BII)/R

    v_h_I = np.cross(v_R_BII,v_V_BII)
    h = np.linalg.norm(v_h_I)

    incl = acos(v_h_I[2] / h)        # Inclination angle
    v_n_I = np.cross([0, 0, 1],v_h_I)     # Line of Nodes
    n = np.linalg.norm(v_n_I)

    if n!= 0:
       RAAN = acos(v_n_I[0] / n)     # Right Ascension Angle of Node
       if v_n_I[1] < 0:
          RAAN = 2 * pi - RAAN
    else:
        RAAN = 0

    v_e_I = (1/MU)*(np.dot((V**2-MU/R),v_R_BII) - np.dot(R*vr,v_V_BII)) # Eccentricity Vector
    e = np.linalg.norm(v_e_I)

    if n != 0:
        if e > tol_e:
           w = acos(np.dot(v_n_I,v_e_I)/(n*e))    # Argument of Perigee
           if v_e_I[2] < 0:
              w = 2*pi - w
        else:
            w = 0
    else:
        w = 0

    if e > tol_e:
        ta = acos(np.dot(v_e_I,v_R_BII)/(e*R))    # True anomaly
        if vr < 0:
            ta = 2*pi - ta

    else:
        cp = np.cross(v_n_I,v_R_BII)
        if cp[2] >= 0:
            ta = acos(np.dot(v_n_I,v_R_BII)/(n*R))
        else:
            ta = 2*pi - acos(np.dot(v_n_I,v_R_BII)/(n*R))

    a = h**2/(MU*(1-e**2))      # semi - major axis

    r_parameters = [h*(10**-10),e,RAAN,incl,w,ta,a*0.001]

    return r_parameters

