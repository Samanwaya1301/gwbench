import numpy as np
import sympy as sp

from gwbench.basic_constants import time_fac, REarth, AU, cLight

cos = sp.cos
sin = sp.sin
exp = sp.exp
PI = np.pi

ap_symbs_string = 'f Mc tc ra dec psi gmst0'
f, Mc, tc, ra, dec, psi, gmst0 = sp.symbols(ap_symbs_string, real=True)

locs = ('H', 'L', 'V', 'K', 'I', 'ET1', 'ET2', 'ET3', 'C', 'N', 'S')

def detector_response_expr(hf_pl,hf_cr,loc,use_rot):
    # input:    loc     location (and implied orientation) of a detector
    #           use_rot  use frequency dependent time due to roa=tation of earth and SPA
    #
    # output:   hf
    return detector_response(f,hf_pl,hf_cr,Mc,tc,ra,dec,psi,gmst0,loc,use_rot)

def antenna_pattern_and_loc_phase_fac_expr(loc,use_rot):
    # input:    loc     location (and implied orientation) of a detector
    #           use_rot  use frequency dependent time due to roa=tation of earth and SPA
    #
    # output:   Fp, Fc
    return antenna_pattern_and_loc_phase_fac(f,Mc,tc,ra,dec,psi,gmst0,loc,use_rot)

def detector_response(f,hf_pl,hf_cr,Mc,tc,ra,dec,psi,gmst0,loc,use_rot):
    # input:    f       frequency domain [Hz]
    #           Mc      chirp Mass [solar mass]
    #           tc      time of coalescence [s]
    #           dec     declination [rad]
    #           ra      right ascencsion [rad]
    #           psi     polarization angle [rad]
    #           gmst0   GreenwichMeanSiderialTime according to LAL
    #           loc     location (and implied orientation) of a detector
    #           use_rot  use frequency dependent time due to rotation of earth and SPA
    #
    # output:   hf      detector strain

    Fp, Fc, Flp = antenna_pattern_and_loc_phase_fac(f,Mc,tc,ra,dec,psi,gmst0,loc,use_rot)
    return Flp * (Fp * hf_pl + Fc * hf_cr)

def antenna_pattern_and_loc_phase_fac(f,Mc,tc,ra,dec,psi,gmst0,loc,use_rot): # theta=dec; phi=ra
    # input:    f       frequency domain [Hz]
    #           Mc      chirp Mass [solar mass]
    #           tc      time of coalescence [s]
    #           dec     declination [rad]
    #           ra      right ascencsion [rad]
    #           psi     polarization angle [rad]
    #           gmst0   GreenwichMeanSiderialTime according to LAL
    #           loc     location (and implied orientation) of a detector
    #           use_rot  use frequency dependent time due to rotation of earth and SPA
    #
    # output:   Fp, Fc

    half_period = 4.32e4
    R = REarth

    D, d = det_ten_and_loc_vec(loc, R)

    if use_rot:
        tf = tc - (5./256.)*(time_fac*Mc)**(-5./3.)*(PI*f)**(-8./3.)
    else:
        tf = 0

    gra = (gmst0 + tf*PI/half_period) - ra
    theta = PI/2. - dec

    r = sp.Matrix([cos(gra) * sin(theta), sin(gra) * sin(theta), cos(theta)])
    XX = sp.Matrix([ -cos(psi)*sin(gra) - sin(psi)*cos(gra)*sin(dec), -cos(psi)*cos(gra) + sin(psi)*sin(gra)*sin(dec), sin(psi)*cos(dec) ])
    YY = sp.Matrix([  sin(psi)*sin(gra) - cos(psi)*cos(gra)*sin(dec),  sin(psi)*cos(gra) + cos(psi)*sin(gra)*sin(dec), cos(psi)*cos(dec) ])
    Fp = 0.5 * (XX.T*D*XX - YY.T*D*YY)
    Fc = 0.5 * (XX.T*D*YY + YY.T*D*XX)

    return Fp[0,0], Fc[0,0], exp(1j * 2*PI * f * r.T*d)[0,0]

def det_ten_and_loc_vec(loc, R):
    i_vec = sp.Matrix([1,0,0])
    j_vec = sp.Matrix([0,1,0])
    k_vec = sp.Matrix([0,0,1])

    et_vec2 = ( i_vec + np.sqrt(3.)*j_vec)/2.
    et_vec3 = (-i_vec + np.sqrt(3.)*j_vec)/2.

    alpha, beta, gamma = det_angles(loc)
    EulerD1 = rot_mat(alpha,'k') * rot_mat(beta,'j') * rot_mat(gamma,'k')

    if loc in   ('ET3','LISA3'):
        eDArm1 = -1 * EulerD1*et_vec2
        eDArm2 = -1 * EulerD1*et_vec3
    elif loc in ('ET2','LISA2'):
        eDArm1 =      EulerD1*et_vec3
        eDArm2 = -1 * EulerD1*i_vec
    elif loc in ('ET1','LISA1'):
        eDArm1 =      EulerD1*i_vec
        eDArm2 =      EulerD1*et_vec2
    else:
        eDArm1 = EulerD1*i_vec
        eDArm2 = EulerD1*j_vec

    return eDArm1*eDArm1.T - eDArm2*eDArm2.T, R/cLight * EulerD1*k_vec

def rot_mat(angle,axis):
    c = np.cos(angle)
    s = np.sin(angle)

    if axis == 'i':
        return sp.Matrix( [ [1,0,0], [0,c,-s], [0,s,c] ] )
    if axis == 'j':
        return sp.Matrix( [ [c,0,s], [0,1,0], [-s,0,c] ] )
    if axis == 'k':
        return sp.Matrix( [ [c,-s,0], [s,c,0], [0,0,1] ] )

def det_angles(loc):
    # return alpha, beta, gamma in radians
    # alpha ... longitude
    # beta  ... pi/2 - latitude
    # gamma ... angle from 'Due East' to y-arm
    if loc == 'H':
        return -2.08406, PI/2.-0.810795, PI-5.65488
    elif loc == 'L':
        return -1.58431, PI/2.-0.533423, PI-4.40318
    elif loc in ('V','ET1','ET2','ET3'):
        return 0.183338, PI/2.-0.761512, PI-0.33916
    elif loc == 'K':
        return 2.3942, PI/2.-0.632682, PI-1.054113
    elif loc == 'I':
        return 1.334013, PI/2.-0.248418, PI-1.570796

    elif loc == 'C':
        return -1.969174, PI/2.-0.764918, 0.
    elif loc == 'N':
        return -1.8584265, PI/2.-0.578751, -PI/3.
    elif loc == 'S':
        return 2.530727, PI/2.+0.593412, PI/4.
