# Copyright (C) 2020  Ssohrab Borhanian
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


import numpy as np

from gwbench.basic_constants import time_fac

PI = np.pi

#-----f_isco from Anuradha-----
def f_isco(M):
    '''
    M is in sec, e.g, M = 10*MTSUN_SI (for 10 solar mass)
    '''
    return 1./6.**(3./2.)/PI/M

def f_isco_Msolar(M):
    '''
    M ... in solar mass
    '''
    # convert to sec
    M = M * time_fac
    '''
    M is in sec, e.g, M = 10*MTSUN_SI (for 10 solar mass)
    '''
    return 1./6.**(3./2.)/PI/M

#-----mass ratio functions-----
def eta_of_q(q):
    return q/np.power(1+q,2)

def delta_of_q(q):
    return np.sqrt(1-4*eta_of_q(q))

def delta_of_eta(eta):
    return np.sqrt(1-4*eta)

def q_of_eta(eta):
    return (1+delta_of_eta(eta))/(1-delta_of_eta(eta))

def M_of_Mc_eta(Mc,eta):
    return Mc/np.power(eta,3./5)

def Mc_of_M_eta(M,eta):
    return M*np.power(eta,3./5)

def m1_m2_of_M_eta(M,eta):
    delta = delta_of_eta(eta)
    return 0.5*M*(1+delta), 0.5*M*(1-delta)

def m1_m2_of_Mc_eta(Mc,eta):
    return m1_m2_of_M_eta(M_of_Mc_eta(Mc,eta),eta)

def M_eta_of_m1_m2(m1,m2,q_gt_1=1):
    if q_gt_1:
        return m1+m2, eta_of_q(m1/m2)
    else:
        return m1+m2, eta_of_q(m2/m1)

def Mc_eta_of_m1_m2(m1,m2,q_gt_1=1):
    if q_gt_1:
        eta = eta_of_q(m1/m2)
    else:
        eta = eta_of_q(m2/m1)
    return Mc_of_M_eta(m1+m2,eta), eta

#-----spin ratio functions-----
def chi_s(chi1,chi2):
    return 0.5*(chi1+chi2)

def chi_a(chi1,chi2):
    return 0.5*(chi1-chi2)

def chi_eff(m1,m2,chi1,chi2):
    return (m1 * chi1 + m2 * chi2) / (m1+m2)

#------Heff defined as functions-----------(added by Samanwaya Mukherjee)
def Heff5(h1, h2, m1, m2, chi1, chi2):
    m = m1 + m2
    return h1 * (m1/m)**3 * (3 * chi1**2 + 1) * chi1 + h2 * (m2/m)**3 * (3 * chi2**2 + 1) * chi2 
def Heff8(h1, h2, m1, m2, chi1, chi2):
    m = m1 + m2
    return 4 * np.pi * Heff5(h1, h2, m1, m2, chi, chi2) + h1 * (m1/m)**4 * (1 + np.sqrt(1 - chi1**2)) * (1 + 3 * chi1**2) + h2 * (m2/m)**4 * (1 + np.sqrt(1 - chi2**2)) * (1 + 3 * chi2**2)
    

#-----derivatives of spin and mass functions-----
def del_Mc_M_of_eta(eta):
    return np.power(eta,-3./5)

def del_eta_M_of_Mc_eta(Mc,eta):
    return -3./5. * Mc * np.power(eta,-8./5)

def del_Mc_m1_of_Mc_eta(Mc,eta):
    delta = delta_of_eta(eta)
    return 1./2 * del_Mc_M_of_eta(eta) * (1 + delta)

def del_eta_m1_of_Mc_eta(Mc,eta):
    M = M_of_Mc_eta(Mc,eta)
    delta = delta_of_eta(eta)
    return 1./2 * del_eta_M_of_Mc_eta(Mc,eta) * (1 + delta) - M/delta

def del_Mc_m2_of_Mc_eta(Mc,eta):
    delta = delta_of_eta(eta)
    return 1./2 * del_Mc_M_of_eta(eta) * (1 - delta)

def del_eta_m2_of_Mc_eta(Mc,eta):
    M = M_of_Mc_eta(Mc,eta)
    delta = delta_of_eta(eta)
    return 1./2 * del_eta_M_of_Mc_eta(Mc,eta) * (1 - delta) + M/delta

def del_Mc_chi_eff(Mc,eta,chi1,chi2):
    M = M_of_Mc_eta(Mc,eta)
    m1, m2 = m1_m2_of_M_eta(M,eta)
    return -1./M * del_Mc_M_of_eta(eta) * chi_eff(m1,m2,chi1,chi2) + 1./M * (del_Mc_m1_of_Mc_eta(Mc,eta) * chi1 + del_Mc_m2_of_Mc_eta(Mc,eta) * chi2)

def del_eta_chi_eff(Mc,eta,chi1,chi2):
    M = M_of_Mc_eta(Mc,eta)
    m1, m2 = m1_m2_of_M_eta(M,eta)
    return -1./M * del_eta_M_of_Mc_eta(Mc,eta) * chi_eff(m1,m2,chi1,chi2) + 1./M * (del_eta_m1_of_Mc_eta(Mc,eta) * chi1 + del_eta_m2_of_Mc_eta(Mc,eta) * chi2)

def del_chi1_chi_eff(Mc,eta,chi1,chi2):
    M = M_of_Mc_eta(Mc,eta)
    m1, m2 = m1_m2_of_M_eta(M,eta)
    return 1./M * (m1 + m2 * chi2)

def del_chi2_chi_eff(Mc,eta,chi1,chi2):
    M = M_of_Mc_eta(Mc,eta)
    m1, m2 = m1_m2_of_M_eta(M,eta)
    return 1./M * (m2 + m1 * chi1)

# tidal parameters

def lam_ts_of_lam_12_eta(lam1,lam2,eta): # from arXiv:1402.5156
#    q = q_of_eta(eta)
#    lam_t = 16./13. * ( (12 + q) * q**4 * lam1 + (12*q + 1) * lam2) / (1 + q)**5
    delta = delta_of_eta(eta)
    lam_t = 8./13. * ( (1. + 7. * eta - 31. * eta**2) * (lam1 + lam2) +
                       delta * (1. + 9. * eta - 11. * eta**2) * (lam1 - lam2) )
    delta_lam_t = 0.5 * ( delta * (1319. - 13272. * eta + 8944. * eta**2) / 1319. * (lam1 + lam2)+
                          (1319. - 15910. * eta + 32850. * eta**2 + 3380. * eta**3) / 1319. * (lam1 - lam2) )
    return lam_t, delta_lam_t

def lam_12_of_lam_ts_eta(lam_t,delta_lam_t,eta):
    delta = delta_of_eta(eta)
    lam1 = ((-(-6.76923076923077*delta_lam_t*delta*(-0.09090909090909091 - 0.8181818181818182*eta + 1.*eta**2) +
            19.076923076923077*delta_lam_t*(-0.03225806451612903 - 0.22580645161290322*eta + 1.*eta**2) +
            3.3904473085670963*delta*(0.1474731663685152 - 1.4838998211091234*eta + 1.*eta**2)*lam_t -
            1.281273692191054*(0.39023668639053255 - 4.707100591715976*eta + 9.718934911242604*eta**2 + 1.*eta**3)*lam_t))/
            (8.881784197001252e-16*eta - 1.4210854715202004e-14*eta**2 + 2.842170943040401e-14*eta**3 + 4.500379075056848*eta**4 - 232.4912812736922*eta**5))
    lam2 = ((delta_lam_t*(-1.5296267736621122e-19 + 3.0592535473242243e-19*delta*eta + 7.342208513578138e-18*eta**2 + 9.789611351437518e-18*eta**3 +
            (0.011550173712335778 - 9.789611351437518e-18*delta)*eta**4 + 0.11646425159938675*eta**5) +
            (-3.8240669341552804e-20 + 3.8240669341552804e-20*delta + (-4.588880320986336e-19 - 6.118507094648449e-19*delta)*eta +
            (9.789611351437518e-18 - 1.2237014189296897e-18*delta)*eta**2 + (-9.789611351437518e-18 + 1.4684417027156276e-17*delta)*eta**3 +
            (-0.0014297928149279568 - 0.007954723326344955*delta)*eta**4 + (0.07386363636363634 - 0.005511061254304498*delta)*eta**5)*lam_t)/
            (eta*(0.0909090909090909 - 0.09090909090909091*delta + (0.6363636363636364 - 0.8181818181818181*delta)*eta +
            (-2.818181818181818 + 1.*delta)*eta**2)*(-3.82026549483612e-18 + 6.112424791737792e-17*eta - 1.2224849583475584e-16*eta**2 -
            0.01935719503287065*eta**3 + 1.*eta**4)))
    return lam1, lam2
