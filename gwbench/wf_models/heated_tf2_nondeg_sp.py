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
import sympy as sp

import gwbench.basic_relations as brs
from gwbench.basic_constants import time_fac, strain_fac

cos = sp.cos
sin = sp.sin
log = sp.log
PI = np.pi

wf_symbs_string = 'f Mc eta chi1z chi2z DL tc phic iota Heff5 Heff8'
f, Mc, eta, chi1z, chi2z, DL, tc, phic, iota, Heff5, Heff8 = sp.symbols(wf_symbs_string, real=True)

#------from Anruadha--------
# defining constants
GammaE = 0.577215664901532

def hfpc(f, Mc, eta, chi1z, chi2z, DL, tc, phic, iota, Heff5, Heff8):
    '''
    Mc ... in solar mass
    DL ... in mega parsec
    '''
    # convert to sec
    Mc = Mc * time_fac
    DL = DL * time_fac/strain_fac

    # get sym and asym chi combinations
    chi_s = brs.chi_s(chi1z,chi2z)
    chi_a = brs.chi_a(chi1z,chi2z)

    '''
    Mc is in sec, e.g., Mc = 10*MTSUN_SI (for 10 solar mass)
    DL is in sec, e.g., DL = 100*1e6*PC_SI/C_SI (for 100 Mpc)
    '''

    M = Mc/eta**(3./5.)
    delta = (1.-4.*eta)**0.5
    v  = (PI*M*f)**(1./3.)
    flso = brs.f_isco(M)
    vlso = (PI*M*flso)**(1./3.)
    A =((5./24.)**0.5/PI**(2./3.))*(Mc**(5./6.)/DL)

    # 3.5PN phasing (point particle limit)
    p0 = 1.

    p1 = 0

    p2 = (3715./756. + (55.*eta)/9.)

    p3 = (-16.*PI + (113.*delta*chi_a)/3. + (113./3. - (76.*eta)/3.)*chi_s)

    p4 = (15293365./508032. + (27145.*eta)/504.+ (3085.*eta**2)/72. + (-405./8. + 200.*eta)*chi_a**2 - (405.*delta*chi_a*chi_s)/4. + (-405./8. + (5.*eta)/2.)*chi_s**2)

    gamma = (732985./2268. - 24260.*eta/81. - 340.*eta**2/9.)*chi_s + (732985./2268. + 140.*eta/9.)*delta*chi_a

    p5 = (38645.*PI/756. - 65.*PI*eta/9. - gamma)

    p5L = (38645.*PI/756. - 65.*PI*eta/9. - gamma)*3*log(v/vlso)

    p6 = (11583231236531./4694215680. - 640./3.*PI**2 - 6848./21.*GammaE + eta*(-15737765635./3048192. + 2255./12.*PI**2) + eta*eta*76055./1728. - eta*eta*eta*127825./1296. \
         - (6848./21.)*log(4.) + PI*(2270.*delta*chi_a/3. + (2270./3. - 520.*eta)*chi_s) + (75515./144. - 8225.*eta/18.)*delta*chi_a*chi_s \
         + (75515./288. - 263245.*eta/252. - 480.*eta**2)*chi_a**2 + (75515./288. - 232415.*eta/504. + 1255.*eta**2/9.)*chi_s**2)

    p6L = -(6848./21.)*log(v)

    p7 = (((77096675.*PI)/254016. + (378515.*PI*eta)/1512.- (74045.*PI*eta**2)/756. + (-25150083775./3048192. + (10566655595.*eta)/762048. - (1042165.*eta**2)/3024. + (5345.*eta**3)/36.
         + (14585./8. - 7270.*eta + 80.*eta**2)*chi_a**2)*chi_s + (14585./24. - (475.*eta)/6. + (100.*eta**2)/3.)*chi_s**3 + delta*((-25150083775./3048192.
         + (26804935.*eta)/6048. - (1985.*eta**2)/48.)*chi_a + (14585./24. - 2380.*eta)*chi_a**3 + (14585./8. - (215.*eta)/2.)*chi_a*chi_s**2)))

    phase = 2*f*PI*tc - phic - PI/4. + (3./(128.*v**5*eta))*(p0 + v*p1 + v**2*p2 + v**3*p3+ v**4*p4 + v**5*(p5+p5L) + v**6*(p6+p6L) + v**7*p7)

    #phase due to tidal heating
    #--------------------------------------------------------------------------------
    psi_so1 = (1/6.)*(-56*eta - 73*np.sqrt(1 - 4*eta) + 73)*chi1z
    psi_so2 = (1/6.)*(-56*eta - 73*np.sqrt(1 - 4*eta) + 73)*chi2z
    psi_so = psi_so1 + psi_so2
    con = (3./(128.*eta))/v**5
    term1 = -(10/9.)*(v**5)*Heff5*(3*np.log(v))
    term2 = -(5/168.)*(v**7)*Heff5*(952*eta + 995)
    term3 = (5/9.)*(v**8)*(3*np.log(v))*(-4*Heff8 + Heff5*psi_so)
    
    heated_phase = con*(term1 + term2 + term3)
    phase += heated_phase
    
    
    hp = 0.5*(1+(cos(iota))**2)*A*f**(-7./6.)*(cos(phase) - 1j*sin(phase))
    hc = -1j*cos(iota)*A*f**(-7./6.)*(cos(phase) - 1j*sin(phase))

    return hp, hc
