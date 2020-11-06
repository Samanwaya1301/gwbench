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


import lalsimulation as lalsim
from numpy import exp, pi

from gwbench.basic_constants import Mpc, Msun
from gwbench.basic_relations import m1_m2_of_M_eta, M_of_Mc_eta

wf_symbs_string = 'f Mc eta chi1x chi1y chi1z chi2x chi2y chi2z DL tc phic iota'

def hfpc(f, Mc, eta, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, DL, tc, phic, iota, approximant, fRef=0., phiRef=0.):
    f_min = f[0]
    delta_f = f[1] - f[0]
    f_max = f[-1] + delta_f

    if not fRef: fRef = f_min

    m1, m2 = m1_m2_of_M_eta(M_of_Mc_eta(Mc,eta),eta)
    m1 *= Msun
    m2 *= Msun
    r = DL * Mpc

    approx = lalsim.GetApproximantFromString(approximant)

    hPlus, hCross = lalsim.SimInspiralChooseFDWaveform(m1=m1, m2=m2, 
                                   S1x = chi1x, S1y = chi1y, S1z = chi1z,
                                   S2x = chi2x, S2y = chi2y, S2z = chi2z,
                                   distance = r, inclination = iota, phiRef = phiRef, 
                                   longAscNodes=0., eccentricity=0., meanPerAno = 0.,
                                   deltaF=delta_f, f_min=f_min, f_max=f_max, f_ref=fRef,
                                   LALpars=None, approximant=approx)

    pf = exp(1j*(2*f*pi*tc - phic))
    i0 = int(round((f_min-hPlus.f0) / delta_f))

    hfp = pf *  hPlus.data.data[i0:i0+len(f)]
    hfc = pf * hCross.data.data[i0:i0+len(f)]

    return hfp, hfc
