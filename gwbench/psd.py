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


import os

import scipy.interpolate as si
from numpy import power, logical_and, inf, pi, exp, tanh, cos, sin, square
from pandas import read_csv

from gwbench.basic_constants import cLight

noise_curves_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),'noise_curves')

tecs = ('aLIGO', 'A+', 'V+', 'K+', 'ET', 'CEwb', 
        'Voyager-CBO', 'Voyager-PMO',
        'LISA-17', 'LISA-Babak17', 'LISA-Robson18',
        'CE1-10-CBO', 'CE1-20-CBO', 'CE1-30-CBO', 'CE1-40-CBO',
        'CE2-10-CBO', 'CE2-20-CBO', 'CE2-30-CBO', 'CE2-40-CBO',
        'CE1-10-PMO', 'CE1-20-PMO', 'CE1-30-PMO', 'CE1-40-PMO',
        'CE2-10-PMO', 'CE2-20-PMO', 'CE2-30-PMO', 'CE2-40-PMO')

def psd(det,f,F_lo=-inf,F_hi=inf,psd_file=None,is_asd=None):

    # limit the PSD to a certain freq range (e.g. if detetor does have such a low PSD)
    f=f[logical_and(f>=F_lo,f<=F_hi)]

    # if a psd file is handed over, reset the det label
    if psd_file is not None: det = ''

    # use analytical PSDs, f_lo, f_hi will be that of freq array f
    if det=='aLIGO':
        f_lo = max(10.,F_lo)
        f_hi = min(2048.,F_hi)
        check_f(det,f,f_lo,f_hi)
        return  psd_aLIGO(f[logical_and(f>=f_lo,f<=f_hi)]), f[logical_and(f>=f_lo,f<=f_hi)]
    elif det == 'CEwb':
        f_lo = max(5.,F_lo)
        f_hi = min(2048.,F_hi)
        check_f(det,f,f_lo,f_hi)
        return  psd_CEwb(f[logical_and(f>=f_lo,f<=f_hi)]), f[logical_and(f>=f_lo,f<=f_hi)]
    elif det == 'LISA-17':
        f_lo = F_lo
        f_hi = F_hi
        check_f(det,f,f_lo,f_hi)
        return psd_LISA_17(f[logical_and(f>=f_lo,f<=f_hi)]), f[logical_and(f>=f_lo,f<=f_hi)]
    elif det == 'LISA-Babak17':
        f_lo = F_lo
        f_hi = F_hi
        check_f(det,f,f_lo,f_hi)
        return psd_LISA_Babak17(f[logical_and(f>=f_lo,f<=f_hi)]), f[logical_and(f>=f_lo,f<=f_hi)]
    elif 'LISA-Robson18' in det:
        f_lo = F_lo
        f_hi = F_hi
        check_f(det,f,f_lo,f_hi)
        curve = det.split('-')[-1]
        if curve == 'gen':
            curve = 0
        elif curve == '6mo':
            curve = 1
        elif curve == '1yr':
            curve = 2
        elif curve == '2yr':
            curve = 3
        elif curve == '4yr':
            curve = 4
        return psd_LISA_Robson18(f[logical_and(f>=f_lo,f<=f_hi)],curve), f[logical_and(f>=f_lo,f<=f_hi)]

    # use PSDs from noise_curves files
    else:
        if det == '':
            filename = psd_file
            asd = is_asd
        else:
            filename, asd = get_filename(det)
            filename = os.path.join(noise_curves_path,filename)
        psd_file = read_csv(filename, sep = " ", header = None, engine = 'python')
        psd_data = psd_file.to_numpy()
        psd = si.interp1d(psd_data[:,0], psd_data[:,1]**(1+asd))

        # find correct limits: file vs user-set limits
        f_lo = max(psd_data[0,0],F_lo)
        f_hi = min(psd_data[-1,0],F_hi)
        check_f(det,f,f_lo,f_hi)

        # return the PSD and corresponding freq array
        return psd(f[logical_and(f>=f_lo,f<=f_hi)]), f[logical_and(f>=f_lo,f<=f_hi)]

def psd_aLIGO(f):
    x = f/245.4
    return 1.e-48 * ( 0.0152 * power(x,-4.) + 0.2935 * power(x,9./4) +
                2.7951 * power(x,3./2) - 6.5080 * power(x,3./4) + 17.7622 )

def psd_CEwb(f):
    return 5.623746655206207e-51 + 6.698419551167371e-50 * power(f,-0.125) + 7.805894950092525e-31 * power(f,-20.) + 4.35400984981997e-43 * power(f,-6.) \
            + 1.630362085130558e-53 * f + 2.445543127695837e-56 * square(f) + 5.456680257125753e-66 * power(f,5)

def psd_LISA_17(f):
    a1 = 8.2047564e-33
    a2 = 3.0292821e-38
    a3 = 1.4990886e-39
    a4 = 1.0216062e-40
    return a1 * (f/10**-4)**-6 + 0.8 * a2 * (f/(10**-3))**-4 + a3 * (f/0.1)**2 + a4

def psd_LISA_Babak17(f):
    L     = 2.5e9
    SnLoc = 2.89e-24
    SnSN  = 7.92e-23
    SnOmn = 4.00e-24
    SnAcc = (9.e-30 + 3.24e-28 * ((3.e-5/f)**10 + (1.e-4/f)**2)) * (2.*pi*f)**-4
    SGal  = 3.266e-44 * f**(-7./3.) * exp(-(f/1.426e-3)**1.183) * 0.5 * (1. + tanh(-(f-2.412e-3)/4.835e-3))
    return SGal + 20./3. * (4 * SnAcc + 2 * SnLoc + SnSN + SnOmn)/L**2 * (1 + (( 2*L*f / (0.41*cLight))**2 ))

def psd_LISA_Robson18(f,curve):
    L     = 2.5e9
    fstar = 19.09e-3
    Poms  = (1.5e-11)**2 * (1. + (2.e-3/f)**4)
    Pacc  = (3.e-15)**2  * (1. + (0.4e-3/f)**2) * (1. + (f/8.e-3)**4)

    if curve == 0:
        Snoc =  10./(3. * L**2) * (Poms + 4./(2. * pi * f)**4 * Pacc)                         * (1. + 0.6 * (f/fstar)**2)
        return  10./(3. * L**2) * (Poms + 2./(2. * pi * f)**4 * (1. + cos(f/fstar)**2) * Pacc) * (1. + 0.6 * (f/fstar)**2)
    else:
        if curve == 1:
            alpha = 0.133
            beta  = 243.
            kappa = 482.
            gamma = 917.
            fk    = 0.00258
        elif curve == 2:
            alpha = 0.171
            beta  = 292.
            kappa = 1020.
            gamma = 1680.
            fk    = 0.00215
        elif curve == 3:
            alpha = 0.165
            beta  = 299.
            kappa = 611.
            gamma = 1340.
            fk    = 0.00173
        elif curve == 4:
            alpha = 0.138
            beta  = -221.
            kappa = 521.
            gamma = 1680.
            fk    = 0.00113
        Snoc =  10./(3. * L**2) * (Poms + 4./(2. * pi * f)**4 * Pacc) * (1. + 0.6 * (f/fstar)**2)
        Sc   = 9.e-45 * f**(-7./3.) * exp(-f**alpha + beta * f * sin(kappa * f)) * (1. + tanh(gamma * (fk - f)))
        return Snoc + Sc

def get_filename(det):
    if det == 'A+':
        filename = 'a_plus.txt'
        asd = 1
    elif det == 'V+':
        filename = 'advirgo_plus.txt'
        asd = 1
    elif det == 'K+':
        filename = 'kagra_plus.txt'
        asd = 1
    elif det == 'Voyager-CBO':
        filename = 'voyager_cb.txt'
        asd = 1
    elif det == 'Voyager-PMO':
        filename = 'voyager_pm.txt'
        asd = 1
    elif det == 'ET':
        filename = 'et.txt'
        asd = 1
    elif det == 'CE1-10-CBO':
        filename = 'ce1_10km_cb.txt'
        asd = 1
    elif det == 'CE1-20-CBO':
        filename = 'ce1_20km_cb.txt'
        asd = 1
    elif det == 'CE1-30-CBO':
        filename = 'ce1_30km_cb.txt'
        asd = 1
    elif det == 'CE1-40-CBO':
        filename = 'ce1_40km_cb.txt'
        asd = 1
    elif det == 'CE2-10-CBO':
        filename = 'ce2_10km_cb.txt'
        asd = 1
    elif det == 'CE2-20-CBO':
        filename = 'ce2_20km_cb.txt'
        asd = 1
    elif det == 'CE2-30-CBO':
        filename = 'ce2_30km_cb.txt'
        asd = 1
    elif det == 'CE2-40-CBO':
        filename = 'ce2_40km_cb.txt'
        asd = 1
    elif det == 'CE1-10-PMO':
        filename = 'ce1_10km_pm.txt'
        asd = 1
    elif det == 'CE1-20-PMO':
        filename = 'ce1_20km_pm.txt'
        asd = 1
    elif det == 'CE1-30-PMO':
        filename = 'ce1_30km_pm.txt'
        asd = 1
    elif det == 'CE1-40-PMO':
        filename = 'ce1_40km_pm.txt'
        asd = 1
    elif det == 'CE2-10-PMO':
        filename = 'ce2_10km_pm.txt'
        asd = 1
    elif det == 'CE2-20-PMO':
        filename = 'ce2_20km_pm.txt'
        asd = 1
    elif det == 'CE2-30-PMO':
        filename = 'ce2_30km_pm.txt'
        asd = 1
    elif det == 'CE2-40-PMO':
        filename = 'ce2_40km_pm.txt'
        asd = 1

    return filename, asd

def check_f(det,f,f_lo,f_hi):
    if f[-1] < f_lo:
        raise ValueError('The maximum frequency is below the low-frequency cutoff of the PSD for '+det,f[-1],f_lo)
    if f[0]  > f_hi:
        raise ValueError('The minimum frequency is above the high-frequency cutoff of the PSD for '+det,f[-1],f_hi)
