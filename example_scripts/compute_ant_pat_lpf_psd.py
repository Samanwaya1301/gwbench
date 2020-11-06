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
from gwbench import network

############################################################################
### User Choices
############################################################################

# include earth's rotation for Fp, Fc, lpf (location phase factor)
# (if set = 1, inj_params['Mc'] and inj_params['tc'] cannot be None)
use_rot = 0

# injection parameters (Mc, tc are only needed, if the motion of earth is relevant)
inj_params = {
    'Mc':    None,
    'tc':    None,
    'ra':    np.pi/4,
    'dec':   np.pi/4,
    'psi':   np.pi/4,
    'gmst0': 0.
    }

# choose frequency array
f_lo = 1.
f_hi = 2.**12
df = 2.**-4
f = np.arange(f_lo,f_hi+df,df)

# choose detectors/network
network_label = 'ECa4cSa4c'

############################################################################
### Initialize network and setup antenna patterns, location phase factors,
### and PSDs
############################################################################

net = network.Network(network_label)
net.set_net_vars(f=f,inj_params=inj_params,use_rot=use_rot)
net.setup_ant_pat_lpf_psds()

############################################################################
### Get antenna patterns, location phase factors, and PSDs
############################################################################

print()

# Method 1 - using the detector keys to get specific detectors
print('Method 1')
print(net.det_keys)
print()
for det_key in net.det_keys:
    print(det_key)
    det = net.get_detector(det_key)
    print(det.f)
    print(det.Fp)
    print(det.Fc)
    print(det.Flp)
    print(det.psd)
    print()

print()

# Method 2 - cycling through the detectors directly
print('Method 2')
print(net.detectors)
print()
for det in net.detectors:
    print(det.det_key)
    print(det.f)
    print(det.Fp)
    print(det.Fc)
    print(det.Flp)
    print(det.psd)
    print()
