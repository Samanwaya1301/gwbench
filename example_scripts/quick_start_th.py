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

# choose the desired detectors
network_spec = ['aLIGO_H','aLIGO_L','aLIGO_V']
# initialize the network with the desired detectors
net = network.Network(network_spec)

# choose the desired waveform 
wf_model_name = 'heated_tf2'
# pass the chosen waveform to the network for initialization
net.set_wf_vars(wf_model_name=wf_model_name)

# pick the desired frequency range
f = np.arange(5.,67.6,2**-4)

# set the injection parameters
inj_params = {
    'Mc':    28.1,
    'eta':   0.247,
    'chi1z': 0,
    'chi2z': 0,
    'DL':    440,
    'tc':    0,
    'phic':  0,
    'iota':  np.pi/4,
    'Heff5': 0,
    'Heff8': 0,
    'ra':    np.pi/4,
    'dec':   np.pi/4,
    'psi':   np.pi/4,
    'gmst0': 0
    }

# assign with respect to which parameters to take derivatives
deriv_symbs_string = 'Mc eta DL tc phic iota ra dec psi'

# assign which parameters to convert to cos or log versions
conv_cos = ('iota','dec')
conv_log = ('Mc','DL')

# choose whether to take Earth's rotation into account
use_rot = 0

# pass all these variables to the network
net.set_net_vars(
    f=f, inj_params=inj_params,
    deriv_symbs_string=deriv_symbs_string,
    conv_cos=conv_cos, conv_log=conv_log,
    use_rot=use_rot
    )

############################################################################
### GW benchmarking
############################################################################

# compute the WF polarizations
net.calc_wf_polarizations()
# compute the WF polarizations and their derivatives
net.calc_wf_polarizations_derivs_num()

# setup antenna patterns, location phase factors, and PSDs
net.setup_ant_pat_lpf_psds()

# compute the detector responses
net.calc_det_responses()
# compute the detector responses and their derivatives
net.calc_det_responses_derivs_num()

# calculate the network and detector SNRs
net.calc_snrs_det_responses()

# calculate the network and detector Fisher matrices, condition numbers,
# covariance matrices, error estimates, and inversion errors
net.calc_errors()

# calculate the 90%-credible sky area (in deg)
net.calc_sky_area_90()

############################################################################
### Print results
############################################################################

# print the contents of the detector objects (inside the network)
net.print_detectors()

# print the contents of the network objects
net.print_network()
