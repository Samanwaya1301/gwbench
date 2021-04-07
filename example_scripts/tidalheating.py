#copyright (C) 2020  Ssohrab Borhanian
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
import gwbench.basic_relations as brs
from gwbench import network

############################################################################
### User Choices
############################################################################

# choose the desired detectors
network_spec = ['CE2-40-CBO_C', 'CE2-40-CBO_N', 'CE2-40-CBO_S']
# initialize the network with the desired detectors
net = network.Network(network_spec)

# choose the desired waveform 
wf_model_name = 'heated_tf2'
# pass the chosen waveform to the network for initialization
net.set_wf_vars(wf_model_name=wf_model_name)

# pick the desired frequency range
#f = np.arange(5.,67.6,2**-4)

#define arrays for errors
err_Log_Mc = []
err_H_eff5 = []
err_H_eff8 = []
chi_val =[]
#input m1, m2 instead of Mc eta and then convert
m1 = 25.
m2 = 2.
chi = 0.
Mc,eta = brs.Mc_eta_of_m1_m2(m1,m2)
#loop for plotting
j = 1

while j>0:
	# set the injection parameters
	inj_params = {
	    'Mc':    Mc,
	    'eta':   eta,
	    'chi1z': chi,
	    'chi2z': chi,
	    'DL':    200,
	    'tc':    0,
	    'phic':  0,
	    'iota':  np.pi/4,
	    'ra':    np.pi/4,
	    'dec':   np.pi/4,
	    'psi':   np.pi/4,
	    'gmst0': 0,
	    'Heff5': 1,
	    'Heff8': 15,
	    'e0': 0.1
	    }
	#calculating isco frequency

	#Mc = inj_params["Mc"]
	#eta = inj_params["eta"]

	M = brs.M_of_Mc_eta(Mc,eta)
	f_isco = brs.f_isco_Msolar(M)


	#check the desired frequency range
	f = np.arange(4.,f_isco,2**-4)



	#wf_other_var_dic = {'Heff5': 0, 'Heff8': 0}
	# assign with respect to which parameters to take derivatives
	deriv_symbs_string = 'Mc Heff5 Heff8'

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
	#net.print_detectors()

	# print the contents of the network objects
	#net.print_network()
	#print(net.errs)
	#err_Log_Mc.append(net.errs["log_Mc"])
	err_H_eff5.append(net.errs["Heff5"])
	err_H_eff8.append(net.errs["Heff8"])
	chi_val.append(chi)
	
	#print injection values of masses and spins
	chi1,chi2 = inj_params["chi1z"],inj_params["chi2z"]
	#m1, m2 = brs.m1_m2_of_Mc_eta(Mc,eta)
	#print("m1 = {}, m2 = {}, chi1 = {}, chi2 = {} ".format(m1, m2, chi1, chi2))
	

	chi += 0.1
	if chi>=0.95: break


#export outputs to a file
#----------------------------

with open ("error_m1_{}_m2_{}.txt".format(m1, m2), "w") as f:
	f.write("!chi\terr_H_eff5\terr_Heff8\n")
	for i in range(0,len(chi_val)):
		f.write("{0}\t{1}\t{2}\n".format(chi_val[i], err_H_eff5[i], err_H_eff8[i]))
