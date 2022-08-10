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
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

############################################################################
### User Choices
############################################################################

# choose the desired detectors
#network_spec = ['aLIGO_H', 'aLIGO_L', 'aLIGO_V']
#network_spec = ['ET_ET1', 'ET_ET2', 'ET_ET3']
#network_spec = ['CE2-40-CBO_C', 'CE2-40-CBO_N', 'CE2-40-CBO_S']
# initialize the network with the desired detectors

# pick the desired frequency range
#f = np.arange(5.,67.6,2**-4)

#define arrays for errors
err_Log_Mc = []

#input m1, m2 instead of Mc eta and then convert
M = 40.
q = 1.1

#Mc,eta = brs.Mc_eta_of_m1_m2(m1,m2)
#loop for plotting
#figure, axs = plt.subplots(2,1)
	
for det in ['ligo','et','ce']:     #  detector loop

	if det=='ligo': network_spec = ['aLIGO_H', 'aLIGO_L', 'aLIGO_V']
	if det=='et': network_spec = ['ET_ET1', 'ET_ET2', 'ET_ET3']
	if det=='ce': network_spec = ['CE2-40-CBO_C', 'CE2-40-CBO_N', 'CE2-40-CBO_S']

	net = network.Network(network_spec)

    # choose the desired waveform 
	wf_model_name = 'heated_tf2_nondeg'
	# pass the chosen waveform to the network for initialization
	net.set_wf_vars(wf_model_name=wf_model_name)

	chi1 = np.linspace(0,1,5)
	chi2 = np.linspace(0,1,5)

	err_H_eff5 = np.zeros((len(chi1),len(chi2)))
	err_H_eff8 = np.zeros((len(chi1),len(chi2)))
	

	#[x,y] = np.meshgrid(chi1,chi2) 

	for i in range(len(chi1)):

		for j in range(len(chi2)):

			# set the injection parameters
			
			m1 = (q/(q+1))*M
			m2 = (1/(q+1))*M
			Mc,eta = brs.Mc_eta_of_m1_m2(m1,m2)
			inj_params = {
			    'Mc':    Mc,
			    'eta':   eta,
			    'chi1z': chi1[i],
			    'chi2z': chi2[j],
			    'DL':    200,
			    'tc':    0,
			    'phic':  0,
			    'iota':  np.pi/4,
			    'ra':    np.pi/4,
			    'dec':   np.pi/4,
			    'psi':   np.pi/4,
			    'gmst0': 0,
			    'Heff5': 0.6,
			    'Heff8': 12.
			    }
			#calculating isco frequency

			#Mc = inj_params["Mc"]
			#eta = inj_params["eta"]

			#M = brs.M_of_Mc_eta(Mc,eta)
			f_isco = brs.f_isco_Msolar_KBH(m1,m2,chi1[j],chi2[j])


			#check the desired frequency range
			if det=='ligo': f = np.arange(10.,f_isco,2**-4)
			else: f = np.arange(4.,f_isco,2**-4)



			#wf_other_var_dic = {'Heff5': 0, 'Heff8': 0}
			# assign with respect to which parameters to take derivatives
			deriv_symbs_string = 'Mc eta DL chi1z chi2z Heff5 Heff8'

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
			net.calc_snrs()

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
			err_H_eff5[i,j] = net.errs["Heff5"]
			err_H_eff8[i,j] = net.errs["Heff8"]
			#chi_val.append(chi)
			
			#print injection values of masses and spins
			#chi1,chi2 = inj_params["chi1z"],inj_params["chi2z"]
			#m1, m2 = brs.m1_m2_of_Mc_eta(Mc,eta)
			#print("m1 = {}, m2 = {}, chi1 = {}, chi2 = {} ".format(m1, m2, chi1, chi2))
	

	np.savetxt("/home/samanwaya/gwbench_data/nondeg/iscoKBH/incl_all/spin-data/spin-h5-1-{}.txt".format(det),err_H_eff5)
	np.savetxt("/home/samanwaya/gwbench_data/nondeg/iscoKBH/incl_all/spin-data/spin-h8-1-{}.txt".format(det),err_H_eff8)


'''

	if det == 'ligo':	
		axs0 = axs[0].contourf(x,y,err_H_eff5)
		axs1 = axs[1].contourf(x,y,err_H_eff8)
		axs0 = axs[0].contour(x,y,err_H_eff5)
		axs1 = axs[1].contour(x,y,err_H_eff8)
	if det == 'et':
		axs0 = axs[0].contour(x,y,err_H_eff5, linestyles='dashed')
		axs1 = axs[1].contour(x,y,err_H_eff8, linestyles='dashed')
	
	axs[0].clabel(axs0,inline=1,fontsize=10)
	axs[1].clabel(axs1,inline=1,fontsize=10)

axs[0].set_xlabel("$\\chi_1$")
axs[0].set_ylabel("$\\chi_2$")
axs[1].set_xlabel("$\\chi_1$")
axs[1].set_ylabel("$\\chi_2$")

plt.show()
'''

