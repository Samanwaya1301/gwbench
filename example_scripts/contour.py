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
import gwbench.basic_relations as brs
from gwbench import network
import matplotlib.pyplot as plt
import scipy.linalg as la

############################################################################
### User Choices
############################################################################

# choose the desired detectors
network_spec = ['aLIGO_H', 'aLIGO_L', 'aLIGO_V']
#network_spec = ['CE2-40-CBO_C', 'CE2-40-CBO_N', 'CE2-40-CBO_S']
# initialize the network with the desired detectors
net = network.Network(network_spec)

# choose the desired waveform 
wf_model_name = 'heated_tf2_nondeg'
# pass the chosen waveform to the network for initialization
net.set_wf_vars(wf_model_name=wf_model_name)


mass_val = [50,60,70,80]
q_val = [4,6]
Heff5 = 0.6
Heff8 = 12.

figure, axis = plt.subplots(len(q_val),len(mass_val))





for i in range(0,len(q_val)):

    for j in range(0,len(mass_val)):

        M = mass_val[j]
        q = q_val[i]
        chi = 0.8
        m1 = (q/(q+1))*M
        m2 = (1/(q+1))*M
        Mc,eta = brs.Mc_eta_of_m1_m2(m1,m2)

        # set the injection parameters
        inj_params = {
            'Mc':    Mc,
            'eta':   eta,
            'chi1z': chi,
            'chi2z': chi,
            'DL':    500,
            'tc':    0,
            'phic':  0,
            'iota':  np.pi/4,
            'ra':    np.pi/4,
            'dec':   np.pi/4,
            'psi':   np.pi/4,
            'gmst0': 0,
            'Heff5': Heff5,
            'Heff8': Heff8
            }


        f_isco = brs.f_isco_Msolar(M)
        f = np.arange(10.,f_isco,2**-4)

        # assign with respect to which parameters to take derivatives
        deriv_symbs_string = 'Heff5 Heff8'

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

        
        fish = net.fisher

        snr = net.snr 

        fish = (1/snr**2)*fish




        eigvals, eigvecs = la.eig(fish)
        eigvals = eigvals.real
        eigvec1 = eigvecs[:,0]
        eigvec2 = eigvecs[:,1]
        print("eigenvalues = ", eigvals)




        #centre of the ellipses (values of Heff5 & Heff8)
        a = Heff5
        b = Heff8
        

        t=np.linspace(0.,360.,360)  # angle parameter for plotting the ellipse
        l = eigvec1[0]
        m = eigvec1[1]
        axs = [1,0]   # x-axis
        vec1 = [l,m]
        vec1 = vec1/np.linalg.norm(vec1)   # normalizing vec1
        vec2 = [-m,l]
        vec2 = vec2/np.linalg.norm(vec2)   # normalizing vec2
        dot = np.dot(axs,vec1)

        theta = np.arctan(m/l)   # angle between vec1 and the x-axis

        lambda1 = np.arccos(np.dot(vec1,vec2))*180/(np.pi)
        theta = theta*180./(np.pi)

        print("theta in degrees = ",theta)
        print("lambda1=",lambda1)
        print("vec1=",vec1)
        print("vec2=",vec2)

        Match = [0.7,0.8,0.9,0.99]
        for k in range(0,len(Match)):  
          #eigenvectors
          l1 = eigvals[0]
          l2 = eigvals[1]

          p1 = np.sqrt((1-Match[k])/l1)
          q1 = np.sqrt((1-Match[k])/l2)

          x = a + p1*np.cos(np.radians(theta))*np.cos(np.radians(t)) - q1*np.sin(np.radians(theta))*np.sin(np.radians(t))
          y = b + p1*np.sin(np.radians(theta))*np.cos(np.radians(t)) + q1*np.cos(np.radians(theta))*np.sin(np.radians(t))

          if i==j==0:
            axis[i,j].plot(x,y,label="{}".format(Match[k]))
          else:
            axis[i,j].plot(x,y)

          

        
        #r=max(p,q)
        axis[i,j].axvline(a,ls='--')
        #plt.axvline(0,ls='--',color='red')
        axis[i,j].axhline(b,ls='--')
        #axis[i,j].axhline(0,ls='--',color='red')
        axis[i,j].grid(color='lightgrey')
        

        w=np.linspace(-500,500,5000)
        z1=(m/l)*(w-a)+b
        z2=(-l/m)*(w-a)+b
        axis[i,j].plot(w,z1,color='black',ls='--')
        axis[i,j].plot(w,z2,color='black',ls='--')
        axis[i,j].set_xlim(-200,200)
        axis[i,j].set_ylim(-200,200)
        axis[i,j].set_title("M={}, q={}".format(M,q))
        #axis[i,j].legend(loc="best")

figure.legend()
for ax in axis.flat:
    ax.set(xlabel='$H_{eff5}$', ylabel='$H_{eff8}$')
    ax.xaxis.label.set_fontsize(18)
    ax.yaxis.label.set_fontsize(18)

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axis.flat:
    ax.label_outer()
plt.show()
