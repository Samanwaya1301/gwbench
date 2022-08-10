import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

figure, axs = plt.subplots(2,2)

chi1 = np.linspace(0,1,5)
chi2 = np.linspace(0,1,5)

[x,y] = np.meshgrid(chi1,chi2)

for det in ['ce','et']:
	err_H_eff5 = np.loadtxt("/home/samanwaya/gwbench_data/nondeg/iscoKBH/incl_all/spin-data/spin-h5-{}.txt".format(det))
	err_H_eff8 = np.loadtxt("/home/samanwaya/gwbench_data/nondeg/iscoKBH/incl_all/spin-data/spin-h8-{}.txt".format(det))
	err_H_eff52 = np.loadtxt("/home/samanwaya/gwbench_data/nondeg/iscoKBH/incl_all/spin-data/spin-h5-1-{}.txt".format(det))
	err_H_eff82 = np.loadtxt("/home/samanwaya/gwbench_data/nondeg/iscoKBH/incl_all/spin-data/spin-h8-1-{}.txt".format(det))


	if det=='et':
		axs0et = axs[1,0].contour(x,y,err_H_eff5,4,colors='k')
		axs1et = axs[1,1].contour(x,y,err_H_eff8,4,colors='k')
		axs[1,0].contourf(x,y,err_H_eff5,4,cmap='viridis',alpha=0.5)
		axs[1,1].contourf(x,y,err_H_eff8,4,cmap='viridis',alpha=0.5)	
		axs[1,0].clabel(axs0et,inline=1, fmt='%1.2f',fontsize=8)
		axs[1,1].clabel(axs1et,inline=1, fmt='%1.1f',fontsize=8)

		axs0et2 = axs[0,0].contour(x,y,err_H_eff52,4,colors='k')
		axs1et2 = axs[0,1].contour(x,y,err_H_eff82,4,colors='k')
		axs[0,0].contourf(x,y,err_H_eff52,4,cmap='viridis',alpha=0.5)
		axs[0,1].contourf(x,y,err_H_eff82,4,cmap='viridis',alpha=0.5)	
		axs[0,0].clabel(axs0et2,inline=1, fmt='%1.2f',fontsize=8)
		axs[0,1].clabel(axs1et2,inline=1, fmt='%1.1f',fontsize=8)

	if det=='ce':
		axs0ce = axs[1,0].contour(x,y,err_H_eff5,4,linestyles='dashed',colors='b')
		axs1ce = axs[1,1].contour(x,y,err_H_eff8,4,linestyles='dashed',colors='b')	
		axs[1,0].clabel(axs0ce,inline=1, fmt='%1.2f',fontsize=8)
		manual_loc2 = [(0.25,0.2),(0.85,0.2)]
		manual_loc = [(0.05,0.6),(0.2,0.2),(0.4,0.2),(0.6,0.2),(0.8,0.2)]
		manual_loc3 = [(0.1,0.5),(0.2,0.2),(0.6,0.3)]
		axs[1,1].clabel(axs1ce,inline=1, manual=manual_loc, fmt='%1.1f',fontsize=8)
		#axs[1,1].clabel(axs1ce,inline=1, fmt='%1.2f',fontsize=8)

		axs0ce2 = axs[0,0].contour(x,y,err_H_eff52,4,linestyles='dashed',colors='b')
		
		axs1ce2 = axs[0,1].contour(x,y,err_H_eff82,4,linestyles='dashed',colors='b')	
		axs[0,0].clabel(axs0ce2,inline=1, manual=manual_loc2, fmt='%1.2f',fontsize=8)
		#manual_loc = [(0.05,0.6),(0.2,0.2),(0.4,0.2),(0.6,0.2),(0.8,0.2),(1,0.2)]
		#axs[1,1].clabel(axs1ligo,inline=1, manual=manual_loc, inline_spacing=10, fmt='%1.1f',fontsize=8)
		axs[0,1].clabel(axs1ce2,inline=1, manual=manual_loc3, fmt='%1.1f',fontsize=8)


labels=['ET','CE']

custom_lines = [Line2D([0], [0], color='k', lw=1.5),
                Line2D([0], [0], color='b', lw=1.5, ls='--')]

xticks = [0,0.2,0.4,0.6,0.8,1]
yticks = [0,0.2,0.4,0.6,0.8,1]

axs[0,0].set_xlabel("$\\chi_1$",fontsize=14)
axs[0,0].set_title("$\\Delta H_{\\rm eff5}$",fontsize=14)
axs[0,1].set_xlabel("$\\chi_1$",fontsize=14)
axs[0,0].set_ylabel("$\\chi_2$",fontsize=14)
#axs[0,1].set_ylabel("$\\chi_2$",fontsize=14)

axs[1,0].set_xlabel("$\\chi_1$",fontsize=14)
axs[1,1].set_xlabel("$\\chi_1$",fontsize=14)
axs[0,1].set_title("$\\Delta H_{\\rm eff8}$",fontsize=14)
axs[1,0].set_ylabel("$\\chi_2$",fontsize=14)
#axs[1,1].set_ylabel("$\\chi_2$",fontsize=14)

for i1 in [0,1]:
	for j1 in [0,1]:
		axs[i1,j1].set_xticklabels(xticks)
		axs[i1,j1].set_yticklabels(yticks) 
		#axs[i1,j1].set_box_aspect(1)


figure.legend(custom_lines,labels,fontsize=10)

#for ax in axs.flat:
#	ax.label_outer()

#axs[1].legend(custom_lines,labels,fontsize=12)
'''
axs[0,0].set_box_aspect(1)
axs[0,1].set_box_aspect(1)
axs[1,0].set_box_aspect(1)
axs[1,1].set_box_aspect(1)
'''
#figure.tight_layout()
#plt.show()
plt.savefig("spin-40-1-3.png",dpi=300,bbox_inches='tight')


