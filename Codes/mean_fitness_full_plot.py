#****************************************************************************
# simNK version 0.1
# Project AevolStatPhys
#****************************************************************************
# mean_max_fitness.py / version 0.1 
#
# Program to calculate the mean maximum fitness values in function of the
# epistatic interaction parameter K, parametrized by the mutations fractions
# M.
#
# Created:  2021-01-29 
# Modified: 2021-02-01
# by Leonardo Trujillo (leonardo.trujillo@gmail.fr)
#
# You may use, share, or modify this file freely
#
#****************************************************************************
# PROGRAM PARAMETERS (Global variables):
#
# N:            The length of the genomes.
# K:            The number of epistatic interactions.
# epi:          The epistatic interactions: ADJ (adjacent) or RND (random).
# M:            Fraction of inversions:
#                   m = 0 <- M = 1.0: inversions
#                   m = 1 <- M = 0.0: point mutationd
#                   m = 3 <- M = 0.5: point mutations and inversions 
		     
# k_min:        Minimum value of the interval of K. 
# k_max:        Maximun value of the interval of K.
# delta_k:      The increments of K.
# realizations: The total number of instances simulated for a given value of
#               the epistatic interaction parameter K
#
#****************************************************************************
# IMPORTAN: The raw data comes from the simulations performed with "nk_walk" 
# saved in the file: "final_fitness.csv".
#
# Before to run the present code, please rename the "final_fitness.csv" file
# for example as: "final_fitness_N50_RND.csv".
# This code assumes that the number of experiments that have been made on 
# each mutation types is equal to realization.
#
#****************************************************************************
from re import A
from statistics import mean
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

N = 100
M = 2
k_min = 0
k_max = N-1
delta_k = 1

realizations = 100 

#  Modify here if you want to use different dots
dot1 = "+"
dot2 = "o"

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

def plotepi(ax, epi,dot):
	'''
	Function plotting the main figure 
	'''
	#-----------------------------------------------------------------------------
	# Arrays initializations
	#-----------------------------------------------------------------------------
	fit = np.zeros((M,N))

	#----------------------------------------------------------------------------
	# Read the data from final_fitness_* files 
	# whose format is: df[0] = N, df[1] = K, df[2] = M, df[3] = max_fitness 
	#-----------------------------------------------------------------------------

	file_name = ("final_fitness_N%s_%s.csv" % (N, epi))
	df = pd.read_csv(file_name, sep = ';',header = None).values
	total_samples = df.shape[0]
	
	#-----------------------------------------------------------------------------
	# Averages 
	#-----------------------------------------------------------------------------
	for i in range(len(df)):
		k = int(df[i,1])
		m = df[i,2]
		m_int = int(m)
		F = df[i,3]
		if m == m_int:
			fit[m_int,k] += F/(realizations)

	#-----------------------------------------------------------------------------
	# Plot 
	#-----------------------------------------------------------------------------
	x1 = np.arange(0, k_max + 1, delta_k)

	for i in range(M):
		if i == 0:
			marker = 'b'+dot
		else:
			marker = 'r'+dot
		for j in range(N):
			#plt.plot(x1[j], fit[i,j], marker, markersize = 2)
			ax.plot(x1[j], fit[i,j], marker, markersize = 4,mfc='none', label = "") #uncomment here
		L = [fit[i,k] for k in range(N)]
		print(marker , "max F(k): " , max(L),"argmax F(k)" ,np.argmax(L))


def plot_all_epi(ax):
	'''
	Function used to define the setup of the main plot, legends, grid, labels , ...
	'''
	red_patch = mpatches.Patch(color='red', label='Inversions')
	blue_patch = mpatches.Patch(color='blue', label='Point mutations')
	purple_patch = mpatches.Patch(color='purple', label=r'$\Delta f_\mathrm{K} = \langle f_I \rangle _{\mathrm{K}} - \langle f_P \rangle _{\mathrm{K}} $')
	adjmark, = ax.plot([], [], 'k'+dot1, markersize = 4,mfc='none', label = "adjacent")
	rndmark, = ax.plot([], [], 'k'+dot2, markersize = 4,mfc='none', label = "random")
	plotepi(ax, "ADJ",dot1)
	plotepi(ax, "RND",dot2)
	ax.legend(handles=[red_patch,blue_patch,purple_patch,adjmark,rndmark],loc = 'lower left')
	ax.grid()
	ax.set_xlabel(r'K', fontsize = 20)
	ax.set_ylabel(r'$\langle f \rangle _{\mathrm{K}} $', fontsize = 20)
	ax.tick_params(labelsize = 12)


fig_file_name = ("fig_final_fitness_N%s_all_epi.png" % (N))


def plotdiff(ax):
	'''
	Function used to plot the difference between 
	'''
	#-----------------------------------------------------------------------------
	# Arrays initializations
	#-----------------------------------------------------------------------------
	diff_fit_adj = [[] for k in range(N)]
	diff_fit_rnd = [[] for k in range(N)]
	means = [[],[]]
	lis_diff = [diff_fit_adj,diff_fit_rnd]
	epi_list = ["ADJ","RND"]
	#----------------------------------------------------------------------------
	# Read the data from final_fitness_* files 
	# whose format is: df[0] = N, df[1] = K, df[2] = M, df[3] = max_fitness 
	#-----------------------------------------------------------------------------
	for epi_ind in range(2):
		epi = epi_list[epi_ind]
		
		file_name = ("final_fitness_N%s_%s.csv" % (N, epi))
		df = pd.read_csv(file_name, sep = ';',header = None).values
		total_samples = df.shape[0]
		#-----------------------------------------------------------------------------
		# Averages 
		#-----------------------------------------------------------------------------
		for i in range(0,len(df),2): #here we assume the data is in shape of inversion then switch on the same env
			k = int(df[i,1])
			F_inv = df[i,3]
			k2 = int(df[i+1,1])
			F_switch = df[i+1,3]
			assert(k == k2)
			lis_diff[epi_ind][k].append(F_inv - F_switch) # here we compute the difference ! 
		means[epi_ind] = [np.mean(lis_diff[epi_ind][k]) for k in range(N)]
		curL = means[epi_ind][0:21]
		print( epi, " K <= 20 :      ", 'max(Delta f_K): ' , max(curL),"; argmax(Delta f_K): " ,np.argmax(curL))
		lowerval = 21
		curL = means[epi_ind][lowerval:100]
		print( epi, " " + str(lowerval-1) +"< K < 100 : ",  'min(Delta f_K): ' , min(curL),"; argmin(Delta f_K): " , lowerval + np.argmin(curL))



	
	#-----------------------------------------------------------------------------
	# Plot 
	#-----------------------------------------------------------------------------
	x1 = np.arange(0, k_max + 1, delta_k)
	

	for epi_ind in range(2):
		for k in range(N):
			for val in [means[epi_ind][k]]: #lis_diff[epi_ind][k]: #
				if epi_ind == 0 :
					marker = dot1
				else:
					marker = dot2
				#plt.plot(x1[j], fit[i,j], marker, markersize = 2)
				#ax.plot(x1[k], val , marker, markersize = 4 ,alpha = 0.01) 
				ax.plot(x1[k], val , color = 'purple', marker =marker, markersize = 4,mfc='none' ,alpha = 1) 
	ax.axhline(y= 0,color = 'orange',linewidth = 1)
	

def plot_ins(ax):
	'''
	Function used to setup the inset, mainly axes labels ...
	'''
	plotdiff(ax)
	ax.set_xlabel(r'K', fontsize = 12)
	ax.set_ylabel(r'$\Delta f_{\mathrm{K}}$', fontsize = 12)
	ax.tick_params(labelsize = 7)



# Finally after all function have been defined, the next part applies them and save the figure
fig, ax1 = plt.subplots()

# These are in unitless percentages of the figure size. (0,0 is bottom left)
left, bottom, width, height = [0.54, 0.51, 0.35, 0.35]
ax2 = fig.add_axes([left, bottom, width, height])

plot_all_epi(ax1)
plot_ins(ax2)
#plt.show()

fig_file_name = ("fig_final_fitness_N%s_all_epi_inset.png" % (N))
plt.savefig(fig_file_name, dpi=300, bbox_inches='tight')

    
#-----------------------------------------------------------------------------
# Bye!  
#-----------------------------------------------------------------------------
