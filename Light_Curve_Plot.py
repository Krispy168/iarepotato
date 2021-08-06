#===================================Pramble====================================#
#Program Description: Plots light curves from list of ids hardcoded in.

#Inputs: Has no inputs to terminal

#Packages: pip install numpy math pdb matplotlib 

#Run: python3 Light_Curve_Plot.py

#Outputs: Saves plots to current directory
#==============================Importing Packages==============================#

import numpy as np
import math
import pdb
import matplotlib.pyplot as plt

#=======================Opening and Reading Catalogue File=====================#

#List of light curves to plot 
ids = ['19884','3126174']

#ids = ["26925", "254251", "282281", "320730", "379582", "1190461", "1473370", "1534907", "2164845", "2305627", "2429787", "2443925", "2857463", "3059882", "3098344", "3126174", "5841454", "6225854", "6226243", "6269796", "6281366", "6294863", "6370555", "6408239", "6460480", "6606804", "6621138", "6639067", "6686288", "6745779", "6865368", "6869137", "6905843", "6926287", "7040060", "7089712", "7317895", "7392340", "7434004", "7444196", "7537266", "7538525", "7562741", "7585442"]  #0<Chi2<1.2

#ids = ["7838", "147204", "228928", "234371", "320730", "347307", "379582", "435346", "473666", "541411", "562018", "647108", "666954", "675083", "701924", "719653", "734665", "756245", "924524", "967088", "980799", "1008531", "1044020", "1116816", "1121819", "1129529", "1129597", "1157897", "1593720", "1633742", "1675922", "1716258", "1731818", "2179775", "2184332", "2192435", "2252824", "2337091", "2362402", "2397698", "2450815", "2476069", "2815198", "2861802", "2998000", "2998960", "3118237", "3120615", "3138327", "3642034", "3813054", "3894418", "3900618", "3969724", "3994311", "4925629", "5382104", "5645121", "6231187", "7538525"]  #1<Chi2<1.5

#================================Function Definition===========================#

def plotfunc(x):

	#------Dictionary for loop use------
	Filters = {'r':{'Time':[],'Mag':[],'Err':[] },
	 'i':{'Time':[],'Mag':[],'Err':[] },
	 'u':{'Time':[],'Mag':[],'Err':[] },
	 'z':{'Time':[],'Mag':[],'Err':[] },
	 'g':{'Time':[],'Mag':[],'Err':[] }}

    #Clearing data variable
	data = None
	with open("LC_files/LC_" + item + ".dat", 'r') as f:
		data = f.read()
		data = data.splitlines()

	#------Splitting Data and Assigning to Dictionaries------
	#Filters Dictionary
	for i in data:
		x = i.split()
		if x[0][0] == "#":
		    pass
		else:
			y = Filters.get(x[1])
			y['Time'].append(float(x[0]))
			y['Mag'].append(float(x[2]))
			y['Err'].append(float(x[3]))

	#------Plotting light curves------

	#Plotting Variables
	rmag = np.array(Filters['r']['Mag'])
	rtime = np.array(Filters['r']['Time'])
	rerr = np.array(Filters['r']['Err'])

	imag = np.array(Filters['i']['Mag'])
	itime = np.array(Filters['i']['Time'])
	ierr = np.array(Filters['i']['Err'])

	umag = np.array(Filters['u']['Mag'])
	utime = np.array(Filters['u']['Time'])
	uerr = np.array(Filters['u']['Err'])

	zmag = np.array(Filters['z']['Mag'])
	ztime = np.array(Filters['z']['Time'])
	zerr = np.array(Filters['z']['Err'])

	gmag = np.array(Filters['g']['Mag'])
	gtime = np.array(Filters['g']['Time'])
	gerr = np.array(Filters['g']['Err'])

	#Scatter Plot
	plt.scatter(rtime,rmag,label='r',color='r')
	plt.scatter(itime,imag,label='i',color='m')
	plt.scatter(utime,umag,label='u',color='blue')
	plt.scatter(ztime,zmag,label='z',color='c')
	plt.scatter(gtime,gmag,label='g',color='green')
	#Plot legend
	plt.legend(scatterpoints=1,loc=2)
	
	'''
	#Plotting error bars
	plt.errorbar(rtime,rmag,yerr=rerr, linestyle="None",color='r')
	plt.errorbar(itime,imag,yerr=ierr, linestyle="None",color='m')
	plt.errorbar(utime,umag,yerr=uerr, linestyle="None",color='blue')
	plt.errorbar(ztime,zmag,yerr=zerr, linestyle="None",color='c')
	plt.errorbar(gtime,gmag,yerr=gerr, linestyle="None",color='green')
    '''

	#Plot labeling
	plt.title('Light Curve for LC_{0}'.format(item))
	plt.xlabel('Time(Julian Days)')
	plt.ylabel('Uncorrected ugriz Magintude')

	#plt.show()
	#Save figure and relavent variables(file) to folder.
	plt.savefig('LC_{0}.png'.format(item)) #Save figure to folder  
	plt.close()

#============================Executing Plot Function===========================#
	
for item in ids: 
	plotfunc(item)
	print(item)
	
#==============================================================================#
#=================================End of File==================================#
#==============================================================================#
