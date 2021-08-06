#===================================Pramble====================================#
#Program Description: This program is used in the initial step of processing 
#   the light curve data. It computes the difference between the filters, for
#   all cases (10 in total). It then computes the Chi-Squared minimum value for 
#   for each filter difference. These values along with the accumulative value
#   of Chi-Squared for each light curve is then saved as output to be used later. 

#Inputs: No inputs to terminal.

#Packages: pip install numpy math pdb matplotlib time pickle

#Run: python3 Major_Project.py

#Outputs: 1) For each object it outputs the id when analysis is complete
#         2) Will save a file 'chi2min.pickle' with results 
#         3) Prints computation time to terminal
#==============================Importing Packages==============================#

import numpy as np
import math
import pdb
import matplotlib.pyplot as plt
import time
import pickle

#=======================Opening and Reading Catalogue File=====================#

f0 = open('stripe82.dat','r')
stripe82 = f0.read()
stripe82 = stripe82.splitlines()

#ID Dictionary 
Obid = []
period = []

#Populating Obid and Array
for i in stripe82:
    x = i.split()
    if x[0][0] != "#":
        per = float(x[3])
        if per > 2:
            Obid.append(x[0])
            period.append(float(x[3]))
    
f0.close()

#===============================Chi-Squared Arrays=============================#

#Chi-Squared min array
chisqmin = []

#==============================Function Definition=============================#

#Function reads in and calculates chi-squared distribution for each light curve
def chi2(item):
    
    #------Dictionaries for loop use------
    Filters = {'r':{'Time':[],'Mag':[],'Err':[] },
               'i':{'Time':[],'Mag':[],'Err':[] },
               'u':{'Time':[],'Mag':[],'Err':[] },
               'z':{'Time':[],'Mag':[],'Err':[] },
               'g':{'Time':[],'Mag':[],'Err':[] }}

    Sigma = {'ug':{'ugtime':[],'ugmag':[],'ugerr':[],'ugchi2':[]},
             'gr':{'grtime':[],'grmag':[],'grerr':[],'grchi2':[]},
             'ri':{'ritime':[],'rimag':[],'rierr':[],'richi2':[]},
             'iz':{'iztime':[],'izmag':[],'izerr':[],'izchi2':[]},
             'zu':{'zutime':[],'zumag':[],'zuerr':[],'zuchi2':[]},
             'ur':{'urtime':[],'urmag':[],'urerr':[],'urchi2':[]},
             'ui':{'uitime':[],'uimag':[],'uierr':[],'uichi2':[]},
             'gi':{'gitime':[],'gimag':[],'gierr':[],'gichi2':[]},
             'gz':{'gztime':[],'gzmag':[],'gzerr':[],'gzchi2':[]},
             'rz':{'rztime':[],'rzmag':[],'rzerr':[],'rzchi2':[]}}

    data = None
    with open("LC_files/LC_" + item + ".dat", 'r') as f:
	    data = f.read()
	    data = data.splitlines()

    #------Splitting Data and Assigning to Dictionaries------

    #Filters Dictionary
    for i in data:
	    x = i.split()
	    if x[0][0] == "#":
		    y = Filters.get(x[1])
		    y['Time'].append(float(x[0][1:]))
		    y['Mag'].append(-99.00) #Changed to None
		    y['Err'].append(-99.00)
	    else:
		    y = Filters.get(x[1])
		    y['Time'].append(float(x[0]))
		    y['Mag'].append(float(x[2]))
		    y['Err'].append(float(x[3]))
    
    #Sigma Dictionary
    #ug filter
    for i in range(len(Filters['g']['Time'])):
	    if Filters['u']['Mag'][i] != -99.00 or Filters['g']['Mag'][i] != -99.00:	
		    Sigma['ug']['ugmag'].append(abs(Filters['u']['Mag'][i]) - abs(Filters['g']['Mag'][i]))	
		    Sigma['ug']['ugtime'].append(float((Filters['u']['Time'][i] + Filters['g']['Time'][i])/2))
		    Sigma['ug']['ugerr'].append(math.sqrt(Filters['u']['Err'][i]**2 + Filters['g']['Err'][i]**2))

    #gr filter
    for i in range(len(Filters['g']['Time'])):
	    if Filters['g']['Mag'][i] != -99.00 or Filters['r']['Mag'][i] != -99.00:
		    Sigma['gr']['grmag'].append(abs(Filters['g']['Mag'][i]) - abs(Filters['r']['Mag'][i]))
		    Sigma['gr']['grtime'].append(float((Filters['g']['Time'][i] + Filters['r']['Time'][i])/2))
		    Sigma['gr']['grerr'].append(math.sqrt(Filters['g']['Err'][i]**2 + Filters['r']['Err'][i]**2))

    #ri filter
    for i in range(len(Filters['i']['Time'])):
	    if Filters['r']['Mag'][i] != -99.00 or Filters['i']['Mag'][i] != -99.00:
		    Sigma['ri']['rimag'].append(abs(Filters['r']['Mag'][i]) - abs(Filters['i']['Mag'][i]))
		    Sigma['ri']['ritime'].append(float((Filters['r']['Time'][i] + Filters['i']['Time'][i])/2))
		    Sigma['ri']['rierr'].append(math.sqrt(Filters['r']['Err'][i]**2 + Filters['i']['Err'][i]**2))

    #iz filter
    for i in range(len(Filters['z']['Time'])):
	    if Filters['i']['Mag'][i] != -99.00 or Filters['z']['Mag'][i] != -99.00:
		    Sigma['iz']['izmag'].append(abs(Filters['i']['Mag'][i]) - abs(Filters['z']['Mag'][i]))
		    Sigma['iz']['iztime'].append(float((Filters['i']['Time'][i] + Filters['z']['Time'][i])/2))
		    Sigma['iz']['izerr'].append(math.sqrt(Filters['i']['Err'][i]**2 + Filters['z']['Err'][i]**2))

    #zu filter
    for i in range(len(Filters['u']['Time'])):
	    if Filters['z']['Mag'][i] != -99.00 or Filters['u']['Mag'][i] != 99.00:
		    Sigma['zu']['zumag'].append(abs(Filters['z']['Mag'][i]) - abs(Filters['u']['Mag'][i]))
		    Sigma['zu']['zutime'].append(float((Filters['z']['Time'][i] + Filters['u']['Time'][i])/2))
		    Sigma['zu']['zuerr'].append(math.sqrt(Filters['z']['Err'][i]**2 + Filters['u']['Err'][i]**2))

    #ur filter
    for i in range(len(Filters['r']['Time'])):
	    if Filters['r']['Mag'][i] != -99.00 or Filters['u']['Mag'][i] != 99.00:
		    Sigma['ur']['urmag'].append(abs(Filters['u']['Mag'][i]) - abs(Filters['r']['Mag'][i]))
		    Sigma['ur']['urtime'].append(float((Filters['u']['Time'][i] + Filters['r']['Time'][i])/2))
		    Sigma['ur']['urerr'].append(math.sqrt(Filters['u']['Err'][i]**2 + Filters['r']['Err'][i]**2))

    #ui filter
    for i in range(len(Filters['i']['Time'])):
	    if Filters['i']['Mag'][i] != -99.00 or Filters['u']['Mag'][i] != 99.00:
		    Sigma['ui']['uimag'].append(abs(Filters['u']['Mag'][i]) - abs(Filters['i']['Mag'][i]))
		    Sigma['ui']['uitime'].append(float((Filters['u']['Time'][i] + Filters['i']['Time'][i])/2))
		    Sigma['ui']['uierr'].append(math.sqrt(Filters['u']['Err'][i]**2 + Filters['i']['Err'][i]**2))

    #gi filter
    for i in range(len(Filters['i']['Time'])):
	    if Filters['g']['Mag'][i] != -99.00 or Filters['i']['Mag'][i] != 99.00:
		    Sigma['gi']['gimag'].append(abs(Filters['g']['Mag'][i]) - abs(Filters['i']['Mag'][i]))
		    Sigma['gi']['gitime'].append(float((Filters['g']['Time'][i] + Filters['i']['Time'][i])/2))
		    Sigma['gi']['gierr'].append(math.sqrt(Filters['g']['Err'][i]**2 + Filters['i']['Err'][i]**2))

    #gz filter
    for i in range(len(Filters['z']['Time'])):
	    if Filters['z']['Mag'][i] != -99.00 or Filters['g']['Mag'][i] != 99.00:
		    Sigma['gz']['gzmag'].append(abs(Filters['g']['Mag'][i]) - abs(Filters['z']['Mag'][i]))
		    Sigma['gz']['gztime'].append(float((Filters['g']['Time'][i] + Filters['z']['Time'][i])/2))
		    Sigma['gz']['gzerr'].append(math.sqrt(Filters['g']['Err'][i]**2 + Filters['g']['Err'][i]**2))

    #rz filter
    for i in range(len(Filters['z']['Time'])):
	    if Filters['z']['Mag'][i] != -99.00 or Filters['r']['Mag'][i] != 99.00:
		    Sigma['rz']['rzmag'].append(abs(Filters['r']['Mag'][i]) - abs(Filters['z']['Mag'][i]))
		    Sigma['rz']['rztime'].append(float((Filters['r']['Time'][i] + Filters['z']['Time'][i])/2))
		    Sigma['rz']['rzerr'].append(math.sqrt(Filters['r']['Err'][i]**2 + Filters['z']['Err'][i]**2))

    #------Chi2 Analysis------
    #Creating straight lines to loop.
    ugrange = np.linspace(min(Sigma['ug']['ugmag']),max(Sigma['ug']['ugmag']),500)
    grrange = np.linspace(min(Sigma['gr']['grmag']),max(Sigma['gr']['grmag']),500)
    rirange = np.linspace(min(Sigma['ri']['rimag']),max(Sigma['ri']['rimag']),500)
    izrange = np.linspace(min(Sigma['iz']['izmag']),max(Sigma['iz']['izmag']),500)
    zurange = np.linspace(min(Sigma['zu']['zumag']),max(Sigma['zu']['zumag']),500)
    urrange = np.linspace(min(Sigma['ur']['urmag']),max(Sigma['ur']['urmag']),500)
    uirange = np.linspace(min(Sigma['ui']['uimag']),max(Sigma['ui']['uimag']),500)
    girange = np.linspace(min(Sigma['gi']['gimag']),max(Sigma['gi']['gimag']),500)
    gzrange = np.linspace(min(Sigma['gz']['gzmag']),max(Sigma['gz']['gzmag']),500)
    rzrange = np.linspace(min(Sigma['rz']['rzmag']),max(Sigma['rz']['rzmag']),500)
	
    #Minimising chi2 for each filter difference
    #ug filter
    ugchi2min = 10E7
    for j in range(len(ugrange)):
        chi2min = np.sum(((ugrange[j]-Sigma['ug']['ugmag'])**2)/np.array(Sigma['ug']['ugerr'])**2)
        if chi2min < ugchi2min:
            ugchi2min = chi2min
        else: 
            pass
    #Dividing by the degrees of freedom (No. of elements - params)
    ugchi2min = ugchi2min/(len(Sigma['ug']['ugmag'])-1)

    #gr filter
    grchi2min = 10E7
    for j in range(len(grrange)):
        chi2min = np.sum(((grrange[j]-Sigma['gr']['grmag'])**2)/np.array(Sigma['gr']['grerr'])**2)
        if chi2min < grchi2min:
            grchi2min = chi2min
        else: 
            pass
    grchi2min = grchi2min/(len(Sigma['gr']['grmag'])-1)

    #ri filter
    richi2min = 10E7
    for j in range(len(rirange)):
	    chi2min = np.sum(((rirange[j]-Sigma['ri']['rimag'])**2)/np.array(Sigma['ri']['rierr'])**2)
	    if chi2min < richi2min:
	        richi2min = chi2min
	    else: 
	        pass
    richi2min = richi2min/(len(Sigma['ri']['rimag'])-1)

    #iz filter
    izchi2min = 10E7
    for j in range(len(izrange)):
	    chi2min = np.sum(((izrange[j]-Sigma['iz']['izmag'])**2)/np.array(Sigma['iz']['izerr'])**2)
	    if chi2min < izchi2min:
	        izchi2min = chi2min
	    else: 
	        pass
    izchi2min = izchi2min/(len(Sigma['iz']['izmag'])-1)

    #zu filter
    zuchi2min = 10E7
    for j in range(len(zurange)):
	    chi2min = np.sum(((zurange[j]-Sigma['zu']['zumag'])**2)/np.array(Sigma['zu']['zuerr'])**2)
	    if chi2min < zuchi2min:
	        zuchi2min = chi2min
	    else: 
	        pass
    zuchi2min = zuchi2min/(len(Sigma['zu']['zumag'])-1)

    #ur filter
    urchi2min = 10E7
    for j in range(len(urrange)):
	    chi2min = np.sum(((urrange[j]-Sigma['ur']['urmag'])**2)/np.array(Sigma['ur']['urerr'])**2)
	    if chi2min < urchi2min:
	        urchi2min = chi2min
	    else: 
	        pass
    urchi2min = urchi2min/(len(Sigma['ur']['urmag'])-1)

    #ui filter
    uichi2min = 10E7
    for j in range(len(uirange)):
	    chi2min = np.sum(((uirange[j]-Sigma['ui']['uimag'])**2)/np.array(Sigma['ui']['uierr'])**2)
	    if chi2min < uichi2min:
	        uichi2min = chi2min
	    else: 
	        pass
    uichi2min = uichi2min/(len(Sigma['ui']['uimag'])-1)

    #gi filter
    gichi2min = 10E7
    for j in range(len(girange)):
	    chi2min = np.sum(((girange[j]-Sigma['gi']['gimag'])**2)/np.array(Sigma['gi']['gierr'])**2)
	    if chi2min < gichi2min:
	        gichi2min = chi2min
	    else: 
	        pass
    gichi2min = gichi2min/(len(Sigma['gi']['gimag'])-1)

    #gz filter
    gzchi2min = 10E7
    for j in range(len(gzrange)):
	    chi2min = np.sum(((gzrange[j]-Sigma['gz']['gzmag'])**2)/np.array(Sigma['gz']['gzerr'])**2)
	    if chi2min < gzchi2min:
	        gzchi2min = chi2min
	    else: 
	        pass
    gzchi2min = gzchi2min/(len(Sigma['gz']['gzmag'])-1)

    #rz filter
    rzchi2min = 10E7
    for j in range(len(rzrange)):
	    chi2min = np.sum(((rzrange[j]-Sigma['rz']['rzmag'])**2)/np.array(Sigma['rz']['rzerr'])**2)
	    if chi2min < rzchi2min:
	        rzchi2min = chi2min
	    else: 
	        pass
    rzchi2min = rzchi2min/(len(Sigma['rz']['rzmag'])-1)

    #Total minimum chi-squared for all filters
    totchi = ugchi2min + grchi2min + richi2min + izchi2min + zuchi2min + urchi2min + uichi2min + gichi2min + gzchi2min + rzchi2min
             
    #Appending the chi2min array
    chisqmin = np.array([ugchi2min, grchi2min, richi2min, izchi2min, zuchi2min,
                 urchi2min, uichi2min, gichi2min, gzchi2min, rzchi2min, totchi])
                         
    return(chisqmin)


#======================================Main Body===============================#

#Starting timer from main body of the code
stime = time.time()
count = 0

#Looping through each id 
for item in Obid['ids']: 
	chimin = chi2(item)
	chisqmin.append(chimin)
	print(item)
	count += 1

#Pickling data and creating the outfile
filename = 'chi2min.pickle'
outfile = open(filename, 'wb')
pickle.dump(chisqmin, outfile)
outfile.close()

#Calculating and printing run time
ftime = time.time()
print(ftime - stime)
'''
#==============================================================================#
#=================================End of File==================================#
#==============================================================================#
