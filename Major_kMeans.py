#===================================Pramble====================================#
#Program Description: This program performs the k-Means algorithm on the data. 
#   It prompts the user for input about the number of clusters, step of k-Means
#   to perform and cut off values for the data. It then prints out information 
#   after each step and a plot of the new clusters. 
#
#Inputs: No inputs to terminal.
#
#Packages: pip install numpy math pdb matplotlib pickle random 
#
#Run: python3 Major_kMeans.py
#
#Outputs: 1) Asks for user inputs in the terminal, prompts given
#         2) Prints initial centeroids
#         3) Prints new centeroid after each iteration of k-Means with the 
#            number of points in each cluster. It also gives the current 
#            iteration of k-Means, the total number of valid points and the 
#            cummulative total change in the centeroid values.
#         4) After the completion of i steps it prints to terminal these total 
#            changes for each cluster in one array. 
#         5) Creates a plot of the cluster which can either be show or save to 
#            file (from within the code).
#==============================Importing Packages==============================#

import numpy as np
import math
import pdb
import matplotlib.pyplot as plt
import pickle
import random as rn

#=======================Opening and Reading Catalogue File=====================#

f0 = open('stripe82.dat','r')
stripe82 = f0.read()
stripe82 = stripe82.splitlines()

#ID Dictionary 
Obid = {'ids':[]}

#Peroid Array
period = []

#Populating Obid and Array
for i in stripe82:
    x = i.split()
    if x[0][0] != "#":
        Obid['ids'].append(x[0])
        period.append(float(x[3]))
            
#Closing Catalogue File      
f0.close()

#Pickle in the Chi-squared data
infile = open('chi2min.pickle', 'rb')
chisqmin = pickle.load(infile)
infile.close()

#Initialising 
totchi2 = []
#Populating array
for i in range(len(period)):
    totchi2.append(chisqmin[i][-1])
    
#==============================Plotting Variables==============================#

#Colour array
colours = ['orange', 'yellow', 'green', 'lime', 'cyan', 'midnightblue', 'blue', 'purple', 'magenta', 'crimson']
#Marker array
markers = ['o', 'v', '^', 's', 'P', 'X', '<', '>', '*', 'D']

#=============================Function Definitions=============================#

#Euclidean distance
def euclid(x1,x2,y1,y2):
    dist = math.sqrt((x2-x1)**2+(y2-y1)**2)
    return dist
    
#Populating Clusters
def pop_clusters(centres,point,Num_of_centres):
    dist = 10E10
    for i in range(len(centres)):
        dist_2_cluster_cen = euclid(centres[i][0],point[0],centres[i][1],point[1])
        if dist_2_cluster_cen < dist:
            dist = dist_2_cluster_cen
            cluster = i
        else:
            pass
    return(int(cluster))
    
#==========================Plotting the Initial Data===========================#

#Plotting Initial zoom
plt.scatter(period, totchi2)
plt.xlabel('Variable Period in Days')
plt.ylabel('Chi-Squared Minimum')
plt.title('Initial Data Structures')
#plt.xlim(0,3500)
#plt.ylim(0,2.5E6)
plt.show()

#=================================User Inputs==================================#

#Number of clusters
print("So, this will break for inputs larger than 10, but it's because I've only ")
print("got 10 colours and markers in my array for plotting. This can be changed ")
print("if needed.")
Num_of_centres = int(input("How many clusters would you like?:"))

#Number of steps of k-Means
k_steps = int(input("How many steps of k-Means would you like to perform?:"))

#Asking for data cuts
#Printing info to the user for variable period cuts
print("For this data, the maximum variability period is about 17,500 days. If you ")
print("would like to see all of the data, I would recommend a cut of 18,000. Because ")
print("of the cadence of these observations there is a natural cut at 3,500 days. I ")
print("would suggest cutting at that point.")
period_cut = float(input("What is the maximum variability period you want (Days)?: "))

#Printing info to the user for chi-square min cuts
print("The Chi-squared minimum is more difficult to judge about where cuts need ")
print("to be made. For this data there seems to be a natural cut-off at approximately ")
print("2E+6. So, I've choosen a cut-off of 2.5E+6 to catch stragglers on the edge. ")
print("I choosen it to be this large and still be valid, as sometimes it's only ")
print("one or two of the Chi-squared values which inflates the total value. The ")
print("maximum value in this dataset is about 1.2E+7. Use this number if you want ")
print("to see all of the data, in this data set.")
chimin_cut = float(input("What is the maximum Chi-squared value you would like to look at?: "))

#=====================Insuring centres are in Choosen Data=====================#

#Initialising centres array
centres = []
k = 0 

#Calculating the Centres
while k < Num_of_centres:
    centre_index = rn.randint(0,len(period)-1)
    if period[centre_index] > period_cut or totchi2[centre_index] > chimin_cut:
        centre_index = rn.randint(0,len(period)-1)
    else:
        centres.append(np.array([period[centre_index],totchi2[centre_index]]))
        k += 1

#Printing the initial centeroids
for i in range(len(centres)):
    print("The initial centre of centeroid ", i+1, " is: ", centres[i])

#==============================Performing k-Means==============================#

#Creating counter for the k-Means steps
counter = 0
#Saving total centre changes
tot_cen_change = []

while counter < k_steps:

    #Creating the array to store the clusters
    #Clusters = np.zeros(Num_of_centres,dtype=object)
    Clusters = []
    for i in range(0,Num_of_centres):
        Clusters.append([])
    
    #Cluster centres before k-Means is starrted
    init_cens = centres
    #print(init_cens)
    centres = []
    #Initialising total points plotted
    tot_points = 0 
    
    #Determining cluster
    for i in range(len(period)):
        if period[i] < period_cut and totchi2[i] < chimin_cut:
            point = np.array([float(period[i]),float(totchi2[i])])
            cluster_num = pop_clusters(init_cens,point,Num_of_centres)
            Clusters[cluster_num].append(np.array([period[i],totchi2[i]]))
        else:
            pass
        
    #Determining the new Centres
    for i in range(len(Clusters)):
        #Initialising new cluster centres
        newx = 0
        newy = 0
        for j in range(len(Clusters[i])):
            newx += Clusters[i][j][0]
            newy += Clusters[i][j][1]
        newx = newx/len(Clusters[i])
        newy = newy/len(Clusters[i])
        new_cens = np.array([newx,newy])
        centres.append(new_cens)
        tot_points += len(Clusters[i])
        print("Cluster ", i+1 , " centeroid, with ", new_cens, " with ", len(Clusters[i]), "data points.")

    #Printing total points plotted
    print("Number of valid points: ", tot_points)
    #Initialising total centeroid differences
    total_dist = 0 

    #Calculating the differences between old and new centres
    for i in range(len(centres)):
        total_dist += euclid(init_cens[i][0],centres[i][0],init_cens[i][1],centres[i][1])
    
    #Informing the user of change in centres  
    print("After ", counter+1, " iterations of the k-Means algorithm")  
    print("Total change in all centres:", round(total_dist,4))
    
    #Saving total change each iteration
    tot_cen_change.append(round(total_dist,2))
    
    #Plotting the new clusters found
    for i in range(len(Clusters)):
        #Resetting x and y co-ordinates
        xcords = []
        ycords = []
        for j in range(len(Clusters[i])):
            xcords.append(Clusters[i][j][0])
            ycords.append(Clusters[i][j][1])
        #Plotting x and y for each cluster
        plt.scatter(xcords,ycords,color=colours[i],marker=markers[i],facecolors='none', label='C{0}'.format(i+1))
        
    #Plotting fluff
    plt.xlabel('Variable Period in Days')
    plt.ylabel('Chi-Squared Minimum')
    plt.title('New Clusters after k-Means Step')
    plt.legend(loc='upper left', bbox_to_anchor=(0.95,1))
    plt.savefig('kMeans_Clusters/ExSmlChi_kMeans_Step_{0}.png'.format(counter+1))
    plt.close()
    #plt.show()
    print('Figure Saved')
    
    #Iterating through loop
    counter += 1

#Printing total change in centres for all iterations of k-Means
print("Total distance change in ceteroids for all iterations of k-Means:")
print(tot_cen_change)

#==============================================================================#
#=================================End of File==================================#
#==============================================================================#

