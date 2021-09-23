#===================================Pramble====================================#
#Program Description: 
#1) Takes the input file locations to the obs_colab function 
#2) This function finds the corresponding files and combines them into one
#3) Within the obs_colab function the air mass correction (correction_add) for 
#   each filter is calculated and added to the new file. The air mass correction 
#   for each filter can be changed within the code, via the CORRECTION_FACTOR 
#   dictionary.
#4) Also within the obs_colab function the min and max values of key variables are
#   added to the meta data of the code (add_meta_data). This is done for each 
#   filter RGB, aswell as the row index for each of these loctions if further
#   information is required at that time. 
#5) Renames the file and saves to the new location given

#Inputs: 1) -v <location for visual data files>  
#        2) -t <location for triangulation files>
#        3) -l <destination for new files>

#Packages: pip install os argparse logging pdb numpy astropy matplotlib pandas

#Run: python3 event_colab.py -v <location for visual data files> -t <location for triangulation files> -l <destination for new files>

#Outputs: Combines and renames the files for each observation 
#==============================Importing Packages==============================#

#Import Reading and Error Modules
import os
import argparse
import logging
import pdb

#Import Science Modules
import numpy as np
from astropy.table import Table, Column, hstack
from astropy.time import Time
import matplotlib.pyplot as plt
import pandas as pd

#==============================Function Definition=============================#
'''
#Calculating the zenith angle
def zenith(alt):
    #Convrting altitude to zenith and degrees to radians
    zenith = (90-alt)*(np.pi/180)
    return zenith
    
#Air mass correction
def air_mass_correction(alt, height):
    #Defining constants and arrays for the function
    rad_E = 6.371E+6 #Ave. Earth Radius (m)
    zen = zenith(alt)
    R_on_high = rad_E/height
    high_on_R = height/rad_E
    cozy = np.cos(zen)
    #Calculating the correction
    AM_correction = R_on_high*np.sqrt((cozy**2)+(2*(high_on_R))+(high_on_R**2))-R_on_high*cozy
    return AM_correction
'''
#Air mass correction 
def air_mass_correction(dist, height):
    AM_correction = dist/height
    return AM_correction
    
#Adding and populating new correction columns 
def correction_add(table):
    #Correction factor per air mass
    CORRECTION_FACTOR = {'R': 0.1,
                         'G': 0.1,
                         'B': 0.1} #Change this for correct values
    #Setting the arrays for calculation
    dist = table['range']
    height = table['height']
    #Calculating the air mass correction
    AM_correction = air_mass_correction(dist,height)
    #Writing columns and data to file
    for key in CORRECTION_FACTOR:
        #Getting the uncorrected magnitudes
        un_mag = table['m_'+key]
        #Setting column name and data
        column_name = str(key+'_airmass_correct')
        column_data = un_mag - (AM_correction*CORRECTION_FACTOR[key])
        table.add_column(column_data, name=column_name)
    #Returns the new table columns and values
    return table #I don't think this is needed

#==============================================================================#

#Adding to the meta data
def add_meta_data(table):
    #Initialising filters 
    FILTERS = {'R', 'G', 'B'}
    #Finding the length of the observation
    table.meta['length_obs'] = (Time(table['datetime'][-1])-Time(table['datetime'][0])).sec
    #Converting the astropy table to a pandas data structure
    df = table.to_pandas()
    #All variables at max brightness for each filter
    for key in FILTERS:
        #Bypassing the NaNs
        no_nan = np.where(np.isnan(df['m_'+key])==False)
        #Setting the min and max index
        min_index = df[df['m_'+key]==max(df['m_'+key][no_nan[0]])].index
        max_index = df[df['m_'+key]==min(df['m_'+key][no_nan[0]])].index
        #Setting the meta data variable names
        time = 'time_max_' + key
        height = 'height_max_' + key
        max_value = 'max_mag_' + key
        min_value = 'min_mag_' + key
        delta_time = 'elapsed_time_' + key
        saturated = 'saturated_' + key
        max_place = 'max_position_' + key
        min_place = 'min_position_' + key
        #Setting the value of the defined variables
        table.meta[time] = df['datetime'][max_index[0]]
        table.meta[height] = df['height'][max_index[0]]
        table.meta[max_value] = df['m_'+key][max_index[0]] #min(table['m_'+key])
        table.meta[min_value] = df['m_'+key][min_index[0]] #max(table['m_'+key])
        table.meta[delta_time] = (Time(df['datetime'][max_index[0]])-Time(df['datetime'][0])).sec
        table.meta[saturated] = df['saturated_flag_'+key][max_index[0]]
        table.meta[max_place] = max_index[0]
        table.meta[min_place] = min_index[0]
    #Message to terminal 
    print("I think this is all you need, but just incase I've included the locations")
    print("in the data where the max and min values occur for each filter") 

#==============================================================================#

#Merging visual and triangulation data for single camera observations
def obs_colab(visual, triangulation, save):
    #Listing all the files
    for f_name in os.listdir(visual):
        if f_name.endswith('.ecsv'):
            #Reading in the tables
            vis_data = os.path.join(visual,f_name)
            vis_table = Table.read(vis_data, format='ascii.ecsv', guess=False, delimiter=',')
            tri_data = os.path.join(triangulation,f_name)
            tri_table = Table.read(tri_data, format='ascii.ecsv', guess=False, delimiter=',')
            #Removing columns with the same names, and therefore data
            tri_table.remove_columns([col for col in tri_table.colnames if col in vis_table.colnames])
            #Stacking the tables horizontally
            try:
                tot_table = hstack([vis_table,tri_table], join_type='exact')
            except Exception as e:
                print(e)
            #Adding the air mass corrections
            correction_add(tot_table)
            #Adding extra meta data 
            add_meta_data(tot_table)
            #Setting the file name, and saving file.
            codename = tot_table.meta['event_codename']
            location = tot_table.meta['location']
            fname = save+tot_table.meta['event_codename']+'_'+tot_table.meta['location']+'.ecsv'
            tot_table.write(fname, format='ascii.ecsv', delimiter=',', overwrite=True)
            #Message to terminal
            print("Table should be merged correctly, don't worry about those errors")
            print("They're more like guidelines anyways")
            
#===============================Execution of Code==============================#

#Setting the logger
logger = logging.getLogger()

#Arguements to pass into the program 
parser = argparse.ArgumentParser(description='Makes comparison figure')
parser.add_argument("-v", "--visual_data", type=str, help="Visual data file location", required=True)
parser.add_argument("-t", "--triangulation", type=str, help="Triangulation data file location", required=True)
parser.add_argument("-l", "--location", type=str, help="Saved data file location", required=True)

args = parser.parse_args()

vis_data = args.visual_data
tri_data = args.triangulation
save_loc = args.location

#idea checks at this stage
obs_colab(vis_data, tri_data, save_loc)

#==============================================================================#
#=================================END OF FILE==================================#
#==============================================================================#
