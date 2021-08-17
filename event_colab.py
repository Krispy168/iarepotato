#===================================Pramble====================================#
#Program Description: 
#1) 

#Inputs: 1) -l <location for data files>  
#        2) -d <destination for figures>

#Packages: pip install os argparse pdb numpy astropy matplotlib

#Run: python3 event_colab.py -l <location for data files> -d <destination for figures>

#Outputs: 
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
#from mpl_toolkit.axes_grid1 import host_subplot
#import datetime

#Local Modules
#import dfn_utils
#import photutils

#Importing files
import fire_light_curve

#==============================Function Definition=============================#

#Air mass correction
def air_mass_correction():
    #This then needs to be different for each filter. 

#Bringing observations together
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
            #Setting the file name, and saving file.
            codename = tot_table.meta['event_codename']
            location = tot_table.meta['location']
            fname = save+tot_table.meta['event_codename']+'_Attempt2_'+tot_table.meta['location']+'.ecsv'
            tot_table.write(fname, format='ascii.ecsv', delimiter=',', overwrite=True)
            #Message to terminal
            print("Table should be merged correctly, don't worry about those errors")
            print("They're more like guidelines anyways")

#==============================================================================#

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
