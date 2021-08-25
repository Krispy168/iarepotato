#===================================Pramble====================================#
#Program Description: 
#1) Reads in the files to create a list of observations
#2) For each observation it creates and saves a figure of the magnitudes for each 
#   filter (RGB) and the colour-indices. 
#3) Prints number of fails and error messages to terminal (shouldn't be any)

#Note: A runtime warning appears when running this code, not sure why. See error below:
#/home/krispy/anaconda3/lib/python3.7/site-packages/astropy/table/column.py:1020: RuntimeWarning: invalid value encountered in less
#  result = getattr(super(), op)(other)
#Note2: Global x-axis time would be great, need to figure out
#Note3: Possibly play with the spacing between the plots

#Inputs: 1) -l <location for data files>  
#        2) -d <destination for figures>

#Packages: pip install os argparse logging pdb numpy astropy matplotlib

#Run: python3 comparison.py -l <location for data files> -d <destination for figures>

#Outputs: A figure for each observation listed in the location file, which is the 
#         side by side of the magnitudes and colour-indices.
#==============================Importing Packages==============================#

#Import Reading and Error Modules
import os
import argparse
import logging
import pdb

#Import Science Modules
import numpy as np
from astropy.table import Table, Column
from astropy.time import Time
import matplotlib.pyplot as plt

#Importing files (Not sure this is completely necessary)
import fire_light_curve

#==============================Function Definition=============================#

#Creating an array of .ecsv file from the directory listed
def ecsv_files(location):
    #Initialising the file path arrays
    file_paths = []
    #Extracting all the .ecsv files in the specified directory
    for f_name in os.listdir(location):
        if f_name.endswith('.ecsv'):
            cur_file = os.path.join(location, f_name)
            file_paths.append(cur_file)
    return(file_paths)

#Plotting all the colour indices 
def colour_indices_compare(input_file, destination):
    logger = logging.getLogger()
    #Reading in the input file
    table = Table.read(input_file, format='ascii.ecsv', guess=False, delimiter=',')
    #Setting the colour filters
    SPECTRAL_BAND_DE_BAYER = {'R': 'red', 'G': 'green', 'B': 'blue'}
    #Setting the different colour indices
    C_INDICES = {'B-R': ['cyan','s'], 'B-G': ['lime','D'], 'G-R': ['orangered','o']}
    #Calculating the x-axis
    ref_time = Time(table['datetime'][0])
    table['reltime'] = (Time(table['datetime']) - ref_time).sec
    #Clearing the previous figure and setting the new one
    plt.close()
    fig, (ax1, ax2) = plt.subplots(1,2)
    
    #Setting the magnitude plot
    for spectral_band in SPECTRAL_BAND_DE_BAYER:
        app_mag_colname = 'm_' + spectral_band
        saturation_colname = 'saturated_flag_' + spectral_band
        snr_colname = 'SNR_' + spectral_band
        if not app_mag_colname in table.colnames:
            logger.warning(f'{spectral_band} not available for plotting')
            continue
        #Marking the saturated points
        sat_mask = (table[saturation_colname] == 1)
        n_tot = np.count_nonzero(~np.isnan(table[saturation_colname]))
        n_sat = np.nansum(table[saturation_colname][sat_mask])
        n_non_sat = n_tot - n_sat
        #Marking the low SNR points 
        low_snr = (table[snr_colname] < 25)
        #Plot unsaturated data points
        ax1.scatter(table['reltime'][~sat_mask], table[app_mag_colname][~sat_mask],
                label=f'{spectral_band} n={n_non_sat:.0f}',
                color=spectral_band.lower())
        #Plot saturated data points
        ax1.scatter(table['reltime'][sat_mask], table[app_mag_colname][sat_mask],
                label=f'{spectral_band} saturated n={n_sat:.0f}',
                marker='x', color=spectral_band.lower())
        #Plot low SNR points 
        ax1.scatter(table['reltime'][low_snr], table[app_mag_colname][low_snr],
                marker='.', color='white', s=75) 
    #Setting the colour-indices plot
    for c_index in C_INDICES:
        c_index_colname = 'colour_index_' + c_index
        if not c_index_colname in table.colnames:
            logger.warning(f'{c_index} not available for plotting')
            continue
        #Plotting the colour indices
        ax2.scatter(table['reltime'], table[c_index_colname], label=f'{c_index}',
                marker=C_INDICES[c_index][1], color=C_INDICES[c_index][0])
    
    #Setting global figure options
    fig.suptitle(table.meta['event_codename']+'_'+table.meta['location']+' Magnitude and Colour-Indices Comparison')
    #Dirty way to note that unfilled circles are low SNR
    ax1.scatter([],[],marker='.',color='white',label='Unfilled point is Low SNR')
    #Sub-figure 1 magnitude comparison plot
    ax1.invert_yaxis()
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Apparent Magnitude')
    ax1.grid()  
    ax1.legend(loc='upper left', bbox_to_anchor=(1,1))
    #Sub-figure 2 colour-indices comparison plot
    ax2.invert_yaxis()
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Colour-Indices')
    ax2.grid()
    ax2.legend(loc='lower right', bbox_to_anchor=(0,0))
    #Printing or saving figure
    #plt.show()
    plt.subplots_adjust(wspace=0.6)
    fig.set_size_inches(15,5)
    fig.savefig(destination+table.meta['event_codename']+'_'+table.meta['location']+'.png', dpi=200)
    
#==========================Plotting the Comparisons============================#

#Setting the logger
logger = logging.getLogger()

#Arguements to pass into the program 
parser = argparse.ArgumentParser(description='Makes comparison figure')
parser.add_argument("-l", "--location", type=str, help="File location", required=True)
parser.add_argument("-d", "--destination", type=str, help="Place for comparison figures", default=None)

args = parser.parse_args()

cur_location = args.location
sav_location = args.destination

#Extracting the .ecsv files to compare
file_paths = ecsv_files(cur_location)
#Error counters
fails = 0
error_type = []
#new_location = '../Figures/15_2015-11-13_172458_DSC_1506-G_DN151113_03_2018-01-02_175311_patrick_nocomment.ecsv'

#Making the graphs
for path in file_paths:
    try:
        colour_indices_compare(path, sav_location)
        print('Comparison figure saved')
    except Exception as e:
        error_type.append(e)
        fails += 1
    #Uncomment for debugging and graph experiments
    #pdb.set_trace()

#Printing errors to terminal
if fails > 0:
    print(fails)
    print(error_type)

#==============================================================================#
#=================================END OF FILE==================================#
#==============================================================================#
