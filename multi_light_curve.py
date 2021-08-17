#===================================Pramble====================================#
#Program Description: 
#1) Program first computes a list of file paths where the initial data can be found,
#   this calls the multi_curves function and cycles through the directory given.
#2) It uses the file_paths array to attempt to compute the fire_light_curve.py 
#   program. This produces a plot of magnitudes and the B-R colour-index, as well 
#   as an updated copy of the data used, and places them in the specified location.
#3) Depending if this can be done, this program produces 4 files with the locations
#   of files that have either passed or failed (and their reason for doing so).
#4) Then closes all open files and prints pass and fail numbers to terminal.

#Inputs: 1) -f <location for initial data> 
#        2) -d <destination for figures and processed data files>

#Packages: pip install os argparse logging pdb numpy astropy datetime

#Run: python3 multi_light_curve.py -f <location for initial data> -d <destination for figures and processed data files>

#Outputs: 1) See fire_light_curve.py. I'll write something better here
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
#import datetime

#Local Modules
import dfn_utils
#import photutils

#Importing files
import fire_light_curve

#==============================Function Definition=============================#

#Finding the most recent .ecsv file
'''
#Old function trying to stream line below, possibly not desired result
def recent_file(time):
	now_time = datetime.datetime(2021,7,1,0,0,0)
	file_time = datetime.datetime(int(time[0:4]), int(time[5:7]), int(time[8:10]), int(time[11:13]), int(time[14:16]), int(time[17:19]))
	return(abs(now_time-file_time))
''' 
def recent_file(last_update_time):
    now_time = Time('2021-07-01')
    now_time.format = 'datetime'
    file_time = Time(last_update_time)
    file_time.format = 'datetime'
    return(abs(now_time-file_time))

#Returns the file paths for all .ecsv files in the directory
def multi_curves(current_dir):
    #Initialising list files and paths
    file_paths = []
    fail_counter = 0
    #Cycling through events in the photometry directory
    with os.scandir(current_dir) as entries:
        for entry in entries:
            #Cycling through cameras that observed each event
            with os.scandir(entry) as obs:
                for ob in obs:
                    #Setting current path and initialising directory list
                    current_cam = os.path.join(ob)
                    files_in_dir = []
                    times_in_dir = []
                    #Finding the .ecsv files for each observation
                    for f_name in os.listdir(current_cam):
                        #Finding the most recent file in directory
                        if f_name.endswith('.ecsv'):
                            current_file = os.path.join(current_cam, f_name)
                            #Checking if file contains information
                            if os.path.getsize(current_file) > 50:
                                table = Table.read(current_file, format='ascii.ecsv', guess=False, delimiter=',')
                                last_update_time = table.meta['point_picking_write_time']
                                times_in_dir.append(recent_file(last_update_time))
                                files_in_dir.append(current_file)
                            else:
                                fail_counter += 1
                                pass
                    #Printing most recent file path to array
                    if len(times_in_dir) < 1:
                        pass
                    else:
                        least_time = times_in_dir.index(min(times_in_dir))
                        file_paths.append(files_in_dir[least_time])
                    #pdb.set_trace()
    #Prints fails to terminal
    print(fail_counter, len(file_paths))
    return(file_paths)
    					
#========================Passing directory into Program========================#

#Passing directory to serach for files
parser = argparse.ArgumentParser(description='Finding data in directory.')
parser.add_argument("-f", "--file", help="File name passed into program", type=str)
parser.add_argument("-d", "--destination", help="Figure destination for output", type=str)

#Setting arguments to pass into the file
args = parser.parse_args()	
#destination = '../Figures'

current_dir = args.file
destination = args.destination

#Creating file paths to cycle through
file_paths = multi_curves(current_dir)

#Creating counters
passes = 0
tfails = 0
ofails = 0
mfails = 0

#opening file to write errors
otherrors = open('otherrors.txt', 'w')
timerrors = open('timerrors.txt', 'w')
filerrors = open('filerrors.txt', 'w')
passfiles = open('passfiles.txt', 'w')

#Computing photometry for events
for i in range(len(file_paths)):
	try:
		fire_light_curve.main_disk_io(file_paths[i], plot=True, destination=destination)
		#passfiles.write('#===========Observation Calculation Complete===========#\n')
		passfiles.write(file_paths[i]+'\n')
		passes += 1
	except dfn_utils.WrongTableTypeExceptions
		#timerrors.write('#===================Not Enouf Timing===================#\n')
		timerrors.write(file_paths[i]+'\n')
		tfails += 1
	except FileNotFoundError:
	    #filerrors.write('#====================Missing a File====================#\n')
	    filerrors.write(file_paths[i]+'\n')
	    mfails += 1
	except Exception as e:
		#otherrors.write('#==============Something Else Went Wrong===============#\n')
		otherrors.write(str(e)+'\n')
		otherrors.write(file_paths[i]+'\n')
		ofails += 1

#Print final numbers to file
otherrors.write('Passes = ' + str(passes) + '\nTiming = ' + str(tfails) + '\nOther = ' + str(ofails) + '\nMissing file = ' + str(mfails))

#Closing open files
otherrors.close()
timerrors.close()
filerrors.close()
passfiles.close()

#Print final numbers to terminal		
print('Passes = ', passes)
print('Timing Failures = ', tfails)
print('Other Failures = ', ofails)
print('There was a missing file', mfails, 'times')

#==============================================================================#
#=================================END OF FILE==================================#
#==============================================================================#
