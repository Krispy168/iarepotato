#!/usr/bin/env python

#===================================Pramble====================================#
#Program Description: 
#1) Reads in ecsv file
#2) Cleans table from ecsv if previously run
#3) Determines if the de Bruijn encoding firmware is usable
#4) Finds the .NEF file for raw data and converts to usable data structure
#5) ??Finds .jpeg for visualisation??
#6) Runs main function
#   6.1) Initialises the de Bruijn sequence
#   6.2) Computes the photometry, different for different apetures
#       6.2.1) Creates light curve from data 
#   6.3) Adds computed data to the .ecsv table
#   6.4) Calibrates photometry ??Should this be done earlier??
#   6.5) Calculates the colour-indices
#       6.5.1) Currently B-R, can be computed for all RGB

#Inputs: Needs at least .ecsv and .NEF I think

#Packages: pip install os sys datetime copy logging argparse numpy imageio astropy photutils

#Install: ./install_still_photometry_tools.sh -d 'location'

#Run: python3 fire_light_curve.py -p 'file_name'.ecsv
#Modified Run: pyhton3 fire_light_curve.py -p ../DFN_photometry

#Outputs: 1) Update input ecsv file
#         2) Light curve figure for RGB of fireball
#         3) B-R colour-index figure, adjustable for all colour-indices
#==============================Importing Packages==============================#

import os
import sys
import datetime
from copy import deepcopy
import logging
import pdb

# science modules
import numpy as np
from imageio import imread
from astropy.table import Table, Column
from astropy import units as u
from astropy.time import Time#, TimeDelta

from photutils import CircularAperture, RectangularAperture#, CircularAnnulus, EllipticalAperture
from photutils import aperture_photometry#, get_cutouts
#import photutils

import matplotlib.pyplot as plt

# local modules
from de_bruijn import de_bruijn_sequence
import extract_photometry_star_cata as cataphot
import dfn_utils
import timing_corrections
import fitsConversion

#================================Authors Information===========================#

__author__ = "Hadrien Devillepoix"
__copyright__ = "Copyright 2016-2021, Desert Fireball Network"
__license__ = "MIT"
__version__ = "3.0"
__scriptName__ = "fire_light_curve.py"

#============================Dictionary Initialisation=========================#

SCRIPT_AND_VERSION =  __scriptName__ + ' ' + __version__
  

LC_SHUTTER_INVERSION_FLAG = 'shutter_inversion'

# list of micro-controller encoding firmwares supported
FIRMWARES_SUPPORTED = ['april2014']


# De Bruijn sequence
DB_SEQ = de_bruijn_sequence(2, 9).de_bruijn_sequence


# dictionary of spectral bands vs de-Bayering methods
SPECTRAL_BAND_DE_BAYER = {'R': 'red',
                          'G': 'green',
                          'B': 'blue'}

# dictionary of correction factors for pixel exposure in different spectral bands
SPECTRAL_BAND_COLLECTING_FRACTION_INVERSE = {'R': 4,
                                        'G': 1,
                                        'B': 4}

#================================Class Definition==============================#

class EncodingInversionError(Exception):
    pass
        

class DBPoint():
    '''
    Class to handle easy navigation in the DB sequence
    '''
    def __init__(self, DB_index, position):
        # de Bruijn index (0-300)
        self.DB_index = int(DB_index)
        # start or end
        self.position = position
    
    def move_forward(self):
        if self.position == 'S':
            self.position = 'E'
        else:
            self.position = 'S'
            self.DB_index += 1
            
    def get_prev(self):
        if self.position == 'E':
            new_position = 'S'
            new_DB_index = self.DB_index
        else:
            new_position = 'E'
            new_DB_index = self.DB_index -1
        return DBPoint(new_DB_index, new_position)
    
    def is_start(self):
        return self.position == 'S'
    
    def __str__(self): 
        return str(self.DB_index) + self.position
    
#================================Function Definition===========================#
 
#Query Hadrien, Finding the bounds of the image??
def get_origins(pos1, pos2):
    min_x, min_y, _, _ = get_bounds(pos1, pos2)
    return (min_x, min_y)


def get_bounds(data):    
    min_x = int(np.min(data['x_image']))
    max_x = int(np.max(data['x_image']))
    min_y = int(np.min(data['y_image']))
    max_y = int(np.max(data['y_image']))
    
    return min_x, min_y, max_x, max_y


#===============Crop of science image around area of interest==================#
def get_bounding_box(scidata, data):
    # crop of science image around area of interest
    min_x, min_y, max_x, max_y = get_bounds(data)
    return scidata[min_y:max_y, min_x:max_x]

#====================Finding saturated pixels in the data======================#
def get_saturated_pixels(scidata, sat_adu=16383):
    
    saturated_pixels = np.where(16383 == scidata)
    
    return saturated_pixels

#==============Mapping errors of saturated pixels and neighbours===============#
def get_error_map(scidata, sat_adu=16383):
    error_map = np.ones(np.shape(scidata))
    sat = np.where(sat_adu == scidata)
    #sat_pix_list = np.transpose(sat)
    #error_map[sat] = sat_adu #*10
    
    # anti-blooming exponential meodel
    # error = sat_adu * exp(beta * N)
    # with N: number of ditect neighbors saturated
    # beta: anti-blooming coefficient: unknown, start with 1
    
    
    for i,j in zip(sat[0],sat[1]):
        # for each saturated pixel, find out how many of its direct neighbors are sturated
        number_sat_neighbors = 0
        try:
            for k in [i-1,i+1]:
                for l in [j-1,j+1]:
                    if error_map[k][l] == sat_adu:
                        number_sat_neighbors += 1
        except IndexError:
            number_sat_neighbors = 0
                    
        # apply empirical correction
        error_map[i][j] = sat_adu * np.exp(1 * number_sat_neighbors)
    
    return error_map


def light_profile(scidata, table, sat_adu=16383):
    '''
    Generate a zero-order light curve
    '''
    x_centers = np.array((table['x_image'][1:] + table['x_image'][:-1])/2).astype(int)
    y_centers = np.array((table['y_image'][1:] + table['y_image'][:-1])/2).astype(int)
    zero_order_brightness = np.asarray([scidata[y][x] for x,y in zip(x_centers, y_centers)])
    return zero_order_brightness


def DB_is_next(row1, row2):
    '''
    Determines if a De Bruijn status change follows another
    row1: supposedely previous DB status change
    row2: supposedely next DB status change
    '''
    db1 = DBPoint(row1['de_bruijn_sequence_element_index'], row1['dash_start_end'])
    db2 = DBPoint(row2['de_bruijn_sequence_element_index'], row2['dash_start_end'])
    db1.move_forward()
    return (db1.DB_index == db2.DB_index) and (db1.position == db2.position)


def table_matches_image(table, scidata):
    '''
    check that the image size matches the image bounds recorded in the data
    '''
    tab_nax1 = int(table.meta['NAXIS1'])
    tab_nax2 = int(table.meta['NAXIS2'])
    max_im_y, max_im_x = np.shape(scidata)
    
    match = ((tab_nax1 == max_im_x) and (tab_nax2 == max_im_y))
    factor = 1.0 * tab_nax1 / max_im_x
    return match, factor

def is_inverted(scidata, table, sat_adu=16383):
    '''
    Determines if a fireball is inverted just by looking at the zero-order light curve
    (pixel values in the middle in the dashes)
    '''
    
    # make a copy
    temp_table = Table(table)
    
    match, factor = table_matches_image(table, scidata)
    if not match:   
        temp_table['x_image'] /= factor
        temp_table['y_image'] /= factor
    
    supposed_to_be_closed = []
    supposed_to_be_open = []
    for r1,r2 in zip(temp_table[:-1],temp_table[1:]):
        if DB_is_next(r1, r2):
            zero_order_light = scidata[int((r1['y_image']+r2['y_image'])/2)][int((r1['x_image']+r2['x_image'])/2)]
            if r1['dash_start_end'] == 'S':
                supposed_to_be_open.append(zero_order_light)
            else:
                supposed_to_be_closed.append(zero_order_light)
    c = np.mean(supposed_to_be_closed)
    o = np.mean(supposed_to_be_open)   
    return c > o, o, c

#==========================Finding shutter open time===========================#
def LC_open_time(table, start_index, end_index, inverted):
    '''
    Determine how long the LC shutter is open for
    '''
    
    leak_factor = 0.04
    
    if table.meta['firmware'] == 'april2014':
        one_S = 0.06
        one_E = 0.1 - one_S
        zero_S = 0.02
        zero_E = 0.1 - zero_S
        if inverted:
            one_S *= leak_factor
            zero_S *= leak_factor
        else:
            one_E *= leak_factor
            zero_E *= leak_factor
            
        DB_iterator = DBPoint(table[start_index]['de_bruijn_sequence_element_index'],
                        table[start_index]['dash_start_end'])    
        last_DB = DBPoint(table[end_index]['de_bruijn_sequence_element_index'],
                        table[end_index]['dash_start_end'])
        
        
        exposure_time = 0.
        while (DB_iterator.DB_index != last_DB.DB_index) or (DB_iterator.position != last_DB.position):
            zero_one = DB_SEQ[DB_iterator.DB_index]
            if DB_iterator.is_start():
                exposure_time += (1-zero_one)*zero_S + zero_one*one_S
            else:
                exposure_time += (1-zero_one)*zero_E + zero_one*one_E
            DB_iterator.move_forward()
    else:
        exposure_time = 0.01
        
    return exposure_time


#def lut_display(image, display_min, display_max, bit_depth=14):
    #'''
    #Look up table for halftoning image
    #'''
    #lut = np.arange(2**16, dtype='uint16')
    #lut = display(lut, display_min, display_max)
    #return np.take(lut, image)

#======================Finding weighted aperture centers=======================#
def aperture_centers_weighted_start_differential(starts):
    # WEIGHTED START DIFFERENTIAL VERSION
    # BEST!

    # weighted centers of apertures
    positions = np.array([(0.7*starts[:-1]['x_image']+0.3*starts[1:]['x_image']),
                        (0.7*starts[:-1]['y_image']+0.3*starts[1:]['y_image'])]).transpose()

    # aperture sizes along the trajectory
    radii = np.sqrt((starts[:-1]['x_image']-starts[1:]['x_image'])**2 + 
                    (starts[:-1]['y_image']-starts[1:]['y_image'])**2)/2

    # not sure it's useful
    #time_deltas = Time(starts[:-1]['datetime']) - Time(starts[1:]['datetime'])

    # thetas: angles for asymetrical apertures definition (eg. rectangles)
    delta_X = starts[:-1]['x_image']-starts[1:]['x_image']
    delta_Y = starts[:-1]['y_image']-starts[1:]['y_image']
    thetas = np.arctan2(delta_Y, delta_X)
    
    return positions, radii, thetas

#===========================Finding apeture centers============================#
def aperture_centers(starts):

    # weighted centers of apertures
    positions = np.array([(starts['x_image']),(starts['y_image'])]).transpose()

    # aperture sizes along the trajectory
    radii = np.ones(len(starts)) * 2.0

    # thetas: angles for asymetrical apertures definition (eg. rectangles)
    delta_X = starts[:-1]['x_image']-starts[1:]['x_image']
    delta_Y = starts[:-1]['y_image']-starts[1:]['y_image']
    thetas = np.arctan2(delta_Y, delta_X)
    
    return positions, radii, thetas

#=======================Initialising the de Bruijn squence=====================#
def ini_sequence_starts(table, scidata, sat_adu=16383):
    '''
    return a subset of table that should be bearing photometric information
    '''
    logger = logging.getLogger()
    
    table['zero_or_one'] = [DB_SEQ[int(i)] for i in table['de_bruijn_sequence_element_index']]

    if table.meta['firmware'] == 'april2014':
        # see if encoding is inverted
        inverted, _, _ = is_inverted(scidata, table, sat_adu=sat_adu)
        if inverted:
            logger.warning('Fireball encoding is inverted')

        # starts variable = actual DB status change to open shutter, even if inverted
        if inverted:
            starts = table[table['dash_start_end'] == 'E']
        else:
            starts = table[table['dash_start_end'] == 'S']
        
        return starts, inverted
    
    elif table.meta['firmware'] == 'june2017':
        # no problem with this on june2017 firmware
        # every astrometric datapoint should lead to a photometric datapoint
        inverted = False
        # so return the entire table
        return table, inverted
    else:
        raise NotImplementedError(f'No idea what to do with firmware {table.meta["firmware"]}')

#==============Newer experimental encoding post june2017 frimware==============#
def symmetric_circular_aperture_photom(scidata, starts, inverted, jpeg='', plot=False, sat_adu = 16383):
    '''
    Method to use for june2017 firmware
    TODO unfinished, still very experimental
    '''
    logger = logging.getLogger()
    
    positions, radii, thetas = aperture_centers(starts)
        
    #===Asymmetric has error map and a co-ord shift for R & B bands(halving)===#
    #error_map = get_error_map(scidata)

    # get median pixel value inside the fireball bounding box
    median_bckgnd = np.median(get_bounding_box(scidata, starts))
    #logger.debug(median_bckgnd)


    if plot:
        plt.close()
        fig = plt.gcf()
        ax = fig.gca()
        plt.imshow(jpeg, cmap=plt.cm.gray)

    #=======In asymmetric, exp_time = nan and aperture_saturated = False=======#

    starts['aperture_sum'] = np.nan
    starts['bckgng_sub_measurements'] = np.nan
    starts['SNR'] = np.nan
    starts['signal'] = np.nan
    starts['exp_time'] = 0.01
    starts['aperture_sum_err'] = np.nan
    starts['aperture_sum_err_plus'] = np.nan
    starts['aperture_sum_err_minus'] = np.nan
    starts['aperture_saturated'] = np.nan

    #===============In asymmetric, shutter open time is calced=================#
    for el, pos, radius, theta in zip(starts, positions, radii, thetas):
        
        zero_order_light = scidata[int(pos[1])][int(pos[0])]
        if zero_order_light >= sat_adu:
            el['aperture_saturated'] = True
        else:
            el['aperture_saturated'] = False
        
        #============Below and above commentout swap for asymmetric============#
        
        #zero_order_saturation_level = float(zero_order_light)/sat_adu
        #=============Switched to RectangularAperture in asymmetric============#
        aperture = CircularAperture(pos, radius*2)
        
        #===In asymmetric, saturation levels assessed here, and mask applied===#
        # do photometry
        aperture_result = aperture_photometry(scidata, aperture)
        
        ## get saturated pixels once photutils 0.4 is available TODO FIXME
        
        el['aperture_sum'] = aperture_result['aperture_sum'].data
        #el['aperture_sum_err'] = aperture_result['aperture_sum_err'].data TODO FIXME
        
        # calculating SNR and error, should bckgng by bckgnd?
        el['bckgng_sub_measurements'] = el['aperture_sum'] -  aperture.area()*median_bckgnd
        el['SNR'] = el['bckgng_sub_measurements'] / np.sqrt(el['bckgng_sub_measurements'] + aperture.area()*median_bckgnd)
        el['aperture_sum_err_plus'] = 1/el['SNR']
        el['aperture_sum_err_minus'] = 1/el['SNR']
        if plot:
            aperture.plot(color='white')

    # plotting something, find out what??
    
    if plot:
        #min_x, min_y, max_x, max_y = get_bounds([p['x_image'] for p in points], pos2)
        margin = 100
        ax.set_xlim([np.min(starts['x_image'])-margin, np.max(starts['x_image'])+margin])
        ax.set_ylim([np.min(starts['y_image'])-margin, np.max(starts['y_image'])+margin])
        
        # saving figure to file
        full_path = starts.meta['self_file_name']
        dirname = os.path.dirname(full_path)
        basename = os.path.basename(full_path).split('.')[0]
        fname = os.path.join(dirname, basename + "_aperture_photometry.jpg")
        #plt.show()
        plt.savefig(fname, dpi=150)

    return starts
    
#========Used for older observations ("April 2014" de Bruijn encoding)=========#
def asymmetric_elongated_aperture_photom(scidata, starts, inverted, jpeg='', plot=False, sat_adu = 16383, spectral_band='G'):
    '''
    Compute photometry using rectangular apertures.
    This is a good method to use for "April 2014" de Bruijn encoding
    '''
    logger = logging.getLogger()
    
    logger.debug(f'science image is dimensions: {np.shape(scidata)}')
    
    match, factor = table_matches_image(starts, scidata)
    if not match:
        logger.warning(f'table data does not match image data. Dividing table coordinates by a factor of {factor:.2f} to table data (spectral band is {spectral_band})')
        starts['x_image'] /= factor
        starts['y_image'] /= factor
    
    positions, radii, thetas = aperture_centers_weighted_start_differential(starts)

    # get median pixel value inside the fireball bounding box
    median_bckgnd = np.median(get_bounding_box(scidata, starts))


    if plot:
        plt.close()
        fig = plt.gcf()
        ax = fig.gca()
        plt.imshow(jpeg, cmap=plt.cm.gray)

    starts['aperture_sum'] = np.nan
    starts['bckgng_sub_measurements'] = np.nan
    starts['SNR'] = np.nan
    starts['signal'] = np.nan
    starts['exp_time'] = np.nan
    starts['aperture_sum_err'] = np.nan
    starts['aperture_sum_err_plus'] = np.nan
    starts['aperture_sum_err_minus'] = np.nan
    starts['aperture_saturated'] = False

    for i in range(len(starts)-1):
        starts[i]['exp_time'] = LC_open_time(starts, i, i+1, inverted)

    for el, pos, radius, theta in zip(starts[:-1], positions, radii, thetas):
        
        #zero_order_light = scidata[int(pos[1])][int(pos[0])]
        
        aperture = RectangularAperture(pos, w=radius*2, h=7.*2, theta=theta) # 7
        
        # determine if any of the pixels in the aperture have reached saturation level,
        # flag accordingly
        aperture_mask = aperture.to_mask('center')
        aperture_values = aperture_mask.multiply(scidata)
        if np.any(aperture_values == sat_adu):
            el['aperture_saturated'] = True
        
        
        # do photometry
        aperture_result = aperture_photometry(scidata, aperture, method='subpixel', subpixels=16)
        
        
        el['aperture_sum'] = aperture_result['aperture_sum'].data
        #el['aperture_sum_err'] = aperture_result['aperture_sum_err'].data TODO FIXME
        
        el['bckgng_sub_measurements'] = el['aperture_sum'] -  aperture.area*median_bckgnd
        el['SNR'] = el['aperture_sum'] / np.sqrt(el['aperture_sum'] + aperture.area*median_bckgnd)
        el['aperture_sum_err_plus'] = 1/el['SNR']
        el['aperture_sum_err_minus'] = 1/el['SNR']
        if plot:
            aperture.plot(color='white')


    
    if plot:
        #min_x, min_y, max_x, max_y = get_bounds([p['x_image'] for p in points], pos2)
        ax.set_xlim([np.min(starts['x_image']), np.max(starts['x_image'])])
        ax.set_ylim([np.min(starts['y_image']), np.max(starts['y_image'])])
        
        full_path = starts.meta['self_file_name']
        dirname = os.path.dirname(full_path)
        basename = os.path.basename(full_path).split('.')[0]
        fname = os.path.join(dirname, basename + "_aperture_photometry_"+spectral_band+".jpg")
        plt.savefig(fname, dpi=150)

    return starts

#Adding different bands to the table, called in add_light_curves_to_observations#
def add_light_curve_band_to_observations(ppTable, starts, spectral_band):
    '''
    Adds photometry results to the table for a particular band
    parameters:
        - ppTable: point picking table (astropy Table)
        - starts: photometry data
        - spectral_band: which spectral band should be added
    '''

    
    # prepare the table: remove colliding columns, initialize then with NaN
    flux_col_ = str('flux_dash_' + spectral_band)
    luminosity_col_ = str('brightness_dash_' + spectral_band)
    saturation_col_ = str('saturated_flag_' + spectral_band)
    snr_col_ = str('SNR_' + spectral_band)
    exptime_col_ = str('aperture_exptime')
 
    #=============What if usable stuff exists in these columns??===============#   
    possible_newcols = [flux_col_, luminosity_col_, exptime_col_, saturation_col_, snr_col_]
    for c in possible_newcols:
        if c in ppTable.colnames:
            ppTable.remove_columns(c)
            
    #Setting initial values of columns to nan
    flux_col = np.ones(len(ppTable)) * np.nan
    luminosity_col = np.ones(len(ppTable)) * np.nan
    exptime_col = np.ones(len(ppTable)) * np.nan
    snr_col = np.ones(len(ppTable)) * np.nan
    saturation_col = np.ones(len(ppTable), dtype=np.bool_) * np.nan
    
    #This is a feature in the astropy package for dealing with data structures
    # add the columns, modify values later
    c1 = Column(name=flux_col_, data=flux_col*u.adu/u.second)
    c2 = Column(name=luminosity_col_, data=luminosity_col*u.adu)
    c3 = Column(name=exptime_col_, data=exptime_col*u.second)
    c4 = Column(name=saturation_col_, data=saturation_col)
    c5 = Column(name=snr_col_, data=snr_col)
    ppTable.add_columns([c2, c1, c3, c4, c5])
    
    
    for dash in starts:
        # store the flux value on the start entry
        # find the corresponding row in the table
        # assume there is only one entry with this time (take the first one)        
        #=======================Is this strictly true???=======================#    
        theRowIndex = np.argwhere(ppTable['datetime'] == dash['datetime'])[0][0]
        #=====================Are these returning numbers?=====================#
        
        ppTable[theRowIndex][luminosity_col_] = dash['bckgng_sub_measurements']
        ppTable[theRowIndex][flux_col_] = dash['bckgng_sub_measurements'] / dash['exp_time']
        ppTable[theRowIndex][exptime_col_] = dash['exp_time']
        ppTable[theRowIndex][saturation_col_] = dash['aperture_saturated']
        ppTable[theRowIndex][snr_col_] = dash['SNR']
        
        
#=====Actually adding to table, using add_light_curve_bands_to_observations====# 
def add_light_curves_to_observations(ppTable, photom_calcs, inverted):
    '''
    Adds photometry results to the table
    parameters:
        - ppTable: point picking table (astropy Table)
    '''
    
    for spectral_band in photom_calcs:
        add_light_curve_band_to_observations(ppTable, photom_calcs[spectral_band], spectral_band)
    
    # meta 
    #================What is the meta, is this the new meta??==================#
    now_timestamp = datetime.datetime.now()   #why timestamp?
    writetime = now_timestamp.strftime('%Y-%m-%dT%H:%M:%S')
    ppTable.meta['photometry_raw_write_time'] = writetime
    ppTable.meta['photometry_raw_software'] = SCRIPT_AND_VERSION
    ppTable.meta['photometry_raw_method'] = 'Rectangular dash aperture photometry'
    ppTable.meta[LC_SHUTTER_INVERSION_FLAG] = str(inverted)
    ppTable.meta['photometric_spectral_band'] = spectral_band
    
#====Performs calibration with background stars, normalises the flux I think===#
def photometric_calibration(table):
    '''
    Perform photometric calibration wrt stars
    Substract zero-point from calculationsparameters:
    parameters:
        - table: point picking table (astropy Table)
    '''
    logger = logging.getLogger()
    
    spectral_bands = [c.split('_')[-1] for c in table.colnames if c.startswith('flux_dash_')]
    
    catafile = os.path.join(os.path.dirname(table.meta['self_file_name']), 'calib', os.path.basename(table.meta['astrometry_catalog']))
    # Keys: zero-point, error, time, dark subtraction
    
    #spectral_band = table.meta['photometric_spectral_band']
    
    if not os.path.exists(catafile):
        logger.error(dfn_utils.warning + 'No reference star catalog found. Skipping zero-point photometric calibration.')
        raise FileNotFoundError()
    
    calibration_header = cataphot.main(catafile)
    
    zero_point = float(str(calibration_header['ZERO-PT']))
    table.meta['photometric_zero_point'] = zero_point
    table.meta['photometric_zero_point_uncertainty'] = float(str(calibration_header['ZPT_ERR']))
    table.meta['photometric_zero_point_datetime'] = str(calibration_header['DATE-OBS'])
    
    for spectral_band in spectral_bands:
        table['m_' + spectral_band] = zero_point - 2.5 * np.log10(table['flux_dash_'+spectral_band] * SPECTRAL_BAND_COLLECTING_FRACTION_INVERSE[spectral_band] / (1 * u.adu/u.second))
    
    return table

#========================Calculating the colour indices========================#
def colour_indices(table):
    '''
    calculate colour indices
    '''
    logger = logging.getLogger()

    for colour_pair in [('B','G'),('G','R'),('B','R')]:
        c1, c2 = colour_pair
        logger.info(f'Calculating {c1}-{c2} colour index')
        c1_mag_name = 'm_' + c1
        c2_mag_name = 'm_' + c2
        c1_sat_name = 'saturated_flag_' + c1
        c2_sat_name = 'saturated_flag_' + c2
        
        # make sure all of required columns are there
        for col in [c1_mag_name, c2_mag_name, c1_sat_name, c2_sat_name]:
            if col not in table.colnames:
                logger.error(f'No column named {c1_mag_name}, therefore connot calculate {c1}-{c2} colour index')
                raise KeyError(col)
        
        
        #Naming column and setting initial value
        colour_index_colname = f'colour_index_{c1}-{c2}'
        table[colour_index_colname] = np.nan
        
        # create a mask that tells us if any of the 2 records used are saturated
        any_saturation = (table[c1_sat_name].astype(np.bool_) | table[c2_sat_name].astype(np.bool_))
        
        # calculate colour index for non saturated records
        table[colour_index_colname][~any_saturation] = table[c1_mag_name][~any_saturation] - table[c2_mag_name][~any_saturation]

#==============================Plotting Functions==============================#

#Plotting the colour indices, can be modified for all graphs, not just one (B-R)#
def plot_colour_index(t, colour_index, outfname):
    '''
    plot light curve
    parameters:
        t: astropy table with the calculated photometry
        index: colour index (eg. "B-R")
    '''
    logger = logging.getLogger()
    
    table = deepcopy(t)
    
    # calulate relative time from first record
    ref_time = Time(table['datetime'][0])
    table['reltime'] = (Time(table['datetime']) - ref_time).sec


    plt.close()
    plt.figure()
    ax = plt.axes()
    
    colname = 'colour_index_' + colour_index

    ax.scatter(table['reltime'], table[colname],
            color='black')

    plt.xlabel(f'time (s) after {ref_time.isot}')
    plt.ylabel(f'{colour_index} colour index')
    plt.grid()
    #plt.legend()
    
    logger.info(f'saving colour index plot to {outfname}')
    plt.savefig(outfname, dpi=100)
    
#Plotting the apparent magnitude of the different filters in the light curve (R,G,B)#
def plot_apparent_mag(t, outfname):
    '''
    plot light curve
    parameters:
        t: astropy table with the calculated photometry
    '''
    logger = logging.getLogger()
    
    table = deepcopy(t)
    
    # calulate relative time from first record
    ref_time = Time(table['datetime'][0])
    table['reltime'] = (Time(table['datetime']) - ref_time).sec


    plt.close()
    plt.figure()
    ax = plt.axes()

    for spectral_band in SPECTRAL_BAND_DE_BAYER:
        app_mag_colname = 'm_' + spectral_band
        saturation_colname = 'saturated_flag_' + spectral_band
        if not app_mag_colname in table.colnames:
            logger.warning(f'{spectral_band} not available for plotting')
            continue
        
        sat_mask = (table[saturation_colname] == 1)
        n_tot = np.count_nonzero(~np.isnan(table[saturation_colname]))
        n_sat = np.nansum(table[saturation_colname][sat_mask])
        n_non_sat = n_tot - n_sat
        
        # plot unsaturated and saturated records with different flags
        ax.scatter(table['reltime'][~sat_mask], table[app_mag_colname][~sat_mask],
                label=f'{spectral_band} n={n_non_sat:.0f}',
                color=spectral_band.lower())
        ax.scatter(table['reltime'][sat_mask], table[app_mag_colname][sat_mask],
                label=f'{spectral_band} saturated n={n_sat:.0f}',
                marker='x', color=spectral_band.lower())

    plt.gca().invert_yaxis()
    plt.xlabel(f'time (s) after {ref_time.isot}')
    plt.ylabel('apparent magnitude')
    plt.grid()
    plt.legend()
    
    logger.info(f'saving light curve plot to {outfname}')
    plt.savefig(outfname, dpi=100)
    
#========================Main function for calculations========================#
def main(raw_hdu, table, aperture_plot=False, jpeg=None):
    '''
    main method for light curve calculations
    parameters:
        - raw_hdu: FITS HDU to work on
        - table: point picking table (astropy Table)
    '''
    logger = logging.getLogger()
    
    # initialize sequence starts
    starts, inverted = ini_sequence_starts(table, np.flipud(raw_hdu.data))
    
    photom_calcs = {}
    for spectral_band in ['R','G','B']:
        logger.info(f'calculating light curve in the {spectral_band} band')
        de_bayer_code = SPECTRAL_BAND_DE_BAYER[spectral_band]
        # de Bayer
        de_bayer_fits_object = fitsConversion.de_Bayer(raw_hdu, color_code=de_bayer_code)
        
        #  flip data array to match pixel coordinates (instead of FITS)
        scidata = np.flipud(de_bayer_fits_object.data)
        
        
        # do photometry
        if table.meta['firmware'] == 'april2014':
            starts = asymmetric_elongated_aperture_photom(scidata, starts, inverted, jpeg=jpeg, plot=aperture_plot, spectral_band=spectral_band)
        else:
            raise NotImplementedError()
            starts = symmetric_circular_aperture_photom(scidata, starts, inverted, jpeg=jpeg, plot=aperture_plot, spectral_band=spectral_band)
            
        photom_calcs[spectral_band] = deepcopy(starts)
    
    # save results to main table
    add_light_curves_to_observations(table, photom_calcs, inverted)
    
    # do photometric calibration at a later step
    try:
        photometric_calibration(table)
    except KeyError as e:
        logger.error(e)
        logger.warning("Cannot perform photometric_calibration, skipping")
        pass
    
    # calculate colour indices
    try:
        colour_indices(table)
    except KeyError as e:
        logger.error(e)
        logger.warning("Cannot calculate colour indices, skipping")
        pass
    
    return table

#===============Cleaning up previous calculations in dictionaries==============#
def clean_up_old_photometry_calculations(table):
    '''
    Clean up old photometry information calculated by previous run
    '''
    logger = logging.getLogger()
    
    cols_to_remove = [c for c in table.colnames if (c.startswith('flux_dash_')
                                                    or c.startswith('brightness_dash_')
                                                    or c.startswith('saturated_flag_')
                                                    or c.startswith('SNR_')
                                                    or (c.startswith('m_') and len(c)==3))]
    if len(cols_to_remove) > 0 :
        table.remove_columns(cols_to_remove)
        logger.info(f'Removed old columns: {cols_to_remove}')
    
    keys_to_remove = [k for k in table.meta if k.startswith('photometry_raw_')]
    for k in keys_to_remove:
        table.meta.pop(k)
        logger.info(f'Removed old key {k}')
    
    

#==================Main function that executes calculations====================#
def main_disk_io(ifile, plot=False, destination=None):
    '''
    main method for light curve calculations
    parameters:
        - ifile: input point picking file
    '''
    #Logs errors and prints to terminal?
    logger = logging.getLogger()
    
    #Initialising file path
    wdir = os.path.dirname(ifile)
    
    #Reading in table from file
    table = Table.read(ifile, format='ascii.ecsv', guess=False, delimiter=',')
    table.meta['self_disk_path'] = ifile
    
    #Raising error if no table exists in file
    if not np.all(dfn_utils.has_reliable_timing(table)):
        raise dfn_utils.WrongTableTypeException('Need timed point picking to do photometry')
    
    #Cleans up old photometry calcs if previously run
    clean_up_old_photometry_calculations(table)
    
    # determine micro-controller firmware version,
    # hence what LC shutter de Bruijn encoding was used for this fireball
    try:
        fw = timing_corrections.determine_firmware_version(table, table_loc='camera')
        if not fw in FIRMWARES_SUPPORTED:
            raise dfn_utils.WrongTableTypeException('Encoding firmware {} is currently not supported for photometry calculations.'.format(fw))
        table.meta['firmware'] = fw
    except timing_corrections.UnknownEncodingError:
        raise dfn_utils.WrongTableTypeException('Cannot determine encoding firmware.')



    # propagate self name for figure saving
    table.meta['self_file_name'] = ifile
    
    image_file_meta = table.meta['image_file']
    # figure out what the RAW file should be called`
    nef_file = os.path.join(os.path.dirname(ifile),
                            "-".join(image_file_meta.split('-')[:-1]) + '.NEF')
        
    # try to locate the RAW image file
    if not os.path.isfile(nef_file):
        logger.error('cannot find expected RAW file {}'.format(nef_file))
        raise FileNotFoundError()
    
    # convert the RAW file to a usable data structure
    raw_hdu = fitsConversion.raw2fits(nef_file)
    
    # jpeg image for visualization
    try:
        jfile = os.path.join(wdir, os.path.splitext(image_file_meta)[0] + '.jpeg')
        jpeg = imread(jfile, as_gray=False, pilmode='P')
        aperture_plot = True
    except FileNotFoundError:
        aperture_plot = False
        jpeg = None
        logger.warning('No jpeg can be found, not plotting apertures')
        # if not jpeg present, forget about plotting apertures
    
    # main processing
    main(raw_hdu, table, aperture_plot=aperture_plot, jpeg=jpeg)
    
    # add MJD column (easier to plot)
    table['MJD'] = Time(table['datetime']).mjd
    
    # plot
    if plot:
        if destination:
            base = os.path.join(destination, os.path.basename(ifile).split('.')[0])
        else:
            base = ifile.split('.')[0]
            
        outfname = base + '_light_curve.png'
        plot_apparent_mag(table, outfname)
        
        c_index = 'B-R'
        outfname = base + '_colour_index_' + c_index + '.png'
        plot_colour_index(table, c_index, outfname)
        
    # if a different destination is specified
    if destination:
        ofile = os.path.join(destination, os.path.basename(ifile))
    else:
        # else overwrite the input file
        ofile = ifile
    # save results    
    table.write(ofile, format='ascii.ecsv', delimiter=',', overwrite=True)

#=========Returns the file paths for all .ecsv files in the directory==========#
'''
def multi_curves(current_dir):
    #Initialising list files and paths
    file_paths = []

    #Cycling through events in the photometry directory
    with os.scandir(current_dir) as enteries:
        for entry in enteries:
            #Cycling through cameras that observed each event
            with os.scandir(entry) as obs:
                for ob in obs:
                    current_cam = os.path.join(ob)
                    #Finding the .ecsv files for each observation
                    for f_name in os.listdir(current_cam):
                        if f_name.endswith('.ecsv'):
                            obs_file = f_name
                            input_file = os.path.join(current_cam, obs_file)
                            file_paths.append(input_file)
                            
    return(file_paths)
''' 
#========================Test function for file location=======================#

def test():
    '''
    bodge test method
    '''
    
    wdir='/home/data/DN150417_01/43_Forrest/'
    pp_file = '43_2015-04-17_200359_DSC_1270-G_DN150417_01_2016-04-13_160324_hadry_nocomment.ecsv'
    ifile = os.path.join(wdir,pp_file)
    main_disk_io(ifile, plot=True)
    
#===================================Code to Compute============================#
if __name__ == '__main__':
    #test_fireball_light_curve()
    #test()
    #exit(1)
    
    # stdout logger
    logger = logging.getLogger()
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s, %(levelname)s, %(module)s, %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.setLevel(logging.INFO)
    
    
    import argparse
    
    parser = argparse.ArgumentParser(description='Compute meteor photometry')
    parser.add_argument('--version', action='version', version=SCRIPT_AND_VERSION)
    parser.add_argument("-p", "--pointsfile", type=str, help="Point picking file", required=True)
    parser.add_argument("-d", "--destination", type=str, help="Different destination folder for output", default=None)
    
    args = parser.parse_args()
    
    pp_file = args.pointsfile
    destination = args.destination
    
    main_disk_io(pp_file, plot=True, destination=destination)
    

#==============================================================================#
#================================END OF FILE===================================#
#==============================================================================#
