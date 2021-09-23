#!/usr/bin/env python

# general imports
import os
import argparse
try:
    from configparser import SafeConfigParser
except ImportError:
    from ConfigParser import SafeConfigParser

# sci packages imports
import numpy as np
from astropy.io import fits
from astropy.table import Table, hstack
#from astropy import units as u
#from astropy.time import Time, TimeDelta
#from astropy.coordinates import EarthLocation
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry
import scipy.spatial
import astropy.units as u

# plotting imports
import matplotlib.pyplot as plt
import seaborn as sns

# local imports
import global_astrometry_error_control as gaec
#import dfn_utils


def get_true_exptime(header):
    '''
    Compute the true exposure time of stiff objects in DFN exposure
    Accounts for the LC shutter operation.
    header: FITS header handle of the image
    '''
    
    exposure_time_main_shutter = float(header['EXPTIME'])
    print('Main shutter exposure time: ' + str(exposure_time_main_shutter))
    
    # TODO parse interval logs to get micro-controller parameters
    number_of_ones = 150
    number_of_zeros = 100
    n_dashes = number_of_ones + number_of_zeros
    one_LC_open_fraction = 0.6
    zero_LC_open_fraction = 0.2
    total_LC_open_fraction = (number_of_ones*one_LC_open_fraction + number_of_zeros*zero_LC_open_fraction)/n_dashes
    print('total_LC_open_fraction: ' + str(total_LC_open_fraction))

    exposure_time_true = exposure_time_main_shutter * total_LC_open_fraction
    print('exposure_time_true: ' + str(exposure_time_true))

    return exposure_time_true


def photometry(refstars, hdu, plot=False):
    '''
    Perform aperture photometry
    refstars: source extraction with astrometric information
    hdu: science image fits handle
    '''

    scidata = hdu.data
    header = hdu.header

    exposure_time_true = get_true_exptime(header)

    positions = np.array(refstars['XC','YC']).tolist()
    apertures = CircularAperture(positions, r=3.5) #3.5
    annulus_apertures = CircularAnnulus(positions, r_in=4., r_out=6.)

    # do aperture photometry on the small circles
    rawflux_table = aperture_photometry(scidata, apertures)
    # do aperture photometry on the large circles
    bkgflux_table = aperture_photometry(scidata, annulus_apertures)
    raw_phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])
    phot_table = hstack([raw_phot_table,refstars], join_type='exact')
    bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area
    bkg_sum = bkg_mean * apertures.area
    # subtract background
    phot_table['residual_aperture_sum'] = phot_table['aperture_sum_raw'] - bkg_sum
    #print phot_table['residual_aperture_sum']

    # first order extinction coefficient in V band
    k_V = 0.2

    # -flux   + mag + mag due to atmos extinction
    # http://www.icq.eps.harvard.edu/ICQExtinct.html
    phot_table['measured_v_mag'] = - 2.5 * np.log10(phot_table['residual_aperture_sum'] / exposure_time_true) 
    phot_table['airmass_extinction'] = k_V / np.cos((90*u.deg - phot_table[gaec.alt_cata_col_]))
    
    phot_table['cata_v_mag'] = phot_table['V_MAG'] + phot_table['airmass_extinction']
    
    phot_table['instrumental_mag'] = phot_table['cata_v_mag'] - phot_table['measured_v_mag']
    
    #from glue import qglue
    #app = qglue(phot_table=phot_table)
    
    if plot:
        plt.close()
        #plt.scatter(phot_table['measured_v_mag'], phot_table['cata_v_mag'], color='r', label='')
        plt.scatter(phot_table['instrumental_mag'] , phot_table['cata_v_mag'] , color='r', label='')
        plt.legend()
        plt.xlabel('visual mag ')
        plt.ylabel('cat mag + atmos extinction')
        plt.show()
    
    
    consolidated_table = phot_table[~np.isnan(phot_table['instrumental_mag'])]
    
    return consolidated_table
    

def cut_off_outlier_stars(cata):
    '''
    Remove dubious reference stars (too far from global fit) and eliminate catalog artifact (multiple star systems)
    cata: input catalog
    '''
     
    # Empirical formula for outlier removal
    good_cata = cata[((cata[gaec.GLOBAL_FIT_DISTANCE] < 220) |
                            (cata['alt_cata'].data < -0.025 * cata[gaec.GLOBAL_FIT_DISTANCE] + 23)) &
                            (cata['V_MAG'] > 1.) & (cata['V_MAG'] < 5.)]
    
    to_remove = set([])
    referenceTree = scipy.spatial.cKDTree(np.column_stack((good_cata['XC'], good_cata['YC'])),
                                          leafsize=100)
    for s in good_cata:
        # get the 2 closest neighbours (including self)
        neighs = referenceTree.query((s['XC'], s['YC']),
                                     k=10,
                                     distance_upper_bound=30)
        # test if there is actually a close neighbour (first result is always self, second is np.inf if there is no close neighbour)
        if not np.isinf(neighs[0][1]):
            to_remove.add(neighs[1][0])
        for i in range(1, 9):
            if not np.isinf(neighs[0][i]):
                to_remove.add(neighs[1][i])
            else:
                break
    # delete the flagged stars
    good_cata.remove_rows(list(to_remove))
    
    return good_cata

def quality_control_plots(cata, cut_cata, overwrite=False):
    '''
    Save plots that show various aspects of the reference catalog
    cata: main reference catalog
    cut_cata: relevant subset of main catalog (eg. outlier cleaned)
    '''
    
    catalog_full_path = cata.meta['self_file_name']
    
    dirname = os.path.dirname(catalog_full_path)
    basename = os.path.basename(catalog_full_path).split('.')[0]
    
    # GENERAL ALTITUDE vs FIT QUALITY PLOT
    fname = os.path.join(dirname, basename + "_instrumental_mag_vs_altitude.png")
    if overwrite or not os.path.isfile(fname):
        plt.close()
        plt.scatter(cata['instrumental_mag'], cata[gaec.alt_cata_col_], color='r', label='all stars')
        plt.scatter(cut_cata['instrumental_mag'], cut_cata[gaec.alt_cata_col_], color='b', label='non dubious stars')
        plt.legend()
        plt.xlabel('Instrumental magnitude ')
        plt.ylabel('Altidude Catalog (degree)')
        plt.title('Used stars: ' + str(int(float(len(cut_cata))/len(cata)*100)) + '%')
        plt.savefig(fname, dpi=150)
    




def calculate_global_instrumental_magnitude(refcatalog, station, hdu, overwrite = False):
    '''
    Main function for ...
    '''
    
    cata = gaec.preprocess_star_cat(refcatalog, station)
    
    photom_table = photometry(cata, hdu)
    
    cut_consolidated_photom_table = cut_off_outlier_stars(photom_table)
    
    #phot_table['visual_mag'], phot_table['cata_v_mag']
    
    #A = np.vstack([cut_consolidated_photom_table['measured_v_mag'], np.ma.asarray(np.ones(len(cut_consolidated_photom_table['measured_v_mag'])))]).T
    #linreg = np.linalg.lstsq(A, cut_consolidated_photom_table['cata_v_mag'])
    
    fit = np.polyfit(cut_consolidated_photom_table['measured_v_mag'],
               cut_consolidated_photom_table['cata_v_mag'],
               1,
               full=True)
    
    #m, c = linreg[0]
    
    coeffs, residuals = fit[0], fit[1]
    
    print('slope: ' + str(coeffs[0]))
    print('offset: ' + str(coeffs[1]))
    print('residuals: ' + str(residuals))
    
    
    i_m_mean = np.mean(cut_consolidated_photom_table['instrumental_mag'])
    i_m_std = np.std(cut_consolidated_photom_table['instrumental_mag'])
    i_m_median = np.median(cut_consolidated_photom_table['instrumental_mag'])
    
    catalog_photometric_efficiency = int(100*float(len(cut_consolidated_photom_table))/len(photom_table))
    
    print('mean instrumental_mag: ' + str(i_m_mean))
    print('stdev instrumental_mag: ' + str(i_m_std))
    print('median instrumental_mag: ' + str(i_m_median))
    print('catalog_photometric_efficiency: ' + str(catalog_photometric_efficiency) + ' %')
    

    hdu.header.set('ZERO-PT', i_m_mean, 'Photometric magnitude zero-point')
    hdu.header.set('ZPT_ERR', i_m_std, '1-SIG ERR photometric magnitude zero-point')
    hdu.header.set('PHOT_EF', catalog_photometric_efficiency, 'Photometric catalog efficiency (%)')
    # (normalised to 1 second of actual exposure time)
    
    quality_control_plots(photom_table, cut_consolidated_photom_table, overwrite = overwrite)
    
    
    
    return cut_consolidated_photom_table, i_m_mean, i_m_std
    
    
def main(refcatalog_file, overwrite = False):
    '''
    High order main that deals with disk I/O
    '''
    
    calib_dirname = os.path.dirname(refcatalog_file)
    fbase = os.path.splitext(os.path.basename(refcatalog_file))[0]
    
    observerCfgFile = os.path.join(calib_dirname, fbase + '_observer.cfg')
    print('observer config file: ' + observerCfgFile)
    
    imfile = os.path.join(calib_dirname, fbase + '.fits')
    print('Image file: ' + imfile)

    
    hdulist = fits.open(imfile, mode='update')
    hdu = hdulist[0]
    
    if not overwrite:
        try:
            # check that the right key words are present
            hdu.header['ZERO-PT']
            hdu.header['ZPT_ERR']
            hdu.header['PHOT_EF']
            
            hdulist.close()
            return hdu.header
        
        except KeyError:
            pass
    
    
    catalog = Table.read(refcatalog_file, format='votable')
    
    # Give the table its own file name as metadata
    catalog.meta['self_file_name'] = refcatalog_file

    station, _ = gaec.parse_observer_cfg_file(observerCfgFile)
    
    cut_consolidated_photom_table, i_m_mean, i_m_std = calculate_global_instrumental_magnitude(catalog, station, hdu, overwrite = overwrite)
    
    
    hdulist.flush()
    hdulist.close()
    
    return hdu.header

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Does photometry extraction on star catalog')
    parser.add_argument("-i", "--inputcatalog", type=str, required=True, help="Input star catalog")

    args = parser.parse_args()
    
    main(args.inputcatalog, overwrite = True)

    

