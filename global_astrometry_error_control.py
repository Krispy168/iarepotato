#!/usr/bin/env python2.7

# general imports
import os
import argparse
import itertools
import warnings

try:
    from configparser import SafeConfigParser
except ImportError:
    from ConfigParser import SafeConfigParser

# sci packages imports
import numpy as np
import astropy.io.votable as votable
from astropy.table import Table, Column, hstack
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation, AltAz, ICRS
import scipy.spatial

# plotting imports
import matplotlib.pyplot as plt
import seaborn as sns

# local imports
import dfn_utils

from astropy.io.fits.verify import VerifyWarning
warnings.simplefilter('ignore', category=VerifyWarning)


alt_cata_col_ = 'alt_cata'
az_cata_col_ = 'az_cata'
alt_fitted_col_ = 'alt_fitted'
az_fitted_col_ = 'az_fitted'
delta_az_cata_col_ = 'delta_az_cata'
delta_alt_cata_col_ = 'delta_alt_cata'
GLOBAL_FIT_DISTANCE = 'DIST_arcsec'


def parse_observer_cfg_file(obs_cfg_file):
    config = SafeConfigParser()
    config.read(obs_cfg_file)
    
    firstsection = 'Observer'
    lat = float(config.get(firstsection,'latitude'))
    lon = float(config.get(firstsection,'longitude'))
    alt = float(config.get(firstsection,'altitude'))
    
    time_str = str(config.get(firstsection,'isodate_mid_obs'))
    time_obs = Time(time_str)
    
    #isodate_mid_obs = 2015-12-02T14:28:40.500
    
    station = EarthLocation(lat=lat*u.deg,
                            lon=lon*u.deg,
                            height=alt*u.meter)

    return station, time_obs

def convert_radec_file_to_altaz(input_file,
                                obs_cfg_file,
                                output_file,
                                alt_key='altitude', az_key='azimuth',
                                ra_key='ra', dec_key='dec'):
    '''
    Convert RA DEC file to Alt-Az
    input_file: File containing equatorial coordinates
    obs_cfg_file: Oberver config file
    output_file: file to save
    keys: column header keywords
    '''
    
    if ra_key.lower() == 'ha':
        ra_unit = u.hour
    else:
        ra_unit = u.deg
    
    input_table = Table.read(input_file, format='ascii.csv', guess=False, delimiter=',')
    
    station, obstime = parse_observer_cfg_file(obs_cfg_file)
    
    equatorial_data = ICRS(ra=input_table[ra_key]*ra_unit, dec=input_table[dec_key]*u.deg)
    horizontal_data = equatorial_data.transform_to(AltAz(obstime=obstime,location=station))
    
    alt_col = Column(data=horizontal_data.alt, name=alt_key)
    az_col = Column(data=horizontal_data.az, name=az_key)
    
    output_table = Table([az_col, alt_col])
    
    output_table.write(output_file, format='ascii.csv', delimiter=',', overwrite=True)
    

def preprocess_star_cat(input_cata, station):
    '''
    Makes sure the star catalog contains all the required information: index alt-az, calculated alt-az, deltas...
    input_cata: astropy Table or votable path
    observer_station: EarthLocation object
    '''
    
    if isinstance(input_cata, Table):
        input_table = input_cata
    else:
        # Read the data
        input_table = Table.read(input_cata, format='votable')
        


    # try to remove alt-az columns if they exist
    possible_newcols = [alt_cata_col_,
                        az_cata_col_,
                        delta_az_cata_col_,
                        delta_alt_cata_col_,
                        alt_fitted_col_,
                        delta_alt_cata_col_]
    
    for c in possible_newcols:
        if c in input_table.colnames:
            input_table.remove_columns(c)
            
    
    times = Time(Time(input_table['jd_mid_obs'], format='jd', scale='utc').isot)
    
    cata_data = ICRS(ra=input_table['RA_cat_deg']*u.deg, dec=input_table['DEC_cat_deg']*u.deg)
    cata_data_altaz = cata_data.transform_to(AltAz(obstime=times,location=station))
    
    fitted_data = ICRS(ra=input_table['RA_deg']*u.deg, dec=input_table['DEC_deg']*u.deg)
    fitted_data_altaz = fitted_data.transform_to(AltAz(obstime=times,location=station))
    
    alt_cata_col = Column(data=cata_data_altaz.alt, name=alt_cata_col_)
    az_cata_col = Column(data=cata_data_altaz.az, name=az_cata_col_)
    alt_fitted_col = Column(data=fitted_data_altaz.alt, name=alt_fitted_col_)
    az_fitted_col = Column(data=fitted_data_altaz.az, name=az_fitted_col_)
    delta_alt_cata_col = Column(data=fitted_data_altaz.alt-cata_data_altaz.alt, name=delta_alt_cata_col_)
    delta_az_cata_col = Column(data=fitted_data_altaz.az-cata_data_altaz.az, name=delta_az_cata_col_)
    
    new_cols = [alt_cata_col,
                az_cata_col,
                delta_az_cata_col,
                delta_alt_cata_col,
                alt_fitted_col,
                az_fitted_col]
    
    
    input_table.add_columns(new_cols)
    
    
    return input_table


def cut_off_outlier_stars(cata):
    '''
    Remove dubious reference stars (too far from global fit) and eliminate catalog artifact (multiple star systems)
    '''
    
    # Empirical formula for outlier removal
    # + formula for removing stars that are within 5 degrees of the coordinate system pole
    good_cata = cata[((cata[GLOBAL_FIT_DISTANCE] < 220) |
                            (cata['alt_cata'] < -0.025 * cata[GLOBAL_FIT_DISTANCE] + 23)) &
                            (cata['alt_cata'] < 85.)]

    to_remove = set([])
    referenceTree = scipy.spatial.cKDTree(np.column_stack((good_cata['XC'], good_cata['YC'])), leafsize=100)
    for s in good_cata:
        # get the 2 closest neighbours (including self)
        neighs = referenceTree.query((s['XC'], s['YC']), k=2, distance_upper_bound=20)
        # test if there is actually a close neighbour (first result is always self, second is np.inf if there is no close neighbour)
        if not np.isinf(neighs[0][1]):
            # find out which one of the 2 is the worst (worst piastro fit), and flag it for deletion
            if good_cata[GLOBAL_FIT_DISTANCE][neighs[1][1]] < good_cata[GLOBAL_FIT_DISTANCE][neighs[1][0]]:
                to_remove.add(neighs[1][0])
            else:
                to_remove.add(neighs[1][1])
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
    fname = os.path.join(dirname, basename + "_global_fit_alt_vs_fitdistance.png")
    print('Attempting to generate astrometric GENERAL ALTITUDE vs FIT QUALITY PLOT: ' + fname)
    if overwrite or not os.path.isfile(fname):
        plt.close()
        plt.scatter(cata[GLOBAL_FIT_DISTANCE], cata[alt_cata_col_], color='r', label='all stars')
        plt.scatter(cut_cata[GLOBAL_FIT_DISTANCE], cut_cata[alt_cata_col_], color='b', label='non dubious stars')
        plt.legend(frameon=True, loc='best', fancybox=True, framealpha=0.5)
        plt.xlabel('Distance from global fit (arcsec)')
        plt.ylabel('Altidude Catalog (degree)')
        plt.title('Used stars: ' + str(int(float(len(cut_cata))/len(cata)*100)) + '%')
        plt.savefig(fname, dpi=150)
    
    # FITTED - CALCULATED plot on both axis
    fname = os.path.join(dirname, basename + "_correlated_horizon_errors.png")
    if overwrite or not os.path.isfile(fname):
        plt.close()
        plt.scatter(3600.*cut_cata[delta_az_cata_col_], 3600.*cut_cata[delta_alt_cata_col_], c=cut_cata[alt_cata_col_], cmap=plt.cm.coolwarm)
        #plt.legend()
        plt.xlabel('Azimuth: Fitted - Catalog (arcsec)')
        plt.ylabel('Altitude: Fitted - Catalog (arcsec)')
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Altidude Catalog (degree)')
        plt.savefig(fname, dpi=150)
    
    # AZIMUTH RESIDUALS vs AZIMUTH
    fname = os.path.join(dirname, basename + "_azimuth_errors.png")
    if overwrite or not os.path.isfile(fname):
        plt.close()
        plt.scatter(3600.*cut_cata[delta_az_cata_col_], cut_cata[az_fitted_col_], color='b', label='non dubious stars')
        plt.legend(frameon=True, loc='best', fancybox=True, framealpha=0.5)
        plt.xlabel('Azimuth: Fitted - Catalog (arcsec)')
        plt.ylabel('Azimuth fitted (degree)')
        plt.title('Used stars: ' + str(int(float(len(cut_cata))/len(cata)*100)) + '%')
        plt.savefig(fname, dpi=150)
    
    
    # ALTITUDE RESIDUALS vs ALTITUDE
    fname = os.path.join(dirname, basename + "_altitude_errors.png")
    print('Attempting to generate astrometric ALTITUDE RESIDUALS vs ALTITUDE: ' + fname)
    if overwrite or not os.path.isfile(fname):
        alt_residuals_fit = np.polyfit(cut_cata[alt_cata_col_], 3600.*cut_cata[delta_alt_cata_col_], 3)
        alt_residuals_poly = np.poly1d(alt_residuals_fit)
        plt.close()
        plt.scatter( 3600.*cut_cata[delta_alt_cata_col_], cut_cata[alt_fitted_col_], color='b', label='non dubious stars')
        yp = np.linspace(0, 90, 100)
        plt.plot(alt_residuals_poly(yp), yp, '-', color='r', label='fit')
        plt.legend(frameon=True, loc='best', fancybox=True, framealpha=0.5)
        plt.xlabel('Altitude: Fitted - Catalog (arcsec)')
        plt.ylabel('Altidude fitted (degree)')
        plt.title('Used stars: ' + str(int(float(len(cut_cata))/len(cata)*100)) + '%')
        plt.savefig(fname, dpi=150)


def correct_altitudes(firetable, cut_cata, station):
    '''
    correct altitude information based on correlated residuals on the global astrometric fit
    firetable: table to correct
    cut_cata: meaningful reference catalog
    station: Observer location
    '''
    
    # TODO FIXME calculate the pixscale
    pixscale = (120 * u.arcsec).to(u.deg).value
    # Default error
    default_error = np.nan
    default_error = 0.1
    
    # fit
    alt_residuals_fit = np.polyfit(cut_cata[alt_cata_col_], cut_cata[delta_alt_cata_col_], 3)
    alt_residuals_poly = np.poly1d(alt_residuals_fit)
    
    #firetable.show_in_browser(jsviewer=True)
    
    firetable['altitude'] = firetable['altitude'] - alt_residuals_poly(firetable['altitude'])
    
    cut_cata['corrected_delta_altitude'] = cut_cata[delta_alt_cata_col_] - alt_residuals_poly(cut_cata[alt_fitted_col_])
    
    times = Time(Time(cut_cata['jd_mid_obs'], format='jd', scale='utc').isot)
    #print(cut_cata[az_cata_col_])
    reference_altaz = AltAz(az=cut_cata[az_cata_col_].data*u.deg,alt=cut_cata[alt_cata_col_].data*u.deg,
                        location=station,
                        obstime=times)

    
    
    for i in range(len(firetable)):
        # choose reference data @ +-30 degrees
        cut = 30.*u.deg
        firepoint = firetable[i]
        firepoint_coord = AltAz(az=firepoint['azimuth']*u.deg, alt=firepoint['altitude']*u.deg,
                        location=station,
                        obstime=Time(firepoint['datetime']))
        ref_subset_az = cut_cata[firepoint_coord.separation(reference_altaz) < cut]
        
        if len(ref_subset_az) < 5:
            firepoint['err_minus_altitude'] = default_error
            firepoint['err_plus_altitude'] = default_error
            firepoint['err_minus_azimuth'] = default_error
            firepoint['err_plus_azimuth'] = default_error
        else:
            
            err_alt = np.std(ref_subset_az['corrected_delta_altitude']) + pixscale * firepoint['err_minus_x_image']
            firepoint['err_minus_altitude'] = err_alt
            firepoint['err_plus_altitude'] = err_alt
            
            err_az = np.std(ref_subset_az[delta_az_cata_col_]) + pixscale * firepoint['err_minus_x_image']
            firepoint['err_minus_azimuth'] = err_az
            firepoint['err_plus_azimuth'] = err_az
            
    
    #firetable.show_in_browser(jsviewer=True)

    
    return firetable

    


def do_corrections_and_compute_errors(refcatalog, station, firetable):
    '''
    Main function for global calibration post-processing
    '''
    
    cata = preprocess_star_cat(refcatalog, station)
    
    cut_cata = cut_off_outlier_stars(cata)
    
    quality_control_plots(cata, cut_cata)
    
    firetable = dfn_utils.add_JD_info(firetable)
    
    correct_altitudes(firetable, cut_cata, station)
    
    
    firetable.meta['astrometry_number_stars'] = len(cut_cata)
    firetable.meta['astrometry_catalog_efficiency'] = int(100*len(cut_cata)/len(cata))
    
    return cut_cata
    
    

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Corrects correlated errors from the global calibration method')
    #parser.add_argument("-r", "--referencefile", type=str, required=True, help="Votable catalog of reference stars")
    #parser.add_argument("-c", "--cfgfile", type=str, required=True, help="Observer config file corresponding to the reference catalog")
    #parser.add_argument("-m", "--multitemporal", action="store_true", default=True, help="Use if catalog is multi-temporal (stars taken from different epochs). In that case, timing from config file will be overridden.")
    parser.add_argument("-f", "--firefile", type=str, required=True, help="Astrometric fireball file")

    args = parser.parse_args()

    #refcatalog = args.referencefile
    #observerCfgFile = args.cfgfile
    firefile = args.firefile
    
    firetable = Table.read(firefile, format='ascii.ecsv', guess=False, delimiter=',')
    
    if not dfn_utils.is_type_pipeline(firetable, 'astrometric'):
        print(firefile + ' does not contain astrometric data.')
        exit(1)
    
        
    if firetable.meta['astrometry_method'] != 'Global calibration':
        print('Global calibration method was not used to convert ' + firefile)
        exit(1)
        
    if 'astrometry_global_correction' in firetable.meta and firetable.meta['astrometry_global_correction'] == True:
        print(firefile + ' astrometric data has already been corrected.')
        exit(1)
        
    calib_dirname = os.path.join(os.path.dirname(firefile), 'calib')
    try:
        calibration_fits = firetable.meta['astrometry_calibration_image']
    except KeyError:
        calibration_fits = firetable.meta['astrometry_catalog']
    
    fbase = os.path.splitext(os.path.basename(calibration_fits))[0]
    
    if not os.path.isfile(calibration_fits):
        # means the catalog is probably found by looking in previous_calib.txt
        previous_calib_file = os.path.join(calib_dirname, 'previous_calib.txt')
        if os.path.isfile(previous_calib_file):
            with open(previous_calib_file) as f:
                lines = f.readlines()
            for l in lines:
                if fbase in l:
                    calib_dirname = os.path.dirname(l)
                    break
        # else let it run itself into a FileNotFoundError
    
    observerCfgFile = os.path.join(calib_dirname, fbase + '_observer.cfg')
    
    refcatalog_file = os.path.join(calib_dirname, fbase + '.xml')
    
    refcatalog = Table.read(refcatalog_file, format='votable')
    
    # Give the table its own file name as metadata
    refcatalog.meta['self_file_name'] = refcatalog_file
    
    station,_ = parse_observer_cfg_file(observerCfgFile)
    
    cut_cata = do_corrections_and_compute_errors(refcatalog, station, firetable)
    
    cut_cata_fname = refcatalog_file.split('.')[0] + '_processed.xml'
    
    cut_cata.write(cut_cata_fname, format='votable', overwrite=True)
    
    firetable.meta['astrometry_global_correction'] = True
    
    if 'astrometry_catalog' in firetable.meta and '.fits' in firetable.meta['astrometry_catalog']:
        firetable.meta['astrometry_calibration_image'] = firetable.meta['astrometry_catalog']
        
    firetable.meta['astrometry_catalog'] = os.path.basename(refcatalog_file)
    
    firetable.write(firefile, format='ascii.ecsv', delimiter=',', overwrite=True)
    
