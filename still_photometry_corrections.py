#!/usr/bin/env python

# general imports
import os
import datetime
import logging
import itertools

# sci packages imports
import numpy as np
from astropy.table import Table
from astropy.time import Time
import astropy.units as u

from scipy.interpolate import interp1d

# plotting imports
import matplotlib.pyplot as plt
import seaborn as sns

# local imports
import dfn_utils


__author__ = "Hadrien Devillepoix"
__copyright__ = "Copyright 2016, Desert Fireball Network"
__license__ = "MIT"
__version__ = "1.0"
__scriptName__ = "still_photometry_corrections.py"
    
    
def correct_photometric_measurements(table):
    '''
    Do corrections on still photometric measures
    
    airmass
    range
    zero-point
    '''
    
    # check raw data is available
    if not dfn_utils.is_type_pipeline(table, 'raw_photometric'):
        raise dfn_utils.WrongTableTypeException('Need raw measurements to be able to compute corrected light curve!')
    
    # check correction data is available
    if not all(key in table.meta for key in ['photometric_zero_point',
                                             'photometric_zero_point_uncertainty',
                                             'photometric_zero_point_datetime']):
        raise dfn_utils.WrongTableTypeException('Need zero point measurements to be able to compute corrected light curve!')
    
    # check range (distance) data is available
    if not dfn_utils.is_type_pipeline(table, 'triangulated'):
        raise dfn_utils.WrongTableTypeException('Need observations ranges to be able to compute corrected light curve!')
    
    
    # do the corrections
    # TODO properly handle spectral band parameters
    #spectral_band = dfn_utils.PROCESSING_FILTER_BAND[processing_filter]
    #flux_col_ = 'flux_dash_' + spectral_band
    flux_col_ = 'flux_dash_V'
    
    photometric_zero_point = table.meta['photometric_zero_point']
    
    table['M_V'] = (photometric_zero_point -
                    2.5 * np.log10(table[flux_col_] / (1. * u.adu/u.second)) -
                    5 * np.log10(table['range']/(100 *u.km)))
    
    # update metadata on what's been done
    now_timestamp = datetime.datetime.now()
    writetime = now_timestamp.strftime('%Y-%m-%dT%H:%M:%S')
    table.meta['photometry_corrections_software'] = __scriptName__ + ' ' + __version__
    table.meta['photometry_corrections_write_time'] = writetime
    
    
    return table


def integrate_light_curve(table, tau=0.005):
    '''
    integrate the absolute light curve
    '''
    
    
    
    # interpolate missing magnitudes
    t = Time(table['datetime'])
    rel_t = (t - t[0]).sec
    # create missing values mask
    missing = np.isnan(table['M_V'])

    # create interpolator using only valid input
    interp_magnitudes = interp1d(rel_t[~missing], table['M_V'][~missing], fill_value='extrapolate')

    filled_table = Table()
    filled_table['relative_time'] = rel_t
    filled_table['M_V'] = interp_magnitudes(rel_t)
    filled_table['velocity'] = table['D_DT_bal_fitted']
    
    # Calculate calibrated intensity (magnitude system back to energetic)
    filled_table['I_V'] = 10**(-0.4 * filled_table['M_V'])
    ceplecha_magic_coefficient = 1.95e10
    filled_table['I'] = filled_table['I_V'] * ceplecha_magic_coefficient
    
    integrated_lightcurve = 0.
    for i in range(len(filled_table)-1):
        integrated_lightcurve += (filled_table['I'][i] /
                                  ((filled_table['velocity'][i]/1000) **2) *
                                  (filled_table['relative_time'][i+1] - filled_table['relative_time'][i]))
    
    integrated_magnitude = -2.5 * np.log10(integrated_lightcurve/ceplecha_magic_coefficient)
    
    photometric_mass = integrated_lightcurve / tau
    
    table.meta['photometry_integrated_M_V'] = float(integrated_magnitude)
    table.meta['photometry_photometric_mass'] = float(photometric_mass)
    table.meta['photometry_luminous_efficiency_used'] = float(tau)
    
    
    
    ## get rid of missing values for the sake of integration
    #valid_photom_measurements = table[~np.isnan(table['I'])]
    ## define a relative time column to feed integrator
    #t = Time(valid_photom_measurements['datetime'].data.tolist())
    #valid_photom_measurements['rel_time'] = Time(t - t[0], format='sec').sec
    
    ## integrate over time
    #integrated_flux = np.trapz(valid_photom_measurements['I'], x=valid_photom_measurements['rel_time'])
    
    #table.meta['photometry_integrated_luminosity'] = integrated_flux.value
    #print('integrated luminosity calculated: {0:.2e}'.format(integrated_flux))
    
    #integrated_magnitude = -2.5 * np.log10(integrated_flux/ceplecha_magic_coefficient)
    #table.meta['photometry_integrated_M_V'] = integrated_magnitude.value
    #print('integrated magnitude calculated: {0:.2f} mag'.format(integrated_magnitude))
    
    

def light_curve_plot(tables, kwargs=None):
    '''
    Plots light curves, as measured from different cameras
    '''
    
    logger = logging.getLogger('trajectory')
    
    phot_tables = [tab for tab in tables if dfn_utils.is_type_pipeline(tab, 'absolute_photometric')]
    
    if len(phot_tables) < 1:
        logger.error('No absolutely calibrated photometry to plot')
        return
    
    logger.debug('Plotting lightcurve using: ' + str([tab.meta['telescope'] for tab in phot_tables]))
    
    plt.clf()
    palette = itertools.cycle(sns.color_palette())
    # Handle datetime axis
    plt.gcf().autofmt_xdate()

    for table in phot_tables:
        non_nan = ~np.isnan(table['M_V'])
        t = Time(table['datetime'].data.tolist())[non_nan]

        plt.plot_date(t.plot_date, table['M_V'][non_nan],
                      xdate=True, ydate=False, color=next(palette),
                      ls='solid', marker='.',
                      label=table.meta['telescope'] + " " + table.meta['location'])

    try:
        event_codename = table.meta['event_codename']
    except:
        event_codename = "fireball"

    plt.gcf().autofmt_xdate()  # orient date labels at a slant
    # inverse y axis (mag plot)
    plt.gca().invert_yaxis()
    plt.draw()
    plt.legend(frameon=True, loc='best', fancybox=True, framealpha=0.5)
    plt.xlabel("Time (t0 + ). t0= " + str(Time(t[0], format='isot', out_subfmt='date')) + " UTC")
    plt.ylabel("M_V")
    #plt.title(event_codename.decode(encoding='UTF-8') + " - Timing consistency check between cameras")
    plt.title(event_codename + " - Light curve")

    dirname = os.path.dirname(phot_tables[0].meta['self_disk_path'])
    #fname = os.path.join(dirname, event_codename.decode(encoding='UTF-8') + "_timing_consistency_check.png")
    fname = os.path.join(dirname, event_codename + "_light_curve_" + kwargs['trajectory_segment'] + ".png")
    plt.savefig(fname)

    
def main_io(inputtable):
    '''
    High order main that deals with disk I/O
    '''
    
    # read
    table = Table.read(inputtable, format='ascii.ecsv', guess=False, delimiter=',')
    
    # do
    correct_photometric_measurements(table)
    
    # write
    table.write(inputtable, format='ascii.ecsv', delimiter=',')
    

def main():
    """
    main
    """
    # parse arguments
    import argparse
    
    parser = argparse.ArgumentParser(description='Does photometry extraction on star catalog')
    parser.add_argument("-i", "--inputtable", type=str, required=True, help="Input observation table")

    args = parser.parse_args()
    
    main_io(args.inputtable)

    
if __name__ == "__main__":
    main()
