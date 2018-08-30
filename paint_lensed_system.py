# -*- coding: utf-8 -*-
"""Demo on using LensPop to "paint" lensed quasar systems

This module demonstrates how to use LensPop to "paint" a lensed quasar system,
i.e. get the magnitudes of the lens and the quasar in multiple filters
(the SDSS filters UGRIZ in this example)

Example
-------
This demo can be run from the root of the LensPop repository by simply typing::

    $ python paint_lensed_system.py

References
----------
    This code was refactored from the db.py script of the 
    OM10 package repository, at https://github.com/drphilmarshall/om10.
    The relevant method in the script is `paint`.

    Use the following citation for OM10:
    Oguri, Masamune, and Philip J. Marshall. "Gravitationally lensed quasars and 
    supernovae in future wide-field optical imaging surveys." 
    Monthly Notices of the Royal Astronomical Society 405.4 (2010): 2579-2593.

"""

from __future__ import division, absolute_import, print_function

def get_sdss_filters():
    """Retrieves the SDSS filter throughputs from files in stellarpop/filters

    Returns
    -------
    dict
        Each key is one of string characters 'u', 'g', 'r', 'i', 'z'
        Each value is a tuple of (Numpy ndarray of wavelengths in Angstroms,
        Numpy ndarray of throughputs) 
    """
    from stellarpop import tools

    filters_dict = {}
    for band in 'ugriz':
        filters_dict[band] = tools.filterfromfile('%s_SDSS' %band)
    return filters_dict

def calculate_rf_quasar_magnitudes(quasar_redshift, apparent_magnitude_i, filters_dict):
    """Calculates the reference-frame quasar magnitudes in multiple filters

    Parameters
    ----------
    quasar_redshift : float
        Redshift of the source quasar
    apparent_magnitude_i : float
        Apparent magnitude of the source AGN in the I-band
    filters_dict : dict
        (See output of `get_sdss_filters` for details)
        Throughputs of various filters

    Returns
    -------
    dict
        Each key is one of string characters 'u', 'g', 'r', 'i', 'z'
        representing the filter
        Each value is the reference-frame apparent magnitude of the quasar
        in the 'key' filter, of type float
    """
    from stellarpop import tools

    quasar_sed = tools.getSED('agn')
    q_offset = apparent_magnitude_i - tools.ABFilterMagnitude(filters_dict['i'], quasar_sed, quasar_redshift)
    rf_quasar_appmag = {}
    if quasar_redshift < 3.9:
        rf_quasar_appmag['u'] = tools.ABFilterMagnitude(filters_dict['u'], quasar_sed, quasar_redshift) + q_offset
    else:
        rf_quasar_appmag['u'] = 99.0
    for band in 'griz':
        rf_quasar_appmag[band] = tools.ABFilterMagnitude(filters_dict[band], quasar_sed, quasar_redshift) + q_offset
    return rf_quasar_appmag

def calculate_rf_lens_magnitudes(lens_redshift, velocity_dispersion, filters_dict):
    """Calculates the reference-frame lens magnitudes in multiple filters

    Parameters
    ----------
    lens_redshift : float
        Redshift of the lens
    velocity_dispersion : float
        Velocity dispersion of the lens
    filters_dict : dict
        (See output of `get_sdss_filters` for details)
        Throughputs of various filters

    Returns
    -------
    dict
        Each key is one of string characters 'u', 'g', 'r', 'i', 'z'
        representing the filter
        Each value is the reference-frame apparent magnitude of the quasar
        in the 'key' filter, of type float
    """
    from stellarpop import tools
    from lenspop import population_functions, distances
    from astropy.cosmology import FlatLambdaCDM

    # Instantiate Distance
    distance = distances.Distance() #TODO: necessary?
    # Instantiate LensPopulation
    lenspop = population_functions.LensPopulation_()
    # Instantiate FlatLambdaCDM cosmology with reasonable parameters
    cosmology = FlatLambdaCDM(H0=70.0, Om0=0.3)

    lens_sed = tools.getSED('BC_Z=1.0_age=9.000gyr')
    velocity_dispersion = velocity_dispersion

    # Absolute --> apparent magnitude conversion in the R-band
    lens_abmag_r = tools.ABFilterMagnitude(filters_dict['r'], lens_sed, lens_redshift)
    distance_modulus = cosmology.distmod(lens_redshift).value
    lens_appmag_r = lens_abmag_r + distance_modulus

    # [Reference frame] Absolute --> apparent magnitude conversion in the R-band
    rf_lens_abmag_r, _ = lenspop.EarlyTypeRelations(velocity_dispersion)
    rf_lens_appmag = {}
    rf_lens_appmag['r'] = rf_lens_abmag_r + distance_modulus

    # Quantity which is added to ~magnitude to convert it into reference-frame ~magnitude
    offset_rf = rf_lens_abmag_r - lens_abmag_r

    # Converting absolute magnitude to reference-frame apparent magnitude
    for band in 'ugiz':
        rf_lens_appmag[band] = tools.ABFilterMagnitude(filters_dict[band], lens_sed, lens_redshift) + offset_rf + distance_modulus
    return rf_lens_appmag

if __name__=="__main__":
    """Grab a random lensed quasar system from the OM10 catalog
    and give the lens and the quasar multi-filter magnitudes

    Notes
    -----
    The following code snippet will grab information about a lensed quasar system
    in dictionary form if you have OM10 installed::
        from om10 import DB 
        lens_catalog = om10.DB()
        example_system = lens_catalog.sample['LENSID'==1077370]

    But since this is just an example, we will simply pre-fetch the information
    of the system relevant to running this demo script.

    Attributes
    ----------
    lens_redshift : float
        Redshift of the lens
    velocity_dispersion : float
        Velocity dispersion of the lens
    quasar_redshift : float
        Redshift of the quasar
    apparent_magnitude_i : float
        Apparent magnitude of the quasar in the I-band

    """
    lens_redshift = 0.276
    velocity_dispersion = 196.7886
    quasar_redshift = 2.33
    apparent_magnitude_i = 23.59

    filters_dict = get_sdss_filters()
    rf_quasar_appmag = calculate_rf_quasar_magnitudes(quasar_redshift, apparent_magnitude_i, filters_dict)
    rf_lens_appmag = calculate_rf_lens_magnitudes(lens_redshift, velocity_dispersion, filters_dict)
    print("Multi-filter quasar magnitudes: ", rf_quasar_appmag)
    print("Multi-filter lens magnitudes: ", rf_lens_appmag)