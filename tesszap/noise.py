"""
Calculate the expected photometric precision for TESS.

(translated from Josh Winn's IDL TESS signal-to-noise calculator
 on the TESS wiki and updated to include calculations
 published with Peter Sullivan's simulation paper).

Be careful:
    -this calculator uses outdated PSF models
    -it pulls some ingredients calculated for Ic magnitudes,
        and treats them immediately as TESS magnitudes,
        thus neglecting a ~0.1 magnitude color term
 """

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
import scipy.interpolate

# create an interpolator to estimate the best number of pixels in a photometric aperture

sullivan_tmag = [0, 5.905, 6.281, 6.903, 7.09, 7.137, 7.459, 7.853, 8.683, 9.764, 10.395, 11.305, 12.616, 13.86, 30]
sullivan_npix = [35, 34.482, 32.564, 29.034, 27.833, 26, 23.965, 22.009, 18.772, 14.383, 11.321, 7.727, 4.723, 3.059, 3]
assert(len(sullivan_tmag) == len(sullivan_npix))
optimal_npix = scipy.interpolate.interp1d(sullivan_tmag, sullivan_npix,
                                            kind='linear', bounds_error=True)

sullivan_enclosednpix = [1.003, 1.285, 2.023, 2.897, 4.745, 7.159, 10.456, 16.958, 24.59, 33.466, 51.732]
sullivan_enclosedfraction = [0.325, 0.396, 0.526, 0.638, 0.79, 0.873, 0.925, 0.982, 1, 1, 1]
enclosed_fraction = scipy.interpolate.interp1d(sullivan_enclosednpix, sullivan_enclosedfraction)

def noise(tmag=10.0, exptime=1800.0, teff=5000.0,
          elon=0.0, elat=30.0, glon=None, glat=None, ra=None, dec=None,
          subexptime=2.0, npix_aper=4, frac_aper=0.76, e_pix_ro=10.0,
          effective_area=73.0, pix_scale=21.1, sys_limit=60.0,
          verbose=False, return_photons=False):
    """Calculate noise, given input TESS magnitude, returing the fractional rms (= 1/snr)

        Mandatory inputs

           tmag,                            apparent TESS magnitude

        Optional inputs

           exptime                           total exposure time in seconds
           teff                              effective temperature in Kelvins
           elon, elat                        ecliptic coordinates in degrees
           subexptime                        subexposure time (n_exp = exptime/subexptime)
           npix_aper                         number of pixels in photometric aperture
           frac_aper                         fraction of flux enclosed in photometric aperture
           e_pix_ro                          rms in no. photons/pixel from readout noise
           effective_area                    geometric collecting area
           pix_scale                         arcsec per pixel
           sys_limit                         minimum uncertainty in 1 hr of data, in ppm
           verbose                           request verbose output
    """


    # pick the optimal number pixels in the photometric aperture
    npix_aper = optimal_npix(tmag)

    # determine the fraction of the stellar flux included
    frac_aper = enclosed_fraction(npix_aper)

    # solid area of a pixel
    omega_pix = pix_scale ** 2.

    # how many subexposures composed this one exposure?
    n_exposures = exptime / subexptime

    # the TESS zeropint
    tmag0 = 1.514e6

    # photoelectrons from the star
    e_star = 10.0 ** (-0.4 * tmag) * tmag0 * effective_area * exptime * frac_aper

    if ra is not None and dec is not None:
        c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
        ecl = c.barycentrictrueecliptic
        elon, elat = ecl.lon.deg, ecl.lat.deg

    if verbose:
        print('tmag = {}'.format(tmag))
        print('tmag = {}'.format(tmag))
        print('tmag0 = {}'.format(tmag0))
        print('exptime = {}'.format(exptime))
        print('teff = {}'.format(teff))
        print('elon = {}'.format(elon))
        print('elat = {}'.format(elat))
        print('npix_aper = {}'.format(npix_aper))
        print('frac_aper = {}'.format(frac_aper))
        print('subexptime = {}'.format(subexptime))
        print('n_exposures = {}'.format(n_exposures))
        print('e_pix_ro = {}'.format(e_pix_ro))
        print('effective_area = {}'.format(effective_area))
        print('pix_scale = {}'.format(pix_scale))
        print('omega_pix = {}'.format(omega_pix))
        print('sys_limit = {}'.format(sys_limit))
        print('e_star = {}'.format(e_star))

    # photoelectrons/pixel from zodiacal light
    dlat = (np.abs(elat) - 90.) / 90.
    vmag_zodi = 23.345 - 1.148 * dlat ** 2.
    e_pix_zodi = 10.0 ** (-0.4 * (vmag_zodi - 22.8)) * 2.39e-3 * effective_area * omega_pix * exptime


    # photoelectrons/pixel from unresolved background stars
    ecl = SkyCoord(lon=elon*u.degree, lat=elat*u.degree, frame='barycentrictrueecliptic')
    glon, glat = ecl.galactic.l.deg, ecl.galactic.b.deg

    glon = np.array([glon])
    glat = np.array([glat])

    if verbose:
        print('vmag_zodi = {}'.format(vmag_zodi))
        print('e_pix_zodi = {}'.format(e_pix_zodi))
        print('glon = {GLON}, glat = {GLAT}'.format(GLON=glon, GLAT=glat))

    dlat = np.abs(glat) / 40.0
    dlon = glon
    q = (dlon > 180.)
    dlon[q] = 360. - dlon[q]
    dlon = np.abs(dlon) / 180.0
    p = [18.9733, 8.833, 4.007, 0.805]
    tmag_bgstars = p[0] + p[1] * dlat + p[2] * dlon ** (p[3])
    e_pix_bgstars = 10.0 ** (-0.4 * tmag_bgstars) * 1.7e6 * effective_area * omega_pix * exptime

    if verbose:
        print('tmag_bgstars = {}'.format(tmag_bgstars))
        print('e_pix_bgstars = {}'.format(e_pix_bgstars))

    noise_star = np.sqrt(e_star) / e_star
    noise_sky = np.sqrt(npix_aper * (e_pix_zodi + e_pix_bgstars)) / e_star
    noise_ro = np.sqrt(npix_aper * n_exposures) * e_pix_ro / e_star
    noise_sys = 0.0 * noise_star + sys_limit / 1e6 / np.sqrt(exptime / 3600.)
    noise = np.sqrt(noise_star ** 2. + noise_sky ** 2. + noise_ro ** 2. + noise_sys ** 2.)

    if verbose:
        print('noise_star [ppm] = ', noise_star*1e6)
        print('noise_sky  [ppm] = ', noise_sky*1e6)
        print('noise_ro   [ppm] = ', noise_ro*1e6)
        print('noise_sys  [ppm] = ', noise_sys*1e6)
        print('noise      [ppm] = ', noise*1e6)

    if return_photons:
        return noise, e_star
    else:
        return noise


def demo(span=27.4, period=12.345678, mean=17, amplitude=1.0):
    """
    Demonstration of the TESS noise calculator,
    by default on a faint and highly-variable star.
    """

    # create times at a half hour spacing
    t = np.arange(0, span, 0.5 / 24.0)
    n = len(t)

    # create a (perfectly smooth) noiseless model
    noiselessmodel = mean + amplitude * np.sin(2 * np.pi * t / period)

    # calculate the per-point photometric uncertainty for each point of that model
    # noinspection PyTypeChecker
    perpointuncertainty = noise(tmag=noiselessmodel)

    # create one random realization of the noise
    noiserealization = np.random.normal(0, 1, n) * perpointuncertainty

    # simulate measurements
    simulated = noiselessmodel + noiserealization

    # create plot showing the demonstration
    plt.ion()
    plt.figure('demonstration', figsize=(10, 3), dpi=200)
    plt.cla()
    # plot the simulated measurements
    plt.errorbar(t, simulated, perpointuncertainty, marker='o', elinewidth=2, linewidth=0, color='black', alpha=0.5)
    # plot the noiseless model
    plt.plot(t, noiselessmodel, color='green', linewidth=2, alpha=0.5)
    # clean up the look of the plot
    # noinspection PyTypeChecker
    plt.ylim(mean + amplitude + np.max(perpointuncertainty) * 5, mean - amplitude - np.min(perpointuncertainty) * 5)
    plt.xlim(np.min(t), np.max(t))
    plt.xlabel('Time (in days)')
    plt.ylabel('Flux (magnitudes)')
    plt.tight_layout()

    return t, simulated, perpointuncertainty
