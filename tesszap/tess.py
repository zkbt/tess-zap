from .timeseries import *
from .noise import noise

class TimeseriesTESS(Timeseries):
    def __init__(self,  tmag, noisekw={}, **kw):
        """
        A wrapper that passes everything through to Timeseries,
        after setting the subcadence and estimating cosmic ray
        parameters.

        Parameters
        ----------

        tmag : float
            TESS magnitude, for estimating the noise
        noisekw : optional, dict
                Additional keywords to pass along to the `noise()` function for
                estimating TESS noise. For example, a custom RA + Dec can be
                supplied for position-dependent sky backgrounds.
        """

        self.tmag = tmag

        # estimate the noise on the subcadence timescale
        kw['subcadence'] = 2.0
        uncertainty, photons = noise(self.tmag,
                                     exptime=kw['subcadence'],
                                     return_photons=True, **noisekw)
        kw['subcadenceuncertainty'] = np.mean(uncertainty)

        # very approximate! these should be updated during commissioning!
        tesscosmicelectrons = 1000.0
        tesscosmicprobability = 0.001

        # calculate how many sigma the cosmic will appear as in the binned cadences
        dilution = kw['cadence']/kw['subcadence']
        fractioninsubcadence = tesscosmicelectrons/photons
        fractionincadence = fractioninsubcadence/dilution
        sigmaincadence = kw['subcadenceuncertainty']/np.sqrt(dilution)
        nsigmaincadence = fractionincadence/sigmaincadence

        kw['addcosmics'] = True
        kw['cosmickw'] = dict(height=nsigmaincadence,
                              probability=tesscosmicprobability)

        # initialize the timeseries with these updated
        Timeseries.__init__(self, **kw)


class StampCadenceTimeseries(TimeseriesTESS):
    def __init__(self, tmin=-0.5, tmax=0.5, model=np.ones_like, tmag=10, noisekw={}):
        '''
        Initialize a short-cadence TESS light curve.

        Parameters
        ----------
        tmin : float
            The minimum time, in *days*.
        tmax : float
            The maximum time, in *days*.
        model : function
            A callable object, that takes time as an input and returns model
            fluxes. For TESS light curves, to get the cosmic ray injection
            approximately right, the model should return values that are
            generally close to 1.
        tmag : float
            TESS magnitude, for estimating the noise
        noisekw : optional, dict
            Additional keywords to pass along to the `noise()` function for
            estimating TESS noise. For example, a custom RA + Dec can be
            supplied for position-dependent sky backgrounds.
        '''

        cadence = 120.0
        TimeseriesTESS.__init__(self, tmin=tmin, tmax=tmax,
                                      model=model, tmag=tmag,
                                      cadence=cadence, noisekw=noisekw)

class FFICadenceTimeseries(TimeseriesTESS):
    def __init__(self, tmin=-0.5, tmax=0.5, model=np.ones_like, tmag=10, noisekw={}):
        '''
        Initialize a long-cadence TESS light curve.
        Note, this is expecting a model that is generally pretty close to 1.

        Parameters
        ----------
        tmin : float
            The minimum time, in *days*.
        tmax : float
            The maximum time, in *days*.
        model : function
            A callable object, that takes time as an input and returns model
            fluxes. For TESS light curves, to get the cosmic ray injection
            approximately right, the model should return values that are
            generally close to 1.
        tmag : float
            TESS magnitude, for estimating the noise
        noisekw : optional, dict
            Additional keywords to pass along to the `noise()` function for
            estimating TESS noise. For example, a custom RA + Dec can be
            supplied for position-dependent sky backgrounds.
        '''

        cadence = 1800.0
        TimeseriesTESS.__init__(self, tmin=tmin, tmax=tmax,
                                      model=model, tmag=tmag,
                                      cadence=cadence, noisekw=noisekw)
