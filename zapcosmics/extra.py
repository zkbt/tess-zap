
    def calculate(self, timeseries):
        '''Populate unbinned and binned arrays, using the filter defined for this Strategy.'''

        # setup the binned and unbinned arrays
        self.timeseries = timeseries

        # set up the unbinned dictionary
        self.unbinned = {'x':self.timeseries.x.flatten(),
                         'flux':self.timeseries.flux.flatten(),
                         'nocosmics':(self.timeseries.flux - self.timeseries.cosmics).flatten()}

        # bin, using the Strategy
        self.binned = self.combine()

        # calculate RMS of the binned light curves, with and without the cosmic ray filter
        self.unmititigated = np.std(self.binned['naive'])
        self.achieved = np.std(self.binned['flux'])

     def combine(self):
        '''Bin the unbinned timeseries, using the appropriate filter.'''

        # define x coordinate for the binned timeseries
        x = np.arange(0.5, self.timeseries.nexposures +0.5, 1)

        # naive timeseries is binned using a simple mean
        naive = self.mean()

        # noirseless binned timeseries
        noiseless = self.noiselessmean()

        # flux timeseries is binned using this Strategy's filter
        flux = self.filter()

        # nocosmics timeseries is binned using a simple mean, but assuming the cosmics were never injected
        nocosmics = np.mean(self.timeseries.flux - self.timeseries.cosmics, 1)

        # return the binned timeseries
        self.binned = {'x':x, 'flux':flux, 'nocosmics':nocosmics, 'naive':naive, 'noiseless':noiseless}
        return self.binned







,
'cosmics':self.subcadencecosmics,
'nocosmics':self.subcadenceflux - self.subcadencecosmics}
