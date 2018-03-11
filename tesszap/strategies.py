'''
Create some cosmic-zapping strategies, using different techniques for
stacking together individual subexposures
'''

from .imports import *

class Strategy():
    '''
    Define structures and methods needed by all cosmic ray rejection strategies.
    '''
    def __init__(self, n=10):
        '''Initialize Strategy object.'''

        # setup basic filtering parameters:
        self.name = 'Nothing'

        # most filters have a single scale (some have more)
        self.n = n

        # define a dictionary of plotting parameters
        self.plotting = {   'nsigma':10,
                            'toplot':96,
                            'flux':'mediumvioletred',
                            'nocosmics':'orange',
                            'naive':'blue'}


    def __call__(self, timeseries):
        raise Exception("Uh-oh! It seems no filtering function was defined!")


class Central(Strategy):
    '''
    Binning Strategy =
            break into subsets, reject the highest and lowest points from each
            and take the mean of the rest, sum these truncated means.
    '''
    def __init__(self, n=None):
        Strategy.__init__(self, n)
        m = None
        if m is None:
          m = self.n - 2
        self.m = m
        self.name = "Central {m} out of {n}".format(m = self.m, n=self.n)

    def __call__(self, timeseries):
        assert(self.m == self.n-2)
        shape = timeseries.shape
        reshapen = timeseries.unbinned['flux'].reshape(shape[0], shape[1]//self.n, self.n)
        sum = np.sum(reshapen, 2)
        max = np.max(reshapen, 2)
        min = np.min(reshapen, 2)
        corrected = np.mean((sum - max - min)/self.m,1)
        return corrected #- np.mean(corrected)


class Mean(Strategy):
    '''
    Binning Strategy =
        the simplest, taking the straight sum of all the subexposures.'''
    def __init__(self):
        Strategy.__init__(self)
        self.name = "Mean"

    def __call__(self, timeseries):
        return np.mean(timeseries.unbinned['flux'], 1)


class Median(Strategy):
    '''
    Binning Strategy =
        break into non-overlapping subsets, take median of each, sum all of
        those medians together.
    '''
    def __init__(self, n=None):
        Strategy.__init__(self, n)
        self.name = "Median of {0}".format(self.n)

    def __call__(self, timeseries):
        shape = timeseries.shape
        reshapen = timeseries.unbinned['flux'].reshape(shape[0], shape[1]//self.n, self.n)
        partial = np.median(reshapen, 2)
        return np.mean(partial, 1)

class ShiftedMedian(Strategy):
    '''
    Binning Strategy =
            break into overlapping subsets, take median of each, sum all of
            those medians together.
    '''
    def __init__(self, n=None):
        Strategy.__init__(self, n)
        self.name = "Shifted Median of {0}".format(self.n)

    def __call__(self, timeseries):
        shape = timeseries.shape
        sumofmedians = np.zeros((shape[0], shape[1]//self.n))
        for i in range(self.n):
          # note, this wraps around to the start at the end of the array; but should work for testing
          reshapen = np.roll(timeseries.unbinned['flux'], i, 1).reshape(shape[0], shape[1]//self.n, self.n)
          sumofmedians += np.median(reshapen, 2)

        return np.mean(sumofmedians/self.n, 1)


class Sullivan(Strategy):
    '''
    Binning Strategy =
        idea Peter and I played around with; didn't work particularly well.
    '''

    def __init__(self, n=None):

        n = 3
        Strategy.__init__(self, n)
        self.name = "Sullivan of {0}".format(self.n)

    def __call__(self, timeseries):
        shape = timeseries.shape
        reshapen = timeseries.unbinned['flux'].reshape(shape[0], shape[1]//self.n, self.n)
        median = np.median(reshapen, 2)
        d = np.sum(reshapen, 2)/3.
        a = (3*d - reshapen[:,:,0])/2.
        b = (3*d - reshapen[:,:,1])/2.
        c = (3*d - reshapen[:,:,2])/2.

        diffs = np.array([np.abs(a-median), np.abs(b - median), np.abs(c -median), np.abs(d - median)])
        mask = diffs == np.min(diffs, 0).reshape(1,diffs.shape[1], diffs.shape[2])*np.ones_like(diffs)

        values = np.array([a,b,c,d])
        partial = np.sum(mask*values,0)/np.sum(mask, 0)
        return np.mean(partial, 1)

class Lowest(Strategy):
    '''
    Binning Strategy =
            break into subsets, reject the highest point from each and take the
            mean of the rest, sum these truncated means.
    '''
    def __init__(self, n=None):
        Strategy.__init__(self, n)
        m = None
        if m is None:
          m = self.n - 1
        self.m = m
        self.name = "Lowest {m} out of {n}".format(m = self.m, n=self.n)

    def __call__(self, timeseries):
        assert(self.m == self.n-1)
        shape = timeseries.shape
        reshapen = timeseries.unbinned['flux'].reshape(shape[0], shape[1]//self.n, self.n)
        sum = np.sum(reshapen, 2)
        max = np.max(reshapen, 2)
        corrected = np.mean((sum - max)/self.m,1)
        return corrected#/ np.mean(corrected)    # this isn't correct -- it should be an empirical correction factor!

class OutlierWithDecay(Strategy):
    '''
    Binning Strategy =
            break into subsets, estimate standard deviation with weighted mean
            of current subset and previous estimate, reject from each and take
            the mean of the rest, sum these truncated means.
    '''
    def __init__(self, n=10, threshold=10.0, memory=0.90, safetybuffer=2.0, diagnostics=False):
        '''Initialize an outlierwith decay Strategy.

        n = the number of subexposures in each "chunk"
        threshold = how many sigma about the noise are required for a point to be an outlier?
        memory = the weight given to previous estimates of standard deviations (best estimate = memory[previous best estimate] + (1 - memory)[most recent chunk])
        safetybuffer = by what factor should we inflate the initial standard deviation to prevent overfitting?'''

        # initialize the basic Strategy class
        Strategy.__init__(self, n)

        # store the parameters of the filter
        self.threshold = threshold
        self.memory = memory
        self.safetybuffer = 2.0

        # define a name for this filter
        self.name = r'Rejecting {threshold}$\sigma$ Outliers; $\sigma$ from chunks of {n};  memory of {memory}%'.format(threshold=self.threshold, memory=self.memory*100, n=self.n)

        # for testing, keep a diagnostics flag to say whether to display the mean + std. estimates
        self.diagnostics = diagnostics

    @property
    def label(self):
        return 'RejectionWith{0:.0f}percentMemory'.format(self.memory*100)


    def strategyprefix(self):
        '''Custom Strategy prefix for this long named Strategy.'''
        name = 'RejectionWith{0:.0f}percentMemory'.format(self.memory*100)
        f = os.path.join(self.directory(), 'crdemo_filter{0}_n{1}'.format(name.replace(' ',''), self.n))
        return f

    def __call__(self, timeseries):
        '''Loop through the unbinned light curve, applying the realtime outlier-rejection filter.'''

        # make sure that the chunk size divides evenly into an exposure
        assert(timeseries.nsubexposures % self.n == 0)
        assert(self.n > 3)

        # define the number of chunks (and other convenience constants)
        nchunks = timeseries.nsubexposures//self.n
        nsubexposures = timeseries.nsubexposures
        nexposures = timeseries.nexposures
        n = self.n

        # create an array to store the final binned timeseries
        finaltimeseries = np.zeros(nexposures)

        # create arrays to store the per-chunk estimates of the mean and the standard deviation (these are arrays simply for plotting purposes)
        running_mean, running_std = np.zeros((nexposures, nchunks)), np.zeros((nexposures, nchunks))

        # initialize these estimates with the first chunk (rejecting no outliers)
        running_mean[0,0] = np.mean(timeseries.unbinned['flux'][0,0:n])
        running_std[0,0] = np.sqrt(np.sum((timeseries.unbinned['flux'][0,0:n] - running_mean[0,0])**2)/(n-1.0))

        # inflate the initial standard deviation measurement to prevent overfitting
        running_std[0,0] *= self.safetybuffer

        # set the first binned point to this mean estimate (there's no other choice)
        finaltimeseries[0] = running_mean[0,0]

        # loop over the binned exposures, and chunks within exposures
        count = 1
        for iexposure in np.arange(nexposures):
          for ichunk in np.arange(nchunks):

            # skip the very first point, because it's already been defined
            if (ichunk == 0)&(iexposure == 0):
              continue

            # pull out the light curve for this little chunk
            flux = timeseries.unbinned['flux'][iexposure, n*ichunk:n*(ichunk+1)]
            # pull out the last mean and standard deviation estimates
            best_mean = running_mean.flatten()[count-1]
            best_std = running_std.flatten()[count-1]

            # determine which points are not outliers, by being less than [threshold] sigma over the mean
            notoutlier = flux < (best_mean + self.threshold*best_std)
            assert(notoutlier.any())
            if notoutlier.any() == False:
              notoulier = np.ones_like(flux)
            # determine the mean and standard deviation of the good points in this chunk
            this_mean = np.mean(flux[notoutlier])
            this_std = np.sqrt(np.sum((flux[notoutlier] - this_mean)**2)/(n - 1.0))

            # store this binned exposure in the final array
            finaltimeseries[iexposure] += this_mean

            # mix this chunk into the running estimates, for the next chunk to use
            running_mean[iexposure, ichunk] = self.memory*best_mean + (1.0 - self.memory)*this_mean
            running_std[iexposure, ichunk] = self.memory*best_std + (1.0 - self.memory)*this_std
            # or would it be better to be doing this in variance space?, such as...
            # running_std[iexposure, ichunk] = np.sqrt(self.memory*best_std**2 + (1.0 - self.memory)*this_std**2)

            # advance the linear counter
            count += 1

        if self.diagnostics:

          # create a plot for showing the mean and standard deviation estimates
          plt.figure('diagnostics for {0}'.format(self.__class__), figsize=(7,4), dpi=150)
          plt.cla()

          # plot the two timeseries
          pargs = dict(linewidth=3, alpha=0.5)
          plt.plot(running_mean.flatten(), linestyle='--', label='mean', **pargs)
          plt.plot(running_std.flatten(), label='standard deviation', **pargs)
          plt.legend()

        # return the binned timeseries
        return finaltimeseries/nchunks
