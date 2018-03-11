
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


















    def plotTimeseries(self, x,y, **kwargs):
        '''Plot a timeseries.'''

        # plot only a subset (for quickness and clarity's sake)
        ok = x <= self.plotting['toplot']

        # plot in the timeseries panel
        ax = self.ax['timeseries']
        ax.plot(x[ok], y[ok], **kwargs)
        ax.set_xlim(0, self.plotting['toplot'])


    def test(self, t, remake=False, niterations=20):
        '''Test Strategy over a range of input cosmic ray amplitudes (using toy model light curves).'''
        self.timeseries = t
        print("testing {0}, with filename of {1}".format(self.name, self.strategyprefix()))
        filename = self.strategyprefix() + '_{0}iterations_{1:03.0f}subexposures.npy'.format(niterations, self.timeseries.nsubexposures)
        try:
          # can we load a preprocessed file?
          assert(remake == False)
          print('trying to load file from ', filename)
          self.amplitudes, self.noise = np.load(filename)
          print(self.strategyprefix() + " loaded from file")
        except:

          # create a grid of relevant cosmic ray amplitudes
          self.amplitudes = np.arange(0, 5, 0.5)
          self.noise = np.zeros((len(self.amplitudes), niterations))

          # loop over iterations, to improve precision of the noise estimates
          for iteration in np.arange(niterations):
            print("  iteration #{0}/{1}".format(iteration, niterations))
            # loop over amplitudes
            for i in range(len(self.amplitudes)):

              # create a new timeseries
              self.timeseries = Timeseries1D(nexposures=t.nexposures, nsubexposures=t.nsubexposures, amplitude=self.amplitudes[i])
              # bin it, using the Strategy
              self.calculate(self.timeseries)

              # create a demo plot of this light curve
              if iteration == 0:
                self.plot()

              # print this timeseries's summary
              print(self.prefix + "{0}".format(self.timeseries))

              # store the achieve noise
              self.noise[i, iteration] = self.achieved/self.timeseries.exposurenoise

          self.noise = np.mean(self.noise, 1)

          # save this calculation, so we can use again if need be
          np.save(filename, (self.amplitudes, self.noise))
          print('   saved results to ', filename)

        # return (x, y) for plotting
        return self.amplitudes, self.noise

    def plot(self,unbinned=False,**kwargs):
        '''Create a pretty plot comparing the binned timeseries and their histograms.'''


        # setup the plots
        includeimage = self.timeseries.toy == False
        self.figure = plt.figure(0, figsize=(9 + 3*includeimage,2), dpi=150);
        self.gs = plt.matplotlib.gridspec.GridSpec(1,6+includeimage,wspace=0,hspace=0,left=0.08, right=0.98, top=0.85, bottom=0.2)
        self.figure = plt.gcf()
        self.figure.clf()
        self.figure.suptitle("Cosmic Ray Rejection with [{0}]".format(self.name))
        self.ax = {}
        self.ax['timeseries'] = plt.subplot(self.gs[0,:-(1 + includeimage)])
        self.ax['timeseries'].set_autoscale_on(False)
        self.ax['timeseries'].set_xlabel(r'Exposure #')
        self.ax['timeseries'].set_ylabel(r'Flux')
        self.ax['histbinned'] = plt.subplot(self.gs[0,-(1 + includeimage)], sharey=self.ax['timeseries'])
        plt.setp(self.ax['histbinned'].get_yticklabels(), visible=False)
        if includeimage:
          self.ax['image'] = plt.subplot(self.gs[0,-1])
          plt.setp(self.ax['image'].get_xticklabels(), visible=False)
          plt.setp(self.ax['image'].get_yticklabels(), visible=False)




        # set up the zero-level for the plot, the width of the plotting windows, and the position of the text labels
        ylevel = np.median(self.timeseries.scale*self.binned['nocosmics'])
        ywidth = self.timeseries.scale*self.plotting['nsigma']*0.8*self.timeseries.exposurenoise
        y =  ylevel - ywidth


        # fiddle with the ylimits of the plots, depending on whether looking at binned or unbinned timeseries
        if unbinned:
          self.ax['histbinned'].set_ylim(np.min(self.timeseries.scale*self.unbinned['flux']), np.max(self.timeseries.scale*self.unbinned['flux']))
          #self.ax['histbinned'].set_ylim(-self.plotting['nsigma']*self.timeseries.subexposurenoise, self.plotting['nsigma']*self.timeseries.subexposurenoise)
        else:
          self.ax['histbinned'].set_ylim(-self.timeseries.scale*self.plotting['nsigma']*self.timeseries.exposurenoise + ylevel, self.timeseries.scale*self.plotting['nsigma']*self.timeseries.exposurenoise + ylevel)

        # plot the unbinned timeseries, if desired
        if unbinned:
          self.plotTimeseries(self.unbinned['x'], self.timeseries.scale*self.unbinned['flux'], linewidth=1, alpha=0.1, color='black')# markersize=2, markeredgewidth=0, markerfacecolor='black', marker='o', linewidth=0, alpha=0.1)
          #self.plotTimeseries(self.unbinned['x'], self.unbinned['nocosmics'], linewidth=1, alpha=0.1, color='orange')# markersize=2, markeredgewidth=0, markerfacecolor='black', marker='o', linewidth=0, alpha=0.1)

        # plot the cosmic ray hits as vertical lines
        cosmichit = ((self.unbinned['flux'] -  self.unbinned['nocosmics']) > 0)*(self.unbinned['x']<=self.plotting['toplot'])
        for x in self.unbinned['x'][cosmichit]:
          self.ax['timeseries'].axvline(x, color='black', alpha=0.1, linewidth=0.5)


        # plot the binned timeseries
        self.plotTimeseries(self.binned['x'], self.timeseries.scale*self.binned['nocosmics'], markersize=4, marker='o', alpha=0.5, markerfacecolor='orange', color=self.plotting['nocosmics'], markeredgewidth=0, linewidth=1)
        self.plotTimeseries(self.binned['x'], self.timeseries.scale*self.binned['naive'], markersize=4, marker='o', alpha=0.5, markerfacecolor=self.plotting['naive'], color=self.plotting['naive'], markeredgewidth=0, linewidth=1)
        self.plotTimeseries(self.binned['x'], self.timeseries.scale*self.binned['flux'], markersize=4, marker='o', alpha=0.5, markerfacecolor=self.plotting['flux'], color=self.plotting['flux'], markeredgewidth=0, linewidth=1)

        # plot the binned histogram
        self.plothistogram(self.timeseries.scale*self.binned['nocosmics'], binwidth=self.timeseries.scale*0.2*self.timeseries.exposurenoise, alpha=0.5, color=self.plotting['nocosmics'])
        self.plothistogram(self.timeseries.scale*self.binned['naive'], binwidth=self.timeseries.scale*0.2*self.timeseries.exposurenoise, alpha=0.5, color=self.plotting['naive'])
        self.plothistogram(self.timeseries.scale*self.binned['flux'], binwidth=self.timeseries.scale*0.2*self.timeseries.exposurenoise, alpha=0.5, color=self.plotting['flux'])

        # plot the image and position, if this isn't just a toy model
        if includeimage:
          self.plotimage(ax=self.ax['image'])

        # labels describing the simulated timeseries
        self.ax['timeseries'].text(self.plotting['toplot']*0.02, y, "{rate:0.2f} cosmic rays per exposure, at {amplitude:0.1f}X noise".format( amplitude=self.timeseries.cosmicsamplitude, rate=self.timeseries.cosmicsperexposure), fontsize=6, color=self.plotting['flux'], alpha=0.7, horizontalalignment='left')
        self.ax['timeseries'].text(self.plotting['toplot']*0.98, y, "{nsubexposures} subexposures".format( nsubexposures=self.timeseries.nsubexposures), fontsize=6, color=self.plotting['flux'], alpha=0.7, horizontalalignment='right')

        # labels showing the unmitigated, achieved, and ideal noise values
        left, right = np.log(self.ax['histbinned'].get_xlim())
        span = right - left
        self.ax['histbinned'].text(np.exp(span*0.2 + left), y, "{:.2f}".format(self.timeseries.scale*self.unmititigated/self.timeseries.exposurenoise), fontsize=6, color=self.plotting['naive'], horizontalalignment='center', alpha=0.7)
        self.ax['histbinned'].text(np.exp(span*0.5 + left), y, "{:.2f}".format(self.timeseries.scale*self.achieved/self.timeseries.exposurenoise), fontsize=6, color=self.plotting['flux'], horizontalalignment='center', alpha=0.7)
        self.ax['histbinned'].text(np.exp(span*0.78 + left), y, "{:.2f}".format(self.timeseries.scale), fontsize=6, color=self.plotting['nocosmics'], horizontalalignment='center', alpha=0.7)
        self.ax['histbinned'].text(np.exp(span*0.2 + left), ylevel - ywidth*1.1, 'unmitigated', fontsize=4, color=self.plotting['naive'], horizontalalignment='center', alpha=0.7)
        self.ax['histbinned'].text(np.exp(span*0.5 + left), ylevel - ywidth*1.1, 'achieved', fontsize=4, color=self.plotting['flux'], horizontalalignment='center', alpha=0.7)
        self.ax['histbinned'].text(np.exp(span*0.8 + left), ylevel - ywidth*1.1, 'perfect', fontsize=4, color=self.plotting['nocosmics'], horizontalalignment='center', alpha=0.7)

        # draw and save the plot
        plt.show();
        plt.savefig(self.fileprefix() + '.pdf')
        print("saved plot to " +  self.fileprefix() + '.pdf')
