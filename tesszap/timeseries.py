'''
Create some 1D timeseries, either with a given S/N or as estimate for a given
TESS magnitude (using a very approximate noise model). For the purposes
of thei
'''

from .imports import *

class TimeseriesWithoutModel():
	"""
	Implement a timeseries that keeps track of both the individual subcadences
	and the actual cadences that have been binned together by a stacker.
	"""

	def __init__(self,  subcadencetime,
						subcadenceflux,
						subcadenceuncertainty,
						subcadence=2.0,
						cadence=120.0,
				):
		"""
		Initialize the TimeseriesWithoutModel object.

		Parameters
		----------
		subcadencetime : 1D array
			Subcadence light curve times, in *days*.
		subcadenceflux : 1D array
			Subcadence light curve flux values.
		subcadenceuncertainty : float
			Subcadence light curve flux uncertainties (currently accepts only a single value).
		subcadence : float, optional
			The duration of each subcadence point, in *seconds*.
		cadence : float, optionl
			The duration of actual cadences, to which the subcadences will be binned, in *seconds*.
		"""

		# this will be updated if its a real
		self.model = False
		self.cosmics = False
		self.unbinned, self.binned = {}, {}


		# how many subcadences will be stacked together?
		self.subcadence = subcadence
		self.cadence = cadence
		if self.cadence % self.subcadence != 0:
			raise ValueError("{} is not evenly divisible by {}".format(self.cadence, self.subcadence))
		self.nsubcadences = np.int(self.cadence/self.subcadence)
		self.ncadence = np.int(len(subcadencetime)/self.nsubcadences)

		# store the subcadence timeseries
		self.unbinned['time'] = subcadencetime
		self.unbinned['flux'] = subcadenceflux
		self.unbinned['unmitigated'] = self.unbinned['flux']

		self.subcadenceuncertainty = subcadenceuncertainty
		self.cadenceuncertainty = self.subcadenceuncertainty/np.sqrt(self.nsubcadences)


	def meancollapse(self, x):
		'''
		Use a simple mean to collapse from subcadences to cadences.

		Returns
		-------

		stackedmean : array
			A stacked array, with the same dimensions as the total number of times.
		'''

		return np.mean(self.twod(x), 1)

	def twod(self, oned):
		'''
		Convert the 1D subcadence arrays into a 2D array, for math.

		Returns
		-------

		reshaped : array
			Simply an array that's been trimmed and reshaped to make it easier to collapse.
		'''
		return oned[:self.ncadence*self.nsubcadences].reshape(self.shape)

	@property
	def shape(self):
		'''
		The shape of a 2D array, with ncadences x nsubcadences

		Returns
		-------

		shape : array
			The shape of a collapsable array.
		'''
		return (self.ncadence, self.nsubcadences)

	def stack(self, strategy):
		'''
		Stack the subcadences into binned cadences, using a given stacking strategy.

		Parameters
		----------

		strategy : a Strategy object
			This will be used for collapsing the subcadences into binned cadences.

		Returns
		-------

		binned : dict
			The cadence-binned light curve, as a dictionary.
		'''

		# store the stacker that's being used
		self.strategy = strategy

		# calculate the binned times (the straight mean assumes we don't know which time was lost)
		self.binned['time'] = self.meancollapse(self.unbinned['time'])

		# apply the stacking filter
		self.binned['flux'] = strategy(self)

		# naive timeseries is binned using a simple mean
		self.binned['unmitigated'] = self.meancollapse(self.unbinned['unmitigated'])

		if self.cosmics:
			self.binned['cosmics'] = self.meancollapse(self.unbinned['cosmics'])
			self.binned['nocosmics'] = self.binned['unmitigated'] - self.binned['cosmics']

		if self.model:
			self.binned['model'] = self.meancollapse(self.unbinned['model'])
		else:
			self.binned['model'] = np.mean(self.binned['flux'])

		# calculate RMS of the binned light curves, with and without the cosmic ray filter
		self.rms = {}
		self.rms['unmitigated'] = np.std(self.binned['unmitigated'] - self.binned['model'])
		self.rms['achieved'] = np.std(self.binned['flux'] - self.binned['model'])
		self.rms['expected'] = self.subcadenceuncertainty/np.sqrt(self.nsubcadences)

		return self.binned











class Timeseries(TimeseriesWithoutModel):
	def __init__(self,  tmin=-0.5, tmax=0.5,
						model=np.ones_like,
						subcadenceuncertainty=0.01,
						subcadence=2.0,
						cadence=120.0,
						addcosmics=True,
						cosmickw=dict(height=1.0, # what's the height of a single cosmic ray?
									  probability=0.001), # what's the probability a subexposure is hit with a cosmic?
				):
		"""
		Initialize the TimeseriesWithoutModel where we know exactly the model that created it

		Parameters
		----------
		tmin : float
			The minimum time, in *days*.
		tmax : float
			The maximum time, in *days*.
		model : function
			A callable object, that takes time as an input and returns model fluxes.
		subcadenceuncertainty : float
			Subcadence light curve flux uncertainties (currently accepts only a single value).
		subcadence : float, optional
			The duration of each subcadence point, in *seconds*.
		cadence : float, optionl
			The duration of actual cadences, to which the subcadences will be binned, in *seconds*.
		addcosmics : bool
			Should we inject some cosmic rays into this timeseries?
		cosmickw : dict
			Keyword arguments that will be passed on to `self.addCosmics(**cosmic)` to define how cosmics should be added.
		"""

		# create a new grid of times
		dtcadence = cadence/24.0/60.0/60.0
		dtsubcadence = subcadence/24.0/60.0/60.0

		cadencetime = np.arange(tmin, tmax, dtcadence)
		nsubcadences = np.int(cadence/subcadence)
		asinglesubcadenceset = np.arange(-dtcadence/2.0, dtcadence/2.0, dtsubcadence)
		subcadencetime = (cadencetime[:, np.newaxis] + asinglesubcadenceset[np.newaxis, :]).flatten()

		# calculate the model
		m =  model(subcadencetime)

		# generate a noise realization
		subcadenceflux = np.random.normal(m, subcadenceuncertainty)

		TimeseriesWithoutModel.__init__(self, subcadencetime,
											subcadenceflux,
											subcadenceuncertainty,
											subcadence,
											cadence)
		self.unbinned['model'] = m
		self.unbinned['residuals'] = subcadenceflux - m

		self.model = model
		if addcosmics:
			self.addCosmics(**cosmickw)

	def addCosmics(self, 	height=1,
							probability=0.001):
		'''
		Add simulated (single-subcadence) cosmic rays into the flux.

		Parameters
		----------
		height : float
			How big will a single cosmic ray hit appear at the binned cadence?
				(this in units of the *per-cadence* RMS)
		probability : float
			How many cosmic ray hits per subcadence?
		'''

		# store the cosmic parameters
		self.cosmics = dict(height=height, probability=probability)

		# how big will the cosmic ray be in the cadences
		diluted = height*self.cadenceuncertainty

		# how big must it have therefore been in the subcadences?
		undiluted = diluted*self.nsubcadences

		# inject the cosmics into the flux timeseries
		N = len(self.unbinned['unmitigated'])
		self.unbinned['cosmics'] =  undiluted*np.random.poisson(probability, N)
		self.unbinned['nocosmics'] = self.unbinned['unmitigated']
		self.unbinned['unmitigated'] =  self.unbinned['nocosmics'] + self.unbinned['cosmics']








	def __str__(self):
		'''Summary of this object.'''
		return "{nexposures} exposures, {nsubexposures} subexposures, cosmic rays {cosmicsamplitude:.2}X noise".format(**self.__dict__)


	def plot(self, xlim=[None, None], unbinned=False):

		# set up the plots
		self.figure = plt.figure('{}'.format(self.model), figsize=(8, 6), dpi=100)
		gs = plt.matplotlib.gridspec.GridSpec(2, 2,
					width_ratios=[1.0, 0.2], height_ratios=[1.0, 0.4],
					wspace=0.0, hspace=0.02)
		self.ax = {}
		self.ax['flux'] = plt.subplot(gs[0,0])
		plt.setp(self.ax['flux'].get_xticklabels(), visible=False)
		plt.ylabel('Relative Flux')
		self.ax['residual'] = plt.subplot(gs[1,0], sharex=self.ax['flux'])
		plt.xlabel('Time (days)')
		plt.ylabel('Residuals\n(in $\sigma$)')
		self.ax['histogramresidual'] = plt.subplot(gs[1,1], sharey=self.ax['residual'])
		self.ax['histogramresidual'].axis('off')

		# set up all of the plotting parameters
		binnedkw = dict(marker='o', markersize=5, alpha=0.5, linewidth=1, markeredgecolor='none')
		unbinnedkw = dict(marker='.', markersize=0.1, alpha=0.3, linewidth=0, markeredgecolor='none')
		histogramkw = dict(linewidth=2, alpha=0.5)
		typekw = {	'flux':dict(color='mediumvioletred', zorder=3),
					'nocosmics':dict(color='darkorange', zorder=1),
					'unmitigated':dict(color='royalblue', zorder=2)}
		modelkw = dict(linewidth=2, color='gray', zorder=-1)

		def plotTimeseries(x, y, **kw):
			'''Helper to plot a timeseries.'''
			if xlim[0]:
				ok = (x >= np.min(xlim)) * (x <= np.max(xlim))
			else:
				ok = np.ones_like(y).astype(np.bool)
			plt.plot(x[ok], y[ok], **kw)

		### plot the actual time series
		plt.sca(self.ax['flux'])
		for k in typekw.keys():
			if unbinned:
				unbinnedkw.update(typekw[k])
				plotTimeseries(self.unbinned['time'], self.unbinned[k], **unbinnedkw)
			binnedkw.update(typekw[k])
			plotTimeseries(self.binned['time'], self.binned[k], **binnedkw)
		plt.plot(self.unbinned['time'], self.unbinned['model'], **modelkw)

		# labels describing the simulated timeseries
		cosmiclabel = "{:0.2f} cosmic rays per subexposure, at {:0.1f}X noise".format(
			self.cosmics['probability']*self.cadence/self.subcadence, self.cosmics['height'])
		plt.text(0, 0, cosmiclabel, transform=plt.gca().transAxes, color='gray', fontsize=4, horizontalalignment='left', alpha=0.7)
		cadencelabel = "{cadence:.0f}s cadence, {subcadence:.0f}s subcadence".format(**self.__dict__)
		plt.text(1, 0, cadencelabel, transform=plt.gca().transAxes, color='gray', fontsize=4, horizontalalignment='right', alpha=0.7)


		# set the scale for plotting the residuals
		scale = self.cadenceuncertainty
		binwidth = np.minimum(100*5.0/self.ncadence, 0.5)

		### plot the residuals
		plt.sca(self.ax['residual'])
		for k in typekw.keys():
			if unbinned:
				unbinnedkw.update(typekw[k])
				plotTimeseries(self.unbinned['time'], (self.unbinned[k] - self.unbinned['model'])/scale, **unbinnedkw)
			binnedkw.update(typekw[k])
			plotTimeseries(self.binned['time'], (self.binned[k] - self.binned['model'])/scale, **binnedkw)
		plt.axhline(0, xmin=np.min(self.unbinned['time']), xmax=np.max(self.unbinned['time']), **modelkw)
		plt.xlim(*xlim)

		# plot the histograms
		plt.sca(self.ax['histogramresidual'])
		def plotHistogram(y, **kwargs):
			'''
			Helper function to plot a histogram,
			rotated clockwise by 90 degrees,
			to represent a projection of a timeseries plot.
			'''

			# create a histogram of the lightcurve values
			b = np.arange(np.min(y) - binwidth, np.max(y) + binwidth, binwidth)
			yhist, edges = np.histogram(y, bins=b, density=True)

			# define the "x"-axis at the centers of the histogram bins
			xhist = (edges[1:] + edges[0:-1])/2.0

			# plot in the histogram panel
			plt.plot(yhist, xhist, **kwargs)
			#plt.xlim(3, np.max(yhist)*1.5)
		for k in typekw.keys():
			histogramkw.update(typekw[k])
			plotHistogram((self.binned[k] - self.binned['model'])/scale, **histogramkw)
		# overplot a model histogram on top of that
		y = np.linspace(*plt.ylim(), num=200)
		plt.plot(np.exp(-0.5*y**2)/np.sqrt(2*np.pi), y, **modelkw)
		plt.xscale('log')
		plt.xlim(2.0/len(self.binned['flux']), None)


		keys = ['expected',  'achieved', 'unmitigated' ]
		kw = [typekw['nocosmics'], typekw['flux'], typekw['unmitigated']]
		textkw = dict(transform=self.ax['histogramresidual'].transAxes,
						horizontalalignment='center', alpha=0.7)

		for i, k in enumerate(keys):
			textkw.update(kw[i])
			plt.text(0.5, 1.3 + i*0.4, k, fontsize=9, **textkw)
			ratio = '{:.2f}'.format(self.rms[k]/self.rms['expected'])
			plt.text(0.5, 1.3 + i*0.4 - 0.12, ratio, fontsize=11, **textkw)
			actual = '({:.0f} ppm)'.format( 1e6*self.rms[k])
			plt.text(0.5, 1.3 + i*0.4 - 0.2, actual, fontsize=8, **textkw)

	'''
	self.ax['histbinned'].text(np.exp(span*0.2 + left), y, "{:.2f}".format(self.timeseries.scale*self.unmititigated/self.timeseries.exposurenoise), fontsize=6, color=self.plotting['naive'], horizontalalignment='center', alpha=0.7)
	self.ax['histbinned'].text(np.exp(span*0.5 + left), y, "{:.2f}".format(self.timeseries.scale*self.achieved/self.timeseries.exposurenoise), fontsize=6, color=self.plotting['flux'], horizontalalignment='center', alpha=0.7)
	self.ax['histbinned'].text(np.exp(span*0.78 + left), y, "{:.2f}".format(self.timeseries.scale), fontsize=6, color=self.plotting['nocosmics'], horizontalalignment='center', alpha=0.7)
	self.ax['histbinned'].text(np.exp(span*0.2 + left), ylevel - ywidth*1.1, 'unmitigated', fontsize=4, color=self.plotting['naive'], horizontalalignment='center', alpha=0.7)
	self.ax['histbinned'].text(np.exp(span*0.5 + left), ylevel - ywidth*1.1, 'achieved', fontsize=4, color=self.plotting['flux'], horizontalalignment='center', alpha=0.7)
	self.ax['histbinned'].text(np.exp(span*0.8 + left), ylevel - ywidth*1.1, 'perfect', fontsize=4, color=self.plotting['nocosmics'], horizontalalignment='center', alpha=0.7)
	'''



		# fiddle with the ylimits of the plots, depending on whether looking at binned or unbinned timeseries
		#nsigma = 5
		#if unbinned:
		#	plt.ylim(None,  None)
		#else:
		#	plt.ylim(-nsigma, nsigma)


	def movie(self, window=1.0, fps=20, bitrate=1800*10, filename='test.mp4'):

		# make a plot with no xlim
		self.plot(xlim=[None, None])
		metadata = dict(title=None, artist=None)
		self.writer = ani.FFMpegWriter(fps=fps, metadata=metadata, bitrate=bitrate)
		with self.writer.saving(self.figure, filename, self.figure.get_dpi()):

			overlap = 0.03
			xlim = plt.xlim()
			span = xlim[1] - xlim[0]
			xmin, xmax = np.min(self.unbinned['time']), np.max(self.unbinned['time'])
			N = (xmax-xmin)/(overlap*span)
			leftedges = np.linspace(xmin, xmax-span, N)

			# loop over spans
			for left in tqdm(leftedges):
				xlim = left, left+span
				self.ax['residual'].set_xlim(*xlim)
				self.writer.grab_frame()

		plt.close()

"""
class Timeseries1D(Timeseries):
	'''1D timeseries, set by a S/N and an amplitude of cosmic rays.'''
	def create(self, 	nexposures=1324,
						nsubexposures=900,
						subexposurecadence = 2.0,
						snr=100,
						amplitude=1.0,
						probabilityofcosmic = 0.001,
						):
		'''
		Initialize a 1D toy model timeseries, given a S/N and an amplitude and rate for cosmics.

		The parameters for a toy model timeseries are:
			nexposures=[how many binned exposures?]
			nsubexposures=[how many subexposures within each exposure?]
			subexposurecadence=[what's the cadence, in seconds, of the subexposures?]
			snr=[what's the S/N of the timeseries?]
			amplitude=[what's the amplitude of cosmic rays relative to binned exposure noise?],
			probabilityofcosmic=[what's the probability of a cosmic ray hit in a given sub-exposure?],
		'''

		# what's the total number of binned exposures (e.g., how many half-hour FFIs?)
		self.nexposures = nexposures

		# what's the number of subexposures per exposure?
		self.nsubexposures = nsubexposures

		# what's the cadence (in real seconds) of those subexposures?
		self.subexposurecadence = subexposurecadence

		# what's the S/N of the subexposure timeseries? (the flux will default to 1)
		self.snr = snr

		# what's the noise of a subexposure
		self.subexposurenoise = 1.0/self.snr

		# is this the right way to set these up???
		self.exposurenoise = self.subexposurenoise/np.sqrt(self.nsubexposures)

		# the cosmic ray amplitude is relative to the binned exposure noise
		self.cosmicsamplitude = amplitude

		# what's the rate of cosmic rays per binned exposure?
		self.cosmicspersubexposure = probabilityofcosmic
		self.cosmicsperexposure = self.cosmicspersubexposure*self.nsubexposures #1.0 - (1.0 - probabilityofhit)**self.nsubexposures
		self.exposurecadence = self.nsubexposures*self.subexposurecadence

		# then create a simple timeseries
		self.createSimple()

	def createSimple(self, cosmics=True, noise=True):
		'''Populate this timeseries with a simple toy model.'''

		# create an array to store flux values
		self.flux = np.ones(self.shape)

		# add noise to the timeseries
		if noise:
			self.addNoise()

		# add cosmic rays to the timeseries
		if cosmics:
			self.addCosmics()

		# define an axis of cadence numbers along the timeseries
		self.x = np.arange(0, self.nexposures, 1.0/self.nsubexposures).reshape(self.shape)

		# note that this is a toy model timeseries
		self.toy = True

	def addNoise(self):
		'''For toy model, include Gaussian noise.'''

		# make a noiseless timeseries
		self.noiseless = self.flux + 0.0

		# add a noise realization to each subexposure
		self.flux += np.random.normal(0,self.subexposurenoise,self.shape)

	def addCosmics(self):
		'''For toy model, include cosmic rays as random impulses with fixed amplitude.'''

		# a flux array with amplitude of the cosmic ray amplitude times the number of cosmics per subexposure
		self.cosmics = self.cosmicsamplitude*self.exposurenoise*self.nsubexposures*np.random.poisson(self.cosmicspersubexposure, self.shape)

		# add these cosmics into the timeseries
		self.flux += self.cosmics


class Timeseries1DTESS(Timeseries1D):
	'''1D cartoon timeseries, simulating a single TESS pixel.'''

	def __init__(self, 	nexposures=1324,
						nsubexposures=900,
						probabilityofcosmic=0.001,
						magnitude=10.0,
						subexposurecadence = 2.0,
						):
		'''
		Initialize a 1D toy model timeseries for a single TESS pixel timeseries,
		given a magnitude for the star, and the photons/cosmic ray.

		The parameters for a toy model timeseries are:
			nexposures=[how many binned exposures?]
			nsubexposures=[how many subexposures within each exposure?]
			probabilityofcosmic=[what's the probability of a cosmic ray hit in a given sub-exposure?],
			magnitude=[the brightness of the star sets the snr of the timeseries (and the amplitude)]
			)
		'''


		# what's the TESS magnitude of the star?
		self.magnitude = magnitude

		# what fraction of the stellar light is contained in that single pixel?
		self.containedflux = 0.3

		# how many stellar photons land on this pixel (per subexposure)?
		self.photonsfromstar = 1.7e6*subexposurecadence*73.0*10**(-0.4*self.magnitude)*self.containedflux

		# what's the noise from the sky, and from readout?
		skynoise = 10.0
		readoutnoise = 10.0

		# total number of stellar photons on the pixel in an exposure (in photoelectrons)
		signal = nsubexposures*self.photonsfromstar

		# calculate the noise (in photoelectrons)
		noise = np.sqrt(nsubexposures*(readoutnoise**2 + skynoise**2 + self.photonsfromstar))
		snr = signal/noise

		# calculate the cosmic ray impact, fractionally, compared to the binned exposure
		photonspercosmic=1000.0

		# this is assuming all cosmic rays are the same in their amplitude
		amplitude = photonspercosmic/noise

		Timeseries1D.__init__(self, nexposures=nexposures,
									nsubexposures=nsubexposures,
									subexposurecadence=subexposurecadence,
									snr=snr,
									amplitude=amplitude,
									probabilityofcosmic=probabilityofcosmic)
"""
