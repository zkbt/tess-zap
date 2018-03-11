'''
A collection of example models that can be fed into zap-cosmics.
'''

from .imports import *

class Model():
    '''
    A base model, defining a handy plot tool.
    '''
    def plot(self, tmin=-0.5, tmax=0.5, n=1000, **plotkw):
        t = np.linspace(tmin, tmax, n)
        plt.plot(t, self(t), label='{}'.format(self), **plotkw)

class Flare(Model):
    '''
    A fast-rising exponential decay.
    '''
    def __init__(self, start=0.02, decay=0.1, amplitude=0.5):
        self.start = start
        self.decay = decay
        self.amplitude = amplitude

    def __call__(self, t):
        m = np.zeros_like(t)
        afterstart = t > self.start
        m[afterstart] += self.amplitude*np.exp(-(t[afterstart]-self.start)/self.decay)
        return m

    def __repr__(arg):
        return '<flare with {:.2}*exp(-t/{:.2})>'.format(self.ampltiude, self.decay)

class ManyFlares(Model):
    '''
    Create lots of flares, all stacked together.
    '''
    def __init__(self, N=25):
        self.N = N
        starts = np.random.uniform(-0.5, 0.5, N)
        decays = np.random.uniform(0, 0.05, N)
        amplitudes = np.random.uniform(0, 0.5, N)

        self.flares = []
        for s, d, a in zip(starts, decays, amplitudes):
            self.flares.append(Flare(start=s, decay=d, amplitude=a))

    def __call__(self, t):
        m = np.ones_like(t)
        for f in self.flares:
            m += f(t)
        return m

    def __repr__(self):
        return '<{} random flares>'.format(self.N)

class Transit(Model):
    '''
    A transit model.
    '''
    def __init__(self, epoch=0.0,
                       period=1.0,
                       rp_over_rs=0.1,
                       a_over_rs=15.0,
                       impact_parameter=0.5,
                       eccentricity=0.0,
                       omega=90.0,
                       limb_coef=[0.1, 0.3],
                       limb_type="quadratic"
                       ):
        '''
        Initialize a transit object and set its parameters.
        '''
        import batman
        self.batman = batman
        self.params = batman.TransitParams()
        self.params.t0 = epoch
        self.params.per = period
        self.params.rp = rp_over_rs
        self.params.a = a_over_rs
        inclination = np.arccos(impact_parameter/a_over_rs)*180.0/np.pi
        self.params.inc = inclination
        self.params.ecc = eccentricity
        self.params.w = omega
        self.params.u = limb_coef
        self.params.limb_dark = limb_type


    def __call__(self, t):
        try:
            assert(np.all(self.t == self.batmant))
        except (AssertionError, AttributeError):
            self.batmant = t
            self.batmanmodel = self.batman.TransitModel(self.params, t)

        return self.batmanmodel.light_curve(self.params)

    def __repr__(self):
        return '<transit with rp_over_rs={:.3}>'.format(self.params.rp)

#def test():
if __name__ == '__main__':
    constant = np.ones
    polynomial = np.poly1d([1,0,0.02])

    t_deep = Transit(rp_over_rs=0.3)
    t_shallow = Transit(rp_over_rs=0.05)

    few = ManyFlares(N=2)
    many = ManyFlares(N=25)

    for model in [t_deep, t_shallow, few, many]:
        model.plot()

    plt.legend(bbox_to_anchor=(1,1), loc='upper left', frameon=False)
