"""Module for mixing functionality.

Provides classes that calculate new proposal/trial values from the history of
previously seen results in ways that accelerate the convergence of the loop to
its fixed point.

Apart from classes implementing different mixing algorithms,
- LinearMixer for linear mixing of the previous result
- DiisMixer for Pulay mixing/DIIS of the result history
- NoMixingMixer for no mixing
there are also the decorator classes InitialMixingDecorator, which
allows switching mixers after a number of iterations, as well as
FlatMixingDecorator and RealMixingDecorator, which resp. flatten the
input and separate its real and imaginary parts before passing it on
to the decorated class, which mainly serves the purpose of bringing
more convoluted input quantities into a form suitable for use with
DiisMixer.

All mixer objects are callable, taking the new values as argument and
returning the next proposal/trial value.
"""
import numpy as np

import w2dyn.auxiliaries.deepflatten as deepflatten

class InitialMixingDecorator(object):
    """
    This mixing decorator is switching from the `init_mixer` to the `mixer`
    after `init_count` calls. Useful if `mixer` has unwanted behavior at the
    beginning.
    """
    def __init__(self, init_count, init_mixer, mixer):
        self.init_count = init_count
        self.init_mixer = init_mixer
        self.mixer = mixer
    def __call__(self, *args):
        if self.init_count > 0:
            self.init_count -= 1
            return self.init_mixer(*args)
        return self.mixer(*args)

class FlatMixingDecorator(object):
    """
    This mixing decorator takes any kind of nested python lists with numbers and
    numpy arrays and calls `mixer` with a one dimensional copy of the data. The
    shape is restored after `mixer` returned.
    """
    def __init__(self, mixer):
        self.mixer = mixer
    def __call__(self, *args):
        if len(args) == 1: args = args[0]
        types = deepflatten.types(args)
        shape = deepflatten.shapes(args)
        x = deepflatten.flatten(args)
        x = self.mixer(x)
        x = deepflatten.restore(x, shape, types)
        return x

class RealMixingDecorator(object):
    """
    This mixing decorator takes an array or list as input and concatenates its
    real and imaginary parts to call `mixer`. A numpy array with complex numbers
    is retruned.
    """
    def __init__(self, mixer):
        self.mixer = mixer
    def __call__(self, x):
        n = x.shape[0]
        x = np.concatenate([np.real(x), np.imag(x)])
        x = self.mixer(x)
        x = x[:n] + 1j*x[n:]
        return x

class NoMixingMixer(object):
    """
    This mixer is just passing the input through and therefore not mixing at
    all.
    """
    def __call__(self, *args):
        return args

class LinearMixer(object):
    """
    Allows (linear) under relaxation of quantities in the DMFT loop.

    This is achieved by mixing into the new self-energy a certain share of the
    previous self-energy (controlled by `oldshare`) every time that `mix()` is
    called:

            A^{mixed}_n = (1-oldshare) * A_n + oldshare * A^{mixed}_{n-1}

    thereby exponentially decreasing the influence of the old iterations by
    `\exp(-n \log oldshare)`. This strategy dampens strong statistical
    fluctuations in the QMC solver and ill-defined chemical potentials in
    insulating cases.
    """
    def __init__(self, old_share=0):
        self.old_share = float(old_share)
        self.old_value = None

    def __call__(self, new_value):
        if self.old_value is None:
            new_trial = new_value
        else:
            new_trial = self.old_share * self.old_value + (1 - self.old_share) * new_value
        self.old_value = new_trial
        return new_trial

class DiisMixer(object):
    """
    This mixing algorithm, also known as Pulay mixing, uses multiple trials and
    results to combine a (hopefully) better trial for the next iteration.
    
    Futher reading
     - https://en.wikipedia.org/wiki/DIIS
     - P. Pulay, Chem. Phys. Lett. 73, 393
       "Convergence acceleration of iterative sequences.
        The case of SCF iteration"
       https://dx.doi.org/10.1016/0009-2614(80)80396-4
     - P. Pulay, J. Comp. Chem. 3, 556
       "Improved SCF convergence acceleration"
       https://dx.doi.org/10.1002/jcc.540030413
     - P. P. Pratapa, P. Suryanarayana, Chem. Phys. Lett. 635, 69
       "Restarted Pulay mixing for efficient and robust acceleration of
        fixed-point iterations"
       https://dx.doi.org/10.1016/j.cplett.2015.06.029
     - A. S. Banerjee, P. Suryanarayana, J. E. Pask, Chem. Phys. Lett. 647, 31
       "Periodic Pulay method for robust and efficient convergence acceleration
        of self-consistent field iterations"
       https://dx.doi.org/10.1016/j.cplett.2016.01.033
     - H. F. Walker, P. Ni, SIAM J. Numer. Anal., 49, 1715
       "Anderson Acceleration for Fixed-Point Iterations"
       https://dx.doi.org/10.1137/10078356X
       https://core.ac.uk/display/47187107
    """
    def __init__(self, old_share, history, period):
        self.alpha = 1 - old_share
        self.history = history
        self.period = period

        self.i = 0
        self.trials = []
        self.residuals = []

    def __call__(self, new_value):
        if self.i <= 0:
            # no history yet
            new_trial = new_value
        else:
            trial = self.trials[-1]
            residual = new_value - trial
            self.residuals.append(residual)

            # trim history
            self.trials = self.trials[-self.history:]
            self.residuals = self.residuals[-self.history:]

            if self.i <= 1 or (self.i % self.period) != 0:
                # linear mixing
                new_trial = trial + self.alpha * residual
            else:
                # Pulay mixing
                R = np.array(self.trials); R = R[1:] - R[:-1]; R = R.T
                F = np.array(self.residuals); F = F[1:] - F[:-1]; F = F.T

                new_trial = trial + self.alpha * residual \
                         - np.linalg.multi_dot([R + self.alpha * F, np.linalg.pinv(F), residual])

        self.i += 1
        self.trials.append(new_trial)
        return new_trial

