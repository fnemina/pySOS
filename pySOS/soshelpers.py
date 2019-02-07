# coding=utf-8

import numpy as np
from scipy.interpolate import interp1d


def RunWavelengths(s, wavelengths=[0.550], angle=0, output="I"):
    """ This method run the simulation for a given pySOS object for a set of
        wavelengths and angles and returns the output from the file up.

        Parameters
        ----------
        s           The pySOS object for which we want to run the simulation.
        wavelengths An interable with the wavelengths in micrometers for which
                    to run the simulations
        angles      The view angle in degrees for which to run the simulation.
        output      The name of the output we want to compute as an string
                    I,Q,U       Stokes parameter at output level Z (in sr-1)
                                normalised to the extraterrestrial solar
                                irradiance
                                (PI * L(z) / Esun)
        """

    if output not in ["I", "Q", "U"]:
        raise(ValueError("Wrong output variable."))

    values = np.array([])

    if type(angle) is int or type(angle) is float:
        angle = np.zeros(np.size(wavelengths))+angle

    for idx, wl in np.ndenumerate(wavelengths):
        # We set the wavelength and run the simulation
        s.wa = wl
        s.run()
        # Convert the output to a directory
        results = vars(s.outputs.up)
        # We interpolte the values and add it to a numpy array
        f = interp1d(results['ang'], results[output])
        values = np.append(values, f(angle[idx[0]]))

    return values
