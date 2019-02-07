# coding=utf-8

import numpy as np
import os
from scipy.io import FortranFile


def ExtractValue(text, reference):
    """ This function extracts the value for a given reference string for the
        text string given and returns it as a float number
        text        Text to search
        reference   Reference string
        """
    for a in text:
        if reference in a:
            res = a.replace(reference, "")

    return float(res)

class PM(object):
    """ Files containing the radiative properties of Aerosols, chlorophyll or
        Mineral-Like particles
        """

    def __init__(self, resroot, filename):
        """ First thirteen header lines provide comments and formatted data on
            the extinction and scattering cross-sections (in μm 2 ), the
            asymmetry factor, the volume of the equivalent mean particle
            (in μm 3 ), the real part of the mean refractive index, the phase
            function truncation coefficient and the single scattering albedo
            (adjusted to the truncation).

            The following lines contain the phase matrix coefficients of the
            development of the Legendre Polynomials of the phase matrix for
            each order k ranging from k= 0 up to the maximum order of
            computations (OS_NB).

            resroot     OSOAA results root directory.
            filename    Filename to look for the results.
            fulltext    Full file text.
            extcs       Extinction cross section (mic^2)
            scacs       Scattering cross section (mic^2)
            asymm       Asymmetry factor (no truncation)
            mpd         Mean particules altitude/depth (m)
            vol         Volume of the mean particule (mic^3)
            rindex      Mean refractive index (real part)
            trunca      Truncation coefficient
            singlesca   Single scattering albedo (truncation)

            alpha       the coefficient alpha Related to the polarized phase
                        functions
            beta        the coefficient beta Related to the polarized phase
                        functions
            gamma       the coefficient gamma Related to the polarized phase
                        functions
            xi          the coefficient xi Related to the polarized phase
                        functions

            These coefficients are adjusted to a phase function truncation if
            applied.
            """
        # We open the file with the corresponding encoding and convert it
        # to a text string.
        with open(resroot+"/SOS/"+filename,
                  encoding="iso-8859-15") as file:
            self.fulltext = file.readlines()
        # Read variables not tabulated
        self.extcs = ExtractValue(self.fulltext,
                                    "EXTINCTION CROSS SECTION (mic^2)     :")
        self.scacs = ExtractValue(self.fulltext,
                                    "SCATTERING CROSS SECTION (mic^2)     :")
        self.asymm = ExtractValue(self.fulltext,
                                    "ASYMMETRY FACTOR (no truncation)     :")
        self.trunca = ExtractValue(self.fulltext,
                                    "TRUNCATION COEFFICIENT               :")
        self.singlesca = ExtractValue(self.fulltext,
                                    "SINGLE SCATTERING ALBEDO (truncation): ")

        # Get header length to skip it
        skipheader = [idx for idx, text in enumerate(self.fulltext)
                      if "ALPHA(K)        BETA11(K)" in text][0]
        self.alpha, self.beta, self.gamma, self.xi = np.genfromtxt(resroot+"/SOS/"+filename,
                                                                   skip_header=skipheader+1, unpack=True,
                                                                   encoding="iso-8859-15")

class RADIANCE(object):
    """ Read the normalized radiances
        """

    def __init__(self, resroot, filename):
        """ Read the normalized radiacnes

            resroot     OSOAA results root directory.
            filename    Filename to look for the results.


            These coefficients are adjusted to a phase function truncation if
            applied.
            """

        self.ang, self.I, self.Q, self.U = np.genfromtxt(resroot+"/SOS/"+filename,
                                                         unpack=True,
                                                         encoding="iso-8859-15")

class TRASMITANCE(object):
    """ Read the transmitances"""

    def __init__(self, resroot, filename):
        """ Read the normalized radiacnes

            resroot     OSOAA results root directory.
            filename    Filename to look for the results.


            These coefficients are adjusted to a phase function truncation if
            applied.
            """
        with open(resroot+"/SOS/"+filename,
                  encoding="iso-8859-15") as file:
            self.fulltext = file.readlines()

        self.solang = ExtractValue(self.fulltext,
                                    "Solar Zenithal Angle  :")
        self.Tdown = ExtractValue(self.fulltext,
                                    "Direct transmission  TOA -> surface :")
        self.tdownang = float(self.fulltext[4].split()[2])
        self.tdownval = float(self.fulltext[4].split()[5])

        ang = np.array([])
        trans = np.array([])
        for s in self.fulltext[7:]:
            ang = np.append(ang, float(s.split()[2]))
            trans = np.append(trans, float(s.split()[5]))

        self.tang = ang
        self.ttrans = trans


class OUTPUTS(object):
    """ This class contains the standard and advanced outputs generated by the
        OSOAA software"""

    def __init__(self, resroot, filenames):
        """ This methods inits the output class with all the avalaible outputs.
            resroot     Results root Directory
            filenames   Object with all filenames
            """

        self.aer = PM(resroot, filenames.aer)
        self.up = RADIANCE(resroot, filenames.advup)
        self.down = RADIANCE(resroot, filenames.advdown)
        self.trans = TRASMITANCE(resroot, filenames.trans)
