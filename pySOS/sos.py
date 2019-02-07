# coding=utf-8

import os
import random
import string
from io import open


class ANG(object):
    """ Angle definitions class."""

    class ANGLES(object):
        """ Angle class to use within object"""

        def __init__(self, nbgauss, userangfile):
            """ Init of the angles class
                nbgaus      Number of gauss angles
                userangfile Filename of the complementary list of user's angles
                """

            self.nbgauss = nbgauss
            self.userangfile = userangfile

    def __init__(self, thetas=30.0, radnb=None, raduser=None,
                 mienb=None, mieuser=None):
        """ Init the angle class
            thetas      Solar zenith angle (degrees)
            radnb       Number of Gauss angles to be used for radiance
                        computations
            raduser     Filename of the complementary list of user's angles to
                        complete the ANG.Rad.NbGauss angles (complete path).
            mienb       Number of Gauss angles to be used for mie computations
            mieuser     Filename of the complementary list of user's angles to
                        complete the ANG.mie.NbGauss angles (complete path).
            """

        self.rad = self.ANGLES(radnb, raduser)
        self.mie = self.ANGLES(mienb, mieuser)
        self.thetas = thetas


class SURFACE(object):
    """ Surface definitions class."""

    def __init__(self, type=0, alb=0):
        """ This method initiates de surface class.
            alb         Surface albedo for wavelength SOS.
            type        Type of surface condition simulated.
                        0 : Lambertian surface
                        1 : Lambertian surface and sunglint over roughness
                            ocean
                        2 : Lambertian surface and sunglint over plane water
                        3 : Lambertian surface with vegetation BRDF
                            (Roujean’s model)
                        4 : Lambertian surface with vegetation BRDF and BPDF
                            (Rondeaux’s model)
                        5 : Lambertian surface with vegetation BRDF
                            (Roujean’s model) and ground polarization
                            BPDF (Breon’s model)
                        6 : Lambertian surface with vegetation BRDF
                            (Roujean’s model) and ground
                            polarization BPDF (Nadal’s model)
            """

        self.type = type
        self.alb = alb

    def setLabertian(self, alb=0):
        """ Lambertian surface configuration
            alb         Surface albedo for wavelength SOS.
            """
        self.type = 0
        self.alb = alb
        self.wind = None
        self.ind = None
        self.file = None

        self.roujeank0 = None
        self.roujeank1 = None
        self.roujeank2 = None

        self.alpha = None
        self.beta = None

    def setLabertianGlintRough(self, alb=0, wind=5, ind=1.34, file="DEFAULT"):
        """ Lambertian surface and sunglint over roughness ocean
            alb         Surface albedo for wavelength SOS.
            wind        Wind velocity (m/s)
            ind         Refractive index air / water for the wavelength of
                        simulation SOS.Wa
            file        Pre-calculated user file (complete path)
                        or “DEFAULT” for code calculation
            """
        self.type = 0
        self.alb = alb
        self.wind = wind
        self.ind = ind
        self.file = file

        self.roujeank0 = None
        self.roujeank1 = None
        self.roujeank2 = None

        self.alpha = None
        self.beta = None

    def setLabertianGlintFlat(self, alb=0, ind=1.34):
        """ Lambertian surface and sunglint over plane water
            alb         Surface albedo for wavelength SOS.
            ind         Refractive index air / water for the wavelength of
                        simulation SOS.Wa
            """
        self.type = 1
        self.alb = alb
        self.wind = None
        self.ind = ind
        self.file = None

        self.roujeank0 = None
        self.roujeank1 = None
        self.roujeank2 = None

        self.alpha = None
        self.beta = None

    def setLabertianRoujean(self, k0, k1, k2, alb=0, file="DEFAULT"):
        """ Lambertian surface with vegetation BRDF (Roujean’s model)
            alb         Surface albedo for wavelength SOS.
            k0,k1,k2    Roujean’s model parameters for wavelength of
                        simulation SOS.Wa
            file        Pre-calculated user file (complete path)
                        or “DEFAULT” for code calculation
            """
        self.type = 3
        self.alb = alb
        self.wind = None
        self.ind = None
        self.file = file

        self.roujeank0 = k0
        self.roujeank1 = k1
        self.roujeank2 = k2

        self.alpha = None
        self.beta = None

    def setLabertianRondeaux(self, k0, k1, k2, ind, alb=0, file="DEFAULT"):
        """ Lambertian surface with vegetation BRDF and BPDF (Rondeaux’s model)
            alb         Surface albedo for wavelength SOS.
            k0,k1,k2    Roujean’s model parameters for wavelength of
                        simulation SOS.Wa
            ind         Refractive index air / ground for wavelength SOS.Wa
            file        Pre-calculated user file (complete path)
                        or “DEFAULT” for code calculation
            """
        self.type = 4
        self.alb = alb
        self.wind = None
        self.ind = ind
        self.file = file

        self.roujeank0 = k0
        self.roujeank1 = k1
        self.roujeank2 = k2

        self.alpha = None
        self.beta = None

    def setLabertianBreon(self, k0, k1, k2, ind, alb=0, file="DEFAULT"):
        """ Lambertian surface with vegetation BRDF (Roujean’s model) and
            ground polarization BPDF (Breon’s model)
            alb         Surface albedo for wavelength SOS.
            k0,k1,k2    Roujean’s model parameters for wavelength of
                        simulation SOS.Wa
            ind         Refractive index air / ground for wavelength SOS.Wa
            file        Pre-calculated user file (complete path)
                        or “DEFAULT” for code calculation
            """
        self.type = 5
        self.alb = alb
        self.wind = None
        self.ind = ind
        self.file = file

        self.roujeank0 = k0
        self.roujeank1 = k1
        self.roujeank2 = k2

        self.alpha = None
        self.beta = None

    def setLabertianNadal(self, k0, k1, k2, ind, alpha, beta,
                          alb=0, file="DEFAULT"):
        """ Lambertian surface with vegetation BRDF (Roujean’s model) and
            ground polarization BPDF (Nadal’s model)
            alb         Surface albedo for wavelength SOS.
            k0,k1,k2    Roujean’s model parameters for wavelength of
                        simulation SOS.Wa
            ind         Refractive index air / ground for wavelength SOS.Wa
            file        Pre-calculated user file (complete path)
                        or “DEFAULT” for code calculation
            alpha       alpha Nadal parameter
            beta        beta Nadal parameter
            """
        self.type = 5
        self.alb = alb
        self.wind = None
        self.ind = ind
        self.file = file

        self.roujeank0 = k0
        self.roujeank1 = k1
        self.roujeank2 = k2

        self.alpha = alpha
        self.beta = beta

class AP(object):
    """ Atmospheric profile parameters object."""

    def __init__(self, mot=None, hr=8.0, ha=2.0):
        """ Init function for the atmospheric profile
            mot         Molecular optical thickness for the wavelength of
                        radiance simulation
                        > 0.0001 : To be considered
                        0        : To ginore
            hr          Molecular heigth scale (km).
            ha          Aerosol heigth scale (km).
            """

        self.mot = mot
        self.hr = hr
        self.ha = ha
        self.type = 1
        self.zmin = None
        self.zmax = None
        self.usefile = None

    def setHeightScale(self, ha):
        """ Profile defined by heights scales
            ha          Aerosol heigth scale (km).
            """
        self.ha = ha
        self.type = 1
        self.zmin = None
        self.zmax = None


    def setMixAltitude(self, zmin, zmax):
        """ Profile for a mixture of molecules and aerosols
            between altitudes Zmin and Zmax (in km)"""
        self.ha = None
        self.type = 2
        self.zmin = zmin
        self.zmax = zmax

    def UserFile(self, usefile):
        """Pre-calculated atmospheric profile file (complete path).
           The user profile is supposed to be a real profile, without
           adjustment to an aerosol phase function truncature.
        """
        self.usefile = usefile
        self.mot = None
        self.hr = None
        self.ha = None
        self.type = None
        self.zmin = None
        self.zmax = None
        self.usefile = None

class SOS(object):
    """ This class creates the SOS object which configures and runs the
        simulation"""

    def __init__(self, wa=0.440, resroot=None, mdf=None,
                 view=1, phi=90, dphi=None, igmax=None,
                 lpolar=0, outputlevel=-1):
        """ This method initiates the OSOAA class
            wa          Wavelength of radiance calculation (microns).
            mdf         Molecular depolarization factor.
            resroot     Working folder for the OSOAA computations (complete
                        path).
            view        Option for field of view representation
                        1 : Viewing in a constant azimuth plan
                        2 : Polar diagram
            phi         Relative azimuth (degrees)
            dphi        Step on the azimuth (degrees, integer value)
            igmax       Maximal order of interaction (scattering & surface
                        reflexion).
            lpolar      Option for polarization turn-off
                        0 : simulation without polarization
            outputlevel Option for output level definitions
                        -1 : Default output : upward radiance at TOA,
                             downward radiance at ground level.
                         N : Specific output : upward and downward
                             radiance at level N of the atmospheric profile
        """

        self.wa = wa
        self.root = os.getenv("OSOAA_ROOT")
        self.view = view
        self.phi = phi
        self.dphi = None
        self.igmax = igmax
        self.lpolar = lpolar

        if resroot is None:
            rnd = ''.join(random.choice(string.ascii_uppercase
                                        + string.ascii_lowercase
                                        + string.digits) for _ in range(16))
            self.resroot = self.root+"/results/"+rnd
        else:
            self.resroot = resroot

    def setConstantView(self, phi):
        """ Set constant view mode for a relative azimuth angles
            phi         Relative azimuth (degrees)
            """
        self.view = 1
        self.phi = phi
        self.dphi = None

    def setPolarDiagram(self, dphi):
        """ Set polar plot for an angle step
            dphi        Step on the azimuth (degrees, integer value)
            """
        self.view = 2
        self.phi = None
