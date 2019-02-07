# coding=utf-8

import os
import random
import string
from io import open


class LOG(object):
    """ Log files for the simulation.
    """
    ang = "log_ang.txt"
    profile = "log_profile.txt"
    aer = "log_aer.txt"
    aermie = "log_aermie.txt"
    surface = "log_surface.txt"
    sos = "log_sos.txt"
    config = "config_sos.txt"


class RESULTS(object):
    """ Result files for the simulation.
        """
    profileatm = None
    aer = None
    angrad = None
    angaer = None
    bin = None
    trans = None
    advup = "resfile_advup.txt"
    advdown = "resfile_advdown.txt"
    advupuser = "resfile_advup_user.txt"
    advdownuser = "resfile_advdown_user.txt"


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


class AEROSOLMODELS(object):
    """ Aerosol models for the AER class."""
    class SF(object):
        """ Shettle and Fenn atmosphere model class."""
        # Shettle and Fenn models
        Tropospheric = 1
        Urban = 2
        Maritime = 3
        Coastal = 4

        def __init__(self, sfmodel, rh):
            """ Init method for the Shettle-Fenn model.
                model       Type of Shettle & Fenn model.
                                1 : Tropospheric S&F model.
                                2 : Urban S&F model.
                                3 : Maritime S&F model.
                                4 : Coastal S&F model.
                rh          Relative humidity (%) for Shettle & Fenn model.
                """
            self.model = sfmodel
            self.rh = rh

    class MM(object):
        """ Mono-modal size distribution"""

        def __init__(self, sdtype):
            """ Init method for the mono-modal size distribution
                sdtype      Type of mono-modal size distribution
                                1 : Log Normal size distribution
                                2 : Junge's law
                lnd         Log normal size distribution
                jd          Junge's law size distribution
                mrwa        Real part of the aerosol refractive index for the
                            wavelength of radiation calculation
                miwa        Imaginary part of the aerosol refractive index for
                            the  wavelength of radiation calculation
                sdradius    Modal radius (um) of the Log-Noprmal size
                            distribution
                sdvar       Standard deviation of the Log-Normal size
                            distribution
                slope       Slope of the Junge's law.
                            Warning: 3 is a singular value.
                rmin        Minimal radius of Junge's law (um)
                rmax        Maximal radius of Junge's law (um)
                mrwaref     Real part of the aerosol refractive index for the
                            reference wavelength of aerosol properties
                            calculation.
                miwaref     Imaginary part of the aerosol refractive index for
                            the reference wavelength of aerosol properties
                            calculation.
                """

            self.sdtype = sdtype
            self.mrwa = None
            self.miwa = None
            self.mrwaref = None
            self.miwaref = None
            if sdtype is 1:
                self.sdradius = None
                self.sdvar = None
            elif sdtype is 2:
                self.slope = None
                self.rmin = None
                self.rmax = None

    class WMO(object):
        """ WMO aerosol models."""

        def __init__(self, wmotype, dl=None, ws=None, oc=None, so=None):
            """ Init method for the WMO aerosol model
                wmotype     Type of WMO model
                                1 : Continental WMO model
                                2 : Maritime WMO model
                                3 : Urban WMO model
                                4 : WMO model by usef definition
                dl          Volume concentration (between 0 and 1) for
                            "Dust like" components
                ws          Volume concentration (between 0 and 1) for
                            "Water soluble" components
                oc          Volume concentration (between 0 and 1) for
                            "Oceanic" components
                so          Volume concentration (between 0 and 1) for
                            "Soot" components
                """

            self.model = wmotype

            if wmotype is 4:
                self.dl = dl
                self.ws = ws
                self.oc = oc
                self.so = so

    class LNB(object):
        """ Log-normal bi-modal aerosol model"""

        def __init__(self, vcdef):
            """ Log-normal bi-modal aerosol model init functions
                vcdef       Choide of the mixture description type
                                1 : Use of predefined volume concentrations.
                                2 : Use of the ratio of aerosol optical
                                    thickness (coarse mode ATO / total AOT)
                coarsevc    User volume concentration of the LND coarse mode
                finevc      User volume concentration of the LND fine mode
                raot        User value of the ration AOT_coarse/AOT_total for
                            the aerosol reference wavelength

                cmrwa       Real part of the aerosol refractive index for the
                            wavelength of radiation calculation for the coarse
                            mode
                cmiwa       Imaginary part of the aerosol refractive index for
                            the  wavelength of radiation calculationfor the
                            coarse mode
                csdradius   Modal radius (um) of the Log-Noprmal size
                            distribution for the coarse mode
                csdvar      Standard deviation of the Log-Normal size
                            distribution for the coarse mode
                cmrwaref    Real part of the aerosol refractive index for the
                            reference wavelength of aerosol properties
                            calculation for the coarse mode
                cmiwaref    Imaginary part of the aerosol refractive index for
                            the reference wavelength of aerosol properties
                            calculation for the coarse mode

                fmrwa       Real part of the aerosol refractive index for the
                            wavelength of radiation calculation for the fine
                            mode
                fmiwa       Imaginary part of the aerosol refractive index for
                            the  wavelength of radiation calculationfor the
                            fine mode
                fsdradius   Modal radius (um) of the Log-Noprmal size
                            distribution for the fine mode
                fsdvar      Standard deviation of the Log-Normal size
                            distribution for the fine mode
                fmrwaref    Real part of the aerosol refractive index for the
                            reference wavelength of aerosol properties
                            calculation for the fine mode
                fmiwaref    Imaginary part of the aerosol refractive index for
                            the reference wavelength of aerosol properties
                            calculation for the fine mode
                """

            self.vcdef = vcdef
            if vcdef is 1:
                self.coarsevc = None
                self.finevc = None
            elif vcdef is 2:
                self.raot = None

            self.cmrwa = None
            self.cmiwa = None
            self.cmrwaref = None
            self.cmiwaref = None
            self.csdradius = None
            self.csdvar = None

            self.fmrwa = None
            self.fmiwa = None
            self.fmrwaref = None
            self.fmiwaref = None
            self.fsdradius = None
            self.fsdvar = None


class AER(object):
    """ This class contains everything related to the aerosol components
        of the atmosphere."""

    def __init__(self, waref=0.550, aotref=0.1, tronca=None, model=2):
        """ Init method for the aerosol componentes class
            waref       Wavelength (microns) for reference aerosol optical
                        thickness.
            aotref      Aerosol optical thickness for the reference wavelength.
                        --> real value, without applied truncation.
            tronca      Option for no aerosol phase function troncation
                        (0 to not apply a troncation). Default is 1.
            model       Type of aerosol model
                            0 : Mono-modal
                            1 : WMO multi-modal
                            2 : Shettle & Fenn bi-modal
                            3 : Log-Normal bi-modal
                            4 : Phase function from an external source
            """

        self.waref = waref
        self.aotref = aotref
        self.tronca = tronca
        self.model = model
        self.sf = AEROSOLMODELS.SF(sfmodel=3, rh=98)
        self.usefile = None

    def SetModel(self, model=2,
                 sdtype=1,
                 sfmodel=3, rh=98,
                 wmotype=1, dl=None, ws=None, oc=None, so=None,
                 vcdef=2,
                 extdata=""):
        """ This methods sets the model for the AER class.

            model       Type of aerosol model
                            0 : Mono-modal
                            1 : WMO multi-modal
                            2 : Shettle & Fenn bi-modal
                            3 : Log-Normal bi-modal
                            4 : Phase function from an external source

            Mono-modal distribution parameters
            ----------------------------------
            mm          Mono-modal model atribute
            sdtype      Type of mono-modal size distribution
                            1 : Log Normal size distribution
                            2 : Junge's law

            WMO model parameters
            -------------------
            wmo         WMO model atribute
            wmotype     Type of WMO model
                            1 : Continental WMO model
                            2 : Maritime WMO model
                            3 : Urban WMO model
                            4 : WMO model by usef definition
            dl          Volume concentration (between 0 and 1) for
                        "Dust like" components
            ws          Volume concentration (between 0 and 1) for
                        "Water soluble" components
            oc          Volume concentration (between 0 and 1) for
                        "Oceanic" components
            so          Volume concentration (between 0 and 1) for
                        "Soot" components

            Shettle and Fenn model parameters
            ---------------------------------
            sf          Shettle and Fenn model atribute
            sfmodel       Type of Shettle & Fenn model.
                            1 : Tropospheric S&F model.
                            2 : Urban S&F model.
                            3 : Maritime S&F model.
                            4 : Coastal S&F model.
            rh          Relative humidity (%) for Shettle & Fenn model.

            Log-Normal bi-modal model parameters
            ------------------------------------
            lnd         Log-Normal bi-modal model atribute
            vcdef       Choide of the mixture description type
                            1 : Use of predefined volume concentrations.
                            2 : Use of the ratio of aerosol optical
                                thickness (coarse mode ATO / total AOT)

            External phase function
            -----------------------
            extdata     Filename (complete path) of user's external phase
                        function data and radiative parameters (extinction and
                        scattering coefficients).
                        The reference aerosol wavelength and the radiance
                        simulation wavelength must be equal
        """
        self.model = model
        if model is 0:
            self.mm = AEROSOLMODELS.MM(sdtype)
            self.wmo = None
            self.sf = None
            self.lnd = None
            self.external = None
        elif model is 1:
            self.mm = None
            self.wmo = AEROSOLMODELS.WMO(wmotype, dl, ws, oc, so)
            self.sf = None
            self.lnd = None
            self.external = None
        elif model is 2:
            self.mm = None
            self.wmo = None
            self.sf = AEROSOLMODELS.SF(sfmodel, rh)
            self.lnd = None
            self.external = None
        elif model is 3:
            self.mm = None
            self.wmo = None
            self.sf = None
            self.lnb = AEROSOLMODELS.LNB(vcdef)
            self.external = None
        elif model is 4:
            self.mm = None
            self.sf = None
            self.wmo = None
            self.lnd = None
            self.extdata = extdata

    def UserFile(self, usefile):
        """Pre-calculated aerosol profile file (complete path).
           The user file may contain aerosol radiative parameters for
            a phase function truncature (user has to check the correct
            agreement with the value of AER.Tronca parameter).
        """
        self.usefile = usefile


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
