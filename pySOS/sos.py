# coding=utf-8

import os
import random
import string
from io import open
from .outputs import OUTPUTS


class LOG(object):
    """ Log files for the simulation.
    """
    ang = "Angle.log"
    profile = "Profile.log"
    aer = "Aerosols.log"
    aermie = "AerMie.log"
    surface = "Surface.log"
    sos = "SOS.Log"
    config = "SOS_config.txt"


class RESULTS(object):
    """ Result files for the simulation.
        """
    profile = "Profile.txt"
    aer = "Aerosols.txt"
    angrad = "SOS_angrad.txt"
    angaer = "SOS_angaer.txt"
    bin = "SOS_Result.bin"
    trans = "SOS_transm.txt"
    advup = "SOS_Up.txt"
    advdown = "SOS_Down.txt"
    advupuser = "SOS_Up_user.txt"
    advdownuser = "SOS_Down_user.txt"


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

    def __init__(self, thetas=32.48, radnb=24, raduser=None,
                 aernb=40, aeruser=None):
        """ Init the angle class
            thetas      Solar zenith angle (degrees)
            radnb       Number of Gauss angles to be used for radiance
                        computations
            raduser     Filename of the complementary list of user's angles to
                        complete the ANG.Rad.NbGauss angles (complete path).
            aernb       Number of Gauss angles to be used for mie computations
            aeruser     Filename of the complementary list of user's angles to
                        complete the ANG.mie.NbGauss angles (complete path).
            """

        self.rad = self.ANGLES(radnb, raduser)
        self.aer = self.ANGLES(aernb, aeruser)
        self.thetas = thetas


class SURFACE(object):
    """ Surface definitions class."""

    def __init__(self, type=1, alb=0.02, ind=1.34, wind=2.0):
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
        self.wind = wind
        self.ind = ind
        self.file = "DEFAULT"

        self.roujeank0 = None
        self.roujeank1 = None
        self.roujeank2 = None

        self.alpha = None
        self.beta = None

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
        self.type = 1
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
        self.type = 2
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
        self.type = 6
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

    def __init__(self, waref=0.550, aotref=0.3, tronca=1, model=1):
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
        self.wmo = AEROSOLMODELS.WMO(wmotype=2)
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

    def __init__(self, mot=0.230, hr=8.0, ha=2.0):
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

        self.ha = 2.0
        self.zmin = None
        self.zmax = None

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

    def __init__(self, wa=0.440, resroot=None, mdf=0.0279,
                 view=1, phi=0, dphi=None, igmax=30,
                 ipolar=0, outputlevel=-1):
        """ This method initiates the SOS class
            wa          Wavelength of radiance calculation (microns).
            mdf         Molecular depolarization factor.
            resroot     Working folder for the SOS computations (complete
                        path).
            view        Option for field of view representation
                        1 : Viewing in a constant azimuth plan
                        2 : Polar diagram
            phi         Relative azimuth (degrees)
            dphi        Step on the azimuth (degrees, integer value)
            igmax       Maximal order of interaction (scattering & surface
                        reflexion).
            ipolar      Option for polarization turn-off
                        0 : simulation without polarization
            outputlevel Option for output level definitions
                        -1 : Default output : upward radiance at TOA,
                             downward radiance at ground level.
                         N : Specific output : upward and downward
                             radiance at level N of the atmospheric profile
        """

        self.wa = wa
        self.root = os.getenv("RACINE")
        self.view = view
        self.phi = phi
        self.dphi = None
        self.igmax = igmax
        self.ipolar = ipolar
        self.outputlevel = outputlevel
        self.mdf = mdf

        if resroot is None:
            rnd = ''.join(random.choice(string.ascii_uppercase
                                        + string.ascii_lowercase
                                        + string.digits) for _ in range(8))
            self.resroot = self.root+"/results/"+rnd
        else:
            self.resroot = resroot

        self.ang = ANG()
        self.surface = SURFACE()
        self.ap = AP()
        self.aer = AER()
        self.results = RESULTS()
        self.log = LOG()

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

    def run(self, root=None):
        """ Run SOS. If no root directory is given for SOS the one
            configured by the system is used.
            """

        if root is not None:
            self.root = root

        ## Previous config
        sc = "export SOS_RACINE=$RACINE"
        # Results
        sc = sc+"\n"+"export SOS_RESULT={}".format(self.resroot)
        # WMO and S&F files
        sc = sc+"\n"+"export SOS_RACINE_FIC=$RACINE/fic"
        # Surface storage
        sc = sc+"\n"+"export dirSUNGLINT=$SOS_RESULT/SURFACE/GLITTER"
        sc = sc+"\n"+"export dirROUJEAN=$SOS_RESULT/SURFACE/ROUJEAN"
        sc = sc+"\n"+"export dirRH=$SOS_RESULT/SURFACE/RH"
        sc = sc+"\n"+"export dirBREON=$SOS_RESULT/SURFACE/BREON"
        sc = sc+"\n"+"export dirNADAL=$SOS_RESULT/SURFACE/NADAL"
        sc = sc+"\n"+"export dirMIE=$SOS_RESULT/MIE"
        sc = sc+"\n"+"export dirLOG=$SOS_RESULT/LOG"
        sc = sc+"\n"+"export dirRESULTS=$SOS_RESULT/SOS"

        sc = sc+"\n"+"{}/exe/main_SOS.ksh \\".format(self.root)
        #   Definition of the working folder : #CHECK!!!
        #   ----------------------------------
        sc = sc+"\n"+"-SOS.ResRoot {} \\".format(self.resroot)
        #
        #   Angles calculation parameters :
        #   --------------------------------
        sc = sc+"\n"+"-ANG.Thetas {} \\".format(self.ang.thetas)
        if self.ang.rad.nbgauss is not None:
            sc = sc+"\n"+"-ANG.Rad.NbGauss {} \\".format(self.ang.rad.nbgauss)
        if self.ang.rad.userangfile is not None:
            sc = sc+"\n"+"-ANG.Rad.UserAngFile {} \\".format(
                self.ang.rad.userangfile)
        if self.results.angrad is not None:
            sc = sc+"\n"+"-ANG.Rad.ResFile ${{dirRESULTS}}/{} \\".format(self.results.angrad)
        if self.ang.aer.nbgauss is not None:
            sc = sc+"\n"+"-ANG.Aer.NbGauss {} \\".format(self.ang.aer.nbgauss)
        if self.ang.aer.userangfile is not None:
            sc = sc+"\n"+"-ANG.Aer.UserAngFile {} \\".format(
                self.ang.aer.userangfile)
        if self.results.angaer is not None:
            sc = sc+"\n"+"-ANG.Aer.ResFile ${{dirRESULTS}}/{} \\".format(self.results.angaer)
        if self.log.ang is not None:
            sc = sc+"\n"+"-ANG.Log ${{dirLOG}}/{} \\".format(self.log.ang)
        #
        #   Radiance calculation parameters :
        #   --------------------------------
        if self.log.sos is not None:
            sc = sc+"\n"+"-SOS.Log ${{dirLOG}}/{} \\".format(self.log.sos)
        sc = sc+"\n"+"-SOS.Wa  {} \\".format(self.wa)
        #
        sc = sc+"\n"+"-SOS.View {} \\".format(self.view)
        if self.view is 1:
            sc = sc+"\n"+"-SOS.View.Phi {} \\".format(self.phi)
        elif self.view is 2:
            sc = sc+"\n"+"-SOS.View.Phi {} \\".format(self.dphi)
        sc = sc+"\n"+"-SOS.View.Level {} \\".format(self.outputlevel)
        #
        if self.log.sos is not None:
            sc = sc+"\n"+"-SOS.Log ${{dirLOG}}/{} \\".format(self.log.sos)
        if self.igmax is not None:
            sc = sc+"\n"+"-SOS.IGmax {} \\".format(self.igmax)

        sc = sc+"\n"+"-SOS.Ipolar {} \\".format(self.ipolar)
        if self.mdf is not None:
                sc = sc+"\n"+"-SOS.MDF {} \\".format(self.mdf)
        if self.results.bin is not None:
            sc = sc+"\n"+"-SOS.ResBin ${{dirRESULTS}}/{} \\".format(self.results.bin)

        if self.results.advup is not None:
            sc = sc+"\n"+"-SOS.ResFileUp ${{dirRESULTS}}/{} \\".format(self.results.advup)
        if self.results.advdown is not None:
            sc = sc+"\n"+"-SOS.ResFileDown ${{dirRESULTS}}/{} \\".format(self.results.advdown)
        if self.results.advupuser is not None:
            sc = sc+"\n"+"-SOS.ResFileUp.UserAng ${{dirRESULTS}}/{} \\".format(self.results.advupuser)
        if self.results.advdown is not None:
            sc = sc+"\n"+"-SOS.ResFileDown.UserAng ${{dirRESULTS}}/{} \\".format(self.results.advdownuser)
        if self.results.trans is not None:
            sc = sc+"\n"+"-SOS.Trans ${{dirRESULTS}}/{} \\".format(self.results.trans)
        if self.log.config is not None:
            sc = sc+"\n"+"-SOS.Config ${{dirRESULTS}}/{} \\".format(self.log.config)
        #
        #   Profile parameters :
        #   -------------------
        if self.ap.usefile is None:
            if self.results.profile is not None:
                sc = sc+"\n"+"-AP.ResFile ${{dirRESULTS}}/{} \\".format(self.results.profile)
            #     Atmospheric Profile parameters
            sc = sc+"\n"+"-AP.MOT {} \\".format(self.ap.mot)
            sc = sc+"\n"+"-AP.HR {} \\".format(self.ap.hr)
            sc = sc+"\n"+"-AP.Type {} \\".format(self.ap.type)
            if self.ap.type is 1:
                sc = sc+"\n"+"-AP.AerHS.HA {} \\".format(self.ap.ha)
            if self.ap.type is 2:
                sc = sc+"\n"+"-AP.AerLayer.Zmax {} \\".format(self.ap.zmax)
                sc = sc+"\n"+"-AP.AerLayer.Zmin {} \\".format(self.ap.zmin)
        elif self.ap.usefile is not None:
            sc = sc+"\n"+"-AP.UseFile {} \\".format(self.ap.usefile)
        if self.log.profile is not None:
            sc = sc+"\n"+"-AP.Log ${{dirLOG}}/{} \\".format(self.log.profile)
        #
        #   Aerosols parameters :
        #   ---------------------
        if self.log.aer is not None:
            sc = sc+"\n"+"-AER.Log ${{dirLOG}}/{} \\".format(self.log.aer)
        if self.aer.usefile is None:
            if self.results.aer is not None:
                sc = sc+"\n"+"-AER.ResFile ${{dirRESULTS}}/{} \\".format(self.results.aer)
            if self.log.aermie is not None:
                sc = sc+"\n"+"-AER.MieLog ${{dirLOG}}/{} \\".format(self.log.aermie)
            sc = sc+"\n"+"-AER.Waref  {} \\".format(self.aer.waref)
            sc = sc+"\n"+"-AER.AOTref {} \\".format(self.aer.aotref)
            sc = sc+"\n"+"-AER.Tronca {} \\".format(self.aer.tronca)
            if self.aer.aotref > 0.0:
                sc = sc+"\n"+"-AER.Model {} \\".format(self.aer.model)
            #     Aerosols parameters for mono-modal models :
            if self.aer.model is 0:
                sc = sc+"\n"+"-AER.MMD.MRwa {} \\".format(self.aer.mm.mrwa)
                sc = sc+"\n"+"-AER.MMD.MIwa {} \\".format(self.aer.mm.miwa)
                if self.wa is not self.aer.waref:
                    sc = sc+"\n"+"-AER.MMD.MRwaref {} \\".format(self.aer.mm.mrwaref)
                    sc = sc+"\n"+"-AER.MMD.MIwaref {} \\".format(self.aer.mm.miwaref)
                sc = sc+"\n"+"-AER.MMD.SDtype {} \\".format(self.aer.mm.sdtype)
                if self.aer.mm.sdtype is 1:
                    sc = sc+"\n"+"-AER.MMD.LNDradius {} \\".format(self.aer.mm.sdradius)
                    sc = sc+"\n"+"-AER.MMD.LNDvar {} \\".format(self.aer.mm.sdvar)
                elif self.aer.mm.sdtype is 2:
                    sc = sc+"\n"+"-AER.MMD.JD.slope {} \\".format(self.aer.mm.slope)
                    sc = sc+"\n"+"-AER.MMD.JD.rmin {} \\".format(self.aer.mm.rmin)
                    sc = sc+"\n"+"-AER.MMD.JD.rmax {} \\".format(self.aer.mm.rmax)
            #     Aerosols parameters for WMO models :
            elif self.aer.model is 1:
                sc = sc+"\n"+"-AER.WMO.Model {} \\".format(self.aer.wmo.model)
                if self.aer.wmo.model is 4:
                    sc = sc+"\n"+"-AER.WMO.DL {} \\".format(self.aer.wmo.dl)
                    sc = sc+"\n"+"-AER.WMO.WS {} \\".format(self.aer.wmo.ws)
                    sc = sc+"\n"+"-AER.WMO.OC {} \\".format(self.aer.wmo.oc)
                    sc = sc+"\n"+"-AER.WMO.SO {} \\".format(self.aer.wmo.so)
            #     Aerosols parameters for Shettle&Fenn models :
            elif self.aer.model is 2:
                sc = sc+"\n"+"-AER.SF.Model {} \\".format(self.aer.sf.model)
                sc = sc+"\n"+"-AER.SF.RH {} \\".format(self.aer.sf.rh)
            #     Aerosols parameters for LND bi-modal models :
            elif self.aer.model is 3:
                sc = sc+"\n"+"-AER.BMD.VCdef {} \\".format(self.aer.lnb.vcdef)
                if self.aer.lnb.vcdef is 1:
                    sc = sc+"\n"+"-AER.BMD.CoarseVC {} \\".format(self.aer.lnb.coarsevc)
                    sc = sc+"\n"+"-AER.BMD.FineVC {} \\".format(self.aer.lnb.finevc)
                elif self.aer.lnb.vcdef is 2:
                    sc = sc+"\n"+"-AER.BMD.RAOT {} \\".format(self.aer.lnb.raot)
                sc = sc+"\n"+"-AER.BMD.CM.MRwa {} \\".format(self.aer.lnb.cmrwa)
                sc = sc+"\n"+"-AER.BMD.CM.MIwa {} \\".format(self.aer.lnb.cmiwa)
                sc = sc+"\n"+"-AER.BMD.CM.MRwaref {} \\".format(self.aer.lnb.cmrwaref)
                sc = sc+"\n"+"-AER.BMD.CM.MIwaref {} \\".format(self.aer.lnb.cmiwaref)
                sc = sc+"\n"+"-AER.BMD.CM.SDradius {} \\".format(self.aer.lnb.csdradius)
                sc = sc+"\n"+"-AER.BMD.CM.SDvar {} \\".format(self.aer.lnb.csdvar)
                sc = sc+"\n"+"-AER.BMD.FM.MRwa {} \\".format(self.aer.lnb.fmrwa)
                sc = sc+"\n"+"-AER.BMD.FM.MIwa {} \\".format(self.aer.lnb.fmiwa)
                sc = sc+"\n"+"-AER.BMD.FM.MRwaref {} \\".format(self.aer.lnb.fmrwaref)
                sc = sc+"\n"+"-AER.BMD.FM.MIwaref {} \\".format(self.aer.lnb.fmiwaref)
                sc = sc+"\n"+"-AER.BMD.FM.SDradius {} \\".format(self.aer.lnb.fsdradius)
                sc = sc+"\n"+"-AER.BMD.FM.SDvar {} \\".format(self.aer.lnb.fsdvar)
            #    Aerosols parameters for external data (phase functions, scattering
            #    and extinction coefficients) :
            elif self.aer.model is 4:
                sc = sc+"\n"+"-AER.ExtData {} \\".format(self.aer.extdata)
        elif self.aer.usefile is not None:
            sc = sc+"\n"+"-AER.UseFile {} \\".format(self.aer.usefile)
        #
        #   Surface :
        #   --------------------------------------

        sc = sc+"\n"+"-SURF.Log ${{dirLOG}}/{} \\".format(self.log.surface)
        sc = sc+"\n"+"-SURF.File {} \\".format(self.surface.file)
        sc = sc+"\n"+"-SURF.Type {} \\".format(self.surface.type)
        sc = sc+"\n"+"-SURF.Alb {} \\".format(self.surface.alb)

        if self.surface.ind is not None:
            sc = sc+"\n"+"-SURF.Ind {} \\".format(self.surface.ind)

        if self.surface.wind is not None:
            sc = sc+"\n"+"-SURF.Glitter.Wind {} \\".format(self.surface.wind)

        if self.surface.roujeank0 is not None:
            sc = sc+"\n"+"-SURF.Roujean.K0 {} \\".format(self.surface.roujeank0)

        if self.surface.roujeank1 is not None:
            sc = sc+"\n"+"-SURF.Roujean.K1 {} \\".format(self.surface.roujeank1)

        if self.surface.roujeank2 is not None:
            sc = sc+"\n"+"-SURF.Roujean.K2 {} \\".format(self.surface.roujeank2)

        if self.surface.alpha is not None:
            sc = sc+"\n"+"-SURF.Nadal.Alpha {} \\".format(self.surface.alpha)

        if self.surface.beta is not None:
            sc = sc+"\n"+"-SURF.Nadal.Beta {} \\".format(self.surface.beta)

        if not os.path.exists(self.resroot):
            os.makedirs(self.resroot)

        # We generate the script
        with open(self.resroot+"/script.kzh", 'w') as file:
            file.write(sc[:-2])

        # Run script with ksh
        os.system("ksh "+self.resroot+"/script.kzh")

        # read OUTPUTS
        self.outputs = OUTPUTS(self.resroot, self.results)

def test():
    s = SOS()
    s.run()
    print("SOS wrapper script by Francisco Nemiña")
    print("Inspired by Py6S wrapper by Robin Wilson")
    print("Using SOS located at {}".format(s.root))
    print("Running SOS using a set of test parameters")
    print("The results are:")
    print("Expected result: 0.346131")
    print("Actual result: {}".format(s.outputs.up.I[0]))
    if (s.outputs.up.I[0]  == 0.346131):
        print("#### Results agree PyOSOAA is working correctly")
    if (s.outputs.up.I[0] != 0.346131):
        print("#### Results do not agree PyOSOAA is not working correctly")
