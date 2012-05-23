import os

from numpy import array, exp, pi, random, sqrt, sum, zeros

from pyculgi import *

from systems import *


__all__ = ["Simulation", "loadpath"]

# All parameters used for these simulations.
parameters = [  "phase", "size", "polymer", "beadvolume", "density", "nchi",
                "kappa", "temperature", "expansion",
                "ca", "cb", "a", "mobility", "population",
                "timestep", "totaltime" ,
                "chmob", "cc", "npname", "sigma" ]
# All instant result types -- no spaces tolerated!
instant_fields_all = "\
SCMBondEnergy,SCMBendingEnergy,SCMTorsionEnergy,SCMNBEnergy,\
SCMElectrostaticsEnergy,SCMExternalPotentialEnergy,SCMPotentialEnergy,\
GCDensityInhomogeneity,GCTotalFreeEnergy,GCIdealFreeEnergy,GCContactFreeEnergy,\
GCCompressibilityFreeEnergy,GCElectrostaticFreeEnergy,CouplingEnergy"

# Maximum number of iterations when solving nonlinear equations.
maxiters = 25000


def getpath(sim):
    """Construct output path for a given simulation object."""

    # First directory just differentiates between phases.
    phase = sim.phase
    dir1 = "phase%i" %phase

    # Second directory: system size, polymer architecture, bead size.
    X,Y,Z = sim.size
    if sim.phase > 14:
        dir2 = "%ix%ix%i_%s_bv%.2f_sigma%.2f" %(X, Y, Z, sim.polymer, sim.beadvolume, sim.sigma)
    else:
        dir2 = "%ix%ix%i_%s_bv%.2f" %(X, Y, Z, sim.polymer, sim.beadvolume)

    # Third directory: system temperature, noise, density, population.
    # From phase 10, add the chain mobility here as a parameter.
    dir3 = "temp%.2f_exp%.2f_den%.1f_pop%i" %(sim.temperature, sim.expansion, sim.density, sim.population)
    if sim.phase > 9:
        dir3 += "_chmob%.2f" %sim.chmob

    # Fourth directory: compressibility, interactions, mobilities.
    # Previously 2 significant digits were used for mobility, now
    #  three or even four or five are needed sometimes,
    #  but we still need to support all the possbilities to be backwards compatible.
    if sim.mobility >= 0.01:
        mobilityformat = "mob%.2f"
    elif sim.mobility >= 0.001:
        mobilityformat = "mob%.3f"
    else:
        mobilityformat = "mob%.5f"
    format = "k%.1f_nchi%.1f"
    args = (sim.kappa, sim.nchi)
    if sim.population > 0:
        if sim.phase >= 13:
            format += "_cc%.1f"
            args += (sim.cc,)
        format += "_ca%.1f_cb%.1f_"+mobilityformat
        args += (sim.ca, sim.cb, sim.mobility)
    if sim.population > 1:
        format += "_a%.1f"
        args += (sim.a,)
    dir4 = format %args

    return "%s/%s/%s/%s" %(dir1, dir2, dir3, dir4)


def getname(sim):
    """Construct output file basename from simulation object."""

    # Until phase 11, this was just the total time, now it contains the time step also.
    # From phase 14, also add the nanoparticle model name.
    name = "tt%i" %sim.totaltime
    if sim.phase > 10:
        name += "_ts%s" %str(sim.timestep)
    if sim.phase > 13:
        name += "_%s" %sim.npname

    return name


def loadpath(path, setup=True, main=False):
    """Create simulation object from an output file."""

    # Separate the path and file name.
    outpath,outfile = os.path.split(path.strip())

    # There will always be four levels  in the path.
    dir1,dir2,dir3,dir4 = outpath.split('/')

    # The first is just the simulation phase.
    phase = int(dir1[5:])

    # The second contain systems size, polymer architecture and bead volume.
    if phase > 14:
        size,pol,bv,sigma = dir2.split("_")
    else:
        size,pol,bv = dir2.split('_')
    size = map(int,size.split('x'))
    bv = float(bv[2:])
    if phase > 14:
        sigma = float(sigma[5:])

    # The third has temperature, exp. factor, density and population.
    # Starting phase 10, this also contains the chain mobility (chmob).
    if phase > 9:
        temp,exp,den,pop,chmob = dir3.split('_')
    else:
        temp,exp,den,pop = dir3.split('_')
    temp = float(temp[4:])
    exp = float(exp[3:])
    den = float(den[3:])
    pop = int(pop[3:])
    if phase > 9:
        chmob = float(chmob[5:])

    # The fourth has kappa, nchi and other things depending on the population.
    # Three cases -- no NPs, one, more than one -- have different path structures
    if pop == 0:
        k,nchi = dir4.split('_')
    elif pop == 1:
        if phase >= 13:
            k,nchi,cc,ca,cb,mob = dir4.split('_')
        else:
            k,nchi,ca,cb,mob = dir4.split('_')
    else:
        if phase >= 13:
            k,nchi,cc,ca,cb,mob,a = dir4.split('_')
        else:
            k,nchi,ca,cb,mob,a = dir4.split('_')
    k = float(k[1:])
    nchi = float(nchi[4:])
    if pop > 0:
        if phase >= 13:
            cc = float(cc[2:])
        ca = float(ca[2:])
        cb = float(cb[2:])
        mob = float(mob[3:])
    else:
        ca = 0.0
        cb = 0.0
        mob = 0.0
    if pop > 1:
        a = float(a[1:])
    else:
        a = 0.0

    # The output file name contain just the total time,
    #  and the time step after phase 10.
    # Before phase 10, the time step is always 0.01
    outname = outfile[:-4]
    if phase > 13:
        totaltime, timestep, npname = outname.split("_")
        totaltime = int(totaltime[2:])
        timestep = float(timestep[2:])
        npname = npname
    elif phase > 10:
        totaltime, timestep = outname.split("_")
        totaltime = int(totaltime[2:])
        timestep = float(timestep[2:])
    else:
        totaltime = int(outname[2:])
        timestep = 0.01

    # Create the simulation object
    # Since chmob appears from phase 10, the call has more arguments
    args = [phase]
    args += [size, pol, bv]
    args += [temp, exp, den, pop]
    args += [k, nchi, ca, cb, a, mob]
    args += [timestep, totaltime]
    if phase > 9:
        args += [chmob]
    if phase > 12:
        args += [cc]
    if phase > 13:
        args += [npname]
    if phase > 14:
        args += [sigma]
    sim = Simulation(*args)

    # Setup the simulation if requested
    if setup:
        sim.setup()

    # If main is set (refers to main gallery), modify the path for the gallery
    if main:
        sim.gallerypath = sim.path

    return sim


class Simulation:

    def __init__(self, phase,
                 size, polymer, beadvolume,
                 temperature, expansion, density, population,
                 kappa, nchi, ca, cb, a, mobility,
                 timestep, totaltime,
                 chmob=None, cc=0.0, npname=None, sigma=0.8):
        """Initialize the simulation.

        This, along with initialization, does the following tasak:
            - sets all relevant parameters passed to the constructor as attributes
            - calculates most importnat derived parameters
            - corrects the effective density for the presence of NPs
            - sets some paths and names

        Everything that takes more than a second or needs some memory should be in setup.
        """

        # Set all the parameters passed.
        for parameter in parameters:
            if parameter != None:
                setattr(self, parameter, eval(parameter))

        # Aliases.
        self.mob = self.mobility
        self.pop = self.population
        self.ts = self.timestep
        self.tt = self.totaltime

        # Number of time steps to run, as well as the frequency of
        #   snapshots and energy saves all follow from the time step and total time.
        # If there are not enough time steps, default frequencies to one.
        self.nsteps = int(self.tt / self.ts)
        self.efreq = int(1.0*self.tt / Nenerg[self.phase-1] / self.ts) or 1
        self.sfreq = int(1.0*self.tt / Nsnaps[self.phase-1] / self.ts) or 1

        # Coupling parameters for the core beads to simplify code later on.
        self.cca = self.ca
        self.ccb = self.cb
        if self.phase > 12:
            self.cca = self.cc
            self.ccb = self.cc

        # Other derived parameters.
        self.dcoupling = self.cb - self.ca
        self.volume = size[0]*size[1]
        self.area = self.volume

        # The nanoparticle prototype.
        # The prototype does not need much memory, so leave it here, but the
        #   box and block copolymer both need significant memory, so create them
        #   in the setup() method, which should be called when we are sure this
        #   simulation is the one we want to run.
        # Until phase 7, the nanoparticle is a single bead, soft core molecule, but
        #   from phase 8 it is a soft core colloid, loaded from a Culgi object file.
        # After phase 13, the model name is variable, determining which file to load
        #   from an assumed, fixed location (and that's the final approach).
        if self.phase <= 7:
            self.np = Palette.CreateSoftCoreMolecule("np", "P")
        else:
            if self.phase <= 13:
                self.npname = "np"
            else:
                self.npname = npname
            self.np = Palette.LoadSoftCoreColloid("phase%i/%s.cof" %(self.phase,self.npname))

        # Number of beads per nanoparticle.
        self.beads_per_np = self.np.GetNumberOfBeads()

        # Estimate nanoparticle volume given the NP model and simulation parameters.
        # Again, until phase 7 things are easy (NPs are single beads), but afterwards
        #   they are made up of a number of beads. At first all those beads have the
        #   the same coupling parameter, but after phase 12, the central bead
        #   can have a different coupling parameter than the shell beads, but that's
        #   still quite simple since there is always only one core bead (for now).
        # In all cases, assume the particle will be in the A domain in the end,
        #   so use the corresponding A-coupling parameter here.
        self.kd1 = self.kappa*self.density + 1
        self.npvolume = self.ca / self.kd1
        if self.phase > 7:
            self.npvolume = self.beads_per_np * self.ca / self.kd1
        if self.phase > 12:
            self.npvolume = self.cc/self.kd1 + (self.beads_per_np-1)*self.ca/self.kd1
        self.dexcluded = self.npvolume

        # Estimate the effective density of the polymer field, with correction as of phase 2.
        self.effective_density = self.density
        if self.phase > 1:
            self.effective_density = self.density - self.pop*self.npvolume/self.volume
        self.effden = self.effective_density

        # Output paths and names used for this run.
        self.path = getpath(self)
        self.name = getname(self)
        self.outname = "%s.out" %self.name
        self.outpath = "%s/%s" %(self.path, self.outname)
        self.gallerypath = self.path.replace("phase%s/" %phase,"")

    def setup(self):
        """Setup the simulation before running.

        This sets up the simulation box and calculator before actually running,
        and generally does things that are too long for the cosntructor.
        """

        # The diblock copolymer model, which we always want to have the same reference density.
        # Use the effective density, which for phase 1 default to the reference density.
        self.bcp = Palette.CreateGaussianChain("bcp", self.polymer, *self.size)
        self.bcp.CreateHomogeneousRelativeDensity(self.effective_density)

        # Pointer to the total density and block densities.
        self.bcp_total = self.bcp.GetRelativeDensityFieldCmds()
        self.bcp_A = self.bcp.GetBeadTypeRelativeDensityFieldCmds("A")
        self.bcp_A.SetDisplayClampLevels(0.0, self.density)
        self.bcp_B = self.bcp.GetBeadTypeRelativeDensityFieldCmds("B")
        self.bcp_B.SetDisplayClampLevels(0.0, self.density)

        # The absolute chi can now be derived based on the BCP model.
        self.chi = self.nchi / self.bcp.GetNumBeads()

        # The simulation box, with BCP and nanoparticles.
        # Note that the nanoparticles are colloids starting from phase 8.
        self.box = Palette.CreateMesoBox("box", *self.size)
        self.box.AddGaussianChain(self.bcp)
        if self.phase <= 7:
            self.box.CopySoftCoreMolecules(self.np, self.pop)
        else:
            self.box.CopySoftCoreColloids(self.np, self.pop)

        # Convenience hooks for nanoparticle commands.
        if self.phase <= 7:
            get_np_cmds = self.box.GetSoftCoreMoleculeCmds
        else:
            get_np_cmds = self.box.GetSoftCoreColloidCmds
        self.nanoparticles = [get_np_cmds("np", i) for i in range(self.pop)]
        self.NPs = self.nanoparticles

        # Initial distribution of nanoparticles inside the simulation box.
        # For phase 5 and 6, use some ordered starting distribution.
        # Starting phase 7, use a random starting distribution again (equilibrate later).
        # Warning! Due to a mistake, phase 7 actually uses the same initial configuration as phase 6.
        # Since this is 2d, the Z coordinates is always equal to 0.5, though.
        # In phase 12, use a sampled initial distribution, to get an experimental-like initial RDF.
        Z = 0.5
        if self.phase == 12 and self.pop != 0:
            dh = 1.6
            ds = 1.33*sqrt(0.9*self.area/self.pop/pi)
            coords = zeros((self.pop, 3))
            norms = lambda ni,r: sqrt(sum((coords[:ni,:2]-r)**2, axis=1))
            hardwalls = lambda ni,r: sum(norms(ni,r) <= dh)
            rdf = lambda n: 1/(1+exp(-(n-dh-ds))) + 1.0*exp(-(n-dh-ds)**2)
            coords[:,2] = Z
            coords[0,:2] = array(self.size[:2]) * random.random(2)
            self.nanoparticles[0].SetCenterOfMass(*coords[0])
            for ni in range(1,self.pop):
                r = coords[ni-1,:2]
                while hardwalls(ni,r) > 0 or random.random() <= exp(-1/rdf(min(norms(ni,r)))):
                    r = array(self.size[:2]) * random.random(2)
                coords[ni,:2] = r
                self.nanoparticles[ni].SetCenterOfMass(*coords[ni])
        elif self.phase in [5,6,7] and self.pop != 0:
            if self.phase == 5:
                rfactor = 1.75
            else:
                rfactor = 3.25
            Nx = int(sqrt(self.pop))
            Ny = int(self.pop // Nx)
            n = self.pop - Nx*Ny
            dX = self.size[0] * 1.0 / Nx
            dY = self.size[1] * 1.0 / (Ny+1)
            for i in range(self.pop):
                dx = dX*MathManager.Rand()/rfactor
                dy = dY*MathManager.Rand()/rfactor
                if i < Nx*Ny:
                    X = dX * (i%Nx + dx)
                    Y = dY * (i//Nx + dy)
                else:
                    X = (i - Nx*Ny + dx) * self.size[0] * 1.0 / n
                    Y = dY * (Ny + dy)
                self.nanoparticles[i].SetCenterOfMass(X, Y, Z)
        else:
            for i in range(self.pop):
                X = self.nanoparticles[i].GetCenterOfMass().GetX()
                Y = self.nanoparticles[i].GetCenterOfMass().GetY()
                self.nanoparticles[i].SetCenterOfMass(X, Y, Z)

        # From phase 16, also orient the nanoparticle randomly.
        if self.phase > 15 and self.pop != 0:
            for i in range(self.pop):
                self.nanoparticles[i].Rotate(0, 0, 1, 180*random.random())

        # The calculator object is a MBF (meso bead-field) hybrid.
        # We will always be using the external potential (entropy field) dynamic scheme.
        self.calc = CalculatorsManager.CreateMBFHybridCalculator()
        self.calc.SetSystem(self.box)
        self.calc.SetTimeStep(self.timestep)
        self.calc.SetFieldDiffusionMethod("EntropyFieldDynamics")

        # The output file settings.
        # Note that the archive plugin needs to be used for this starting Culgi 6.
        self.calc.SetSaveInstantResultsOn(instant_fields_all)
        self.calc.SetInstantResultsSaveFrequency(self.efreq)
        if Culgi.GetVersionID()[0] == "6":
                self.archive = PluginsManager.CreateArchivePlugin()
                self.calc.AddPlugin(self.archive)
                self.archive.SetFieldsArchiveOn()
                self.archive.SetParticlesArchiveOn()
                self.archive.SetFrequency(self.sfreq)
        if Culgi.GetVersionID()[0] == "5":
        	self.calc.SetSaveRelativeDensityFieldsOn()
        	self.calc.SetSaveCoordinatesOn()
        	self.calc.SetCoordinateSaveFrequency(self.sfreq)
        	self.calc.SetRelativeDensityFieldSaveFrequency(self.sfreq)

        # Calculator parameter commands and simulation engine.
        self.params = self.calc.GetParametersCmds()
        self.params_GC = self.params.GetGCParametersCmds()
        self.params_SC = self.params.GetSCMParametersCmds()
        self.engine = self.calc.GetDDFTEngineCmds()
        self.engine.SetMaxIterations(maxiters)

        # Parameters pertaining to mobility.
        # Until phase 7, nanoparticle are single beads (SC molecules with 'P' beads),
        #   so the mobility has the same meaning as in standard Culgi.
        # From phase 8, nanoparticles are colloids with two different kinds of beads,
        #  but the diffusion constant is set for the colloid as a whole.
        # Phase 8 assumes zero angular diffusion, but phase 9 and higher do not.
        # From phase 10, the GC chain mobilities can also change.
        self.calc.SetTemperature(self.temperature)
        self.params_GC.SetExpansionParameter("A", self.expansion)
        self.params_GC.SetExpansionParameter("B", self.expansion)
        if (self.phase <= 7):
            self.params_SC.SetDiffusionFactor("P", self.mobility)
        else:
            for i in range(self.pop):
                self.nanoparticles[i].SetConstVelocity('Z', 0.0)
                self.nanoparticles[i].SetDiffusionFactor(self.mobility)
                if self.phase == 8:
                    self.nanoparticles[i].SetConstAngularVelocity(0.0, 0.0, 0.0)
                else:
                    self.nanoparticles[i].SetConstAngularVelocity('X', 0.0)
                    self.nanoparticles[i].SetConstAngularVelocity('Y', 0.0)
                    self.nanoparticles[i].SetRotationalFactors(0.0, 0.0, self.mob)
        if self.phase > 9:
            self.params_GC.SetDiffusionFactor("A", self.chmob)
            self.params_GC.SetDiffusionFactor("B", self.chmob)

        # Parameters concerning interactions.
        # From phase 8, nanoparticle are colloids with two different types of beads.
        self.calc.SetKappa(self.kappa)
        self.params_GC.SetBeadVolume("A", self.beadvolume)
        self.params_GC.SetBeadVolume("B", self.beadvolume)
        self.params_GC.SetChi("A", "B", self.chi)
        if self.phase <= 7:
            self.params.SetBeadFieldCoupling("P", "A", self.ca)
            self.params.SetBeadFieldCoupling("P", "B", self.cb)
            self.params_SC.SetA("P", "P", self.a)
        else:
            self.params.SetBeadFieldCoupling("C", "A", self.cca)
            self.params.SetBeadFieldCoupling("C", "B", self.ccb)
            self.params.SetBeadFieldCoupling("S", "A", self.ca)
            self.params.SetBeadFieldCoupling("S", "B", self.cb)
            self.params_SC.SetA("C", "C", self.a)
            self.params_SC.SetA("C", "S", self.a)
            self.params_SC.SetA("S", "S", self.a)
        if self.phase > 14:
            self.params_SC.SetSigma("C", self.sigma)
            self.params_SC.SetSigma("S", self.sigma)

        # Output path and file name
        self.calc.SetOutputFileName(self.name)


    def getcenterofmass(self, name, num=0):
        """Return the center of mass of a soft core molecule."""

        com = self.box.GetSoftCoreMoleculeCmds(name, num).GetCenterOfMass()
        return com.GetX(), com.GetY(), com.GetZ()

    def getorderparameter(self, x, y):
        """Return the order parameter at coordinates x,y."""

        a = self.bcp_A.GetGaussianInterpolatedValue(x, y, 0.0)
        b = self.bcp_B.GetGaussianInterpolatedValue(x, y, 0.0)
        return a-b

    def gettotalfield(self, x, y):
        """Return total field density at coordinates x,y."""

        return self.bcp_total.GetGaussianInterpolatedValue(x, y, 0.0)

    def run(self):
        """Run the simulation.

        The whole simulation ultimately consists of three steps:
            - equilibration of the NPs without the field
            - equilibration of the field around nanoparticles
            - the simulation proper
        """

        # Equilibrations stat from phase 4
        if self.phase >= 4:

            # Turn off all archives
            self.calc.SetSaveInstantResultsOff()
            if Culgi.GetVersionID()[0] == "5":
                self.calc.SetSaveRelativeDensityFieldsOff()
                self.calc.SetSaveCoordinatesOff()
            if Culgi.GetVersionID()[0] == "6":
                self.calc.RemovePlugin(self.archive)

            # No phase separation at all during equilibration
            self.params_GC.SetChi("A", "B", 0.0)

            # From phase 7, equilibrate the nanoparticles first,
            #   without them fealing the field at all.
            if self.phase >= 7:
                if self.phase == 7:
                    self.params.SetBeadFieldCoupling("P", "A", 0.0)
                    self.params.SetBeadFieldCoupling("P", "B", 0.0)
                else:
                    self.params.SetBeadFieldCoupling("C", "A", 0.0)
                    self.params.SetBeadFieldCoupling("C", "B", 0.0)
                    self.params.SetBeadFieldCoupling("S", "A", 0.0)
                    self.params.SetBeadFieldCoupling("S", "B", 0.0)
                self.calc.SetFieldDiffusionOff()
                self.calc.Run(int(Teq_np[self.phase-1]/self.timestep))

            # Now equilibrate fields with only base coupling,
            #   keeping the nanoparticle positions fixed.
            if self.phase <= 7:
                self.params_SC.SetDiffusionFactor("P", 0)
                self.params.SetBeadFieldCoupling("P", "A", self.ca)
                self.params.SetBeadFieldCoupling("P", "B", self.ca)
            else:
                for i in range(self.pop):
                    self.nanoparticles[i].SetDiffusionFactor(0.0)
                self.params.SetBeadFieldCoupling("C", "A", self.cca)
                self.params.SetBeadFieldCoupling("C", "B", self.cca)
                self.params.SetBeadFieldCoupling("S", "A", self.ca)
                self.params.SetBeadFieldCoupling("S", "B", self.ca)
            self.calc.SetFieldDiffusionOn()
            self.calc.SetBeadDiffusionOff()
            self.calc.Run(int(Teq_field[self.phase-1]/self.timestep))

            # Set relevant parameters back to the original values
            self.params_GC.SetChi("A", "B", self.chi)
            if self.phase <= 7:
                self.params_SC.SetDiffusionFactor("P", self.mobility)
                self.params.SetBeadFieldCoupling("P", "A", self.ca)
                self.params.SetBeadFieldCoupling("P", "B", self.cb)
            else:
                for i in range(self.pop):
                    self.nanoparticles[i].SetDiffusionFactor(self.mobility)
                self.params.SetBeadFieldCoupling("C", "A", self.cca)
                self.params.SetBeadFieldCoupling("C", "B", self.ccb)               
                self.params.SetBeadFieldCoupling("S", "A", self.ca)
                self.params.SetBeadFieldCoupling("S", "B", self.cb)               
            self.calc.SetBeadDiffusionOn()

            # Turn archiving back on
            if Culgi.GetVersionID()[0] == "5":
                self.calc.SetSaveRelativeDensityFieldsOn()
                self.calc.SetSaveCoordinatesOn()
            if Culgi.GetVersionID()[0] == "6":
                self.calc.AddPlugin(self.archive)
            self.calc.SetSaveInstantResultsOn(instant_fields_all)

        # Run the simulation proper.
        self.calc.Run(self.nsteps)

