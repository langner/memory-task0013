"""Contains main functions and class for task0013."""


import os

from numpy import array, average, exp, pi, random
from numpy import round, sqrt, std, sum, zeros

from pyculgi import *

from systems import *


__all__ = ["Simulation", "loadpath"]


# All parameters used for these simulations.
parameters = [  "phase", "size", "polymer", "beadvolume", "density", "nchi",
                "kappa", "temperature", "expansion",
                "ca", "cb", "a", "mobility", "population",
                "timestep", "totaltime" ,
                "chmob", "cc", "npname", "sigma",
                "disp", "poly", "stencilsize" ]

# All instant result types -- no spaces tolerated!
instant_fields_all = "\
SCMBondEnergy,SCMBendingEnergy,SCMTorsionEnergy,SCMNBEnergy,\
SCMElectrostaticsEnergy,SCMExternalPotentialEnergy,SCMPotentialEnergy,\
GCDensityInhomogeneity,GCTotalFreeEnergy,GCIdealFreeEnergy,GCContactFreeEnergy,\
GCCompressibilityFreeEnergy,GCElectrostaticFreeEnergy,CouplingEnergy"

# Maximum number of iterations when solving nonlinear equations.
# This used to be 25k, but that turned out to be to small occasionally.
maxiters = 50000


def getpath(sim):
    """Construct output path for a given simulation object."""

    # First directory just differentiates between phases.
    dir1 = "phase%i" %sim.phase

    # Second directory: system size, polymer architecture, bead size.
    # After phase 14, the sigma parameter is appended to this directory name.
    # After phase 19, the coupling stencil can also be changed.
    dir2 = "%ix%ix%i_%s_bv%.2f" %(sim.size[0], sim.size[1], sim.size[2], sim.polymer, sim.beadvolume)
    if sim.phase > 14:
        dir2 += "_sigma%.2f" %sim.sigma
    if sim.phase > 19:
        dir2 += "_ss%i" %sim.stencilsize

    # Third directory: system temperature, noise, density, population.
    # From phase 10, add the chain mobility here as a parameter.
    # From phase 17 we want to consider very low temperatures, so adjust the format accordingly.
    if (sim.temperature == 0.00) or (sim.temperature >= 0.01):
        tempformat = "temp%.2f"
    else:
        tempformat = "temp%.4f"
    dir3 = (tempformat + "_exp%.2f_den%.1f_pop%i") %(sim.temperature, sim.expansion, sim.density, sim.population)
    if sim.phase > 9:
        dir3 += "_chmob%.2f" %sim.chmob

    # Fourth directory: compressibility, interactions, mobilities.
    # The NP mobility and coupling parameters are added only if the population is larger than zero,
    #   and interparticle interaction parameters are added when the population is larger than one.
    # At first, two significant digits were used for mobility, now three to five digits
    #   are needed, but we still need to support just two for backwards compatibility.
    # After phase 13, this directory also hold the coupling parameter for core beads, which
    #   were differentiated from the coupling parameters of shell beads from phase 14.
    dir4 = "k%.1f_nchi%.1f" %(sim.kappa, sim.nchi)
    if sim.population > 0:
        if sim.phase >= 13:
            dir4 += "_cc%.1f" %sim.cc
        if sim.mobility >= 0.01:
            mobilityformat = "mob%.2f"
        elif sim.mobility >= 0.001:
            mobilityformat = "mob%.3f"
        else:
            mobilityformat = "mob%.5f"
        dir4 += ("_ca%.1f_cb%.1f_"+mobilityformat) %(sim.ca, sim.cb, sim.mobility)
    if sim.population > 1:
        dir4 += "_a%.1f" %sim.a

    return "%s/%s/%s/%s" %(dir1, dir2, dir3, dir4)


def getname(sim):
    """Construct output file basename from simulation object."""

    # Up to phase 10, this was just the total time, afterwards it also contains the time step.
    # After phase 13, also add the nanoparticle model name.
    # After phase 17, also add the initial dispersion factor.
    # After phase 18, also add the dispersion factor.
    name = "tt%i" %sim.totaltime
    if sim.phase > 10:
        name += "_ts%s" %str(sim.timestep)
    if sim.phase > 13 and sim.pop > 0:
        name += "_%s" %sim.npname
    if sim.phase > 17 and sim.pop > 1:
        name += "_disp%.1f" %sim.disp
    if sim.phase > 18 and sim.pop > 1:
        name += "_poly%.2f" %sim.poly

    return name


def loadpath(path, setup=True, main=False):
    """Create simulation object from an output file."""

    # Separate the path and file name.
    outpath,outfile = os.path.split(path.strip())

    # There will always be four levels  in the path.
    dir1,dir2,dir3,dir4 = outpath.split('/')

    # The first is just the simulation phase.
    phase = int(dir1[5:])

    # The second contains the systems size, polymer architecture and bead volume.
    # After phase 14, the sigma parameter is in this directory name, too.
    # After phase 19, the stencil size is also in this directory name.
    if phase > 19:
        size,pol,bv,sigma,stencilsize = dir2.split("_")
    elif phase > 14:
        size,pol,bv,sigma = dir2.split("_")
    else:
        size,pol,bv = dir2.split('_')
    size = map(int,size.split('x'))
    bv = float(bv[2:])
    if phase > 14:
        sigma = float(sigma[5:])
    if phase > 19:
        stencilsize = int(stencilsize[2:])

    # The third has temperature, exp. factor, density and population.
    # After phase 9, this also contains the chain mobility (chmob).
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

    # The fourth has kappa, nchi and other things, depending on the population.
    # The three cases -- no NPs, one, more than one -- have different path structures.
    # And after phase 12, this directory name can also contain the core bead coupling parameter.
    if pop == 0:
        k,nchi = dir4.split('_')
    elif pop == 1:
        if phase > 12:
            k,nchi,cc,ca,cb,mob = dir4.split('_')
        else:
            k,nchi,ca,cb,mob = dir4.split('_')
    else:
        if phase > 12:
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

    # The output file name contains just the total time, and the time step
    #   after phase 10 (before that the time step is always 0.01).
    # After phase 13, the name also contains the NP model name (if there are any NPs).
    # After phase 17, the name can also contain the dispersion factor (if there are any NPs).
    # After phase 18, the name can also contain the spread for polydispersity.
    outname = outfile[:-4]
    if phase < 11:
        totaltime = int(outname[2:])
        timestep = 0.01
    else:
        if pop == 0 or phase < 14:
            totaltime, timestep = outname.split("_")
            npname = None
            disp = 0.0
        elif phase < 18:
            totaltime, timestep, npname = outname.split("_")
        elif phase < 19:
            totaltime, timestep, npname, disp = outname.split("_")
            disp = float(disp[4:])
        else:
            totaltime, timestep, npname, disp, poly = outname.split("_")
            disp = float(disp[4:])
            poly = float(poly[4:])
        totaltime = int(totaltime[2:])
        timestep = float(timestep[2:])

    # Create the simulation object now.
    # After phase 9, the parameter list starts to grow.
    args = [phase]
    args += [temp, mob, a, ca, cb]
    args += [pol, bv, den, exp, k, nchi]
    args += [size, pop]
    args += [timestep, totaltime]
    if phase > 9:
        args += [chmob]
    if phase > 12:
        args += [cc]
    if phase > 13:
        args += [npname]
    if phase > 14:
        args += [sigma]
    if phase > 17:
        args += [disp]
    if phase > 18:
        args += [poly]
    if phase > 19:
        args += [stencilsize]
    sim = Simulation(*args)

    # Setup the simulation if requested.
    if setup:
        sim.setup()

    # If main is set (refers to main gallery), modify the path to the gallery.
    if main:
        sim.gallerypath = sim.path

    return sim


def restrict_colloid_z(colloid):
    """Restrict all movements of a colloids to the XY plane."""

    colloid.SetConstVelocity("Z", 0.0)
    colloid.SetConstAngularVelocity("X", 0.0)
    colloid.SetConstAngularVelocity("Y", 0.0)


def gen_shell_beads():
    """Generator for all types of polydisperse shell beads.

    Will return the digits (int), sign (int),
    value of the delta (float) and name of bead (string).
    """

    for digits in range(1,100):
        for sign in -1,1:
            delta = sign * digits / 100.0
            name = "S%s%s" %("p"*(sign>0) or "m", str(digits).zfill(2))
            yield digits, sign, delta, name


def set_shell_beads_A(sim, a=None):
    """Set A parameter between all pairs of shell bead types."""

    # For any given shell bead type, there will be interactions between
    #   core beads (C), normal shell beads (S), the same bead type,
    #   as well as all other shell bead types.
    # This below is a bit redundant, but it's still fast.
    # If no value provided, use that of simulation object.
    a = a or sim.a
    for digits1,sign1,delta1,name1 in gen_shell_beads():
        sim.params_SC.SetA("C", name1, a)
        sim.params_SC.SetA("S", name1, a)
        sim.params_SC.SetA(name1, name1, a)
        for digits2,sign2,delta2,name2 in gen_shell_beads():
            sim.params_SC.SetA(name1, name2, a)
            sim.params_SC.SetA(name2, name2, a)


class Simulation:

    def __init__(self, phase,
                temperature, mobility, a, ca, cb,
                polymer, beadvolume, density, expansion, kappa, nchi,
                size, population,
                timestep, totaltime,
                chmob=None, cc=0.0, npname=None, sigma=0.8,
                disp=4.0, poly=0.0, stencilsize=3):
        """Initialize the simulation.

        This, along with initialization, does the following tasks:
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

        # Aliases
        self.arch = self.polymer
        self.bv = self.beadvolume
        self.mob = self.mobility
        self.pop = self.population
        self.ss = self.stencilsize
        self.ts = self.timestep
        self.tt = self.totaltime

        # Number of time steps to run, as well as the frequency of
        #   snapshots and energy saves all follow from the time step and total time.
        # If there are not enough time steps, default frequencies are one (every time step).
        self.nsteps = int(self.tt / self.ts)
        self.efreq = int(1.0*self.tt / phases_energs[self.phase-1] / self.ts) or 1
        self.sfreq = int(1.0*self.tt / phases_snaps[self.phase-1] / self.ts) or 1

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
        if self.phase < 8:
            self.np = Palette.CreateSoftCoreMolecule("np", "P")
            self.npname = "np"
        else:
            if self.phase < 14:
                self.npname = "np"
            else:
                self.npname = npname
            self.np = Palette.LoadSoftCoreColloid("phase%i/%s.cof" %(self.phase,self.npname))

        # If there are no nanoparticles, set the NP model name to None (orders the gallery nicely).
        if self.pop == 0:
            self.npname = None

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

        # Estimate the effective density of the polymer field, with a correction as of phase 2.
        self.effective_density = self.density
        if self.phase > 1:
            self.effective_density = self.density - self.pop*self.npvolume/self.volume
        self.effden = self.effective_density

        # Output paths and names used for this run.
        self.path = getpath(self)
        self.name = getname(self)
        self.outroot = "%s/%s" %(self.path, self.name)
        self.outname = "%s.out" %self.name
        self.outpath = "%s/%s" %(self.path, self.outname)
        self.gallerypath = self.path.replace("phase%s/" %phase,"")

    def setup(self):
        """Setup the simulation before running.

        This sets up the simulation box and calculator before actually running,
        and generally does things that are too long for the cosntructor.
        """

        # The diblock copolymer model, which we always want to have the same reference density.
        # Use the effective density, which for phase 1 defaults to the reference density.
        self.bcp = Palette.CreateGaussianChain("bcp", self.polymer, *self.size)
        self.bcp.CreateHomogeneousRelativeDensity(self.effective_density)

        # Pointers to the total density and block densities.
        self.bcp_total = self.bcp.GetRelativeDensityFieldCmds()
        self.bcp_A = self.bcp.GetBeadTypeRelativeDensityFieldCmds("A")
        self.bcp_A.SetDisplayClampLevels(0.0, self.density)
        self.bcp_B = self.bcp.GetBeadTypeRelativeDensityFieldCmds("B")
        self.bcp_B.SetDisplayClampLevels(0.0, self.density)

        # The absolute chi can now be derived based on the BCP model.
        self.chi = self.nchi / self.bcp.GetNumBeads()

        # Create the simulation box, with the BCP and nanoparticles.
        # Note that the nanoparticles are colloids starting from phase 8.
        self.box = Palette.CreateMesoBox("box", *self.size)
        self.box.AddGaussianChain(self.bcp)
        if self.phase < 8:
            self.box.CopySoftCoreMolecules(self.np, self.pop)
        else:
            self.box.CopySoftCoreColloids(self.np, self.pop)

        # Convenience hooks for nanoparticle commands.
        if self.phase < 8:
            get_np_cmds = self.box.GetSoftCoreMoleculeCmds
        else:
            get_np_cmds = self.box.GetSoftCoreColloidCmds
        self.nanoparticles = [get_np_cmds("np", i) for i in range(self.pop)]
        self.NPs = self.nanoparticles

        # Initial distribution of nanoparticles inside the simulation box.
        # Since this is 2D, the Z coordinats is always equal to 0.5.
        # For phases 5 and 6, use some ordered starting distribution.
        # After phase 6, use a random starting distribution again (equilibrate later). Due to a mistake,
        #   phase 7 actually initially used the same initial configuration as phase 6 (now fixed). So in this
        #   case, all we need to do is set the Z coordinates, because Culgi randomizes by default.
        # In phase 12, use a kind-of sampled initial distribution, to get an experimental-like initial RDF.
        # After phase 16, the distribution is generated by equilibrating in a smaller box,
        #   but that is done later in during the actual simulation run (so this here is irrelevant).
        Z = 0.5
        if self.pop != 0:
            if self.phase == 12:
                dh = 1.6
                ds = 1.33 * sqrt(0.9*self.area/self.pop/pi)
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
            elif self.phase in [5,6,7]:
                rfactor = 1.75*(self.phase == 5) or 3.25
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
            elif self.phase < 17:
                for i in range(self.pop):
                    X = self.nanoparticles[i].GetCenterOfMass().GetX()
                    Y = self.nanoparticles[i].GetCenterOfMass().GetY()
                    self.nanoparticles[i].SetCenterOfMass(X, Y, Z)

        # After phase 15, also orient the nanoparticle randomly.
        if self.phase > 15 and self.pop != 0:
            for i in range(self.pop):
                self.nanoparticles[i].Rotate(0, 0, 1, 180*random.random())

        # The calculator object is a MBF (meso bead-field) hybrid.
        # We will always be using the external potential (entropy field) dynamic scheme.
        # After phase 19, we might be changing the stencil size.
        self.calc = CalculatorsManager.CreateMBFHybridCalculator()
        self.calc.SetFieldDiffusionMethod("EntropyFieldDynamics")
        self.calc.SetTimeStep(self.timestep)
        self.calc.SetSystem(self.box)
        if self.phase > 19:
            self.calc.SetStencilSize(self.stencilsize)

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

        # Calculator parameter commands and simulation engine commands.
        self.params = self.calc.GetParametersCmds()
        self.params_GC = self.params.GetGCParametersCmds()
        self.params_SC = self.params.GetSCMParametersCmds()
        self.engine = self.calc.GetDDFTEngineCmds()
        self.engine.SetMaxIterations(maxiters)

        # Parameters pertaining to mobility.
        # Up to phase 7, nanoparticle are single beads (SC molecules with 'P' beads),
        #   so the mobility has the same meaning as in standard Culgi Brownian motion.
        # After phase 7, nanoparticles are colloids with two different kinds of beads,
        #   but the diffusion constant is set for the colloid as a whole (custom code).
        # Phase 8 assumes zero angular diffusion, but phase 9 and higher phases do not,
        #   and the rotational factor about the Z axis is the same as the translational
        #   mobility. Not that restrict_colloid_z takes care of fixing the Z axis.
        # After phase 9, the GC chain mobilities can also change.
        self.calc.SetTemperature(self.temperature)
        self.params_GC.SetExpansionParameter("A", self.expansion)
        self.params_GC.SetExpansionParameter("B", self.expansion)
        if (self.phase < 8):
            self.params_SC.SetDiffusionFactor("P", self.mobility)
        else:
            for i in range(self.pop):
                restrict_colloid_z(self.nanoparticles[i])
                self.nanoparticles[i].SetDiffusionFactor(self.mobility)
                if self.phase == 8:
                    self.nanoparticles[i].SetConstAngularVelocity(0.0, 0.0, 0.0)
                else:
                    self.nanoparticles[i].SetRotationalFactors(0.0, 0.0, self.mob)
        if self.phase > 9:
            self.params_GC.SetDiffusionFactor("A", self.chmob)
            self.params_GC.SetDiffusionFactor("B", self.chmob)

        # Parameters related directly to interactions.
        # After phase 7, nanoparticle are colloids with two different types of beads.
        # After phase 16, coupling and repulsion will be set later, after equilibration.
        self.calc.SetKappa(self.kappa)
        self.params_GC.SetBeadVolume("A", self.beadvolume)
        self.params_GC.SetBeadVolume("B", self.beadvolume)
        self.params_GC.SetChi("A", "B", self.chi)
        if self.phase < 8:
            self.params.SetBeadFieldCoupling("P", "A", self.ca)
            self.params.SetBeadFieldCoupling("P", "B", self.cb)
            self.params_SC.SetA("P", "P", self.a)
        elif self.phase < 17:
            self.params.SetBeadFieldCoupling("C", "A", self.cca)
            self.params.SetBeadFieldCoupling("C", "B", self.ccb)
            self.params.SetBeadFieldCoupling("S", "A", self.ca)
            self.params.SetBeadFieldCoupling("S", "B", self.cb)
            self.params_SC.SetA("C", "C", self.a)
            self.params_SC.SetA("C", "S", self.a)
            self.params_SC.SetA("S", "S", self.a)

        # After phase 14, the sigma parameter can also change.
        if self.phase > 14:
            self.params_SC.SetSigma("C", self.sigma)
            self.params_SC.SetSigma("S", self.sigma)

        # After phase 18, we also have some notion of polydispersity.
        # This is implemented by modifying the coupling parameter of shell beads,
        #   or more precisely by changing their bead types to new bead types with
        #   slightly different coupling parameters. These new beat types
        #   have names with the pattern (S[mp][0-9][0-9]), where [mp] means
        #   plus or minus and the digits are the digits found in the delta,
        #   that is the difference between the normal couplind and that new
        #   shell bead type in parts per hundred. That makes 2*99 new bead types
        # We want a normal distribution of coupling parameters, so the coupling
        #   is first taken from a distribution and then rounded and converted
        #   to a bead type with an appropriate name. Bead types are defined
        #   implicitely, so it will be enough to just set the coupling parameter
        #   and then the beat type for each shell bead. If zero is sampled,
        #   the bead type of course does not change.
        # Note that Culgi requires parentheses around the bead type name when
        #   it is longer than one character.
        # Overall, the average coupling does not change, and so the effective
        #   density in the field does not need to be corrected. What should be
        #   seen is a population of NPs with varying shapes.
        # Remember that the interaction will need to be set again after equilibration.
        if self.phase > 18 and self.poly > 0.0:

            # Set the coupling parameters and sigmas for all new shell bead types.
            for digits,sign,delta,name in gen_shell_beads():
                self.params.SetBeadFieldCoupling(name, "A", self.ca + delta)
                self.params.SetBeadFieldCoupling(name, "B", self.cb + delta)
                self.params_SC.SetSigma(name, self.sigma)

            # We also need to set repulsions for each pair of new beads.
            # This will will be redone for equilibration, so there is a separate procedure.
            set_shell_beads_A(self)

            # Generate all deltas in one go, and check their mean and standard
            #   deviation (and bail out if not appropriate). The core bead is not
            #   to be changed, so there are population*(beads_per_np-1) deltas.
            deltas = random.normal(0.0, self.poly, self.pop*(self.beads_per_np-1))
            deltas_avg = average(deltas)
            deltas_std = std(deltas)
            limit = 0.05
            assert abs(deltas_avg)/self.ca < limit, \
                "Average delta (%.3f) is more than %.2f from zero." %(deltas_avg, limit)
            assert (abs(deltas_std) - self.poly)/self.poly < limit, \
                "Deviation of deltas (%.3f) differs from choice (%.3f) by more than %.2f." %(deltas_std, self.poly, limit)

            # Round the deltas to two significant digits and make sure no values
            #   are smaller than -0.99 and none larger than 0.99.
            deltas = round(deltas, 2)
            deltas[deltas < -0.99] = -0.99
            deltas[deltas > 0.99] = 0.99

            # Now set new bead types of all shell beads.
            # Remember that the core bead does not change!
            for npi in range(self.pop):
                for bi in range(self.beads_per_np-1):
                    d = deltas[npi*(self.beads_per_np-1) + bi]
                    if d != 0.00:
                        name = "S%s%s" %("p"*(d>0) or "m", str(int(100*abs(d))).zfill(2))
                        self.nanoparticles[npi].GetBeadCmds(1+bi).SetElementType(name)

        # Output path and file name.
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

        # Equilibrations start after phase 3, and consist of two phases,
        #   one where the particles are equilibrated without any interactions
        #   with the field, and another where the density field is equilibrated with
        #   static nanoparticles and base interactions (the lower ones).
        # After phase 16, we do something less straightforward, namely we equilibrate
        #   the nanoparticles in a box where the DPD cutoff is comparable to the average
        #   distance between nanoparticles. The temperature is arbitrarily chosen as 0.25,
        #   which gave reasonably ordered, fairly noisy distribution for test runs.
        if self.phase > 3:

            # Turn off all archives.
            self.calc.SetSaveInstantResultsOff()
            if Culgi.GetVersionID()[0] == "5":
                self.calc.SetSaveRelativeDensityFieldsOff()
                self.calc.SetSaveCoordinatesOff()
            if Culgi.GetVersionID()[0] == "6":
                self.calc.RemovePlugin(self.archive)

            # No phase separation at all during equilibration.
            self.params_GC.SetChi("A", "B", 0.0)

            # After phase 6, first equilibrate nanoparticles not interacting with the field.
            if self.phase > 6:

                # Set all coupling parameters to zero.
                # For phase 7, nanoparticles are still molecules, but after that
                #   they are colloids and so the bead names differ.
                # After phase 18, we have many types of polydipserse beads.
                if self.phase == 7:
                    self.params.SetBeadFieldCoupling("P", "A", 0.0)
                    self.params.SetBeadFieldCoupling("P", "B", 0.0)
                else:
                    self.params.SetBeadFieldCoupling("C", "A", 0.0)
                    self.params.SetBeadFieldCoupling("C", "B", 0.0)
                    self.params.SetBeadFieldCoupling("S", "A", 0.0)
                    self.params.SetBeadFieldCoupling("S", "B", 0.0)
                if self.phase > 18:
                    for digits,sign,delta,name in gen_shell_beads():
                        self.params.SetBeadFieldCoupling(name, "A", 0.0)
                        self.params.SetBeadFieldCoupling(name, "B", 0.0)

                # No field diffusion in this part of the simulation.
                self.calc.SetFieldDiffusionOff()

                # After phase 16, we will create a temporary, smaller box to euqilibrate
                #   the particles at a higher density, so as to induce some degree
                #   of order. The size of the temporary box can be estimated by
                #   assuming an interparticle distance of 1.0 (counting from the shell)
                #   and calculating the appropriate area needed to achieve that.
                if self.phase > 16:

                    # Packing fraction for maximal circles in two dimensions.
                    packing = 0.91

                    # Distance between touching maximal circles in this box.
                    dmax = 2*sqrt(packing*self.area/(self.pop*pi))

                    # Approximate observed distance between nanoparticles before simulation.
                    # This is entirely phenomenological, although based on experimental data.
                    # The calibration is done so as not to go much over 3.0-4.0, which is just above
                    #   the practical limit of depletion interactions in these simulations, and to
                    #   approach the optimum 1.5-1.8 rather slowly (for higher concentrations).
                    # Since this callibration depends on other parameters, we allow it to be
                    #   adjusted by an additional dispersion parameter that describes the initial
                    #   extent of NP dispersion (lower disp parameter implies more dispersion).
                    dobs = 4.0*(1 - exp(-dmax/self.disp))

                    # Distance between core and shell beads in the nanoparticle model.
                    b0 = self.np.GetBeadCmds(0).GetCoordinates()
                    b1 = self.np.GetBeadCmds(1).GetCoordinates()
                    dx = b0.GetX() - b1.GetX()
                    dy = b0.GetY() - b1.GetY()
                    rnp = sqrt(dx**2 + dy**2)

                    # Now we want a temporary box that will scale dobs to 2*rnp+1,
                    #   because that is always the approximate final distance when
                    #   the DPD cutoff is 1. Add one to dimension for some excess space.
                    f = (2*rnp+1) / dobs
                    X = int(f*self.size[0])+1
                    Y = int(f*self.size[1])+1
                    box = Palette.CreateMesoBox("tmp", X, Y, 1)

                    # Add some nanoparticle to this box, too.
                    # After phase 18, we want to change the bead types of shell beads,
                    #   to correspond to the bead types in the original box.
                    box.CopySoftCoreColloids(self.np, self.pop)
                    if self.phase > 18:
                        for npi,np in enumerate(self.nanoparticles):
                            for bi in range(self.beads_per_np-1):
                                bt = np.GetBeadCmds(1+bi).GetElementType()
                                box.GetSoftCoreColloidCmds(npi).GetBeadCmds(1+bi).SetElementType(bt)

                    # Make this box the system of the simulation calculator, temporarily.
                    self.calc.SetSystem(box)

                    # These settings were tested for simpler DPD/Brownian systems.
                    self.calc.SetTemperature(0.1)
                    self.calc.SetTimeStep(0.01)
                    self.params_SC.SetA("C", "C", 10)
                    self.params_SC.SetA("C", "S", 10)
                    self.params_SC.SetA("S", "S", 10)
                    if self.phase > 18:
                        set_shell_beads_A(self, 10)
                    for i in range(self.pop):
                        np = box.GetSoftCoreColloidCmds(i)
                        np.SetCenterOfMass(X*random.random(), Y*random.random(), 0.5)
                        np.SetDiffusionFactor(0.01)
                        np.SetRotationalFactors(0.0, 0.0, 0.01)
                        restrict_colloid_z(np)

                # Run the equilibration.
                self.calc.Run(int(phases_eqtime_np[self.phase-1]/self.timestep))

                # We need to revert some of the simualtion parameters and transfer coordinates
                #   after equilibrating in phases after 16. The temporary box should be deleted.
                # Note that the NP coordinates need to be scaled back to the simulation box size.
                if self.phase > 16:

                    self.calc.SetSystem(self.box)
                    self.calc.SetTemperature(self.temperature)
                    self.calc.SetTimeStep(self.timestep)
                    self.params_SC.SetA("C", "C", self.a)
                    self.params_SC.SetA("C", "S", self.a)
                    self.params_SC.SetA("S", "S", self.a)
                    if self.phase > 18:
                        set_shell_beads_A(self)
                    for i in range(self.pop):
                        com = box.GetSoftCoreColloidCmds(i).GetCenterOfMass()
                        x = com.GetX() * self.size[0] / X
                        y = com.GetY() * self.size[1] / Y
                        self.nanoparticles[i].SetCenterOfMass(x, y, 0.5)

                    Palette.DeleteObject(box)

            # Now equilibrate fields with only base coupling, keeping nanoparticles fixed.
            # Remember that after phase 7 nanoparticles are colloids.
            if self.phase < 8:
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
                if self.phase > 18:
                    for digits,ssign,delta,name in gen_shell_beads():
                        self.params.SetBeadFieldCoupling(name, "A", self.ca)
                        self.params.SetBeadFieldCoupling(name, "B", self.ca)
            self.calc.SetFieldDiffusionOn()
            self.calc.SetBeadDiffusionOff()

            self.calc.Run(int(phases_eqtime_field[self.phase-1]/self.timestep))

            # Set relevant parameters back to their target simulation values.
            self.params_GC.SetChi("A", "B", self.chi)
            if self.phase < 8:
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
                if self.phase > 18:
                    for digits,sign,delta,name in gen_shell_beads():
                        self.params.SetBeadFieldCoupling(name, "A", self.ca+delta)
                        self.params.SetBeadFieldCoupling(name, "B", self.cb+delta)
            self.calc.SetBeadDiffusionOn()

            # Turn archiving back on.
            if Culgi.GetVersionID()[0] == "5":
                self.calc.SetSaveRelativeDensityFieldsOn()
                self.calc.SetSaveCoordinatesOn()
            if Culgi.GetVersionID()[0] == "6":
                self.calc.AddPlugin(self.archive)
            self.calc.SetSaveInstantResultsOn(instant_fields_all)

        # Run the target simulation.
        self.calc.Run(self.nsteps)

        # Save the box after running for reference.
        self.box.SaveAs(self.name)
