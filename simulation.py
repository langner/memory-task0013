import os

from numpy import sqrt

from pyculgi import *


Nsnaps = 11000
Nenerg = 110000
maxiters = 25000
parameters = [  "phase", "size", "polymer", "beadvolume", "density", "nchi",
                "kappa", "temperature", "expansion",
                "ca", "cb", "a", "mobility", "population",
                "timestep", "totaltime" ]
Teq_np = 1000
Teq_field = 100

def printnow(fname, content, mode='a'):
    """ Print to a file without delay, so make sure to flush and close. """

    f = open(fname, mode)
    print >>f, content
    f.flush()
    f.close()

def getpath(sim):
    """ Construct output path based on simulation object. """

    # First directory just differentiates between phases.
    phase = sim.phase
    dir1 = "phase%i" %phase

    # Second directory: system size, polymer architecture, bead size.
    X,Y,Z = sim.size
    dir2 = "%ix%ix%i_%s_bv%.2f" %(X, Y, Z, sim.polymer, sim.beadvolume)

    # Third directory: system temperature, noise, density, population.
    dir3 = "temp%.2f_exp%.2f_den%.1f_pop%i" %(sim.temperature, sim.expansion, sim.density, sim.population)

    # Fourth directory: compressibility, interactions, mobilities.
    # Previously 2 significant digits were used for mobility, now
    #  three or even four of five are needed sometimes,
    #  but we still need to support all the possbilities.
    if sim.mobility >= 0.01:
        mobilityformat = "mob%.2f"
    elif sim.mobility >= 0.001:
        mobilityformat = "mob%.3f"
    else:
        mobilityformat = "mob%.5f"
    if sim.population == 0:
        format = "k%.1f_nchi%.1f"
        dir4 = format %(sim.kappa, sim.nchi)
    elif sim.population == 1:
        format = "k%.1f_nchi%.1f_ca%.1f_cb%.1f_"+mobilityformat
        dir4 = format %(sim.kappa, sim.nchi, sim.ca, sim.cb, sim.mobility)
    else:
        format = "k%.1f_nchi%.1f_ca%.1f_cb%.1f_"+mobilityformat+"_a%.1f"
        dir4 = format %(sim.kappa, sim.nchi, sim.ca, sim.cb, sim.mobility, sim.a)

    return "%s/%s/%s/%s" %(dir1, dir2, dir3, dir4)

def getname(sim):
    """ Construct output file basename from simulation object. """

    return "tt%i" %sim.totaltime  

def loadpath(path, setup=True, main=False):
    """ Create simulation object from an output file. """

    outpath,outfile = os.path.split(path.strip())
    dir1,dir2,dir3,dir4 = outpath.split('/')

    phase = int(dir1[5:])
    size,pol,bv = dir2.split('_')
    size = map(int,size.split('x'))
    bv = float(bv[2:])
    temp,exp,den,pop = dir3.split('_')
    temp = float(temp[4:])
    exp = float(exp[3:])
    den = float(den[3:])
    pop = int(pop[3:])

    if pop == 0:
        k,nchi = dir4.split('_')
        k = float(k[1:])
        nchi = float(nchi[4:])
        ca = 0.0
        cb = 0.0
        mob = 0.0
        a = 0.0
    elif pop ==1:
        k,nchi,ca,cb,mob = dir4.split('_')
        k = float(k[1:])
        nchi = float(nchi[4:])
        ca = float(ca[2:])
        cb = float(cb[2:])
        mob = float(mob[3:])
        a = 0.0
    else:
        k,nchi,ca,cb,mob,a = dir4.split('_')
        k = float(k[1:])
        nchi = float(nchi[4:])
        ca = float(ca[2:])
        cb = float(cb[2:])
        mob = float(mob[3:])
        a = float(a[1:])

    outname = outfile[:-4]
    totaltime = int(outname[2:])

    sim = simulation(phase,size,pol,bv,temp,exp,den,pop,k,nchi,ca,cb,a,mob,0.01,totaltime)
    sim.outpath = outpath
    if not main:
        sim.gallerypath = sim.outpath.replace("phase%s/" %phase,"")
    else:
        sim.gallerypath = sim.outpath
    sim.outfile = outfile
    sim.outname = outname
    sim.phase = phase
    if setup:
        sim.setup()

    return sim


class simulation:

    def __init__(self, phase, size, polymer, beadvolume, temperature, expansion, density, population, kappa, nchi, ca, cb, a, mobility, timestep, totaltime):

        # Set all the parameters
        for parameter in parameters:
            setattr(self, parameter, eval(parameter))

        # Some derived parameters
        self.dcoupling = self.cb - self.ca
        self.volume = size[0]*size[1]
        self.npvolume = self.ca/(self.kappa*self.density+1)
        if self.phase > 1:
            self.effective_density = self.density - self.population*self.npvolume/self.volume
        else:
            self.effective_density = self.density

        # The nanoparticle prototype and initial NP distribution
        # Until phase 7, the nanoparticle is a single bead soft core molecule
        # From phase 8, it is a soft core colloid loaded from a Culgi object file
        if self.phase <= 7:
            self.np = Palette.CreateSoftCoreMolecule("np", "P")
        else:
            self.np = Palette.LoadSoftCoreColloid("phase%i/np.cof" %self.phase)

        # We must correct the effective density from phase 8, because nanoparticles
        # are colloids with more than one bead.
        if self.phase > 7:
            self.npvolume = self.np.GetNumberOfBeads()*self.ca/(self.kappa*self.density+1)
            self.effective_density = self.density - self.population*self.npvolume/self.volume

        # Use this for organizing results (not necessary for actual simultion)
        # The difference cb-ca is not representative, use the excluded volume
        self.dexcluded = self.npvolume

        # The diblock copolymer
        self.bcp = Palette.CreateGaussianChain("bcp", self.polymer, *self.size)
        self.bcp.CreateHomogeneousRelativeDensity(self.effective_density)
        self.bcp_total = self.bcp.GetRelativeDensityFieldCmds()
        self.bcp_A = self.bcp.GetBeadTypeRelativeDensityFieldCmds("A")
        self.bcp_A.SetDisplayClampLevels(0.0, self.density)
        self.bcp_B = self.bcp.GetBeadTypeRelativeDensityFieldCmds("B")
        self.bcp_B.SetDisplayClampLevels(0.0, self.density)

        # The simulation box
        # Note again that the nanoparticles are colloids starting from phase 8
        self.box = Palette.CreateMesoBox("box", *self.size)
        self.box.AddGaussianChain(self.bcp)
        if self.phase <= 7:
            self.box.CopySoftCoreMolecules(self.np, self.population)
        else:
            self.box.CopySoftCoreColloids(self.np, self.population)

        # The nanoparticle hooks
        if self.phase <= 7:
            self.nanoparticles = [self.box.GetSoftCoreMoleculeCmds("np", i) for i in range(self.population)]
        else:
            self.nanoparticles = [self.box.GetSoftCoreColloidCmds("np", i) for i in range(self.population)]

    def setup(self):

        # These follow from the time step
        self.nsteps = int(self.totaltime/self.timestep)
        self.efreq = int(1.0*self.totaltime/Nenerg/self.timestep)
        self.sfreq = int(1.0*self.totaltime/Nsnaps/self.timestep)

        # Initial distribution of nanoparticles inside the simulation box
        # For phase 5 and 6, use some ordered starting distribution
        # Starting phase 7, use a random starting distribution again (equilibrate later)
        # (due to a mistake, phase 7 actually uses the same initial configuration as phase 6)
        if self.phase >= 5 and self.phase <=7 and self.population != 0:
            if self.phase == 5:
                rfactor = 1.75
            if self.phase >= 6:
                rfactor = 3.25
            Nx = int(sqrt(self.population))
            Ny = int(self.population // Nx)
            n = self.population - Nx*Ny
            dX = self.size[0] * 1.0 / Nx
            dY = self.size[1] * 1.0 / (Ny+1)
            Z = 0.5
            for i in range(self.population):
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
            for i in range(self.population):
                Z = 0.5
                X = self.nanoparticles[i].GetCenterOfMass().GetX()
                Y = self.nanoparticles[i].GetCenterOfMass().GetY()
                self.nanoparticles[i].SetCenterOfMass(X, Y, Z)

        # The calculator
        self.calc = CalculatorsManager.CreateMBFHybridCalculator()
        self.calc.SetSystem(self.box)
        self.calc.SetTimeStep(self.timestep)
        self.calc.SetFieldDiffusionMethod("EntropyFieldDynamics")

        # The output settings
        # The archive plugin is used in Culgi 6
        self.calc.SetSaveInstantResultsOn("SCMBondEnergy,SCMBendingEnergy,SCMTorsionEnergy,SCMNBEnergy,SCMElectrostaticsEnergy,SCMExternalPotentialEnergy,SCMPotentialEnergy,GCDensityInhomogeneity,GCTotalFreeEnergy,GCIdealFreeEnergy,GCContactFreeEnergy,GCCompressibilityFreeEnergy,GCElectrostaticFreeEnergy,CouplingEnergy")
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

        # Calculator parameters and simulation engine
        self.params = self.calc.GetParametersCmds()
        self.params_GC = self.params.GetGCParametersCmds()
        self.params_SC = self.params.GetSCMParametersCmds()
        self.engine = self.calc.GetDDFTEngineCmds()
        self.engine.SetMaxIterations(maxiters)

        # Parameters concerning mobility
        # Before phase 7 nanoparticle are single beads (SC molecules with 'P' beads)
        # From phase 8, nanoparticles are colloids with two different kinds of beads,
        #  but the diffusion constant is set for the whole colloid
        # Phase 8 assumes zero angular diffusion, but phase 9 already does not
        # From phase 9 the colloids can rotate
        self.calc.SetTemperature(self.temperature)
        self.params_GC.SetExpansionParameter("A", self.expansion)
        self.params_GC.SetExpansionParameter("B", self.expansion)
        if (self.phase <= 7):
            self.params_SC.SetDiffusionFactor("P", self.mobility)
        else:
            for i in range(self.population):
                self.nanoparticles[i].SetConstVelocity('Z', 0.0)
                self.nanoparticles[i].SetDiffusionFactor(self.mobility)
                if self.phase == 8:
                    self.nanoparticles[i].SetConstAngularVelocity(0.0, 0.0, 0.0)
                else:
                    self.nanoparticles[i].SetConstAngularVelocity('X', 0.0)
                    self.nanoparticles[i].SetConstAngularVelocity('Y', 0.0)
                    self.nanoparticles[i].SetRotationalFactors(0.0, 0.0, self.mobility)

        # Parameters concerning interactions
        # From phase 8, nanoparticle are colloids with two different types of beads
        self.calc.SetKappa(self.kappa)
        self.params_GC.SetBeadVolume("A", self.beadvolume)
        self.params_GC.SetBeadVolume("B", self.beadvolume)
        self.params_GC.SetChi("A", "B", self.nchi/self.bcp.GetNumBeads())
        if self.phase <= 7:
            self.params.SetBeadFieldCoupling("P", "A", self.ca)
            self.params.SetBeadFieldCoupling("P", "B", self.cb)
            self.params_SC.SetA("P", "P", self.a)
        else:
            self.params.SetBeadFieldCoupling("C", "A", self.ca)
            self.params.SetBeadFieldCoupling("C", "A", self.ca)
            self.params.SetBeadFieldCoupling("S", "B", self.cb)
            self.params.SetBeadFieldCoupling("S", "B", self.cb)
            self.params_SC.SetA("C", "C", self.a)
            self.params_SC.SetA("C", "S", self.a)
            self.params_SC.SetA("S", "S", self.a)

        # Output path and file name
        self.path = getpath(self)
        self.name = getname(self)
        self.calc.SetOutputFileName(self.name)

    def getcenterofmass(self, name, num=0):

        com = self.box.GetSoftCoreMoleculeCmds(name, num).GetCenterOfMass()
        return com.GetX(), com.GetY(), com.GetZ()

    def getorderparameter(self, x, y):

        return self.bcp_A.GetGaussianInterpolatedValue(x,y,0.0) - self.bcp_B.GetGaussianInterpolatedValue(x,y,0.0)

    def gettotalfield(self, x, y):

        return self.bcp_total.GetGaussianInterpolatedValue(x,y,0.0)

    def run(self):

        # From phase 4, we will be doing some equilibration before the actual simulation
        if self.phase >= 4:

            # Turn off all archives
            if Culgi.GetVersionID()[0] == "5":
                self.calc.SetSaveRelativeDensityFieldsOff()
                self.calc.SetSaveCoordinatesOff()
            if Culgi.GetVersionID()[0] == "6":
                self.calc.RemovePlugin(self.archive)
            self.calc.SetSaveInstantResultsOff()

            # No phase separation at all during equilibration
            self.params_GC.SetChi("A", "B", 0.0)

            # From phase 7, equilibrate the nanoparticles, too, first
            # Nanoparticles should not feel the field here at all, which is not updated
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
                self.calc.Run(int(Teq_np/self.timestep))

            # Now equilibrate fields with only base coupling, keeping NP positions fixed.
            if self.phase <= 7:
                self.params_SC.SetDiffusionFactor("P", 0)
                self.params.SetBeadFieldCoupling("P", "A", self.ca)
                self.params.SetBeadFieldCoupling("P", "B", self.ca)
            else:
                for i in range(self.population):
                    self.nanoparticles[i].SetDiffusionFactor(0.0)
                self.params.SetBeadFieldCoupling("C", "A", self.ca)
                self.params.SetBeadFieldCoupling("C", "B", self.ca)
                self.params.SetBeadFieldCoupling("S", "A", self.ca)
                self.params.SetBeadFieldCoupling("S", "B", self.ca)
            self.calc.SetFieldDiffusionOn()
            self.calc.SetBeadDiffusionOff()
            self.calc.Run(int(Teq_field/self.timestep))

            # Set relevant parameters back to the original values
            self.params_GC.SetChi("A", "B", self.nchi/self.bcp.GetNumBeads())
            if self.phase <= 7:
                self.params_SC.SetDiffusionFactor("P", self.mobility)
                self.params.SetBeadFieldCoupling("P", "A", self.ca)
                self.params.SetBeadFieldCoupling("P", "B", self.cb)
            else:
                for i in range(self.population):
                    self.nanoparticles[i].SetDiffusionFactor(self.mobility)
                self.params.SetBeadFieldCoupling("C", "A", self.ca)
                self.params.SetBeadFieldCoupling("C", "B", self.cb)               
                self.params.SetBeadFieldCoupling("S", "A", self.ca)
                self.params.SetBeadFieldCoupling("S", "B", self.cb)               
            self.calc.SetBeadDiffusionOn()

            # Turn archiving back on
            if Culgi.GetVersionID()[0] == "5":
                self.calc.SetSaveRelativeDensityFieldsOn()
                self.calc.SetSaveCoordinatesOn()
            if Culgi.GetVersionID()[0] == "6":
                self.calc.AddPlugin(self.archive)
            self.calc.SetSaveInstantResultsOn("SCMBondEnergy,SCMBendingEnergy,SCMTorsionEnergy,SCMNBEnergy,SCMElectrostaticsEnergy,SCMExternalPotentialEnergy,SCMPotentialEnergy,GCDensityInhomogeneity,GCTotalFreeEnergy,GCIdealFreeEnergy,GCContactFreeEnergy,GCCompressibilityFreeEnergy,GCElectrostaticFreeEnergy,CouplingEnergy")

        self.calc.Run(self.nsteps)

