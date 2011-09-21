from common import *
from systems import *


snapsfreq = 1.0
energfreq = 0.1
maxiters = 25000
parameters = [  "size", "polymer", "beadvolume", "density", "nchi",
                "kappa", "temperature", "expansion",
                "ca", "cb", "a", "mobility", "population",
                "timestep", "totaltime" ]


def printnow(fname, content, mode='a'):
    f = open(fname, mode)
    print >>f, content
    f.flush()
    f.close()

def getpath(sim):
    X,Y,Z = sim.size
    pol = sim.polymer
    bv = sim.beadvolume
    den = sim.density
    pop = sim.population
    temp = sim.temperature
    exp = sim.expansion
    return "%ix%ix%i_%s_bv%.2f/temp%.2f_exp%.2f_den%.1f_pop%i" %(X,Y,Z,pol,bv,temp,exp,den,pop)

def getname(sim):
    k = sim.kappa
    nchi = sim.nchi
    ca = sim.ca
    cb = sim.cb
    a = sim.a
    mob = sim.mobility
    if sim.population == 0:
        name = "k%.1f_nchi%.1f" %(k,nchi)
    elif sim.population == 1:
        name = "k%.1f_nchi%.1f_ca%.1f_cb%.1f_mob%.2f" %(k,nchi,ca,cb,mob)
    else:
        name = "k%.1f_nchi%.1f_ca%.1f_cb%.1f_mob%.2f_a%.1f" %(k,nchi,ca,cb,mob,a)
    return name    

def loadpath(path):

    outpath,outfile = os.path.split(path.strip())
    system,model = outpath.split('/')

    size,pol,bv = system.split('_')
    size = map(int,size.split('x'))
    bv = float(bv[2:])

    temp,exp,den,pop = model.split('_')
    temp = float(temp[4:])
    exp = float(exp[3:])
    den = float(den[3:])
    pop = int(pop[3:])

    outname = outfile[:-4]
    if pop == 0:
        k,nchi = outname.split('_')
        k = float(k[1:])
        nchi = float(nchi[4:])
        ca = 0.0
        cb = 0.0
        mob = 0.0
        a = 0.0
    elif pop ==1:
        k,nchi,ca,cb,mob = outname.split('_')
        k = float(k[1:])
        nchi = float(nchi[4:])
        ca = float(ca[2:])
        cb = float(cb[2:])
        mob = float(mob[3:])
        a = 0.0
    else:
        k,nchi,ca,cb,mob,a = outname.split('_')
        k = float(k[1:])
        nchi = float(nchi[4:])
        ca = float(ca[2:])
        cb = float(cb[2:])
        mob = float(mob[3:])
        a = float(a[1:])

    sim = simulation(size,pol,bv,temp,exp,den,pop,k,nchi,ca,cb,a,mob,0.01,10)
    sim.outpath = outpath
    sim.outfile = outfile
    sim.outname = outname
    sim.setup()

    return sim


class simulation:

    def __init__(self, size, polymer, beadvolume, temperature, expansion, density, population, kappa, nchi, ca, cb, a, mobility, timestep, totaltime):

        # Set all the parameters.
        for parameter in parameters:
            setattr(self, parameter, eval(parameter))

        # Some derived parameters.
        self.dexcluded = self.ca - self.kappa
        self.dcoupling = self.cb - self.ca
        self.volume = size[0]*size[1]
        self.npvolume = self.ca/(self.kappa*self.density+1)
        self.effective_density = self.density - self.population*self.npvolume/self.volume

    def setup(self):

        # These follow from the time step.
        self.nsteps = int(self.totaltime/self.timestep)
        self.efreq = int(energfreq/self.timestep)
        self.sfreq = int(snapsfreq/self.timestep)

        # The diblock copolymer.
        self.bcp = Palette.CreateGaussianChain("bcp", self.polymer, *self.size)
        self.bcp.CreateHomogeneousRelativeDensity(self.density)
        self.bcp_total = self.bcp.GetRelativeDensityFieldCmds()
        self.bcp_A = self.bcp.GetBeadTypeRelativeDensityFieldCmds("A")
        self.bcp_A.SetDisplayClampLevels(0.0, self.density)
        self.bcp_B = self.bcp.GetBeadTypeRelativeDensityFieldCmds("B")
        self.bcp_B.SetDisplayClampLevels(0.0, self.density)

        # The nanoparticle prototype.
        self.np = Palette.CreateSoftCoreMolecule("np", "P")

        # The simulation box.
        self.box = Palette.CreateMesoBox("box", *self.size)
        self.box.AddGaussianChain(self.bcp)
        self.box.CopySoftCoreMolecules(self.np, self.population)

        # The actual nanoparticles inside the simulation box.
        self.nanoparticles = []
        for i in range(self.population):
            self.nanoparticles.append(self.box.GetSoftCoreMoleculeCmds("np", i))
            X = self.nanoparticles[i].GetCenterOfMass().GetX()
            Y = self.nanoparticles[i].GetCenterOfMass().GetY()
            Z = 0.5
            self.nanoparticles[i].SetCenterOfMass(X, Y, Z)

        # The calculator.
        self.calc = CalculatorsManager.CreateMBFHybridCalculator()
        self.calc.SetSystem(self.box)
        self.calc.SetTimeStep(self.timestep)
        self.calc.SetFieldDiffusionMethod("EntropyFieldDynamics")

        # The output settings.
        self.calc.SetSaveRelativeDensityFieldsOn()
        self.calc.SetSaveCoordinatesOn()
        self.calc.SetSaveInstantResultsOn("SCMBondEnergy,SCMBendingEnergy,SCMTorsionEnergy,SCMNBEnergy,SCMElectrostaticsEnergy,SCMExternalPotentialEnergy,SCMPotentialEnergy,GCDensityInhomogeneity,GCTotalFreeEnergy,GCIdealFreeEnergy,GCContactFreeEnergy,GCCompressibilityFreeEnergy,GCElectrostaticFreeEnergy,CouplingEnergy")
        self.calc.SetInstantResultsSaveFrequency(self.efreq)
        self.calc.SetRelativeDensityFieldSaveFrequency(self.sfreq)
        self.calc.SetCoordinateSaveFrequency(self.sfreq)

        # Calculator parameters and simulation engine.
        self.params = self.calc.GetParametersCmds()
        self.params_GC = self.params.GetGCParametersCmds()
        self.params_SC = self.params.GetSCMParametersCmds()
        self.engine = self.calc.GetDDFTEngineCmds()
        self.engine.SetMaxIterations(maxiters)

        # Parameters concerning mobility.
        self.calc.SetTemperature(self.temperature)
        self.params_GC.SetExpansionParameter("A", self.expansion)
        self.params_GC.SetExpansionParameter("B", self.expansion)
        self.params_SC.SetDiffusionFactor("P", self.mobility)

        # Parameters concerning interactions.
        self.calc.SetKappa(self.kappa)
        self.params_GC.SetBeadVolume("A", self.beadvolume)
        self.params_GC.SetBeadVolume("B", self.beadvolume)
        self.params_GC.SetChi("A", "B", self.nchi/self.bcp.GetNumBeads())
        self.params.SetBeadFieldCoupling("P", "A", self.ca)
        self.params.SetBeadFieldCoupling("P", "B", self.cb)
        self.params_SC.SetA("P", "P", self.a)

        # Output path and file name.
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

        self.calc.Run(self.nsteps)
