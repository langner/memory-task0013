import sys
import time

import numpy

from pyculgi import *

from simulation import *


if __name__ == "__main__":

    t = time.time()

    # Load the simulation object
    fout = sys.argv[1]
    sim = loadpath(fout, setup=False)

    # Construct the other file names
    fname = fout[:-4]
    fcga = fname + ".cga"
    fctf = fname + "_Inst.ctf"
    if sim.population > 0:
        fcsa = fname + ".csa"

    # Make sure that simulation has finished
    # Note that from phase 4, there are two runs for each simulation
    if not "Time used" in open(fout).read().strip().split('\n')[-1]:
        print "This simulation has not finished."
        sys.exit(1)

    # Create the archive object
    archive = UtilitiesManager.CreateArchiveReader()
    if sim.population > 0:
        archive.SetSystem(sim.box, fcga, fcsa)
    else:
        archive.SetSystem(sim.box, fcga)
    nframes = archive.GetNumberOfFrames()

    # Bead command objects, sorted per bead index
    # Parameter Nbeads is number of beads per nanoparticle
    if sim.population > 0:
        Nbeads = sim.nanoparticles[0].GetNumberOfBeads()
        beads = [[np.GetBeadCmds(i) for np in sim.nanoparticles] for i in range(Nbeads)]

    # Create the arrays beforehand, which should save memory
    # There are always two kinds of fields, but not number of beads
    cga = numpy.zeros((nframes,2,64,64))
    if sim.population > 0:
        csa = numpy.zeros((nframes,Nbeads,sim.population,3))

    print "init:", time.time() - t

    T1, T2, T3 = 0.0, 0.0, 0.0
    for nframe in range(nframes):

        t = time.time()
        archive.LoadFrame(nframe+1)
        T1 += time.time() - t

        t = time.time()

        # Transform NP coordinates, only if there are any
        # Second index is bead type
        if sim.population > 0:
            coords = [[bb.GetCoordinates() for bb in b] for b in beads]
            for i in range(Nbeads):
                csa[nframe,i] = [[c.GetX(),c.GetY(),c.GetZ()] for c in coords[i]]

        T2 += time.time() - t

        t = time.time()
        A = [[sim.bcp_A.GetValue(i,j,0) for j in range(64)] for i in range(64)]
        B = [[sim.bcp_B.GetValue(i,j,0) for j in range(64)] for i in range(64)]
        cga[nframe][0] = A
        cga[nframe][1] = B
        T3 += time.time() - t

    print "archive:", T1, "coords:", T2, "fields:", T3

    t = time.time()
    ctf = numpy.loadtxt(fctf, skiprows=22)
    print "energy:", time.time() - t

    t = time.time()

    # Save the arrays
    numpy.save(fname+".ctf.npy", ctf)
    numpy.save(fcga+".npy", cga)
    if sim.population > 0:
        numpy.save(fcsa+".npy", csa)

    print "save:", time.time() - t