import sys
import time

import numpy

from simulation import *


if __name__ == "__main__":

    t = time.time()

    fout = sys.argv[1]
    fname = fout[:-4]
    fcga = fname + ".cga"
    fcsa = fname + ".csa"
    fctf = fname + "_Inst.ctf"

    # Make sure that simulations has finished.
    # Note that from phase 4, there are two runs for each simulation.
    if not "Time used" in open(fout).read().strip().split('\n')[-1]:
        print "This simulation has not finished."
        sys.exit(1)

    sim = loadpath(fout)

    archive = UtilitiesManager.CreateArchiveReader()
    archive.SetSystem(sim.box, fcga, fcsa)
    nframes = archive.GetNumberOfFrames()
    beads = [sim.box.GetSoftCoreMoleculeCmds("np",np).GetBeadCmds(0) for np in range(sim.population)]

    # Create the arrays beforehand, which should save memory
    cga = numpy.zeros((nframes,2,64,64))
    csa = numpy.zeros((nframes,1,sim.population,3))

    print "init:", time.time() - t

    T1, T2, T3 = 0.0, 0.0, 0.0
    for nframe in range(nframes):

        t = time.time()
        archive.LoadFrame(nframe+1)
        T1 += time.time() - t

        t = time.time()

        # Transform NP coordinates, only if there are any
        if sim.population > 0:
            coords = [b.GetCoordinates() for b in beads]
            csa[nframe,0] = [[c.GetX(),c.GetY(),c.GetZ()] for c in coords]

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
    numpy.save(fcsa+".npy", csa)

    print "save:", time.time() - t