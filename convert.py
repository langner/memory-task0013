"""Script to convert Culgi archive files to npy archives."""


import sys
import time

import numpy

from pyculgi import *
from task0013 import *


if __name__ == "__main__":

    t = time.time()

    # First argument must be the output path we want to load as a simulation object.
    fout = sys.argv[1]

    # Make sure that simulation has finished.
    # Note that from phase 4, there are two runs for each simulation.
    if not "Time used" in open(fout).read().strip().split('\n')[-1]:
        print "This simulation has not finished."
        sys.exit(1)

    # Load the simulation object now, and do the setup step, too.
    sim = loadpath(fout, setup=True)

    # Construct the other file names, but the csa file only if there are nanoparticles.
    fname = fout[:-4]
    fcga = fname + ".cga"
    fctf = fname + "_Inst.ctf"
    if sim.population > 0:
        fcsa = fname + ".csa"

    # Create the archive object.
    archive = UtilitiesManager.CreateArchiveReader()
    if sim.population > 0:
        archive.SetSystem(sim.box, fcga, fcsa)
    else:
        archive.SetSystem(sim.box, fcga)

    # Get the number of frames in this archive.
    nframes = archive.GetNumberOfFrames()

    # Bead command objects, sorted per bead index, where Nbeads is number of beads per nanoparticle.
    if sim.population > 0:
        Nbeads = sim.nanoparticles[0].GetNumberOfBeads()
        beads = [[np.GetBeadCmds(i) for np in sim.nanoparticles] for i in range(Nbeads)]

    # Create the csa and cga arrays beforehand, which should save memory in the end.
    # There are always two kinds of fields, but not number of beads. Do not lump beads
    #   together, so there are always the same number of beads for each index and
    #   consequently the array is nicely shaped.
    cga = numpy.zeros((nframes,2, sim.size[0],sim.size[1]))
    if sim.population > 0:
        csa = numpy.zeros((nframes,Nbeads,sim.population,3))

    print "init:", time.time() - t

    T1, T2, T3 = 0.0, 0.0, 0.0
    for nframe in range(nframes):

        t = time.time()
        archive.LoadFrame(nframe+1)
        T1 += time.time() - t

        # Set the bead coordinates for this frame.
        # The second index is bead type.
        t = time.time()
        if sim.pop > 0:
            coords = [[bb.GetCoordinates() for bb in b] for b in beads]
            for i in range(Nbeads):
                csa[nframe,i] = [[c.GetX(),c.GetY(),c.GetZ()] for c in coords[i]]
        T2 += time.time() - t

        # Set the field density values for this frame.
        t = time.time()
        A = [[sim.bcp_A.GetValue(i,j,0) for j in range(sim.size[1])] for i in range(sim.size[0])]
        B = [[sim.bcp_B.GetValue(i,j,0) for j in range(sim.size[1])] for i in range(sim.size[0])]
        cga[nframe][0] = A
        cga[nframe][1] = B
        T3 += time.time() - t

    print "archive:", T1, "coords:", T2, "fields:", T3

    # Get the energy values from the ctf file.
    t = time.time()
    ctf = numpy.loadtxt(fctf, skiprows=22)
    print "energy:", time.time() - t

    # Save the arrays.
    t = time.time()
    numpy.save(fname+".ctf.npy", ctf)
    numpy.save(fcga+".npy", cga)
    if sim.population > 0:
        numpy.save(fcsa+".npy", csa)
    print "save:", time.time() - t