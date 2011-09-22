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

    if not "Time used" in open(fout).read():
        print "This simulation has not finished."
        sys.exit(1)

    sim = loadpath(fout)

    archive = UtilitiesManager.CreateArchiveReader()
    archive.SetSystem(sim.box, fcga, fcsa)
    nframes = archive.GetNumberOfFrames()
    beads = [sim.box.GetSoftCoreMoleculeCmds("np",np).GetBeadCmds(0) for np in range(sim.population)]

    cga = []
    csa = []

    print "init:", time.time() - t

    T1, T2, T3 = 0.0, 0.0, 0.0
    for nframe in range(1,nframes+1):

        t = time.time()
        archive.LoadFrame(nframe)
        T1 += time.time() - t

        t = time.time()
        coords = [b.GetCoordinates() for b in beads]
        coords = [[c.GetX(),c.GetY(),c.GetZ()] for c in coords]
        csa.append([coords])
        T2 += time.time() - t

        t = time.time()
        A = [[sim.bcp_A.GetValue(i,j,0) for j in range(64)] for i in range(64)]
        B = [[sim.bcp_B.GetValue(i,j,0) for j in range(64)] for i in range(64)]
        cga.append([A,B])
        T3 += time.time() - t

    print "archive:", T1, "coords:", T2, "fields:", T3

    t = time.time()
    ctf = numpy.loadtxt(fctf, skiprows=22)
    print "energy:", time.time() - t

    t = time.time()
    numpy.save(fname+".ctf.npy", ctf)
    numpy.save(fcsa+".npy", csa)
    numpy.save(fcga+".npy", cga)
    print "save:", time.time() - t