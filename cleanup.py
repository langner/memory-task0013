import glob
import os
import sys

from simulation import *
from systems import *


if __name__ == "__main__":

    if len(sys.argv) == 1:
        phases_to_clean = phases
    else:
        phases_to_clean = [int(p) for p in sys.argv[1:]]

    for p in [p for p in phases_to_clean if os.path.exists("phase%i" %p)]:

        print "Checking phase %i..." %p

        sims_to_keep = [Simulation(p, *s) for s in systems[p]]
        outs_to_keep = [s.outpath for s in sims_to_keep]

        outpattern = "phase%i/*x*x*_A*B*_bv*/temp*_exp?.??_den?.?_pop*/k*_nchi*/t*.out" %p
        outs = glob.glob(outpattern)
        sims_ondisk = [loadpath(out) for out in outs]
        for sim in sims_ondisk:
            out = sim.outpath
            if not out in outs_to_keep:
                outroot = out[:-4]
                to_remove = glob.glob(outroot+"*")
                print "Removing %s* (%i files)" %(outroot, len(to_remove))
                for fn in to_remove:
                    os.remove(fn)
                try:
                    os.removedirs(sim.path)
                except OSError:
                    print "Simulation directory remains (not empty)."

            