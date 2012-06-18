"""Script to clean up simulation files of runs that were removed from the system definitions.

The idea is that some runs will be removed after they are finished, either because
they did not turn out as nice as thought, or becuase they're redundant, or for any other reason.
"""


import glob
import os
import sys

from systems import *
from task0013 import *


if __name__ == "__main__":

    # Take phases to clean as arguments, or use all phases.
    if len(sys.argv) == 1:
        phases_to_clean = phases
    else:
        phases_to_clean = [int(p) for p in sys.argv[1:]]

    # Loop over phases (if the directory exists)
    for p in [p for p in phases_to_clean if os.path.exists("phase%i" %p)]:

        print "Checking phase %i..." %p

        # These are the simulation and output files we want to keep.
        sims_to_keep = [Simulation(p, *s) for s in systems[p]]
        outs_to_keep = [s.outpath for s in sims_to_keep]

        # Find all output files in this phase and load them as simulation objects.
        outpattern = "phase%i/*x*x*_A*B*_bv*/temp*_exp?.??_den?.?_pop*/k*_nchi*/t*.out" %p
        outs = glob.glob(outpattern)
        sims_ondisk = [loadpath(out, setup=False) for out in outs]

        # Loop over all simulation found on disk.
        for sim in sims_ondisk:

            # Remove all associated files the output file is not in the keep list.
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

            