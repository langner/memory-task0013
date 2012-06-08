"""Script to run task0013 simulations."""


import os

from pyculgi import *

from systems import *
from task0013 import *


def printnow(fname, content, mode='a'):
    """Print to a file without delay, so make sure to flush and close right now."""

    f = open(fname, mode)
    print >>f, content
    f.flush()
    f.close()


if __name__ == "__main__":

    # Set log file for this run job.
    flogname, flogpath = None, None
    if os.environ.has_key('PBS_JOBID'):
        flogname = 'run_%s.log' %os.environ['PBS_JOBID']

    # Loop over all phases flagged to run.
    for p in phases_to_run:

        # Loop over all parameter sets in this phase.
        for s in systems[p]:

            # Remember the current base directory we start from.
            basedir = os.path.abspath(os.curdir)

            # Create the simulation object.
            sim = Simulation(p, *s)

            # Create the output directory if needed.
            outdir = "%s/%s" %(basedir, sim.path)
            if Culgi.GetProcRank() == 0 and not os.path.exists(outdir):
                os.makedirs(outdir)

            # Move to the output directory.
            Culgi.MpiBarrier()
            os.chdir(outdir)

            # Skip this parameter set if output file already exists.
            # Must return to base directory here, too!
            if os.path.exists("%s.out" %sim.name):
                os.chdir(basedir)
                continue

            # Full path to log file.
            if flogname:
                flogpath = "%s/%s" %(basedir,flogname)

            # Print a log message that the simulation is starting.
            line = "Path: %s Output file: %s" %(sim.path,sim.name)
            if Culgi.GetProcRank() == 0:
                if flogpath:
                    printnow(flogpath, line)
                else:
                    print line

            # Setup simulation now (this can take some time), no earlier.
            sim.setup()

            # Run the simulation.
            Culgi.MpiBarrier()
            sim.run()

            # Change back to base directory.
            os.chdir(basedir)
