from simulation import *
from systems import *


# What we currently want to run.
phases = [ 4, 5 ]

if __name__ == "__main__":

    # Set log file.
    flogname, flogpath = None, None
    if os.environ.has_key('PBS_JOBID'):
        flogname = 'run_%s.log' %os.environ['PBS_JOBID']

    for phase in phases:

        for system in systems[phase]:

            # Create the simulation object.
            sim = simulation(phase, *system)
            sim.setup()

            # Base and output path. 
            basedir = os.path.abspath(os.curdir)
            if flogname:
                flogpath = "%s/%s" %(basedir,flogname)
            outpath = "%s/%s" %(basedir,sim.path)
            if Culgi.GetProcRank() == 0 and not os.path.exists(outpath):
                os.makedirs(outpath)

            # Change to directory for this model.
            Culgi.MpiBarrier()
            os.chdir(outpath)

            # Skip this parameter set if output file already exists.
            # Must return to base directory here, too!
            if os.path.exists("%s.out" %sim.name):
                os.chdir(basedir)
                continue

            # Print log message.
            line = "Path: %s Output file: %s" %(sim.path,sim.name)
            if Culgi.GetProcRank() == 0:
                if flogpath:
                    printnow(flogpath, line)
                else:
                    print line

            # Run it.
            Culgi.MpiBarrier()
            sim.run()

            # Change back to base directory.
            os.chdir(basedir)
