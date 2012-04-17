import glob
import shutil

from run import phases_to_run
from simulation import loadpath, simulation
from systems import systems


phases_to_clean = range(1,2)


if __name__ == "__main__":

    for p in phases_to_clean:

        sims_to_keep = [simulation(p, *s) for s in systems[p]]
        for s in sims_to_keep:
            s.setup()
        outs_to_keep = [s.outpath for s in sims_to_keep]

        outpattern = "phase%i/*x*x*_A*B*_bv?.??/temp?.??_exp?.??_den?.?_pop*/k*_nchi*/t*.out" %p
        outs = glob.glob(outpattern)
        sims_ondisk = [loadpath(out) for out in outs]
        for sim in sims_ondisk:

            if not sim.outpath in outs_to_keep:
                print "Removing %s..." %sim.path
                shutil.rmtree(sim.path)

            