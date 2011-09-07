import gzip
import sys
import time

from numpy import arange
from numpy import array
from numpy import average
from numpy import histogram
from numpy import hstack
from numpy import load
from numpy import round
from numpy import savez
from numpy import sqrt
from numpy import std
from numpy import vstack
from numpy.linalg import norm
from scipy.interpolate import RectBivariateSpline
from scipy.spatial.distance import pdist


# Order of ctf columns.
energies = ['nonbonded', 'inhomo', 'ideal', 'contact', 'compress', 'coupling']
columns = dict(zip(energies,[4, 8, 10, 11, 12, 14]))

if __name__ == "__main__":

    # #######
    # Initial
    # #######

    T = time.time()

    # Filenames.
    fraw = sys.argv[1]
    fname = fraw[:-13]
    fcga = fname + ".cga"
    fcsa = fname + ".csa"
    fout = fname + ".out"

    # Load NumPy archive.
    npz = load(fraw)
    ctf = npz['ctf']
    cga = npz['cga'].swapaxes(0,1)
    csa = npz['csa'].swapaxes(0,1)

    print "init:", time.time() - T

    # ###############
    # Energy analysis
    # ###############

    T = time.time()

    # Number of points to keep.
    npoints = 1100
    freq = ctf.shape[0]/npoints
    indices = ctf[::freq,0]

    # Instantaneous energies as found in ctf file.
    instants = [ctf[::freq,columns[en]] for en in energies]

    # Average between points kept.
    averages = [[average(ctf[i*freq:(i+1)*freq,columns[en]]) for i in range(npoints+1)] for en in energies]
    averages = [array(avg) for avg in averages]

    # Standard deviation between points kept.
    deviations = [[std(ctf[i*freq:(i+1)*freq,columns[en]]) for i in range(npoints+1)] for en in energies]
    deviations = [array(dev) for dev in deviations]

    # Keep the indices together with the energies.
    energy = vstack([indices]+instants+averages+deviations).transpose()

    print 'energy:', time.time() - T

    # ################
    # Archive analysis
    # ################

    T = time.time()

    # Generate list of frames to analyze, with context for averaging.
    baseframes = [1, 11, 101, 1001, 10001]
    nsamples = 10
    frames = hstack([0]+[range(bf,bf+nsamples) for bf in baseframes])

    # Coordinates need to be shifted to coincide with the field.
    coords = (csa[frames][:,0,:,:2] - 0.5) % 64.0

    # Disregard PBC, we are interested in short distances mostly.
    pds = [pdist(c) for c in coords]
    hists_radials = array([histogram(pd, bins=1000, range=[0.0,32*sqrt(2.0)])[0] for pd in pds])

    # Field totals and order parameters at particle positions.
    fields = cga[frames]
    edge = range(64)
    splines = [[RectBivariateSpline(edge,edge,ff,kx=2,ky=2) for ff in f] for f in fields]
    values = array([[s.ev(c[:,0],c[:,1]) for s in splines[i]] for i,c in enumerate(coords)])
    hists_totals = array([histogram(t, bins=256) for t in values[:,0,:]+values[:,1,:]])
    hists_orders = array([histogram(t, bins=256) for t in values[:,0,:]-values[:,1,:]])

    print 'archive:', time.time() - T

    # ##################
    # Save analyzed data
    # ##################

    T = time.time()

    savez(fname+".data-analyzed.npz", energy=energy, hist_frames=frames,
        hists_radials=hists_radials, hists_totals=hists_totals, hists_orders=hists_orders)

    print "save:", time.time() - T
