import gzip
import os
import sys
import time

from numpy import all
from numpy import arange
from numpy import array
from numpy import average
from numpy import histogram
from numpy import hstack
from numpy import load
from numpy import max
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
    cga = npz['cga']
    csa = npz['csa']

    # Sanity check for bead Z coordinates after task0013.
    if not (os.path.abspath(os.curdir).split('/')[-1] == "task0013" or all(csa[:,:,:,2] == 0.5)):
        print "Not all Z coordinates are 0.5."
        sys.exit(1)

    print "init:", time.time() - T

    # ###############
    # Energy analysis
    # ###############

    T = time.time()

    # Number of points to keep.
    # Note that npoints might be correct by +/-1.
    npoints = 1100
    freq = ctf.shape[0]/npoints
    indices = ctf[::freq,0]
    npoints = ctf.shape[0]/freq

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
    nsamples = 16
    frames = hstack([0]+[range(bf,bf+nsamples) for bf in baseframes])

    # Leave only the csa and cga frames we want.
    csa = csa[frames]
    cga = cga[frames]

    # Histograms of field values and order parameters.
    # Chose some arbitrary range for the bins, so to catch all values.
    nbins = 256
    fmax = 1.5
    totals = cga[:,0] + cga[:,1]
    orders = cga[:,0] - cga[:,1]
    hists_totals = array([histogram(f.flatten(), bins=nbins, range=(0.0,fmax))[0] for f in totals])
    hists_orders = array([histogram(f.flatten(), bins=nbins, range=(-fmax,fmax))[0] for f in orders])

    # Coordinates need to be shifted to coincide with the field.
    coords = (csa[:,0,:,:2] - 0.5) % 64.0

    # Radial distribution histograms.
    # Disregard PBC, we are interested in short distances mostly.
    # Therefore, generate the histogram only up to a distance of 10-20 units or so.
    nbins = 1024
    rmax = 16.0
    pds = [pdist(c) for c in coords]
    hists_radials = array([histogram(pd, bins=nbins, range=[0.0,rmax])[0] for pd in pds])

    # Residual field totals and order parameters at particle positions.
    nbins = 256
    edge = range(64)
    splines = [[RectBivariateSpline(edge,edge,ff,kx=1,ky=1) for ff in f] for f in cga]
    values_res = array([[s.ev(c[:,0],c[:,1]) for s in splines[i]] for i,c in enumerate(coords)])
    totals_res = values_res[:,0,:] + values_res[:,1,:]
    orders_res = values_res[:,0,:] - values_res[:,1,:]
    hists_totals_res = array([histogram(t, bins=nbins, range=(0.0,fmax))[0] for t in totals_res])
    hists_orders_res = array([histogram(t, bins=nbins, range=(-fmax,fmax))[0] for t in orders_res])

    print 'archive:', time.time() - T

    # ##################
    # Save analyzed data
    # ##################

    T = time.time()

    savez(fname+".data-analyzed.npz", energy=energy, hist_frames=frames,
        hists_totals=hists_totals, hists_orders=hists_orders,
        hists_radials=hists_radials, hists_totals_res=hists_totals_res, hists_orders_res=hists_orders_res)

    print "save:", time.time() - T
