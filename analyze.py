import os
import sys
import time

from numpy import all
from numpy import arange
from numpy import array
from numpy import average
from numpy import concatenate
from numpy import histogram
from numpy import hstack
from numpy import load
from numpy import max
from numpy import round
from numpy import save
from numpy import sqrt
from numpy import std
from numpy import vstack
from numpy.linalg import norm
from scipy.interpolate import RectBivariateSpline
from scipy.spatial.distance import pdist


# Base frames and number of samples to analyze, by phase.
PhaseFrames = {
    "phase1" : ( [1, 11, 101, 1001, 10001], 16 ),
    "phase2" : ( [1, 11, 101, 1001, 10001], 16 ),
    "phase3" : ( [1, 11, 101, 1001, 10001], 16 ),
    "phase4" : ( [1, 11, 101, 1001, 10001], 16 ),
    "phase5" : ( [1, 11, 101, 1001, 10001], 16 ),
    "phase6" : ( [1, 11, 51, 101, 501, 1001, 5001, 10001, 50001], 16 ),
}

# Order of ctf columns.
energies = ['nonbonded', 'inhomo', 'ideal', 'contact', 'compress', 'coupling']
columns = dict(zip(energies,[4, 8, 10, 11, 12, 14]))

if __name__ == "__main__":

    fout = sys.argv[1]
    fname = fout[:-4]
    phase = fname.split('/')[0]
    fctf = fname + ".ctf.npy"
    fcga = fname + ".cga.npy"
    fcsa = fname + ".csa.npy"

    # Make sure that simulations has finished.
    # Note that from phase 4, there are two runs for each simulation.
    if not "Time used" in open(fout).read().strip().split('\n')[-1]:
        print "This simulation has not finished."
        sys.exit(1)

    # ###############
    # Energy analysis
    # ###############

    T = time.time()

    # Load NumPy archive.
    ctf = load(fctf)

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
    baseframes, nsamples = PhaseFrames[phase]
    frames = hstack([0]+[range(bf,bf+nsamples) for bf in baseframes])
    nframes = len(frames)

    # Load NumPy archives.
    csa = load(fcsa)
    cga = load(fcga)

    # Due to a stupid mistake in the convert script, the shape sometimes needs fixing.
    # Insert missing dimension, which represents the type of bead.
    if len(csa.shape) == 3:
        csa = csa.reshape((csa.shape[0],1,csa.shape[1],csa.shape[2]))

    # Leave only the csa and cga frames we want.
    csa = csa[frames]
    cga = cga[frames]

    # Sanity check for bead Z coordinates after task0013.
    if not (phase == "phase1" or all(csa[:,:,:,2] == 0.5)):
        print "Not all Z coordinates are 0.5."
        sys.exit(1)

    print 'archive (load):', time.time() - T

    T = time.time()

    # Histograms of field values and order parameters.
    # Chose some arbitrary range for the bins, so to catch all values.
    nbins = 256
    fmax = 1.5
    totals = cga[:,0] + cga[:,1]
    orders = cga[:,0] - cga[:,1]
    hist_field_total = array([histogram(f.flatten(), bins=nbins, range=(0.0,fmax))[0] for f in totals])
    hist_field_order = array([histogram(f.flatten(), bins=nbins, range=(-fmax,fmax))[0] for f in orders])
    hist_field_shape = (nframes,1,nbins)

    # Coordinates need to be shifted to coincide with the field.
    coords = (csa[:,0,:,:2] - 0.5) % 64.0

    # Radial distribution histograms.
    # Disregard PBC, we are interested in short distances mostly.
    # Therefore, generate the histogram only up to a distance of 10-20 units or so.
    nbins = 1024
    rmax = 16.0
    pds = [pdist(c) for c in coords]
    hist_radial = array([histogram(pd, bins=nbins, range=[0.0,rmax])[0] for pd in pds])
    hist_radial_shape = (nframes,1,nbins)

    # Residual field totals and order parameters at particle positions.
    nbins = 256
    edge = range(64)
    splines = [[RectBivariateSpline(edge,edge,ff,kx=1,ky=1) for ff in f] for f in cga]
    values_res = array([[s.ev(c[:,0],c[:,1]) for s in splines[i]] for i,c in enumerate(coords)])
    totals_res = values_res[:,0,:] + values_res[:,1,:]
    orders_res = values_res[:,0,:] - values_res[:,1,:]
    hist_residual_total = array([histogram(t, bins=nbins, range=(0.0,fmax))[0] for t in totals_res])
    hist_residual_order = array([histogram(t, bins=nbins, range=(-fmax,fmax))[0] for t in orders_res])
    hist_residual_shape = (nframes,1,nbins)

    print 'archive (other):', time.time() - T

    # ##################
    # Save analyzed data
    # ##################

    T = time.time()

    save(fname+".energy.npy", energy)
    save(fname+".hist-frames.npy", frames)
    
    save_field = concatenate([hist_field_total.reshape(hist_field_shape), hist_field_order.reshape(hist_field_shape)], axis=1)
    save(fname+".hist-field.npy", save_field)

    save(fname+".hist-radial.npy", array(hist_radial))

    save_residual = concatenate([hist_residual_total.reshape(hist_residual_shape), hist_residual_order.reshape(hist_residual_shape)], axis=1)
    save(fname+".hist-residual.npy", save_residual)

    print "save:", time.time() - T
