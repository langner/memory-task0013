import gzip
import sys
import time

import numpy
from scipy import histogram
from scipy.interpolate import RectBivariateSpline
from scipy.spatial.distance import pdist


energies = ['nonbonded', 'inhomo', 'ideal', 'contact', 'compress', 'coupling']
columns = dict(zip(energies,[4, 8, 10, 11, 12, 14]))

if __name__ == "__main__":

    T = time.time()
    # Filenames.
    fraw = sys.argv[1]
    fname = fraw[:-13]
    fcga = fname + ".cga"
    fcsa = fname + ".csa"
    fout = fname + ".out"
    # Load NumPy archive.
    npz = numpy.load(fnpz)
    ctf = npz['ctf']
    cga = npz['cga'].swapaxes(0,1)
    csa = npz['csa'].swapaxes(0,1)
    print "init:", time.time() - T

    # Energy analysis.
    T = time.time()
    npoints = 1100
    freq = ctf.shape[0]/npoints
    indices = ctf[::freq,0]
    instants = [ctf[::freq,columns[en]] for en in energies]
    averages = [[numpy.average(ctf[i*freq:(i+1)*freq,columns[en]]) for i in range(npoints+1)] for en in energies]
    averages = [numpy.array(avg) for avg in averages]
    deviations = [[numpy.std(ctf[i*freq:(i+1)*freq,columns[en]]) for i in range(npoints+1)] for en in energies]
    deviations = [numpy.array(dev) for dev in deviations]
    energy = numpy.vstack([indices]+instants+averages+deviations).transpose()
    print 'energy:', time.time() - T

    # Archive analysis.
    T = time.time()
    coords = (csa[:,:,:,:2] - 0.5) / 64.0
    N = coords.shape[2]
    nbins = 1000
    frames = [0, 1, 2, 10, 20, 100, 200, 1000, 2000, 10000, 11000]
    hists_radials = numpy.array([histogram(pdist(coords[iframe,0], "euclidean"), bins=nbins)[0] for iframe in frames])
    print 'archive:', time.time() - T

    # Save analyzed data.
    T = time.time()
    numpy.savez(fname+".data-analyzed.npz", energy=energy, hists_radials=hists_radials)
    print "save:", time.time() - T
