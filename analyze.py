import gzip
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

from scipy import mean
from scipy import var
from scipy.interpolate import RectBivariateSpline
from scipy.spatial.distance import pdist
from scipy.stats import kurtosis
from scipy.stats import skew


# Base frames and number of samples to analyze, by phase
PhaseFrames = {
    1 : ( [1, 11, 101, 1001, 10001], 16 ),
    2 : ( [1, 11, 101, 1001, 10001], 16 ),
    3 : ( [1, 11, 101, 1001, 10001], 16 ),
    4 : ( [1, 11, 101, 1001, 10001], 16 ),
    5 : ( [1, 11, 101, 1001, 10001], 16 ),
    6 : ( [1, 11, 51, 101, 501, 1001, 5001, 10001, 50001], 16 ),
    7 : ( [1, 11, 21, 51, 101, 201, 501, 1001, 2001, 5001, 10001], 10 ),
}

# Order of ctf columns
# Need also a version for neat systems
# We choose the version later, once population is known
energies_nps = ['nonbonded', 'inhomo', 'ideal', 'contact', 'compress', 'coupling']
columns_nps = dict(zip(energies_nps, [4, 8, 10, 11, 12, 14]))
energies_neat = ['inhomo', 'ideal', 'contact', 'compress']
columns_neat = dict(zip(energies_neat, [1, 3, 4, 5]))


class Analysis():

    def __init__(self, fout):
        """ Load files and initialization tasks """

        # Remove extension from Culgi output file
        self.fout = fout
        self.fname = fout[:-4]

        # Extract phase from the path
        self.phase = int(self.fname.split('/')[0][5:])

        # Extract NP population from the path
        self.pop = int(self.fname.split('/')[2].split('_')[-1][3:])

        # Choose columns/energie to use for ctf file
        if self.pop == 0:
            self.energies = energies_neat
            self.columns = columns_neat
        else:
            self.energies = energies_nps
            self.columns = columns_nps

        # Build ctf and archive file paths from the root
        # Empty csa file if no NPs in the system
        self.fctf = self.fname + ".ctf.npy.gz"
        self.fcga = self.fname + ".cga.npy.gz"
        self.fcsa = (self.fname + ".csa.npy.gz")*(self.pop > 0)

        # Make sure that simulations has finished
        # Note that from phase 4, there are two runs for each simulation
        if not "Time used" in open(self.fout).read().strip().split('\n')[-1]:
            print "This simulation has not finished."
            sys.exit(1)

        # Load NumPy archives
        self.ctf = load(gzip.open(self.fctf))
        self.cga = load(gzip.open(self.fcga))

        # Load csa archive only if there are NPs
        if self.pop > 0:

            self.csa = load(gzip.open(self.fcsa))

            # Due to a stupid mistake in the convert script, the shape sometimes needs fixing
            # Insert missing dimension, which represents the type of bead
            if len(self.csa.shape) == 3:
                self.csa = self.csa.reshape((self.csa.shape[0],1,self.csa.shape[1],self.csa.shape[2]))

            # Sanity check for bead Z coordinates after phase1
            if not (self.phase == 1 or all(self.csa[:,:,:,2] == 0.5)):
                print "Not all Z coordinates are 0.5"
                sys.exit(1)

            # Also coordinates need to be shifted to coincide with the field
            # And proceed only with two coordinatse (since Z is constant)
            self.csa = (self.csa[:,:,:,:2] - 0.5) % 64.0

    def analyze_energy(self):
        """ Analyze energy and other scalars """

        T = time.time()

        # Number of points to keep
        # Note that npoints might be correct to within +/-1
        self.npoints = 1100
        self.freq = self.ctf.shape[0]/self.npoints
        if self.pop > 0:
            self.freq_csa = self.csa.shape[0]/self.npoints
        self.indices = self.ctf[::self.freq,0]
        self.npoints = self.ctf.shape[0]/self.freq
        self.nrange = range(self.npoints)

        # Instantaneous energies as found in ctf file
        self.instants = [self.ctf[::self.freq,self.columns[en]] for en in self.energies]

        # Average and standard deviation between points kept
        cols = [self.columns[en] for en in self.energies]
        E = [array([self.ctf[i*self.freq:(i+1)*self.freq,c] for i in self.nrange]) for c in cols]
        self.averages = [average(e, axis=1) for e in E]
        self.deviations = [std(e, axis=1) for e in E]

        # Statistics of NPs on the grid
        # We only need look at the X and Y coordinates, since Z is fixed
        # Mean - from 0 to 1, where 0.5 is the grid point
        # Variance - should ideally be 0.085
        # Skewness - should ideally be 0
        # Kurtosis - should ideally be -1.2
        if self.pop > 0:
            self.offsets = (self.csa+0.5 - (self.csa+0.5).astype(int))[::self.freq_csa]
            self.offsets = [[f(o, axis=None) for o in self.offsets[:,0]] for f in mean,var,skew,kurtosis]
            self.offsets = [array(o) for o in self.offsets]

        return time.time() - T

    def analyze_histograms(self):
        """ Analyze data that is expressed as histograms """

        T = time.time()

        # Generate list of frames to analyze, with context for averaging
        self.baseframes, self.nsamples = PhaseFrames[self.phase]
        self.frames = hstack([0]+[range(bf,bf+self.nsamples) for bf in self.baseframes])
        self.nframes = len(self.frames)

        # Leave only the cga frames we want
        # Do the same for coordinates if there are in fact NPs
        # Note that since the cga/csa attributes are changed,
        #  we should do this ONLY after energy analysis is finished,
        #  because it has uses more frames
        self.cga = self.cga[self.frames]
        if self.pop > 0:
            self.csa = self.csa[self.frames]

        # Maximum values, ranges and bins to be used in histograms
        self.fmax = 1.5
        self.rmax = 16.0
        self.totalrange = (0.0,self.fmax)
        self.orderrange = (-self.fmax,self.fmax)
        self.rrange = [0.0,self.rmax]
        self.nbins_f = 256
        self.nbins_rad = 1024
        self.nbins_res = 256

        # Histograms of field values and order parameters
        # Chose some arbitrary range for the bins, so to catch all values
        self.totals = self.cga[:,0] + self.cga[:,1]
        self.orders = self.cga[:,0] - self.cga[:,1]
        self.hist_field_total = array([histogram(f.flatten(), bins=self.nbins_f, range=self.totalrange)[0] for f in self.totals])
        self.hist_field_order = array([histogram(f.flatten(), bins=self.nbins_f, range=self.orderrange)[0] for f in self.orders])
        self.hist_field_shape = (self.nframes,1,self.nbins_f)

        # Radial distribution histograms (only if there are NPs)
        # Disregard PBC, we are interested in short distances mostly
        # Therefore, generate the histogram only up to a distance of 10-20 units
        if self.pop > 0:
            self.pds = [pdist(c) for c in self.csa[:,0]]
            self.hist_radial = array([histogram(pd, bins=self.nbins_rad, range=self.rrange)[0] for pd in self.pds])
            self.hist_radial_shape = (self.nframes,1,self.nbins_rad)

        # Residual field totals and order parameters at particle positions
        if self.pop > 0:
            self.edge = range(64)
            self.splines = [[RectBivariateSpline(self.edge,self.edge,ff,kx=1,ky=1) for ff in f] for f in self.cga]
            self.values_res = array([[s.ev(c[:,0],c[:,1]) for s in self.splines[i]] for i,c in enumerate(self.csa[:,0])])
            self.totals_res = self.values_res[:,0,:] + self.values_res[:,1,:]
            self.orders_res = self.values_res[:,0,:] - self.values_res[:,1,:]
            self.hist_residual_total = array([histogram(t, bins=self.nbins_res, range=self.totalrange)[0] for t in self.totals_res])
            self.hist_residual_order = array([histogram(t, bins=self.nbins_res, range=self.orderrange)[0] for t in self.orders_res])
            self.hist_residual_shape = (self.nframes,1,self.nbins_res)

        return time.time() - T

    def save(self):
        """ Save the analyzed data to archives """

        T = time.time()

        if hasattr(self, 'instants'):

            # Stack vectors of energies and other instantaneous scalar results
            # Keep the indices together with all these, for reference
            # We might need to drop the last elements in instantaneous results
            # Note that the saved array is smaller when np NPs are in the system
            if len(self.instants[0]) > len(self.averages[0]):
                n = len(self.averages[0])
                self.indices = self.indices[:-1]
                self.instants = [i[:-1] for i in self.instants]
                if self.pop > 0:
                    self.offsets = [o[:-1] for o in self.offsets]
            if self.pop > 0:
                energy = vstack([self.indices]+self.instants+self.averages+self.deviations+self.offsets).transpose()
            else:
                energy = vstack([self.indices]+self.instants+self.averages+self.deviations).transpose()

            # Save the instananeous results as energies
            save(self.fname+".energy.npy", energy)

        if hasattr(self, 'frames'):

            # Save also the frames used in histograms for reference
            save(self.fname+".hist-frames.npy", self.frames)
            
            # Save histograms of field values
            total = self.hist_field_total.reshape(self.hist_field_shape)
            order = self.hist_field_order.reshape(self.hist_field_shape)
            save_field = concatenate([total,order], axis=1)
            save(self.fname+".hist-field.npy", save_field)

            # Save histograms of radial distributions and residual fields
            if self.pop > 0:
                save(self.fname+".hist-radial.npy", array(self.hist_radial))
                total = self.hist_residual_total.reshape(self.hist_residual_shape)
                order = self.hist_residual_order.reshape(self.hist_residual_shape)
                save_residual = concatenate([total, order], axis=1)
                save(self.fname+".hist-residual.npy", save_residual)

        return time.time() - T


if __name__ == "__main__":

    # First argument must be the output file
    fout = sys.argv[1]

    # Load data and initiate
    # Here we must time manually
    T = time.time()
    a = Analysis(fout)
    print "init: %.2fs" %(time.time()-T)

    # Do energy analysis
    if "energy" in sys.argv:
        print "energy: %.2fs" %a.analyze_energy()

    # Do histogram analysis
    if "histograms" in sys.argv:
        print "histograms: %.2fs" %a.analyze_histograms()

    # Save archives
    if "energy" in sys.argv or "histgrams" in sys.argv:
        print "save:, %.2fs" %a.save()
