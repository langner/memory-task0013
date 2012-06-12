"""Script to analyze npy archives from simulations.

The npy archives were originally converted from
Culgi archives (text files) by another script.
The idea is that we only need this processed data for
presenting results, plus it then takes up much less space
than even the converted npy archives.
"""


import gzip
import os
import sys
import time

from numpy import all, arange, arctan2, array, average
from numpy import concatenate, histogram, hstack, load
from numpy import max, round, save, sqrt, std, vstack
from numpy.linalg import norm

from scipy import mean, var
from scipy.interpolate import RectBivariateSpline
from scipy.spatial.distance import pdist
from scipy.stats import kurtosis
from scipy.stats import skew

from systems import *


# Order of ctf columns for neat systems and system with nanoparticles.
energies_nps = ['nonbonded', 'inhomo', 'ideal', 'contact', 'compress', 'coupling']
columns_nps = dict(zip(energies_nps, [4, 8, 10, 11, 12, 14]))
energies_neat = ['inhomo', 'ideal', 'contact', 'compress']
columns_neat = dict(zip(energies_neat, [1, 3, 4, 5]))


class Analysis():

    def __init__(self, fout):
        """Load files and initialization certain data."""

        # Remove extension from Culgi output file.
        self.fout = fout
        self.fname = fout[:-4]

        # Extract phase and populations from the path.
        # We could load a simulation object, but this way we skip the dependency
        #   on the Culgi Python module.
        self.phase = int(self.fname.split('/')[0][5:])
        self.pop = int(self.fname.split('/')[2].split('_')[3][3:])

        # Choose columns/energies to use for ctf file.
        if self.pop == 0:
            self.energies = energies_neat
            self.columns = columns_neat
        else:
            self.energies = energies_nps
            self.columns = columns_nps

        # Build ctf and archive file paths from the root.
        # Empty csa file if no NPs in the system.
        self.fctf = self.fname + ".ctf.npy.gz"
        self.fcga = self.fname + ".cga.npy.gz"
        self.fcsa = (self.fname + ".csa.npy.gz")*(self.pop > 0)

        # Make sure that simulation has finished.
        # Note that from phase 4, there are two runs for each simulation.
        if not "Time used" in open(self.fout).read().strip().split('\n')[-1]:
            print "This simulation has not finished."
            sys.exit(1)

        # Load the NumPy archives as arrays.
        self.ctf = load(gzip.open(self.fctf))
        self.cga = load(gzip.open(self.fcga))

        # Load csa archive only if there are NPs.
        if self.pop > 0:

            self.csa = load(gzip.open(self.fcsa))

            # Due to a stupid mistake in the convert script, the shape sometimes needs fixing.
            # Namely, we need to insert a missing dimension, which represents the type of bead.
            if len(self.csa.shape) == 3:
                self.csa = self.csa.reshape((self.csa.shape[0],1,self.csa.shape[1],self.csa.shape[2]))

            # Sanity check for bead Z coordinates after phase 1.
            if self.phase > 1 and not all(self.csa[:,:,:,2] == 0.5):
                print "Not all Z coordinates are 0.5 (and they should)"
                sys.exit(1)

            # Coordinates need to be shifted to coincide with the field.
            # Also, proceed only with two coordinates (since Z is irrelevant).
            # Acutally, in phase 1 the third coordinates is relevant, but
            #   the later phases are of much more interest and this will save memory.
            self.csa = (self.csa[:,:,:,:2] - 0.5)

            # Apply PBC to bead coordinates.
            # After phase 7, we have some internal structure, so the modulo of a core bead
            #   implies the modulo of all related shell beads, and shell beads should not
            #   be subject to the modulo independently of core beads (thus the tripple loop).
            if self.phase <= 7:
                self.csa %= 64.0
            else:
                for snap in self.csa:
                    for icb, corebead in enumerate(snap[0]):
                        for ix in range(2):
                            if corebead[ix] >= 64.0 or corebead[ix] < 0.0:
                                snap[:,icb,ix] %= 64.0

        # Indices of beads to use for analysis.
        # After phase 7, there are both core and shell beads, and we want
        #   to analyze those separately.
        if self.phase <= 7:
            self.bead_ind = [[0]]
        else:
            self.bead_ind = [[0],range(1,9)]

    def analyze_energy(self):
        """Analyze energy and other scalars."""

        T = time.time()

        # Number of points to keep and the frequency with which to sample data.
        # Note that npoints might be correct to within +/-1.
        # If npoints is larger than number of frames, default the frequency to 1.
        self.npoints = 1100
        self.freq = self.ctf.shape[0]/self.npoints or 1
        if self.pop > 0:
            self.freq_csa = self.csa.shape[0]/self.npoints or 1
        self.indices = self.ctf[::self.freq,0]
        self.npoints = self.ctf.shape[0]/self.freq
        self.nrange = range(self.npoints)

        # Instantaneous energies as found in ctf file.
        self.instants = [self.ctf[::self.freq,self.columns[en]] for en in self.energies]

        # Average and standard deviation between points kept.
        cols = [self.columns[en] for en in self.energies]
        E = [array([self.ctf[i*self.freq:(i+1)*self.freq,c] for i in self.nrange]) for c in cols]
        self.averages = [average(e, axis=1) for e in E]
        self.deviations = [std(e, axis=1) for e in E]

        # Statistics of NPs on the grid (offsets).
        #   Mean        - from 0 to 1, where 0.5 is the grid point
        #   Variance    - should ideally be 1/12
        #   Skewness    - should ideally be 0
        #   Kurtosis    - should ideally be -6/5
        # From phase 8, nanoparticles are colloids, so offsets can be calculated for both the core
        #   and shell beads separately (all shell beads combined, and stacked right after core offsets).
        funcs = [mean, var, skew, kurtosis]
        if self.pop > 0:
            self.offsets = (self.csa+0.5 - (self.csa+0.5).astype(int))[::self.freq_csa]
            self.offsets = [[f(o, axis=None) for o in self.offsets[:,bi]] for f in funcs for bi in self.bead_ind]
            self.offsets = [array(o) for o in self.offsets]

        # Evaluate the rotation of nanoparticles and calculate statistics on the angles.
        # This makes sense starting from phase 9, and use bead with index 0 as reference, assuming
        #   that it started in the "north" position.
        # Note: we take the absolute value of the angle, to make the function continuous => angle goes up to pi/2.
        if self.phase > 8 and self.pop > 0:
            icenter = self.bead_ind[0][0]
            iref = self.bead_ind[1][0]
            self.angles = [abs(arctan2(snap[iref,:,1]-snap[icenter,:,1],snap[iref,:,0]-snap[icenter,:,0])) for snap in self.csa[::self.freq_csa]]
            self.offsets_angles = [[f(a, axis=None) for a in self.angles] for f in funcs]

        return time.time() - T

    def analyze_histograms(self):
        """Analyze data that we want to express as histograms."""

        T = time.time()

        # Generate list of frames to analyze, with context for averaging.
        self.baseframes, self.nsamples = phases_frames[self.phase]
        self.frames = hstack([0]+[range(bf,bf+self.nsamples) for bf in self.baseframes])
        self.nframes = len(self.frames)

        # Leave only the cga frames we want. Do the same for coordinates if there are in fact NPs.
        # Note that since the cga/csa attributes are changed, we should do this ONLY after
        #   energy analysis is finished, because there we use more frames.
        self.cga = self.cga[self.frames]
        if self.pop > 0:
            self.csa = self.csa[self.frames]

        # Maximum values, ranges and bins to be used in histograms.
        self.fmax = 1.5
        self.rmax = 16.0
        self.totalrange = (0.0,self.fmax)
        self.orderrange = (-self.fmax,self.fmax)
        self.rrange = [0.0,self.rmax]
        self.nbins_f = 256
        self.nbins_rad = 1024
        self.nbins_res = 256
        self.nbins_ang = 128

        # Histograms of field values and order parameters.
        self.totals = self.cga[:,0] + self.cga[:,1]
        self.orders = self.cga[:,0] - self.cga[:,1]
        self.hist_field_total = array([histogram(f.flatten(), bins=self.nbins_f, range=self.totalrange)[0] for f in self.totals])
        self.hist_field_order = array([histogram(f.flatten(), bins=self.nbins_f, range=self.orderrange)[0] for f in self.orders])
        self.hist_field_shape = (self.nframes,1,self.nbins_f)

        # Radial distribution histograms (only if there are NPs).
        # Disregard PBC, because we are interested in short distances. Therefore, generate the histogram
        #   only up to a distance of 10-20 grid points. After phase 7, we can also generate distributions
        #   between all shell beads, but this is handled automatically by bead indexing.
        # Build pair distance lists sequentially for each time frame and throw away the smaller distances to save memory.
        # Note that these histograms are nor normlized, so they need to be normalized later.
        if self.pop > 0:
            self.pds = []
            for c in self.csa[:]:
                self.pds.append([])
                for bi in self.bead_ind:
                    pd = pdist(c[bi].reshape((self.pop*len(bi),2)))
                    pd = pd[pd<self.rmax]
                    self.pds[-1].append(pd)
            self.hist_radial = array([[histogram(pd, bins=self.nbins_rad, range=self.rrange)[0] for pd in pds] for pds in self.pds])
            self.hist_radial_shape = (self.nframes,len(self.bead_ind),self.nbins_rad)

        # Residual total fields and order parameters at bead positions (core and shell).
        if self.pop > 0:
            self.edge = range(64)
            self.splines = [[RectBivariateSpline(self.edge,self.edge,ff,kx=1,ky=1) for ff in f] for f in self.cga]
            self.values_res = array([[[s.ev(cc[:,0],cc[:,1]) for cc in c] for s in self.splines[i]] for i,c in enumerate(self.csa[:])])
            self.totals_res = self.values_res[:,0,:,:] + self.values_res[:,1,:,:]
            self.orders_res = self.values_res[:,0,:,:] - self.values_res[:,1,:,:]
            hist_total = lambda data: histogram(data, bins=self.nbins_res, range=self.totalrange)[0]
            hist_order = lambda data: histogram(data, bins=self.nbins_res, range=self.orderrange)[0]
            self.hist_res_total = array([[hist_total(t[bi].reshape((self.pop*len(bi),))) for bi in self.bead_ind] for t in self.totals_res])
            self.hist_res_order = array([[hist_order(t[bi].reshape((self.pop*len(bi),))) for bi in self.bead_ind] for t in self.orders_res])
            self.hist_res_shape = (self.nframes,len(self.bead_ind),self.nbins_res)

        # Nanoparticle rotation angle distributions, which makes sense after phase 8.
        if self.phase > 8 and self.pop > 0:
            self.hist_angles = [histogram(angles, bins=self.nbins_ang)[0] for angles in self.angles]

        return time.time() - T

    def save(self):
        """Save the analyzed data to archives."""

        T = time.time()

        if hasattr(self, 'instants'):

            # Stack vectors of energies and other instantaneous scalar results. Keep the indices together
            #   with all these, for reference. We might need to drop the last elements in instantaneous results.
            # Note that the saved array is smaller when no NPs are in the system (no offsets).
            if len(self.instants[0]) > len(self.averages[0]):
                n = len(self.averages[0])
                self.indices = self.indices[:-1]
                self.instants = [i[:-1] for i in self.instants]
                if self.pop > 0:
                    self.offsets = [o[:-1] for o in self.offsets]
                    if self.phase > 8:
                        self.offsets_angles = [o[:-1] for o in self.offsets_angles]
            if self.pop > 0:
                if self.phase > 8:
                    energy = vstack([self.indices]+self.instants+self.averages+self.deviations+self.offsets+self.offsets_angles).transpose()
                else:
                    energy = vstack([self.indices]+self.instants+self.averages+self.deviations+self.offsets).transpose()
            else:
                energy = vstack([self.indices]+self.instants+self.averages+self.deviations).transpose()

            # Save all instananeous results in the energy archive.
            save(self.fname+".energy.npy", energy)

        if hasattr(self, 'frames'):

            # Save also the frames used in histograms for reference, in a separate archive.
            save(self.fname+".hist-frames.npy", self.frames)
            
            # Save histograms of field values.
            total = self.hist_field_total.reshape(self.hist_field_shape)
            order = self.hist_field_order.reshape(self.hist_field_shape)
            save_field = concatenate([total,order], axis=1)
            save(self.fname+".hist-field.npy", save_field)

            # Save histograms of radial distributions and residual fields.
            # Note that the residuals for core and shell beads are stacked after phase 7.
            if self.pop > 0:
                save(self.fname+".hist-radial.npy", array(self.hist_radial))
                total = self.hist_res_total.reshape(self.hist_res_shape)
                order = self.hist_res_order.reshape(self.hist_res_shape)
                save_residual = concatenate([total, order], axis=1)
                save(self.fname+".hist-residual.npy", save_residual)

            # Save histograms of angular distributions, which is available after phase 8.
            if self.phase > 8 and self.pop > 0:
                save(self.fname+".hist-ang.npy", array(self.hist_angles))

        return time.time() - T


if __name__ == "__main__":

    # First argument must be the output file.
    fout = sys.argv[1]

    # Load data and initiate analysis object.
    T = time.time()
    a = Analysis(fout)
    print "init: %.2fs" %(time.time()-T)

    # Do energy analysis.
    if "energy" in sys.argv:
        print "energy: %.2fs" %a.analyze_energy()

    # Do histogram analysis.
    if "histograms" in sys.argv:
        print "histograms: %.2fs" %a.analyze_histograms()

    # Save archives.
    if "energy" in sys.argv or "histograms" in sys.argv:
        print "save:, %.2fs" %a.save()
