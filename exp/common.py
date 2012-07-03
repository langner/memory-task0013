"""
These are quite general tools for to analyzing these types of SEM images.
In particular, they are used here for both benchmarking and real analysis.
"""


import bz2
import os
import random
import sys

import numpy
import scipy

from scipy import misc
from scipy import ndimage
from scipy import optimize
from scipy import spatial

import mahotas
import pymorph


def getscalebarlength(im, xstart, ystart):
    """Get the length of a scale bar in an image.

    This function will measure a scale bar that is placed on the right
    hand side of the position (xstart,ystart), assuming the background is
    close to black (<100) and that the scale bar is close to white (>200)
    at laest in the first channel."""

    # Start at the given coordinates. Their order should be reversed before
    #   this function call, if needed.
    i = xstart
    j = ystart

    # Increase second index until something not dark is encountered.
    while im[i,j] < 100:
        j += 1

    # Now increase the measure until something quite bright is encountered.
    length = 0
    while im[i,j] > 200:
        length += 1
        j += 1

    return length


def stretchimage(img):
    """Stretch out the contrast of an image."""

    stretched = img - img.min()
    stretched = 255.0 * stretched / float(stretched.ptp())
    return stretched.astype(numpy.uint8)


def balanceimage(img, r, R):
    """Balance the brightness of an image by leveling.

    This is achieved here by applying a minimum filter over radius r
    and a uniform filter over radius R, and substracting the minimum
    of the two from the original image."""

    img_min = ndimage.minimum_filter(img, r)
    img_uni = ndimage.uniform_filter(img, R)
    return img - numpy.minimum(img_min, img_uni)


def rdfcorrection(X, Y, R):
    """Correction to radial distribution function due to periodic
    boundary conditions, without actually considering them."""

    pXY = numpy.pi * X * Y
    return 1.0 - 2*R*(X+Y)/pXY + R*R*(4.0-numpy.pi)/pXY


def normalize_rdf(hist, X, Y, bins, N, npairs=None):
    """Normalize a radial distribution function."""
    
    # Area of strips at all distances (assume constant spacing).
    dx = bins[1] - bins[0]
    strips = 2 * numpy.pi * bins * rdfcorrection(X, Y, bins) * dx

    # The radial scaling factor.
    npairs = npairs or N*(N-1)/2
    factor = X*Y / (strips * npairs)

    return factor * hist


def clusterparticles(coms, threshold=1.0):
    """Cluster the coordinates of nanoparticles in an image.

    Use the DBSCAN algorithm with minimum sample count of 1."""

    # Import this just now, because it seems to drag in dependecy in matplotlib.
    from sklearn import cluster

    # Setup the analysis object and do the clustering
    dbscan = cluster.DBSCAN(eps=threshold, min_samples=1)
    dbscan.fit(coms)

    # Return the labels, but incremented by one.
    return dbscan.labels_ + 1


def radialdistribution(coords, img, ps):
    """Radial ditribution function for points in an image"""

    # Number of coordinates and image shape.
    # Now, the RDF correction is effective to no more than about a half of the
    #   smallest dimension, so limit the histogram range to that,
    #   or 256 pixels, whichever is smaller.
    # The width of a bin should always be at most 1nm, so the number of bin will vary.
    N = len(coords)
    X = img.shape[0]
    Y = img.shape[1]
    maxr = min(256,min(X,Y)/2)

    pds = spatial.distance.pdist(coords)
    nbins = max(maxr,int(maxr*ps))
    hist,bins = scipy.histogram(pds, bins=nbins, range=[0,maxr])
    bins = bins[:-1] + (bins[1]-bins[0])/2.0
    return normalize_rdf(hist, X, Y, bins, N), bins


def rdfmodel(X, p, m, l, a, b, t, gd):
    """Approximate Radial distribution function

    Coded according to Matteoli and Mansoori, JPC 1995
    All parameters have same meaning, but the first peak position, p, is here additionally
    """
    Y = X/p
    Ym1 = Y-1
    left = gd*numpy.exp(-t*Ym1*Ym1)
    right = 1.0 + numpy.power(Y,-m)*(gd-1.0-l) + numpy.exp(-a*Ym1)*numpy.cos(b*Ym1)*(Ym1+l)/Y
    return (Y<1)*left + (Y>=1)*right


class SEMAnalysis():

    def __init__(self, fpath=None, cropx=None, cropy=None, scalebarstart=None):

        if fpath != None:
            self.loadimage(fpath, cropx=cropx, cropy=cropy, scalebarstart=scalebarstart)
    
    def loadimage(self, fpath, cropx=None, cropy=None, scale=None, scalebarstart=None):
        """Load image with cropping, set filename derivatives and extract scale if possible"""

        # Load the image as an array, flattened.
        self.fpath = fpath
        self.img = misc.imread(self.fpath, flatten=1)
        self.img = self.img.astype(numpy.uint8)

        # If no cropping limits chosen, use entire image.
        if cropx == None:
            cropx = (0,-1)
        if cropy == None:
            cropy = (0,-1)

        # Cropy image saving the cropped out portion for reading the scale bar.
        self.img_scale = self.img[cropy[1]:,:]
        self.img = self.img[cropy[0]:cropy[1],cropx[0]:cropx[1]]

        # Only now can we take the image dimensions.
        self.imgwidth = self.img.shape[1]
        self.imgheight = self.img.shape[0]

        # Various name components.
        # Replace both 'sem' and 'homo', but these should exclude each other.
        self.dirname, self.fname = os.path.split(self.fpath)
        self.name, self.fext = os.path.splitext(self.fname)
        self.dirname = self.dirname.replace("sem/","sem-analyzed/")
        self.dirname = self.dirname.replace("homo/","homo-analyzed/")
        self.outname = "%s/%s" %(self.dirname,self.name)
        if not os.path.exists(self.dirname):
            os.makedirs(self.dirname)

        # Various related filenames
        self.ffiltered = "%s-filtered" %self.outname
        self.fthreshold = "%s-threshold" %self.outname
        self.fseeds = "%s-seeds" %self.outname
        self.fdist = "%s-dist" %self.outname
        self.fnps = "%s-nps" %self.outname
        self.fcoms = "%s-coms" %self.outname
        self.fclusters = "%s-clusters" %self.outname
        self.frdf = "%s-rdf" %self.outname

        # This will be useful, since there are some difference for the benchmark
        self.isbench = self.fpath[:5] == "bench"

        # The scale can usually be extracted from filename
        if scale:
            self.scale = scale
        if not scale and not self.isbench:
            self.scale = int(self.fname.split('-')[2][:-4])

        # Extract the ligand type from the file name (not for benchmarks).
        if not self.isbench:
            self.ligand_type = self.dirname.split("/")[1]

        # Extract the length of the scale bar, if starting position provided
        # Assume that the bar is in the portion cropped away
        if scalebarstart:
            xstart = scalebarstart[1] - cropy[1] + (cropy[1] == -1) - 1
            ystart = scalebarstart[0] - cropx[1] + (cropx[1] == -1) - 1
            self.scalebarlength = getscalebarlength(self.img_scale, xstart, ystart)

        # If both above are available, derive the pixel size, and area
        if hasattr(self, "scale") and hasattr(self, "scalebarlength"):
            self.pixelsize = 1.0 * self.scale / self.scalebarlength
            self.ps = self.pixelsize
            self.area = self.imgwidth*self.imgheight*self.ps**2

        # A time can be extracted from the filename (not for benchmarks)
        if not self.isbench:
            self.time = int(self.fname.split('-')[0])

    def load_archive(self, archive):
        fname = "f"+archive
        fpath = getattr(self, fname) + ".npy.bz2"
        if os.path.exists(fpath):
            data = numpy.load(bz2.BZ2File(fpath))
            if archive == "rdf":
                setattr(self, "bins", data[0])
                setattr(self, "rdf", data[1])
            else:
                setattr(self, archive, data)
            return True
        else:
            return False

    def filterimage(self, gaussian=1, balance=[]):
        """Try to balance brightness and smooth a bit"""

        if self.load_archive("filtered"):
            return

        self.filtered = ndimage.gaussian_filter(self.img, gaussian)
        for b in balance:
            if b:
                self.filtered = balanceimage(self.filtered, b[0], b[1])
        self.filtered = stretchimage(self.filtered)

    def thresholdimage(self, factor=1.0, erode=True):
        """Thresholding with optional erosion"""

        if not hasattr(self, "filtered"):
            self.filterimage()

        if self.load_archive("threshold"):
            return

        limit = factor*mahotas.otsu(self.filtered)
        self.threshold = self.filtered > limit
        if erode:
            self.threshold = pymorph.erode(self.threshold)

    def labelimage(self, threshold=None):
        """Label regions in image corresponding to regional maxima"""

        if threshold is None:
            threshold = self.threshold

        self.seeds = pymorph.regmax(self.filtered) * threshold
        self.labels, self.nlabels = ndimage.label(self.seeds, structure=numpy.ones((3,3)))

    def distimage(self):
        """Distance transform of the Boolean threashold"""

        if not hasattr(self, "threshold"):
            self.thresholdimage()

        self.dist = ndimage.distance_transform_edt(self.threshold)
        self.dist = self.dist.max() - self.dist
        self.dist = stretchimage(self.dist)

    def watershedimage(self):
        """Watershed and extract region centers of mass"""

        # Generate distribution transform first.
        self.distimage()

        if not hasattr(self, "labels"):
            self.labelimage()

        # Load regions if available (takes a bit).
        if not self.load_archive("nps"):
            self.nps = pymorph.cwatershed(self.dist, self.labels, Bc=numpy.ones((3,3), dtype=bool))

        if not self.load_archive("coms"):
            self.coms = numpy.array(ndimage.center_of_mass(self.filtered, self.nps, range(1,self.nlabels+1)))

        self.regcom = numpy.zeros(self.img.shape, dtype="uint8")
        for icom,com in enumerate(self.coms):
            self.regcom[round(com[0]),round(com[1])] = 1

        self.npcount = len(self.coms)
        self.npconc = 10.0**6*self.npcount/self.area

    def clusterparticles(self):

        # Use a threshold of 35nm for 5k and 45nm for 20k ligands.
        # Benchmarks should use their own thresholds (we can try some default, though).
        if not self.isbench:
            self.cluster_threshold = (32.0 + (self.ligand_type == "20k")*10)/self.ps
        else:
            self.cluster_threshold = 10.0

        if not hasattr(self, "coms"):
            self.watershedimage()

        # Load the cluster data if possible.
        if self.load_archive("clusters"):

            # We still need to restore the cluster labels.
            self.cluster_labels =  [self.clusters[com[0],com[1]] for com in self.coms]

        else:

            # Do the clustering in a separate function and save the labels.
            self.cluster_labels = clusterparticles(self.coms, self.cluster_threshold)

            # Now use the cluster labels to color nanoparticle clusters.
            self.clusters = numpy.copy(self.nps)
            for ci,cn in enumerate(self.cluster_labels):
                self.clusters[self.nps == ci+1] = cn

        # Some statistics.
        self.cluster_count = len(set(self.cluster_labels))
        self.cluster_sizes = numpy.array([numpy.sum(self.cluster_labels==l) for l in set(self.cluster_labels)])
        self.cluster_size = 1.0 * self.npcount / self.cluster_count
        self.cluster_size_inv = 1.0 / self.cluster_size
        self.clusters_per_micron = 10.0**6 * self.cluster_count / self.area
        self.clusters_per_micron_inv = 1.0 / self.clusters_per_micron

    def radialdistribution(self):
        """Calculate radial distribution function of regions centers"""

        if not hasattr(self, "coms"):
            self.watershedimage()

        # Load RDFs if available.
        if not self.load_archive("rdf"):
            self.rdf, self.bins = radialdistribution(self.coms, self.img, self.ps)

        # Estimate the bin width.
        self.dbin = self.ps*(self.bins[1] - self.bins[0])

    def fitrdf(self):
        """Fit a simple equation to the radial distribution function"""

        if not hasattr(self, "rdf"):
            self.radialdistribution()

        self.X = self.bins*self.ps
        self.Y = self.rdf

        # We need to discard the zeros before the peak
        istart = numpy.where(self.Y>0)[0][0]
        self.X = self.X[istart:]
        self.Y = self.Y[istart:]

        # It is also nice to have an apprixmate peak position to start from,
        #   but just start from 20 nanometers if the first peak is not highest.
        imax = numpy.argmax(self.Y)
        imax = imax*(self.X[imax] < 100) or 50

        # Do the actual curve fitting now.
        P0 = [self.X[imax], 5, 0.5, 2, 5, 100, 2]
        self.popt, self.pcov = optimize.curve_fit(rdfmodel, self.X, self.Y, p0=P0)

        # An RDF model equation based on this fit.
        self.rdfmodel = lambda X: rdfmodel(X, *self.popt)

        # Some handy aliases for analysis and plotting.
        self.rdfmax = self.popt[0]
        self.rdfheight = self.popt[-1]
        self.rdfmax_var = self.pcov[0,0]
        self.rdfheight_var = self.pcov[1,1]

    def plot(self, pylab, xmax=100, areaunit="micron"):
        """Helper routine for plotting in this task"""

        if not hasattr(self, "clusters"):
            self.clusterparticles()
        if not hasattr(self, "rdf"):
            self.radialdistribution()

        # Various stages of the image analysis are straightforward
        # 1 - filtered images
        # 2 - threshold with seeds
        # 3 - nanoparticle regions
        # 4 - centers of mass on to of original image
        # 5 - colored clusters of nanoparticles
        # The list figure_params contains parameter tuples to pass to imshow,
        #   with the mandatory and optional parameters (in a dict) as elements
        figure_params = [ (self.filtered, dict()),
                          (pymorph.overlay(self.threshold, self.seeds), dict()),
                          ((self.nps % 20 + 5)*(self.nps>0), dict(cmap=pylab.cm.spectral)),
                          (pymorph.overlay(self.img, self.regcom), dict()),
                          ((self.clusters % 20 + 5)*(self.clusters>0), dict(cmap=pylab.cm.spectral)),
                        ]
        for ip,params in enumerate(figure_params):
            pylab.figure(ip+1)
            pylab.imshow(params[0], **params[1])
            pylab.xticks([])
            pylab.yticks([])

        # Now the RDF plot (zoomed in to the first few peaks)
        # Remember to scale everything by the pixel size
        # Add lines and labels at multiples of the first peak
        # Also add a larger portion of the RDF in an inset
        pylab.figure(len(figure_params)+1)
        ymax = 1.01*max(2.5,max(self.rdf))
        rdfxunits = "nm"*hasattr(self,"ps") or "pixels"
        expconc = self.npconc*(areaunit=="micron") or self.npconc / 10.0**6
        explabel = "experimental %.1f/micron$^2$"*(areaunit=="micron") or "experiment (%.1f/nm$^2$)"
        pylab.plot(self.bins*self.ps, self.rdf, label=explabel %expconc)
        if hasattr(self, "rdfmodel"):
            pylab.plot(self.bins*self.ps, self.rdfmodel(self.bins*self.ps), label="hard sphere liquid model")
            for peak,label in [(1,"$d=%.1f\mathrm{nm}$" %self.rdfmax), (numpy.sqrt(3),"$\sqrt{3}$"), (2,"$2$"), (numpy.sqrt(7),"$\sqrt{7}$")]:
                x = peak*self.rdfmax
                pylab.axvline(x=x, ymin=0, ymax=ymax, linestyle='--', color='gray')
                pylab.text(x, ymax*1.02, label, fontsize=15, horizontalalignment="center")
        pylab.xlabel("distance [%s]" %rdfxunits)
        pylab.ylabel("g(r)")
        pylab.xlim([0, xmax])
        pylab.ylim([0.0, ymax])
        pylab.legend()
        pylab.grid()
        a = pylab.axes([0.60, 0.50, 0.25, 0.25], axisbg='y')
        pylab.plot(self.bins*self.ps, self.rdf)
        pylab.setp(a, xticks=map(int,[xmax/2.0,xmax,1.5*xmax]), yticks=[1.0])
        pylab.xlim([0,2*xmax])
        pylab.grid()

    def savearchives(self):
        """Helper routine for saving NumPy archives in this task"""

        if not os.path.exists(self.ffiltered+".npy.bz2"):
            numpy.save(self.ffiltered+".npy", self.filtered)
        if not os.path.exists(self.fthreshold+".npy.bz2"):
            numpy.save(self.fthreshold+".npy", self.threshold)
        if not os.path.exists(self.fnps+".npy.bz2"):
            numpy.save(self.fnps+".npy", self.nps)
        if not os.path.exists(self.fcoms+".npy.bz2"):
            numpy.save(self.fcoms+".npy", self.coms)
        if not os.path.exists(self.fclusters+".npy.bz2"):
            numpy.save(self.fclusters+".npy", self.clusters)
        if not os.path.exists(self.frdf+".npy.bz2"):
            numpy.save(self.frdf+".npy", [self.bins, self.rdf])

    def savefigures(self, pylab):
        """Helper routine for saving figures from pylab in this task"""

        if pylab.get_fignums() == []:
            self.plot(pylab)

        if not os.path.exists(self.ffiltered+".png"):
            pylab.figure(1)
            pylab.savefig(self.ffiltered+".png", bbox_inches="tight")
        if not os.path.exists(self.fthreshold+".png"):
            pylab.figure(2)
            pylab.savefig(self.fthreshold+".png", bbox_inches="tight")
        if not os.path.exists(self.fnps+".png"):
            pylab.figure(3)
            pylab.savefig(self.fnps+".png", bbox_inches="tight")
        if not os.path.exists(self.fcoms+".png"):
            pylab.figure(4)
            pylab.savefig(self.fcoms+".png", bbox_inches="tight")
        if not os.path.exists(self.fclusters+".png"):
            pylab.figure(5)
            pylab.savefig(self.fclusters+".png", bbox_inches="tight")
        pylab.figure(6)
        pylab.savefig(self.frdf+".png")
