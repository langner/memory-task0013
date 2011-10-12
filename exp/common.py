import bz2
import os
import sys

import numpy
import pylab
import scipy
from scipy import misc
from scipy import ndimage
from scipy import spatial

import mahotas
import pymorph


def stretchimage(img):
    """ Stretch out the contrast of an image. """

    stretched = img - img.min()
    stretched = 255.0 * stretched / float(stretched.ptp())
    stretched = stretched.astype(numpy.uint8)
    return stretched

def balanceimage(img, r, R):
    """ Balance the brightness of an image by leveling. """

    img_min = ndimage.minimum_filter(img, r)
    img_uni = ndimage.uniform_filter(img, R)
    return img - numpy.minimum(img_min, img_uni)

def rdfcorrection(X,Y,R):
    """ Correction to radial distribution function due to periodic
    boundary conditions, without actually considering them. """

    ang1 = 2 * scipy.pi
    ang2 = 2 * scipy.pi - 2.0
    ang3 = 2 * scipy.pi - 2.0
    ang4 = 3 * scipy.pi / 2.0 - 2.0
    V1 = (0.5*X - R)*(0.5*Y - R)
    V2 = R*(0.5*X - R)
    V3 = R*(0.5*Y - R)
    V4 = R**2
    ref = (V1+V2+V3+V4) * 2 * scipy.pi
    return (ang1*V1 + ang2*V2 + ang3*V3 + ang4*V4) / ref

def radialdistribution(coords, img):
    """ Radial ditribution function for points in an image."""

    N = len(coords)
    X = img.shape[0]
    Y = img.shape[1]

    pds = spatial.distance.pdist(coords)
    hist,bins = scipy.histogram(pds, bins=512, range=[0,256])
    dr = bins[1] - bins[0]
    bins = bins[:-1] + dr/2.0
    hist = 1.0 * hist
    hist *= X*Y * 2.0 / (N*(N-1))
    hist /= 2 * scipy.pi * rdfcorrection(X,Y,bins) * bins * dr
    return hist, bins


class SEMAnalysis():

    def __init__(self):

        pass
    
    def loadimage(self, fpath, cropx=(0,-1), cropy=(0,-1)):
        """ Load image with cropping and set various filename serivatives. """

        self.fpath = fpath

        self.img = misc.imread(self.fpath, flatten=1)
        self.img = self.img[cropy[0]:cropy[1],cropx[0]:cropx[1]]
        self.img = self.img.astype(numpy.uint8)

        self.dirname, self.fname = os.path.split(self.fpath)
        self.name, self.fext = os.path.splitext(self.fname)
        self.dirname = self.dirname.replace("sem/","sem-analyzed/")
        self.outname = "%s/%s" %(self.dirname,self.name)
        if not os.path.exists(self.dirname):
            os.makedirs(self.dirname)

        self.ffiltered = "%s-filtered" %self.outname
        self.fthreshold = "%s-threshold" %self.outname
        self.fseeds = "%s-seeds" %self.outname
        self.fdist = "%s-dist" %self.outname
        self.fnps = "%s-nps" %self.outname
        self.fcoms = "%s-coms" %self.outname
        self.frdf = "%s-rdf" %self.outname

    def filterimage(self, gaussian1=1, gaussian2=1, balance=[(5,75),(3,25)]):
        """ Try to balance brightness and smooth a bit. """

        self.filtered = ndimage.gaussian_filter(self.img, gaussian1)
        for b in balance:
            self.filtered = balanceimage(self.filtered, b[0], b[1])
        self.filtered = ndimage.gaussian_filter(self.filtered, gaussian2)
        self.filtered = stretchimage(self.filtered)

    def thresholdimage(self, factor=1.1, erode=True):
        """ Thresholding with optional erosion. """

        if not hasattr(self, "filtered"):
            self.filterimage()

        self.threshold = self.filtered > factor*mahotas.otsu(self.filtered)
        if erode:
            self.threshold = pymorph.erode(self.threshold)

    def labelimage(self, threshold=None):
        """ Label regions in image corresponding to regional maxima. """

        if threshold is None:
            threshold = self.threshold
        self.seeds = pymorph.regmax(self.filtered) * threshold
        self.labels, self.nlabels = ndimage.label(self.seeds, structure=numpy.ones((3,3)))

    def distimage(self):
        """ Distance transform of the Boolean threashold. """

        if not hasattr(self, "threshold"):
            self.thresholdimage()

        self.dist = ndimage.distance_transform_edt(self.threshold)
        self.dist = self.dist.max() - self.dist
        self.dist = stretchimage(self.dist)

    def watershedimage(self):

        # Generate distrivution transform first.
        self.distimage()

        if not hasattr(self, "labels"):
            self.labelimage()

        # Load regions if available (takes a bit).
        if os.path.exists(self.fnps+".npy.bz2"):
            self.nps = numpy.load(bz2.BZ2File(self.fnps+".npy.bz2"))
        else:
            self.nps = pymorph.cwatershed(self.dist, self.labels, Bc=numpy.ones((3,3), dtype=bool))

        self.coms = numpy.array(ndimage.center_of_mass(self.filtered, self.nps, range(1,self.nlabels+1)))
        self.regcom = numpy.zeros(self.img.shape, dtype="uint8")
        for icom,com in enumerate(self.coms):
            self.regcom[round(com[0]),round(com[1])] = 1

    def radialdistribution(self):
        
        if not hasattr(self, "coms"):
            self.watershedimage()
        self.rdf, self.bins = radialdistribution(self.coms, self.img)

    def plot(self):
        """ Helper routine for plotting in this task. """

        pylab.figure(1)
        pylab.imshow(self.filtered)

        pylab.figure(2)
        pylab.imshow(pymorph.overlay(self.threshold, self.seeds))

        pylab.figure(3)
        pylab.imshow(self.nps % 55)

        pylab.figure(4)
        pylab.imshow(pymorph.overlay(self.img, self.regcom))

        pylab.figure(5)
        pylab.plot(self.bins, self.rdf)
        pylab.ylim([0.0, max(3.0,max(self.rdf))])
        pylab.xlim([0, 50])
        pylab.xlabel("distance [pixels]")
        pylab.grid()
        a = pylab.axes([0.60, 0.60, 0.25, 0.25], axisbg='y')
        pylab.setp(a, xlim=(0,.2), xticks=[50,100,200], yticks=[1.0])
        pylab.xlim([0,256])
        pylab.plot(self.bins, self.rdf)
        pylab.grid()

    def savearchives(self):
        """ Helper routine for saving NumPy archives in this task. """

        numpy.save(self.ffiltered+".npy", self.filtered)
        numpy.save(self.fthreshold+".npy", self.threshold)
        numpy.save(self.fnps+".npy", self.nps)
        numpy.save(self.fcoms+".npy", self.coms)
        numpy.save(self.frdf+".npy", [self.bins, self.rdf])

    def savefigures(self):
        """ Helper routine for saving figures from pylab in this task. """

        pylab.figure(1)
        pylab.savefig(self.ffiltered+".png")
        pylab.figure(2)
        pylab.savefig(self.fthreshold+".png")
        pylab.figure(3)
        pylab.savefig(self.fnps+".png")
        pylab.figure(4)
        pylab.savefig(self.fcoms+".png")
        pylab.figure(5)
        pylab.savefig(self.frdf+".png")
