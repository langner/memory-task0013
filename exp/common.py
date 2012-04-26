import bz2
import os
import sys

import numpy
import pylab
import scipy
from scipy import misc
from scipy import ndimage
from scipy import optimize
from scipy import spatial

import mahotas
import pymorph


def getscalefromname(fn):
    """Get task-specific scale from file name."""

    return int(fn.split('-')[2].split('.')[0])


def getscalebarlength(im, xstart, ystart):
    """Get the length of a scale bar in a an image

    This function will measure a scale bar that is placed on the right
    hand side of the position `start`, assuming the background is black (<100)
    and the scale bar is white (>200) in the first channel."""

    # Start at the given coordinates
    # Their order should be reversed before this function call, if needed
    i = xstart
    j = ystart

    # Increase second index until white encountered
    while im[i,j] < 100:
        j += 1

    # Now increase the measure until white ends
    length = 0
    j += 1
    while im[i,j] > 200:
        length += 1
        j += 1

    return length


def stretchimage(img):
    """Stretch out the contrast of an image"""

    stretched = img - img.min()
    stretched = 255.0 * stretched / float(stretched.ptp())
    stretched = stretched.astype(numpy.uint8)
    return stretched


def balanceimage(img, r, R):
    """Balance the brightness of an image by leveling"""

    img_min = ndimage.minimum_filter(img, r)
    img_uni = ndimage.uniform_filter(img, R)
    return img - numpy.minimum(img_min, img_uni)


def rdfcorrection(X,Y,R):
    """Correction to radial distribution function due to periodic
    boundary conditions, without actually considering them"""

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
    """Radial ditribution function for points in an image"""

    # Number of coordinates and image shape.
    # Now, the RDF correction is effective to no more than about a half of the
    #   smallest dimension, so limit the histogram range to that,
    #   or 256 pixels, whichever is smaller.
    # The width of a bin should be maximum half a pixel.
    N = len(coords)
    X = img.shape[0]
    Y = img.shape[1]
    maxr = min(256,min(X,Y)/2)

    pds = spatial.distance.pdist(coords)
    hist,bins = scipy.histogram(pds, bins=maxr*2, range=[0,maxr])
    dr = bins[1] - bins[0]
    bins = bins[:-1] + dr/2.0
    hist = 1.0 * hist
    hist *= X*Y * 2.0 / (N*(N-1))
    hist /= 2 * scipy.pi * rdfcorrection(X,Y,bins) * bins * dr
    return hist, bins


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

        # Load the image as an array, flattened
        self.fpath = fpath
        self.img = misc.imread(self.fpath, flatten=1)
        self.img = self.img.astype(numpy.uint8)
        self.imgwidth = self.img.shape[1]
        self.imgheight = self.img.shape[0]

        # If no cropping limits chosen, use entire image
        if cropx == None:
            cropx = (0,-1)
        if cropy == None:
            cropy = (0,-1)

        # Cropy image saving the cropped out portion for reading the scale bar
        self.img_scale = self.img[cropy[1]:,:]
        self.img = self.img[cropy[0]:cropy[1],cropx[0]:cropx[1]]

        # Various name components
        self.dirname, self.fname = os.path.split(self.fpath)
        self.name, self.fext = os.path.splitext(self.fname)
        self.dirname = self.dirname.replace("sem/","sem-analyzed/")
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
        self.frdf = "%s-rdf" %self.outname

        # This will be useful, since there are some difference for the benchmark
        self.isbench = self.fpath[:5] == "bench"

        # The scale can usually be extracted from filename
        if scale:
            self.scale = scale
        if not scale and not self.isbench:
            self.scale = int(self.fname.split('-')[2][:-4])

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

    def filterimage(self, gaussian=1, balance=[]):
        """Try to balance brightness and smooth a bit"""

        if os.path.exists(self.ffiltered+".npy.bz2"):
            self.filtered = numpy.load(bz2.BZ2File(self.ffiltered+".npy.bz2"))
        else:
            self.filtered = ndimage.gaussian_filter(self.img, gaussian)
            for b in balance:
                if b:
                    self.filtered = balanceimage(self.filtered, b[0], b[1])
            self.filtered = stretchimage(self.filtered)

    def thresholdimage(self, factor=1.0, erode=True):
        """Thresholding with optional erosion"""

        if not hasattr(self, "filtered"):
            self.filterimage()

        if os.path.exists(self.fthreshold+".npy.bz2"):
            self.threshold = numpy.load(bz2.BZ2File(self.fthreshold+".npy.bz2"))
        else:
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

        # Generate distrivution transform first
        self.distimage()

        if not hasattr(self, "labels"):
            self.labelimage()

        # Load regions if available (takes a bit)
        if os.path.exists(self.fnps+".npy.bz2"):
            self.nps = numpy.load(bz2.BZ2File(self.fnps+".npy.bz2"))
        else:
            self.nps = pymorph.cwatershed(self.dist, self.labels, Bc=numpy.ones((3,3), dtype=bool))

        if os.path.exists(self.fcoms+".npy.bz2"):
            self.coms = numpy.load(bz2.BZ2File(self.fcoms+".npy.bz2"))
        else:
            self.coms = numpy.array(ndimage.center_of_mass(self.filtered, self.nps, range(1,self.nlabels+1)))

        self.regcom = numpy.zeros(self.img.shape, dtype="uint8")
        for icom,com in enumerate(self.coms):
            self.regcom[round(com[0]),round(com[1])] = 1

        self.npcount = len(self.coms)
        self.npconc = 10.0**6*self.npcount/self.area

    def radialdistribution(self):
        """Calculate radial distribution function of regions centers"""

        if not hasattr(self, "coms"):
            self.watershedimage()

        # Load RDFs if available
        if os.path.exists(self.frdf+".npy.bz2"):
            self.bins, self.rdf = numpy.load(bz2.BZ2File(self.frdf+".npy.bz2"))
        else:
            self.rdf, self.bins = radialdistribution(self.coms, self.img)

        # Estimate the bin width
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
        imax = imax*(self.X[imax] < 60) or 20
        P0 = [self.X[imax], 5, 0.5, 2, 5, 100, 2]
        self.popt, self.pcov = optimize.curve_fit(rdfmodel, self.X, self.Y, p0=P0)
        self.rdfmodel = lambda X: rdfmodel(X, *self.popt)
        self.rdfmax = self.popt[0]

    def plot(self, xmax=100, areaunit="micron"):
        """Helper routine for plotting in this task"""

        # Various stages of the image analysis are straightforward
        # 1 - filtered images
        # 2 - threshold with seeds
        # 3 - nanoparticle regions
        # 4 - centers of mass on to of original image
        # The list figure_params contains parameter tuples to pass to imshow,
        #   with the mandatory and optional parameters (in a dict) as elements
        figure_params = [ (self.filtered, dict()),
                          (pymorph.overlay(self.threshold, self.seeds), dict()),
                          ((self.nps % 20 + 5)*(self.nps>0), dict(cmap=pylab.cm.spectral)),
                          (pymorph.overlay(self.img, self.regcom), dict()),
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
        if not os.path.exists(self.frdf+".npy.bz2"):
            numpy.save(self.frdf+".npy", [self.bins, self.rdf])

    def savefigures(self):
        """Helper routine for saving figures from pylab in this task"""

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
        pylab.figure(5)
        pylab.savefig(self.frdf+".png")
