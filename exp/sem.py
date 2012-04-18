import sys

from common import pylab, SEMAnalysis


# The size of images in pixels
# Note the the height is already the cropped value
imagewidth = 645
imageheight = 430

# Physical image size in nanometers (calculated from one image)
imagesize = 10**6 * 124.03846153846154

# Magnifications (for image scales)
magnifications = {
    100: 200 * 10**3,
    200: 100 * 10**3,
    500: 50 * 10**3,
    1000: 25 * 10**3,
    2000: 12.5 * 10**3,
    5000: 3.125 * 10**3
}

# Pixel size calculation
pixelsize = lambda s: imagesize / magnifications[s] / imagewidth

# Balance filters for images
balance_default = {
    500 : [(2,25),(1,5)],
}
balance_custom = {
    "sem/5k/C3/0-Image5-1000.jpg":  [(1,3), (2,25)],
    "sem/5k/C3/2-Image7-1000.jpg":  [(1,3), (2,25)],
    "sem/5k/C8/0-Image27-500.jpg":  [(5,75), (2,25)],
}

# Threshold values for images
threshold_default = {
    500 : 0.75
}
threshold_custom = {
    "sem/5k/C1/2-Image2-200.jpg":           1.0,
    "sem/5k/C3/0-Image5-1000.jpg":          0.9,
    "sem/5k/C3/2-Image7-1000.jpg":          0.8,
    "sem/20k/C1/14-Image22-500.jpg":        0.9,
    "sem/20k/C5/14-Image33-500.jpg":        1.0,
}


if __name__ == "__main__":

    # Create object and load image
    fn = sys.argv[1]
    analysis = SEMAnalysis()
    analysis.loadimage(fn, cropy=(0,imageheight))

    # Do the analysis
    # In hand picked cases, apply custom filter parameters
    balance = balance_custom.get(fn) or balance_default[analysis.magnification]
    threshold = threshold_custom.get(fn) or threshold_default[analysis.magnification]
    analysis.filterimage(balance=balance)
    analysis.thresholdimage(factor=threshold)
    analysis.labelimage()
    analysis.distimage()
    analysis.watershedimage()
    analysis.radialdistribution()
    analysis.fitrdf()

    # Print some data
    N = len(analysis.coms)
    D = 10.0**6*N/(imageheight*imagewidth*analysis.ps**2)
    print "**************************"
    print "Number of particles: %i" %N
    print "Particle concentration: %.1f/micron^2" %D
    print "First peak: %.1fnm" %analysis.popt[0]
    print "**************************"

    # Do the plotting
    ps = pixelsize(int(analysis.name.split('-')[-1]))
    analysis.plot(rdfpixelsize=ps, rdfxunits="nm")

    if "save" in sys.argv:
        analysis.savearchives()
        analysis.savefigures()
    else:
        pylab.show()
