import sys

from common import pylab, SEMAnalysis


# The cropping height should be identical in all images
cropy = 430

# A position just before the scale bar
# This needs to be the same for all images
scalebarstart = (238,441)

# Balance filters for images
# The default for each resolution should be OK,
#   with custom ones defined for chosen images
balance_default = {
    200 : [(2,25), (1,5)],
    500 : [(2,25), (1,5)],
}
balance_custom = {
    "sem/5k/C3/0-Image5-1000.jpg":  [(1,3), (2,25)],
    "sem/5k/C3/2-Image7-1000.jpg":  [(1,3), (2,25)],
    "sem/5k/C8/0-Image27-500.jpg":  [(5,75), (2,25)],
}

# Threshold values for images
# Again, the default should yield acceptable nanoparticles
#   in most cases, but this can be tuned as needed
threshold_default = {
    200 : 0.75,
    500 : 0.75,
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
    analysis.loadimage(fn, cropy=(0,430), scalebarstart=scalebarstart)

    # Do the analysis
    balance = balance_custom.get(fn) or balance_default[analysis.scale]
    threshold = threshold_custom.get(fn) or threshold_default[analysis.scale]
    analysis.filterimage(balance=balance)
    analysis.thresholdimage(factor=threshold)
    analysis.labelimage()
    analysis.distimage()
    analysis.watershedimage()
    analysis.radialdistribution()
    analysis.fitrdf()

    # Print some data
    N = len(analysis.coms)
    D = 10.0**6*N/analysis.area
    print "**************************"
    print "Number of particles: %i" %N
    print "Particle concentration: %.1f/micron^2" %D
    print "First peak: %.1fnm" %analysis.popt[0]
    print "Pixel size: %.2fnm" %analysis.ps
    print "**************************"

    # Do the plotting
    analysis.plot(xmax=150)

    # Save or show
    if "save" in sys.argv:
        analysis.savearchives()
        analysis.savefigures()
    else:
        pylab.show()
