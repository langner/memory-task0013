import sys

from common import SEMAnalysis, getscalefromname


# The cropping height should be identical in all images
cropy = 430

# A position just before the scale bar
# This needs to be the same for all images,
#   but appears to change with different magnifications
scalebarstart = {
    200 : (244,441),
    500 : (238,441),
   1000 : (236,441),
}

# Gaussian filter radii
gaussian_default = {
    200 : 2,
    500 : 1,
   1000 : 0,
}
gaussian_custom = {
    "sem/5k/C6/48-Image40-200.jpg" : 3,
    "sem/20k/C3/0-Image22-500.jpg" : 2,
    "sem/20k/C3/0-Image23-500.jpg" : 2,
    "sem/20k/C3/0-Image38-500.jpg" : 2,
    "sem/10k/C3/48-Image7-200.jpg" : 3,
    "sem/20k/C5/0-Image4-200.jpg" : 3,
    "sem/20k/C5/0-Image6-1000.jpg" : 1,
    "sem/20k/C5/0-Image24-500.jpg" : 2,
    "sem/20k/C5/0-Image25-1000.jpg" : 1,
}

# Balance filters for images
# The default for each resolution should be OK,
#   with custom ones defined for chosen images
balance_default = {
    200 : [(25,50), (1,10)],
    500 : [(13,25), (1,5)],
   1000 : [(10,13), (1,2)],
}
balance_custom = {
    "sem/20k/C1/14-Image20-1000.jpg" : [(10,20), (1,3)],
    "sem/20k/C3/0-Image22-500.jpg": [(20,20), (1,13)],
    "sem/20k/C3/0-Image23-500.jpg": [(20,20), (1,13)],
    "sem/20k/C3/0-Image38-500.jpg": [(20,20), (1,13)],
    "sem/20k/C3/0-Image39-1000.jpg" : [(10,10), (1,5)],
    "sem/20k/C5/48-Image24-500.jpg": [(20,20), (1,13)],
}

# Threshold values for images
# Again, the default should yield acceptable nanoparticles
#   in most cases, but this can be tuned as needed
threshold_default = {
    200 : 0.9,
    500 : 0.8,
   1000 : 0.7,
}
threshold_custom = {
    "sem/5k/C3/14-Image3-200.jpg" : 0.7,
    "sem/5k/C3/14-Image6-500.jpg" : 0.9,
    "sem/5k/C6/0-Image20-1000.jpg": 0.9,
    "sem/5k/C6/0-Image23-1000.jpg": 0.9,
    "sem/5k/C6/2-Image3-1000.jpg" : 1.5,
    "sem/5k/C6/48-Image40-200.jpg" : 1.4,
    "sem/5k/C8/48-Image18-200.jpg" : 1.0,
    "sem/5k/C8/48-Image20-1000.jpg" : 1.1,
    "sem/20k/C1/0-Image37-1000.jpg" : 1.4,
    "sem/20k/C1/14-Image20-1000.jpg" : 1.1,
    "sem/20k/C1/14-Image23-1000.jpg" : 1.1,
    "sem/20k/C3/0-Image22-500.jpg" : 1.6,
    "sem/20k/C3/0-Image23-500.jpg" : 1.6,
    "sem/20k/C3/0-Image38-500.jpg" : 1.8,
    "sem/20k/C3/0-Image39-1000.jpg" : 1.5,
    "sem/20k/C3/48-Image6-1000.jpg" : 1.0,
    "sem/10k/C3/48-Image7-200.jpg" : 1.6,
    "sem/20k/C5/0-Image4-200.jpg" : 1.8,
    "sem/20k/C5/0-Image6-1000.jpg" : 1.6,
    "sem/20k/C5/0-Image25-1000.jpg" : 1.0,
    "sem/20k/C5/2-Image18-1000.jpg" : 0.8,
    "sem/20k/C5/14-Image32-1000.jpg" : 1.0,
    "sem/20k/C5/48-Image18-200.jpg" : 2.5,
    "sem/20k/C5/48-Image23-500.jpg" : 1.2,
    "sem/20k/C5/48-Image24-500.jpg" : 1.5,
    "sem/20k/C5/48-Image25-200.jpg" : 2.5,
    "sem/20k/C5/48-Image27-500.jpg" : 2.5,
}

# Erosion on the thresholded image
# Normally we don't want to erode, but not if the resolution is too low
erode_default = {
    200 : False,
    500 : False,
   1000 : False,
}
erode_custom = {
}


if __name__ == "__main__":

    # Import pylab appropriately.
    import matplotlib
    if "save" in sys.argv:
        matplotlib.use("Agg")
    import pylab

    # Create object and load image.
    fn = sys.argv[1]
    sbstart = scalebarstart[getscalefromname(fn)]
    analysis = SEMAnalysis(fn, cropy=(0,430), scalebarstart=sbstart)

    # Do the analysis.
    gaussian = gaussian_custom.get(fn) or gaussian_default[analysis.scale]
    balance = balance_custom.get(fn) or balance_default[analysis.scale]
    threshold = threshold_custom.get(fn) or threshold_default[analysis.scale]
    erode = erode_custom.get(fn) or erode_default[analysis.scale]
    analysis.filterimage(gaussian=gaussian, balance=balance)
    analysis.thresholdimage(factor=threshold, erode=erode)
    analysis.labelimage()
    analysis.distimage()
    analysis.watershedimage()
    analysis.radialdistribution()
    try:
        analysis.fitrdf()
    except:
        "Error fitting to RDF"

    # Print some data.
    print "**************************"
    print "Number of particles: %i" %analysis.npcount
    print "Particle concentration: %.1f/micron^2" %analysis.npconc
    print "First peak: %.1fnm" %analysis.popt[0]
    print "Pixel size: %.2fnm" %analysis.ps
    print "RDF bin size: %.2fnm" %analysis.dbin
    print "**************************"

    # Do the plotting.
    analysis.plot(pylab, xmax=150)

    # Save or show.
    if "save" in sys.argv:
        analysis.savearchives()
        analysis.savefigures(pylab)
    else:
        pylab.show()
