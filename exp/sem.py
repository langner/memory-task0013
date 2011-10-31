import sys

from common import pylab, SEMAnalysis


imagewidth = 645

# Physical image size in nanometers
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


if __name__ == "__main__":

    # Create object and load image
    analysis = SEMAnalysis()
    analysis.loadimage(sys.argv[1], cropy=(0,430))

    # Do the analysis
    analysis.filterimage()
    analysis.thresholdimage()
    analysis.labelimage()
    analysis.distimage()
    analysis.watershedimage()
    analysis.radialdistribution()

    # Do the plotting
    ps = pixelsize(int(analysis.name.split('-')[-1]))
    analysis.plot(rdfpixelsize=ps, rdfxunits="nm")

    if "save" in sys.argv:
        analysis.savearchives()
        analysis.savefigures()
    else:
        pylab.show()
