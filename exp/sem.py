import sys

from common import pylab, SEMAnalysis


if __name__ == "__main__":

    # Create object and load image.
    analysis = SEMAnalysis()
    analysis.loadimage(sys.argv[1], cropy=(0,430))

    # Do the analysis.
    analysis.filterimage()
    analysis.thresholdimage()
    analysis.labelimage()
    analysis.distimage()
    analysis.watershedimage()
    analysis.radialdistribution()

    # Do the plotting.
    analysis.plot()

    if "save" in sys.argv:
        analysis.savearchives()
        analysis.savefigures()
    else:
        pylab.show()
