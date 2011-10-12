import sys

from common import pylab, SEMAnalysis


if __name__ == "__main__":

    # Create object and load image.
    analysis = SEMAnalysis()
    analysis.loadimage(sys.argv[1], cropy=(0,850))

    # Do the analysis.
    analysis.radialdistribution()

    # Do the plotting.
    analysis.plot()

    if "save" in sys.argv:
        analysis.savearchives()
        analysis.savefigures()
    else:
        pylab.show()
