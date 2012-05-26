import sys

from common import SEMAnalysis


# There are more of these than for the batch SEM jobs,
#   because each benchmark is from a different source.
benchmarks_size = {
    "benchmark/F20-atom_res.preview.png": (980,980),
    "benchmark/a06fig02a.png": (288,267),
    "benchmark/r6.png": (288,288),
    "benchmark/r7.png": (371,372),
}
benchmarks_cropy = {
    "benchmark/F20-atom_res.preview.png": (0,850),
    "benchmark/a06fig02a.png": (0,228),
    "benchmark/r6.png": (0,270),
    "benchmark/r7.png": (0,350),
}
benchmarks_scalebarstart = {
    "benchmark/F20-atom_res.preview.png": (29, 944),
    "benchmark/a06fig02a.png": (85,245),
    "benchmark/r6.png": (6,280),
    "benchmark/r7.png": (7,361),
}
benchmarks_scale = {
    "benchmark/F20-atom_res.preview.png": 2.0,
    "benchmark/a06fig02a.png": 20.0,
    "benchmark/r6.png": 100.0,
    "benchmark/r7.png": 50.0,
}
benchmarks_areaunits = {
    "benchmark/F20-atom_res.preview.png": "nm",
    "benchmark/a06fig02a.png": "micron",
    "benchmark/r6.png": "micron",
    "benchmark/r7.png": "micron",
}
benchmarks_xmax = {
    "benchmark/F20-atom_res.preview.png": 2.0,
    "benchmark/a06fig02a.png": 33.0,
    "benchmark/r6.png": 100.0,
    "benchmark/r7.png": 50.0,
}

# These, in turn, should have fewer exceptions than the batch job
gaussian_default = 1
gaussian_custom = {
    "benchmark/a06fig02a.png": 2
}
balance_default = [(25,25), (1,5)],
balance_custom = {
}
threshold_default = 0.7
threshold_custom = {
}

erode_default = True
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
    analysis = SEMAnalysis()
    cropy = benchmarks_cropy[fn]
    start = benchmarks_scalebarstart[fn]
    scale = benchmarks_scale[fn]
    analysis.loadimage(sys.argv[1], cropy=cropy, scalebarstart=start, scale=scale)

    # Do the analysis.
    gaussian = gaussian_custom.get(fn) or gaussian_default
    balance = balance_custom.get(fn) or balance_default
    threshold = threshold_custom.get(fn) or threshold_default
    erode = erode_custom.get(fn) or erode_default
    analysis.filterimage(gaussian=gaussian, balance=balance)
    analysis.thresholdimage(factor=threshold, erode=erode)
    analysis.labelimage()
    analysis.distimage()
    analysis.watershedimage()
    analysis.radialdistribution()
    analysis.fitrdf()

    # Print some data.
    print "**************************"
    print "Number of particles: %i" %analysis.npcount
    print "Particle concentration: %.1f/micron^2" %analysis.npconc
    print "First peak: %.1fnm" %analysis.popt[0]
    print "Pixel size: %.2fnm" %analysis.ps
    print "RDF bin size: %.2fnm" %analysis.dbin
    print "**************************"

    # Do the plotting.
    xmax = benchmarks_xmax[fn]
    units = benchmarks_areaunits[fn]
    analysis.plot(pylab, xmax=xmax, areaunit=units)

    # Save or show.
    if "save" in sys.argv:
        analysis.savearchives()
        analysis.savefigures(pylab)
    else:
        pylab.show()
