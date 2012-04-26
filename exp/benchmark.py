import sys

from common import pylab, SEMAnalysis


benchmarks_size = {
    "benchmark/F20-atom_res.preview.png": (980,980),
}
benchmarks_cropy = {
    "benchmark/F20-atom_res.preview.png": (0,850),
}
benchmarks_scalebarstart = {
    "benchmark/F20-atom_res.preview.png": (29, 944),
}
benchmarks_scale = {
    "benchmark/F20-atom_res.preview.png": 2.0,
}


if __name__ == "__main__":

    # Create object and load image
    fn = sys.argv[1]
    analysis = SEMAnalysis()
    cropy = benchmarks_cropy[fn]
    start = benchmarks_scalebarstart[fn]
    scale = benchmarks_scale[fn]
    analysis.loadimage(sys.argv[1], cropy=cropy, scalebarstart=start, scale=scale)

    # Do the analysis
    analysis.radialdistribution()
    analysis.fitrdf()

    # Do the plotting
    analysis.plot(xmax=scale)

    # Save or show
    if "save" in sys.argv:
        analysis.savearchives()
        analysis.savefigures()
    else:
        pylab.show()
