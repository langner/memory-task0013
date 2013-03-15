import glob
import itertools
import os
import sys

import numpy

from common import SEMAnalysis


array = numpy.array

# Ignore magnifications and concentrations other than these.
# There are ten concentration/series which are called
# simply C1-C10, but denote (usually) a collection of images
# taken from a single sample.
magnifications = [200, 500, 1000]
concentrations = ["C%i" %i for i in range(1,11)]
ligand_lengths = [5, 20]
times = [0, 2, 14, 48]

# The cropping height should be identical in all images.
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
    "sem/5k/C2/48-Image15-200.jpg" : 3,
    "sem/5k/C2/48-Image23-200.jpg" : 3,
    "sem/5k/C5/48-Image35-200.jpg" : 3,
    "sem/5k/C6/48-Image40-200.jpg" : 3,
    "sem/5k/C9/2-Image12-200.jpg" : 3,
    "sem/5k/C9/2-Image15-500.jpg" : 2,
    "sem/5k/C10/2-Image36-200.jpg" : 3,
    "sem/20k/C2/48-Image4-200.jpg" : 3,
    "sem/20k/C3/0-Image22-500.jpg" : 2,
    "sem/20k/C3/0-Image23-500.jpg" : 2,
    "sem/20k/C3/0-Image38-500.jpg" : 2,
    "sem/20k/C4/0-Image2-500.jpg" : 2,
    "sem/20k/C5/0-Image4-200.jpg" : 3,
    "sem/20k/C5/0-Image6-1000.jpg" : 1,
    "sem/20k/C5/0-Image24-500.jpg" : 2,
    "sem/20k/C5/0-Image25-1000.jpg" : 2,
    "sem/20k/C5/48-Image23-500.jpg" : 2,
    "sem/20k/C6/0-Image1-500.jpg": 2,
    "sem/20k/C6/2-Image35-500.jpg": 2,
    "sem/20k/C9/0-Image15-1000.jpg" : 1,
    "sem/20k/C9/14-Image13-1000.jpg" : 1,
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
    "sem/20k/C6/0-Image1-500.jpg": [(10,20), (1,10)],
    "sem/20k/C6/2-Image35-500.jpg": [(10,20), (1,10)],
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
    "sem/5k/C1/0-Image16-200.jpg" : 1.1,
    "sem/5k/C1/14-Image16-500.jpg" : 1.3,
    "sem/5k/C1/14-Image17-1000.jpg" : 1.3,
    "sem/5k/C2/0-Image2-200.jpg" : 1.1,
    "sem/5k/C2/14-Image2-200.jpg" : 1.1,
    "sem/5k/C2/48-Image15-200.jpg" : 1.1,
    "sem/5k/C2/48-Image23-200.jpg" : 1.1,
    "sem/5k/C3/14-Image3-200.jpg" : 0.7,
    "sem/5k/C3/14-Image6-500.jpg" : 0.9,
    "sem/5k/C4/0-Image6-1000.jpg" : 1.0,
    "sem/5k/C5/0-Image9-1000.jpg" : 1.0,
    "sem/5k/C5/48-Image35-200.jpg" : 1.1,
    "sem/5k/C6/0-Image20-1000.jpg": 0.9,
    "sem/5k/C6/0-Image23-1000.jpg": 0.9,
    "sem/5k/C6/2-Image3-1000.jpg" : 1.5,
    "sem/5k/C6/48-Image40-200.jpg" : 1.4,
    "sem/5k/C7/2-Image7-1000.jpg": 0.9,
    "sem/5k/C7/2-Image8-200.jpg": 1.3,
    "sem/5k/C7/48-Image14-1000.jpg": 1.0,
    "sem/5k/C8/48-Image18-200.jpg" : 1.0,
    "sem/5k/C8/48-Image20-1000.jpg" : 1.1,
    "sem/5k/C9/2-Image12-200.jpg" : 1.0,
    "sem/5k/C10/0-Image34-1000.jpg" : 1.2,
    "sem/5k/C10/0-Image36-200.jpg" : 1.5,
    "sem/5k/C10/48-Image25-1000.jpg" : 1.0,
    "sem/20k/C1/0-Image37-1000.jpg" : 1.4,
    "sem/20k/C1/14-Image20-1000.jpg" : 1.1,
    "sem/20k/C1/14-Image23-1000.jpg" : 1.1,
    "sem/20k/C2/2-Image31-1000.jpg" : 1.1,
    "sem/20k/C2/48-Image4-200.jpg" : 1.2,
    "sem/20k/C3/0-Image22-500.jpg" : 1.6,
    "sem/20k/C3/0-Image23-500.jpg" : 1.6,
    "sem/20k/C3/0-Image38-500.jpg" : 1.8,
    "sem/20k/C3/0-Image39-1000.jpg" : 1.5,
    "sem/20k/C3/48-Image6-1000.jpg" : 1.0,
    "sem/20k/C4/0-Image3-1000.jpg" : 1.1,
    "sem/20k/C4/2-Image33-1000.jpg" : 1.1,
    "sem/20k/C4/48-Image16-1000.jpg" : 1.1,
    "sem/20k/C5/0-Image4-200.jpg" : 1.8,
    "sem/20k/C5/0-Image6-1000.jpg" : 1.6,
    "sem/20k/C5/0-Image25-1000.jpg" : 1.0,
    "sem/20k/C5/2-Image18-1000.jpg" : 0.8,
    "sem/20k/C5/14-Image32-1000.jpg" : 1.0,
    "sem/20k/C5/48-Image18-200.jpg" : 2.5,
    "sem/20k/C5/48-Image23-500.jpg" : 1.2,
    "sem/20k/C5/48-Image25-200.jpg" : 2.5,
    "sem/20k/C5/48-Image27-500.jpg" : 2.5,
    "sem/20k/C6/0-Image1-500.jpg" : 1.6,
    "sem/20k/C6/0-Image2-1000.jpg" : 1.3,
    "sem/20k/C6/2-Image35-500.jpg" : 1.6,
    "sem/20k/C6/2-Image36-1000.jpg" : 1.3,
    "sem/20k/C6/2-Image37-200.jpg" : 1.3,
    "sem/20k/C6/14-Image1-200.jpg" : 1.3,
    "sem/20k/C6/48-Image6-1000.jpg" : 1.3,
    "sem/20k/C8/0-Image11-500.jpg" : 1.6,
    "sem/20k/C8/0-Image9-500.jpg" : 1.6,
    "sem/20k/C9/0-Image15-1000.jpg" : 2.0,
    "sem/20k/C9/14-Image13-1000.jpg" : 2.0,
    "sem/20k/C10/0-Image11-500.jpg" : 1.4,
    "sem/20k/C10/0-Image15-1000.jpg" : 1.4,
    "sem/20k/C10/2-Image33-500.jpg" : 1.4,
    "sem/20k/C10/48-Image25-500.jpg" : 1.4,
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

# Certain expensive analysis results are calculated in a separate execution branch
# and saved for later use in figures.
fname_results = "%s/figs.results.txt" %os.path.dirname(os.path.abspath(__file__))
def retrieve_results():
    lines = open(fname_results).readlines()
    isresults = lambda l: not ("ligand" in l or "conc" in l)
    results = [map(float,l.split()[-4:]) for l in lines if isresults(l)]
    return numpy.array(results).reshape((2,10,11,4))

def getscalefromname(fn):
    """Get task-specific scale from file name."""
    return int(fn.split('-')[2].split('.')[0])

def loadsem(fn, scalebarstart=scalebarstart):
    print fn
    sbs = scalebarstart[getscalefromname(fn)]
    return SEMAnalysis(fpath=fn, cropy=(0,430), scalebarstart=sbs)


if __name__ == "__main__":

    # ########################
    # Compute and save results
    # ########################

    # This double loop runs over all ligand lengths and concentration/series,
    # and does intensive analyses on all images available.
    # We save these averages/uncertainties to a text file that can be read
    # later on when drawing figures.
    if "compute" in sys.argv:

        root = os.path.dirname(fname_results)
        f = open(fname_results, "w")

        chain = lambda l: list(itertools.chain(*l))

        # This helper function uses globs to find all filenames that match
        # a triple of ligand length (k), concentration (s) and time (t).
        fmt = root + "/sem/%ik/%s/%i-*-%i.jpg"
        get_kst = lambda k, s, t: chain([glob.glob(fmt % (k, s, t, m)) for m in magnifications])

        for length in ligand_lengths:

            print >>f, "%ik ligands" %length

            for conc in concentrations:

                print >>f, "conc %s" %conc

                # File names and images are loaded into nested lists,
                # and we need to analyze all images.
                fnames = [get_kst(length, conc, t) for t in times]
                sems = [[loadsem(fn) for fn in t] for t in fnames]
                for s in chain(sems):
                    s.clusterparticles()
                    s.radialdistribution()
                    s.fitrdf()

                # These are helper functions for extracting the numerical
                # attributes of nested lists of images, and calculating
                # their averages/devation. The areal densities of nanoparticles
                # are used as weights(or something proportional).
                # When no data is available, insert a negative value. 
                counts = [len(s) for s in sems]
                weights = [[1.0*im.npcount/im.ps for im in s] or 1 for s in sems]
                getall = lambda p: [[getattr(im,p) for im in s] for s in sems]
                getavg = lambda p: [numpy.average(s or -1, weights=w) for s,w in zip(getall(p),weights)]
                getstd = lambda p: [numpy.std(s, ddof=1) for s in getall(p)]

                npconc_avg = getavg("npconc")
                npconc_std = getstd("npconc")

                # These quantities we want to have as functions of time.
                #   pp - peak positions
                #   ph - average peak height
                #   cs - inverse cluster size
                #   cd - cluster density per squared micron
                pp_avg = getavg("rdfmax")
                ph_avg = getavg("rdfheight")
                pp_err_std = getstd("rdfmax")
                pp_err_fit = getavg("rdfmax_var")
                ph_err_std = getstd("rdfheight")
                ph_err_fit = getavg("rdfheight_var")
                cs_avg = getavg("cluster_size_inv")
                cd_avg = getavg("clusters_per_micron_inv")
                cs_std = getstd("cluster_size_inv")
                cd_std = getstd("clusters_per_micron_inv")

                # The total deviation we estimate from the sum of the fitting
                # error (which is a variance) and the standard deviation of the
                # fitting parameter.
                pp_err_total = numpy.sqrt(array(pp_err_std)**2 + array(pp_err_fit))
                ph_err_total = numpy.sqrt(array(ph_err_std)**2 + array(ph_err_fit))

                print >>f, "Weights                 %s " %(" ".join(map(str, [numpy.sum(w) for w in weights])))
                print >>f, "Concentration average   %s " %(" ".join(map(str, npconc_avg)))
                print >>f, "Concentration error     %s " %(" ".join(map(str, npconc_std)))
                print >>f, "Peak position average   %s " %(" ".join(map(str, pp_avg)))
                print >>f, "Peak position error     %s " %(" ".join(map(str, pp_err_total)))
                print >>f, "Peak height average     %s " %(" ".join(map(str, ph_avg)))
                print >>f, "Peak height error       %s " %(" ".join(map(str, ph_err_total)))
                print >>f, "Cluster size average    %s " %(" ".join(map(str, cs_avg)))
                print >>f, "Cluster size error      %s " %(" ".join(map(str, cs_std)))
                print >>f, "Cluster density average %s " %(" ".join(map(str, cd_avg)))
                print >>f, "Cluster density error   %s " %(" ".join(map(str, cd_std)))

        f.close()

    else:

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
        analysis.clusterparticles()
        analysis.radialdistribution()
        try:
            analysis.fitrdf()
        except:
            print "Error fitting to RDF"

        # Print some data.
        print "**************************"
        print "Number of particles: %i" %analysis.npcount
        print "Particle concentration: %.1f/micron^2" %analysis.npconc
        print "**************************"
        print "RDF bin size: %.2fnm" %analysis.dbin
        print "First RDF peak: %.1fnm" %analysis.rdfmax
        print "Pixel size: %.2fnm" %analysis.ps
        print "**************************"
        print "Cluster threshold: %.3f pixels" %analysis.cluster_threshold
        print "Cluster threshold: %.3f nm" %(analysis.cluster_threshold*analysis.ps)
        print "Number of clusters: %i" %analysis.cluster_count
        print "Average cluster size: %.2f NPs" %analysis.cluster_size
        print "Cluster per micron^2: %.2f" %analysis.clusters_per_micron
        print "**************************"

        # Do the plotting.
        analysis.plot(pylab, xmax=150)

        # Save or show.
        if "save" in sys.argv:
            analysis.savearchives()
            analysis.savefigures(pylab)
        else:
            pylab.show()
