"""Plot results of task0013 simulations."""


import bz2
import os
import sys

import numpy

from systems import phases_frames
from task0013 import loadpath

sys.path.append("exp/")
from common import rdfcorrection


def rdfmodel(X, p, m, l, a, b, t, gd):
    """Approximate radial distribution function for a hard sphere liquid.

    Coded according to Matteoli and Mansoori, JPC 1995.
    All parameters have same meaning, but I also added the first peak position, p.
    """
    Y = X/p
    Ym1 = Y-1
    left = gd*numpy.exp(-t*Ym1*Ym1)
    right = 1.0 + numpy.power(Y,-m)*(gd-1.0-l) + numpy.exp(-a*Ym1)*numpy.cos(b*Ym1)*(Ym1+l)/Y
    return (Y<1)*left + (Y>=1)*right


if __name__ == "__main__":

    # Do not use X backend in matplotlib if saving.
    # Backend must be changed before pylab is imported.
    import matplotlib
    if "save" in sys.argv:
        matplotlib.use("Agg")
    import pylab

    # Makefile passes the path to some archive file, from which we should
    #   know what kind of plot to make. But we don't know beforehand what
    #   type of file this will be exactly, or its full name. In hindsight,
    #   that was not a good idea! In any case, it will be a bziped npy archive,
    #   so we can load the data right away.
    fpath = sys.argv[1].strip()
    data = numpy.load(bz2.BZ2File(fpath))

    # Change the path passed to the output file of the simulation,
    #   and load a simulation object from that. This is a lot easier
    #   than writing all the code to decrypt the path over again.
    out = ".".join(fpath.split('.')[:-3]) + ".out"
    sim = loadpath(out, setup=False)

    # Archives with energy data are used to produce energy plots (three different types),
    #   or grid offset plots. The energy plots are quite similar, so we can treat them
    #   ina a single clause, but the offsets option has many additional details.
    # Grid offsets plots can additionally be for angles.
    isenergy = "energy.npy.bz2" in fpath
    if isenergy:
        istotal = "total" in sys.argv
        isfield = "field" in sys.argv
        iscoupling = "coupl" in sys.argv
        isoffsets = "offsets" in sys.argv
        if isoffsets:
            isangles = "angles" in sys.argv
        isinterface = "interface" in sys.argv

    # Histogram plots can be of several kinds. Field histograms can show the distribution
    #   of total field values or the distirbution of order parameters. Radial histograms
    #   are essentialy RDFs, for core or shell beads, full range or zoomed. Residual histograms
    #   show the values of the field or the order parameter at bead coordinates, for cores
    #   or shells. Angular histogram plot the distriubtion of rotation angles for nanoparticles.
    ishist = "hist" in fpath
    if ishist:
        isfield = "hist-field" in fpath
        isradial = "hist-radial" in fpath
        isresidual = "hist-residual" in fpath
        isangles = "hist-ang" in fpath
        if isfield or isresidual:
            istotal = "total" in sys.argv
            isorder = "order" in sys.argv
        if isradial:
            iszoom = "zoom" in sys.argv
        if not isangles:
            isshell = "shell" in sys.argv

    # Flag to save depending on whether appropriate argument is passed.
    issave = "save" in sys.argv

    # ############
    # ENERGY PLOTS
    # ############

    if isenergy and not isoffsets and not isinterface:

        # Things that are specific to each type. Remember that coupling requires
        #   there to be at least one nanoparticle. If no type, then bail out here.
        if istotal:
            instant = data[:,1] + data[:,3] + data[:,4]  + data[:,5]  + data[:,6]
            average = data[:,7] + data[:,9] + data[:,10] + data[:,11] + data[:,12]
            plotfname = "%s.energy-total.png" %sim.outroot
            ylabel = "total free energy"
        elif isfield:
            instant = data[:,3] + data[:,4]  + data[:,5]
            average = data[:,9] + data[:,10] + data[:,11]
            plotfname = "%s.energy-field.png" %sim.outroot
            ylabel = "field free energy"
        elif iscoupling and sim.pop > 0:
            instant = data[:,6]
            average = data[:,12]
            plotfname = "%s.energy-coupl.png" %sim.outroot
            ylabel = "coupling free energy"
        else:
            print "Energy plots must specify total, field or coupling."
            sys.exit(0)

        # Things all types have in common.
        datas = [ "instant.", "average" ]
        dataset = { "instant.": instant, "average": average }

        # With the above set, plotting is now straightforward.
        for ds in datas:
            pylab.plot(data[:,0], dataset[ds], label=ds)
        pylab.xlabel("time step")
        pylab.ylabel(ylabel)

        # Set the X axis range a bit further than the last point.
        xmin = data[0,0]
        xmax = data[-1,0]*1.1

        # Now play around with the range for the Y axis.
        ymin = min([min(dataset[ds][100:]) for ds in dataset])
        ymax = max([max(dataset[ds][100:]) for ds in dataset])
        dy = ymax - ymin
        if isfield:
            ymin = ymax - dy/1.05
            ymax = ymax + dy/100
        else:
            ymin = ymin - dy/100
            ymax = ymin + dy/1.05

        # Finally set the ranges for axes.
        pylab.xlim(xmin,xmax)
        pylab.ylim(ymin,ymax)

        # Let pylab place the legend on its own.
        pylab.legend(loc="best")

    # ############
    # OFFSET PLOTS
    # ############

    if isenergy and isoffsets:

        # These are possible after phase 6 and when there are nanoparticles.
        if sim.phase < 7 or sim.pop == 0:
            print "Offsets plots available after phase 6 and with at aleast one NP."
            sys.exit(0)

        # Angles are quite the same as the initial grid offset plots, but the latter
        #   have more intricacies , since shell beads are involved.
        if isangles:
            mean = data[:,27]
            var = data[:,28]
            skew = data[:,29]
            kurtosis = data[:,30]
            datas = [ "mean", "variance", "skewness", "kurtosis" ]
            dataset = { "mean": mean, "variance": var, "skewness": skew, "kurtosis": kurtosis }
            plotfname = "%s.offsets-ang.png" %sim.outroot
            ylabel = "offset from zero rotation"
        elif sim.phase == 7:
            mean = data[:,19]
            var = data[:,20]
            skew = data[:,21]
            kurtosis = data[:,22]
            datas = [ "mean", "variance", "skewness", "kurtosis" ]
            dataset = { "mean": mean, "variance": var, "skewness": skew, "kurtosis": kurtosis }
            plotfname = "%s.offsets.png" %sim.outroot
            ylabel = "offset relative to grid midpoint"
        else:
            mean = data[:,19]
            mean_shell = data[:,20]
            var = data[:,21]
            var_shell = data[:,22]
            skew = data[:,23]
            skew_shell = data[:,24]
            kurtosis = data[:,25]
            kurtosis_shell = data[:,26]
            datas = [ "mean", "mean (shell)", "variance", "variance (shell)", "skewness", "skewness (shell)", "kurtosis", "kurtosis (shell)" ]
            dataset = { "mean": mean, "mean (shell)": mean_shell, "variance": var, "variance (shell)": var_shell,
                        "skewness": skew, "skewness (shell)": skew_shell, "kurtosis": kurtosis, "kurtosis (shell)": kurtosis_shell }
            plotfname = "%s.offsets.png" %sim.outroot
            ylabel = "offset relative to grid midpoint"

        # With the above set, plotting is now straightforward.
        for ds in datas:
            pylab.plot(data[:,0], dataset[ds], label=ds)
        pylab.xlabel("time step")
        pylab.ylabel(ylabel)

        # Set the X axis range a bit further than the last point.
        xmin = data[0,0]
        xmax = data[-1,0]*1.1

        # Y axis range is quite straightforward.
        ymin = -1.3
        ymax = 0.6
        if isangles:
            ymin = -2.0
            ymax = 4.0

        # Finally set the ranges for axes.
        pylab.xlim(xmin,xmax)
        pylab.ylim(ymin,ymax)

        # Let pylab place the legend on its own.
        pylab.legend(loc="best")

        # For offset plots, we want to add a number of extra eye candy.
        y_mean = 0.5
        y_var = 1.0 / 12
        y_skew = 0.0
        y_kurt = - 6.0 / 5
        str_mean = r"$\frac{1}{2}$"
        str_var = r"$\frac{1}{12}$"
        str_skew = r"$0.0$"
        str_kurt = r"-$\frac{6}{5}$"
        if isangles:
            y_mean = numpy.pi / 2
            y_var = numpy.pi**2 / 12
            str_mean = r"$\frac{\pi}{2}$"
            str_var = r"$\frac{\pi^2}{12}$"

        y_values = (y_mean, y_var, y_skew, y_kurt)
        y_strings = (str_mean, str_var, str_skew, str_kurt)
        for y in y_values:
            pylab.axhline(y=y, xmin=data[0,0], xmax=data[-1,0], linestyle='--', color='gray')
        for y,s in zip(y_values,y_strings):
            pylab.text(xmax*1.01, y, s, fontsize=12, horizontalalignment="left", verticalalignment="center")

    # #######################
    # DISTANCE FROM INTERFACE
    # #######################

    if isenergy and isinterface:

        plotfname = "%s.interface.png" %sim.outroot

        dist_core = data[:,31]
        dist_shell = data[:,32]

        pylab.plot(100*dist_core, label="core beads")
        pylab.plot(100*dist_shell, label="shell beads")

        pylab.xlabel("time step")
        pylab.ylabel(r"NPs at interface [%]")

        pylab.legend(loc="best")

    # ######################
    # SETUP FIELD HISTOGRAMS
    # ######################

    if ishist and isfield:

        if istotal:
            plotfname = "%s.hist-field-total.png" %sim.outroot
            xlabel = "overall field density"
            xmin, xmax = 0.0, 1.5
            data = data[:,0]
        if isorder:
            plotfname = "%s.hist-field-order.png" %sim.outroot
            xlabel = "overall order parameter"
            xmin, xmax = -1.5, 1.5
            data = data[:,1]

    # #######################
    # SETUP RADIAL HISTOGRAMS
    # #######################

    if ishist and isradial:

        # Generate the plot file name.
        # RDFs for shell beads are available after phase 7.
        if isshell:
            if sim.phase > 7:
                data = data[:,1,:]
                plotfname = "%s.hist-radial-shell" %sim.outroot
            else:
                print "RDFs for shell beads are available after phase 7."
                sys.exit(0)
        else:
            data = data[:,0,:]
            plotfname = "%s.hist-radial" %sim.outroot
        if iszoom:
            plotfname += ".zoom.png"
        else:
            plotfname += ".png"

        # The X axis is fixed from analysis, but we will change the plotted scale later.
        xlabel = "NP-NP distance"
        xmin = 0.0
        xmax = 16.0

    # #########################
    # SETUP RESIDUAL HISTOGRAMS
    # #########################
    if ishist and isresidual:

        if isshell and sim.phase < 8:
            print "Shell beads available only after phase 7."
            sys.exit(0)

        if not (istotal or isorder):
            print "Must specify total or order residual histogram type."
            sys.exit(0)

        # Set the plot file name tersely.
        s1 = "total"*istotal + "order"*isorder
        s2 = "-shell"*isshell
        plotfname = "%s.hist-residual-%s%s.png" %(sim.outroot, s1, s2)

        # X axis label and limits. Set the limits to 1.5 for clarity.
        xlabel = "residual field density"*istotal + "residual order parameter"*isorder
        xmin = 0.0*istotal + -1.5*isorder
        xmax = 1.5

        # Data to plot. This is tricky for the order parameter, since there are
        #   additional rows in the array after phase 7 when shell beads are involved.
        if istotal:
            if isshell:
                data = data[:,1]
            else:
                data = data[:,0]
        if isorder:
            if sim.phase > 7:
                if isshell:
                    data = data[:,3]
                else:
                    data = data[:,2]
            else:
                data = data[:,1]

    # ########################
    # SETUP ANGULAR HISTOGRAMS
    # ########################

    if ishist and isangles:

        plotfname = "%s.hist-ang.png" %sim.outroot
        xlabel = "absolute rotation angle [rad]"
        xmin = 0.0
        xmax = 3.3

    # ###############
    # PLOT HISTOGRAMS
    # ###############

    if ishist:

        # Frames are stored in a separate archive file.
        frames = numpy.load(bz2.BZ2File("%s.hist-frames.npy.bz2" %sim.outroot)).tolist()

        # All histograms are averages over some number of samples.
        nsamples = phases_frames[sim.phase][1]

        # The number of bins is implicit.
        nbins = data.shape[1]

        # Bin width and list of bin center points.
        dx = (xmax-xmin)/nbins
        xrange = numpy.arange(xmin+dx/2.0,xmax+dx/2.0,dx)

        # Number of pair interactions in this simulation.
        Npairs = sim.pop*(sim.pop-1)/2

        # Calculate the normalization factors for histograms.
        if isradial:

            # Cut out the strong peaks due to the internal structure of nanoparticles.
            # Must also adjust the number of pair interactions.
            # This very specific hack works for nanoparticle models np4 and np8.
            # Unfortunately, this is still somewhat approximate (but good enough).
            if isshell:
                if sim.npname == "np4":
                    a, b = 0.2828, 0.4
                if sim.npname == "np8":
                    a, b = 0.5657, 0.8
                r1 = numpy.sqrt( (b-a)**2 + a**2 )
                r2 = numpy.sqrt(2.0)*b
                r3 = numpy.sqrt( (b+a)**2 + a**2 )
                r4 = 2*b
                b1 = int(r1//dx)
                b2 = int(r2//dx)
                b3 = int(r3//dx)
                b4 = int(r4//dx)
                data[:,b1] -= 8*sim.pop
                data[:,b2] -= 8*sim.pop
                data[:,b3] -= 8*sim.pop
                data[:,b4] -= 8*sim.pop/2
                Npairs = 8*sim.pop*(8*sim.pop-1) / 2

            
            strips = rdfcorrection(sim.size[0],sim.size[1],xrange) * 2 * numpy.pi * xrange * dx
            factor = sim.size[0] * sim.size[1] / strips / Npairs

        elif isfield or (isresidual and istotal):
            factor = 1.0 / (dx * sum(data[0]))

        elif isresidual and isorder:
            factor = 2.0 / (dx * sum(data[0]))

        elif isangles:
            factor = numpy.pi / (dx * sum(data[0]))

        # This is a quick hack to plot the trend in the position of the largest histogram peak.
        # Should be turned into a proper plot in the future.
        if "trend" in sys.argv:

            D = xrange[numpy.argmax(data[:,:250], axis=1)]
            d = []
            for fi,f in enumerate(phases_frames[13][0]):
                i = phases_frames[13][1]
                d.append(numpy.average(D[fi*i:(fi+1)*i]))
            pylab.plot(phases_frames[13][0], d)

        else:

            frames_to_plot = [1, 11, 21, 101, 201, 1001]
            for f in frames_to_plot:
                hist = factor * numpy.sum([data[frames.index(f+i)] for i in range(nsamples)], axis=0) / nsamples
                pylab.plot(xrange, hist, label="frames %i-%i" %(f,f+nsamples))

            pylab.xlabel(xlabel)
            pylab.ylabel("probability density")

        # Turn on grid and legend.
        pylab.grid()
        pylab.legend()

        # We are most interested in small radial distances for histograms.
        if isradial and not "trend" in sys.argv:
            xmin = 0.0
            if iszoom:
                xmin = 1.35
                if sim.phase > 7:
                    xmax = 3.7
                else:
                    xmax = 2.2
            else:
                xmax = 16.0
            pylab.xlim(xmin,xmax)

    # Finally, save or show the plot.
    if issave:
        pylab.savefig(plotfname)
    else:
        pylab.show()
