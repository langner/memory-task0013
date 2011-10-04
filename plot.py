import bz2
import sys

import numpy
import pylab


if __name__ == "__main__":

    fpath = sys.argv[1].strip()
    data = numpy.load(bz2.BZ2File(fpath))

    if "energy.npy.bz2" in fpath:

        froot = fpath.replace(".energy.npy.bz2","")
        E = data

        if "total" in sys.argv:
            instant = data[:,1] + data[:,3] + data[:,4]  + data[:,5]  + data[:,6]
            average = data[:,7] + data[:,9] + data[:,10] + data[:,11] + data[:,12]
            plotfname = "%s.energy-total.png" %froot
            ylabel = "total free energy"
        if "field" in sys.argv:
            instant = data[:,3] + data[:,4]  + data[:,5]
            average = data[:,9] + data[:,10] + data[:,11]
            plotfname = "%s.energy-field.png" %froot
            ylabel = "field free energy"
        if "coupl" in sys.argv:
            instant = data[:,6]
            average = data[:,12]
            plotfname = "%s.energy-coupl.png" %froot
            ylabel = "coupling free energy"

        pylab.plot(data[:,0], instant, label="instant.")
        pylab.plot(data[:,0], average, label="average")

        ymin = min(instant[100:])
        ymax = max(instant[100:])
        dy = ymax - ymin
        if "field" in sys.argv:
            ymin = ymax - dy/1.05
            ymax = ymax + dy/250
        else:
            ymin = ymin - dy/250
            ymax = ymin + dy/1.05

        pylab.xlabel("time step")
        pylab.ylabel(ylabel)

        pylab.ylim(ymin,ymax)
        pylab.legend()

    if "hist-field.npy.bz2" in fpath:
        froot = fpath.replace(".hist-field.npy.bz2","")

        if "total" in sys.argv:
            plotfname = froot+".hist-field-total.png"
            xlabel = "overall field density"
            xmin, xmax = 0.0, 1.5
            data = data[:,0]

        if "order" in sys.argv:
            plotfname = froot+".hist-field-order.png"
            xlabel = "overall order parameter"
            xmin, xmax = -1.5, 1.5
            data = data[:,1]

    if "hist-radial.npy.bz2" in fpath:
        froot = fpath.replace(".hist-radial.npy.bz2","")
        plotfname = froot+".hist-radial.png"
        xlabel = "NP-NP distance"
        xmin, xmax = 0.0, 16.0

    if "hist-residual.npy.bz2" in fpath:
        froot = fpath.replace(".hist-residual.npy.bz2","")

        if "total" in sys.argv:
            plotfname = froot+".hist-residual-total.png"
            xlabel = "residual field density"
            xmin, xmax = 0.0, 1.5
            data = data[:,0]

        if "order" in sys.argv:
            plotfname = froot+".hist-residual-order.png"
            xlabel = "residual order parameter"
            xmin, xmax = -1.5, 1.5
            data = data[:,1]


    if "hist" in fpath:

        S = 1.0*sum(data[0])
        nbins = data.shape[1]
        dx = (xmax-xmin)/nbins
        xrange = numpy.arange(xmin+dx/2.0,xmax+dx/2.0,dx)
        frames = numpy.load(bz2.BZ2File(froot+".hist-frames.npy.bz2")).tolist()

        if "hist-radial.npy" in fpath:
            pylab.plot(xrange, data[frames.index(0)]/S, label="frame 0")

        hist = numpy.sum([data[frames.index(1+i)] for i in range(16)], axis=0)/S/16.0
        pylab.plot(xrange, hist, label="frames 1-16")
        hist = numpy.sum([data[frames.index(101+i)] for i in range(16)], axis=0)/S/16.0
        pylab.plot(xrange, hist, label="frames 101-116")
        hist = numpy.sum([data[frames.index(1001+i)] for i in range(16)], axis=0)/S/16.0
        pylab.plot(xrange, hist, label="frames 1001-1016")
        hist = numpy.sum([data[frames.index(10001+i)] for i in range(16)], axis=0)/S/16.0
        pylab.plot(xrange, hist, label="frames 10001-10016")
        if 50001 in frames:
            hist = numpy.sum([data[frames.index(50001+i)] for i in range(16)], axis=0)/S/16.0
            pylab.plot(xrange, hist, label="frames 50001-50016")
        pylab.xlabel(xlabel)
        pylab.ylabel("probability density")

        pylab.grid()
        pylab.legend()

        pylab.xlim(xmin,xmax)
        if "hist-radial.npy" in fpath and "zoom" in sys.argv:
            plotfname = froot+".hist-radial-zoom.png"
            pylab.xlim(0.0, 2.2)

    if "save" in sys.argv:
        pylab.savefig(plotfname)
    else:
        pylab.show()