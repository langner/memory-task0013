import sys

import numpy
import pylab


if __name__ == "__main__":

    fpath = sys.argv[1].strip()
    data = numpy.load(fpath)
    froot = fpath.replace(".data-analyzed.npz","")

    if "energy" in sys.argv:

        E = data['energy']

        if "total" in sys.argv:
            instant = E[:,1]+E[:,3]+E[:,4]+E[:,5]+E[:,6]
            average = E[:,7]+E[:,9]+E[:,10]+E[:,11]+E[:,12]
            plotfname = "%s.energy-total.png" %froot
            ylabel = "total free energy"
        if "field" in sys.argv:
            instant = E[:,3]+E[:,4]+E[:,5]
            average = E[:,9]+E[:,10]+E[:,11]
            plotfname = "%s.energy-field.png" %froot
            ylabel = "field free energy"
        if "coupl" in sys.argv:
            instant = E[:,6]
            average = E[:,12]
            plotfname = "%s.energy-coupl.png" %froot
            ylabel = "coupling free energy"

        pylab.plot(E[:,0], instant, label="instant.")
        pylab.plot(E[:,0], average, label="average")

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

    if "hist" in sys.argv:

        frames = data["hist_frames"].tolist()

        if "totals" in sys.argv:

            H = data["hists_totals"]
            nbins = H.shape[1]
            plotfname = froot+".hist-totals.png"
            dmax = numpy.sqrt(2.0)*32
            dx = dmax/nbins
            xrange = numpy.arange(dx,dmax+dx,dx)

        if "radials" in sys.argv:

            H = data["hists_radials"]
            nbins = H.shape[1]
            plotfname = froot+".hist-radials.png"
            dmax = numpy.sqrt(2.0)*32
            dx = dmax/nbins
            xrange = numpy.arange(dx,dmax+dx,dx)
            S = 1.0*sum(H[0])
            pylab.plot(xrange, H[frames.index(0)]/S, label="frame 0")
            hist = numpy.sum([H[frames.index(1+i)] for i in range(10)], axis=0)/S/10.0
            pylab.plot(xrange, hist, label="frames 1-10")
            hist = numpy.sum([H[frames.index(101+i)] for i in range(10)], axis=0)/S/10.0
            pylab.plot(xrange, hist, label="frames 101-110")
            hist = numpy.sum([H[frames.index(1001+i)] for i in range(10)], axis=0)/S/10.0
            pylab.plot(xrange, hist, label="frames 1001-1010")
            hist = numpy.sum([H[frames.index(10001+i)] for i in range(10)], axis=0)/S/10.0
            pylab.plot(xrange, hist, label="frames 10001-10010")
            pylab.xlabel("NP-NP distance")
            pylab.ylabel("number of pairs")
            pylab.xlim(0.0,4.0)
            pylab.grid()
            pylab.legend()

    if "save" in sys.argv:
        pylab.savefig(plotfname)
    else:
        pylab.show()