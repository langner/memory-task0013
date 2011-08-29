import bz2
import sys

import numpy
import pylab


if __name__ == "__main__":

    fpath = sys.argv[1].strip()
    data = numpy.load(bz2.BZ2File(fpath))

    etype = sys.argv[2].strip()
    if etype == "total":
        instant = data[:,1]+data[:,3]+data[:,4]+data[:,5]+data[:,6]
        average = data[:,7]+data[:,9]+data[:,10]+data[:,11]+data[:,12]
        plotfname = "%s.energy-total.png" %fpath[:-15]
    if etype == "field":
        instant = data[:,3]+data[:,4]+data[:,5]
        average = data[:,9]+data[:,10]+data[:,11]
        plotfname = "%s.energy-field.png" %fpath[:-15]
    if etype == "coupling":
        instant = data[:,6]
        average = data[:,12]
        plotfname = "%s.energy-coupling.png" %fpath[:-15]

    pylab.plot(data[:,0], instant, label="instant.")
    pylab.plot(data[:,0], average, label="average")

    ymin = min(instant[100:])
    ymax = max(instant[100:])
    dy = ymax - ymin
    if etype in ("total", "coupling"):
        ymin = ymin - dy/250
        ymax = ymin + dy/2.5
    if etype == "field":
        ymin = ymax - dy/1.5
        ymax = ymax + dy/250
    pylab.ylim(ymin,ymax)
    pylab.legend()

    if "save" in sys.argv:
        pylab.savefig(plotfname)
    else:
        pylab.show()