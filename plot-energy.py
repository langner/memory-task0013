import bz2
import sys

import numpy
import pylab


if __name__ == "__main__":

    fpath = sys.argv[1].strip()
    data = numpy.load(bz2.BZ2File(fpath))

    instant = data[:,1]+data[:,3]+data[:,4]+data[:,5]+data[:,6]
    average = data[:,7]+data[:,9]+data[:,10]+data[:,11]+data[:,12]
    pylab.plot(data[:,0], instant, label="instant.")
    pylab.plot(data[:,0], average, label="average")

    dy = max(instant) - min(instant)
    ymin = min(instant) - dy/250
    ymax = min(instant) + dy/10
    pylab.ylim(ymin,ymax)
    pylab.legend()

    if "save" in sys.argv:
        pylab.savefig("%s.energy.png" %fpath[:-15])
    else:
        pylab.show()