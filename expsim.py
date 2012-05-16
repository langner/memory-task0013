import bz2

import numpy
import pylab

from exp.common import SEMAnalysis, getscalefromname
from exp.sem import scalebarstart


sbs = scalebarstart
gsf = getscalefromname
loadsem = lambda fn: SEMAnalysis(fpath=fn, cropy=(0,430), scalebarstart=sbs[gsf(fn)])

chosen_exp = "exp/sem/5k/C6/48-Image5-500.jpg"
chosem_sim = "phase11/64x64x1_A20B16_bv1.00/temp0.01_exp0.10_den1.0_pop512_chmob1.00/k15.0_nchi24.0_ca4.0_cb8.0_mob0.001_a10.0/tt1100_ts0.01.hist-radial.npy.bz2"

sem = loadsem(chosen_exp)
sem.fitrdf()

radial = 1.0 * numpy.load(bz2.BZ2File(chosem_sim))[:,0,:]

def correction(X,Y,R):
    ang1 = 2 * numpy.pi
    ang2 = 2 * numpy.pi - 2.0
    ang3 = 2 * numpy.pi - 2.0
    ang4 = 3 * numpy.pi / 2.0 - 2.0
    V1 = (0.5*X - R)*(0.5*Y - R)
    V2 = R*(0.5*X - R)
    V3 = R*(0.5*Y - R)
    V4 = R**2
    ref = (V1+V2+V3+V4) * 2 * numpy.pi
    return (ang1*V1 + ang2*V2 + ang3*V3 + ang4*V4) / ref

nbins = radial.shape[1]
size = (64,64)
xmin, xmax = 0.0, 16.0
dx = (xmax-xmin)/nbins
xrange = numpy.arange(xmin+dx/2.0,xmax+dx/2.0,dx)
N = 512
Npairs = N*(N-1)/2

strips = 1.0 * correction(size[0],size[1],xrange) * 2 * numpy.pi * xrange * dx
factor = 1.0 * size[0] * size[1] / strips / Npairs
radial *= factor

iplot = -1
xopt = xrange[radial[iplot].argmax()-1]
xrange *= sem.rdfmax / xopt

pylab.plot(sem.X, sem.Y, color="red", label="experimental")
pylab.plot(sem.X, sem.rdfmodel(sem.X), color="blue", label="exp. fit")

#pylab.plot(xrange, numpy.sum(radial[0:10], axis=0)/10, color="green")
pylab.plot(xrange, numpy.sum(radial[iplot-25:iplot], axis=0)/25, color="gray", label="simulation")

pylab.xlim([5,150])
pylab.ylim([0,9])

pylab.xlabel("NP-NP distance [nm]")
pylab.ylabel("g(r)")

pylab.legend(loc=[0.6,0.25])
pylab.grid()

ax = pylab.axes([0.50, 0.50, 0.35, 0.35], axisbg='y')
pylab.plot(sem.X, sem.Y, color="black")
pylab.xticks([100, 500])

pylab.show()