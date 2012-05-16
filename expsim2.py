import bz2

import numpy
import pylab

from analyze import PhaseFrames
from exp.common import SEMAnalysis, getscalefromname
from exp.plot import fig2
from exp.sem import scalebarstart


sbs = scalebarstart
gsf = getscalefromname
loadsem = lambda fn: SEMAnalysis(fpath=fn, cropy=(0,430), scalebarstart=sbs[gsf(fn)])

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

chosen_exp = ["exp/"+fn for fn in fig2]
chosen_sim = [    #"phase11/64x64x1_A20B16_bv1.00/temp0.01_exp0.10_den1.0_pop128_chmob1.00/k15.0_nchi24.0_ca4.0_cb8.0_mob0.001_a10.0/tt1100_ts0.01.hist-radial.npy.bz2",
"phase11/64x64x1_A20B16_bv1.00/temp0.01_exp0.10_den1.0_pop256_chmob1.00/k15.0_nchi24.0_ca4.0_cb8.0_mob0.001_a10.0/tt1100_ts0.01.hist-radial.npy.bz2",
"phase11/64x64x1_A20B16_bv1.00/temp0.01_exp0.10_den1.0_pop512_chmob1.00/k15.0_nchi24.0_ca4.0_cb8.0_mob0.001_a10.0/tt1100_ts0.01.hist-radial.npy.bz2",
"phase11/64x64x1_A20B16_bv1.00/temp0.01_exp0.10_den1.0_pop768_chmob1.00/k15.0_nchi24.0_ca4.0_cb8.0_mob0.001_a10.0/tt1100_ts0.01.hist-radial.npy.bz2",
]

sems = [loadsem(s) for s in chosen_exp]
radials = [1.0*numpy.load(bz2.BZ2File(fn))[:,0,:] for fn in chosen_sim]

for s in sems:
    s.fitrdf()

nbins = radials[0].shape[1]
size = (64,64)
xmin = 0.0
xmax = 16.0
dx = (xmax-xmin)/nbins
xrange = numpy.arange(xmin+dx/2.0,xmax+dx/2.0,dx)

ns = [256, 512, 768]
npairs = [n*(n-1)/2 for n in ns]
strips = 1.0 * correction(size[0],size[1],xrange) * 2 * numpy.pi * xrange * dx
factors = numpy.array([1.0 * size[0] * size[1] / strips / np for np in npairs])

for ni,n in enumerate(ns):
    for fi,f in enumerate(factors[ni]):
        radials[ni][:,fi] *= f

pylab.figure(1)

exp_d = [s.rdfmax for s in sems]
pylab.plot([0,2,14,48], exp_d[:4], "d-", label="C3")
pylab.plot([0,2,14,48], exp_d[4:8], "d-", label="C6")
pylab.plot([0,2,14,48], exp_d[8:12], "d-", label="C8")

D = [xrange[numpy.argmax(radial[:,:250], axis=1)] for radial in radials]
for di,d in enumerate(D):
    newd = []
    for fi,f in enumerate(PhaseFrames[11][0]):
        i = PhaseFrames[11][1]
        newd.append(numpy.average(d[fi*i:(fi+1)*i]))
        D[di] = newd

X = numpy.array(PhaseFrames[11][0])
X *= 48.0 / 200.0

for di,d in enumerate(D):
    d = numpy.array(d)
    d *= exp_d[3]/d[-1]
    pylab.plot(X, d, "d-", label="%i" %ns[di])

pylab.grid()
pylab.legend()

pylab.figure(2)

exp_g = [s.popt[-1] for s in sems]
pylab.plot([0,2,14,48], exp_g[:4], "d-", label="exp. C3")
pylab.plot([0,2,14,48], exp_g[4:8], "d-", label="exp. C6")
pylab.plot([0,2,14,48], exp_g[8:12], "d-", label="exp. C8")

G = [numpy.max(radial[:,:250], axis=1) for radial in radials]
for gi,g in enumerate(G):
    newg = []
    for fi,f in enumerate(PhaseFrames[11][0]):
        i = PhaseFrames[11][1]
        newg.append(numpy.average(g[fi*i:(fi+1)*i]))
        G[gi] = newg

ns = [256, 512, 768]
for gi,g in enumerate(G):
    g = numpy.array(g)
    pylab.plot(X, g, "d-", label="sim. %i" %ns[gi])

pylab.grid()
pylab.legend(loc="upper left")

pylab.show()