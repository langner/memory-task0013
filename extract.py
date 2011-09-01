import sys
import time

import numpy

from scipy import histogram
from scipy.spatial.distance import pdist

from simulation import *

energies = ['nonbonded', 'inhomo', 'ideal', 'contact', 'compress', 'coupling']
columns = dict(zip(energies,[4, 8, 10, 11, 12, 14]))


if __name__ == "__main__":

    fname = sys.argv[1]
    fext = fname[-4:]

    if fext == ".ctf":
        fout = fname[:-9] + ".out"
    else:
        fout = fname[:-4] + ".out"

    if not "Time used" in open(fout).read():

        print "This simulation has not finished."
        sys.exit(1)

    if fext == ".ctf":

        data = numpy.loadtxt(fname, skiprows=22)

        nsamples = 256
        npoints = data.shape[0]/nsamples
        ifreq = data.shape[0]/npoints

        indices = data[::ifreq,0]
        instants = [data[::ifreq,columns[en]] for en in energies]
        averages = [[numpy.average(data[i*ifreq:(i+1)*ifreq,columns[en]]) for i in range(npoints+1)] for en in energies]
        averages = [numpy.array(avg) for avg in averages]
        deviations = [[numpy.std(data[i*ifreq:(i+1)*ifreq,columns[en]]) for i in range(npoints+1)] for en in energies]
        deviations = [numpy.array(dev) for dev in deviations]

        tosave = numpy.vstack([indices]+instants+averages+deviations).transpose()
        numpy.save(fout.replace(".out", ".energy.npy"), tosave)

    if "histograms" in sys.argv:

        t = time.time()

        sim = loadpath(fout)

        cga = "%s/%s.cga" %(sim.outpath,sim.outname)
        csa = "%s/%s.csa" %(sim.outpath,sim.outname)

        archive = UtilitiesManager.CreateArchiveReader()
        archive.SetSystem(sim.box, cga, csa)

        nframes = archive.GetNumberOfFrames()
        nsamples = 16
        npoints = nframes/nsamples
        ifreq = nframes/npoints

        pairs = sim.population*(sim.population-1)/2
        hists_totalfields, totalfields = [], []
        hists_orderparams, orderparams = [], []
        hists_radialdists, radialdists = [], []

        npi1, npi2 = [], []
        for i1 in range(1,sim.population):
            for i2 in range(i1):
                npi1.append(i1)
                npi2.append(i2)
        #print 'init:', time.time()-t

        T1, T2, T3 = 0.0, 0.0, 0.0
        frames = [0] + [2**n for n in range(100) if 2**n <= nframes-nsamples-1] + [nframes-nsamples-1]
        #print 'frames:', len(frames)
        for nframe in frames:

            for iframe in range(nframe,nframe+nsamples):

                t = time.time()
                archive.LoadFrame(iframe+1)
                coordinates = numpy.array([sim.getcenterofmass("np", i)[:2] for i in range(sim.population)])
                coordinates %= 64.0
                T1 += time.time() - t

                t = time.time()
                totalfields.append([sim.gettotalfield(*c) for c in coordinates])
                orderparams.append([sim.getorderparameter(*c) for c in coordinates])
                #drs = (coordinates.take(npi1, axis=0) - coordinates.take(npi2, axis=0)) % 32.0
                #radials = numpy.array([numpy.linalg.norm(dr) for dr in drs])
                #radials = numpy.array([numpy.sqrt(dr[0]*dr[0] + dr[1]*dr[1]) for dr in drs])
                #radials = numpy.array([numpy.sqrt(numpy.dot(dr,dr)) for dr in drs])
                #radials = numpy.sqrt(numpy.array([numpy.dot(dr,dr) for dr in drs]))
                #radials = numpy.sqrt((drs*drs).sum(axis=1))
                #radialdists.extend(numpy.sqrt((drs*drs).sum(axis=1)))
                radialdists.append(pdist(coordinates, 'euclidean'))
                T2 += time.time() - t

            t = time.time()
            totalfields = numpy.array(totalfields).flatten()
            orderparams = numpy.array(orderparams).flatten()
            radialdists = numpy.array(radialdists).flatten()
            hists_totalfields.append(histogram(totalfields, bins=64, range=[0.0,1.0])[0])
            hists_orderparams.append(histogram(orderparams, bins=64, range=[-1.0,1.0])[0])
            hists_radialdists.append(histogram(radialdists, bins=64, range=[0.0,48.0])[0])
            totalfields = []
            orderparams = []
            radialdists = []
            #filter = SelectorsManager.CreateSoftCoreBeadsFilter()
            #filter.AddName("P", "INCLUDE")
            #analyzer = AnalyzersManager.CreateSCMAnalyzer()
            #analyzer.SetSystem(sim.box, csa)
            #analyzer.SetFrames("%i-%i" %(nframe-10,nframe))
            #PC = analyzer.PairCorrelation(filter, filter, 0.1)
            T3 += time.time() - t

        #print 'coordinates:', T1
        #print 'indexing:', T2
        #print 'histograms:', T3

        t = time.time()
        tosave = numpy.array([hists_totalfields, hists_orderparams, hists_radialdists])
        numpy.save(fout.replace(".out", ".histograms.npy"), tosave)
        #print 'save:', time.time() - t
