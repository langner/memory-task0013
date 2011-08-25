import sys

import numpy

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

    if "archives" in sys.argv:

        sim = loadpath(out)

        cga = "%s/%s.cga" %(sim.outpath,sim.outname)
        csa = "%s/%s.csa" %(sim.outpath,sim.outname)

        archive = UtilitiesManager.CreateArchiveReader()
        archive.SetSystem(sim.box, cga, csa)

        nframes = archive.GetNumberOfFrames()
        nsamples = 128
        npoints = nframes/nsamples
        ifreq = nframes/npoints

        hists_totalfields, totalfields = [], []
        hists_orderparams, orderparams = [], []
        hists_radialdists, radialdists = [], []
        for nframe in range(nframes):

            archive.LoadFrame(nframe+1)
            coordinates = numpy.zeros((sim.population,2))
            for np in range(sim.population):
                x,y,z = sim.getcenterofmass("np",np)
                coordinates[np,0] = x % 64
                coordinates[np,1] = y % 64
                totalfields.append(sim.gettotalfield(*coordinates[np]))
                orderparams.append(sim.getorderparameter(*coordinates[np]))
                for npp in range(np):
                    dr = coordinates[npp]-coordinates[np]
                    radialdists.append(numpy.linalg.norm(dr%32.0))

            if nframe == 0:
                pass
            elif nframe % ifreq == 0:
                hists_totalfields.append(numpy.histogram(totalfields, bins=64, range=[0.0,1.0])[0])
                hists_orderparams.append(numpy.histogram(orderparams, bins=64, range=[-1.0,1.0])[0])
                hists_radialdists.append(numpy.histogram(radialdists, bins=64, range=[0.0,48.0])[0])
                totalfields = []
                orderparams = []
                radialdists = []

        tosave = numpy.array([hists_totalfields, hists_orderparams, hists_radialdists])
        numpy.save(out.replace(".out", ".histograms.npy"), tosave)
