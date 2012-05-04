#Fig.2

import glob
import string
import sys

import numpy

from common import pylab, SEMAnalysis, getscalefromname
from sem import cropy, scalebarstart


fig1 = [
        "sem/5k/C6/14-Image35-500.jpg",
        "sem/5k/C6/14-Image37-200.jpg",
        "sem/20k/C3/14-Image25-500.jpg",
        "sem/20k/C3/14-Image24-200.jpg",
]

fig2_conc = [3, 6, 8]
fig2 = [
        "sem/5k/C3/0-Image4-500.jpg",
        "sem/5k/C3/2-Image5-500.jpg",
        "sem/5k/C3/14-Image6-500.jpg",
        "sem/5k/C3/48-Image26-500.jpg",
        "sem/5k/C6/0-Image21-500.jpg",
        "sem/5k/C6/2-Image1-500.jpg",
        "sem/5k/C6/14-Image35-500.jpg",
        "sem/5k/C6/48-Image5-500.jpg",
        "sem/5k/C8/0-Image27-500.jpg",
        "sem/5k/C8/2-Image29-500.jpg",
        "sem/5k/C8/14-Image56-500.jpg",
        "sem/5k/C8/48-Image17-500.jpg",
]

fig3_conc = [1, 3, 5]
fig3 = [
        "sem/20k/C1/0-Image36-500.jpg",
        "sem/20k/C1/2-Image9-500.jpg",
        "sem/20k/C1/14-Image22-500.jpg",
        "sem/20k/C1/48-Image43-500.jpg",
        "sem/20k/C3/0-Image23-500.jpg",
        "sem/20k/C3/2-Image13-500.jpg",
        "sem/20k/C3/14-Image24-500.jpg",
        "sem/20k/C3/48-Image10-500.jpg",
        "sem/20k/C5/0-Image24-500.jpg",
        "sem/20k/C5/2-Image17-500.jpg",
        "sem/20k/C5/14-Image33-500.jpg",
        "sem/20k/C5/48-Image27-500.jpg",
]

fig5 = [
        "sem/5k/C1/2-Image2-200.jpg",
        "sem/5k/C1/48-Image5-200.jpg",
]

# Times at which SEMS images were registered
times = [0, 2, 14, 48]

# Convenience function for loading SEM analyses
sbs, gsf = scalebarstart, getscalefromname
loadsem = lambda fn: SEMAnalysis(fpath=fn, cropy=(0,430), scalebarstart=sbs[gsf(fn)])

# Other convenience functions
array = numpy.array
sqrt = numpy.sqrt
avg = numpy.average
std = numpy.std


if __name__ == "__main__":

    # Fig.2 and Fig.3
    isfig2 = "fig2" in sys.argv
    isfig3 = "fig3" in sys.argv
    if isfig2 or isfig3:

        # Number of columns and rows
        ncols = 3
        nrows = 4

        # Load and analyze all images
        to_plot = fig2*isfig2 + fig3*isfig3
        F = [loadsem(fn) for fn in to_plot]
        for f in F:
            f.radialdistribution()
            f.fitrdf()

        # Notice how things are parametrized by number of columns and rows
        pylab.figure(figsize=(4.8*ncols,3.25*nrows))
        for fi,f in enumerate(F):

            # Compute the coordinates
            # Note that we want top-down order, whereas matplotlib does left-right
            x = (1.0/ncols)*(fi//nrows)
            y = 1.0*(nrows-1)/nrows - 1.0*(fi%nrows)/nrows

            # Show image, without ticks
            a = pylab.axes([x+0.005, y+0.005, 1.0/ncols-0.01, 1.0/nrows-0.005])
            pylab.setp(a, xticks=[], yticks=[])
            pylab.imshow(f.img, cmap='gray')

            # Determine scaling (same in all plots)
            ymax = max(3.0, max([max(ff.rdf) for ff in F if hasattr(ff, 'rdf')]))

            # Plot the RDF in the upper right corner
            # Plot the fitted function over the experimental data
            a = pylab.axes([x+1.0/ncols*2/3, y+1.0/nrows*0.6, 1.0/ncols/3-0.01, 1.0/nrows*0.4-0.01])
            pylab.setp(a, xticks=[50], xticklabels=[""], yticks=[], xlim=[0, 30*f.ps], ylim=[0, ymax])
            pylab.plot(f.X, f.Y, label="%ih" %f.time, color="black", linewidth=1.00)
            pylab.plot(f.X, f.rdfmodel(f.X), color="black", linewidth=1.75)

            txt_label = string.ascii_lowercase[fi]
            txt_sigma = "$\mathrm{\sigma}$=%i/$\mathrm{\mu}$m$^2$" %f.npconc
            txt_peak = "d=%.1fnm" %f.popt[0]
            pylab.text(22*f.ps, 0.75*ymax, txt_label, fontsize=22)
            pylab.text(16*f.ps, 0.50*ymax, txt_sigma, fontsize=8)
            pylab.text(16*f.ps, 0.35*ymax, txt_peak, fontsize=8)
            pylab.grid()

        pylab.savefig("plot-fig%i.png" %(2*isfig2+3*isfig3))

    # Fig.4 -- single images
    if "fig4" in sys.argv and "single" in sys.argv:

        sem_fig2 = [fig2[:4], fig2[4:8], fig2[8:]]
        sem_fig3 = [fig3[:4], fig3[4:8], fig3[8:]]
        sem_5k = [[loadsem(fn) for fn in fnames] for fnames in sem_fig2]
        sem_20k = [[loadsem(fn) for fn in fnames] for fnames in sem_fig3]
        for sems in sem_5k+sem_20k:
            for s in sems:
                s.radialdistribution()
                s.fitrdf()
        peaks_5k = [[s.rdfmax for s in sems] for sems in sem_5k]
        peaks_20k = [[s.rdfmax for s in sems] for sems in sem_20k]
        npconc_5k = [[s.npconc for s in sems] for sems in sem_5k]
        npconc_20k = [[s.npconc for s in sems] for sems in sem_20k]

        pylab.figure(figsize=(10,6))
        styles = ["o-", "^--", "s:"]
        kwargs = { "color": "black", "linewidth": 1.50 }

        pylab.subplot(121)
        labels = ["5k PEG, %i$\pm$%i/$\mu$m$^2$" %(avg(c),std(c, ddof=1)) for c in npconc_5k]
        for i in range(3):
            pylab.plot(times, peaks_5k[i], styles[i], label=labels[i], **kwargs)
        pylab.legend()
        pylab.grid()
        pylab.xlim([-1,50])
        pylab.ylim([20,55])

        pylab.subplot(122)
        labels = ["20k PEG, %i$\pm$%i/$\mu$m$^2$" %(avg(c),std(c, ddof=1)) for c in npconc_20k]
        for i in range(3):
            pylab.plot(times, peaks_5k[i], styles[i], label=labels[i], **kwargs)
        pylab.legend()
        pylab.grid()
        pylab.xlim([-1,50])
        pylab.ylim([23,53])

        pylab.savefig("plot-fig4-single.png")

    # Fig.4 -- with statistics
    if "fig4" in sys.argv and "statistics" in sys.argv:

        # Load and analyze all images
        # Note that we first organize the file names according to times
        fnames_5k_C3 = glob.glob("sem/5k/C3/*-200.jpg") + glob.glob("sem/5k/C3/*-500.jpg")
        fnames_5k_C6 = glob.glob("sem/5k/C6/*-200.jpg") + glob.glob("sem/5k/C6/*-500.jpg")
        fnames_5k_C8 = glob.glob("sem/5k/C8/*-200.jpg") + glob.glob("sem/5k/C8/*-500.jpg")
        fnames_20k_C1 = glob.glob("sem/20k/C1/*-200.jpg") + glob.glob("sem/20k/C1/*-500.jpg")
        fnames_20k_C3 = glob.glob("sem/20k/C3/*-200.jpg") + glob.glob("sem/20k/C3/*-500.jpg")
        fnames_20k_C5 = glob.glob("sem/20k/C5/*-200.jpg") + glob.glob("sem/20k/C5/*-500.jpg")
        fnames_5k_C3 = [[fn for fn in fnames_5k_C3 if "%i-Image" %t in fn] for t in times]
        fnames_5k_C6 = [[fn for fn in fnames_5k_C6 if "%i-Image" %t in fn] for t in times]
        fnames_5k_C8 = [[fn for fn in fnames_5k_C8 if "%i-Image" %t in fn] for t in times]
        fnames_20k_C1 = [[fn for fn in fnames_20k_C1 if "%i-Image" %t in fn] for t in times]
        fnames_20k_C3 = [[fn for fn in fnames_20k_C3 if "%i-Image" %t in fn] for t in times]
        fnames_20k_C5 = [[fn for fn in fnames_20k_C5 if "%i-Image" %t in fn] for t in times]
        sem_5k_C3 = [[loadsem(fn) for fn in fnames] for fnames in fnames_5k_C3]
        sem_5k_C6 = [[loadsem(fn) for fn in fnames] for fnames in fnames_5k_C6]
        sem_5k_C8 = [[loadsem(fn) for fn in fnames] for fnames in fnames_5k_C8]
        sem_20k_C1 = [[loadsem(fn) for fn in fnames] for fnames in fnames_20k_C1]
        sem_20k_C3 = [[loadsem(fn) for fn in fnames] for fnames in fnames_20k_C3]
        sem_20k_C5 = [[loadsem(fn) for fn in fnames] for fnames in fnames_20k_C5]
        sem_5k = [sem_5k_C3, sem_5k_C6, sem_5k_C8]
        sem_20k = [sem_20k_C1, sem_20k_C3, sem_20k_C5]
        for sems in sem_5k+sem_20k:
            for s in sems:
                for st in s:
                    st.radialdistribution()
                    st.fitrdf()
        peaks_5k = [[[st.rdfmax for st in s] for s in sems] for sems in sem_5k]
        peaks_20k = [[[st.rdfmax for st in s] for s in sems] for sems in sem_20k]
        npconc_5k = [[[st.npconc for st in s] for s in sems] for sems in sem_5k]
        npconc_20k = [[[st.npconc for st in s] for s in sems] for sems in sem_20k]

        pylab.figure(figsize=(10,6))
        styles = ["o-", "^--", "s:"]
        kwargs = { "color": "black", "linewidth": 1.50 }
        times = numpy.array(times)
        d = numpy.array([[0,0,0,0],[0.5,0.5,0.5,0.5],[1,1,1,1]])

        pylab.subplot(121)
        label_avg = [avg([c[0]+c[1]+c[2]]) for c in npconc_5k]
        label_std = [std([c[0]+c[1]+c[2]], ddof=1) for c in npconc_5k]
        labels = ["5k PEG, %i$\pm$%i/$\mathrm{\mu m^2}$" %(la,ls) for la,ls in zip(label_avg,label_std)]
        peaks_avg = [[avg(pp) for pp in p] for p in peaks_5k]
        err_std = array([[std(pp, ddof=1) for pp in p] for p in peaks_5k])
        err_fit = array([[avg([f.pcov[0,0] for f in pp]) for pp in p] for p in sem_5k])
        err_bin = array([[avg([f.dbin/2.0 for f in pp]) for pp in p] for p in sem_5k])
        err_total =  sqrt(err_std**2 + err_fit**2 + err_bin**2)
        for i in range(3):
            pylab.plot(times+d[i], peaks_avg[i], styles[i], label=labels[i], **kwargs)
            kwargs2 = { "color": "black" }
            for it,t in enumerate(times):
                pylab.errorbar(t+d[i][it], peaks_avg[i][it], yerr=err_total[i][it], **kwargs2)
        pylab.legend()
        pylab.grid()
        pylab.xlim([-1,50])
        pylab.ylim([20,55])

        pylab.subplot(122)
        label_avg = [avg([c[0]+c[1]+c[2]]) for c in npconc_20k]
        label_std = [std([c[0]+c[1]+c[2]], ddof=1) for c in npconc_20k]
        lformat = "20k PEG, %i$\pm$%i/$\mathrm{\mu m^2}$"
        labels = [lformat %(la,ls) for la,ls in zip(label_avg,label_std)]
        peaks_avg = [[avg(pp) for pp in p] for p in peaks_20k]
        err_std = array([[std(pp, ddof=1) for pp in p] for p in peaks_20k])
        err_fit = array([[avg([f.pcov[0,0] for f in pp]) for pp in p] for p in sem_20k])
        err_bin = array([[avg([f.dbin/2.0 for f in pp]) for pp in p] for p in sem_20k])
        err_total =  sqrt(err_std**2 + err_fit**2 + err_bin**2)
        for i in range(3):
            pylab.plot(times+d[i], peaks_avg[i], styles[i], label=labels[i], **kwargs)
            kwargs2 = { "color": "black" }
            for it,t in enumerate(times):
                pylab.errorbar(t+d[i][it], peaks_avg[i][it], yerr=err_total[i][it], **kwargs2)
        pylab.legend(loc="lower right")
        pylab.grid()
        pylab.xlim([-1,50])
        pylab.ylim([23,53])

        pylab.savefig("plot-fig4-statistics.png")
