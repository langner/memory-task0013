#Fig.2

import bz2
import glob
import itertools
import string
import sys

import numpy
import pylab

sys.path.append("exp/")
from common import SEMAnalysis, getscalefromname
from sem import cropy, scalebarstart


fig1 = [
        "5k/C6/14-Image35-500.jpg",
        "5k/C6/14-Image37-200.jpg",
        "20k/C3/14-Image25-500.jpg",
        "20k/C3/14-Image24-200.jpg",
]

fig2_conc = [3, 6, 8]
fig2 = [
        "5k/C3/0-Image4-500.jpg",
        "5k/C3/2-Image5-500.jpg",
        "5k/C3/14-Image6-500.jpg",
        "5k/C3/48-Image26-500.jpg",
        "5k/C6/0-Image21-500.jpg",
        "5k/C6/2-Image1-500.jpg",
        "5k/C6/14-Image35-500.jpg",
        "5k/C6/48-Image5-500.jpg",
        "5k/C8/0-Image27-500.jpg",
        "5k/C8/2-Image29-500.jpg",
        "5k/C8/14-Image56-500.jpg",
        "5k/C8/48-Image17-500.jpg",
]

fig3_conc = [1, 3, 5]
fig3 = [
        "20k/C1/0-Image36-500.jpg",
        "20k/C1/2-Image9-500.jpg",
        "20k/C1/14-Image22-500.jpg",
        "20k/C1/48-Image43-500.jpg",
        "20k/C3/0-Image23-500.jpg",
        "20k/C3/2-Image13-500.jpg",
        "20k/C3/14-Image24-500.jpg",
        "20k/C3/48-Image10-500.jpg",
        "20k/C5/0-Image24-500.jpg",
        "20k/C5/2-Image17-500.jpg",
        "20k/C5/14-Image33-500.jpg",
        "20k/C5/48-Image27-500.jpg",
]

fig8 = [
        "5k/C1/2-Image2-200.jpg",
        "5k/C1/48-Image5-200.jpg",
]

# Times at which SEMS images were registered.
times = [0, 2, 14, 48]

# Aliases.
sbs, gsf = scalebarstart, getscalefromname
loadsem = lambda fn: SEMAnalysis(fpath=fn, cropy=(0,430), scalebarstart=sbs[gsf(fn)])
array = numpy.array
sqrt = numpy.sqrt
avg = numpy.average
std = numpy.std


if __name__ == "__main__":

    # Fig.2 and Fig.3
    isfig2 = "fig2" in sys.argv
    isfig3 = "fig3" in sys.argv
    if isfig2 or isfig3:

        pylab.rc('text', usetex=True)

        # Load and analyze all images
        to_plot = fig2*isfig2 + fig3*isfig3
        F = [loadsem("exp/sem/"+fn) for fn in to_plot]
        for f in F:
            f.radialdistribution()
            f.fitrdf()

        # Notice how things are parametrized by number of columns and rows.
        # We have decided to remove the bottom row (48h) from Figure 3.
        ncols, nrows = 3, 4
        if isfig3:
            nrows = 3
            F = [f for i,f in enumerate(F) if i%4 != 3]
        margin_x, margin_y = 0.03, 0.03
        effx, effy = 1.0 - margin_x, 1.0 - margin_y
        img_width, img_height = 4.8, 3.25
        fig_width = ncols*img_width*(1.0+margin_x)
        fig_height = nrows*img_height*(1.0+margin_x)
        ax_main = pylab.figure(figsize=(fig_width,fig_height))

        # Now fill in the table with images and insets.
        for fi,f in enumerate(F):

            # Compute the starting coordinates for this image.
            # Note that we want top-down order, whereas matplotlib does left-right.
            x = margin_x + (effx/ncols)*(fi//nrows) + 0.01
            y = effy*(nrows-1)/nrows - effy*(fi%nrows)/nrows + 0.01

            # Show image, without any ticks.
            a = pylab.axes([x, y, effx/ncols-0.02, effy/nrows-0.005])
            pylab.setp(a, xticks=[], yticks=[])
            pylab.imshow(f.img, cmap='gray')

            # Add a scale bar in the lower left of the image (500nm base on pixel size).
            scale_bar = pylab.Rectangle((f.img.shape[0]*0.03,cropy-20), width=500.0/f.ps, height=10, color="white")
            pylab.gca().add_patch(scale_bar)
            pylab.text(f.img.shape[0]*0.05, cropy-25, "500 nm", fontsize=11, color="white", weight="bold")

            # Determine scaling (same in all plots)
            ymax = max(3.0, max([max(ff.rdf) for ff in F if hasattr(ff, 'rdf')]))

            # Plot the RDF in the upper left corner
            # Plot the fitted function over the experimental data
            a = pylab.axes([x+0.01, y+effy/nrows*0.6-0.01, effx/ncols/3-0.01, effy/nrows*0.4-0.01])
            if fi%nrows == 0:
                peak_x = f.rdfmax
            pylab.setp(a, xticks=[peak_x], xticklabels=[""], yticks=[], xlim=[0, 30*f.ps], ylim=[0, ymax])
            pylab.plot(f.X, f.Y, label="%ih" %f.time, color="black", linewidth=1.00)
            pylab.plot(f.X, f.rdfmodel(f.X), color="black", linewidth=1.75)
            txt_label = string.ascii_lowercase[fi]
            txt_peak = r"\boldmath d=%.1fnm" %f.popt[0]
            pylab.text(28*f.ps, 0.75*ymax, txt_label, fontsize=20, weight="bold", ha="right")
            pylab.text(14*f.ps, 0.55*ymax, txt_peak, fontsize=12, weight="bold")
            pylab.grid()

        # Add times in text to side of image table.
        times = ["0h", "2h", "14h", "48h"]
        for i in range(nrows):
            ax_main.text(margin_x/2, (1-(1.0*i+0.5)/nrows)*effy, times[i], rotation=90, fontsize=20, weight="bold")

        # Add concentrations in text to top of image table.
        sigma_avg = [round(numpy.average([F[3*i+j].npconc for j in range(nrows)]), -1) for i in range(3)]
        sigma_std = [round(numpy.std([F[3*i+j].npconc for j in range(nrows)]), -1) for i in range(3)]
        sigma_dat = zip(sigma_avg,sigma_std)
        sigma_fmt = r"\boldmath $\sigma$ = %i$\pm$%i NP/$\mu$m$^2$"
        for i in range(ncols):
            ax_main.text(margin_x+(i+0.25)*effx/ncols, effy+margin_y/4, sigma_fmt %sigma_dat[i], fontsize=20, weight="bold")

        pylab.savefig("article0015-fig%i.png" %(2*isfig2+3*isfig3))

    # Fig.4 -- single images
    if "fig4" in sys.argv and "single" in sys.argv:

        sem_fig2 = [fig2[:4], fig2[4:8], fig2[8:]]
        sem_fig3 = [fig3[:4], fig3[4:8], fig3[8:]]
        sem_5k = [[loadsem("exp/sem/"+fn) for fn in fnames] for fnames in sem_fig2]
        sem_20k = [[loadsem("exp/sem/"+fn) for fn in fnames] for fnames in sem_fig3]
        for sems in sem_5k+sem_20k:
            for s in sems:
                s.radialdistribution()
                s.fitrdf()
        peaks_5k = [[s.rdfmax for s in sems] for sems in sem_5k]
        peaks_20k = [[s.rdfmax for s in sems] for sems in sem_20k]
        npconc_5k = [[s.npconc for s in sems] for sems in sem_5k]
        npconc_20k = [[s.npconc for s in sems] for sems in sem_20k]

        ax = pylab.figure(figsize=(10,6))
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

        # Labels are common, so create them on the main figure axes.
        ax.text(0.48, 0.02, "time [h]")
        ax.text(0.05, 0.70, "nearest neighbor distance ($d$) [nm]", rotation=90)

        pylab.savefig("article0015-fig4-single.png")

    # Fig.4 -- with statistics
    if "fig4" in sys.argv and "single" not in sys.argv:

        # Load and analyze all images
        # Note that we first organize the file names according to times
        def find_images(*pats):
            return list(itertools.chain(*[glob.glob("exp/sem/"+p) for p in pats]))
        fnames_5k_C3 = find_images("5k/C3/*-200.jpg", "5k/C3/*-500.jpg", "5k/C3/*-1000.jpg")
        fnames_5k_C6 = find_images("5k/C6/*-200.jpg", "5k/C6/*-500.jpg", "5k/C6/*-1000.jpg")
        fnames_5k_C8 = find_images("k/C8/*-200.jpg", "5k/C8/*-500.jpg", "5k/C8/*-1000.jpg")
        fnames_20k_C1 = find_images("20k/C1/*-200.jpg", "20k/C1/*-500.jpg", "20k/C1/*-1000.jpg")
        fnames_20k_C3 = find_images("20k/C3/*-200.jpg", "20k/C3/*-500.jpg", "20k/C3/*-1000.jpg")
        fnames_20k_C5 = find_images("20k/C5/*-200.jpg", "20k/C5/*-500.jpg", "20k/C5/*-1000.jpg")
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
        npconc_5k = [[[st.npconc for st in s] for s in sems] for sems in sem_5k]
        npconc_20k = [[[st.npconc for st in s] for s in sems] for sems in sem_20k]

        ax = pylab.figure(figsize=(10,6))
        styles = ["o-", "^--", "s:"]
        kwargs = { "color": "black", "linewidth": 1.50 }
        times = numpy.array(times)
        d = numpy.array([[0,0,0,0],[0.5,0.5,0.5,0.5],[1,1,1,1]])

        pylab.subplot(121)
        label_avg = numpy.round([avg([c[0]+c[1]+c[2]]) for c in npconc_5k])
        label_std = numpy.round([std([c[0]+c[1]+c[2]], ddof=1) for c in npconc_5k])
        labels = ["5k PEG, %i$\pm$%i/$\mathrm{\mu m^2}$" %(la,ls) for la,ls in zip(label_avg,label_std)]
        peaks_avg = array([[avg([st.rdfmax for st in s], weights=[st.scale for st in s]) for s in sems] for sems in sem_5k])
        # add weights to std!
        # binning does no contribute, since we are taking the fit result
        err_std = array([[std([st.rdfmax for st in s], ddof=1) for s in sems] for sems in sem_5k])
        err_fit = array([[avg([f.pcov[0,0] for f in pp], weights=[f.scale for f in pp]) for pp in p] for p in sem_5k])
        err_bin = array([[avg([f.dbin/2.0 for f in pp], weights=[f.scale for f in pp]) for pp in p] for p in sem_5k])
        err_total =  sqrt(err_std**2 + err_fit**2)
        print "Error for 5k case..."
        print "Standard deviations:\n", err_std
        print "Errors from fitting\n", err_fit
        print "Total errors\n", err_total
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
        peaks_avg = array([[avg([st.rdfmax for st in s], weights=[st.scale for st in s]) for s in sems] for sems in sem_20k])
        # add weights to std!
        # binning does no contribute, since we are taking the fit result
        err_std = array([[std([st.rdfmax for st in s], ddof=1) for s in sems] for sems in sem_20k])
        err_fit = array([[avg([f.pcov[0,0] for f in pp], weights=[f.scale for f in pp]) for pp in p] for p in sem_20k])
        err_bin = array([[avg([f.dbin/2.0 for f in pp], weights=[f.scale for f in pp]) for pp in p] for p in sem_20k])
        err_total =  sqrt(err_std**2 + err_fit**2)
        print "Error for 20k case..."
        print "Standard deviations:\n", err_std
        print "Errors from fitting\n", err_fit
        print "Total errors\n", err_total
        for i in range(3):
            pylab.plot(times+d[i], peaks_avg[i], styles[i], label=labels[i], **kwargs)
            kwargs2 = { "color": "black" }
            for it,t in enumerate(times):
                pylab.errorbar(t+d[i][it], peaks_avg[i][it], yerr=err_total[i][it], **kwargs2)
        pylab.legend(loc="lower right")
        pylab.grid()
        pylab.xlim([-1,50])
        pylab.ylim([23,53])

        # Labels are common, so create them on the main figure axes.
        ax.text(0.48, 0.02, "time [h]")
        ax.text(0.05, 0.70, "nearest neighbor distance ($d$) [nm]", rotation=90)

        pylab.savefig("article0015-fig4.png")

    if "fig5" in sys.argv:

        pass

    if "fig6" in sys.argv:

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

    if "fig7" in sys.argv:

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
