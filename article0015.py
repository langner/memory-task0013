#Fig.2

import bz2
import glob
import itertools
import string
import sys

import numpy
import pylab

from scipy import misc

from systems import PhaseFrames

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
            # However, only add it for the most lower left image, since all bars are the same.
            if fi == ncols*nrows-1:
                scale_bar = pylab.Rectangle((f.img.shape[0]*1.1,cropy-20), width=500.0/f.ps, height=10, color="white")
                pylab.gca().add_patch(scale_bar)
                pylab.text(f.img.shape[0]*1.25, cropy-25, r"\textbf{500 nm}", fontsize=16, color="white", ha="center")

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
            txt_label = r"\textbf{%s}" %string.ascii_lowercase[fi]
            txt_peak = r"\textbf{d=%.0fnm}" %f.popt[0]
            pylab.text(3*f.ps, 0.75*ymax, txt_label, fontsize=13, ha="left")
            pylab.text(14*f.ps, 0.75*ymax, txt_peak, fontsize=13)
            pylab.grid()

        # Add times in text to side of image table.
        times = ["0h", "2h", "14h", "48h"]
        for i in range(nrows):
            ax_main.text(margin_x/2, (1-(1.0*i+0.5)/nrows)*effy, r"\textbf{%s}" %times[i], rotation=90, fontsize=20)

        # Add concentrations in text to top of image table (read from file).
        lines = open("article0015.results.txt").readlines()
        nline_avg = isfig2*1 + isfig3*6
        nline_std = isfig2*2 + isfig3*7
        sigma_avg = map(int,lines[nline_avg].split()) #[round(numpy.average([F[3*i+j].npconc for j in range(nrows)]), -1) for i in range(3)]
        sigma_std = map(int,lines[nline_std].split()) #[round(numpy.std([F[3*i+j].npconc for j in range(nrows)]), -1) for i in range(3)]
        sigma_dat = zip(sigma_avg,sigma_std)
        sigma_fmt = r"\boldmath $\sigma=%i\pm%i\,\mathrm{NP}/\mu\mathrm{m}^2$"
        for i in range(ncols):
            ax_main.text(margin_x+(i+0.25)*effx/ncols, effy+margin_y/4, sigma_fmt %sigma_dat[i], fontsize=18, weight="bold")

        pylab.savefig("article0015-fig%i.png" %(2*isfig2+3*isfig3))

    # Fig.4 -- single images
    if "fig4" in sys.argv and "single" in sys.argv:

        sem_fig2 = [fig2[:4], fig2[4:8], fig2[8:]]
        sem_fig3 = [fig3[:4], fig3[4:8], fig3[8:]]
        sem_5k = [[loadsem("exp/sem/"+fn) for fn in fnames] for fnames in sem_fig2]
        sem_20k = [[loadsem("exp/sem/"+fn) for fn in fnames] for fnames in sem_fig3]

        for sems in sem_5k+sem_20k:
            for s in sems:
                print "Analyzing %s..." %s.name
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
            pylab.plot(times, peaks_20k[i], styles[i], label=labels[i], **kwargs)
        pylab.legend()
        pylab.grid()
        pylab.xlim([-1,50])
        pylab.ylim([20,55])

        # Labels are common, so create them on the main figure axes.
        ax.text(0.48, 0.02, "time [h]")
        ax.text(0.05, 0.70, "nearest neighbor distance ($d$) [nm]", rotation=90)

        pylab.savefig("article0015-fig4-single.png")

    # Fig.4 -- with statistics
    if "fig4" in sys.argv and "single" not in sys.argv:

        # Load and analyze all images.
        # Note that we first organize the file names according to times.
        def find_images(*pats):
            return list(itertools.chain(*[glob.glob("exp/sem/"+p) for p in pats]))
        fnames_5k_C3 = find_images("5k/C3/*-200.jpg", "5k/C3/*-500.jpg", "5k/C3/*-1000.jpg")
        fnames_5k_C6 = find_images("5k/C6/*-200.jpg", "5k/C6/*-500.jpg", "5k/C6/*-1000.jpg")
        fnames_5k_C8 = find_images("5k/C8/*-200.jpg", "5k/C8/*-500.jpg", "5k/C8/*-1000.jpg")
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
                    print "Analyzing %s..." %st.outname
                    st.radialdistribution()
                    st.fitrdf()

        npconc_5k = [[[st.npconc for st in s] for s in sems] for sems in sem_5k]
        npconc_20k = [[[st.npconc for st in s] for s in sems] for sems in sem_20k]

        ax = pylab.figure(figsize=(10,6))
        styles = ["o-", "^--", "s:"]
        kwargs = { "color": "black", "linewidth": 1.50 }
        times = numpy.array(times)
        d = numpy.array([[0,0,0,0],[0.5,0.5,0.5,0.5],[1,1,1,1]])

        print "Error for 5k case..."

        pylab.subplot(121)
        label_avg_5k = numpy.round([avg([c[0]+c[1]+c[2]]) for c in npconc_5k], -1)
        label_std_5k = numpy.round([std([c[0]+c[1]+c[2]], ddof=1) for c in npconc_5k], -1)
        labels_5k = ["5k PEG, %i$\pm$%i/$\mathrm{\mu m^2}$" %(la,ls) for la,ls in zip(label_avg_5k,label_std_5k)]
        weights_5k = [[[1.0*st.npcount/st.ps for st in s] for s in sems] for sems in sem_5k]
        peaks_avg_5k = array([[avg([st.rdfmax for st in s], weights=weights_5k[isems][si]) for si,s in enumerate(sems)] for isems,sems in enumerate(sem_5k)])
        # add weights to std!
        # binning does not contribute, since we are taking the fit result
        err_std_5k = array([[std([st.rdfmax for st in s], ddof=1) for s in sems] for sems in sem_5k])
        err_fit_5k = array([[avg([f.pcov[0,0] for f in s], weights=weights_5k[isems][si]) for si,s in enumerate(sems)] for isems,sems in enumerate(sem_5k)])
        err_bin_5k = array([[avg([f.dbin/2.0 for f in s], weights=weights_5k[isems][si]) for si,s in enumerate(sems)] for isems,sems in enumerate(sem_5k)])
        err_total_5k =  sqrt(err_std_5k**2 + err_fit_5k**2)
        print "Standard deviations:\n", err_std_5k
        print "Errors from fitting\n", err_fit_5k
        print "Total errors\n", err_total_5k
        for i in range(3):
            pylab.plot(times+d[i], peaks_avg_5k[i], styles[i], label=labels_5k[i], **kwargs)
            kwargs2 = { "color": "black" }
            for it,t in enumerate(times):
                pylab.errorbar(t+d[i][it], peaks_avg_5k[i][it], yerr=err_total_5k[i][it], **kwargs2)
        pylab.legend()
        pylab.grid()
        pylab.xlim([-1,50])
        pylab.ylim([20,55])

        pylab.subplot(122)
        label_avg_20k = numpy.round([avg([c[0]+c[1]+c[2]]) for c in npconc_20k], -1)
        label_std_20k = numpy.round([std([c[0]+c[1]+c[2]], ddof=1) for c in npconc_20k], -1)
        labels_20k = ["20k PEG, %i$\pm$%i/$\mathrm{\mu m^2}$" %(la,ls) for la,ls in zip(label_avg_20k,label_std_20k)]
        weights_20k = [[[1.0*st.npcount/st.ps for st in s] for s in sems] for sems in sem_20k]
        peaks_avg_20k = array([[avg([st.rdfmax for st in s], weights=weights_20k[isems][si]) for si,s in enumerate(sems)] for isems,sems in enumerate(sem_20k)])
        # add weights to std!
        # binning does no contribute, since we are taking the fit result
        err_std_20k = array([[std([st.rdfmax for st in s], ddof=1) for s in sems] for sems in sem_20k])
        err_fit_20k = array([[avg([f.pcov[0,0] for f in s], weights=weights_20k[isems][si]) for si,s in enumerate(sems)] for isems,sems in enumerate(sem_20k)])
        err_bin_20k = array([[avg([f.dbin/2.0 for f in s], weights=weights_20k[isems][si]) for si,s in enumerate(sems)] for isems,sems in enumerate(sem_20k)])
        err_total_20k =  sqrt(err_std_20k**2 + err_fit_20k**2)
        print "Error for 20k case..."
        print "Standard deviations:\n", err_std_20k
        print "Errors from fitting\n", err_fit_20k
        print "Total errors\n", err_total_20k
        for i in range(3):
            pylab.plot(times+d[i], peaks_avg_20k[i], styles[i], label=labels_20k[i], **kwargs)
            kwargs2 = { "color": "black" }
            for it,t in enumerate(times):
                pylab.errorbar(t+d[i][it], peaks_avg_20k[i][it], yerr=err_total_20k[i][it], **kwargs2)
        pylab.legend(loc="lower right")
        pylab.grid()
        pylab.xlim([-1,50])
        pylab.ylim([20,55])

        # Save some results for reuse (notably concentrations in Figs 2 and 3).
        f = open("article0015.results.txt", "w")
        print >>f, "5k"
        print >>f, " ".join(["%i" %n for n in label_avg_5k])
        print >>f, " ".join(["%i" %n for n in label_std_5k])
        print >>f, " ".join(["%.1f" %x for x in itertools.chain(*peaks_avg_5k)])
        print >>f, " ".join(["%.1f" %x for x in itertools.chain(*err_total_5k)])
        print >>f, "20k"
        print >>f, " ".join(["%i" %n for n in label_avg_20k])
        print >>f, " ".join(["%i" %n for n in label_std_20k])
        print >>f, " ".join(["%.1f" %x for x in itertools.chain(*peaks_avg_20k)])
        print >>f, " ".join(["%.1f" %x for x in itertools.chain(*err_total_20k)])
        f.close()

        # Labels are common, so create them on the main figure axes.
        ax.text(0.48, 0.02, "time [h]")
        ax.text(0.05, 0.70, "nearest neighbor distance ($d$) [nm]", rotation=90)

        pylab.savefig("article0015-fig4.png")

    if "fig5" in sys.argv:

        pylab.rc('text', usetex=True)

        from simulation import loadpath

        ncols, nrows = 3, 4

        chosen = []
        for i in range(ncols):
            chosen.append([None, None, None, None])
        chosen[0][0] = "phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop700_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.frame0001.jpg"
        chosen[0][1] = "phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop700_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.frame0045.jpg"
        chosen[0][2] = "phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop700_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.frame0321.jpg"
        chosen[0][3] = "phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop700_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.frame1101.jpg"
        chosen[1][0] = "phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop400_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.frame0001.jpg"
        chosen[1][1] = "phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop400_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.frame0045.jpg"
        chosen[1][2] = "phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop400_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.frame0321.jpg"
        chosen[1][3] = "phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop400_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.frame1101.jpg"
        chosen[2][0] = "phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop200_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.frame0001.jpg"
        chosen[2][1] = "phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop200_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.frame0045.jpg"
        chosen[2][2] = "phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop200_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.frame0321.jpg"
        chosen[2][3] = "phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop200_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.frame1101.jpg"
        imgs = [[misc.imread(fn)[25:-24,25:-22] for fn in ch] for ch in chosen]
        rdfs = [[numpy.load(bz2.BZ2File(fn.replace(fn[-14:],".hist-radial.npy.bz2"))) for fn in ch] for ch in chosen]

        margin_x, margin_y = 0.04, 0.03
        effx, effy = 1.0 - margin_x, 1.0 - margin_y
        img_width, img_height = 3.25, 3.25
        fig_width = ncols*img_width*(1.0+margin_x)
        fig_height = nrows*img_height*(1.0+margin_x)
        ax_main = pylab.figure(figsize=(fig_width,fig_height))

        for fi in range(ncols*nrows):

            x = margin_x + (effx/ncols)*(fi//nrows)
            y = effy*(nrows-1)/nrows - effy*(fi%nrows)/nrows + 0.01

            # Show image, without any ticks.
            a = pylab.axes([x, y, effx/ncols-0.01, effy/nrows-0.005])
            pylab.setp(a, xticks=[], yticks=[])
            pylab.imshow(imgs[fi//nrows][fi%nrows], cmap='gray')

        # Add times in text to side of image table.
        times = ["0", "45", "320", "1100"]
        for i in range(nrows):
            ax_main.text(margin_x/3, (1-(1.0*i+0.5)/nrows)*effy, times[i], rotation=90, fontsize=20, weight="bold")

        # Add concentrations in text to top of image table.
        conc = [700, 400, 200]
        for i in range(ncols):
            ax_main.text(margin_x+(i+0.35)*effx/ncols, effy+margin_y/4, "%i NPs" %conc[i], fontsize=20, weight="bold")

        pylab.savefig("article0015-fig5.png")

    if "fig6" in sys.argv:

        chosen_exp = "exp/sem/5k/C6/48-Image5-500.jpg"
        chosem_sim = "phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop400_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.hist-radial.npy.bz2"

        sem = loadsem(chosen_exp)
        sem.fitrdf()

        radial = 1.0 * numpy.load(bz2.BZ2File(chosem_sim))[:,0,:]

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

        pylab.plot(sem.X, sem.Y, color="black", label="experimental", linewidth=1.5)
        pylab.plot(sem.X, sem.rdfmodel(sem.X), color="black", label="exp. fit", linewidth=1.5)

        pylab.plot(xrange, numpy.sum(radial[iplot-25:iplot], axis=0)/25, color="gray", label="simulation", linewidth=1.5)

        pylab.xlim([5,150])
        pylab.ylim([0,9])

        pylab.xlabel("NP-NP distance [nm]")
        pylab.ylabel("g(r)")

        #pylab.legend(loc=[0.20,0.70])
        pylab.grid()

        ax = pylab.axes([0.45, 0.45, 0.4, 0.4], axisbg='lightgray')
        pylab.plot(sem.X, sem.Y, color="black")
        pylab.xticks([100, 500])
        pylab.yticks([1, 3 ,5])
        pylab.grid()

        pylab.savefig("article0015-fig6.png")

    if "fig7" in sys.argv:

        pylab.rc('text', usetex=True)

        chosen_exp = ["exp/sem/"+fn for fn in fig2]
        chosen_sim = ["phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop700_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.hist-radial.npy.bz2",
"phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop400_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.hist-radial.npy.bz2",
"phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop200_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.hist-radial.npy.bz2",
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

        ns = [700, 400, 200]
        npairs = [n*(n-1)/2 for n in ns]
        strips = 1.0 * correction(size[0],size[1],xrange) * 2 * numpy.pi * xrange * dx
        factors = numpy.array([1.0 * size[0] * size[1] / strips / np for np in npairs])

        for ni,n in enumerate(ns):
            for fi,f in enumerate(factors[ni]):
                radials[ni][:,fi] *= f

        pylab.subplot(121)

        exp_d = [s.rdfmax for s in sems]

        D = [xrange[numpy.argmax(radial[:,:250], axis=1)] for radial in radials]
        for di,d in enumerate(D):
            newd = []
            for fi,f in enumerate(PhaseFrames[11][0]):
                i = PhaseFrames[11][1]
                newd.append(numpy.average(d[fi*i:(fi+1)*i]))
                D[di] = newd

        X = numpy.array(PhaseFrames[11][0])
        X *= 48.0 / 200.0

        styles = ["ko-", "ks-", "k^-"]
        labels = ["a-d", "e-h", "i-l"]
        for di,d in enumerate(D):
            d = numpy.array(d)
            d *= exp_d[3]/d[-1]
            pylab.plot(X, d, styles[di], label= "sim. %i NPs" %(ns[di]))

        pylab.plot([0,2,14,48], exp_d[:4], "ko-.", mfc='None', label="exp. 760 NP/$\mu$m$^2$")
        pylab.plot([0,2,14,48], exp_d[4:8], "ks-.", mfc='None', label="exp. 510 NP/$\mu$m$^2$")
        pylab.plot([0,2,14,48], exp_d[8:12], "k^-.", mfc='None', label="exp. 300 NP/$\mu$m$^2$")

        pylab.xlim([0, 50])
        pylab.ylim([23, 47])

        pylab.xlabel("time [h]")
        pylab.ylabel("interparticle distance $d$ [nm]")
        pylab.grid()
        pylab.legend(prop={'size':10})

        pylab.subplot(122)

        exp_g = [s.popt[-1] for s in sems]
        pylab.plot([0,2,14,48], exp_g[:4], "ko-.", mfc='None')
        pylab.plot([0,2,14,48], exp_g[4:8], "ks-.", mfc='None')
        pylab.plot([0,2,14,48], exp_g[8:12], "k^-.", mfc='None')

        G = [numpy.max(radial[:,:250], axis=1) for radial in radials]
        for gi,g in enumerate(G):
            newg = []
            for fi,f in enumerate(PhaseFrames[11][0]):
                i = PhaseFrames[11][1]
                newg.append(numpy.average(g[fi*i:(fi+1)*i]))
                G[gi] = newg

        for gi,g in enumerate(G):
            g = numpy.array(g)
            pylab.plot(X, g, styles[gi])

        pylab.xlim([0, 50])
        pylab.ylim([0, 40])

        pylab.xlabel("time [h]")
        pylab.ylabel("first $g(r)$ peak height")

        pylab.grid()

        pylab.savefig("article0015-fig7.png")
