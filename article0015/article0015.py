"""Script that generates plots for articl00013."""


import bz2
import glob
import itertools
import string
import sys

import numpy
import pylab

from scipy import misc

sys.path.append("../")
from systems import phases_frames

sys.path.append("../exp/")
from common import normalize_rdf, rdfcorrection
from sem import cropy, loadsem


# Some aliases.
array = numpy.array
average = numpy.average

# These are the images that are used for figures.
# For Fig.2 and Fig.3, the three different are stacked sequentially.
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
# This is the same for all series/concentrations.
times = [0, 2, 14, 48]


if __name__ == "__main__":

    # Fig.2 and Fig.3
    isfig2 = "fig2" in sys.argv
    isfig3 = "fig3" in sys.argv
    if isfig2 or isfig3:

        # Interpret all text as TeX.
        pylab.rc('text', usetex=True)

        # Load and analyze all images.
        to_plot = fig2*isfig2 + fig3*isfig3
        F = [loadsem("../exp/sem/"+fn) for fn in to_plot]
        print "Analyzing images..."
        for f in F:
            f.radialdistribution()
            f.fitrdf()

        # Most things are parametrized by number of columns and rows.
        ncols, nrows = 3, 4
        if isfig3:
            nrows = 3

        # We have decided to remove the bottom row (48h) from Figure 3.
        if isfig3:
            F = [f for i,f in enumerate(F) if i%4 != 3]

        # Margins around plots.
        margin_x, margin_y = 0.03, 0.03

        # Effective relative size of entire plots.
        effx, effy = 1.0 - margin_x, 1.0 - margin_y

        # Width and height of images in these plots.
        img_width, img_height = 4.8, 3.25

        # Width and height of the entire figure.
        fig_width = ncols*img_width*(1.0+margin_x)
        fig_height = nrows*img_height*(1.0+margin_x)

        # The main axis of the entire figure.
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
        nline_avg = isfig2*1 + isfig3*8
        nline_std = isfig2*2 + isfig3*9
        sigma_avg = map(int,lines[nline_avg].split())
        sigma_std = map(int,lines[nline_std].split())
        sigma_dat = zip(sigma_avg,sigma_std)
        sigma_fmt = r"\boldmath $\sigma=%i\pm%i\,\mathrm{NP}/\mu\mathrm{m}^2$"
        for i in range(ncols):
            ax_main.text(margin_x+(i+0.25)*effx/ncols, effy+margin_y/4, sigma_fmt %sigma_dat[i], fontsize=18, weight="bold")

        # Save the figure.
        pylab.savefig("article0015-fig%i.png" %(2*isfig2+3*isfig3))

    # Fig.4 -- single images
    if "fig4" in sys.argv and "single" in sys.argv:

        # Divide the images chosen for Figs.2-3 into sublists.
        sem_fig2 = [fig2[:4], fig2[4:8], fig2[8:]]
        sem_fig3 = [fig3[:4], fig3[4:8], fig3[8:]]

        # Load the file names into SEM objects.
        sem_5k = [[loadsem("../exp/sem/"+fn) for fn in fnames] for fnames in sem_fig2]
        sem_20k = [[loadsem("../exp/sem/"+fn) for fn in fnames] for fnames in sem_fig3]

        # Analyze all images (redundant, but I don't feel like saving the results).
        print "Analyzing images..."
        for sems in sem_5k+sem_20k:
            for s in sems:
                s.radialdistribution()
                s.fitrdf()

        # We will use peak positions and concentrations repeatedly.
        peaks_5k = [[s.rdfmax for s in sems] for sems in sem_5k]
        peaks_20k = [[s.rdfmax for s in sems] for sems in sem_20k]
        peaks = [peaks_5k, peaks_20k]
        npconc_5k = [[s.npconc for s in sems] for sems in sem_5k]
        npconc_20k = [[s.npconc for s in sems] for sems in sem_20k]
        npconc = [npconc_5k, npconc_20k]

        # Axes of the root figure.
        ax = pylab.figure(figsize=(10,6))

        # Styles, colors and linewidths we will use across figures.
        styles = ["o-", "^--", "s:"]
        kwargs = { "color": "black", "linewidth": 1.50 }

        # Text that will be in all legends.
        legend_format = r"%ik PEG, %i$\pm$%i/$\mu$m$^2$"

        # Plot the left and right figures now.
        nligand = [5, 20]
        for fi in range(2):
            pylab.subplot(120+fi+1)
            labels = [legend_format %(nligand[fi], numpy.average(c), numpy.std(c, ddof=1)) for c in npconc[fi]]
            for i in range(3):
                pylab.plot(times, peaks[fi][i], styles[i], label=labels[i], **kwargs)
            pylab.xlim([-1,50])
            pylab.ylim([20,55])
            pylab.legend()
            pylab.grid()

        # Labels are common, so create them on the main figure axes.
        ax.text(0.48, 0.02, "time [h]")
        ax.text(0.05, 0.70, "nearest neighbor distance ($d$) [nm]", rotation=90)

        pylab.savefig("article0015-fig4-single.png")

    # Fig.4 -- with statistics
    if "fig4" in sys.argv and "single" not in sys.argv:

        # Load and analyze all images.
        # Note that we first organize the file names according to times.
        def find_images(*pats):
            return list(itertools.chain(*[glob.glob("../exp/sem/"+p) for p in pats]))
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
        sems = [sem_5k, sem_20k]

        # Analyze all images (redundant, but anyway).
        print "Analyzing images..."
        for sem in itertools.chain(*sems):
            for st in itertools.chain(*sem):
                st.radialdistribution()
                st.fitrdf()

        # Concentrations of images.
        npconc_5k = [[[st.npconc for st in s] for s in sem] for sem in sem_5k]
        npconc_20k = [[[st.npconc for st in s] for s in sem] for sem in sem_20k]
        npconc = [npconc_5k, npconc_20k]

        # Axes of the root figure.
        ax = pylab.figure(figsize=(10,6))

        # Styles, color, linewidth and legend position to be used.
        styles = ["o-", "^--", "s:"]
        kwargs = { "color": "black", "linewidth": 1.50 }
        lpos = ['upper right', 'lower right']

        # Shift some points on the plots to the right, so the error bars do not overlap.
        # The times need to be an array to be nicely used later on.
        times = array(times)
        d = array([[0,0,0,0],[0.5,0.5,0.5,0.5],[1,1,1,1]])

        # Save some results for reuse after this, mostly for using the average concentrations in Figs 2 and 3.
        f = open("article0015.results.txt", "w")

        # Do the plotting for 5k and 20k with the same code.
        nligand = [5, 20]
        for ni,n in enumerate(nligand):

            print "Plotting for %ik case..." %n
            pylab.subplot(120+ni+1)

            # Round the average concentration and standard deviation to tens.
            label_avg = numpy.round([numpy.average([c[0]+c[1]+c[2]]) for c in npconc[ni]], -1)
            label_std = numpy.round([numpy.std([c[0]+c[1]+c[2]], ddof=1) for c in npconc[ni]], -1)

            # Labels for the legends.
            fmt = "%ik PEG, %i$\pm$%i/$\mathrm{\mu m^2}$"
            labels = [fmt %(n,la,ls) for la,ls in zip(label_avg,label_std)]

            # Calculate weights for averages -- use population divided by pixel size.
            weights = [[[1.0*ss.npcount/ss.ps for ss in s] for s in sem] for sem in sems[ni]]

            # Get the average peak positions (pp) and peak hieghts (ph), using above weights.        
            pp_avg = [[average([ss.rdfmax for ss in s], weights=w) for w,s in zip(wgt,sem)] for wgt,sem in zip(weights,sems[ni])]
            ph_avg = [[average([ss.popt[-1] for ss in s], weights=w) for w,s in zip(wgt,sem)] for wgt,sem in zip(weights,sems[ni])]

            # The total variance comprises the variance for the images available for this sample
            #   and the variance of the fitting parameter. 
            pp_err_std = [[numpy.std([ss.rdfmax for ss in s], ddof=1) for s in sem] for sem in sems[ni]]
            pp_err_fit = [[average([ss.pcov[0,0] for ss in s], weights=w) for w,s in zip(wgt,sem)] for wgt,sem in zip(weights,sems[ni])]
            pp_err_total = numpy.sqrt(array(pp_err_std)**2 + array(pp_err_fit))
            ph_err_std = [[numpy.std([ss.popt[-1] for ss in s], ddof=1) for s in sem] for sem in sems[ni]]
            ph_err_fit = [[average([ss.pcov[-1,-1] for ss in s], weights=w) for w,s in zip(wgt,sem)] for wgt,sem in zip(weights,sems[ni])]
            ph_err_total = numpy.sqrt(array(ph_err_std)**2 + array(ph_err_fit))
            err_kwargs = { "color": "black" }

            # Now plot the trends and error bars.
            for i in range(3):
                pylab.plot(times+d[i], pp_avg[i], styles[i], label=labels[i], **kwargs)
                for ti,t in enumerate(times):
                    pylab.errorbar(t+d[i][ti], pp_avg[i][ti], yerr=pp_err_total[i][ti], **err_kwargs)

            pylab.xlim([-1,50])
            pylab.ylim([20,55])
            pylab.legend(loc=lpos[ni])
            pylab.grid()

            # Write the concenctrations and positions to the results file.
            print >>f, "%ik" %n
            print >>f, " ".join(["%i" %n for n in label_avg])
            print >>f, " ".join(["%i" %n for n in label_std])
            print >>f, " ".join(["%.1f" %x for x in itertools.chain(*pp_avg)])
            print >>f, " ".join(["%.1f" %x for x in itertools.chain(*pp_err_total)])
            print >>f, " ".join(["%.1f" %x for x in itertools.chain(*ph_avg)])
            print >>f, " ".join(["%.1f" %x for x in itertools.chain(*ph_err_total)])

        # Close the file with results.
        f.close()

        # Labels are common, so create them on the main figure axes.
        ax.text(0.48, 0.02, "time [h]")
        ax.text(0.05, 0.70, "nearest neighbor distance ($d$) [nm]", rotation=90)

        pylab.savefig("article0015-fig4.png")

    # Fig.5: simulation snapshots.
    if "fig5" in sys.argv:

        # Use TeX formatting throughout.
        pylab.rc('text', usetex=True)

        # We will be using the Culgi module here, that's why we import only in this clause.
        from task0013 import loadpath

        # The number of columns and rows is important, like in Figs.2-3.
        ncols, nrows = 3, 4

        # These are the snapshots chosen for this figure.
        # This could be smarter.
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

        # Load the images chosen.
        imgs = [[misc.imread("../"+fn)[25:-24,25:-22] for fn in ch] for ch in chosen]

        # Load the RDFs for the se images also (will they be used, though?).
        rdfs = [[numpy.load(bz2.BZ2File("../"+fn.replace(fn[-14:],".hist-radial.npy.bz2"))) for fn in ch] for ch in chosen]

        # Relative margins for this figure.
        margin_x, margin_y = 0.04, 0.03

        # Effective relative size of the root figure.
        effx, effy = 1.0 - margin_x, 1.0 - margin_y

        # Image width and height.
        img_width, img_height = 3.25, 3.25

        # Root figure dimensions.
        fig_width = ncols*img_width*(1.0+margin_x)
        fig_height = nrows*img_height*(1.0+margin_x)

        # Axes for the root figure.
        ax_main = pylab.figure(figsize=(fig_width,fig_height))

        # Plot the images in their corresponding slots now.
        # Like in Figs.2-3, without any ticks.
        for fi in range(ncols*nrows):
            x = margin_x + (effx/ncols)*(fi//nrows)
            y = effy*(nrows-1)/nrows - effy*(fi%nrows)/nrows + 0.01
            a = pylab.axes([x, y, effx/ncols-0.01, effy/nrows-0.005])
            pylab.setp(a, xticks=[], yticks=[])
            pylab.imshow(imgs[fi//nrows][fi%nrows], cmap='gray')

        # Add times in text to the left side of image table.
        times = ["0", "45", "320", "1100"]
        for i in range(nrows):
            ax_main.text(margin_x/3, (1-(1.0*i+0.5)/nrows)*effy, times[i], rotation=90, fontsize=20, weight="bold")

        # Add concentrations in text to top of image table.
        conc = [700, 400, 200]
        for i in range(ncols):
            ax_main.text(margin_x+(i+0.35)*effx/ncols, effy+margin_y/4, "%i NPs" %conc[i], fontsize=20, weight="bold")

        pylab.savefig("article0015-fig5.png")

    # Fig.6: comparison between experimental and simulation RDFs.
    if "fig6" in sys.argv:

        # The chosen experimental images and simulation archive.
        chosen_exp = "../exp/sem/5k/C6/48-Image5-500.jpg"
        chosem_sim = "../phase17/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop400_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4.hist-radial.npy.bz2"

        # Load the SEM object and analyze it.
        sem = loadsem(chosen_exp)
        sem.fitrdf()

        # Get the RDF and normalize it. Some parameters are hardcoded here, because there is only one case.
        radial = 1.0 * numpy.load(bz2.BZ2File(chosem_sim))[:,0,:]
        N = 512
        size = (64,64)
        nbins = radial.shape[1]
        xmin, xmax = 0.0, 16.0
        dx = (xmax-xmin) / nbins
        xrange = numpy.arange(xmin+dx/2.0,xmax+dx/2.0,dx)
        radial = normalize_rdf(radial, size[0], size[1], xrange, N)

        # Which RDF to actually plot from simulation?
        iplot = -1

        # Match the peak position from simulation and experiment.
        xopt = xrange[radial[iplot].argmax()-1]
        xrange *= sem.rdfmax / xopt

        # Plot the experimental RDF and the corresponding fit.
        pylab.plot(sem.X, sem.Y, color="black", label="experimental", linewidth=1.5)
        pylab.plot(sem.X, sem.rdfmodel(sem.X), color="black", label="exp. fit", linewidth=1.5)

        # Plot the simulation RDF (already matched in the domain).
        pylab.plot(xrange, numpy.sum(radial[iplot-25:iplot], axis=0)/25, color="gray", label="simulation", linewidth=1.5)

        # Limits, labels and grid.
        pylab.xlabel("NP-NP distance [nm]")
        pylab.ylabel("g(r)")
        pylab.xlim([5,150])
        pylab.ylim([0,9])
        pylab.grid()

        # Plot the full experimental RDF inside an inset.
        ax = pylab.axes([0.45, 0.45, 0.4, 0.4], axisbg='lightgray')
        pylab.plot(sem.X, sem.Y, color="black")
        pylab.xticks([100, 500])
        pylab.yticks([1, 3 ,5])
        pylab.grid()

        pylab.savefig("article0015-fig6.png")

    # Fig.7: comparison of RDF peak trends for experiments and simulations.
    if "fig7" in sys.argv:

        # Always interpret text as TeX.
        pylab.rc('text', usetex=True)

        # The experimental results can be read from file after doing Fig.4,
        #   using the results for Fig.2 images (first part).
        lines = open("article0015.results.txt").readlines()
        sigma_avg = map(int,lines[1].split())
        sigma_std = map(int,lines[2].split())
        pp_avg = map(float,lines[3].split())
        pp_std = map(float,lines[4].split())
        ph_avg = map(float,lines[5].split())
        ph_std = map(float,lines[6].split())

        # Chosen simulations to compare against.
        phase = 17
        chosen_sim = [  "phase18/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop700_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4_disp3.5.hist-radial.npy.bz2",
                        "phase18/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop400_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4_disp3.5.hist-radial.npy.bz2",
                        "phase18/64x64x1_A20B16_bv1.00_sigma1.00/temp0.01_exp0.01_den1.0_pop200_chmob1.00/k15.0_nchi24.0_cc10.0_ca1.0_cb2.0_mob0.001_a10.0/tt11000_ts0.01_np4_disp3.5.hist-radial.npy.bz2",
        ]

        # Load the RDF and normalize it.
        # We could use the normalize_rdf method, but the factors will be the same
        #   for many RDFs, so that would be a bit of an overkill.
        radials = [1.0*numpy.load(bz2.BZ2File("../"+fn))[:,0,:] for fn in chosen_sim]
        nbins = radials[0].shape[1]
        size = (64,64)
        xmin = 0.0
        xmax = 16.0
        dx = (xmax-xmin)/nbins
        xrange = numpy.arange(xmin+dx/2.0,xmax+dx/2.0,dx)
        ns = [700, 400, 200]
        npairs = [n*(n-1)/2 for n in ns]
        strips = 1.0 * rdfcorrection(size[0],size[1],xrange) * 2 * numpy.pi * xrange * dx
        factors = numpy.array([1.0 * size[0] * size[1] / strips / np for np in npairs])
        for ni,n in enumerate(ns):
            radials[ni][:] *= factors[ni]

        # Frames and number of context frames.
        frames = phases_frames[phase-1][0]
        nf = phases_frames[phase-1][1]

        # Extract the peak position as the index of the maximum value (quite crude).
        # Once known, average the peak position for the context frames (becomes less crude).
        # Finally, match the final experimental/simulation peak positions.
        D = [xrange[numpy.argmax(r[:,:250], axis=1)] for r in radials]
        D = array([[average(d[fi*nf:(fi+1)*nf]) for fi,f in enumerate(frames)] for di,d in enumerate(D)])
        D *= average(pp_avg[3::4]) / average(D[:,-1])

        # Extract the peak height as the maximum value (quite crude).
        # Again, average for context frames (less crude).
        # Do not match this with experiment -- leave raw.
        G = [numpy.max(radial[:,:250], axis=1) for radial in radials]
        G = array([[average(g[fi*nf:(fi+1)*nf]) for fi,f in enumerate(frames)] for gi,g in enumerate(G)])

        # Now set and scale the domain for the simulation data.
        X = numpy.array(phases_frames[11][0])
        X *= 48.0 / 200.0

        # Styles and labels to use.
        styles = ["ko-", "ks-", "k^-"]
        labels = ["a-d", "e-h", "i-l"]

        # The left hand side plot will be of peak position (pp).
        pylab.subplot(121)
        pylab.plot([0,2,14,48], pp_avg[:4], "ko-.", mfc='None', label="exp. %i NP/$\mu$m$^2$" %sigma_avg[0])
        pylab.plot([0,2,14,48], pp_avg[4:8], "ks-.", mfc='None', label="exp. %i NP/$\mu$m$^2$" %sigma_avg[1])
        pylab.plot([0,2,14,48], pp_avg[8:12], "k^-.", mfc='None', label="exp. %i NP/$\mu$m$^2$" %sigma_avg[2])
        for di,d in enumerate(D):
            pylab.plot(X, d, styles[di], label= "sim. %i NPs" %(ns[di]))
        pylab.xlim([0, 50])
        pylab.ylim([23, 47])
        pylab.xlabel("time [h]")
        pylab.ylabel("interparticle distance $d$ [nm]")
        pylab.legend(prop={'size':10})
        pylab.grid()

        # The right hand side plot will be of peak height (ph).
        pylab.subplot(122)
        pylab.plot([0,2,14,48], ph_avg[:4], "ko-.", mfc='None')
        pylab.plot([0,2,14,48], ph_avg[4:8], "ks-.", mfc='None')
        pylab.plot([0,2,14,48], ph_avg[8:12], "k^-.", mfc='None')
        for gi,g in enumerate(G):
            pylab.plot(X, g, styles[gi])
        pylab.xlim([0, 50])
        pylab.ylim([0, 40])
        pylab.xlabel("time [h]")
        pylab.ylabel("first $g(r)$ peak height")
        pylab.grid()

        pylab.savefig("article0015-fig7.png")
